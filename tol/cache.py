"""Polite, cached HTTP access to the upstream APIs.

Every request is keyed by (method, url, body). The first time a key is seen
the response is fetched (rate-limited, with retries) and written to the cache
directory; afterwards it is always read from the cache. In offline mode a
missing key is an error, so an offline build can only ever use cached inputs.
"""

import hashlib
import json
import os
import ssl
import time
import urllib.error
import urllib.request

USER_AGENT = "tree-of-life-data-build/0.1"
# Minimum seconds between two requests to the same host.
MIN_INTERVAL = {
    "api.opentreeoflife.org": 1.0,
    "paleobiodb.org": 1.0,
    "eutils.ncbi.nlm.nih.gov": 0.5,  # NCBI asks for <= 3 req/s without a key
}
DEFAULT_INTERVAL = 1.0


class OfflineMiss(RuntimeError):
    pass


def _key(method, url, body):
    h = hashlib.sha256()
    h.update(method.encode())
    h.update(b"\n")
    h.update(url.encode())
    h.update(b"\n")
    if body is not None:
        h.update(json.dumps(body, sort_keys=True, separators=(",", ":")).encode())
    return h.hexdigest()[:40]


def _host(url):
    return url.split("/")[2]


class Cache:
    def __init__(self, root, offline=False, verbose=True):
        self.root = root
        self.offline = offline
        self.verbose = verbose
        self._last = {}
        self.hits = 0
        self.fetched = 0
        self.used = set()  # cache files touched by this build

    def path_for(self, method, url, body):
        return os.path.join(self.root, _host(url), _key(method, url, body) + ".json")

    def get_json(self, url, body=None, ok_errors=()):
        """GET (body is None) or POST JSON (body is a dict). Returns parsed JSON."""
        method = "GET" if body is None else "POST"
        path = self.path_for(method, url, body)
        rel = os.path.relpath(path, self.root)
        if os.path.exists(path):
            self.hits += 1
            self.used.add(rel)
            with open(path, encoding="utf-8") as f:
                rec = json.load(f)
            return json.loads(rec["body"])
        if self.offline:
            raise OfflineMiss("not in cache (offline build): %s %s %s" % (method, url, body))
        text = self._fetch(method, url, body, ok_errors)
        data = json.loads(text)  # validate before caching
        os.makedirs(os.path.dirname(path), exist_ok=True)
        rec = {
            "request": {"method": method, "url": url, "body": body},
            "fetched": time.strftime("%Y-%m-%d", time.gmtime()),
            "body": text,
        }
        tmp = path + ".tmp"
        with open(tmp, "w", encoding="utf-8") as f:
            json.dump(rec, f, sort_keys=True, ensure_ascii=False)
        os.replace(tmp, path)
        self.fetched += 1
        self.used.add(rel)
        return data

    def _wait(self, host):
        gap = MIN_INTERVAL.get(host, DEFAULT_INTERVAL)
        last = self._last.get(host)
        if last is not None:
            delay = gap - (time.monotonic() - last)
            if delay > 0:
                time.sleep(delay)
        self._last[host] = time.monotonic()

    def _fetch(self, method, url, body, ok_errors=()):
        host = _host(url)
        data = None
        headers = {"User-Agent": USER_AGENT, "Accept": "application/json"}
        if body is not None:
            data = json.dumps(body, sort_keys=True).encode()
            headers["Content-Type"] = "application/json"
        ctx = ssl.create_default_context()
        last_err = None
        for attempt in range(5):
            self._wait(host)
            if self.verbose:
                print("  fetch %s %s%s" % (method, url[:160], " " + json.dumps(body)[:80] if body else ""))
            req = urllib.request.Request(url, data=data, headers=headers, method=method)
            try:
                with urllib.request.urlopen(req, timeout=180, context=ctx) as resp:
                    return resp.read().decode("utf-8")
            except urllib.error.HTTPError as e:
                last_err = e
                if e.code in ok_errors:
                    # a well-formed API error that is part of the answer (cached like a success)
                    return e.read().decode("utf-8")
                if e.code < 500 and e.code != 429:
                    raise
            except (urllib.error.URLError, TimeoutError, ConnectionError) as e:
                last_err = e
            time.sleep(5 * (attempt + 1))
        raise RuntimeError("giving up on %s: %s" % (url, last_err))
