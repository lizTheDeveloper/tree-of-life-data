"""Thin wrappers around the three upstream APIs. All access goes through Cache."""

import urllib.parse

OT = "https://api.opentreeoflife.org/v3"
PBDB = "https://paleobiodb.org/data1.2"
NCBI = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


# ---------------------------------------------------------------- Open Tree
def ot_about(cache):
    return cache.get_json(OT + "/tree_of_life/about", {})


def ott_about(cache):
    return cache.get_json(OT + "/taxonomy/about", {})


def ott_taxon(cache, ott_id, children=False, lineage=False):
    body = {"ott_id": int(ott_id)}
    if children:
        body["include_children"] = True
    if lineage:
        body["include_lineage"] = True
    return cache.get_json(OT + "/taxonomy/taxon_info", body)


def ot_induced_subtree(cache, ott_ids):
    """Induced synthetic subtree. Ids absent from the synthesis (pruned, extinct,
    unplaced...) make the API answer 400 with an 'unknown' map; drop them and retry.
    Returns (response, {ott_id: reason} for dropped ids)."""
    ids = sorted({int(i) for i in ott_ids})
    dropped = {}
    for _ in range(10):
        body = {"ott_ids": ids, "label_format": "id"}
        res = cache.get_json(OT + "/tree_of_life/induced_subtree", body, ok_errors=(400,))
        unknown = res.get("unknown")
        if "newick" in res:
            return res, dropped
        if not unknown:
            raise RuntimeError("induced_subtree failed: %s" % res.get("message"))
        for k, v in unknown.items():
            dropped[int(k[3:])] = v
        ids = [i for i in ids if i not in dropped]
    raise RuntimeError("induced_subtree: too many retries")


# ---------------------------------------------------------------- PBDB
def _pbdb_url(endpoint, params):
    q = dict(params)
    q.setdefault("vocab", "pbdb")
    # sorted, so the cache key does not depend on dict order
    return "%s/%s.json?%s" % (PBDB, endpoint, urllib.parse.urlencode(sorted(q.items()), safe=",:"))


def pbdb(cache, endpoint, params):
    data = cache.get_json(_pbdb_url(endpoint, params))
    if "errors" in data:
        raise RuntimeError("PBDB error for %s %s: %s" % (endpoint, params, data["errors"]))
    return data


def pbdb_paged(cache, endpoint, params, page=1000):
    """Fetch every record of a listing in small pages (offset/limit)."""
    out = []
    offset = 0
    while True:
        p = dict(params)
        p.update({"limit": page, "offset": offset})
        data = pbdb(cache, endpoint, p)
        recs = data.get("records", [])
        out.extend(recs)
        if len(recs) < page:
            return out
        offset += page


def pbdb_datainfo(cache):
    """One small request with datainfo, to record the licence PBDB declares."""
    return pbdb(cache, "taxa/single", {"id": 1, "datainfo": ""})


# ---------------------------------------------------------------- NCBI
def ncbi_summaries(cache, ids):
    ids = sorted({str(int(i)) for i in ids})
    out = {}
    for i in range(0, len(ids), 150):
        chunk = ids[i:i + 150]
        url = "%s/esummary.fcgi?%s" % (NCBI, urllib.parse.urlencode(
            [("db", "taxonomy"), ("id", ",".join(chunk)), ("retmode", "json"), ("tool", "tree-of-life-data-build")], safe=","))
        res = cache.get_json(url).get("result", {})
        for u in res.get("uids", []):
            out[u] = res[u]
    return out
