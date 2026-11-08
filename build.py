#!/usr/bin/env python3
"""Fetch (into cache/) and build (into data/) the tree-of-life dataset.

    python3 -I build.py              # fetch anything missing from cache, then build
    python3 -I build.py --offline    # build from cache only; fails if anything is missing
    python3 -I build.py --out /tmp/x # write somewhere else (used by the determinism check)
"""
import argparse
import os
import sys

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, ROOT)

from tol.build import build  # noqa: E402
from tol.cache import Cache  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cache", default=os.path.join(ROOT, "cache"))
    ap.add_argument("--out", default=os.path.join(ROOT, "data"))
    ap.add_argument("--offline", action="store_true", help="never touch the network")
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args()
    cache = Cache(a.cache, offline=a.offline, verbose=not a.quiet)
    log = (lambda *x: None) if a.quiet else print
    build(cache, ROOT, a.out, log=log)
    if not a.quiet:
        print("cache: %d hits, %d fetched" % (cache.hits, cache.fetched))


if __name__ == "__main__":
    main()
