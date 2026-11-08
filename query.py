#!/usr/bin/env python3
"""Query the built dataset.

  python3 -I query.py --time 500 --size-mm 0.5 --habitat marine
  python3 -I query.py --time 407 --habitat terrestrial --limit 20
  python3 -I query.py --clade Tardigrada --depth 2        # tree around a clade
  python3 -I query.py --slice 300                         # representative genera of the slice at 300 Ma
Add --json for JSON-lines output.
"""
import argparse
import json
import os
import sys

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, ROOT)

from tol.query import Dataset  # noqa: E402


def fmt_size(s):
    if not s:
        return "size ?"
    a, b = s.get("min_mm"), s.get("max_mm")
    if a == b or b is None:
        return "%s mm" % a
    if a is None:
        return "<=%s mm" % b
    return "%s-%s mm" % (a, b)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", default=os.path.join(ROOT, "data"))
    ap.add_argument("--time", type=float, help="time T in Ma (millions of years ago)")
    ap.add_argument("--size-mm", type=float, help="body size S in mm")
    ap.add_argument("--decades", type=float, default=1.0, help="size tolerance in log10 decades (default 1)")
    ap.add_argument("--habitat", help="marine | brackish | freshwater | terrestrial")
    ap.add_argument("--include-unknown-size", action="store_true")
    ap.add_argument("--include-unknown-habitat", action="store_true")
    ap.add_argument("--kind", choices=["clade", "genus", "both"], default="both")
    ap.add_argument("--limit", type=int, default=25)
    ap.add_argument("--clade", help="name or id (ott:NNN); prints lineage and subtree")
    ap.add_argument("--depth", type=int, default=2)
    ap.add_argument("--slice", type=float, help="print the time slice containing this age (Ma)")
    ap.add_argument("--json", action="store_true")
    a = ap.parse_args()
    ds = Dataset(a.data)

    if a.clade:
        hits = ds.find(a.clade)
        if not hits:
            sys.exit("no clade named %r" % a.clade)
        for r in hits:
            if a.json:
                print(json.dumps({"clade": r, "lineage": [x["id"] for x in ds.lineage(r["id"])],
                                  "subtree": ds.subtree(r["id"], a.depth)}, sort_keys=True))
                continue
            print("%s  %s (%s)  %s" % (r["id"], r["name"], r["rank"], (r.get("common_name") or {}).get("en") or ""))
            print("  lineage: " + " < ".join(x["name"] for x in ds.lineage(r["id"])))
            f = r.get("fossil")
            if f:
                print("  fossil (PBDB): first %s-%s Ma (%s), last %s-%s Ma, %d occs [%s]" % (
                    f["first_max_ma"], f["first_min_ma"], f.get("early_interval"), f["last_max_ma"],
                    f["last_min_ma"], f["n_occs"], f["basis"]))
            for d in r.get("divergence") or []:
                print("  divergence: %s %s (%s-%s) [%s]" % (d["kind"], d.get("ma"), d.get("min_ma"), d.get("max_ma"), d["src"]))
            print("  size: %s; habitat: %s" % (fmt_size(r.get("size_mm")), ",".join(r.get("habitat") or []) or "?"))

            def show(node, ind):
                print("  " * ind + "- %s (%s)" % (node["name"], node["rank"]))
                for c in node.get("children", []):
                    show(c, ind + 1)
            show(ds.subtree(r["id"], a.depth), 1)
        return

    if a.slice is not None:
        s = ds.slice_at(a.slice)
        if not s:
            sys.exit("no slice covers %s Ma" % a.slice)
        if a.json:
            print(json.dumps(s, sort_keys=True))
            return
        print("%s (%s), %s-%s Ma, %d genera with fossils in scope" % (s["name"], s["type"], s["start_ma"], s["end_ma"], s["n_genera"]))
        for g, ids in s["representatives"].items():
            grp = ds.by_id.get(g, {}).get("name", g)
            print("  %s: %s" % (grp, ", ".join(ds.by_id[i]["name"] for i in ids)))
        for dec, ids in s["by_size_decade"].items():
            print("  size ~10^%s mm: %s" % (dec, ", ".join(ds.by_id[i]["name"] for i in ids)))
        return

    if a.time is None:
        ap.error("give --time (and optionally --size-mm/--habitat), --clade or --slice")
    kinds = ("clade", "genus") if a.kind == "both" else (a.kind,)
    res = ds.around(a.time, a.size_mm, a.habitat, a.decades, a.include_unknown_size,
                    a.include_unknown_habitat, kinds, a.limit)
    for r in res:
        if a.json:
            print(json.dumps(r, sort_keys=True))
        else:
            print("%-6s %-22s %-28s %-12s %-16s %s" % (
                r["kind"], r["id"], r["name"], r["rank"], fmt_size(r["size_mm"]), ",".join(r["habitat"])))
    if not a.json:
        print("(%d shown)" % len(res))


if __name__ == "__main__":
    main()
