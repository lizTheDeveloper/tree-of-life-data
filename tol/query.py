"""Query helpers over the built JSON-lines files (stdlib only)."""

import json
import math
import os


def load_jsonl(path):
    with open(path, encoding="utf-8") as f:
        return [json.loads(line) for line in f if line.strip()]


class Dataset:
    def __init__(self, data_dir):
        self.dir = data_dir
        self.clades = load_jsonl(os.path.join(data_dir, "clades.jsonl"))
        self.genera = load_jsonl(os.path.join(data_dir, "genera.jsonl"))
        self.slices = load_jsonl(os.path.join(data_dir, "slices.jsonl"))
        self.by_id = {r["id"]: r for r in self.clades}
        self.by_id.update({r["id"]: r for r in self.genera})

    # ------------------------------------------------------------ around()
    def around(self, time_ma, size_mm=None, habitat=None, decades=1.0,
               include_unknown_size=False, include_unknown_habitat=False,
               kinds=("clade", "genus"), limit=None):
        """What lived around time T (Ma), at body size S (mm), in habitat H.

        A record matches when
          * its time_range covers T (start_ma >= T >= end_ma);
          * if size_mm is given: its size range overlaps [S/10^d, S*10^d]
            (d = decades); records without a size match only with include_unknown_size;
          * if habitat is given: it carries that tag; records without habitat
            tags match only with include_unknown_habitat.
        Results are sorted: genera first, then clades from the most to the least
        nested, then by number of PBDB occurrences, then name.
        """
        out = []
        rows = []
        if "clade" in kinds:
            rows += [("clade", r) for r in self.clades]
        if "genus" in kinds:
            rows += [("genus", r) for r in self.genera]
        lo = hi = None
        if size_mm is not None:
            lo, hi = size_mm / 10 ** decades, size_mm * 10 ** decades
        for kind, r in rows:
            tr = r.get("time_range")
            if not tr or not (tr["start_ma"] >= time_ma >= tr["end_ma"]):
                continue
            size_known = None
            if size_mm is not None:
                s = r.get("size_mm")
                if not s or (s.get("min_mm") is None and s.get("max_mm") is None):
                    if not include_unknown_size:
                        continue
                else:
                    smin = s["min_mm"] if s.get("min_mm") is not None else s["max_mm"]
                    smax = s["max_mm"] if s.get("max_mm") is not None else s["min_mm"]
                    if smax < lo or smin > hi:
                        continue
                    size_known = True
            if habitat:
                tags = r.get("habitat") or []
                if not tags:
                    if not include_unknown_habitat:
                        continue
                elif habitat not in tags:
                    continue
            occ = (r.get("fossil") or {}).get("n_occs") or 0
            out.append({
                "kind": kind, "id": r["id"], "name": r["name"], "rank": r["rank"],
                "common_name": (r.get("common_name") or {}).get("en"),
                "time_range": tr, "size_mm": r.get("size_mm"), "habitat": r.get("habitat") or [],
                "n_occs": occ, "size_known": bool(size_known),
                "_sort": (0 if kind == "genus" else 1, -(r.get("depth") or 0), -occ, r["name"], r["id"]),
            })
        out.sort(key=lambda x: x["_sort"])
        for x in out:
            del x["_sort"]
        return out[:limit] if limit else out

    # ------------------------------------------------------------ tree helpers
    def find(self, name_or_id):
        if name_or_id in self.by_id:
            return [self.by_id[name_or_id]]
        return [r for r in self.clades if name_or_id in (r["name"], r["unique_name"])]

    def lineage(self, cid):
        out = []
        r = self.by_id.get(cid)
        while r is not None and r.get("parent"):
            r = self.by_id.get(r["parent"])
            if r:
                out.append(r)
        return out

    def subtree(self, cid, depth=2):
        r = self.by_id[cid]
        node = {"id": r["id"], "name": r["name"], "rank": r["rank"]}
        if depth > 0 and r.get("children"):
            node["children"] = [self.subtree(c, depth - 1) for c in r["children"] if c in self.by_id]
        return node

    def slice_at(self, time_ma):
        hits = [s for s in self.slices if s["start_ma"] >= time_ma >= s["end_ma"]]
        hits.sort(key=lambda s: s["start_ma"] - s["end_ma"])
        return hits[0] if hits else None


def size_decade(mm):
    return int(math.floor(math.log10(mm)))
