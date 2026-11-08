"""Assemble the dataset from cached API responses.

Pipeline (each step only reads the cache through tol.sources):
  1. Open Tree taxonomy (OTT): expand the focal groups, add the lineage of every
     focal group and seed taxon up to the root, plus siblings of those ancestors.
  2. Open Tree synthetic tree: induced subtree of the selected OTT ids.
  3. PBDB: taxa above genus and genera under pinned anchor ids; match OTT clades
     to PBDB taxa by exact accepted name *within an id-anchored scope* plus a rank
     check (never by a global name search); fossil ranges and environment.
  4. PBDB specimen measurements: body length for genera (fossil specimens).
  5. NCBI Taxonomy: English common names (via OTT's ncbi: source ids).
  6. Curated, cited tables: divergence dates, sizes, habitats.
  7. Time slices (PBDB/ICS intervals) with representative genera.
"""

import hashlib
import json
import math
import os
import re
from collections import defaultdict

from . import sources

PACKAGE_VERSION = "0.1.0"

RANKS = [
    "domain", "superkingdom", "kingdom", "subkingdom", "infrakingdom",
    "superphylum", "phylum", "subphylum", "infraphylum",
    "superclass", "class", "subclass", "infraclass", "megacohort", "supercohort", "cohort", "subcohort",
    "superorder", "order", "suborder", "infraorder", "parvorder",
    "superfamily", "family", "subfamily", "tribe", "subtribe",
    "genus", "subgenus", "section", "subsection", "series", "species group", "species subgroup",
    "species", "subspecies", "variety", "forma",
]
RANK_INDEX = {r: i for i, r in enumerate(RANKS)}
UNRANKED = {"no rank", "unranked clade", "unranked", "informal", "clade", "", None}

ENV_TAGS = {
    "marine": ["marine"],
    "brackish": ["brackish"],
    "freshwater": ["freshwater"],
    "terrestrial": ["terrestrial"],
    "hypersaline": ["hypersaline"],
    "marine,freshwater": ["freshwater", "marine"],
}
ROLE_PRIORITY = ["ancestor", "focal", "seed", "seed_path", "context"]
PLACEHOLDER = re.compile(r"^(uncultured|unclassified|unidentified|environmental)\b", re.I)
BINOMIAL = re.compile(r"^[A-Z][a-z-]+ (\([A-Z][a-z]+\) )?[a-z]")
BODY_PARTS = {"body", "whole body", "entire body", "total body", "complete body", "whole"}


# ------------------------------------------------------------------ helpers
def rank_check(ott_rank, pbdb_rank):
    """Return 'exact', 'near', 'unranked' or None (reject)."""
    a = (ott_rank or "").strip().lower()
    b = (pbdb_rank or "").strip().lower()
    if a in UNRANKED or b in UNRANKED:
        return "unranked"
    if a == b:
        return "exact"
    if a in RANK_INDEX and b in RANK_INDEX and abs(RANK_INDEX[a] - RANK_INDEX[b]) <= 3:
        return "near"
    return None


def num(x):
    if x is None or x == "":
        return None
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    if math.isnan(v):
        return None
    return int(v) if v == int(v) and abs(v) < 1e15 else round(v, 6)


def r3(x):
    if x is None:
        return None
    v = float("%.4g" % x)
    return int(v) if v == int(v) else v


def env_tags(raw):
    if not raw:
        return []
    raw = raw.strip().lower()
    if raw in ENV_TAGS:
        return list(ENV_TAGS[raw])
    tags = set()
    for part in re.split(r"[,/;]| and ", raw):
        p = part.strip()
        for key in ("marine", "brackish", "freshwater", "terrestrial", "hypersaline"):
            if key in p:
                tags.add(key)
    return sorted(tags)


def pbdb_fossil(rec):
    if rec is None or not rec.get("firstapp_max_ma") and rec.get("firstapp_max_ma") != 0:
        return None
    return {
        "first_max_ma": num(rec.get("firstapp_max_ma")),
        "first_min_ma": num(rec.get("firstapp_min_ma")),
        "last_max_ma": num(rec.get("lastapp_max_ma")),
        "last_min_ma": num(rec.get("lastapp_min_ma")),
        "early_interval": rec.get("early_interval"),
        "late_interval": rec.get("late_interval") or rec.get("early_interval"),
        "n_occs": int(rec.get("n_occs") or 0),
        "basis": "pbdb_taxon",
        "src": "PBDB",
    }


def _agg(fn, rows, key):
    vals = [r[key] for r in rows if r.get(key) is not None]
    return fn(vals) if vals else None


def dumps(obj):
    return json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":"))


def write_jsonl(path, rows):
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        for r in rows:
            f.write(dumps(r))
            f.write("\n")


def write_json(path, obj):
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(json.dumps(obj, sort_keys=True, ensure_ascii=False, indent=1))
        f.write("\n")


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        h.update(f.read())
    return h.hexdigest()


def load_json(path):
    with open(path, encoding="utf-8") as f:
        return json.load(f)


# ------------------------------------------------------------------ step 1: OTT
class Taxa:
    def __init__(self):
        self.nodes = {}  # ott_id -> dict

    def add(self, t, parent, role):
        oid = int(t["ott_id"])
        n = self.nodes.get(oid)
        if n is None:
            n = {
                "ott_id": oid,
                "name": t["name"],
                "unique_name": t.get("unique_name") or t["name"],
                "rank": t.get("rank") or "no rank",
                "flags": sorted(t.get("flags") or []),
                "tax_sources": sorted(t.get("tax_sources") or []),
                "parent": parent,
                "roles": set(),
            }
            self.nodes[oid] = n
        elif n["parent"] is None and parent is not None:
            n["parent"] = parent
        n["roles"].add(role)
        return n


def collect_taxa(cache, cfg, filtered_flags, log):
    taxa = Taxa()
    skip_ranks = set(cfg["skip_ranks"])
    synth_flags = set(filtered_flags)
    skip_flags = synth_flags | set(cfg.get("skip_flags_extra", []))
    root_id = cfg["root"]["ott_id"]

    def keep(child, context=False):
        if child.get("rank") in skip_ranks:
            return False
        # unranked binomials / trinomials are species-level names: skip them
        if (child.get("rank") or "no rank") == "no rank" and BINOMIAL.match(child["name"]):
            return False
        if (synth_flags if context else skip_flags) & set(child.get("flags") or []):
            return False
        return True

    def add_lineage(info, role):
        """info: taxon_info with lineage. Adds the path root..info."""
        lin = info["lineage"]  # nearest parent first
        ids = [int(x["ott_id"]) for x in lin]
        if root_id not in ids:
            raise RuntimeError("%s (ott %s) is not under the root" % (info["name"], info["ott_id"]))
        path = lin[: ids.index(root_id) + 1][::-1]  # root first
        parent = None
        for x in path:
            taxa.add(x, parent, "ancestor")
            parent = int(x["ott_id"])
        taxa.add(info, parent, role)
        return [int(x["ott_id"]) for x in path]

    ancestors = set()
    focal_roots = {}
    for g in cfg["focal_groups"]:
        info = sources.ott_taxon(cache, g["ott_id"], lineage=True)
        if info["name"] != g["name"]:
            raise RuntimeError("OTT %s is %r, expected %r" % (g["ott_id"], info["name"], g["name"]))
        ancestors.update(add_lineage(info, "focal"))
        focal_roots[int(g["ott_id"])] = g

    # depth-limited expansion of each focal group
    for g in cfg["focal_groups"]:
        stop = RANK_INDEX[g["stop_rank"]]
        frontier = [(int(g["ott_id"]), 0)]
        while frontier:
            oid, depth = frontier.pop(0)
            if depth >= g["max_depth"]:
                continue
            info = sources.ott_taxon(cache, oid, children=True)
            for ch in sorted(info.get("children", []), key=lambda c: int(c["ott_id"])):
                if not keep(ch) or PLACEHOLDER.match(ch["name"]):
                    continue
                r = ch.get("rank") or "no rank"
                ri = RANK_INDEX.get(r)
                if ri is not None and ri > stop:
                    continue
                taxa.add(ch, oid, "focal")
                if ri is None or ri < stop:
                    frontier.append((int(ch["ott_id"]), depth + 1))
        log("  ott: %s expanded, %d taxa so far" % (g["name"], len(taxa.nodes)))

    # seeds (representative species) and their lineages
    for s in cfg["seed_taxa"]:
        info = sources.ott_taxon(cache, s["ott_id"], lineage=True)
        if info["name"] != s["name"]:
            raise RuntimeError("OTT %s is %r, expected %r" % (s["ott_id"], info["name"], s["name"]))
        for oid in add_lineage(info, "seed"):
            n = taxa.nodes[oid]
            if oid not in ancestors and "focal" not in n["roles"]:
                n["roles"].add("seed_path")

    # siblings of the ancestors of the focal groups (one level, no expansion)
    if cfg.get("include_ancestor_siblings"):
        for oid in sorted(ancestors):
            if oid in focal_roots:
                continue
            info = sources.ott_taxon(cache, oid, children=True)
            for ch in sorted(info.get("children", []), key=lambda c: int(c["ott_id"])):
                if keep(ch, context=True) and not PLACEHOLDER.match(ch["name"]):
                    taxa.add(ch, oid, "context")
        log("  ott: ancestor siblings added, %d taxa" % len(taxa.nodes))

    for n in taxa.nodes.values():
        if "ancestor" in n["roles"] and int(n["ott_id"]) not in ancestors:
            # only on a seed path
            n["roles"].discard("ancestor")
            n["roles"].add("seed_path")
        if int(n["ott_id"]) in ancestors and int(n["ott_id"]) not in focal_roots:
            n["roles"] = {"ancestor"} | (n["roles"] & {"seed"})
    if len(taxa.nodes) > cfg["max_clades"]:
        raise RuntimeError("scope too large: %d clades > max_clades %d" % (len(taxa.nodes), cfg["max_clades"]))
    return taxa, focal_roots, ancestors


# ------------------------------------------------------------------ step 2: synthetic tree
def parse_newick(s):
    """Minimal Newick parser for unquoted labels. Returns {label: parent_label}."""
    s = s.strip().rstrip(";")
    parents = {}
    stack = []
    i = 0
    pending_children = []  # list of lists of labels of nodes whose parent is not yet known
    pending_children.append([])
    label_re = re.compile(r"[^(),:;\[\]]*")
    while i < len(s):
        c = s[i]
        if c == "(":
            stack.append(pending_children)
            pending_children = [[]]
            i += 1
        elif c == ",":
            i += 1
        elif c == ")":
            i += 1
            m = label_re.match(s, i)
            label = m.group(0).strip()
            i = m.end()
            if i < len(s) and s[i] == ":":
                m2 = re.compile(r":[^(),;]*").match(s, i)
                i = m2.end()
            for ch in pending_children[0]:
                parents[ch] = label
            pending_children = stack.pop()
            pending_children[0].append(label)
        else:
            m = label_re.match(s, i)
            label = m.group(0).strip()
            i = m.end()
            if i < len(s) and s[i] == ":":
                m2 = re.compile(r":[^(),;]*").match(s, i)
                i = m2.end()
            if label:
                pending_children[0].append(label)
            elif m.end() == m.start():
                i += 1
    for ch in pending_children[0]:
        parents.setdefault(ch, None)
    return parents


def synth_relations(cache, taxa, log):
    ids = sorted(taxa.nodes)
    res, dropped = sources.ot_induced_subtree(cache, ids)
    newick = res["newick"]
    broken = {int(k[3:]): v for k, v in (res.get("broken") or {}).items()}
    parents = parse_newick(newick)
    in_tree = {int(l[3:]) for l in parents if re.fullmatch(r"ott\d+", l)}
    rel = {}
    for oid in ids:
        if oid in in_tree:
            label, status = "ott%d" % oid, "in_tree"
        elif oid in broken:
            label, status = broken[oid], "broken"
        else:
            rel[oid] = {"status": "not_in_tree", "synth_parent": None, "node": None,
                        "reason": dropped.get(oid)}
            continue
        p = parents.get(label)
        while p is not None and not (re.fullmatch(r"ott\d+", p) and int(p[3:]) in taxa.nodes):
            p = parents.get(p)
        rel[oid] = {
            "status": status,
            "synth_parent": ("ott:" + p[3:]) if p else None,
            "node": label,
        }
    log("  synth: %d in tree, %d broken, %d not in tree" % (
        sum(1 for r in rel.values() if r["status"] == "in_tree"),
        sum(1 for r in rel.values() if r["status"] == "broken"),
        sum(1 for r in rel.values() if r["status"] == "not_in_tree")))
    return rel, newick, res.get("broken") or {}


# ------------------------------------------------------------------ step 3/4: PBDB
class PbdbIndex:
    """PBDB taxa from one id-anchored scope, indexed by accepted id and accepted name."""

    def __init__(self):
        self.by_id = {}  # accepted_no -> canonical record
        self.alias = {}  # any taxon_no / orig_no -> accepted_no
        self.by_name = defaultdict(set)

    def add(self, rec):
        acc = int(rec["accepted_no"])
        for k in ("taxon_no", "orig_no"):
            if rec.get(k):
                self.alias[int(rec[k])] = acc
        self.alias[acc] = acc
        canonical = rec.get("taxon_name", rec.get("accepted_name")) == rec.get("accepted_name") and not rec.get("difference")
        cur = self.by_id.get(acc)
        if cur is None or (canonical and not cur["_canonical"]) or (
                canonical == cur["_canonical"] and int(rec["taxon_no"]) < int(cur["taxon_no"])):
            r = dict(rec)
            r["_canonical"] = canonical
            self.by_id[acc] = r
        self.by_name[rec["accepted_name"]].add(acc)

    def candidates(self, name):
        return sorted(self.by_name.get(name, ()))

    def resolve(self, any_no):
        return self.alias.get(int(any_no))


def fetch_pbdb(cache, cfg, log):
    show = "app,ecospace,etbasis,classext"
    anchors = {}
    scopes = {}
    genera = {}
    measurements = {}
    for g in cfg["focal_groups"]:
        a = g["pbdb_anchor"]
        aid = int(a["id"])
        if aid not in anchors:
            rec = sources.pbdb(cache, "taxa/list", {"id": aid, "show": show})["records"][0]
            if rec["accepted_name"] != a["name"] or int(rec["accepted_no"]) != aid:
                raise RuntimeError("PBDB anchor %s is %r, expected %r" % (aid, rec["accepted_name"], a["name"]))
            anchors[aid] = rec
            idx = PbdbIndex()
            for r in sources.pbdb_paged(cache, "taxa/list", {"base_id": aid, "rank": "above_genus", "show": show}):
                idx.add(r)
            scopes[aid] = idx
        if g.get("pbdb_genera") and aid not in genera:
            genera[aid] = sources.pbdb_paged(cache, "taxa/list", {"base_id": aid, "rank": "genus", "show": show})
            for r in genera[aid]:
                scopes[aid].add(r)
        if g.get("pbdb_measurements") and aid not in measurements:
            measurements[aid] = sources.pbdb_paged(
                cache, "specs/measurements", {"base_id": aid, "show": "spec,classext"}, page=2500)
        log("  pbdb: %s anchor %s: %d scope taxa, %d genera, %d measurements" % (
            g["name"], aid, len(scopes[aid].by_id), len(genera.get(aid, [])), len(measurements.get(aid, []))))

    # the scope for ancestors / context: every PBDB parent of the anchors, and their children
    ids = ",".join(str(i) for i in sorted(anchors))
    parents = sources.pbdb(cache, "taxa/list", {"id": ids, "rel": "all_parents", "show": show})["records"]
    parent_ids = sorted({int(r["accepted_no"]) for r in parents})
    lineage_idx = PbdbIndex()
    for r in parents:
        lineage_idx.add(r)
    for i in range(0, len(parent_ids), 40):
        chunk = ",".join(str(x) for x in parent_ids[i:i + 40])
        for r in sources.pbdb(cache, "taxa/list", {"id": chunk, "rel": "children", "show": show})["records"]:
            lineage_idx.add(r)
    # pinned ids
    pins = {}
    for k, v in cfg.get("pbdb_pins", {}).items():
        if k.startswith("_"):
            continue
        rec = sources.pbdb(cache, "taxa/list", {"id": v["id"], "show": show})["records"][0]
        if rec["accepted_name"] != v["name"]:
            raise RuntimeError("PBDB pin %s is %r, expected %r" % (v["id"], rec["accepted_name"], v["name"]))
        pins[int(k)] = rec
    log("  pbdb: lineage scope %d taxa, %d pins" % (len(lineage_idx.by_id), len(pins)))
    return anchors, scopes, genera, measurements, lineage_idx, pins


def match_pbdb(taxa, focal_roots, cfg, scopes, lineage_idx, pins):
    """OTT clade -> PBDB record, by exact accepted name within an id-anchored scope + rank check."""
    group_anchor = {int(g["ott_id"]): int(g["pbdb_anchor"]["id"]) for g in cfg["focal_groups"]}
    matches = {}
    rejected = []
    for oid in sorted(taxa.nodes):
        n = taxa.nodes[oid]
        if oid in pins:
            rec = pins[oid]
            matches[oid] = (rec, "pinned", rank_check(n["rank"], rec.get("accepted_rank")) or "pinned")
            continue
        scope_list = []
        g = n.get("focal_group")
        if g is not None:
            scope_list.append(scopes[group_anchor[g]])
        scope_list.append(lineage_idx)
        for idx in scope_list:
            cands = idx.candidates(n["name"])
            if not cands:
                continue
            if len(cands) > 1:
                rejected.append({"ott_id": oid, "name": n["name"], "reason": "ambiguous in scope", "pbdb": cands})
                break
            rec = idx.by_id[cands[0]]
            rc = rank_check(n["rank"], rec.get("accepted_rank"))
            if rc is None:
                rejected.append({"ott_id": oid, "name": n["name"], "reason": "rank mismatch %s vs %s" % (
                    n["rank"], rec.get("accepted_rank")), "pbdb": cands})
                break
            matches[oid] = (rec, "exact_name_in_scope", rc)
            break
    return matches, rejected


# ------------------------------------------------------------------ main build
def build(cache, root, out_dir, log=print):
    cfg = load_json(os.path.join(root, "config", "scope.json"))
    refs = load_json(os.path.join(root, "curated", "references.json"))
    cur_div = load_json(os.path.join(root, "curated", "divergence.json"))
    cur_size = load_json(os.path.join(root, "curated", "sizes.json"))
    cur_hab = load_json(os.path.join(root, "curated", "habitat.json"))
    corrections = load_json(os.path.join(root, "curated", "pbdb_corrections.json"))["genus_first_appearance"]
    bad_sigs = {(c["bad"]["first_max_ma"], c["bad"]["early_interval"]) for c in corrections.values()}

    about = sources.ot_about(cache)
    tax_about = sources.ott_about(cache)
    log("open tree synthesis %s, taxonomy %s" % (about["synth_id"], tax_about["source"]))

    taxa, focal_roots, ancestors = collect_taxa(cache, cfg, about.get("filtered_flags", []), log)

    # children lists and focal group membership (nearest focal root on the path)
    children = defaultdict(list)
    for oid, n in taxa.nodes.items():
        if n["parent"] is not None:
            children[n["parent"]].append(oid)
    for oid, n in taxa.nodes.items():
        g, p = None, oid
        while p is not None:
            if p in focal_roots:
                g = p
                break
            p = taxa.nodes[p]["parent"]
        n["focal_group"] = g

    def depth(oid):
        d = 0
        while taxa.nodes[oid]["parent"] is not None:
            oid = taxa.nodes[oid]["parent"]
            d += 1
        return d

    synth, newick, broken = synth_relations(cache, taxa, log)

    anchors, scopes, genera_raw, meas_raw, lineage_idx, pins = fetch_pbdb(cache, cfg, log)
    matches, rejected = match_pbdb(taxa, focal_roots, cfg, scopes, lineage_idx, pins)
    pbdb_to_ott = defaultdict(set)
    for oid, (rec, how, rc) in matches.items():
        pbdb_to_ott[int(rec["accepted_no"])].add(oid)

    # ---- genera: assign each to the nearest matched clade by walking PBDB parents
    group_anchor = {int(g["ott_id"]): int(g["pbdb_anchor"]["id"]) for g in cfg["focal_groups"]}
    anchor_root = {}
    for g in cfg["focal_groups"]:
        a = int(g["pbdb_anchor"]["id"])
        if g.get("pbdb_genera"):
            anchor_root.setdefault(a, int(g["ott_id"]))

    def clade_for_pbdb(aid, rec):
        idx = scopes[aid]
        seen = set()
        p = rec.get("parent_no")
        while p:
            acc = idx.resolve(p)
            if acc is None or acc in seen:
                break
            seen.add(acc)
            hits = [o for o in pbdb_to_ott.get(acc, ()) if group_anchor.get(taxa.nodes[o]["focal_group"]) == aid]
            if hits:
                return max(hits, key=lambda o: (depth(o), -o))
            nxt = idx.by_id.get(acc, {}).get("parent_no")
            p = nxt
        return anchor_root[aid]

    min_occs = cfg["genera"]["min_occs"]
    genera = {}
    for aid in sorted(genera_raw):
        idx = scopes[aid]
        for acc in sorted({int(r["accepted_no"]) for r in genera_raw[aid]}):
            rec = idx.by_id[acc]
            if rec.get("accepted_rank") != "genus" or rec["accepted_name"] != rec.get("taxon_name", rec["accepted_name"]):
                continue
            if int(rec.get("n_occs") or 0) < min_occs:
                continue
            if str(acc) in corrections:
                c = corrections[str(acc)]
                if (num(rec.get("firstapp_max_ma")), rec.get("early_interval")) != (c["bad"]["first_max_ma"], c["bad"]["early_interval"]):
                    raise RuntimeError("PBDB correction for %s no longer applies; review curated/pbdb_corrections.json" % acc)
                rec = dict(rec, firstapp_max_ma=c["corrected"]["first_max_ma"], firstapp_min_ma=c["corrected"]["first_min_ma"],
                           early_interval=c["corrected"]["early_interval"], _corrected=True)
            clade = clade_for_pbdb(aid, rec)
            genera[acc] = {"rec": rec, "clade": clade, "anchor": aid}

    # ---- body length from PBDB specimen measurements, per genus
    gsize = defaultdict(list)
    parts_seen = defaultdict(int)
    for aid in sorted(meas_raw):
        idx = scopes[aid]
        for m in meas_raw[aid]:
            if (m.get("measurement_type") or "").lower() != "length":
                continue
            part = (m.get("specimen_part") or "").strip().lower()
            parts_seen[part] += 1
            if part not in BODY_PARTS:
                continue
            gno = m.get("genus_no")
            if not gno:
                continue
            acc = idx.resolve(gno)
            if acc is None or acc not in genera:
                continue
            for k in ("min", "average", "max"):
                v = num(m.get(k))
                if v is not None and v > 0:
                    gsize[acc].append(float(v))
    for acc, vals in gsize.items():
        vals.sort()
        med = vals[len(vals) // 2] if len(vals) % 2 else (vals[len(vals) // 2 - 1] + vals[len(vals) // 2]) / 2
        genera[acc]["size"] = {
            "min_mm": r3(vals[0]), "max_mm": r3(vals[-1]), "median_mm": r3(med), "n": len(vals),
            "basis": "pbdb_measurements", "src": ["PBDB"],
            "note": "Body length of fossil specimens (PBDB specimen measurements, specimen_part=body). Units: mm.",
        }

    # ---- NCBI common names
    ncbi_ids = {}
    for oid, n in taxa.nodes.items():
        for s in n["tax_sources"]:
            if s.startswith("ncbi:"):
                ncbi_ids[oid] = s[5:]
                break
    summaries = sources.ncbi_summaries(cache, ncbi_ids.values())

    # ---- assemble clade records
    clades = {}
    for oid in sorted(taxa.nodes):
        n = taxa.nodes[oid]
        rec = {
            "id": "ott:%d" % oid,
            "ott_id": oid,
            "name": n["name"],
            "unique_name": n["unique_name"],
            "rank": n["rank"],
            "parent": ("ott:%d" % n["parent"]) if n["parent"] is not None else None,
            "children": sorted("ott:%d" % c for c in children.get(oid, [])),
            "role": next(r for r in ROLE_PRIORITY if r in n["roles"]),
            "is_focal_root": oid in focal_roots,
            "focal_group": ("ott:%d" % n["focal_group"]) if n["focal_group"] is not None else None,
            "depth": depth(oid),
            "ott_flags": n["flags"],
            "extinct": bool({"extinct", "extinct_inherited"} & set(n["flags"])),
            "xref": {"ncbi": ncbi_ids.get(oid)},
            "synth": synth[oid],
            "common_name": None,
            "pbdb": None,
            "fossil": None,
            "habitat": [],
            "habitat_src": None,
            "ecology": None,
            "divergence": [],
            "size_mm": None,
        }
        s = summaries.get(ncbi_ids.get(oid) or "")
        if s and s.get("commonname") and s.get("scientificname") == n["name"]:
            rec["common_name"] = {"en": s["commonname"], "src": "NCBI Taxonomy"}
        m = matches.get(oid)
        if m:
            prec, how, rc = m
            rec["pbdb"] = {"id": "txn:%s" % prec["accepted_no"], "name": prec["accepted_name"],
                           "rank": prec.get("accepted_rank"), "match": how, "rank_check": rc}
            rec["fossil"] = pbdb_fossil(prec)
            env = prec.get("taxon_environment")
            if env:
                rec["habitat"] = env_tags(env)
                rec["habitat_src"] = {"src": "PBDB", "pbdb_environment": env,
                                      "inherited_from": prec.get("environment_basis")}
            eco = {k: prec[k] for k in ("motility", "life_habit", "diet") if prec.get(k)}
            rec["ecology"] = ({"src": "PBDB", **eco}) if eco else None
        key = str(oid)
        if key in cur_div:
            rec["divergence"] = sorted(cur_div[key], key=lambda e: (e["kind"], e["src"]))
        if key in cur_hab and not rec["habitat"]:
            h = cur_hab[key]
            rec["habitat"] = sorted(h["tags"])
            rec["habitat_src"] = {"src": "curated", "ref": h["src"], "check": h["check"]}
        if key in cur_size:
            c = cur_size[key]
            rec["size_mm"] = {"min_mm": c["min_mm"], "max_mm": c["max_mm"], "basis": "curated",
                              "src": [c["src"]], "check": c["check"], "note": c.get("note")}
        clades[oid] = rec

    # ---- bottom-up derivations: fossil range and size from descendants and assigned genera
    genera_by_clade = defaultdict(list)
    for acc, g in genera.items():
        genera_by_clade[g["clade"]].append(acc)
    order = sorted(clades, key=lambda o: (-depth(o), o))
    for oid in order:  # habitat: union of members' tags where PBDB/curation gives none
        if clades[oid]["habitat"] or clades[oid]["focal_group"] is None:
            continue  # never generalise from this narrow scope up to ancestors / context
        tags = set()
        for acc in genera_by_clade.get(oid, []):
            tags.update(env_tags(genera[acc]["rec"].get("taxon_environment")))
        for c in children.get(oid, []):
            tags.update(clades[c]["habitat"])
        if tags:
            clades[oid]["habitat"] = sorted(tags)
            clades[oid]["habitat_src"] = {"src": "derived_from_descendants"}
    agg_f = {}
    agg_s = {}
    for oid in order:
        fs, ss = [], []
        for acc in genera_by_clade.get(oid, []):
            f = pbdb_fossil(genera[acc]["rec"])
            if f:
                fs.append(f)
            if genera[acc].get("size"):
                ss.append(genera[acc]["size"])
        for c in children.get(oid, []):
            if agg_f.get(c):
                fs.append(agg_f[c])
            if agg_s.get(c):
                ss.append(agg_s[c])
        own_f = clades[oid]["fossil"]
        if own_f and (own_f["first_max_ma"], own_f["early_interval"]) in bad_sigs and fs:
            firsts = [f for f in fs if f.get("first_max_ma") is not None]
            if firsts:
                top = max(firsts, key=lambda f: (f["first_max_ma"], f.get("first_min_ma") or 0))
                own_f.update({"first_max_ma": top["first_max_ma"], "first_min_ma": top["first_min_ma"],
                              "early_interval": top.get("early_interval"), "basis": "pbdb_taxon_adjusted",
                              "note": "PBDB first appearance came from a known misassigned occurrence "
                                      "(curated/pbdb_corrections.json); recomputed from genera and child clades."})
        if own_f is None and fs:
            clades[oid]["fossil"] = {
                "first_max_ma": _agg(max, fs, "first_max_ma"),
                "first_min_ma": _agg(max, fs, "first_min_ma"),
                "last_max_ma": _agg(min, fs, "last_max_ma"),
                "last_min_ma": _agg(min, fs, "last_min_ma"),
                "early_interval": None, "late_interval": None,
                "n_occs": sum(f["n_occs"] for f in fs),
                "basis": "derived_from_descendants", "src": "PBDB",
            }
        agg_f[oid] = clades[oid]["fossil"]
        ss = [s for s in ss if s.get("basis") != "curated" or s.get("min_mm") is not None]
        if clades[oid]["size_mm"] is None and ss and clades[oid]["focal_group"] is not None:
            lo = [s["min_mm"] for s in ss if s.get("min_mm") is not None]
            hi = [s["max_mm"] for s in ss if s.get("max_mm") is not None]
            if lo and hi:
                clades[oid]["size_mm"] = {
                    "min_mm": min(lo), "max_mm": max(hi), "n": sum(s.get("n", 0) for s in ss),
                    "basis": "derived_from_descendants", "src": ["PBDB"],
                    "note": "Range of fossil body lengths of measured genera in this clade; not a typical size.",
                }
        agg_s[oid] = clades[oid]["size_mm"] if (clades[oid]["size_mm"] or {}).get("basis") != "curated" else None

    # ---- time range used by the query
    def time_range(fossil, divergence, extinct, child_starts):
        """start: oldest of (PBDB fossil first appearance [focal clades only],
        curated divergence estimate, oldest child time_range start)."""
        cands = []
        if fossil and fossil.get("first_max_ma") is not None:
            cands.append((fossil["first_max_ma"], "fossil"))
        for e in divergence:
            d = e.get("ma") if e.get("ma") is not None else (e.get("max_ma") if e.get("max_ma") is not None else e.get("min_ma"))
            if d is not None:
                cands.append((d, "divergence"))
        cands += [(c, "descendants") for c in child_starts]
        if not cands:
            return None
        start = max(c[0] for c in cands)
        basis = sorted({b for v, b in cands if v == start})
        if not extinct:
            end = 0
            basis.append("extant")
        elif fossil and fossil.get("last_min_ma") is not None:
            end = fossil["last_min_ma"]
        else:
            return None
        return {"start_ma": start, "end_ma": end, "basis": "+".join(basis)}

    caveat = ("PBDB taxon-level range for a high-level clade outside the focal groups; it can include problematic "
              "or misassigned occurrences (PBDB files Paleoproterozoic 'Problematica' under Animalia, for example). "
              "Not used for time_range.")
    for oid in order:  # deepest first
        rec = clades[oid]
        focal = rec["focal_group"] is not None
        if rec["fossil"] and not focal:
            rec["fossil"]["caveat"] = caveat
        child_starts = [clades[c]["time_range"]["start_ma"] for c in children.get(oid, [])
                        if clades[c].get("time_range") and (clades[c]["focal_group"] is not None or clades[c]["role"] == "ancestor"
                                                            or clades[c]["role"].startswith("seed"))]
        rec["time_range"] = time_range(rec["fossil"] if focal else None, rec["divergence"], rec["extinct"],
                                       child_starts if not focal else [])

    # ---- genus records
    genus_rows = []
    for acc in sorted(genera):
        g = genera[acc]
        rec = g["rec"]
        env = rec.get("taxon_environment")
        f = pbdb_fossil(rec)
        extant = rec.get("is_extant") == "extant"
        tr = None
        if f:
            tr = {"start_ma": f["first_max_ma"], "end_ma": 0 if extant else f["last_min_ma"],
                  "basis": "fossil" + ("+extant" if extant else "")}
        eco = {k: rec[k] for k in ("motility", "life_habit", "diet") if rec.get(k)}
        if f:
            f = {k: v for k, v in f.items() if k not in ("basis", "src")}
        size = g.get("size")
        if size:
            size = {k: v for k, v in size.items() if k not in ("note", "src")}
        genus_rows.append({
            "id": "pbdb:txn:%d" % acc,
            "pbdb_id": acc,
            "name": rec["accepted_name"],
            "rank": "genus",
            "clade": "ott:%d" % g["clade"],
            "focal_group": clades[g["clade"]]["focal_group"],
            "extant": extant if rec.get("is_extant") else None,
            "fossil": f,
            "time_range": tr,
            "habitat": env_tags(env),
            "pbdb_environment": env,
            "pbdb_environment_basis": rec.get("environment_basis"),
            "pbdb_corrected": bool(rec.get("_corrected")),
            "ecology": eco or None,
            "size_mm": size,
            "higher": {k: rec.get(k) for k in ("phylum", "class", "order", "family") if rec.get(k)},
        })

    # ---- time slices
    sc = cfg["slices"]
    ints = sources.pbdb(cache, "intervals/list", {"scale_id": 1})["records"]
    slices = []
    for it in ints:
        b, t = num(it["b_age"]), num(it["t_age"])
        if it["type"] == sc["phanerozoic_type"] or (it["type"] == sc["precambrian_type"] and t >= 538.8 and b <= sc["oldest_ma"]):
            slices.append({"id": "pbdb:int:%s" % it["interval_no"], "name": it["interval_name"], "type": it["type"],
                           "start_ma": b, "end_ma": t})
    slices.sort(key=lambda s: (-s["start_ma"], s["id"]))
    k_group = sc["representatives_per_group"]
    k_size = sc["representatives_per_size_decade"]
    for s in slices:
        live = [g for g in genus_rows if g["fossil"] and g["fossil"]["first_max_ma"] is not None
                and g["fossil"]["last_min_ma"] is not None
                and g["fossil"]["first_max_ma"] > s["end_ma"] and g["fossil"]["last_min_ma"] < s["start_ma"]]
        rank_key = lambda g: (-g["fossil"]["n_occs"], g["name"], g["pbdb_id"])
        by_group = defaultdict(list)
        for g in live:
            by_group[g["focal_group"] or g["clade"]].append(g)
        s["representatives"] = {k: [g["id"] for g in sorted(v, key=rank_key)[:k_group]] for k, v in sorted(by_group.items())}
        by_size = defaultdict(list)
        for g in live:
            if g["size_mm"] and g["size_mm"].get("median_mm"):
                by_size[str(int(math.floor(math.log10(g["size_mm"]["median_mm"]))))].append(g)
        s["by_size_decade"] = {k: [g["id"] for g in sorted(v, key=rank_key)[:k_size]] for k, v in sorted(by_size.items())}
        s["n_genera"] = len(live)
        s["clades_with_fossils"] = sorted(
            c["id"] for c in clades.values()
            if c["fossil"] and c["focal_group"] and c["fossil"]["basis"].startswith("pbdb_taxon")
            and c["fossil"]["first_max_ma"] is not None
            and c["fossil"]["last_min_ma"] is not None
            and c["fossil"]["first_max_ma"] > s["end_ma"] and c["fossil"]["last_min_ma"] < s["start_ma"])
        s["src"] = "PBDB intervals (ICS); genera ranked by PBDB occurrence count"

    # ---- write
    os.makedirs(out_dir, exist_ok=True)
    clade_rows = [clades[o] for o in sorted(clades)]
    write_jsonl(os.path.join(out_dir, "clades.jsonl"), clade_rows)
    write_jsonl(os.path.join(out_dir, "genera.jsonl"), genus_rows)
    write_jsonl(os.path.join(out_dir, "slices.jsonl"), slices)
    with open(os.path.join(out_dir, "synthetic_subtree.nwk"), "w", encoding="utf-8", newline="\n") as f:
        f.write(newick.strip() + "\n")
    used_refs = sorted({e["src"] for v in cur_div.values() if isinstance(v, list) for e in v}
                       | {v["src"] for k, v in cur_size.items() if not k.startswith("_")}
                       | {v["src"] for k, v in cur_hab.items() if not k.startswith("_")})
    write_json(os.path.join(out_dir, "references.json"), {k: refs[k] for k in used_refs})

    pbdb_info = sources.pbdb_datainfo(cache)
    report = {
        "pbdb_match_rejections": sorted(rejected, key=lambda r: r["ott_id"]),
        "pbdb_specimen_parts_length": dict(sorted(parts_seen.items())),
        "synth_broken": dict(sorted(broken.items())),
    }
    write_json(os.path.join(out_dir, "build_report.json"), report)

    fetched = sorted({load_json(os.path.join(cache.root, p))["fetched"] for p in cache.used})
    manifest = {
        "package": "tree-of-life-data",
        "version": PACKAGE_VERSION,
        "counts": {
            "clades": len(clade_rows),
            "clades_with_pbdb_match": sum(1 for c in clade_rows if c["pbdb"]),
            "clades_with_fossil_range": sum(1 for c in clade_rows if c["fossil"]),
            "clades_with_size": sum(1 for c in clade_rows if c["size_mm"]),
            "clades_with_habitat": sum(1 for c in clade_rows if c["habitat"]),
            "clades_with_divergence": sum(1 for c in clade_rows if c["divergence"]),
            "clades_with_common_name": sum(1 for c in clade_rows if c["common_name"]),
            "clades_in_synthetic_tree": sum(1 for c in clade_rows if c["synth"]["status"] == "in_tree"),
            "genera": len(genus_rows),
            "genera_with_size": sum(1 for g in genus_rows if g["size_mm"]),
            "slices": len(slices),
        },
        "sources": {
            "open_tree": {"synth_id": about["synth_id"], "taxonomy": tax_about["source"],
                          "licence": "CC0 1.0", "api": sources.OT},
            "pbdb": {"api": sources.PBDB, "declared_licence": pbdb_info.get("data_license"),
                     "declared_licence_url": pbdb_info.get("license_url"),
                     "licence_applied": "CC BY 4.0 (attribution given; see DATA-LICENSE.md)"},
            "ncbi_taxonomy": {"api": sources.NCBI, "fields": ["common_name"], "licence": "public domain (US government work)"},
            "curated": {"files": ["curated/divergence.json", "curated/sizes.json", "curated/habitat.json"],
                        "licence": "CC BY 4.0"},
        },
        "fetched_dates": [fetched[0], fetched[-1]] if fetched else [],
        "n_cached_requests_used": len(cache.used),
        "files": {},
    }
    manifest["sources"]["curated"]["files"].append("curated/pbdb_corrections.json")
    for fn in ("clades.jsonl", "genera.jsonl", "slices.jsonl", "synthetic_subtree.nwk", "references.json", "build_report.json"):
        manifest["files"][fn] = {"sha256": sha256_file(os.path.join(out_dir, fn)),
                                 "bytes": os.path.getsize(os.path.join(out_dir, fn))}
    write_json(os.path.join(out_dir, "manifest.json"), manifest)
    log("wrote %d clades, %d genera, %d slices to %s" % (len(clade_rows), len(genus_rows), len(slices), out_dir))
    return manifest
