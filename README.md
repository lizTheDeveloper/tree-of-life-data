# tree-of-life-data

A small, reproducible JSON-lines dataset that joins **phylogeny** (Open Tree of
Life), **fossil time ranges and environments** (Paleobiology Database),
**common names** (NCBI Taxonomy) and a **cited table of divergence dates,
sizes and habitats**, so that you can ask:

> What lived around time *T*, at body size *S*, in habitat *H*?
> What does the tree look like around clade *X*?

Version 0.1 covers a focused slice of the eukaryote tree (see [Scope](#scope)).
It is built by a stdlib-only Python 3 script, and the same cached inputs
always give byte-identical outputs.

## Quick start

```sh
python3 -I query.py --time 500 --size-mm 0.5 --habitat marine
# clade  ott:111438  Tardigrada  phylum  0.05-1.2 mm  freshwater,marine,terrestrial

python3 -I query.py --time 310 --size-mm 100 --habitat terrestrial --limit 5   # Carboniferous giants
python3 -I query.py --time 407 --habitat terrestrial --include-unknown-size     # Rhynie-chert time
python3 -I query.py --slice 407                    # representative genera, Early Devonian
python3 -I query.py --clade Tardigrada --depth 2   # lineage, dates, size, subtree
python3 -I query.py --clade ott:150772 --json      # machine-readable
```

From Python:

```python
import sys; sys.path.insert(0, "path/to/tree-of-life-data")
from tol.query import Dataset
ds = Dataset("path/to/tree-of-life-data/data")
ds.around(407, size_mm=0.5, habitat="terrestrial")   # list of dicts
ds.lineage("ott:111438"); ds.subtree("ott:150772", depth=3); ds.slice_at(300)
```

How the query matches a record:

- **Time:** its `time_range` covers *T* (`start_ma >= T >= end_ma`).
- **Size:** if `--size-mm` is given, its `size_mm` range overlaps *S* ÷ 10^d to
  *S* × 10^d, where d is `--decades` (default 1). Records with no size match
  only with `--include-unknown-size`.
- **Habitat:** if `--habitat` is given, it carries that tag. Records with no
  habitat match only with `--include-unknown-habitat`.
- **Order:** genera first, then clades from the most to the least nested, then
  by PBDB occurrence count.

## Files (`data/`)

| File | What | Rows / size (v0.1) |
|---|---|---|
| `clades.jsonl` | one OTT clade per line | 1,138 clades, ~0.76 MB |
| `genera.jsonl` | one PBDB genus per line, with at least one fossil occurrence, in the focal groups | 15,409 genera, ~10.4 MB |
| `slices.jsonl` | ICS epochs (Phanerozoic) and periods (Proterozoic to 2500 Ma) with representative genera | 44 slices |
| `synthetic_subtree.nwk` | Open Tree synthetic tree induced on the clades (labels `ottNNN` / `mrcaottAottB`) | |
| `references.json` | full citations for every `src` key used by curated values | |
| `manifest.json` | counts, source versions, fetch dates, SHA-256 of every output | |
| `build_report.json` | PBDB matches that were rejected, and why; specimen-part counts; taxa that are broken in the synthesis | |

Every line is canonical JSON: sorted keys, no extra whitespace, UTF-8. Rows
are sorted by id.

### Clade record (`clades.jsonl`)

| Field | Meaning |
|---|---|
| `id` | stable id, `ott:<OTT id>` |
| `ott_id`, `name`, `unique_name`, `rank` | from OTT 3.7 |
| `parent`, `children` | taxonomy (OTT) parent and the children included in this package |
| `depth` | number of steps from Eukaryota |
| `synth` | where the clade sits in the Open Tree synthetic tree: `status` is `in_tree`, `broken` (not monophyletic in the synthesis; `node` is the MRCA node) or `not_in_tree` (with `reason`); `synth_parent` is the nearest included ancestor in the synthetic tree, which may differ from the taxonomy parent |
| `role` | `ancestor` (on the path from Eukaryota to a focal group), `focal` (inside a focal group), `seed` / `seed_path` (representative species and their lineage), `context` (sibling of an ancestor, one level, not expanded). `is_focal_root` marks the focal groups themselves |
| `focal_group` | nearest focal-group root above or at this clade |
| `extinct`, `ott_flags` | from OTT flags |
| `common_name` | `{en, src}` from NCBI Taxonomy when NCBI's scientific name equals the OTT name |
| `pbdb` | linked PBDB taxon: `id` (`txn:NNN`), `name`, `rank`, `match` (`exact_name_in_scope` or `pinned`), `rank_check` (`exact`, `near` (≤3 rank steps apart) or `unranked`) |
| `fossil` | PBDB first/last appearance: `first_max_ma`, `first_min_ma`, `last_max_ma`, `last_min_ma`, intervals, `n_occs`, `basis` (`pbdb_taxon`, `pbdb_taxon_adjusted` (see corrections) or `derived_from_descendants`). High-level clades outside the focal groups carry a `caveat` |
| `divergence` | list of cited estimates: `kind` (`clock`, `fossil_minimum`, `biomarker`), `ma`, `min_ma`, `max_ma`, `src` (key into `references.json`), `check`, `note` |
| `size_mm` | `{min_mm, max_mm, basis, src, ...}`; `basis` is `curated` (cited), or `derived_from_descendants` (the range of PBDB body-length measurements in the clade). It is a range, not a typical size. `null` when no licensed source exists |
| `habitat`, `habitat_src` | tags from `marine`, `brackish`, `freshwater`, `terrestrial`; source is PBDB `taxon_environment` (with `inherited_from`, the PBDB taxon the value is inherited from), a curated citation, or the union of members' tags (`derived_from_descendants`, inside focal groups only) |
| `ecology` | PBDB `motility`, `life_habit`, `diet` when present |
| `time_range` | what the query uses: `start_ma`, `end_ma`, `basis`. The start is the oldest of: the PBDB first appearance (focal clades only), a curated divergence estimate (`ma`, else the interval's older end), and the oldest child's start. The end is 0 for extant clades |
| `xref.ncbi` | NCBI taxonomy id (via OTT) |

### Genus record (`genera.jsonl`)

`id` (`pbdb:txn:NNN`), `pbdb_id`, `name`, `rank`, `clade` (the nearest
included OTT clade, found by walking PBDB parent ids to a matched clade),
`focal_group`, `extant`, `fossil`, `time_range`, `habitat`,
`pbdb_environment`, `pbdb_environment_basis` (the PBDB taxon the value is
inherited from), `ecology`, `higher` (PBDB phylum/class/order/family
names), `pbdb_corrected`, and `size_mm`. `size_mm` is the body length of
fossil specimens from PBDB specimen measurements (`specimen_part` = `body`,
`measurement_type` = `length`): `min_mm`, `max_mm`, `median_mm`, and `n`
values. Genera have no OTT id in v0.1.

### Slice record (`slices.jsonl`)

`id` (`pbdb:int:NNN`), `name`, `type`, `start_ma`, `end_ma`, `n_genera` (genera
in scope with a fossil range that overlaps the slice), `representatives` (the
top 5 genera per focal group, by PBDB occurrence count), `by_size_decade` (the
top 5 genera per log10 mm decade of median body length, so `"0"` means 1–10 mm),
and `clades_with_fossils`.

## Scope

v0.1 is Eukaryota cut down to the groups an early-life setting needs:

| Focal group (OTT id) | Depth limit | Most specific rank kept |
|---|---|---|
| Choanoflagellida (202765) | 4 | genus |
| Volvocaceae (391142) | 2 | genus |
| Tardigrada (111438) | 5 | genus |
| Nematoda (395057) | 3 | order |
| Arachnida (511967) | 3 | order |
| Acari (150772) | 6 | superfamily |
| Insecta (1062253) | 7 | order |
| Hemiptera (603650) | 6 | family |
| Aphididae (155403) | 2 | subfamily |
| Bryophyta (246594) | 5 | order |

The package also includes:

- 12 seed species (for example *Volvox carteri*, *Hypsibius exemplaris*,
  *C. elegans*, *Acyrthosiphon pisum*, *Physcomitrella patens* and the fairyfly
  *Dicopomorpha echmepterygis*), with their lineages;
- every ancestor up to Eukaryota;
- one level of siblings of those ancestors, as context: sponges, cnidarians,
  deuterostomes, Trilobitomorpha and others.

Excluded from the focal groups:

- species-rank names, and unranked binomials;
- OTT taxa flagged as unplaced or *incertae sedis*. Most are fossil genera,
  and those are covered by `genera.jsonl`;
- OTT placeholders such as "unclassified ..." and "environmental samples";
- every flag that the Open Tree synthesis itself filters (`barren`,
  `not_otu`, `major_rank_conflict`, ...).

### Widening the scope

Everything is in `config/scope.json`:

1. **Add a focal group:** add an entry with its OTT id and name, a
   `max_depth`, a `stop_rank`, and a **PBDB anchor**. The anchor is a
   hand-checked PBDB taxon id whose subtree is used for PBDB matching. It
   should be the same taxon or the nearest PBDB ancestor. Set
   `pbdb_genera` / `pbdb_measurements` to pull its genera and body lengths.
   Find OTT ids with Open Tree `tnrs/match_names` and PBDB ids with
   `taxa/list?id=...&rel=all_children`. Always confirm a candidate's lineage;
   never trust a name lookup on its own.
2. **Go deeper:** raise `max_depth` or move `stop_rank` toward genus. The OTT
   expansion makes one API call per expanded node. `max_clades` (5000) is a
   safety stop.
3. **Add a seed species:** add it to `seed_taxa`, which pulls in its lineage.
4. **Fix a name mismatch between OTT and PBDB:** add a `pbdb_pins` entry (OTT
   id → PBDB id, both checked by hand).
5. Rebuild: `python3 -I build.py`. Only missing requests are fetched.

## Rebuilding

```sh
python3 -I build.py            # fetch missing responses into cache/, then build data/
python3 -I build.py --offline  # build only from cache/; fails if a response is missing
python3 -I -m unittest discover -s tests -v
```

- **Requirements:** Python 3.8 or later, standard library only.
- **First run:** a cold build makes about 290 API requests and takes roughly
  10 minutes, because it is rate-limited to 1 request per second per host
  (NCBI: 2 per second). Most of the time goes on large PBDB pages. The cache
  is about 82 MB and is not committed.
- **Politeness:** every response is cached. Requests use small pages (1,000
  PBDB records, 2,500 for measurements), the User-Agent is
  `tree-of-life-data-build/0.1`, and failed requests are retried with backoff.
- **Determinism:** the build reads only the cache, sorts every list, writes
  canonical JSON, and contains no timestamps except the fetch dates recorded
  in the cache. `tests/` rebuilds twice offline and checks the results against
  `data/` byte for byte. It also checks `manifest.json`'s SHA-256 values.
- **Fresh fetches:** a fresh fetch reflects the live upstream databases, so
  `data/` can change. The source versions are in `manifest.json`.

### How OTT clades are linked to PBDB (and why not by name)

Many PBDB names are homonyms. Searching PBDB for "Tardigrada" returns a
junior name for **sloths** (Folivora), not the water bears. So a clade is
linked only if all of the following hold:

- the PBDB candidate lies inside an **id-anchored scope**: the subtree of the
  focal group's pinned PBDB anchor, or, for ancestors, the PBDB parents of the
  anchors and their children;
- its accepted name equals the OTT name exactly;
- it is the **only** candidate with that name in the scope;
- it passes a **rank check** (same rank, within 3 rank steps, or unranked on
  either side).

Candidates that fail a check are listed in `build_report.json`. Tardigrada
→ PBDB `txn:67137` (phylum, first appearance 506.5–497 Ma, Wuliuan) is
covered by a test.

### Corrections to upstream data

`curated/pbdb_corrections.json` lists known PBDB errors that would mislead.

- In v0.1 it has one entry. The Cretaceous wasp genus *Humiryssus* has one
  Ordovician occurrence (471 Ma) attached to it. That made Insecta,
  Pterygota and Hymenoptera appear about 470 Ma.
- The genus's first appearance is corrected from its other occurrences.
- Clades whose PBDB first appearance equals the bad value are recomputed
  from their genera and child clades (`basis: pbdb_taxon_adjusted`).
- The build fails if PBDB fixes the record, so the correction can then be
  removed.

## Sources, licences and attribution

The data licence is in [DATA-LICENSE.md](DATA-LICENSE.md). In short:

- **Code:** MIT ([LICENSE](LICENSE)).
- **Combined data:** CC BY 4.0.
- **Open Tree of Life** (topology, names, ranks, OTT ids): CC0. Synthesis
  `opentree16.1`, taxonomy `ott3.7draft3`, via
  <https://api.opentreeoflife.org/v3/>.
- **Paleobiology Database** (fossil ranges, environments, genera, specimen
  measurements, ICS intervals): treated as CC BY 4.0 and attributed, via
  <https://paleobiodb.org/data1.2/>. *Note:* on 2026-11-08 the API's
  `datainfo` declared CC0; see open questions.
- **NCBI Taxonomy** (English common names): public domain, US government work.
- **Curated tables** (`curated/`): divergence dates, sizes and habitats,
  restated with citations from the papers in `data/references.json`; CC BY 4.0.

Attribution text to reproduce:

> Tree-of-life data package. Topology and taxonomy from the Open Tree of Life
> (synthetic tree opentree16.1, OTT 3.7; CC0). Fossil ranges, environments,
> genera and specimen measurements from the Paleobiology Database
> (paleobiodb.org; CC BY 4.0). Common names from NCBI Taxonomy. Divergence
> dates, sizes and habitats from the published literature cited in the data.
> Combined data licensed CC BY 4.0.

### Deliberately not used

- **TimeTree:** its terms forbid redistribution of its data and
  transformations. No TimeTree number appears anywhere, not even as a
  fallback.
- **Heim et al. 2015** marine body sizes: no open licence could be verified
  for the Stanford Digital Repository deposit, so it is not shipped.
- **PBDB taxon "size":** in PBDB this field counts subtaxa. It is not body
  size, and it is ignored.

## Known limitations (v0.1)

- **Sizes are sparse.** Seven clades have a curated, cited size. Everything
  else comes from PBDB fossil body lengths, which are mostly insects and a
  few arachnids. Derived clade sizes are min–max ranges over measured genera,
  not typical sizes. Measurement units are not stated in the PBDB API. They
  are taken to be mm, and spot checks agree (fairyflies ~0.3 mm, fossil
  odonates ~100 mm).
- **Curated values are restated from secondary summaries.** Every curated
  date, size and habitat has `check: "secondary"`: it was restated from the
  paper's abstract or a review, and has not been re-read against the full
  text.
- **PBDB ranges for high-level clades are unreliable.** For example, PBDB
  files Paleoproterozoic "Problematica" under Animalia, giving a first
  appearance of 2050 Ma. They are shipped with a `caveat` and are not used
  for `time_range` outside the focal groups. Context clades without a
  curated date have no `time_range`, so time queries never return them.
- **PBDB environments are mostly inherited.** Almost all genera inherit
  their environment from a class or order. For 12,937 of them it is simply
  "terrestrial" from Insecta. So the Devonian sea scorpion *Palaeoscorpius*
  is tagged terrestrial through Arachnida. Use `pbdb_environment_basis` to
  tell an assigned value from an inherited one.
- **Microbial and soft-bodied groups have no fossils in PBDB:**
  choanoflagellates and volvocine algae have none, and tardigrades have 6
  occurrences. Their `time_range` comes from curated divergence dates
  where one exists.
- **Genera have no OTT ids.**
- **Synthetic tree:** 67 OTT taxa are "broken" in the synthesis, which means
  they are not monophyletic there. Acari is one of them.

## Open questions before publishing

- **PBDB licence:** CC BY 4.0, as treated here, or CC0, as the API now
  declares? The safe choice is kept.
- **Curated values:** verify each against the full text of its paper.
- **Inputs:** ship the cache as a release artefact, so that others can
  rebuild offline from exactly these inputs?
