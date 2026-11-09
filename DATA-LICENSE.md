# Data licence

The files in `data/` are a combination of several sources. **The combined
dataset is distributed under the Creative Commons Attribution 4.0
International licence (CC BY 4.0)**, <https://creativecommons.org/licenses/by/4.0/>.
If you redistribute it, or a derivative of it, keep the attribution text
below.

The code that builds it is MIT-licensed (see `LICENSE`).

## Parts and their sources

| Part | Fields | Source | Source licence |
|---|---|---|---|
| Topology, names, ranks, OTT ids | `id`, `ott_id`, `name`, `unique_name`, `rank`, `parent`, `children`, `ott_flags`, `extinct`, `xref.ncbi`, `synth.*`, `synthetic_subtree.nwk` | Open Tree of Life: synthetic tree `opentree16.1`, taxonomy OTT 3.7 (draft 3) | CC0 1.0 (public domain dedication), "where not limited by pre-existing terms of use". OTT aggregates other taxonomies (NCBI, GBIF, WoRMS, IRMNG, ...); only names, ranks, ids and parent links are used here. |
| Fossil ranges, environment, ecology, genera, specimen body lengths, time slices | `pbdb`, `fossil`, `habitat` (when `habitat_src.src` is `PBDB`), `ecology`, all of `genera.jsonl`, `slices.jsonl` | The Paleobiology Database (PBDB), data service 1.2 | Treated here as **CC BY 4.0** (attribution required). Note: on 2026-11-08 the PBDB API's `datainfo` header declared `Creative Commons CC0` for the taxa, occurrence and measurement endpoints; this package nevertheless applies the stricter CC BY 4.0 terms and gives attribution. |
| English common names | `common_name` | NCBI Taxonomy (E-utilities) | US government work, no restrictions placed by NCBI. |
| Curated tables | `divergence`, `size_mm` (when `basis` is `curated`), `habitat` (when `habitat_src.src` is `curated`), `references.json` | Compiled for this package from the published literature cited in each record | CC BY 4.0. The facts are restated with citations; no text, figures or tables are copied from the papers. |

Fields derived by the build (`time_range`, `size_mm` with basis
`derived_from_descendants` or `pbdb_measurements`, `fossil` with basis
`derived_from_descendants`) are computations over the PBDB data above and
carry the same terms.

## Attribution text

> Tree-of-life data package. Topology and taxonomy from the Open Tree of Life
> (synthetic tree opentree16.1, OTT 3.7; CC0). Fossil ranges, environments,
> genera and specimen measurements from the Paleobiology Database
> (paleobiodb.org; CC BY 4.0). Common names from NCBI Taxonomy. Divergence
> dates, sizes and habitats from the published literature cited in the data.
> Combined data licensed CC BY 4.0.

Please also cite the sources directly where you can:

- Open Tree of Life: Hinchliff CE, et al. (2015) Synthesis of phylogeny and
  taxonomy into a comprehensive tree of life. *PNAS* 112: 12764-12769.
  doi:10.1073/pnas.1423041112. Taxonomy: Rees JA, Cranston K (2017)
  Automated assembly of a reference taxonomy for phylogenetic data synthesis.
  *Biodiversity Data Journal* 5: e12581. doi:10.3897/BDJ.5.e12581.
- Paleobiology Database: "Data were downloaded from the Paleobiology
  Database (paleobiodb.org) on <date in data/manifest.json>." PBDB asks users
  to acknowledge the database and, where possible, the contributors of the
  records used.
- The papers listed in `data/references.json` for curated values.

## Not included, and why

- **TimeTree** divergence times: TimeTree's terms of use permit personal
  research and teaching use and state that redistribution of TimeTree data and
  its transformations is not permitted. No TimeTree values are used, even as
  defaults; divergence dates come only from the cited papers.
- **Heim et al. 2015 body-size data** (Stanford Digital Repository,
  purl.stanford.edu/rf761bx8302): the catalogue record (checked 2026-11-08 via
  SearchWorks) carries only the repository's generic use statement and
  "may be subject to additional license and use restrictions applied by the
  depositor"; no open licence is stated. Redistribution rights could not be
  verified, so none of it is shipped.
