# Canterbury CPT provenance and depth-limited ground change

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23308668.svg)](https://doi.org/10.5281/zenodo.23308668) [![GitHub release](https://img.shields.io/github/v/release/Johnsonlijian/canterbury-cpt-provenance-ground-change)](https://github.com/Johnsonlijian/canterbury-cpt-provenance-ground-change/releases/tag/v1.0.0)

Public reproducibility package for the study:

> Test dates expose investigation provenance and depth-limited ground change in the Canterbury CPT liquefaction case histories

This release contains the analysis scripts, non-identifying aggregate outputs,
figure exports and data-source links needed to document the Paper B evidence
chain. It deliberately excludes active manuscript files, cover letters,
reviewer material, raw third-party data, per-sounding source tables and NZGD
exports.

## Main result encoded by this release

The Canterbury CPT archive contains 5,667 decoded test dates among 5,668
soundings. The same unrefitted liquefaction index discriminates 2010 Darfield
manifestation more strongly in soundings pushed before the next earthquake
(AUC 0.75–0.78) than in 2012–2013 soundings (AUC 0.59–0.61). At matched
locations, saturated sand-like layers at 1–5 m show about a 50% increase in
cone resistance in 2012–2014, whereas deeper sand, the unsaturated crust and
friction ratio show no detectable change. After property-scale location
control, the data do not support a monotonic elapsed-time recovery law.

## Contents

- `src/` — selected Paper B analysis and figure-generation scripts.
- `outputs/` — aggregate JSON summaries only; no per-sounding data.
- `figures/` — aggregate B2–B7 figure masters/previews and the graphical abstract. Point-level B1/B8 exports remain local because they reveal exact sounding locations or site identifiers.
- `DATASETS_AND_LINKS.csv` — external data sources and access boundaries.
- `REPRODUCIBLE_RUNBOOK.md` — environment, data placement and run order.
- `CITATION.cff` — citation metadata for this release.

## Data boundary

The DesignSafe PRJ-2937 dataset and the New Zealand Geotechnical Database are
external sources. Raw files and per-sounding derived tables are not included.
Users must obtain them directly and comply with their current licences and
terms of use. Aggregate outputs in this repository are provided for method
checking and do not authorize redistribution of source data.

The public figure set intentionally omits the point-level timeline map and
repeat-site panel. Their local versions contain exact sounding locations or
site identifiers; the aggregate results and graphical abstract preserve the
release-level evidence without exposing those fields.

## Release status

Version `v1.0.0` is a public code-and-aggregate-output release. It is separate
from the local submission package and contains no active manuscript PDF,
cover letter or reviewer correspondence.

## License

Original scripts and documentation are released under the MIT License. External
datasets and their derived records remain governed by their source terms.

