# Reproducible runbook

This public repository documents the analysis route without bundling raw
DesignSafe or NZGD inputs. A full rerun therefore requires the user to obtain
the source data directly and place local paths in a private configuration.

## Environment

Python 3.11 with `numpy`, `scipy`, `pandas`, `pyarrow`, `scikit-learn`,
`statsmodels`, `shapely`, `pyproj` and `matplotlib`. A LaTeX installation with
`elsarticle` is needed only for local manuscript compilation; manuscript source
files are intentionally outside this public release.

## Analysis route

1. Obtain PRJ-2937 and NZGD data under their current terms.
2. Decode CPT test dates and construct the case-history table.
3. Link investigation programmes and summarize provenance.
4. Reproduce the reported critical-point resistance metric.
5. Run location-controlled strata, elapsed-time, negative-control and repeat-sounding analyses.
6. Regenerate aggregate B2–B7 figures with `python -I src/make_figures.py <release_root>`; this writes to `figures/regenerated/` and uses only aggregate JSON summaries.
7. Use `python -I src/final_numbers.py <private_outputs_dir>` for the Paper B headline-number route. Use `python -I src/prepare_geo.py <private_soundings.parquet> <private_soundings_geo.parquet>` when a local B1 map is needed. Use `python -I src/nzgd_raw_profiles.py <private_outputs_dir> <private_nzgd_zip_dir>...` to parse local NZGD AGS4 archives.
8. Compare regenerated aggregate outputs with `outputs/`.

The scripts in `src/` accept project-specific input/output paths where the
original analysis route supported them. The default public figure route is
aggregate-only. The private `--full` route of `make_figures.py` requires
point-level inputs and is not part of this release. The scripts are provided as
research code, not as a turnkey data downloader.

