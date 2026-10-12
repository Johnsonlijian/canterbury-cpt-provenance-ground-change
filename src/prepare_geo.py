"""Prepare minimal geospatial aggregate for Paper B figures.

Usage: python -I prepare_geo.py <soundings.parquet> <soundings_geo.parquet>

The input must contain WGS84 ``lon`` and ``lat`` columns. The output is a
non-identifying table with ``rec_index``, ``X2193`` and ``Y2193`` only.
"""
from pathlib import Path
import sys
import pandas as pd
from pyproj import Transformer


def main(soundings_path, output_path):
    s = pd.read_parquet(soundings_path)
    required = {"rec_index", "lon", "lat"}
    missing = required - set(s.columns)
    if missing:
        raise ValueError(f"missing required columns: {sorted(missing)}")
    tr = Transformer.from_crs(4326, 2193, always_xy=True)
    x, y = tr.transform(s["lon"].to_numpy(), s["lat"].to_numpy())
    g = pd.DataFrame({"rec_index": s["rec_index"], "X2193": x, "Y2193": y})
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    g.to_parquet(output_path, index=False)
    print(f"wrote {len(g):,} rows to {output_path}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("Usage: python -I prepare_geo.py <soundings.parquet> <soundings_geo.parquet>")
    main(sys.argv[1], sys.argv[2])
