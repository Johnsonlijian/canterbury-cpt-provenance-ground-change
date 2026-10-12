"""Headline Paper B numbers with uncertainty.

Computes per-event AUC by investigation period with block-bootstrap intervals,
difference-in-differences of the 2010 AUC across February 2011 label strata,
and pooled versus within-event AUC. The output contains only the ``paperB`` key.

Usage: python -I final_numbers.py <outputs_dir>
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

FEB11 = pd.Timestamp("2011-02-22")
JUN11 = pd.Timestamp("2011-06-13")
DEC11 = pd.Timestamp("2011-12-23")
RNG = np.random.default_rng(42)


def auc(y, s):
    y = np.asarray(y)
    return float(roc_auc_score(y, s)) if len(np.unique(y)) > 1 else np.nan


def period(d):
    if pd.isna(d):
        return "nodate"
    if d < FEB11:
        return "P0"
    if d < JUN11:
        return "P1"
    if d < DEC11:
        return "P2"
    if d < pd.Timestamp("2014-01-01"):
        return "P3"
    return "P4"


def block_boot(df, fn, n_boot=2000):
    blocks = df["block"].values
    ub = np.unique(blocks)
    rows = {b: np.where(blocks == b)[0] for b in ub}
    vals = []
    for _ in range(n_boot):
        pick = RNG.choice(ub, size=len(ub), replace=True)
        idx = np.concatenate([rows[b] for b in pick])
        v = fn(df.iloc[idx])
        if np.isfinite(v):
            vals.append(v)
    v = np.array(vals)
    return [float(np.quantile(v, .025)), float(np.quantile(v, .975))]


def main(out_dir):
    out = Path(out_dir)
    c = pd.read_parquet(out / "cases.parquet")
    c = c[(c["variant"] == "inv") & c["y"].notna()].copy()
    c["p0"] = (c["date"] < FEB11).astype(int)
    c["period"] = c["date"].map(period)
    c["block"] = (np.floor(c["E"] / 1000).astype(int) * 100000 + np.floor(c["N"] / 1000).astype(int))
    res = {"paperB": {}}

    # Paper B: per-event AUC by period with CI
    tab = []
    for ev in ["2010", "2011", "2016"]:
        d = c[c["event"] == ev]
        if ev == "2016":
            d = d[~d["cpt_after_event"].astype(bool)]
        for per, g in d.groupby("period"):
            if len(g) < 30 or g["y"].nunique() < 2:
                continue
            row = {"event": ev, "period": per, "n": int(len(g)), "n_pos": int(g["y"].sum()), "base": float(g["y"].mean())}
            for ix in ["LSN", "LPI", "LPIish"]:
                row[f"AUC_{ix}"] = auc(g["y"], g[ix])
                row[f"AUC_{ix}_ci95"] = block_boot(g, lambda s, ix=ix: auc(s["y"], s[ix]), n_boot=1000)
            tab.append(row)
    res["paperB"]["auc_by_event_period"] = tab

    # DiD on 2010 AUC: (P0 - P3) within Feb-2011-manifested minus within not-manifested
    w = c.pivot_table(index="rec_index", columns="event", values="label")
    d = c[(c["event"] == "2010") & c["period"].isin(["P0", "P3"])].merge(
        w["2011"].rename("lab11"), left_on="rec_index", right_index=True)
    d["m11"] = d["lab11"].between(1, 5).astype(int)
    did = {}
    for ix in ["LSN", "LPI", "LPIish"]:
        def did_fn(s, ix=ix):
            a = {}
            for per in ["P0", "P3"]:
                for m in [0, 1]:
                    g = s[(s["period"] == per) & (s["m11"] == m)]
                    a[(per, m)] = auc(g["y"], g[ix])
            return (a[("P0", 1)] - a[("P3", 1)]) - (a[("P0", 0)] - a[("P3", 0)])
        point = {}
        for per in ["P0", "P3"]:
            for m in [0, 1]:
                g = d[(d["period"] == per) & (d["m11"] == m)]
                point[f"AUC_{per}_m11={m}"] = auc(g["y"], g[ix])
                point[f"n_{per}_m11={m}"] = int(len(g))
                point[f"AUC_{per}_m11={m}_ci95"] = block_boot(g, lambda s, ix=ix: auc(s["y"], s[ix]), n_boot=1000)
        point["DiD"] = did_fn(d)
        point["DiD_ci95"] = block_boot(d, did_fn, n_boot=1000)
        did[ix] = point
    res["paperB"]["did_2010_by_feb2011_label"] = did

    # pooled vs within-event AUC (all soundings, inverse)
    pw = {}
    for ix in ["LSN", "LPI", "LPIish"]:
        pooled = auc(c["y"], c[ix])
        within = [auc(g["y"], g[ix]) for _, g in c.groupby("event")]
        ns = [len(g) for _, g in c.groupby("event")]
        pw[ix] = {"pooled": pooled, "within_weighted": float(np.average(within, weights=ns)), "within": within}
    res["paperB"]["pooled_vs_within"] = pw
    (out / "final_numbers.json").write_text(json.dumps(res, indent=2, default=float))
    print(json.dumps(res["paperB"]["did_2010_by_feb2011_label"], indent=1, default=float))


if __name__ == "__main__":
    main(sys.argv[1])
