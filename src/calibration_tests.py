"""Do later (post-event) soundings under-predict later manifestation? And what drives the
within-window rises in resistance?

1. Timing bins within windows. For P0 (days since 4 Sep 2010) and P2 (days since 13 Jun 2011),
   the metric is demeaned within 250-m cells and averaged in elapsed-time bins, separately for
   sites that did and did not manifest in the window's opening earthquake.
2. Calibration by investigation period. For the 2011 and 2016 outcomes, a logistic model
   y ~ ln(1 + index) + period dummies is fitted on soundings within r of a P0 sounding (shared
   footprint), with 1-km block cluster-robust errors. A positive period coefficient means that, at
   the same index value, later soundings' sites manifested more often, i.e. the later soundings
   under-state the hazard (the bias claimed by Yang et al. 2026). Conditional logit with 500-m
   cell strata is a robustness check. Indices from the inverse-filtered profiles, as in Paper B.

Usage: python -I calibration_tests.py <outputs_dir>
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.spatial import cKDTree
from statsmodels.discrete.conditional_models import ConditionalLogit

sys.path.insert(0, str(Path(__file__).resolve().parent))
from recovery_tests import load  # noqa: E402


def timing_bins(D, out):
    res = {}
    specs = {"P0": ("2010-09-04", "2011-02-22", "label_2010", [0, 30, 60, 90, 151]),
             "P2": ("2011-06-13", "2011-12-23", "label_2011", [0, 7, 30, 90, 173])}
    for per, (t0, t1, lab, edges) in specs.items():
        w = D[(D["date"] > pd.Timestamp(t0)) & (D["date"] < pd.Timestamp(t1))].copy()
        w["days"] = (w["date"] - pd.Timestamp(t0)).dt.days
        w["bin"] = pd.cut(w["days"], edges, right=True)
        w["man"] = np.where(w[lab].between(1, 5), "manifested", np.where(w[lab] == 0, "none", "unknown"))
        res[per] = {"n_by_bin": {str(k): int(v) for k, v in w.groupby("bin", observed=True).size().items()},
                    "programme_by_bin": {str(k): v for k, v in w.groupby("bin", observed=True)["programme"].agg(lambda x: x.value_counts().to_dict()).items()}}
        for m in ["Y_rep", "p10", "med"]:
            x = w.dropna(subset=[f"ln_{m}"])
            x = x[x.groupby("cell")[f"ln_{m}"].transform("size") > 1]
            x["dm"] = x[f"ln_{m}"] - x.groupby("cell")[f"ln_{m}"].transform("mean")
            res[per][m] = {}
            for man, g in x.groupby("man"):
                res[per][m][man] = {str(k): [round(float(v.mean()), 3), int(v.size)] for k, v in g.groupby("bin", observed=True)["dm"]}
    return {str(k): v for k, v in res.items()}


def calibration(out):
    C = pd.read_parquet(out / "cases.parquet")
    Geo = pd.read_parquet(out / "soundings_geo.parquet")[["rec_index", "X2193", "Y2193"]]
    Pp = pd.read_parquet(out / "programme_classes.parquet")[["rec_index", "period", "programme"]]
    res = {}
    for ev in ["2011", "2016"]:
        d = C[(C["variant"] == "inv") & (C["event"] == ev) & C["y"].notna()].merge(Geo, on="rec_index").merge(Pp, on="rec_index")
        d = d[d["period"].notna()]
        d["block"] = np.floor(d["X2193"] / 1000).astype(int).astype(str) + "_" + np.floor(d["Y2193"] / 1000).astype(int).astype(str)
        d["cell500"] = np.floor(d["X2193"] / 500).astype(int).astype(str) + "_" + np.floor(d["Y2193"] / 500).astype(int).astype(str)
        p0 = d[d["period"] == "P0"]
        tree = cKDTree(p0[["X2193", "Y2193"]].values)
        dist = tree.query(d[["X2193", "Y2193"]].values)[0]
        d["dist_p0"] = np.where(d["period"] == "P0", 0.0, dist)
        res[ev] = {}
        for r in [100, 200, 400]:
            s = d[(d["dist_p0"] <= r) & d["period"].isin(["P0", "P1", "P2", "P3"])].copy()
            for p in ["P1", "P2", "P3"]:
                s[f"is_{p}"] = (s["period"] == p).astype(float)
            out_r = {"n": int(len(s)), "n_by_period": s["period"].value_counts().to_dict(),
                     "base_rate_by_period": s.groupby("period", observed=True)["y"].mean().round(3).to_dict()}
            for idx in ["LSN", "LPI", "LPIish"]:
                s["x"] = np.log1p(s[idx])
                X = sm.add_constant(s[["x", "is_P1", "is_P2", "is_P3"]])
                try:
                    fit = sm.Logit(s["y"], X).fit(disp=0, cov_type="cluster", cov_kwds={"groups": pd.factorize(s["block"])[0]})
                    ci = fit.conf_int()
                    out_r[idx] = {k: [float(fit.params[k]), float(ci.loc[k, 0]), float(ci.loc[k, 1])] for k in ["x", "is_P1", "is_P2", "is_P3"]}
                except Exception as e:  # noqa: BLE001
                    out_r[idx] = str(e)
                # calibration in the large: model fitted on P0 only, applied to later periods
                f0 = sm.Logit(s.loc[s["period"] == "P0", "y"], sm.add_constant(s.loc[s["period"] == "P0", ["x"]])).fit(disp=0)
                s["p_hat0"] = f0.predict(sm.add_constant(s[["x"]]))
                out_r[idx + "_calib_from_P0"] = {p: [round(float(g["y"].mean()), 3), round(float(g["p_hat0"].mean()), 3), int(len(g))]
                                                 for p, g in s.groupby("period", observed=True)}
                # conditional logit within 500-m cells
                try:
                    cl = ConditionalLogit(s["y"].values, s[["x", "is_P1", "is_P2", "is_P3"]].values, groups=s["cell500"].values).fit(disp=0)
                    cci = cl.conf_int()
                    out_r[idx + "_condlogit500"] = {k: [float(cl.params[i]), float(cci[i, 0]), float(cci[i, 1])]
                                                    for i, k in enumerate(["x", "is_P1", "is_P2", "is_P3"])}
                except Exception as e:  # noqa: BLE001
                    out_r[idx + "_condlogit500"] = str(e)
                out_r[idx + "_median_index_by_period_manifested"] = s[s["y"] == 1].groupby("period", observed=True)[idx].median().round(2).to_dict()
            res[ev][f"r{r}"] = out_r
    return res


def main(out_dir):
    out = Path(out_dir)
    D = load(out)
    for m in ["Y_rep", "p10", "med"]:
        D[f"ln_{m}"] = np.log(D[m])
    res = {"timing_bins": timing_bins(D, out), "calibration": calibration(out)}
    (out / "calibration_tests.json").write_text(json.dumps(res, indent=1, default=str))
    print(json.dumps(res, indent=1, default=str))


if __name__ == "__main__":
    main(sys.argv[1])
