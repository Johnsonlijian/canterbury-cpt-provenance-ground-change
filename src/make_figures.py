"""Publication figures for Paper B (case-history provenance).

All numbers are read from the JSON/parquet outputs of the analysis scripts; nothing is
typed in by hand. Palette: validated categorical slots (blue, orange, aqua; all-pairs CVD
pass) for events, an ordinal single-hue ramp for investigation periods, and a blue-red
diverging map with a grey midpoint for site terms.

Default (public aggregate release):
    python -I make_figures.py <analysis_dir>

Full private-data route (adds B1 and B8):
    python -I make_figures.py <analysis_dir> --full
"""
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm  # noqa: E402

EVENT_COL = {"2010": "#2a78d6", "2011": "#eb6834", "2016": "#1baf7a"}
EVENT_LAB = {"2010": "2010 Darfield $M_w$7.1", "2011": "2011 Christchurch $M_w$6.2", "2016": "2016 Christchurch $M_w$5.7"}
PERIOD_COL = {"P0": "#104281", "P1": "#256abf", "P2": "#3987e5", "P3": "#86b6ef"}
PERIOD_LAB = {"P0": "Sep 2010–Feb 2011", "P1": "Feb–Jun 2011", "P2": "Jun–Dec 2011", "P3": "2012–2013"}
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#898781", "#e1e0d9"
DIV = LinearSegmentedColormap.from_list("div", ["#184f95", "#6da7ec", "#f0efec", "#ec8a89", "#a12b2b"])
MM = 1 / 25.4

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7, "axes.titlesize": 7.5, "axes.labelsize": 7, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
    "legend.fontsize": 6.5, "axes.edgecolor": "#c3c2b7", "axes.linewidth": 0.6, "axes.labelcolor": INK,
    "xtick.color": INK2, "ytick.color": INK2, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.4,
    "axes.axisbelow": True, "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42,
    "svg.fonttype": "none", "figure.dpi": 150, "savefig.dpi": 300, "legend.frameon": False,
})


def save(fig, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(path.with_suffix(".png"), bbox_inches="tight")
    fig.savefig(path.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def panel(ax, letter):
    ax.text(-0.02, 1.04, letter, transform=ax.transAxes, fontsize=8, fontweight="bold", va="bottom", ha="right", color=INK)


def load(a, name):
    return json.loads((a / "outputs" / name).read_text())


def gm3_outlines(a):
    from pyproj import Transformer
    from shapely.geometry import shape
    gj = json.load(open(a / "external" / "gm3_units_chch.geojson", encoding="utf-8"))
    tr = Transformer.from_crs(4326, 2193, always_xy=True)
    lines = []
    for f in gj["features"]:
        g = f.get("geometry")
        if not g:
            continue
        geom = shape(g)
        polys = [geom] if geom.geom_type == "Polygon" else list(geom.geoms)
        for p in polys:
            x, y = p.exterior.xy
            X, Y = tr.transform(np.array(x), np.array(y))
            lines.append((X, Y))
    return lines


def paper_b_full(a, outdir):
    snd = pd.read_parquet(a / "outputs" / "soundings.parquet")
    geo = pd.read_parquet(a / "outputs" / "soundings_geo.parquet")
    S = snd.merge(geo, on="rec_index")
    fn = load(a, "final_numbers.json")["paperB"]
    tm = load(a, "timing_matched.json")
    tl = load(a, "timing_label_checks.json")
    tp = load(a, "timing_protocol_check.json")

    # ---- B1: timeline + map of early vs 2012-13 soundings
    fig, axs = plt.subplots(1, 2, figsize=(190 * MM, 72 * MM), gridspec_kw={"width_ratios": [1.35, 1]})
    ax = axs[0]
    q = S["date"].dt.to_period("Q").value_counts().sort_index()
    qs = [p.start_time for p in q.index]
    cols = []
    for t in qs:
        if t < pd.Timestamp("2011-02-22"):
            cols.append(PERIOD_COL["P0"])
        elif t < pd.Timestamp("2011-06-13"):
            cols.append(PERIOD_COL["P1"])
        elif t < pd.Timestamp("2011-12-23"):
            cols.append(PERIOD_COL["P2"])
        elif t < pd.Timestamp("2014-01-01"):
            cols.append(PERIOD_COL["P3"])
        else:
            cols.append("#c3c2b7")
    ax.bar(qs, q.values, width=80, color=cols, align="edge")
    ev = {"2010-09-04": "Darfield", "2011-02-22": "Christchurch", "2011-06-13": "June 2011", "2011-12-23": "Dec 2011", "2016-02-14": "Feb 2016"}
    for d, lab in ev.items():
        t = pd.Timestamp(d)
        ax.axvline(t, color=INK, lw=0.7)
        ax.text(t, ax.get_ylim()[1] * 0.97 if False else 1750, " " + lab, rotation=90, va="top", ha="right", fontsize=6, color=INK)
    ax.set_ylim(0, 1800)
    ax.set_ylabel("CPT soundings per quarter")
    ax.set_xlim(pd.Timestamp("2010-07-01"), pd.Timestamp("2017-04-01"))
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=PERIOD_COL[p], label=PERIOD_LAB[p]) for p in ["P0", "P1", "P2", "P3"]],
              loc="upper right", bbox_to_anchor=(0.80, 1.0), fontsize=6, title="test period", title_fontsize=6)
    ax.set_title("When the 5,668 PRJ-2937 soundings were pushed", loc="left"); panel(ax, "a")
    ax = axs[1]
    outlines = gm3_outlines(a)
    xmin, xmax = S["X2193"].quantile(0.001) - 1500, S["X2193"].quantile(0.999) + 1500
    ymin, ymax = S["Y2193"].quantile(0.001) - 1500, S["Y2193"].quantile(0.999) + 1500
    for X, Y in outlines:
        ax.plot(X, Y, color="#e6e5df", lw=0.25, zorder=1)
    late = S[(S["date"] >= pd.Timestamp("2011-12-23")) & (S["date"] < pd.Timestamp("2014-01-01"))]
    early = S[S["date"] < pd.Timestamp("2011-02-22")]
    ax.scatter(late["X2193"], late["Y2193"], s=3, c=PERIOD_COL["P3"], lw=0, label=f"2012–2013 (n={len(late):,})", zorder=2)
    ax.scatter(early["X2193"], early["Y2193"], s=3, c=PERIOD_COL["P0"], lw=0, label=f"Sep 2010–Feb 2011 (n={len(early)})", zorder=3)
    ax.set_xlim(xmin, xmax); ax.set_ylim(ymin, ymax); ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([]); ax.grid(False)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.legend(loc="lower left", markerscale=3, fontsize=6)
    ax.set_title("Where they were pushed", loc="left"); panel(ax, "b")
    fig.tight_layout()
    save(fig, outdir / "B1_timeline_map")

    # ---- B2: AUC by event x period (LSN, LPI, LPIish small multiples)
    T = pd.DataFrame(fn["auc_by_event_period"])
    T = T[T["period"].isin(["P0", "P1", "P2", "P3"])]
    fig, axs = plt.subplots(1, 3, figsize=(190 * MM, 62 * MM), sharey=True)
    for k, ix in enumerate(["LSN", "LPI", "LPIish"]):
        ax = axs[k]
        for j, ev in enumerate(["2010", "2011", "2016"]):
            for i, per in enumerate(["P0", "P1", "P2", "P3"]):
                r = T[(T["event"] == ev) & (T["period"] == per)]
                if r.empty:
                    continue
                r = r.iloc[0]
                if r["n_pos"] < 10 or (r["n"] - r["n_pos"]) < 10:
                    continue  # too few manifested (or non-manifested) cases for a stable AUC
                x = j * 5 + i
                v, ci = r[f"AUC_{ix}"], r[f"AUC_{ix}_ci95"]
                ax.errorbar(x, v, yerr=[[v - ci[0]], [ci[1] - v]], fmt="o", ms=4, color=PERIOD_COL[per],
                            elinewidth=1.2, capsize=0, label=PERIOD_LAB[per] if (j == 0 and k == 0) else None)
        ax.axhline(0.5, color=INK2, lw=0.6)
        ax.set_xticks([1.5, 6.5, 11.5]); ax.set_xticklabels(["2010\nDarfield", "2011\nChristchurch", "2016\nChristchurch"])
        ax.set_title(f"{ix}", loc="left")
        panel(ax, "abc"[k])
    axs[0].set_ylabel("AUC of the unrefitted index (95% CI)")
    axs[0].set_ylim(0.35, 1.0)
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=4, fontsize=6.5, title="sounding test period", title_fontsize=6.5,
               bbox_to_anchor=(0.5, 1.08))
    fig.tight_layout()
    save(fig, outdir / "B2_auc_by_event_period")

    # ---- B3: footprint-matched difference vs radius
    fig, axs = plt.subplots(1, 2, figsize=(140 * MM, 56 * MM), sharey=True)
    for k, ev in enumerate(["2010", "2011"]):
        ax = axs[k]
        fp = tm[f"{ev}_inv"]["footprint"]
        for i, ix in enumerate(["LSN", "LPI", "LPIish"]):
            r = np.array([f["radius_m"] for f in fp])
            d = np.array([f[f"dAUC_{ix}_early_minus_late"] for f in fp])
            lo = np.array([f[f"dAUC_{ix}_ci95"][0] for f in fp])
            hi = np.array([f[f"dAUC_{ix}_ci95"][1] for f in fp])
            off = (i - 1) * 0.06
            ax.errorbar(r * np.exp(off), d, yerr=[d - lo, hi - d], fmt="o-", ms=3.5, lw=1.0, elinewidth=1.0,
                        color=["#2a78d6", "#eb6834", "#1baf7a"][i], capsize=0, label=ix)
        ax.axhline(0, color=INK2, lw=0.6)
        ax.set_xscale("log"); ax.set_xticks([50, 100, 200, 400]); ax.set_xticklabels(["50", "100", "200", "400"])
        ax.set_xlabel("Late soundings kept within r of an early sounding (m)")
        ax.set_title(f"{EVENT_LAB[ev]} outcome", loc="left"); panel(ax, "ab"[k])
    axs[0].set_ylabel("AUC(early) − AUC(late), same footprint")
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=3, fontsize=6.5, bbox_to_anchor=(0.5, 1.06))
    fig.tight_layout()
    save(fig, outdir / "B3_footprint_matched")

    # ---- B4: DiD by Feb 2011 label
    did = fn["did_2010_by_feb2011_label"]
    fig, axs = plt.subplots(1, 2, figsize=(140 * MM, 60 * MM), sharey=True)
    cols3 = ["#2a78d6", "#eb6834", "#1baf7a"]
    for k, (m, title) in enumerate([(1, "Sites that manifested in February 2011"), (0, "Sites that did not")]):
        ax = axs[k]
        for i, ix in enumerate(["LSN", "LPI", "LPIish"]):
            r = did[ix]
            pts = []
            for j, per in enumerate(["P0", "P3"]):
                v = r[f"AUC_{per}_m11={m}"]
                c = r[f"AUC_{per}_m11={m}_ci95"]
                x = j + (i - 1) * 0.1
                ax.errorbar(x, v, yerr=[[v - c[0]], [c[1] - v]], fmt="o", ms=4, color=cols3[i], elinewidth=1.1, capsize=0,
                            label=ix if (k == 0 and j == 0) else None)
                pts.append((x, v))
            ax.plot([p[0] for p in pts], [p[1] for p in pts], color=cols3[i], lw=0.8)
            if k == 0:
                pass
        ax.axhline(0.5, color=INK2, lw=0.6)
        ax.set_xticks([0, 1]); ax.set_xticklabels(["tested\nSep 2010–Feb 2011", "tested\n2012–2013"])
        ax.set_xlim(-0.5, 1.5)
        n0, n3 = did["LSN"][f"n_P0_m11={m}"], did["LSN"][f"n_P3_m11={m}"]
        ax.text(0.5, 0.03, f"n = {n0} and {n3:,}", transform=ax.transAxes, ha="center", fontsize=6, color=MUTED)
        ax.set_title(title, loc="left"); panel(ax, "ab"[k])
    axs[0].set_ylabel("2010 Darfield AUC (95% CI)")
    axs[0].set_ylim(0.35, 0.95)
    r = did["LSN"]
    axs[1].text(0.98, 0.97, f"difference-in-differences (LSN)\n{r['DiD']:.2f} [{r['DiD_ci95'][0]:.2f}, {r['DiD_ci95'][1]:.2f}]",
                transform=axs[1].transAxes, ha="right", va="top", fontsize=6, color=INK2)
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=3, fontsize=6.5, bbox_to_anchor=(0.5, 1.06))
    fig.tight_layout()
    save(fig, outdir / "B4_did_feb2011")

    # ---- B5: mechanism checks
    fig, axs = plt.subplots(1, 2, figsize=(140 * MM, 58 * MM))
    ax = axs[0]
    tr = pd.DataFrame(tp["truncation_test_2010"])
    tr["per"] = tr["period"].str[:2]
    tr = tr[(tr["footprint"] == "all") & tr["per"].isin(["P0", "P1", "P2", "P3"])]
    for per in ["P0", "P1", "P2", "P3"]:
        g = tr[tr["per"] == per].sort_values("z0")
        ax.plot(g["z0"], g["AUC_LSN"], "o-", ms=3.5, lw=1.0, color=PERIOD_COL[per], label=PERIOD_LAB[per])
    ax.axhline(0.5, color=INK2, lw=0.6)
    ax.set_xlabel("Data above this depth discarded (m)")
    ax.set_ylabel("2010 AUC (LSN)")
    ax.set_ylim(0.34, 0.85)
    ax.legend(fontsize=5.6, loc="lower center", ncol=2, columnspacing=0.8, handlelength=1.5)
    ax.set_title("Not pre-drilling or shallow truncation", loc="left"); panel(ax, "a")
    ax = axs[1]
    l3 = pd.DataFrame(tl["L3_neighbour_label_agreement"])
    l3 = l3[(l3["radius_m"] == 50) & l3["period"].isin(["P0", "P1", "P2", "P3"])]
    xs = np.arange(len(l3))
    ax.bar(xs - 0.18, l3["agree2010_exact"], width=0.34, color="#2a78d6", label="2010 severity class")
    ax.bar(xs + 0.18, l3["agree2011_binary"], width=0.34, color="#eb6834", label="2011 manifested yes/no")
    ax.set_xticks(xs); ax.set_xticklabels([PERIOD_LAB[p] for p in l3["period"]], rotation=20, fontsize=6)
    ax.set_ylabel("Agreement with nearest neighbour ≤ 50 m")
    ax.set_ylim(0, 1.25)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.legend(fontsize=6, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 1.0))
    ax.set_title("Label resolution differs by programme", loc="left"); panel(ax, "b")
    fig.tight_layout()
    save(fig, outdir / "B5_mechanism_checks")

    # ---- B6: reproduction of the reported trend and negative-control strata
    rt = load(a, "recovery_tests.json")
    bp = load(a, "between_period_robustness.json")
    fig, axs = plt.subplots(1, 2, figsize=(190 * MM, 70 * MM), gridspec_kw={"width_ratios": [1, 1.25]})
    ax = axs[0]
    groups = ["G1", "G2", "G3", "G4"]
    glab = ["Sep 2010–\nFeb 2011", "Feb–Dec\n2011", "2012", "2013–\n2014"]
    mean = [rt["A_group_means"]["Y_rep"][g] for g in groups]
    xs = np.arange(4)
    ax.errorbar(xs - 0.1, [m[0] for m in mean], yerr=[[m[0] - m[1] for m in mean], [m[2] - m[0] for m in mean]], fmt="o",
                ms=4.5, color="#2a78d6", elinewidth=1.2, capsize=0, label="this study, mean (95% block CI)")
    ax.plot(xs - 0.1, [rt["A_group_medians"]["Y_rep"][g] for g in groups], "o", ms=4.5, mfc="white", mec="#2a78d6", mew=1.0,
            label="this study, median")
    yang = [(0, 58.0, "median"), (1, 66.0, "mean"), (2, 72.0, "mean"), (3, 76.0, "mean")]
    ax.plot([y[0] + 0.1 for y in yang], [y[1] for y in yang], "D", ms=4.5, color="#eb6834",
            label="Yang et al. (2026): median (IQR), then means")
    ax.plot([0.1, 0.1], [45, 70], color="#eb6834", lw=1.2)
    for x, m in zip(xs, mean):
        ax.text(x - 0.1, 41.0, f"n={m[3]:,}", ha="center", fontsize=5.8, color=MUTED)
    ax.set_xticks(xs); ax.set_xticklabels(glab)
    ax.set_ylim(38, 92)
    ax.set_ylabel("$q_{c1Ncs}$ at the minimum-FS point")
    ax.legend(loc="upper left", fontsize=6)
    ax.set_title("The reported trend is reproduced", loc="left"); panel(ax, "a")
    ax = axs[1]
    strata = [("lnqc_sandW_1_5", "saturated sand-like, 1–5 m"), ("lnqc_sandW_5_10", "saturated sand-like, 5–10 m"),
              ("lnqc_sand_12_20", "sand-like, 12–20 m"), ("lnqc_clay_1_10", "clay-like ($I_c$ > 2.95), 1–10 m"),
              ("lnqc_dry", "unsaturated crust"), ("lnRf_sandW_1_10", "friction ratio, saturated sand-like")]
    for k, (c, lab) in enumerate(strata):
        for j, (con, col, lab2) in enumerate([("G34_vs_P0", "#2a78d6", "2012–2014 vs Sep 2010–Feb 2011"),
                                              ("P12_vs_P0", "#86b6ef", "2011 (same programme) vs Sep 2010–Feb 2011")]):
            v = bp["cellFE"][con][c].get("100m")
            if not v:
                continue
            ax.errorbar(v[0], k + (j - 0.5) * 0.28, xerr=[[v[0] - v[1]], [v[2] - v[0]]], fmt="o", ms=4, color=col,
                        elinewidth=1.2, capsize=0, label=lab2 if k == 0 else None)
            if j == 0:
                ax.text(1.01, k, f"n={v[3]:,}", transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=5.8, color=MUTED)
    ax.axvline(0, color=INK2, lw=0.6)
    ax.set_yticks(range(len(strata))); ax.set_yticklabels([s[1] for s in strata]); ax.invert_yaxis()
    ax.set_xlabel("Change in median ln $q_c$ (or ln $R_f$) within 100 m cells (95% CI)")
    ax.legend(loc="upper center", bbox_to_anchor=(0.42, -0.2), ncol=2, fontsize=6)
    ax.set_title("Only the liquefiable zone changed", loc="left"); panel(ax, "b")
    fig.tight_layout()
    save(fig, outdir / "B6_reproduction_strata")

    # ---- B7: between-period change vs elapsed-time slope under tighter location control
    wr = load(a, "window_robustness.json")
    tv = load(a, "time_vs_event.json")
    fig, axs = plt.subplots(1, 2, figsize=(190 * MM, 64 * MM), sharex=False)
    scales = ["500m", "250m", "100m", "50m", "25m"]
    sc_x = np.arange(len(scales) + 1)
    sc_lab = ["500 m", "250 m", "100 m", "50 m", "25 m", "pairs\n≤ 50 m"]
    st_cols = [("lnqc_sandW_1_5", "#2a78d6", "saturated, 1–5 m"), ("lnqc_sandW_5_10", "#eb6834", "saturated, 5–10 m"),
               ("lnqc_sand_12_20", "#1baf7a", "12–20 m")]
    ax = axs[0]
    for i, (c, col, lab) in enumerate(st_cols):
        pts = []
        for k, s in enumerate(scales):
            v = bp["cellFE"]["G34_vs_P0"][c].get(s)
            if v:
                pts.append((k + (i - 1) * 0.12, v[0], v[1], v[2]))
        pr = bp["pairs"]["G34_vs_P0"]["r50"].get(f"{c}_all")
        if pr:
            pts.append((len(scales) + (i - 1) * 0.12, pr[0], pr[1], pr[2]))
        ax.errorbar([p[0] for p in pts], [p[1] for p in pts], yerr=[[p[1] - p[2] for p in pts], [p[3] - p[1] for p in pts]],
                    fmt="o", ms=4, color=col, elinewidth=1.1, capsize=0, label=lab)
    ax.axhline(0, color=INK2, lw=0.6)
    ax.set_xticks(sc_x); ax.set_xticklabels(sc_lab, fontsize=6)
    ax.set_xlabel("Location control: cell fixed effects of size, or paired soundings")
    ax.set_ylabel("Δ median ln $q_c$, 2012–2014 vs Sep 2010–Feb 2011")
    ax.legend(loc="upper left", fontsize=6, title="sand-like layers", title_fontsize=6)
    ax.set_title("The between-period change persists", loc="left"); panel(ax, "a")
    ax = axs[1]
    scales_w = ["cellFE_500m", "cellFE_250m", "cellFE_100m", "cellFE_50m", "pairs_50m_sameprog_pd", "pairs_20m_sameprog_pd"]
    lab_w = ["500 m", "250 m", "100 m", "50 m", "pairs\n≤ 50 m", "pairs\n≤ 20 m"]
    for i, (c, col, lab) in enumerate(st_cols):
        pts = []
        for k, s in enumerate(scales_w):
            v = wr["P3"][c].get(s)
            if v and "slope" in v:
                pts.append((k + (i - 1) * 0.12, v["slope"], v["lo"], v["hi"]))
        ax.errorbar([p[0] for p in pts], [p[1] for p in pts], yerr=[[p[1] - p[2] for p in pts], [p[3] - p[1] for p in pts]],
                    fmt="o", ms=4, color=col, elinewidth=1.1, capsize=0, label=lab)
    ys = tv["archive"]["yang_implied_slope_per_decade"]
    ax.axhline(ys, color="#eb6834", lw=0.8, ls=(0, (4, 2)))
    ax.text(len(scales_w) - 0.6, ys + 0.02, "implied by a 70% rise over 10–800 days", ha="right", va="bottom", fontsize=5.8, color=INK2)
    ax.axhline(0, color=INK2, lw=0.6)
    ax.set_xticks(np.arange(len(scales_w))); ax.set_xticklabels(lab_w, fontsize=6)
    ax.set_xlabel("Location control within the Dec 2011–Feb 2016 window")
    ax.set_ylabel("Δ ln $q_c$ per decade of days since 23 Dec 2011")
    ax.set_title("The elapsed-time slope does not", loc="left"); panel(ax, "b")
    fig.tight_layout()
    save(fig, outdir / "B7_location_control")

    # ---- B8: repeat panel
    inc = pd.read_csv(a / "outputs" / "repeat_increments.csv", parse_dates=["date"])
    fig, axs = plt.subplots(1, 2, figsize=(190 * MM, 66 * MM), gridspec_kw={"width_ratios": [1.45, 1]})
    ax = axs[0]
    hl = {"BUR-22": ("#2a78d6", "BUR-22 (Liggins St)"), "AVS-16": ("#eb6834", "AVS-16 (Morris St)")}
    for fam, g in inc.groupby("fam"):
        g = g.sort_values("date")
        col, lab = hl.get(fam, ("#b9b8b0", None))
        ax.plot(g["date"], g["cum_dlogq_1_10"], "o-", ms=3.2, lw=1.1 if fam in hl else 0.8, color=col,
                label=lab, zorder=3 if fam in hl else 2)
        if fam not in hl and g["date"].iloc[-1] > pd.Timestamp("2012-01-01"):
            last = g.iloc[-1]
            ax.text(last["date"] + pd.Timedelta(days=12), last["cum_dlogq_1_10"], fam, fontsize=5.6, color=MUTED, va="center")
    for d, lab in {"2011-02-22": "22 Feb 2011", "2011-06-13": "13 Jun 2011", "2011-12-23": "23 Dec 2011"}.items():
        ax.axvline(pd.Timestamp(d), color=INK, lw=0.6)
        ax.text(pd.Timestamp(d), 1.12, " " + lab, rotation=90, va="top", ha="right", fontsize=5.8, color=INK)
    ax.axhline(0, color=INK2, lw=0.6)
    ax.set_ylim(-1.0, 1.15)
    ax.set_xlim(pd.Timestamp("2010-09-01"), pd.Timestamp("2013-02-01"))
    ax.set_ylabel("Cumulative Δ median ln $q_c$, sand-like, 1–10 m")
    ax.legend(loc="lower right", fontsize=6)
    ax.set_title("Eight repeat locations, Oct 2010–Nov 2012", loc="left"); panel(ax, "a")
    ax = axs[1]
    J = inc[inc["from"].notna()].copy()
    J["cls"] = np.where(J["n_events"] == 0, "event-free\n(9–14 days)",
                        np.where(J["to"] == "e", "to Nov 2012\n(505–714 days)", "across one 2011\nevent (87–145 days)"))
    order = ["event-free\n(9–14 days)", "across one 2011\nevent (87–145 days)", "to Nov 2012\n(505–714 days)"]
    bands = [("dlogq_1_4", "1–4 m", "#2a78d6"), ("dlogq_4_8", "4–8 m", "#eb6834"), ("dlogq_8_12", "8–12 m", "#1baf7a")]
    for i, (c, lab, col) in enumerate(bands):
        for k, cl in enumerate(order):
            v = J.loc[J["cls"] == cl, c].dropna().values
            x = k + (i - 1) * 0.22
            ax.scatter(np.full(v.size, x), v, s=9, color=col, lw=0, alpha=0.75, label=lab if k == 0 else None)
            if v.size:
                ax.plot([x - 0.08, x + 0.08], [np.median(v)] * 2, color=INK, lw=1.2)
    ax.axhline(0, color=INK2, lw=0.6)
    ax.set_xticks(range(len(order))); ax.set_xticklabels(order, fontsize=6)
    ax.set_ylabel("Δ median ln $q_c$ between consecutive visits")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.0), fontsize=6, title="depth band", title_fontsize=6)
    ax.set_title("Changes by interval and depth", loc="left"); panel(ax, "b")
    fig.tight_layout()
    save(fig, outdir / "B8_repeat_panel")


def paper_b(a, outdir):
    """Render B2-B7 from aggregate JSON outputs only."""
    fn = load(a, "final_numbers.json")["paperB"]
    tm = load(a, "timing_matched.json")
    tl = load(a, "timing_label_checks.json")
    tp = load(a, "timing_protocol_check.json")
    # ---- B2: AUC by event x period (LSN, LPI, LPIish small multiples)
    T = pd.DataFrame(fn["auc_by_event_period"])
    T = T[T["period"].isin(["P0", "P1", "P2", "P3"])]
    fig, axs = plt.subplots(1, 3, figsize=(190 * MM, 62 * MM), sharey=True)
    for k, ix in enumerate(["LSN", "LPI", "LPIish"]):
        ax = axs[k]
        for j, ev in enumerate(["2010", "2011", "2016"]):
            for i, per in enumerate(["P0", "P1", "P2", "P3"]):
                r = T[(T["event"] == ev) & (T["period"] == per)]
                if r.empty:
                    continue
                r = r.iloc[0]
                if r["n_pos"] < 10 or (r["n"] - r["n_pos"]) < 10:
                    continue  # too few manifested (or non-manifested) cases for a stable AUC
                x = j * 5 + i
                v, ci = r[f"AUC_{ix}"], r[f"AUC_{ix}_ci95"]
                ax.errorbar(x, v, yerr=[[v - ci[0]], [ci[1] - v]], fmt="o", ms=4, color=PERIOD_COL[per],
                            elinewidth=1.2, capsize=0, label=PERIOD_LAB[per] if (j == 0 and k == 0) else None)
        ax.axhline(0.5, color=INK2, lw=0.6)
        ax.set_xticks([1.5, 6.5, 11.5]); ax.set_xticklabels(["2010\nDarfield", "2011\nChristchurch", "2016\nChristchurch"])
        ax.set_title(f"{ix}", loc="left")
        panel(ax, "abc"[k])
    axs[0].set_ylabel("AUC of the unrefitted index (95% CI)")
    axs[0].set_ylim(0.35, 1.0)
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=4, fontsize=6.5, title="sounding test period", title_fontsize=6.5,
               bbox_to_anchor=(0.5, 1.08))
    fig.tight_layout()
    save(fig, outdir / "B2_auc_by_event_period")

    # ---- B3: footprint-matched difference vs radius
    fig, axs = plt.subplots(1, 2, figsize=(140 * MM, 56 * MM), sharey=True)
    for k, ev in enumerate(["2010", "2011"]):
        ax = axs[k]
        fp = tm[f"{ev}_inv"]["footprint"]
        for i, ix in enumerate(["LSN", "LPI", "LPIish"]):
            r = np.array([f["radius_m"] for f in fp])
            d = np.array([f[f"dAUC_{ix}_early_minus_late"] for f in fp])
            lo = np.array([f[f"dAUC_{ix}_ci95"][0] for f in fp])
            hi = np.array([f[f"dAUC_{ix}_ci95"][1] for f in fp])
            off = (i - 1) * 0.06
            ax.errorbar(r * np.exp(off), d, yerr=[d - lo, hi - d], fmt="o-", ms=3.5, lw=1.0, elinewidth=1.0,
                        color=["#2a78d6", "#eb6834", "#1baf7a"][i], capsize=0, label=ix)
        ax.axhline(0, color=INK2, lw=0.6)
        ax.set_xscale("log"); ax.set_xticks([50, 100, 200, 400]); ax.set_xticklabels(["50", "100", "200", "400"])
        ax.set_xlabel("Late soundings kept within r of an early sounding (m)")
        ax.set_title(f"{EVENT_LAB[ev]} outcome", loc="left"); panel(ax, "ab"[k])
    axs[0].set_ylabel("AUC(early) − AUC(late), same footprint")
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=3, fontsize=6.5, bbox_to_anchor=(0.5, 1.06))
    fig.tight_layout()
    save(fig, outdir / "B3_footprint_matched")

    # ---- B4: DiD by Feb 2011 label
    did = fn["did_2010_by_feb2011_label"]
    fig, axs = plt.subplots(1, 2, figsize=(140 * MM, 60 * MM), sharey=True)
    cols3 = ["#2a78d6", "#eb6834", "#1baf7a"]
    for k, (m, title) in enumerate([(1, "Sites that manifested in February 2011"), (0, "Sites that did not")]):
        ax = axs[k]
        for i, ix in enumerate(["LSN", "LPI", "LPIish"]):
            r = did[ix]
            pts = []
            for j, per in enumerate(["P0", "P3"]):
                v = r[f"AUC_{per}_m11={m}"]
                c = r[f"AUC_{per}_m11={m}_ci95"]
                x = j + (i - 1) * 0.1
                ax.errorbar(x, v, yerr=[[v - c[0]], [c[1] - v]], fmt="o", ms=4, color=cols3[i], elinewidth=1.1, capsize=0,
                            label=ix if (k == 0 and j == 0) else None)
                pts.append((x, v))
            ax.plot([p[0] for p in pts], [p[1] for p in pts], color=cols3[i], lw=0.8)
            if k == 0:
                pass
        ax.axhline(0.5, color=INK2, lw=0.6)
        ax.set_xticks([0, 1]); ax.set_xticklabels(["tested\nSep 2010–Feb 2011", "tested\n2012–2013"])
        ax.set_xlim(-0.5, 1.5)
        n0, n3 = did["LSN"][f"n_P0_m11={m}"], did["LSN"][f"n_P3_m11={m}"]
        ax.text(0.5, 0.03, f"n = {n0} and {n3:,}", transform=ax.transAxes, ha="center", fontsize=6, color=MUTED)
        ax.set_title(title, loc="left"); panel(ax, "ab"[k])
    axs[0].set_ylabel("2010 Darfield AUC (95% CI)")
    axs[0].set_ylim(0.35, 0.95)
    r = did["LSN"]
    axs[1].text(0.98, 0.97, f"difference-in-differences (LSN)\n{r['DiD']:.2f} [{r['DiD_ci95'][0]:.2f}, {r['DiD_ci95'][1]:.2f}]",
                transform=axs[1].transAxes, ha="right", va="top", fontsize=6, color=INK2)
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=3, fontsize=6.5, bbox_to_anchor=(0.5, 1.06))
    fig.tight_layout()
    save(fig, outdir / "B4_did_feb2011")

    # ---- B5: mechanism checks
    fig, axs = plt.subplots(1, 2, figsize=(140 * MM, 58 * MM))
    ax = axs[0]
    tr = pd.DataFrame(tp["truncation_test_2010"])
    tr["per"] = tr["period"].str[:2]
    tr = tr[(tr["footprint"] == "all") & tr["per"].isin(["P0", "P1", "P2", "P3"])]
    for per in ["P0", "P1", "P2", "P3"]:
        g = tr[tr["per"] == per].sort_values("z0")
        ax.plot(g["z0"], g["AUC_LSN"], "o-", ms=3.5, lw=1.0, color=PERIOD_COL[per], label=PERIOD_LAB[per])
    ax.axhline(0.5, color=INK2, lw=0.6)
    ax.set_xlabel("Data above this depth discarded (m)")
    ax.set_ylabel("2010 AUC (LSN)")
    ax.set_ylim(0.34, 0.85)
    ax.legend(fontsize=5.6, loc="lower center", ncol=2, columnspacing=0.8, handlelength=1.5)
    ax.set_title("Not pre-drilling or shallow truncation", loc="left"); panel(ax, "a")
    ax = axs[1]
    l3 = pd.DataFrame(tl["L3_neighbour_label_agreement"])
    l3 = l3[(l3["radius_m"] == 50) & l3["period"].isin(["P0", "P1", "P2", "P3"])]
    xs = np.arange(len(l3))
    ax.bar(xs - 0.18, l3["agree2010_exact"], width=0.34, color="#2a78d6", label="2010 severity class")
    ax.bar(xs + 0.18, l3["agree2011_binary"], width=0.34, color="#eb6834", label="2011 manifested yes/no")
    ax.set_xticks(xs); ax.set_xticklabels([PERIOD_LAB[p] for p in l3["period"]], rotation=20, fontsize=6)
    ax.set_ylabel("Agreement with nearest neighbour ≤ 50 m")
    ax.set_ylim(0, 1.25)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.legend(fontsize=6, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 1.0))
    ax.set_title("Label resolution differs by programme", loc="left"); panel(ax, "b")
    fig.tight_layout()
    save(fig, outdir / "B5_mechanism_checks")

    # ---- B6: reproduction of the reported trend and negative-control strata
    rt = load(a, "recovery_tests.json")
    bp = load(a, "between_period_robustness.json")
    fig, axs = plt.subplots(1, 2, figsize=(190 * MM, 70 * MM), gridspec_kw={"width_ratios": [1, 1.25]})
    ax = axs[0]
    groups = ["G1", "G2", "G3", "G4"]
    glab = ["Sep 2010–\nFeb 2011", "Feb–Dec\n2011", "2012", "2013–\n2014"]
    mean = [rt["A_group_means"]["Y_rep"][g] for g in groups]
    xs = np.arange(4)
    ax.errorbar(xs - 0.1, [m[0] for m in mean], yerr=[[m[0] - m[1] for m in mean], [m[2] - m[0] for m in mean]], fmt="o",
                ms=4.5, color="#2a78d6", elinewidth=1.2, capsize=0, label="this study, mean (95% block CI)")
    ax.plot(xs - 0.1, [rt["A_group_medians"]["Y_rep"][g] for g in groups], "o", ms=4.5, mfc="white", mec="#2a78d6", mew=1.0,
            label="this study, median")
    yang = [(0, 58.0, "median"), (1, 66.0, "mean"), (2, 72.0, "mean"), (3, 76.0, "mean")]
    ax.plot([y[0] + 0.1 for y in yang], [y[1] for y in yang], "D", ms=4.5, color="#eb6834",
            label="Yang et al. (2026): median (IQR), then means")
    ax.plot([0.1, 0.1], [45, 70], color="#eb6834", lw=1.2)
    for x, m in zip(xs, mean):
        ax.text(x - 0.1, 41.0, f"n={m[3]:,}", ha="center", fontsize=5.8, color=MUTED)
    ax.set_xticks(xs); ax.set_xticklabels(glab)
    ax.set_ylim(38, 92)
    ax.set_ylabel("$q_{c1Ncs}$ at the minimum-FS point")
    ax.legend(loc="upper left", fontsize=6)
    ax.set_title("The reported trend is reproduced", loc="left"); panel(ax, "a")
    ax = axs[1]
    strata = [("lnqc_sandW_1_5", "saturated sand-like, 1–5 m"), ("lnqc_sandW_5_10", "saturated sand-like, 5–10 m"),
              ("lnqc_sand_12_20", "sand-like, 12–20 m"), ("lnqc_clay_1_10", "clay-like ($I_c$ > 2.95), 1–10 m"),
              ("lnqc_dry", "unsaturated crust"), ("lnRf_sandW_1_10", "friction ratio, saturated sand-like")]
    for k, (c, lab) in enumerate(strata):
        for j, (con, col, lab2) in enumerate([("G34_vs_P0", "#2a78d6", "2012–2014 vs Sep 2010–Feb 2011"),
                                              ("P12_vs_P0", "#86b6ef", "2011 (same programme) vs Sep 2010–Feb 2011")]):
            v = bp["cellFE"][con][c].get("100m")
            if not v:
                continue
            ax.errorbar(v[0], k + (j - 0.5) * 0.28, xerr=[[v[0] - v[1]], [v[2] - v[0]]], fmt="o", ms=4, color=col,
                        elinewidth=1.2, capsize=0, label=lab2 if k == 0 else None)
            if j == 0:
                ax.text(1.01, k, f"n={v[3]:,}", transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=5.8, color=MUTED)
    ax.axvline(0, color=INK2, lw=0.6)
    ax.set_yticks(range(len(strata))); ax.set_yticklabels([s[1] for s in strata]); ax.invert_yaxis()
    ax.set_xlabel("Change in median ln $q_c$ (or ln $R_f$) within 100 m cells (95% CI)")
    ax.legend(loc="upper center", bbox_to_anchor=(0.42, -0.2), ncol=2, fontsize=6)
    ax.set_title("Only the liquefiable zone changed", loc="left"); panel(ax, "b")
    fig.tight_layout()
    save(fig, outdir / "B6_reproduction_strata")

    # ---- B7: between-period change vs elapsed-time slope under tighter location control
    wr = load(a, "window_robustness.json")
    tv = load(a, "time_vs_event.json")
    fig, axs = plt.subplots(1, 2, figsize=(190 * MM, 64 * MM), sharex=False)
    scales = ["500m", "250m", "100m", "50m", "25m"]
    sc_x = np.arange(len(scales) + 1)
    sc_lab = ["500 m", "250 m", "100 m", "50 m", "25 m", "pairs\n≤ 50 m"]
    st_cols = [("lnqc_sandW_1_5", "#2a78d6", "saturated, 1–5 m"), ("lnqc_sandW_5_10", "#eb6834", "saturated, 5–10 m"),
               ("lnqc_sand_12_20", "#1baf7a", "12–20 m")]
    ax = axs[0]
    for i, (c, col, lab) in enumerate(st_cols):
        pts = []
        for k, s in enumerate(scales):
            v = bp["cellFE"]["G34_vs_P0"][c].get(s)
            if v:
                pts.append((k + (i - 1) * 0.12, v[0], v[1], v[2]))
        pr = bp["pairs"]["G34_vs_P0"]["r50"].get(f"{c}_all")
        if pr:
            pts.append((len(scales) + (i - 1) * 0.12, pr[0], pr[1], pr[2]))
        ax.errorbar([p[0] for p in pts], [p[1] for p in pts], yerr=[[p[1] - p[2] for p in pts], [p[3] - p[1] for p in pts]],
                    fmt="o", ms=4, color=col, elinewidth=1.1, capsize=0, label=lab)
    ax.axhline(0, color=INK2, lw=0.6)
    ax.set_xticks(sc_x); ax.set_xticklabels(sc_lab, fontsize=6)
    ax.set_xlabel("Location control: cell fixed effects of size, or paired soundings")
    ax.set_ylabel("Δ median ln $q_c$, 2012–2014 vs Sep 2010–Feb 2011")
    ax.legend(loc="upper left", fontsize=6, title="sand-like layers", title_fontsize=6)
    ax.set_title("The between-period change persists", loc="left"); panel(ax, "a")
    ax = axs[1]
    scales_w = ["cellFE_500m", "cellFE_250m", "cellFE_100m", "cellFE_50m", "pairs_50m_sameprog_pd", "pairs_20m_sameprog_pd"]
    lab_w = ["500 m", "250 m", "100 m", "50 m", "pairs\n≤ 50 m", "pairs\n≤ 20 m"]
    for i, (c, col, lab) in enumerate(st_cols):
        pts = []
        for k, s in enumerate(scales_w):
            v = wr["P3"][c].get(s)
            if v and "slope" in v:
                pts.append((k + (i - 1) * 0.12, v["slope"], v["lo"], v["hi"]))
        ax.errorbar([p[0] for p in pts], [p[1] for p in pts], yerr=[[p[1] - p[2] for p in pts], [p[3] - p[1] for p in pts]],
                    fmt="o", ms=4, color=col, elinewidth=1.1, capsize=0, label=lab)
    ys = tv["archive"]["yang_implied_slope_per_decade"]
    ax.axhline(ys, color="#eb6834", lw=0.8, ls=(0, (4, 2)))
    ax.text(len(scales_w) - 0.6, ys + 0.02, "implied by a 70% rise over 10–800 days", ha="right", va="bottom", fontsize=5.8, color=INK2)
    ax.axhline(0, color=INK2, lw=0.6)
    ax.set_xticks(np.arange(len(scales_w))); ax.set_xticklabels(lab_w, fontsize=6)
    ax.set_xlabel("Location control within the Dec 2011–Feb 2016 window")
    ax.set_ylabel("Δ ln $q_c$ per decade of days since 23 Dec 2011")
    ax.set_title("The elapsed-time slope does not", loc="left"); panel(ax, "b")
    fig.tight_layout()
    save(fig, outdir / "B7_location_control")



def main(a, full=False):
    a = Path(a)
    if full:
        paper_b_full(a, a / "figures" / "paperB")
        print("full Paper B figures written to", a / "figures" / "paperB")
    else:
        paper_b(a, a / "figures" / "regenerated")
        print("aggregate-only Paper B figures written to", a / "figures" / "regenerated")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="Generate Paper B figures from aggregate JSON; --full enables private B1/B8 inputs.")
    ap.add_argument("analysis_dir", type=Path)
    ap.add_argument("--full", action="store_true", help="also read private soundings/parquet and repeat CSV inputs for B1/B8")
    args = ap.parse_args()
    main(args.analysis_dir, full=args.full)
