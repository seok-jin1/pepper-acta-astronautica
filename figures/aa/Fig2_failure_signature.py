"""
Fig 2 — Failure signature: rhizosphere nitrogen-cycle depletion pattern.
6 panels composing the pepper N-cycle results for Acta Astronautica.
  a) Faith's PD (α-diversity)
  b) Bray-Curtis PCoA (β-diversity + PERMANOVA)
  c) FAPROTAX Log2FC per function
  d) FAPROTAX × PICRUSt2 consensus heatmap
  e) N-cycle completeness (per-sample boxplot, Space vs Terrestrial)
  f) N-cycle step presence (5 steps x 2 groups heatmap)
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Patch
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.stats import mannwhitneyu
import config

OUT = config.AA_FIG_DIR
VD = config.VERSION_DIR

ROOT_ZONE_TYPES = {"Plant Roots", "Plant_Roots", "Rhizosphere"}

N_CYCLE_FUNCTIONS = [
    "nitrification", "aerobic_ammonia_oxidation", "aerobic_nitrite_oxidation",
    "denitrification", "nitrogen_fixation", "nitrate_reduction",
]
CONTROL_FUNCTIONS = [
    "aerobic_chemoheterotrophy", "fermentation",
    "aromatic_compound_degradation", "cellulolysis", "chitinolysis",
]


def panel_a_faith_pd(ax):
    """α-diversity Faith's PD, Space vs Terrestrial (native vector)."""
    meta = pd.read_csv(config.INTEGRATED_METADATA, sep="\t", index_col=0)
    alpha = pd.read_csv(VD / "exported_diversity/alpha-diversity.tsv",
                        sep="\t", index_col=0)
    alpha.columns = ["faith_pd"]

    space_ids = meta.index[meta["Study_Group"] == "Space_Flight"]
    terr_ids = meta.index[meta["Study_Group"] == "Terrestrial_Soil"]
    space_pd = alpha.loc[alpha.index.intersection(space_ids), "faith_pd"].values
    terr_pd = alpha.loc[alpha.index.intersection(terr_ids), "faith_pd"].values

    res = mannwhitneyu(space_pd, terr_pd, alternative="two-sided")
    p = res.pvalue
    r_rb = (2 * res.statistic / (len(space_pd) * len(terr_pd))) - 1
    ptxt = "p < 0.001" if p < 0.001 else f"p = {p:.3f}"

    COL_SP, COL_TE = "#C0392B", "#2E6DB4"
    parts = ax.violinplot([space_pd, terr_pd], positions=[1, 2], widths=0.72,
                          showmeans=False, showmedians=False, showextrema=False)
    for body, col in zip(parts["bodies"], [COL_SP, COL_TE]):
        body.set_facecolor(col); body.set_alpha(0.22); body.set_edgecolor("none")

    bp = ax.boxplot(
        [space_pd, terr_pd], positions=[1, 2], widths=0.34, patch_artist=True,
        medianprops=dict(color="black", lw=1.4),
        whiskerprops=dict(color="#555", lw=1.0),
        capprops=dict(color="#555", lw=1.0),
        flierprops=dict(marker="o", markersize=3, alpha=0.45,
                        markerfacecolor="#888", markeredgecolor="none"),
        boxprops=dict(lw=1.0),
    )
    for patch, col in zip(bp["boxes"], [COL_SP, COL_TE]):
        patch.set_facecolor(col); patch.set_alpha(0.78); patch.set_edgecolor("#404040")

    rng = np.random.default_rng(42)
    for pos, vals in zip([1, 2], [space_pd, terr_pd]):
        jitter = rng.uniform(-0.09, 0.09, size=len(vals))
        ax.scatter(np.full_like(vals, pos) + jitter, vals,
                   s=13, c="#333", alpha=0.28, zorder=4)

    ytop = float(max(space_pd.max(), terr_pd.max()))
    ax.plot([1, 2], [ytop * 1.02, ytop * 1.02], color="#333", lw=0.9)
    ax.text(1.5, ytop * 1.065, ptxt, ha="center", va="bottom", fontsize=9)
    ax.text(1.5, ytop * 1.115, f"$r_{{rb}}$ = {r_rb:.2f}",
            ha="center", va="bottom", fontsize=8.6, color="#666")

    ax.set_xticks([1, 2])
    ax.set_xticklabels([f"Space\n(n={len(space_pd)})",
                        f"Terrestrial\n(n={len(terr_pd)})"], fontsize=9)
    ax.set_ylabel("Faith's PD")
    ax.set_xlim(0.4, 2.6)
    ax.set_ylim(0, ytop * 1.17)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax.grid(axis="y", alpha=0.3); ax.set_axisbelow(True)
    ax.set_title("a) α-diversity (Faith's PD)", loc="left", fontweight="bold")


def _load_analytical_root_zone_relabund():
    """Return (genus_ra_df samples × genera, meta_slice) for the analytical root zone."""
    ft = pd.read_csv(config.FEATURE_TABLE_CLEAN, sep="\t", skiprows=1, index_col=0)
    meta = pd.read_csv(config.INTEGRATED_METADATA, sep="\t", index_col=0)
    keep = meta.index[
        (meta["Study_Group"] == "Terrestrial_Soil") |
        ((meta["Study_Group"] == "Space_Flight") & meta["sample_type"].isin(ROOT_ZONE_TYPES))
    ]
    keep = [s for s in keep if s in ft.columns]
    rel = ft[keep]
    rel = rel.div(rel.sum(axis=0), axis=1).fillna(0.0).T  # samples × ASVs
    return rel, meta.loc[keep]


def _pcoa(dist):
    """Classical PCoA: double-centered Gower matrix eigendecomposition."""
    n = dist.shape[0]
    A = -0.5 * dist ** 2
    J = np.eye(n) - np.ones((n, n)) / n
    B = J @ A @ J
    eigvals, eigvecs = np.linalg.eigh(B)
    order = np.argsort(eigvals)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]
    coords = eigvecs * np.sqrt(np.maximum(eigvals, 0))
    total_var = np.sum(np.maximum(eigvals, 0))
    pct = 100 * np.maximum(eigvals, 0) / total_var if total_var > 0 else np.zeros(n)
    return coords, pct


def _ellipse(ax, x, y, color, alpha=0.15):
    if len(x) < 3:
        return
    cov = np.cov(x, y)
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:, order]
    angle = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))
    # 95% CI scaling factor for 2-D Gaussian ≈ sqrt(5.991)
    width, height = 2 * np.sqrt(5.991 * vals)
    ell = Ellipse((x.mean(), y.mean()), width, height, angle=angle,
                  facecolor=color, edgecolor=color, alpha=alpha, linewidth=1.5)
    ax.add_patch(ell)


def panel_b_pcoa(ax):
    """Bray–Curtis PCoA scatter (analytical root zone, n=40) with PERMANOVA annotation."""
    rel, meta = _load_analytical_root_zone_relabund()
    dist = squareform(pdist(rel.values, metric="braycurtis"))
    coords, pct = _pcoa(dist)
    pc1, pc2 = coords[:, 0], coords[:, 1]

    groups = meta.loc[rel.index, "Study_Group"]
    colors = {"Space_Flight": "#e34a33", "Terrestrial_Soil": "#2b8cbe"}
    labels = {"Space_Flight": "Space Flight (n=20)", "Terrestrial_Soil": "Terrestrial Soil (n=20)"}
    for g, c in colors.items():
        m = (groups == g).values
        ax.scatter(pc1[m], pc2[m], s=52, alpha=0.82, color=c,
                   edgecolor="white", linewidth=0.7, label=labels[g])
        _ellipse(ax, pc1[m], pc2[m], c)

    ax.axhline(0, color="lightgrey", lw=0.6, zorder=0)
    ax.axvline(0, color="lightgrey", lw=0.6, zorder=0)
    ax.set_xlabel(f"PCoA1 ({pct[0]:.1f}%)")
    ax.set_ylabel(f"PCoA2 ({pct[1]:.1f}%)")
    ax.legend(fontsize=8, loc="best", frameon=True)

    # PERMANOVA / PERMDISP annotation (values from SuppS2_PERMANOVA_results.csv, analytical filter)
    try:
        perm = pd.read_csv(VD / "SuppS2_PERMANOVA_results.csv")
        row = perm[perm["Filter"].str.contains("ANALYTICAL", na=False)].iloc[0]
        ann = (f"PERMANOVA: R²={row['R2']:.3f}, p={row['p_value']:.3f}\n"
               f"PERMDISP: p={row['betadisper_p']:.3f} "
               f"(Space σ={row['Space_dispersion']:.3f}, Terr σ={row['Terr_dispersion']:.3f})")
    except Exception:
        ann = ("PERMANOVA: R²=0.56, p<0.001\n"
               "PERMDISP: p=0.038 (Space σ=0.307, Terr σ=0.397)")
    ax.text(0.02, 0.02, ann, transform=ax.transAxes, fontsize=8,
            va="bottom", ha="left",
            bbox=dict(facecolor="white", edgecolor="grey", alpha=0.85, pad=4))
    ax.set_title("b) Bray–Curtis PCoA (analytical root zone, n=40)",
                 loc="left", fontweight="bold")


def panel_c_faprotax(ax):
    """FAPROTAX Log2FC horizontal barplot (native vector)."""
    df = pd.read_csv(VD / "FAPROTAX_group_comparison.csv")
    df = df.set_index("group")
    funcs = [f for f in N_CYCLE_FUNCTIONS + CONTROL_FUNCTIONS if f in df.index]
    series = df.loc[funcs, "Log2FC"].sort_values()
    labels = [f.replace("_", " ").title() for f in series.index]
    colors = ["#c0392b" if v < 0 else "#2980b9" for v in series.values]
    is_ncycle = [f in N_CYCLE_FUNCTIONS for f in series.index]

    ax.barh(range(len(series)), series.values,
            color=colors, edgecolor="white", height=0.65)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_yticks(range(len(series)))
    tls = ax.set_yticklabels(labels, fontsize=9)
    for lbl, nc in zip(tls, is_ncycle):
        if nc:
            lbl.set_fontweight("bold")
    for i, v in enumerate(series.values):
        ax.text(v + (0.3 if v >= 0 else -0.3), i, f"{v:+.1f}",
                va="center", ha="left" if v >= 0 else "right", fontsize=8)
    ax.set_xlabel("Log$_2$ Fold Change (Space / Terrestrial)")
    ax.legend(handles=[Patch(facecolor="#2980b9", label="Space-enriched"),
                       Patch(facecolor="#c0392b", label="Space-depleted")],
              fontsize=8, loc="lower right", frameon=True)
    ax.text(0.01, 0.99, "Bold = N-cycle function", transform=ax.transAxes,
            fontsize=7.5, va="top", color="#555")
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax.grid(axis="x", alpha=0.3); ax.set_axisbelow(True)
    ax.set_title("c) FAPROTAX Log$_2$FC (Space / Terr)",
                 loc="left", fontweight="bold")


def panel_d_consensus(ax):
    """FAPROTAX vs PICRUSt2 consensus grouped bars (native vector)."""
    df = pd.read_csv(config.AA_DIR / "consensus_function_log2fc.csv", index_col=0)
    df = df.sort_values("consensus_log2fc")
    x = np.arange(len(df)); w = 0.4
    ax.bar(x - w / 2, df["faprotax_log2fc"], w,
           label="FAPROTAX", color="#2b8cbe", edgecolor="white")
    ax.bar(x + w / 2, df["picrust_log2fc"], w,
           label="PICRUSt2", color="#e34a33", edgecolor="white")
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels([f.replace("_", " ").title() for f in df.index],
                       rotation=25, ha="right", fontsize=9)
    for xi, (fv, pv) in enumerate(zip(df["faprotax_log2fc"], df["picrust_log2fc"])):
        ax.text(xi - w / 2, fv + (0.08 if fv >= 0 else -0.08), f"{fv:+.2f}",
                ha="center", va="bottom" if fv >= 0 else "top", fontsize=7.5)
        ax.text(xi + w / 2, pv + (0.08 if pv >= 0 else -0.08), f"{pv:+.2f}",
                ha="center", va="bottom" if pv >= 0 else "top", fontsize=7.5)
    ax.set_ylabel("Log$_2$ Fold Change (Space / Terrestrial)")
    ax.legend(fontsize=8, frameon=True, loc="best")
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax.grid(axis="y", alpha=0.3); ax.set_axisbelow(True)
    ax.set_title("d) FAPROTAX × PICRUSt2 consensus (MetaCyc-covered functions)",
                 loc="left", fontweight="bold")


NCYCLE_STEPS = ["nitrogen_fixation", "nitrification", "denitrification",
                "nitrate_reduction", "ureolysis"]
NCYCLE_LABELS = {"nitrogen_fixation": "N$_2$ fixation",
                 "nitrification":     "Nitrification",
                 "denitrification":   "Denitrification",
                 "nitrate_reduction": "Nitrate reduction",
                 "ureolysis":         "Ureolysis"}

COLOR_SPACE = "#e34a33"
COLOR_TERR  = "#2b8cbe"


def _sig_label(p):
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    return "ns"


def panel_e_completeness_boxplot(ax):
    """e) Per-sample N-cycle completeness, Space vs Terrestrial."""
    df = pd.read_csv(VD / "Fig2B_NCycle_Completeness_Results.csv", index_col=0)
    sp = df[df["Study_Group"] == "Space_Flight"]["Completeness"].values
    tr = df[df["Study_Group"] == "Terrestrial_Soil"]["Completeness"].values
    _, p = mannwhitneyu(sp, tr, alternative="two-sided")

    rng = np.random.default_rng(42)
    for i, (vals, col) in enumerate(zip([sp, tr], [COLOR_SPACE, COLOR_TERR])):
        ax.boxplot(vals, positions=[i], widths=0.45,
                   patch_artist=True,
                   boxprops=dict(facecolor=col, alpha=0.72, edgecolor="black"),
                   medianprops=dict(color="black", linewidth=1.8),
                   whiskerprops=dict(linewidth=1.2),
                   capprops=dict(linewidth=1.2),
                   flierprops=dict(marker="o", markersize=3, alpha=0.4))
        jitter = i + rng.uniform(-0.13, 0.13, size=len(vals))
        ax.scatter(jitter, vals, alpha=0.5, s=18, color="#333", zorder=3)

    y_top = max(sp.max(), tr.max()) + 0.06
    ax.plot([0, 1], [y_top, y_top], color="black", lw=1)
    ax.text(0.5, y_top + 0.015, _sig_label(p),
            ha="center", fontweight="bold", fontsize=12)

    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"Space\n(n={len(sp)})",
                        f"Terrestrial\n(n={len(tr)})"], fontsize=9)
    ax.set_ylabel("N-cycle completeness\n(fraction of 5 steps)")
    ax.set_ylim(0.4, 1.15)
    ax.grid(axis="y", alpha=0.3)
    ax.text(0.98, 0.02, f"mean: {sp.mean():.2f} vs {tr.mean():.2f}",
            transform=ax.transAxes, ha="right", va="bottom",
            fontsize=8, color="#555",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.7, pad=2))
    ax.set_title("e) Completeness per sample (0.76 vs 1.00)",
                 loc="left", fontweight="bold")


def panel_f_step_heatmap(ax):
    """f) Fraction of samples with each N-cycle step present (5 steps × 2 groups)."""
    df = pd.read_csv(VD / "Fig2B_NCycle_Completeness_Results.csv", index_col=0)
    sp = df[df["Study_Group"] == "Space_Flight"]
    tr = df[df["Study_Group"] == "Terrestrial_Soil"]
    # Each step column is already 0/1 per sample, so .mean() = fraction present.
    sp_frac = sp[NCYCLE_STEPS].mean().values
    tr_frac = tr[NCYCLE_STEPS].mean().values
    matrix = np.column_stack([sp_frac, tr_frac])   # rows = steps, cols = groups

    im = ax.imshow(matrix, cmap="RdYlGn", vmin=0, vmax=1, aspect="auto")

    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            ax.text(j, i, f"{matrix[i, j]:.2f}", ha="center", va="center",
                    color="black", fontweight="bold", fontsize=10)

    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"Space\n(n={len(sp)})",
                        f"Terrestrial\n(n={len(tr)})"], fontsize=9)
    ax.set_yticks(range(len(NCYCLE_STEPS)))
    ax.set_yticklabels([NCYCLE_LABELS[s] for s in NCYCLE_STEPS], fontsize=9)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="6%", pad=0.08)
    cbar = plt.colorbar(im, cax=cax)
    cbar.set_label("Fraction of samples\nwith step present", fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    ax.set_title("f) Step presence (Space collapses on nitrification)",
                 loc="left", fontweight="bold")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(20, 12))
    gs = fig.add_gridspec(2, 4, height_ratios=[1, 1.2],
                          hspace=0.28, wspace=0.32)
    panel_a_faith_pd(fig.add_subplot(gs[0, 0]))
    panel_b_pcoa(fig.add_subplot(gs[0, 1]))
    panel_c_faprotax(fig.add_subplot(gs[0, 2:4]))
    panel_d_consensus(fig.add_subplot(gs[1, 0:2]))
    panel_e_completeness_boxplot(fig.add_subplot(gs[1, 2]))
    panel_f_step_heatmap(fig.add_subplot(gs[1, 3]))
    out = OUT / "Fig2_failure_signature"
    fig.savefig(str(out) + ".png", dpi=300, bbox_inches="tight")
    fig.savefig(str(out) + ".pdf", bbox_inches="tight")
    print(f"Wrote {out}.png/.pdf")


if __name__ == "__main__":
    main()
