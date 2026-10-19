"""
38j_lettuce_diversity_pcoa_figures.py
======================================
Publication-quality figures for lettuce α-diversity and β-diversity (PCoA).

Inputs (from QIIME2 core-metrics-phylogenetic depth=100):
  - faith_pd_vector/alpha-diversity.tsv
  - shannon_vector/alpha-diversity.tsv
  - evenness_vector/alpha-diversity.tsv
  - observed_features_vector/alpha-diversity.tsv
  - bray_curtis_pcoa_results/ordination.txt
  - weighted_unifrac_pcoa_results/ordination.txt
  - unweighted_unifrac_pcoa_results/ordination.txt

Outputs (acta_astronautica/lettuce/):
  - Supp_lettuce_alpha_diversity.{pdf,png}  — 4-panel violin+strip
  - Supp_lettuce_pcoa.{pdf,png}              — 3-panel PCoA (BC/wUF/uUF)
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import mannwhitneyu

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config

LETTUCE_DIR = config.VEGGIE_LETTUCE_DIR
CORE_DIR = LETTUCE_DIR / "core-metrics-d100"
OUT_DIR = config.AA_DIR / "lettuce"

COLORS = {"Space_Flight": "#3498db", "Terrestrial_Soil": "#2ecc71"}
MARKERS = {"OSD-267": "o", "OSD-269": "s"}
GROUP_ORDER = ["Space_Flight", "Terrestrial_Soil"]


def cliffs_delta(a, b):
    a = np.asarray(a); b = np.asarray(b)
    if len(a) == 0 or len(b) == 0:
        return np.nan
    dom = np.sum(a[:, None] > b[None, :]) - np.sum(a[:, None] < b[None, :])
    return dom / (len(a) * len(b))


def hedges_g(a, b):
    a = np.asarray(a); b = np.asarray(b)
    if len(a) < 2 or len(b) < 2:
        return np.nan
    na, nb = len(a), len(b)
    sa2, sb2 = np.var(a, ddof=1), np.var(b, ddof=1)
    sp = np.sqrt(((na - 1) * sa2 + (nb - 1) * sb2) / (na + nb - 2))
    if sp == 0:
        return np.nan
    d = (np.mean(a) - np.mean(b)) / sp
    j = 1 - 3 / (4 * (na + nb) - 9)
    return d * j


def load_alpha(metric_dir, metric_name):
    """QIIME2 alpha-diversity exports a TSV with header 'sample-id\t{metric}'."""
    path = CORE_DIR / metric_dir / "alpha-diversity.tsv"
    df = pd.read_csv(path, sep="\t", index_col=0)
    col = df.columns[0]
    return df[col].rename(metric_name)


def load_pcoa(pcoa_dir):
    """Parse QIIME2 ordination.txt (skoot format). Returns DataFrame of coords + proportion explained."""
    path = CORE_DIR / pcoa_dir / "ordination.txt"
    lines = path.read_text().splitlines()

    # Parse proportion explained
    prop_idx = lines.index("Proportion explained\t" + lines[lines.index("Proportion explained\t"
                                    ) + 0].split("\t", 1)[-1]) if False else None
    # Simpler: scan for sections
    sections = {}
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith("Proportion explained"):
            _, n = line.split("\t")
            n = int(n)
            sections["proportion"] = list(map(float, lines[i+1].strip().split("\t")))
            i += 2
        elif line.startswith("Site"):
            _, nrow, ncol = line.split("\t")
            nrow = int(nrow); ncol = int(ncol)
            rows = []
            for j in range(i + 1, i + 1 + nrow):
                parts = lines[j].strip().split("\t")
                rows.append([parts[0]] + list(map(float, parts[1:])))
            coords = pd.DataFrame(rows, columns=["sample_id"] + [f"PC{k+1}"
                                                                  for k in range(ncol)])
            coords.set_index("sample_id", inplace=True)
            sections["coords"] = coords
            i = i + 1 + nrow
        else:
            i += 1
    return sections


def annotate_effect(ax, sp, tr, ylim_top=None):
    """Place Cliff's δ + Hedges g + MWU p below title."""
    if len(sp) == 0 or len(tr) == 0:
        return
    cd = cliffs_delta(sp, tr)
    hg = hedges_g(sp, tr)
    try:
        _, p = mannwhitneyu(sp, tr, alternative="two-sided")
        p_str = f"p={p:.3g}"
    except ValueError:
        p_str = "p=NA"
    txt = f"Cliff's δ={cd:.2f}, g={hg:.2f}, {p_str}"
    ax.text(0.5, 0.98, txt, transform=ax.transAxes,
            ha="center", va="top", fontsize=8,
            bbox=dict(boxstyle="round,pad=0.3",
                      facecolor="white", edgecolor="gray", alpha=0.8))


def alpha_figure(meta):
    metrics = [
        ("shannon_vector", "Shannon entropy", "Shannon"),
        ("observed_features_vector", "Observed ASVs", "Observed ASVs"),
        ("faith_pd_vector", "Faith's phylogenetic diversity", "Faith's PD"),
        ("evenness_vector", "Pielou's evenness", "Evenness"),
    ]

    data = meta.copy()
    for met_dir, met_name, col in metrics:
        try:
            s = load_alpha(met_dir, col)
            data[col] = s.reindex(data.index)
        except FileNotFoundError:
            print(f"  Missing: {met_dir}")
            data[col] = np.nan

    # 4-panel
    fig, axes = plt.subplots(1, 4, figsize=(15, 4.5), constrained_layout=True)
    for ax, (_, met_name, col) in zip(axes, metrics):
        if col not in data.columns or data[col].isna().all():
            ax.axis("off"); continue
        d = data.dropna(subset=[col])
        # violin
        sns.violinplot(data=d, x="Study_Group", y=col, order=GROUP_ORDER,
                       hue="Study_Group", palette=COLORS, ax=ax,
                       inner=None, alpha=0.5, legend=False, cut=0)
        # strip with mission shapes
        for mission, marker in MARKERS.items():
            sub = d[d["OSD"] == mission]
            if len(sub) == 0:
                continue
            for grp in GROUP_ORDER:
                gg = sub[sub["Study_Group"] == grp]
                if len(gg) == 0:
                    continue
                x_pos = GROUP_ORDER.index(grp) + np.random.uniform(-0.1, 0.1, len(gg))
                ax.scatter(x_pos, gg[col], marker=marker,
                           s=70, edgecolor="black", linewidth=0.7,
                           facecolor=COLORS[grp], alpha=0.9,
                           label=f"{mission}" if grp == GROUP_ORDER[0] else None)
        # effect size annotation
        sp = d.loc[d["Study_Group"] == "Space_Flight", col].values
        tr = d.loc[d["Study_Group"] == "Terrestrial_Soil", col].values
        annotate_effect(ax, sp, tr)
        ax.set_title(met_name, fontsize=10, fontweight="bold")
        ax.set_xlabel("")
        ax.set_ylabel(met_name, fontsize=9)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Space\n(ISS Veggie)", "Ground\n(Veggie control)"], fontsize=9)
        ax.grid(axis="y", alpha=0.3)

    # Mission legend on far right
    mission_handles = [plt.Line2D([0], [0], marker=MARKERS[m], color="w",
                                  markerfacecolor="gray", markeredgecolor="black",
                                  markersize=10, label=f"{m} ({'VEG-01A' if '267' in m else 'VEG-03A'})")
                       for m in ["OSD-267", "OSD-269"]]
    axes[-1].legend(handles=mission_handles, loc="lower right", fontsize=8,
                    title="Mission", frameon=True)

    fig.suptitle("Lettuce α-diversity (QIIME2 core-metrics, rarefaction depth 100, n=4+4)",
                 fontsize=11, fontweight="bold", y=1.02)
    for ext in ("pdf", "png"):
        fig.savefig(OUT_DIR / f"Supp_lettuce_alpha_diversity.{ext}",
                    dpi=300 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: Supp_lettuce_alpha_diversity.pdf/.png")
    return data


def pcoa_figure(meta, beta_df):
    """Three-panel PCoA (BC, weighted UniFrac, unweighted UniFrac)."""
    metrics = [
        ("bray_curtis_pcoa_results", "Bray–Curtis"),
        ("weighted_unifrac_pcoa_results", "Weighted UniFrac"),
        ("unweighted_unifrac_pcoa_results", "Unweighted UniFrac"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.3), constrained_layout=True)
    for ax, (dir_, title) in zip(axes, metrics):
        try:
            sect = load_pcoa(dir_)
        except Exception as e:
            print(f"  Skip {dir_}: {e}")
            ax.axis("off"); continue
        coords = sect["coords"]
        prop = sect.get("proportion", [0.0, 0.0])
        merged = coords.join(meta, how="inner")
        for grp in GROUP_ORDER:
            for mission, marker in MARKERS.items():
                sub = merged[(merged["Study_Group"] == grp) & (merged["OSD"] == mission)]
                if len(sub) == 0:
                    continue
                ax.scatter(sub["PC1"], sub["PC2"],
                           color=COLORS[grp], marker=marker, s=110,
                           edgecolor="black", linewidth=1.0, alpha=0.9,
                           label=f"{grp.split('_')[0]} ({mission})")
        ax.axhline(0, color="gray", lw=0.5, alpha=0.5)
        ax.axvline(0, color="gray", lw=0.5, alpha=0.5)
        ax.set_xlabel(f"PC1 ({prop[0]*100:.1f}%)", fontsize=9)
        ax.set_ylabel(f"PC2 ({prop[1]*100:.1f}%)", fontsize=9)
        ax.set_title(title, fontsize=10, fontweight="bold")
        # PERMANOVA annotation if available
        if beta_df is not None:
            row = beta_df[beta_df["Subset"] == "aggregate_n4x2"]
            if len(row):
                r2 = row.iloc[0].get("PERMANOVA_R2", "NA")
                p = row.iloc[0].get("PERMANOVA_p", "NA")
                bd = row.iloc[0].get("Betadisper_p", "NA")
                ax.text(0.02, 0.98, f"PERMANOVA\nR²={r2}, p={p}\nbetadisper p={bd}",
                        transform=ax.transAxes, ha="left", va="top", fontsize=7.5,
                        bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                                  edgecolor="gray", alpha=0.85))
        ax.grid(alpha=0.25)

    # Combined legend
    all_h = []
    for grp in GROUP_ORDER:
        for mission, marker in MARKERS.items():
            all_h.append(plt.Line2D([0], [0], marker=marker, linestyle="",
                                    markerfacecolor=COLORS[grp], markeredgecolor="black",
                                    markersize=10,
                                    label=f"{grp.split('_')[0]} — {mission}"))
    axes[-1].legend(handles=all_h, loc="lower right", fontsize=7.5,
                    frameon=True, ncol=1)

    fig.suptitle("Lettuce β-diversity PCoA (QIIME2 core-metrics, rarefaction depth 100, n=4+4)",
                 fontsize=11, fontweight="bold", y=1.03)
    for ext in ("pdf", "png"):
        fig.savefig(OUT_DIR / f"Supp_lettuce_pcoa.{ext}",
                    dpi=300 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: Supp_lettuce_pcoa.pdf/.png")


def main():
    print("=" * 60)
    print("38j: Lettuce α-diversity + PCoA figures")
    print("=" * 60)

    meta = pd.read_csv(LETTUCE_DIR / "metadata_harmonized.tsv", sep="\t", index_col=0)
    print(f"  Metadata: {len(meta)} samples")

    beta_df = None
    try:
        beta_df = pd.read_csv(OUT_DIR / "lettuce_beta_permanova.csv")
    except FileNotFoundError:
        pass

    print("\n[1] α-diversity figure...")
    alpha_data = alpha_figure(meta)
    alpha_data.to_csv(OUT_DIR / "lettuce_alpha_all_metrics.csv")

    print("\n[2] PCoA figure...")
    pcoa_figure(meta, beta_df)

    print("\nDone.")


if __name__ == "__main__":
    main()
