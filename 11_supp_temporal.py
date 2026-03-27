#!/usr/bin/env python3
"""
Supplementary Figure: Q1→Q4 Temporal Community Dynamics
ISS Chile Pepper Rhizosphere Microbiome

Panels:
  A — Faith's Phylogenetic Diversity (Plant Roots + Rhizosphere, Q1-Q4)
  B — Shannon Entropy              (Plant Roots + Rhizosphere, Q1-Q4)
  C — PCoA trajectory              (all Space_Flight, Q1-Q4, centroid arrows)
  D — Phylum composition           (Plant Roots + Rhizosphere, Q1-Q4 stacked bar)

Output: version-2_integrated/SuppS3_Temporal_Q1_Q4.pdf / .png
Corresponds to: Supplementary Figure S2 in the manuscript

Run from repository root:
  conda activate qiime2-amplicon
  python 14_supp_temporal.py
"""

from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
from scipy.stats import kruskal

# ── Paths ──────────────────────────────────────────────────────────────────
BASE   = str(Path(__file__).parent / 'version-2_integrated')
META   = f"{BASE}/integrated_metadata.tsv"
ALPHA  = f"{BASE}/exported_diversity/alpha-diversity.tsv"   # faith_pd
SHAN   = f"{BASE}/exported_diversity/shannon.tsv"
ORDIN  = f"{BASE}/exported_diversity/ordination.txt"
FTAB   = f"{BASE}/exported_table_clean/feature-table.tsv"
TAX    = f"{BASE}/exported_taxonomy/taxonomy.tsv"
OUTDIR = str(Path(__file__).parent / 'version-2_integrated')
Path(OUTDIR).mkdir(parents=True, exist_ok=True)

Q_ORDER  = ["Q1", "Q2", "Q3", "Q4"]
Q_COLORS = {"Q1": "#2196F3", "Q2": "#4CAF50", "Q3": "#FF9800", "Q4": "#E91E63"}
RHIZO_TYPES = {"Plant Roots", "Rhizosphere"}   # n=5 per Q cycle


# ══════════════════════════════════════════════════════════════════════════
# Data loading
# ══════════════════════════════════════════════════════════════════════════

def load_metadata():
    meta = pd.read_csv(META, sep="\t", index_col=0)
    space = meta[meta["Study_Group"] == "Space_Flight"].copy()
    rhizo = space[space["sample_type"].isin(RHIZO_TYPES)].copy()
    return space, rhizo


def load_alpha(rhizo_meta):
    faith = pd.read_csv(ALPHA, sep="\t", index_col=0)
    faith.columns = ["faith_pd"]
    shan  = pd.read_csv(SHAN,  sep="\t", index_col=0)
    shan.columns  = ["shannon"]
    alpha = faith.join(shan, how="outer")
    return rhizo_meta.join(alpha, how="inner")


def parse_ordination(path):
    """Parse QIIME2 Bray-Curtis PCoA ordination.txt → (DataFrame, proportions)."""
    with open(path) as f:
        lines = f.readlines()
    prop_idx = next(i for i, l in enumerate(lines) if l.startswith("Proportion explained"))
    props = list(map(float, lines[prop_idx + 1].strip().split("\t")))
    site_idx = next(i for i, l in enumerate(lines) if l.startswith("Site\t"))
    n = int(lines[site_idx].split("\t")[1])
    coords = {}
    for line in lines[site_idx + 1: site_idx + 1 + n]:
        parts = line.strip().split("\t")
        coords[parts[0]] = list(map(float, parts[1:]))
    df = pd.DataFrame(coords).T
    df.columns = [f"PC{i+1}" for i in range(df.shape[1])]
    return df, props


def load_phylum_by_q(rhizo_meta):
    """Return mean phylum relative abundance per Q cycle (top 8 + Others)."""
    tax = pd.read_csv(TAX, sep="\t", index_col=0)
    tax["Phylum"] = tax["Taxon"].str.extract(r"p__([\w\-]+)")
    tax["Phylum"] = tax["Phylum"].fillna("Unassigned")

    ftab = pd.read_csv(FTAB, sep="\t", skiprows=1, index_col=0)
    rhizo_ids = [s for s in rhizo_meta.index if s in ftab.columns]
    rel = ftab[rhizo_ids].div(ftab[rhizo_ids].sum(), axis=1) * 100
    rel = rel.join(tax["Phylum"])
    phylum = rel.groupby("Phylum").sum()

    top8 = phylum.sum(axis=1).nlargest(8).index.tolist()
    result = phylum.loc[top8].copy()
    result.loc["Others"] = phylum.loc[~phylum.index.isin(top8)].sum()

    q_means = {}
    for q in Q_ORDER:
        ids = [s for s in rhizo_meta[rhizo_meta["location"] == q].index if s in result.columns]
        if ids:
            q_means[q] = result[ids].mean(axis=1)
    return pd.DataFrame(q_means)[Q_ORDER]


def kw_label(groups):
    """Return Kruskal-Wallis p-value string from list of arrays."""
    valid = [g for g in groups if len(g) >= 2]
    if len(valid) < 2:
        return "n.s."
    stat, p = kruskal(*valid)
    if p < 0.001:
        return "KW p < 0.001"
    elif p < 0.05:
        return f"KW p = {p:.3f}"
    else:
        return f"KW p = {p:.3f} (n.s.)"


# ══════════════════════════════════════════════════════════════════════════
# Plotting
# ══════════════════════════════════════════════════════════════════════════

PHYLUM_COLORS = [
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
    "#FF7F00", "#A65628", "#F781BF", "#999999", "#CCCCCC"
]
FONT  = dict(fontsize=9)
TFONT = dict(fontsize=10, fontweight="bold")


def add_boxpanel(ax, alpha_df, col, ylabel, title_label):
    groups = {q: alpha_df.loc[alpha_df["location"] == q, col].dropna().values
              for q in Q_ORDER}
    bp = ax.boxplot(
        [groups[q] for q in Q_ORDER],
        positions=range(len(Q_ORDER)),
        patch_artist=True, widths=0.5,
        medianprops=dict(color="black", linewidth=1.5),
        whiskerprops=dict(linewidth=1), capprops=dict(linewidth=1),
        flierprops=dict(marker="o", markersize=3, alpha=0.5)
    )
    for patch, q in zip(bp["boxes"], Q_ORDER):
        patch.set_facecolor(Q_COLORS[q])
        patch.set_alpha(0.8)
    # Overlay individual points
    for i, q in enumerate(Q_ORDER):
        jitter = np.random.uniform(-0.1, 0.1, size=len(groups[q]))
        ax.scatter(i + jitter, groups[q], color=Q_COLORS[q],
                   s=18, zorder=5, alpha=0.7, edgecolors="white", linewidths=0.3)
    ax.text(0.97, 0.96, kw_label(list(groups.values())),
            transform=ax.transAxes, ha="right", va="top", **FONT)
    ax.set_xticks(range(len(Q_ORDER)))
    ax.set_xticklabels(Q_ORDER, **FONT)
    ax.set_ylabel(ylabel, **FONT)
    ax.set_title(title_label, loc="left", **TFONT)
    ax.spines[["top", "right"]].set_visible(False)


def plot_pcoa(ax, space_meta, pcoa_df, props):
    space_pcoa = space_meta.join(pcoa_df[["PC1", "PC2"]], how="inner")
    for q in Q_ORDER:
        pts = space_pcoa[space_pcoa["location"] == q]
        ax.scatter(pts["PC1"], pts["PC2"],
                   c=Q_COLORS[q], label=q, s=25, alpha=0.65, edgecolors="none")
    # Centroid trajectory
    centroids = {
        q: space_pcoa[space_pcoa["location"] == q][["PC1", "PC2"]].mean()
        for q in Q_ORDER
    }
    for i in range(len(Q_ORDER) - 1):
        q1, q2 = Q_ORDER[i], Q_ORDER[i + 1]
        if q1 in centroids and q2 in centroids:
            c1, c2 = centroids[q1], centroids[q2]
            ax.annotate("", xy=(c2["PC1"], c2["PC2"]), xytext=(c1["PC1"], c1["PC2"]),
                        arrowprops=dict(arrowstyle="-|>", color="gray",
                                        lw=1.5, mutation_scale=12))
    ax.set_xlabel(f"PC1 ({props[0]*100:.1f}%)", **FONT)
    ax.set_ylabel(f"PC2 ({props[1]*100:.1f}%)", **FONT)
    ax.set_title("C   PCoA Trajectory (All Space Flight)", loc="left", **TFONT)
    ax.legend(title="Cycle", fontsize=8, title_fontsize=8,
              frameon=False, handletextpad=0.3)
    ax.spines[["top", "right"]].set_visible(False)


def plot_phylum(ax, phylum_q):
    bottoms = np.zeros(len(Q_ORDER))
    for i, phylum in enumerate(phylum_q.index):
        vals = phylum_q.loc[phylum].values
        ax.bar(range(len(Q_ORDER)), vals, bottom=bottoms,
               color=PHYLUM_COLORS[i % len(PHYLUM_COLORS)],
               label=phylum, width=0.65)
        bottoms += vals
    ax.set_xticks(range(len(Q_ORDER)))
    ax.set_xticklabels(Q_ORDER, **FONT)
    ax.set_ylabel("Mean Relative Abundance (%)", **FONT)
    ax.set_ylim(0, 105)
    ax.set_title("D   Phylum Composition (Rhizosphere)", loc="left", **TFONT)
    ax.legend(loc="upper left", fontsize=7, frameon=False,
              bbox_to_anchor=(1.01, 1.0))
    ax.spines[["top", "right"]].set_visible(False)


# ══════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════

def main():
    print("Loading data...")
    space_meta, rhizo_meta = load_metadata()
    alpha_df  = load_alpha(rhizo_meta)
    pcoa_df, props = parse_ordination(ORDIN)
    phylum_q  = load_phylum_by_q(rhizo_meta)

    print(f"  Rhizosphere samples: {len(rhizo_meta)} "
          f"({dict(rhizo_meta['location'].value_counts())})")
    print(f"  All Space_Flight:    {len(space_meta)}")

    np.random.seed(42)
    fig = plt.figure(figsize=(15, 9))
    fig.patch.set_facecolor("white")
    gs = GridSpec(2, 2, figure=fig, hspace=0.40, wspace=0.38)

    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    add_boxpanel(ax_a, alpha_df, "faith_pd",
                 "Faith's Phylogenetic Diversity",
                 "A   Phylogenetic Diversity (Rhizosphere)")

    add_boxpanel(ax_b, alpha_df, "shannon",
                 "Shannon Entropy",
                 "B   Alpha Diversity (Rhizosphere)")

    plot_pcoa(ax_c, space_meta, pcoa_df, props)

    plot_phylum(ax_d, phylum_q)

    fig.suptitle(
        "Temporal dynamics of the ISS rhizosphere microbiome across planting cycles (Q1–Q4)",
        fontsize=11, fontweight="bold", y=0.98
    )

    for fmt in ("pdf", "png"):
        path = f"{OUTDIR}/SuppS3_Temporal_Q1_Q4.{fmt}"
        fig.savefig(path, dpi=300, bbox_inches="tight")
        print(f"Saved: {path}")

    # Current Biology supplementary figure output
    import config
    cb_dir = config.CURRENT_BIOLOGY_DIR
    cb_dir.mkdir(exist_ok=True)
    fig.savefig(cb_dir / 'SuppS1d_Temporal.png', dpi=300, bbox_inches='tight')
    fig.savefig(cb_dir / 'SuppS1d_Temporal.pdf', bbox_inches='tight')
    print(f"Saved Current Biology: {cb_dir / 'SuppS1d_Temporal.png'}")


if __name__ == "__main__":
    main()
