"""
Fig 3 — Root cause: functional redundancy loss with genome-level evidence.
(a) contributing-genera count per function
(b) GTDB N-cycle gene presence/absence matrix (script 36)
(c) PGP Shannon entropy
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import pandas as pd
import seaborn as sns
import config

OUT = config.AA_FIG_DIR
VD = config.VERSION_DIR


def panel_a_redundancy(ax):
    ax.imshow(mpimg.imread(VD / "Fig1d_FunctionalRedundancy.png"))
    ax.set_axis_off()
    ax.set_title("a) Contributing genera per function (pepper)",
                 loc="left", fontweight="bold")


def panel_b_gtdb_matrix(ax):
    mat = pd.read_csv(config.AA_GTDB_MATRIX_CSV, index_col=0)
    # Order: nitrifiers first, then denitrifiers, then generalists
    nitrif = ["Nitrosomonas", "Nitrosospira", "Nitrospira", "Nitrobacter"]
    denitrif = ["Rhodanobacter", "Pseudomonas", "Burkholderia", "Variovorax",
                "Sphingomonas", "Azospirillum", "Paenibacillus"]
    others = ["Bradyrhizobium", "Bacillus"]
    order = [g for g in nitrif + denitrif + others if g in mat.index]
    mat = mat.loc[order]

    cmap = sns.color_palette(["#f5f5f5", "#08519c"], as_cmap=False)
    sns.heatmap(mat, ax=ax, cmap=["#f5f5f5", "#08519c"], cbar=False,
                linewidths=0.6, linecolor="grey",
                annot=True, fmt="d", annot_kws={"fontsize": 10})
    ax.set_xlabel("N-cycle gene (product-description)")
    ax.set_ylabel("")
    ax.set_title("b) GTDB/NCBI genome-level gene evidence",
                 loc="left", fontweight="bold")
    for tick in ax.get_yticklabels():
        tick.set_fontstyle("italic")


def panel_c_pgp(ax):
    img = VD / "Fig1c_NCycle_Heatmap.png"
    # The existing NCycle heatmap already bundles PGP entropy; but for clarity
    # we show the PGP entropy CSV as a bar plot.
    try:
        df = pd.read_csv(config.VERSION_DIR / "Fig1c_PGP_Specificity_Results.csv")
        if "Study_Group" in df.columns and "Shannon_Entropy" in df.columns:
            means = df.groupby("Study_Group")["Shannon_Entropy"].mean()
            ax.bar(means.index, means.values,
                   color=["#e34a33", "#2b8cbe"])
            ax.set_ylabel("Shannon entropy (PGP specificity)")
            ax.set_title("c) PGP specificity entropy",
                         loc="left", fontweight="bold")
            for i, v in enumerate(means.values):
                ax.text(i, v, f"{v:.2f}", ha="center", va="bottom")
        else:
            # Fallback — show any column ending with entropy
            num_cols = df.select_dtypes("number").columns.tolist()
            if num_cols:
                df.boxplot(column=num_cols[0], ax=ax)
                ax.set_title("c) PGP specificity", loc="left", fontweight="bold")
            else:
                ax.set_axis_off()
    except FileNotFoundError:
        ax.set_axis_off()


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(16, 6))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.1, 1.5, 0.9], wspace=0.35)
    panel_a_redundancy(fig.add_subplot(gs[0, 0]))
    panel_b_gtdb_matrix(fig.add_subplot(gs[0, 1]))
    panel_c_pgp(fig.add_subplot(gs[0, 2]))
    fig.suptitle("Figure 3. Functional redundancy loss with genome-level gene evidence",
                 fontsize=13, fontweight="bold", y=1.02)
    out = OUT / "Fig3_root_cause"
    fig.savefig(str(out) + ".png", dpi=300, bbox_inches="tight")
    fig.savefig(str(out) + ".pdf", bbox_inches="tight")
    print(f"Wrote {out}.png/.pdf")


if __name__ == "__main__":
    main()
