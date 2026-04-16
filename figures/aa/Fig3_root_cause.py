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
import numpy as np
import pandas as pd
import seaborn as sns
import config

OUT = config.AA_FIG_DIR
VD = config.VERSION_DIR


def panel_a_redundancy(ax):
    df = pd.read_csv(VD / "Fig2C_Functional_Redundancy_Summary.csv")
    key_funcs = ["nitrification", "nitrogen_fixation", "denitrification",
                 "nitrate_reduction", "fermentation", "aerobic_chemoheterotrophy"]
    df = df[df["Function"].isin(key_funcs)].set_index("Function").reindex(key_funcs)

    y = np.arange(len(df))
    h = 0.38
    ax.barh(y - h/2, df["Terr_Mean"], h, xerr=df["Terr_SD"],
            label="Terrestrial Soil", color="#2b8cbe",
            error_kw={"elinewidth": 0.8, "capsize": 2})
    ax.barh(y + h/2, df["Space_Mean"], h, xerr=df["Space_SD"],
            label="Space Flight", color="#e34a33",
            error_kw={"elinewidth": 0.8, "capsize": 2})
    ax.set_yticks(y); ax.set_yticklabels(df.index, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("# contributing genera (mean ± SD)", fontsize=9)
    # Significance stars
    for i, (_, row) in enumerate(df.iterrows()):
        sig = str(row.get("sig", "")).strip()
        if sig and sig != "nan":
            xmax = max(float(row["Terr_Mean"] + row["Terr_SD"]),
                       float(row["Space_Mean"] + row["Space_SD"]))
            ax.text(xmax + 0.3, i, sig, va="center", fontsize=10, color="red")
    ax.legend(fontsize=8, loc="lower right")
    ax.set_title("a) Contributing genera per function (pepper, BH q < 0.05)",
                 loc="left", fontweight="bold", fontsize=11)


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
    summary = pd.read_csv(VD / "Fig1c_PGP_Specificity_Summary.csv")
    results = pd.read_csv(VD / "Fig1c_PGP_Specificity_Results.csv")
    # Box + jittered points per group
    groups = ["Terrestrial_Soil", "Space_Flight"]
    data = [results.loc[results["Study_Group"] == g, "Specificity_H"].values
            for g in groups]
    colors = ["#2b8cbe", "#e34a33"]
    bp = ax.boxplot(data, patch_artist=True, widths=0.55)
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c); patch.set_alpha(0.5)
    for line in bp["medians"]:
        line.set_color("black")
    # jitter points
    import numpy as np
    rng = np.random.default_rng(0)
    for i, (g, c) in enumerate(zip(groups, colors)):
        vals = data[i]
        x = rng.normal(i + 1, 0.06, size=len(vals))
        ax.scatter(x, vals, s=14, color=c, edgecolor="black", linewidth=0.3, alpha=0.7)
    ax.set_xticks([1, 2])
    ax.set_xticklabels(["Terrestrial", "Space"], fontsize=9)
    ax.set_ylabel("PGP Shannon entropy (H)", fontsize=9)
    # Annotate means from summary
    if "H_Mean" in summary.columns:
        means = summary.set_index("Group")["H_Mean"]
        for i, g in enumerate(groups):
            if g in means.index:
                ax.text(i + 1, ax.get_ylim()[1] * 0.95,
                        f"{means[g]:.2f}", ha="center", fontsize=8)
    # MWU p from summary
    p_text = ""
    if "MWU_p" in summary.columns:
        p_vals = summary["MWU_p"].dropna()
        if len(p_vals) > 0:
            p_text = f"MWU p = {p_vals.iloc[0]:.1e}"
    if p_text:
        ax.text(0.5, 1.02, p_text, transform=ax.transAxes,
                fontsize=8, color="dimgrey", ha="center")
    ax.set_title("c) PGP specificity entropy",
                 loc="left", fontweight="bold", fontsize=11)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(17, 6.5))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.3, 1.6, 0.7], wspace=0.4)
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
