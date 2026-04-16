"""
Fig 4 — System-level consequence: network collapse and nitrifier-free ISS source pool.

(a) Network metrics Space vs Terrestrial (from Fig1e_Network_Genus_Metrics.csv)
(b) Threshold sensitivity — Rhodanobacter hub degree across |r| thresholds
(c) ISS env overall overlap (genera shared vs unique)
(d) N-cycle genus presence/absence in ISS env vs pepper Space (KEY ARGUMENT)
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import config

OUT = config.AA_FIG_DIR


def panel_a_network_metrics(ax):
    df = pd.read_csv(config.VERSION_DIR / "Fig2E_Network_Genus_Metrics.csv")
    metrics = ["Nodes", "Edges", "Density", "Avg_Degree", "Negative_Edges"]
    df_plot = df.set_index("Group")[metrics].T
    x = np.arange(len(metrics))
    w = 0.38
    ax.bar(x - w/2, df_plot["Space Flight"].values, w,
           label="Space Flight", color="#e34a33")
    ax.bar(x + w/2, df_plot["Terrestrial Soil"].values, w,
           label="Terrestrial Soil", color="#2b8cbe")
    ax.set_xticks(x); ax.set_xticklabels(metrics, rotation=25, ha="right")
    ax.legend(fontsize=8, loc="upper right")
    ax.set_title("a) Co-occurrence network metrics (|r|>0.4)",
                 loc="left", fontweight="bold")


def panel_b_threshold_sensitivity(ax):
    df = pd.read_csv(config.VERSION_DIR / "SuppS2b_ThresholdSensitivity_Results.csv")
    groups = df["Group"].unique()
    colors = {"Space_Flight": "#e34a33", "Terrestrial_Soil": "#2b8cbe"}
    for g in groups:
        sub = df[df["Group"] == g].sort_values("r_threshold")
        ax.plot(sub["r_threshold"], sub["Rhodanobacter_Degree"],
                marker="o", label=g.replace("_", " "),
                color=colors.get(g, "grey"), linewidth=2)
    ax.set_xlabel("Spearman |r| threshold")
    ax.set_ylabel("Rhodanobacter degree")
    ax.legend(fontsize=8)
    ax.set_title("b) Rhodanobacter hub robustness across thresholds",
                 loc="left", fontweight="bold")


def panel_c_overall_overlap(ax):
    s = pd.read_csv(config.AA_ISS_ENV_SOURCES_CSV, index_col=0).iloc[:, 0]
    # Visualize 3 categories
    labels = ["Pepper only", "Shared\n(ISS env + pepper)", "ISS env only"]
    values = [int(s["n_pepper_only"]), int(s["n_overlap"]), int(s["n_iss_only"])]
    colors = ["#e34a33", "#9ecae1", "#2b8cbe"]
    ax.pie(values, labels=labels, colors=colors, autopct="%1.0f",
           startangle=90, textprops={"fontsize": 9})
    ax.set_title(f"c) Genus pool overlap\n({values[1]}/{values[0]+values[1]} = "
                 f"{s['pct_pepper_recoverable_from_iss_env']:.0f}% of pepper Space recoverable)",
                 loc="left", fontweight="bold", fontsize=10)


def panel_d_ncycle_iss_env(ax):
    df = pd.read_csv(config.AA_DIR / "iss_env_ncycle_genera.csv", index_col=0)
    # Categorize rows for color: nitrifiers vs denitrifiers vs other
    nitrif = {"Nitrosomonas", "Nitrosospira", "Nitrospira", "Nitrobacter"}
    denitrif = {"Rhodanobacter", "Pseudomonas", "Burkholderia", "Azospirillum"}
    categories = []
    for g in df.index:
        if g in nitrif: categories.append("Nitrifier")
        elif g == "Rhodanobacter": categories.append("Keystone denitrifier")
        elif g in denitrif: categories.append("Denitrifier")
        else: categories.append("Generalist")
    df["category"] = categories
    df = df.sort_values(["category", "pepper_space_count"], ascending=[True, False])

    # Build 2-column binary matrix
    mat = df[["in_iss_env", "in_pepper_space"]].astype(int)
    import matplotlib.colors as mcolors
    colmap = {"Nitrifier": "#d94801",
              "Keystone denitrifier": "#993404",
              "Denitrifier": "#2171b5",
              "Generalist": "#bdbdbd"}

    y = np.arange(len(mat))
    for i, (idx, row) in enumerate(mat.iterrows()):
        cat = df.loc[idx, "category"]
        ax.scatter(0, i, s=260, marker="s",
                   color="white" if row["in_iss_env"] == 0 else colmap[cat],
                   edgecolor=colmap[cat], linewidth=1.6)
        ax.scatter(1, i, s=260, marker="s",
                   color="white" if row["in_pepper_space"] == 0 else colmap[cat],
                   edgecolor=colmap[cat], linewidth=1.6)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["ISS env\n(OSD-65)", "Pepper Space\nrhizosphere"])
    ax.set_yticks(y); ax.set_yticklabels(mat.index, fontstyle="italic", fontsize=9)
    ax.set_xlim(-0.5, 1.5); ax.set_ylim(-0.5, len(mat) - 0.5)
    ax.invert_yaxis()
    # Legend
    from matplotlib.patches import Patch
    handles = [Patch(facecolor=c, edgecolor=c, label=k) for k, c in colmap.items()]
    ax.legend(handles=handles, loc="center right", bbox_to_anchor=(1.55, 0.5),
              fontsize=8, frameon=False)
    ax.set_title("d) N-cycle genera: all nitrifiers absent from ISS env",
                 loc="left", fontweight="bold", fontsize=10)
    ax.grid(axis="y", color="lightgrey", alpha=0.4, linestyle=":")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(14, 10),
                             gridspec_kw={"height_ratios": [1, 1.25],
                                          "width_ratios": [1, 1.2]})
    panel_a_network_metrics(axes[0, 0])
    panel_b_threshold_sensitivity(axes[0, 1])
    panel_c_overall_overlap(axes[1, 0])
    panel_d_ncycle_iss_env(axes[1, 1])
    plt.tight_layout()
    fig.suptitle("Figure 4. Network collapse and nitrifier-free ISS source pool",
                 fontsize=13, fontweight="bold", y=1.01)
    out = OUT / "Fig4_system_consequence"
    fig.savefig(str(out) + ".png", dpi=300, bbox_inches="tight")
    fig.savefig(str(out) + ".pdf", bbox_inches="tight")
    print(f"Wrote {out}.png/.pdf")


if __name__ == "__main__":
    main()
