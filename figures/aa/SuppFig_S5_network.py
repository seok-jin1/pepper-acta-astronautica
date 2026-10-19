"""
Supplementary Fig. S5 — Exploratory co-occurrence network context.

(a) Network metrics Space vs Terrestrial (from Fig2E_Network_Genus_Metrics.csv)
(b) Threshold sensitivity — Rhodanobacter hub degree across |r| thresholds
    (with BH-FDR correction; at n=20 the FDR cutoff is stricter than the
    magnitude thresholds tested, so all three thresholds yield identical
    edge sets)

Moved from main Fig 4 (a, b) because the network analysis is explicitly
framed as supplementary exploratory context (Methods §2.5, Results §3.5).
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
    ax.set_title("a) Co-occurrence network metrics (|r|>0.4, BH q<0.05)",
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


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    panel_a_network_metrics(axes[0])
    panel_b_threshold_sensitivity(axes[1])
    plt.tight_layout()
    out = OUT / "SuppFig_S5_network"
    fig.savefig(str(out) + ".png", dpi=300, bbox_inches="tight")
    fig.savefig(str(out) + ".pdf", bbox_inches="tight")
    print(f"Wrote {out}.png/.pdf")


if __name__ == "__main__":
    main()
