"""
Fig 4 — ISS interior surface pool does not contain nitrifier genera.

(a) ISS env overall overlap (genera shared vs unique)
(b) N-cycle genus presence/absence in ISS env vs pepper Space (KEY ARGUMENT)

Network panels (Rhodanobacter hub, threshold sensitivity) moved to
Supplementary Fig. S5 per manuscript positioning (see SuppFig_S5_network.py).
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import config

OUT = config.AA_FIG_DIR


def panel_a_overall_overlap(ax):
    s = pd.read_csv(config.AA_ISS_ENV_SOURCES_CSV, index_col=0).iloc[:, 0]
    labels = ["Pepper only", "Shared\n(ISS env + pepper)", "ISS env only"]
    values = [int(s["n_pepper_only"]), int(s["n_overlap"]), int(s["n_iss_only"])]
    colors = ["#e34a33", "#9ecae1", "#2b8cbe"]
    ax.pie(values, labels=labels, colors=colors, autopct="%1.0f",
           startangle=90, textprops={"fontsize": 9})
    ax.set_title(f"a) Genus pool overlap\n({values[1]}/{values[0]+values[1]} = "
                 f"{s['pct_pepper_recoverable_from_iss_env']:.0f}% of pepper Space recoverable)",
                 loc="left", fontweight="bold", fontsize=10)


def panel_b_ncycle_iss_env(ax):
    df = pd.read_csv(config.AA_DIR / "iss_env_ncycle_genera.csv", index_col=0)
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

    mat = df[["in_iss_env", "in_pepper_space"]].astype(int)
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
    from matplotlib.patches import Patch
    handles = [Patch(facecolor=c, edgecolor=c, label=k) for k, c in colmap.items()]
    ax.legend(handles=handles, loc="center right", bbox_to_anchor=(1.55, 0.5),
              fontsize=8, frameon=False)
    ax.set_title("b) N-cycle genera: all nitrifiers absent from ISS env",
                 loc="left", fontweight="bold", fontsize=10)
    ax.grid(axis="y", color="lightgrey", alpha=0.4, linestyle=":")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(13, 7),
                             gridspec_kw={"width_ratios": [1, 1.4]})
    panel_a_overall_overlap(axes[0])
    panel_b_ncycle_iss_env(axes[1])
    plt.tight_layout()
    out = OUT / "Fig4_system_consequence"
    fig.savefig(str(out) + ".png", dpi=300, bbox_inches="tight")
    fig.savefig(str(out) + ".pdf", bbox_inches="tight")
    print(f"Wrote {out}.png/.pdf")


if __name__ == "__main__":
    main()
