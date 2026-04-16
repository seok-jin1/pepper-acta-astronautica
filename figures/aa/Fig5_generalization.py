"""
Fig 5 — Crop-specific failure mode: divergent lettuce vs pepper signatures.
Reframed from 'replication' to evidence for crop-specific BLSS risk models.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import config

OUT = config.AA_FIG_DIR
KEY_FUNCS = ["nitrification", "denitrification", "nitrogen_fixation", "chemoheterotrophy"]


def panel_a_log2fc_sidebyside(ax):
    df = pd.read_csv(config.AA_VEGGIE_RESULTS_CSV)
    # pepper entries use FAPROTAX group names; lettuce uses the same.
    # Keep only key functions and latest-crop-per-function.
    df = df[df["function"].isin(KEY_FUNCS)]
    pivot = df.pivot_table(index="function", columns="crop", values="log2fc",
                           aggfunc="mean").reindex(KEY_FUNCS)

    x = np.arange(len(pivot.index)); w = 0.38
    ax.bar(x - w/2, pivot["pepper"].values, w, label="Pepper (OSD-772)",
           color="#d94801")
    ax.bar(x + w/2, pivot["lettuce"].values, w, label="Lettuce (OSD-267+269)",
           color="#2171b5")
    ax.axhline(0, color="black", lw=0.6)
    ax.set_xticks(x); ax.set_xticklabels(pivot.index, rotation=20, ha="right")
    ax.set_ylabel("FAPROTAX Log₂FC (Space / Terrestrial)")
    ax.legend(fontsize=9, loc="upper left")
    ax.set_title("a) Divergent crop responses: pepper nitrification collapses, "
                 "lettuce enriches", loc="left", fontweight="bold", fontsize=10)
    # Annotate values
    for i, fn in enumerate(pivot.index):
        for j, col in enumerate(["pepper", "lettuce"]):
            val = pivot.loc[fn, col]
            off = -0.19 + j * 0.38
            y_text = val + (0.4 if val >= 0 else -0.8)
            ax.text(i + off, y_text, f"{val:+.1f}", ha="center",
                    va="bottom" if val >= 0 else "top", fontsize=8)


def panel_b_ncycle_completeness(ax):
    pepper_summary = pd.read_csv(config.VERSION_DIR / "Fig1c_NCycle_Completeness_Summary.csv")
    # Pepper summary columns: Group, Completeness_Mean (fractional), ...
    # Convert fractional to percentage for parity with lettuce
    p_mean = (pepper_summary.set_index("Group")["Completeness_Mean"] * 100)
    lettuce_nc = pd.read_csv(config.AA_DIR / "veggie_ncycle_completeness.csv")
    l_mean = lettuce_nc.groupby("Study_Group")["completeness"].mean()

    cats = ["Terrestrial_Soil", "Space_Flight"]
    x = np.arange(len(cats)); w = 0.38
    ax.bar(x - w/2, [p_mean.get(c, 0) for c in cats], w,
           label="Pepper", color="#d94801")
    ax.bar(x + w/2, [l_mean.get(c, 0) for c in cats], w,
           label="Lettuce", color="#2171b5")
    ax.set_xticks(x); ax.set_xticklabels([c.replace("_", " ") for c in cats])
    ax.set_ylabel("N-cycle completeness (%)")
    ax.set_ylim(0, 105)
    for i, c in enumerate(cats):
        for j, (mean, label) in enumerate([(p_mean, "P"), (l_mean, "L")]):
            v = mean.get(c, 0)
            ax.text(i - 0.19 + j * 0.38, v + 1.5, f"{v:.0f}%",
                    ha="center", fontsize=8)
    ax.legend(fontsize=9)
    ax.set_title("b) N-cycle completeness diverges between crops",
                 loc="left", fontweight="bold", fontsize=10)


def panel_c_redundancy_sidebyside(ax):
    pepper_r = pd.read_csv(config.VERSION_DIR / "Fig1d_Functional_Redundancy_Summary.csv")
    lettuce_r = pd.read_csv(config.AA_DIR / "veggie_functional_redundancy.csv")

    def get_pepper_value(func, grp):
        # Try common column shapes
        if "function" in pepper_r.columns and "Study_Group" in pepper_r.columns:
            m = pepper_r[(pepper_r["function"] == func) &
                         (pepper_r["Study_Group"] == grp)]
            for col in ["n_contributing_genera", "Mean_N_Genera", "N_Genera"]:
                if col in m.columns and len(m):
                    return float(m[col].iloc[0])
        return 0.0

    def get_lettuce_value(func, grp):
        m = lettuce_r[(lettuce_r["function"] == func) &
                      (lettuce_r["Study_Group"] == grp)]
        return float(m["n_contributing_genera"].iloc[0]) if len(m) else 0.0

    funcs = ["nitrification", "denitrification", "nitrogen_fixation"]
    x = np.arange(len(funcs))

    vals_pepper_terr = [get_pepper_value(f, "Terrestrial_Soil") for f in funcs]
    vals_pepper_sp = [get_pepper_value(f, "Space_Flight") for f in funcs]
    vals_lettuce_terr = [get_lettuce_value(f, "Terrestrial_Soil") for f in funcs]
    vals_lettuce_sp = [get_lettuce_value(f, "Space_Flight") for f in funcs]

    w = 0.2
    ax.bar(x - 1.5*w, vals_pepper_terr, w, label="Pepper Terr",
           color="#fdd0a2", edgecolor="#d94801")
    ax.bar(x - 0.5*w, vals_pepper_sp, w, label="Pepper Space",
           color="#d94801")
    ax.bar(x + 0.5*w, vals_lettuce_terr, w, label="Lettuce Terr",
           color="#c6dbef", edgecolor="#2171b5")
    ax.bar(x + 1.5*w, vals_lettuce_sp, w, label="Lettuce Space",
           color="#2171b5")
    ax.set_xticks(x); ax.set_xticklabels(funcs, rotation=15, ha="right")
    ax.set_ylabel("# contributing genera")
    ax.legend(fontsize=7, ncol=2)
    ax.set_title("c) Redundancy: pepper nitrifiers collapse (4.8→0.2), "
                 "lettuce pattern differs",
                 loc="left", fontweight="bold", fontsize=10)


def panel_d_message(ax):
    """Summary message box — Crop-specific failure mode implication."""
    ax.axis("off")
    msg = (
        "Crop-specific failure modes:\n\n"
        "• Pepper (OSD-772) shows catastrophic nitrification\n"
        "  collapse (Log₂FC = −10.68, 4.8 → 0.2 genera).\n\n"
        "• Lettuce (OSD-267+269) shows the opposite:\n"
        "  nitrification Log₂FC = +11.45 in Space.\n\n"
        "• Conclusion: BLSS biological risk models built\n"
        "  on a single crop cannot be generalized.\n"
        "  Per-crop microbiome profiling is required.\n"
    )
    ax.text(0.02, 0.98, msg, transform=ax.transAxes, fontsize=10,
            va="top", family="sans-serif",
            bbox=dict(boxstyle="round,pad=0.6",
                      facecolor="#fff5eb", edgecolor="#d94801"))
    ax.set_title("d) Engineering implication", loc="left",
                 fontweight="bold", fontsize=10)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    panel_a_log2fc_sidebyside(axes[0, 0])
    panel_b_ncycle_completeness(axes[0, 1])
    panel_c_redundancy_sidebyside(axes[1, 0])
    panel_d_message(axes[1, 1])
    plt.tight_layout()
    fig.suptitle("Figure 5. Crop-specific failure mode — pepper and lettuce diverge",
                 fontsize=13, fontweight="bold", y=1.01)
    out = OUT / "Fig5_crop_specific"
    fig.savefig(str(out) + ".png", dpi=300, bbox_inches="tight")
    fig.savefig(str(out) + ".pdf", bbox_inches="tight")
    print(f"Wrote {out}.png/.pdf")


if __name__ == "__main__":
    main()
