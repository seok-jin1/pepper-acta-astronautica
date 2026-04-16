"""
Fig 6 — Assembly mechanism (βNTI/RCbray) + BLSS failure-mode synthesis.

(a) βNTI vs RCbray quadrant scatter with process classification.
(b) Failure-mode cascade: environment → redundancy loss → N-cycle failure
    → BLSS risk → mitigations.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import pandas as pd
import numpy as np
import config

OUT = config.AA_FIG_DIR


# ---------------------------------------------------------------------------
# Panel a — βNTI/RCbray quadrant
# ---------------------------------------------------------------------------
def panel_a_bnti_rcbray(ax):
    df = pd.read_csv(config.VERSION_DIR / "SuppS4_BNTI_RCbray_results.csv")
    # Focus on within-group comparisons
    within = df[df["group1"] == df["group2"]].copy()
    colors = {"Space_Flight": "#e34a33", "Terrestrial_Soil": "#2b8cbe"}
    for grp, sub in within.groupby("group1"):
        ax.scatter(sub["bnti"], sub["rcbray"], alpha=0.35, s=22,
                   color=colors.get(grp, "grey"),
                   label=f"{grp.replace('_', ' ')} (n={len(sub)})")
    # Quadrant lines
    ax.axvline(-2, color="grey", ls="--", lw=0.7)
    ax.axvline(2, color="grey", ls="--", lw=0.7)
    ax.axhline(0.95, color="grey", ls=":", lw=0.7)
    ax.axhline(-0.95, color="grey", ls=":", lw=0.7)
    ax.set_xlabel("βNTI")
    ax.set_ylabel("RCbray")
    ax.legend(fontsize=9, loc="upper right")
    ax.text(-5, 0.9, "Homogeneous\nselection", fontsize=7, ha="left",
            va="top", color="dimgrey")
    ax.text(5, 0.9, "Variable\nselection", fontsize=7, ha="right",
            va="top", color="dimgrey")
    ax.text(0, -1, "Homogenizing\ndispersal", fontsize=7, ha="center",
            va="top", color="dimgrey")
    ax.set_title("a) Community assembly: Terrestrial is deterministic, "
                 "Space becomes stochastic",
                 loc="left", fontweight="bold", fontsize=10)


# ---------------------------------------------------------------------------
# Panel b — Failure-mode cascade
# ---------------------------------------------------------------------------
def panel_b_cascade(ax):
    ax.set_xlim(0, 10); ax.set_ylim(0, 10)
    ax.axis("off")

    boxes = [
        # (x, y, w, h, text, color)
        (0.2, 7.8, 2.8, 1.6,
         "Spaceflight environmental filter\n(microgravity, closed system,\nlimited substrate microbiome)",
         "#fcbba1"),
        (3.6, 7.8, 2.8, 1.6,
         "Preferential loss of\nnitrifier genera\n(Nitrosomonas, Nitrosospira,\nNitrospira, Nitrobacter)",
         "#fc9272"),
        (7.0, 7.8, 2.8, 1.6,
         "GTDB-verified\nstructural non-redundancy:\namoA in only 2/13 genera",
         "#fb6a4a"),
        (3.6, 5.2, 2.8, 1.6,
         "N-cycle completeness\ncollapse 100 % → 73 %\n(Log₂FC nitrif = −10.68)",
         "#de2d26"),
        (7.0, 5.2, 2.8, 1.6,
         "Network keystone\nRhodanobacter hub lost\n(degree 8 → 0)",
         "#a50f15"),
        (3.6, 2.6, 2.8, 1.6,
         "BLSS single-point-of-failure:\nno exogenous nitrifier source\nin ISS interior (OSD-65 null)",
         "#67000d"),
        (0.2, 0.0, 4.6, 2.0,
         "Mitigation 1: Exogenous nitrifier inoculation\n"
         "(load-bearing, supported by OSD-65 data)",
         "#c7e9c0"),
        (5.2, 0.0, 4.6, 2.0,
         "Mitigation 2: Per-crop risk profiling\n+ Real-time amoA biosensor monitoring\n"
         "(supported by Fig 5 crop divergence)",
         "#a1d99b"),
    ]
    for x, y, w, h, t, c in boxes:
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.06",
                                    facecolor=c, edgecolor="#333"))
        ax.text(x + w/2, y + h/2, t, ha="center", va="center", fontsize=8.5,
                color="black" if c != "#67000d" else "white",
                fontweight="bold" if c in ("#67000d", "#a50f15") else "normal")

    arrows = [
        # (xa, ya, xb, yb)
        (3.0, 8.6, 3.6, 8.6),
        (6.4, 8.6, 7.0, 8.6),
        (5.0, 7.8, 5.0, 6.8),
        (8.4, 7.8, 8.0, 6.8),
        (5.0, 5.2, 5.0, 4.2),
        (8.0, 5.2, 6.4, 3.4),
        (4.0, 2.6, 2.5, 2.0),
        (6.0, 2.6, 7.5, 2.0),
    ]
    for xa, ya, xb, yb in arrows:
        ax.add_patch(FancyArrowPatch((xa, ya), (xb, yb),
                                     arrowstyle="->",
                                     mutation_scale=15, color="#333",
                                     linewidth=1.2))
    ax.set_title("b) Failure-mode cascade and BLSS mitigations",
                 loc="left", fontweight="bold", fontsize=10)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(16, 7),
                             gridspec_kw={"width_ratios": [1, 1.3]})
    panel_a_bnti_rcbray(axes[0])
    panel_b_cascade(axes[1])
    plt.tight_layout()
    fig.suptitle("Figure 6. Community assembly mechanism and BLSS failure-mode synthesis",
                 fontsize=13, fontweight="bold", y=1.01)
    out = OUT / "Fig6_mechanism_synthesis"
    fig.savefig(str(out) + ".png", dpi=300, bbox_inches="tight")
    fig.savefig(str(out) + ".pdf", bbox_inches="tight")
    print(f"Wrote {out}.png/.pdf")


if __name__ == "__main__":
    main()
