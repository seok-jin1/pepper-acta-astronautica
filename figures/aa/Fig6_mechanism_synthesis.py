"""
Fig 6 — Evidence synthesis: structural bottleneck and BLSS design considerations.

Evidence cascade: environmental filter → nitrifier loss → structural
bottleneck → design considerations.

The βNTI/RCbray assembly proxy (previously panel a) has been moved to
Supplementary Fig. S4 because Methods §2.5 and Results §3.5 frame the
assembly analysis as exploratory supplementary context (99 null-model
randomisations, sample-level metric averaged to pairwise).
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


def panel_cascade(ax):
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
         "RefSeq/PGAP-verified\nstructural non-redundancy:\namoA in only 2/13 genera",
         "#fb6a4a"),
        (3.6, 5.2, 2.8, 1.6,
         "N-cycle completeness\ndepletion 100 % → 76 %\n(Log₂FC nitrif = −13.08)",
         "#de2d26"),
        (3.6, 2.6, 2.8, 1.6,
         "No nitrifier genera detected\nin ISS interior surface pool\n(OSD-65, n=196)",
         "#67000d"),
        (0.2, 0.0, 4.6, 2.0,
         "Design consideration 1:\nExogenous nitrifier inoculation\n(consistent with OSD-65 null)",
         "#c7e9c0"),
        (5.2, 0.0, 4.6, 2.0,
         "Design consideration 2:\nPer-crop characterisation +\namoA-targeted monitoring\n(consistent with Fig 5 crop divergence)",
         "#a1d99b"),
    ]
    for x, y, w, h, t, c in boxes:
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.06",
                                    facecolor=c, edgecolor="#333"))
        ax.text(x + w/2, y + h/2, t, ha="center", va="center", fontsize=8.5,
                color="black" if c != "#67000d" else "white",
                fontweight="bold" if c == "#67000d" else "normal")

    arrows = [
        # (xa, ya, xb, yb)
        (3.0, 8.6, 3.6, 8.6),
        (6.4, 8.6, 7.0, 8.6),
        (5.0, 7.8, 5.0, 6.8),
        (8.4, 7.8, 6.4, 4.2),
        (5.0, 5.2, 5.0, 4.2),
        (4.0, 2.6, 2.5, 2.0),
        (6.0, 2.6, 7.5, 2.0),
    ]
    for xa, ya, xb, yb in arrows:
        ax.add_patch(FancyArrowPatch((xa, ya), (xb, yb),
                                     arrowstyle="->",
                                     mutation_scale=15, color="#333",
                                     linewidth=1.2))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(1, 1, figsize=(11, 8))
    panel_cascade(ax)
    plt.tight_layout()
    out = OUT / "Fig6_mechanism_synthesis"
    fig.savefig(str(out) + ".png", dpi=300, bbox_inches="tight")
    fig.savefig(str(out) + ".pdf", bbox_inches="tight")
    print(f"Wrote {out}.png/.pdf")


if __name__ == "__main__":
    main()
