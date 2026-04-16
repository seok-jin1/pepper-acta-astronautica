"""
Fig 2 — Failure signature: rhizosphere functional collapse.
5 panels composing existing pepper results for Acta Astronautica.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import pandas as pd
import numpy as np
import config

OUT = config.AA_FIG_DIR
VD = config.VERSION_DIR


def panel_a_faith_pd(ax):
    """α-diversity Faith's PD, Space vs Terrestrial."""
    img = VD / "SuppS1a_FaithPD.png"
    ax.imshow(mpimg.imread(img))
    ax.set_axis_off()
    ax.set_title("a) α-diversity (Faith's PD)", loc="left", fontweight="bold")


def panel_b_permanova(ax):
    img = VD / "SuppS1c_PERMANOVA.png"
    ax.imshow(mpimg.imread(img))
    ax.set_axis_off()
    ax.set_title("b) β-diversity PERMANOVA (bootstrap)", loc="left", fontweight="bold")


def panel_c_faprotax(ax):
    img = VD / "Fig2A_FAPROTAX_Log2FC.png"
    ax.imshow(mpimg.imread(img))
    ax.set_axis_off()
    ax.set_title("c) FAPROTAX Log₂FC (Space / Terr)", loc="left", fontweight="bold")


def panel_d_consensus(ax):
    img = OUT / "Fig2d_consensus.png"
    ax.imshow(mpimg.imread(img))
    ax.set_axis_off()
    ax.set_title("d) FAPROTAX × PICRUSt2 consensus", loc="left", fontweight="bold")


def panel_e_ncycle_heatmap(ax):
    img = config.CURRENT_BIOLOGY_DIR / "Fig1c_NCycle_Heatmap.png"
    ax.imshow(mpimg.imread(img))
    ax.set_axis_off()
    ax.set_title("e) N-cycle completeness (Terr 100 % → Space 76 %)", loc="left", fontweight="bold")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(18, 13))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.25],
                          hspace=0.22, wspace=0.12)
    panel_a_faith_pd(fig.add_subplot(gs[0, 0]))
    panel_b_permanova(fig.add_subplot(gs[0, 1]))
    panel_c_faprotax(fig.add_subplot(gs[0, 2]))
    panel_d_consensus(fig.add_subplot(gs[1, 0:2]))
    panel_e_ncycle_heatmap(fig.add_subplot(gs[1, 2]))
    fig.suptitle("Figure 2. Rhizosphere functional collapse in ISS-grown Capsicum annuum",
                 fontsize=13, fontweight="bold", y=0.995)
    out = OUT / "Fig2_failure_signature"
    fig.savefig(str(out) + ".png", dpi=300, bbox_inches="tight")
    fig.savefig(str(out) + ".pdf", bbox_inches="tight")
    print(f"Wrote {out}.png/.pdf")


if __name__ == "__main__":
    main()
