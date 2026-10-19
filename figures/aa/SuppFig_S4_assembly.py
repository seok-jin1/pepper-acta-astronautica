"""
Supplementary Fig. S4 — Exploratory community-assembly proxy.

βNTI vs RCbray quadrant scatter for Space vs Terrestrial within-group pairs.

Methodological caveats (see Methods §2.5):
  - ses.mntd() is a sample-level metric; pairwise values derived by
    averaging the two sample-level values. This is not a true pairwise
    null model and introduces bias relative to a proper pairwise βNTI
    implementation.
  - 99 null-model randomisations (not 999) due to runtime with
    11,345 ASVs × 40 samples. 99 randomisations are insufficient to
    reliably resolve individual sample null distributions at |z| > 2.

Results presented as exploratory directional context only, not as
quantitatively precise process-fraction estimates.

Moved from main Fig 6 panel a per manuscript positioning.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import matplotlib.pyplot as plt
import pandas as pd
import config

OUT = config.AA_FIG_DIR


def panel_bnti_rcbray(ax):
    df = pd.read_csv(config.VERSION_DIR / "SuppS4_BNTI_RCbray_results.csv")
    within = df[df["group1"] == df["group2"]].copy()
    colors = {"Space_Flight": "#e34a33", "Terrestrial_Soil": "#2b8cbe"}
    for grp, sub in within.groupby("group1"):
        ax.scatter(sub["bnti"], sub["rcbray"], alpha=0.35, s=22,
                   color=colors.get(grp, "grey"),
                   label=f"{grp.replace('_', ' ')} (n={len(sub)})")
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
    ax.set_title("Exploratory assembly proxy (sample-level ses.mntd, "
                 "99 null randomisations; see Methods)",
                 loc="left", fontweight="bold", fontsize=10)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(1, 1, figsize=(7, 6))
    panel_bnti_rcbray(ax)
    plt.tight_layout()
    out = OUT / "SuppFig_S4_assembly"
    fig.savefig(str(out) + ".png", dpi=300, bbox_inches="tight")
    fig.savefig(str(out) + ".pdf", bbox_inches="tight")
    print(f"Wrote {out}.png/.pdf")


if __name__ == "__main__":
    main()
