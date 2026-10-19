"""
48_rbd_simulation.py
====================
Reliability Block Diagram (RBD) simulation of the ISS rhizosphere N-cycle,
Acta Astronautica engineering framing — MELiSSA Compartment 2 analogue.

Each N-cycle biochemical step is a BLOCK. Within a block, the multiple taxa
that carry the relevant gene are a PARALLEL redundancy. Across blocks, the
conversion chain is SERIES (all steps must work to deliver plant-available
NO3- from Osmocote NH4+).

Parameters come from the GTDB/PGAP matrix (script 36):
  - amoA:        2/13 canonical rhizosphere genera (Nitrosomonas, Nitrosospira)
  - nxrA-NOB:    2 canonical NOB lineages (Nitrobacter + Nitrospira
                 literature-confirmed despite generic PGAP annotation)
  - nirK (denit loss): 10/13 (parallel, but here a "leak", not a supply step)

Reliability math (standard RBD):
  Within a block with n_parallel taxa, each with independent loss probability p
  after g generations:
      P(single taxon survives g gen)     = (1 - p)^g
      P(block fails at g)                = (1 - (1-p)^g)^n_parallel
      Block reliability                  = 1 - P(block fails)
  System is SERIES:
      R_system(g) = ∏ R_block(g)
  And FAILURE timeline:
      g_fail(x) = min g such that 1 - R_system(g) ≥ x

Outputs:
  - acta_astronautica/simulation/rbd_blocks.csv
  - acta_astronautica/simulation/rbd_reliability_grid.csv
  - acta_astronautica/simulation/rbd_diagram.pdf/.png
  - acta_astronautica/simulation/rbd_reliability_curves.pdf/.png
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config

OUT_DIR = config.AA_DIR / "simulation"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# -------------------------------------------------------------------
# Block definitions (from GTDB/PGAP + canonical literature)
# -------------------------------------------------------------------
BLOCKS = [
    {
        "id": "OSMOCOTE",
        "label": "NH\u2084\u207A supply\n(Osmocote CRF)",
        "role": "source",
        "n_parallel": 1,       # chemistry release, no taxonomic redundancy
        "genes": "—",
    },
    {
        "id": "AMO",
        "label": "Ammonia oxidation\n(amoA, NH\u2084\u207A \u2192 NO\u2082\u207B)",
        "role": "series",
        "n_parallel": 2,       # Nitrosomonas, Nitrosospira
        "genes": "amoA",
        "taxa": ["Nitrosomonas", "Nitrosospira"],
    },
    {
        "id": "NXR",
        "label": "Nitrite oxidation\n(nxrA, NO\u2082\u207B \u2192 NO\u2083\u207B)",
        "role": "series",
        "n_parallel": 2,       # Nitrobacter + Nitrospira (literature)
        "genes": "nxrA",
        "taxa": ["Nitrobacter", "Nitrospira"],
    },
    {
        "id": "PLANT",
        "label": "Plant NO\u2083\u207B uptake",
        "role": "sink",
        "n_parallel": 1,
        "genes": "—",
    },
]

# Leak: denitrification removes NO3- (BLSS undesirable)
LEAK_BLOCK = {
    "id": "DENIT",
    "label": "Denitrification\n(nirK leak)",
    "role": "leak",
    "n_parallel": 10,          # 10/13 carry nirK
    "genes": "nirK",
}


def block_reliability(n_parallel: int, p_loss: float, g: int) -> float:
    """Probability a block still has ≥1 functional taxon after g generations."""
    surv = (1.0 - p_loss) ** g
    block_fails = (1.0 - surv) ** n_parallel
    return 1.0 - block_fails


def system_reliability(p_loss: float, g: int) -> float:
    """Series: NH4 supply * AMO * NXR * Plant (source/sink treated as 1.0)."""
    r = 1.0
    for blk in BLOCKS:
        if blk["role"] in ("source", "sink"):
            continue
        r *= block_reliability(blk["n_parallel"], p_loss, g)
    return r


def generations_to_failure(p_loss: float, threshold: float = 0.5,
                           max_gen: int = 200) -> int:
    """Smallest g such that cumulative failure probability ≥ threshold."""
    for g in range(1, max_gen + 1):
        if 1.0 - system_reliability(p_loss, g) >= threshold:
            return g
    return -1  # not reached within horizon


# -------------------------------------------------------------------
# 1. Emit block-definition table
# -------------------------------------------------------------------

def write_block_table():
    rows = []
    for blk in BLOCKS:
        rows.append({
            "block_id": blk["id"],
            "label": blk["label"].replace("\n", " "),
            "role": blk["role"],
            "n_parallel_taxa": blk["n_parallel"],
            "genes": blk["genes"],
            "taxa": ";".join(blk.get("taxa", [])) if "taxa" in blk else "",
        })
    rows.append({
        "block_id": LEAK_BLOCK["id"],
        "label": LEAK_BLOCK["label"].replace("\n", " "),
        "role": LEAK_BLOCK["role"],
        "n_parallel_taxa": LEAK_BLOCK["n_parallel"],
        "genes": LEAK_BLOCK["genes"],
        "taxa": "",
    })
    df = pd.DataFrame(rows)
    df.to_csv(OUT_DIR / "rbd_blocks.csv", index=False)
    print(df.to_string(index=False))
    return df


# -------------------------------------------------------------------
# 2. Parametric sweep: p_loss × generations → system reliability
# -------------------------------------------------------------------

def parametric_sweep():
    p_values = [0.02, 0.05, 0.10, 0.15, 0.20, 0.30]
    g_values = np.arange(0, 31)
    rows = []
    for p in p_values:
        for g in g_values:
            rows.append({
                "p_loss_per_generation": p,
                "generation": int(g),
                "R_AMO":   block_reliability(2, p, g),
                "R_NXR":   block_reliability(2, p, g),
                "R_DENIT": block_reliability(10, p, g),
                "R_system": system_reliability(p, g),
                "P_failure": 1.0 - system_reliability(p, g),
            })
    df = pd.DataFrame(rows)
    df.to_csv(OUT_DIR / "rbd_reliability_grid.csv", index=False)

    # Summary: p_loss → generation at 50% cumulative failure
    summary_rows = []
    for p in p_values:
        g50 = generations_to_failure(p, 0.5)
        g25 = generations_to_failure(p, 0.25)
        g75 = generations_to_failure(p, 0.75)
        summary_rows.append({
            "p_loss_per_generation": p,
            "g_to_25pct_failure": g25,
            "g_to_50pct_failure": g50,
            "g_to_75pct_failure": g75,
        })
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(OUT_DIR / "rbd_failure_thresholds.csv", index=False)
    print("\nFailure thresholds (generations until cumulative failure probability ≥ X):")
    print(summary.to_string(index=False))
    return df, summary


# -------------------------------------------------------------------
# 3. Block-diagram figure
# -------------------------------------------------------------------

def draw_block_diagram():
    fig, ax = plt.subplots(figsize=(13, 6))

    # Layout: four horizontal boxes for main chain + one leak below
    block_x = [0.05, 0.25, 0.50, 0.78]
    block_y = 0.55
    box_w = 0.17
    box_h = 0.28

    colors = {
        "source": "#b3cde3",
        "series": "#fbb4ae",
        "sink": "#ccebc5",
        "leak": "#decbe4",
    }

    # Draw main chain — source/sink with dashed outline to signal "boundary term"
    for i, blk in enumerate(BLOCKS):
        x = block_x[i]
        face = colors[blk["role"]]
        ls = "--" if blk["role"] in ("source", "sink") else "-"
        rect = mpatches.FancyBboxPatch((x, block_y), box_w, box_h,
                                        boxstyle="round,pad=0.01",
                                        linewidth=2, edgecolor="black",
                                        facecolor=face, linestyle=ls)
        ax.add_patch(rect)
        ax.text(x + box_w / 2, block_y + box_h * 0.75, blk["label"],
                ha="center", va="center", fontsize=10, fontweight="bold")
        if blk["role"] == "series":
            ax.text(x + box_w / 2, block_y + box_h * 0.3,
                    f"n = {blk['n_parallel']} parallel taxa\n"
                    f"({', '.join(blk.get('taxa', []))})",
                    ha="center", va="center", fontsize=8, fontstyle="italic")
        elif blk["role"] in ("source", "sink"):
            ax.text(x + box_w / 2, block_y + box_h * 0.3,
                    "[boundary term;\nexcluded from $R_\\mathrm{system}$]",
                    ha="center", va="center", fontsize=7.5, fontstyle="italic",
                    color="#555")

        # Arrow to next (with NO3- label on the arrow into plant uptake)
        if i < len(BLOCKS) - 1:
            ax.annotate("", xy=(block_x[i + 1], block_y + box_h / 2),
                        xytext=(x + box_w, block_y + box_h / 2),
                        arrowprops=dict(arrowstyle="->", lw=2.5, color="black"))
            # Label the NO3- flow between NXR (index 2) and Plant (index 3)
            if i == 2:
                mid_x = (x + box_w + block_x[i + 1]) / 2
                ax.text(mid_x, block_y + box_h / 2 + 0.04,
                        "NO$_3^-$ pool",
                        ha="center", va="bottom", fontsize=8, color="black")

    # Draw leak (denitrification) below — branch from the NXR → plant NO3- arrow,
    # NOT from the NXR block itself, because denitrification removes NO3- from
    # the NO3- pool generated by NXR, not directly from NXR's output internally.
    leak_x = 0.58   # sit under the NO3- pool arrow
    leak_y = 0.08
    leak_rect = mpatches.FancyBboxPatch((leak_x, leak_y), box_w, box_h,
                                         boxstyle="round,pad=0.01",
                                         linewidth=2, edgecolor="black",
                                         facecolor=colors["leak"], linestyle="--")
    ax.add_patch(leak_rect)
    ax.text(leak_x + box_w / 2, leak_y + box_h * 0.75, LEAK_BLOCK["label"],
            ha="center", va="center", fontsize=10, fontweight="bold")
    ax.text(leak_x + box_w / 2, leak_y + box_h * 0.3,
            f"n = {LEAK_BLOCK['n_parallel']} parallel taxa\n(architectural contrast)",
            ha="center", va="center", fontsize=8, fontstyle="italic")
    # Branch from the NO3- pool arrow midpoint (NXR → plant), going down to leak block
    nxr_end_x   = block_x[2] + box_w
    plant_start = block_x[3]
    branch_x    = (nxr_end_x + plant_start) / 2  # midpoint of NO3- arrow
    arrow_y     = block_y + box_h / 2
    ax.annotate("", xy=(leak_x + box_w / 2, leak_y + box_h),
                xytext=(branch_x, arrow_y),
                arrowprops=dict(arrowstyle="->", lw=2, color="gray",
                                linestyle="--"))
    ax.text(leak_x + box_w / 2 + 0.12, (arrow_y + leak_y + box_h) / 2,
            "NO$_3^-$ $\\to$ N$_2$ gas\n(BLSS-undesirable leak;\nnot part of $R_\\mathrm{system}$)",
            ha="left", va="center", fontsize=8, color="gray", fontstyle="italic")

    # Low-redundancy bottleneck annotation (replaces SPOF)
    for i in (1, 2):  # AMO, NXR
        x = block_x[i] + box_w / 2
        ax.annotate("low-redundancy bottleneck ($n=2$)",
                    xy=(x, block_y + box_h + 0.005),
                    xytext=(x, block_y + box_h + 0.05),
                    ha="center", fontsize=8, color="#c0392b", fontweight="bold",
                    arrowprops=dict(arrowstyle="->", color="#c0392b", lw=1.5))

    ax.set_xlim(0, 1)
    ax.set_ylim(-0.02, 1.05)
    ax.axis("off")
    ax.set_title("Reliability Block Diagram — ISS rhizosphere N-cycle\n"
                 "(series chain of biochemical blocks; parallel taxa within each block)",
                 fontsize=11, fontweight="bold", pad=20)

    legend = [
        mpatches.Patch(color=colors["source"], label="Boundary term (excluded from $R_\\mathrm{system}$)"),
        mpatches.Patch(color=colors["series"], label="Redundant biochemical step (low-redundancy bottleneck)"),
        mpatches.Patch(color=colors["sink"],   label="Boundary term / plant uptake"),
        mpatches.Patch(color=colors["leak"],   label="Architectural contrast: denitrification leak"),
    ]
    ax.legend(handles=legend, loc="lower left", fontsize=9, bbox_to_anchor=(0.02, -0.02))

    for ext in ("pdf", "png"):
        fig.savefig(OUT_DIR / f"rbd_diagram.{ext}",
                    dpi=300 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: rbd_diagram.pdf/.png")


# -------------------------------------------------------------------
# 4. Reliability curves figure
# -------------------------------------------------------------------

def draw_reliability_curves(df: pd.DataFrame):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)

    palette = plt.cm.viridis(np.linspace(0.15, 0.85, df["p_loss_per_generation"].nunique()))
    for (p, grp), color in zip(df.groupby("p_loss_per_generation"), palette):
        ax1.plot(grp["generation"], grp["R_system"],
                 label=f"p = {p:.2f}/gen", color=color, linewidth=2)

    ax1.axhline(0.5, color="red", linestyle=":", linewidth=1.3,
                label="50% reliability threshold")
    ax1.set_xlabel("Generations (crop cycles)", fontsize=11)
    ax1.set_ylabel("System reliability R$_{\\mathrm{system}}(g)$", fontsize=11)
    ax1.set_title("N-cycle system reliability vs crop generations\n"
                  "(AMO × NXR series chain, taxa-parallel within blocks)",
                  fontsize=11, fontweight="bold")
    ax1.legend(fontsize=9, loc="lower left")
    ax1.grid(alpha=0.3)
    ax1.set_xlim(0, 30); ax1.set_ylim(0, 1.02)

    # Per-block detail at p = 0.10
    sub = df[df["p_loss_per_generation"] == 0.10]
    ax2.plot(sub["generation"], sub["R_AMO"],
             label="AMO block (2 parallel taxa)", linewidth=2.5, color="#d62728")
    ax2.plot(sub["generation"], sub["R_NXR"],
             label="NXR block (2 parallel taxa)", linewidth=2.5, color="#ff7f0e")
    ax2.plot(sub["generation"], sub["R_DENIT"],
             label="Denit block (10 parallel taxa)\n[shown for comparison]",
             linewidth=2.5, color="#2ca02c", linestyle="--")
    ax2.plot(sub["generation"], sub["R_system"],
             label="Full N-cycle system\n(AMO × NXR)",
             linewidth=3, color="black")
    ax2.axhline(0.5, color="red", linestyle=":", linewidth=1.3)
    ax2.set_xlabel("Generations (crop cycles)", fontsize=11)
    ax2.set_ylabel("Block reliability", fontsize=11)
    ax2.set_title("Per-block reliability decomposition\n(at p_loss = 0.10/gen/taxon)",
                  fontsize=11, fontweight="bold")
    ax2.legend(fontsize=9, loc="lower left")
    ax2.grid(alpha=0.3)
    ax2.set_xlim(0, 30); ax2.set_ylim(0, 1.02)

    fig.suptitle("MELiSSA-style reliability analysis of the ISS pepper N-cycle",
                 fontsize=13, fontweight="bold")
    for ext in ("pdf", "png"):
        fig.savefig(OUT_DIR / f"rbd_reliability_curves.{ext}",
                    dpi=300 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: rbd_reliability_curves.pdf/.png")


def main():
    print("=" * 60)
    print("48: Reliability Block Diagram simulation")
    print("=" * 60)

    print("\n[1] Block definitions...")
    write_block_table()

    print("\n[2] Parametric sweep...")
    df, summary = parametric_sweep()

    print("\n[3] Drawing block diagram...")
    draw_block_diagram()

    print("\n[4] Drawing reliability curves...")
    draw_reliability_curves(df)

    print(f"\nDone. Outputs in {OUT_DIR}/")


if __name__ == "__main__":
    main()
