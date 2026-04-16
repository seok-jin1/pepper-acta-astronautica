"""
Fig 1 — Study design and BLSS system context.
(a) ISS Veggie hardware schematic (stylized layered view).
(b) Sample lineage flowchart.
(c) Analysis pipeline DAG.
All panels drawn in matplotlib for reproducibility.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, FancyArrowPatch, Ellipse, Wedge
import matplotlib.lines as mlines
import config

OUT = config.AA_FIG_DIR


# ---------------------------------------------------------------------------
# Panel a — Veggie hardware cross-section
# ---------------------------------------------------------------------------
def panel_a_veggie_hardware(ax):
    ax.set_xlim(0, 10); ax.set_ylim(0, 10)
    ax.axis("off")

    # LED canopy
    led_y = 8.0
    for xx in range(2, 9):
        ax.add_patch(Rectangle((xx - 0.25, led_y), 0.5, 0.4,
                               facecolor="#ffeb3b", edgecolor="#f57f17"))
    ax.text(5, led_y + 0.7, "LED canopy (red / blue / green)",
            ha="center", fontsize=8, color="#e65100", fontweight="bold")

    # Airflow arrows
    for xx in (1.2, 8.8):
        ax.add_patch(FancyArrowPatch((xx, 7.5), (xx, 3.0),
                                     arrowstyle="->",
                                     mutation_scale=15, color="#1976d2",
                                     linewidth=1.2, alpha=0.65))
    ax.text(0.8, 5.2, "air\nflow", fontsize=7, color="#1976d2",
            ha="center", fontstyle="italic")
    ax.text(9.2, 5.2, "air\nflow", fontsize=7, color="#1976d2",
            ha="center", fontstyle="italic")

    # Capsicum plant (stem + leaves)
    ax.plot([5, 5], [3.5, 7.2], color="#4caf50", linewidth=2)
    for dy, dx in [(0.6, 0.9), (1.5, -0.9), (2.4, 0.9), (3.2, -0.6)]:
        ax.add_patch(Ellipse((5 + dx, 3.5 + dy), 1.2, 0.45,
                             facecolor="#66bb6a", edgecolor="#2e7d32",
                             angle=20 * (1 if dx > 0 else -1)))
    # Fruit
    ax.add_patch(Ellipse((5.5, 5.5), 0.4, 0.7, facecolor="#e53935",
                         edgecolor="#b71c1c"))

    # Plant pillow (arcillite substrate)
    ax.add_patch(FancyBboxPatch((2.2, 2.2), 5.6, 1.2,
                                boxstyle="round,pad=0.05",
                                facecolor="#8d6e63", edgecolor="#3e2723"))
    ax.text(5, 2.8, "Arcillite substrate (plant pillow)",
            ha="center", fontsize=8, color="white", fontweight="bold")

    # Wick
    ax.add_patch(Rectangle((4.7, 1.0), 0.6, 1.3,
                           facecolor="#ffcc80", edgecolor="#e65100"))
    ax.text(4.4, 1.6, "wick", fontsize=7, color="#e65100", ha="right")

    # Water reservoir
    ax.add_patch(FancyBboxPatch((2.2, 0.3), 5.6, 0.8,
                                boxstyle="round,pad=0.05",
                                facecolor="#4fc3f7", edgecolor="#01579b"))
    ax.text(5, 0.7, "Water reservoir",
            ha="center", fontsize=8, color="#01579b", fontweight="bold")

    # Chamber outline
    ax.add_patch(Rectangle((1, 0), 8, 9, fill=False,
                           edgecolor="#424242", linewidth=1.5))
    ax.text(5, 9.4, "ISS Veggie growth chamber",
            ha="center", fontsize=10, fontweight="bold")


# ---------------------------------------------------------------------------
# Panel b — Sample lineage flowchart
# ---------------------------------------------------------------------------
def panel_b_sample_lineage(ax):
    ax.set_xlim(0, 10); ax.set_ylim(0, 10)
    ax.axis("off")
    ax.set_title("b) Sample lineage", loc="left", fontweight="bold", fontsize=10)

    nodes = [
        # (x, y, text, colour, kind)
        (5.0, 9.0, "Seed\n(pre-flight)", "#fff9c4", "source"),
        (2.0, 7.0, "Arcillite\nsubstrate", "#bcaaa4", "source"),
        (5.0, 7.0, "Wick", "#ffcc80", "source"),
        (8.0, 7.0, "Water\nsystem", "#81d4fa", "source"),
        (3.5, 4.5, "Foam /\nPillow surface", "#d7ccc8", "source"),
        (6.5, 4.5, "Swab /\ncargo surface", "#cfd8dc", "source"),
        (3.5, 2.2, "Plant roots", "#8bc34a", "sink"),
        (6.5, 2.2, "Rhizosphere\n(pepper, lettuce)", "#558b2f", "sink"),
    ]
    for x, y, t, c, kind in nodes:
        shape = "round,pad=0.1"
        ax.add_patch(FancyBboxPatch((x - 0.9, y - 0.55), 1.8, 1.1,
                                    boxstyle=shape, facecolor=c,
                                    edgecolor="#424242",
                                    linewidth=1.2 if kind == "sink" else 0.8))
        ax.text(x, y, t, ha="center", va="center", fontsize=7.5,
                fontweight="bold" if kind == "sink" else "normal")

    edges = [
        (5.0, 8.5, 5.0, 7.5),    # seed -> wick
        (5.0, 8.5, 2.0, 7.5),    # seed -> substrate
        (2.0, 6.5, 3.3, 5.0),    # substrate -> roots
        (5.0, 6.5, 3.7, 5.0),    # wick -> roots
        (8.0, 6.5, 6.8, 5.0),    # water -> rhizo
        (3.5, 4.0, 3.5, 2.7),    # foam -> roots
        (6.5, 4.0, 6.5, 2.7),    # swab -> rhizo
        (3.5, 2.0, 6.3, 2.0),    # roots -> rhizo
    ]
    for xa, ya, xb, yb in edges:
        ax.add_patch(FancyArrowPatch((xa, ya), (xb, yb),
                                     arrowstyle="->",
                                     mutation_scale=10, color="#616161",
                                     linewidth=0.9))
    ax.text(5.0, 0.8, "FEAST source-tracking (scripts 22–24, 39)",
            ha="center", fontsize=8, style="italic", color="#424242")


# ---------------------------------------------------------------------------
# Panel c — Analysis pipeline DAG
# ---------------------------------------------------------------------------
def panel_c_pipeline_dag(ax):
    ax.set_xlim(0, 10); ax.set_ylim(0, 10)
    ax.axis("off")
    ax.set_title("c) Analysis pipeline", loc="left", fontweight="bold", fontsize=10)

    nodes = [
        (5.0, 9.2, "16S rRNA FASTQ\n(pepper + lettuce + ISS env)", "#e1f5fe"),
        (5.0, 7.8, "DADA2 denoise + SILVA 138", "#b3e5fc"),
        (1.8, 6.2, "FAPROTAX\n(script 12)", "#c5e1a5"),
        (4.2, 6.2, "PICRUSt2\n(script 07)", "#dcedc8"),
        (6.6, 6.2, "GTDB / NCBI\nN-cycle genes\n(script 36)", "#aed581"),
        (9.0, 6.2, "βNTI / RCbray\n(script 13)", "#ffe082"),
        (3.0, 4.2, "Consensus\n(script 35)", "#a5d6a7"),
        (5.8, 4.2, "Genome matrix\nFig 3b", "#81c784"),
        (8.2, 4.2, "Assembly\nFig 6a", "#ffd54f"),
        (1.5, 2.4, "Network\n(scripts 9, 17)", "#fdd0a2"),
        (4.0, 2.4, "FEAST\n(scripts 22–24, 39)", "#ffab91"),
        (6.5, 2.4, "Redundancy\n(script 14)", "#d7ccc8"),
        (9.0, 2.4, "N-cycle\n(script 16)", "#b0bec5"),
        (5.0, 0.6, "BLSS failure-mode synthesis (Fig 6b)", "#e57373"),
    ]
    for x, y, t, c in nodes:
        ax.add_patch(FancyBboxPatch((x - 1.2, y - 0.55), 2.4, 1.1,
                                    boxstyle="round,pad=0.08",
                                    facecolor=c, edgecolor="#424242"))
        ax.text(x, y, t, ha="center", va="center", fontsize=7.5)

    edges = [
        (5.0, 8.6, 5.0, 8.3),
        (5.0, 7.3, 1.8, 6.7),
        (5.0, 7.3, 4.2, 6.7),
        (5.0, 7.3, 6.6, 6.7),
        (5.0, 7.3, 9.0, 6.7),
        (1.8, 5.7, 3.0, 4.7),
        (4.2, 5.7, 3.0, 4.7),
        (6.6, 5.7, 5.8, 4.7),
        (9.0, 5.7, 8.2, 4.7),
        (1.8, 5.7, 1.5, 2.9),
        (4.2, 5.7, 4.0, 2.9),
        (1.8, 5.7, 6.5, 2.9),
        (1.8, 5.7, 9.0, 2.9),
        (3.0, 3.7, 5.0, 1.1),
        (5.8, 3.7, 5.0, 1.1),
        (8.2, 3.7, 5.0, 1.1),
    ]
    for xa, ya, xb, yb in edges:
        ax.add_patch(FancyArrowPatch((xa, ya), (xb, yb),
                                     arrowstyle="->",
                                     mutation_scale=8, color="#616161",
                                     linewidth=0.7))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(22, 9),
                             gridspec_kw={"width_ratios": [0.9, 1.1, 1.4]})
    panel_a_veggie_hardware(axes[0])
    axes[0].set_title("a) ISS Veggie growth chamber",
                      loc="left", fontweight="bold", fontsize=10)
    panel_b_sample_lineage(axes[1])
    panel_c_pipeline_dag(axes[2])
    plt.tight_layout()
    fig.suptitle("Figure 1. Study design and BLSS system context",
                 fontsize=13, fontweight="bold", y=1.02)
    out = OUT / "Fig1_system_context"
    fig.savefig(str(out) + ".png", dpi=300, bbox_inches="tight")
    fig.savefig(str(out) + ".pdf", bbox_inches="tight")
    print(f"Wrote {out}.png/.pdf")


if __name__ == "__main__":
    main()
