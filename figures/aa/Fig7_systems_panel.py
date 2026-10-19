"""
Fig7_systems_panel.py
======================
Systems panel: observation → risk → monitoring indicator → mitigation response.

Flow:
  [Nitrifier depletion detected]
      → [amoA structural bottleneck risk]
          → [amoA/qPCR monitoring indicator]
              → [Substrate-targeted inoculation contingency]
  [ISS surface reservoir absent]
      → [Passive recolonization implausible]
          → (same monitoring / mitigation)

Output: figures/aa/Fig7_systems_panel.pdf / .png
"""

from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

OUTPUT_DIR = Path(__file__).parent

# ── Color palette ────────────────────────────────────────────────────────────
C_OBS     = '#2C3E50'   # dark navy  — observation
C_RISK    = '#C0392B'   # red        — risk assessment
C_MONITOR = '#E67E22'   # orange     — monitoring indicator
C_MITIG   = '#27AE60'   # green      — mitigation / design response
C_ARROW   = '#7F8C8D'   # grey       — arrows
C_BG      = '#FDFEFE'


def draw_box(ax, x, y, w, h, text, color, fontsize=8.5, text_color='white', style='round,pad=0.1'):
    box = FancyBboxPatch(
        (x - w/2, y - h/2), w, h,
        boxstyle=style,
        facecolor=color, edgecolor='white', linewidth=1.5,
        zorder=3
    )
    ax.add_patch(box)
    ax.text(x, y, text, ha='center', va='center',
            fontsize=fontsize, color=text_color,
            fontweight='bold', wrap=True,
            multialignment='center', zorder=4)


def draw_arrow(ax, x1, y1, x2, y2, label=''):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(
                    arrowstyle='->', color=C_ARROW,
                    lw=1.8, connectionstyle='arc3,rad=0.0'
                ), zorder=2)
    if label:
        mx, my = (x1+x2)/2 + 0.02, (y1+y2)/2
        ax.text(mx, my, label, ha='left', va='center',
                fontsize=7, color=C_ARROW, style='italic')


fig, ax = plt.subplots(figsize=(12, 7))
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis('off')
ax.set_facecolor(C_BG)
fig.patch.set_facecolor(C_BG)

# ── Column headers ────────────────────────────────────────────────────────────
col_x  = [0.13, 0.37, 0.63, 0.87]
col_labels = ['OBSERVATION', 'RISK ASSESSMENT', 'MONITORING INDICATOR', 'MITIGATION / DESIGN RESPONSE']
col_colors = [C_OBS, C_RISK, C_MONITOR, C_MITIG]

for cx, cl, cc in zip(col_x, col_labels, col_colors):
    ax.text(cx, 0.95, cl, ha='center', va='center',
            fontsize=8.5, fontweight='bold', color=cc,
            bbox=dict(boxstyle='round,pad=0.3', facecolor=cc, alpha=0.15, edgecolor=cc))

# ── Row 1: Main nitrifier depletion chain ────────────────────────────────────
y1 = 0.72
bw, bh = 0.20, 0.14

draw_box(ax, col_x[0], y1, bw, bh,
         'Nitrifier genera\n(Nitrosospira, Nitrospira)\nundetected: 0/20 ISS\nvs 19/20 Terrestrial',
         C_OBS, fontsize=7.8)

draw_box(ax, col_x[1], y1, bw, bh,
         'amoA structural\nbottleneck:\nonly 2/13 genera\ncarry amoA (GTDB)',
         C_RISK, fontsize=7.8)

draw_box(ax, col_x[2], y1, bw, bh,
         'amoA abundance\nor transcript\n(qPCR / metatranscriptome)\nin rhizosphere',
         C_MONITOR, fontsize=7.8)

draw_box(ax, col_x[3], y1, bw, bh,
         'Substrate-embedded\nNitrosomonas +\nNitrospira consortium\ninoculation',
         C_MITIG, fontsize=7.8)

for i in range(3):
    draw_arrow(ax, col_x[i]+bw/2+0.01, y1, col_x[i+1]-bw/2-0.01, y1)

# ── Row 2: ISS reservoir absence chain ──────────────────────────────────────
y2 = 0.47
draw_box(ax, col_x[0], y2, bw, bh,
         'ISS-MO surface pool\n(OSD-65, n=196)\ncontains zero\nnitrifier genera',
         C_OBS, fontsize=7.8)

draw_box(ax, col_x[1], y2, bw, bh,
         'Passive recolonization\nfrom ISS environment\nimplausible:\nno source taxa present',
         C_RISK, fontsize=7.8)

draw_box(ax, col_x[2], y2, bw, bh,
         'Longitudinal\nnitrifier prevalence\ntracking across\nmission duration',
         C_MONITOR, fontsize=7.8)

draw_box(ax, col_x[3], y2, bw, bh,
         'Pre-flight substrate\nsterilization +\ndefined inoculum\n(closed-loop)',
         C_MITIG, fontsize=7.8)

for i in range(3):
    draw_arrow(ax, col_x[i]+bw/2+0.01, y2, col_x[i+1]-bw/2-0.01, y2)

# ── Row 3: Crop specificity ──────────────────────────────────────────────────
y3 = 0.22
draw_box(ax, col_x[0], y3, bw, bh,
         'Nitrification\nLog₂FC diverges:\npeper −13.08\nvs lettuce +11.45',
         C_OBS, fontsize=7.8)

draw_box(ax, col_x[1], y3, bw, bh,
         'Nitrogen-cycle\nrisk is crop-specific;\nsingle-species\ncharacterisation insufficient',
         C_RISK, fontsize=7.8)

draw_box(ax, col_x[2], y3, bw, bh,
         'Per-crop nitrifier\nprevalence / amoA\nprofiling before\nmission deployment',
         C_MONITOR, fontsize=7.8)

draw_box(ax, col_x[3], y3, bw, bh,
         'Mandatory per-crop\nbiological risk\nassessment protocol\nfor BLSS deployment',
         C_MITIG, fontsize=7.8)

for i in range(3):
    draw_arrow(ax, col_x[i]+bw/2+0.01, y3, col_x[i+1]-bw/2-0.01, y3)

# ── Vertical connectors between rows ────────────────────────────────────────
for cx in col_x:
    draw_arrow(ax, cx, y1 - bh/2 - 0.01, cx, y2 + bh/2 + 0.01)
    draw_arrow(ax, cx, y2 - bh/2 - 0.01, cx, y3 + bh/2 + 0.01)

# ── Title and note ────────────────────────────────────────────────────────────
ax.set_title(
    'BLSS nitrogen-cycle risk framework: observation → risk → monitoring → mitigation',
    fontsize=11, fontweight='bold', pad=10, color=C_OBS
)
ax.text(0.5, 0.04,
        'Evidence strength: Row 1 (direct ASV detection + genome-level gene evidence) is the strongest tier.\n'
        'Rows 2–3 provide supporting context and generalization scope. All mitigation implications are \n'
        'consistent with but not proven by the present observational data.',
        ha='center', va='center', fontsize=7.5, color='#5D6D7E', style='italic',
        transform=ax.transAxes)

fig.tight_layout()

for ext in ('pdf', 'png'):
    p = OUTPUT_DIR / f'Fig7_systems_panel.{ext}'
    fig.savefig(p, dpi=300, bbox_inches='tight')
    print(f'Saved: {p}')

plt.close(fig)
print('Done.')
