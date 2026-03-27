#!/usr/bin/env python3
"""
Supplementary Figure S1: Faith's PD violin + boxplot + jitter
(Space Flight vs Terrestrial Soil)

Output:
  version-2_integrated/SuppS1_FaithPD.pdf
  version-2_integrated/SuppS1_FaithPD.png
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
import config

BASE    = Path(__file__).parent / 'version-2_integrated'
OUT_DIR = config.VERSION_DIR
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Load data ──────────────────────────────────────────────────────────────
# Only metadata + alpha-diversity needed; feature-table not required here.
meta  = pd.read_csv(BASE / 'integrated_metadata.tsv', sep='\t', index_col=0)
alpha = pd.read_csv(BASE / 'exported_diversity/alpha-diversity.tsv',
                    sep='\t', index_col=0)
alpha.columns = ['faith_pd']

space_ids = meta.index[meta['Study_Group'] == 'Space_Flight'].tolist()
terr_ids  = meta.index[meta['Study_Group'] == 'Terrestrial_Soil'].tolist()

# ── Group assignment ───────────────────────────────────────────────────────
pd_df = alpha.copy()
pd_df['group'] = None
pd_df.loc[pd_df.index.isin(space_ids), 'group'] = 'Space Flight'
pd_df.loc[pd_df.index.isin(terr_ids),  'group'] = 'Terrestrial'
pd_df = pd_df.dropna(subset=['group'])

space_pd = pd_df.loc[pd_df['group'] == 'Space Flight', 'faith_pd'].values
terr_pd  = pd_df.loc[pd_df['group'] == 'Terrestrial',  'faith_pd'].values

n_sp = len(space_pd)
n_te = len(terr_pd)
print(f"n  Space: {n_sp}, Terrestrial: {n_te}")
print(f"   Space  mean ± SD: {space_pd.mean():.2f} ± {space_pd.std(ddof=1):.2f}")
print(f"   Terr.  mean ± SD: {terr_pd.mean():.2f} ± {terr_pd.std(ddof=1):.2f}")
print(f"   (Note: Space Flight n={n_sp} of 106 total; {106 - n_sp} excluded — "
      "not present in alpha-diversity output after rarefaction/depth filtering)")

# ── Statistics ────────────────────────────────────────────────────────────
res  = mannwhitneyu(space_pd, terr_pd, alternative='two-sided')
pval = res.pvalue
W    = res.statistic

# Rank-biserial correlation: (2W / n1*n2) - 1
# Negative when Space (group1) tends to be LOWER than Terrestrial (group2)
r_rb = (2 * W / (n_sp * n_te)) - 1

ptxt = 'p < 0.001' if pval < 0.001 else f'p = {pval:.3f}'
rtxt = f'$r_{{rb}}$ = {r_rb:.2f}'
print(f"   W = {W:.0f}, p = {pval:.2e}, r_rb = {r_rb:.2f}")

# ── Plot ───────────────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family':     'DejaVu Serif',
    'axes.titlesize':  11,
    'axes.labelsize':  10,
    'xtick.labelsize': 9.5,
    'ytick.labelsize': 9,
})

fig, ax = plt.subplots(figsize=(4.4, 5.8))

COL_SP = '#C0392B'   # space  — red
COL_TE = '#2E6DB4'   # terr.  — blue

# Violin
parts = ax.violinplot([space_pd, terr_pd], positions=[1, 2], widths=0.72,
                      showmeans=False, showmedians=False, showextrema=False)
for body, col in zip(parts['bodies'], [COL_SP, COL_TE]):
    body.set_facecolor(col)
    body.set_alpha(0.22)
    body.set_edgecolor('none')

# Boxplot
bp = ax.boxplot(
    [space_pd, terr_pd], positions=[1, 2], widths=0.34,
    patch_artist=True,
    medianprops=dict(color='black', lw=1.4),
    whiskerprops=dict(color='#555555', lw=1.0),
    capprops=dict(color='#555555', lw=1.0),
    flierprops=dict(marker='o', markersize=3.0, alpha=0.45,
                    markerfacecolor='#888888', markeredgecolor='none'),
    boxprops=dict(lw=1.0),
)
for patch, col in zip(bp['boxes'], [COL_SP, COL_TE]):
    patch.set_facecolor(col)
    patch.set_alpha(0.78)
    patch.set_edgecolor('#404040')

# Jitter
rng = np.random.default_rng(42)
for pos, vals, col in zip([1, 2], [space_pd, terr_pd], [COL_SP, COL_TE]):
    jitter = rng.uniform(-0.09, 0.09, size=len(vals))
    ax.scatter(np.full_like(vals, pos) + jitter, vals,
               s=13, c='#333333', alpha=0.28, zorder=4)

# Significance bracket + p-value + r_rb
ytop = float(max(pd_df['faith_pd']))
y_line = ytop * 1.02
y_p    = ytop * 1.065
y_r    = ytop * 1.115

ax.plot([1, 2], [y_line, y_line], color='#333333', lw=0.9)
ax.text(1.5, y_p, ptxt, ha='center', va='bottom', fontsize=9.2, color='#222222')
ax.text(1.5, y_r, rtxt, ha='center', va='bottom', fontsize=8.6, color='#666666')

# x-axis labels with n=
ax.set_xticks([1, 2])
ax.set_xticklabels([f'Space\nFlight\n(n={n_sp})', f'Terrestrial\n(n={n_te})'],
                   fontsize=9.2)

ax.set_ylabel("Faith's Phylogenetic Diversity", fontsize=10)
ax.set_title("Faith's PD distribution by habitat", fontsize=11, fontweight='bold', pad=8)

# Axis style
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.grid(axis='y', color='#EAEAEA', lw=0.8)
ax.set_axisbelow(True)
ax.set_xlim(0.4, 2.6)
ax.set_ylim(bottom=max(0, float(min(pd_df['faith_pd'])) - 1),
            top=ytop * 1.17)

fig.tight_layout()

out_pdf = OUT_DIR / 'SuppS1_FaithPD.pdf'
out_png = OUT_DIR / 'SuppS1_FaithPD.png'
fig.savefig(out_pdf, bbox_inches='tight')
fig.savefig(out_png, dpi=320, bbox_inches='tight')
print(f"Saved:\n  {out_pdf}\n  {out_png}")

# Current Biology supplementary figure output
cb_dir = config.CURRENT_BIOLOGY_DIR
cb_dir.mkdir(exist_ok=True)
fig.savefig(cb_dir / 'SuppS1a_FaithPD.png', dpi=320, bbox_inches='tight')
fig.savefig(cb_dir / 'SuppS1a_FaithPD.pdf', bbox_inches='tight')
print(f"Saved Current Biology:\n  {cb_dir / 'SuppS1a_FaithPD.png'}\n  {cb_dir / 'SuppS1a_FaithPD.pdf'}")
