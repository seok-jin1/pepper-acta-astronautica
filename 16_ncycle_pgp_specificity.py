"""
N-Cycle Completeness Score + PGP Functional Specificity Index
=============================================================
Two complementary metrics that quantify functional capacity per sample:

1. N-Cycle Completeness Score
   - 5 steps: N2-fixation → nitrification → denitrification →
              nitrate_reduction → ureolysis
   - Score per sample = fraction of steps with ≥1 genus present (0–1)
   - Compared between Space Flight and Terrestrial Soil

2. PGP Functional Specificity Index
   - Measures diversity of PGP functional categories present in a sample
   - = Shannon entropy (H') across 5 PGP categories weighted by abundance
   - H' = 0 → only one PGP category present
   - H' = log(5) → all 5 categories equally represented
   - Captures whether PGP capacity is broad-spectrum or narrowly specialized

Ecological rationale:
  An N-cycle completeness score near 1.0 indicates the microbial community can
  sustain the full nitrogen cycle. Space samples losing multiple N-cycle steps
  is strong evidence of functional network collapse. PGP specificity captures
  whether any remaining PGP capacity is well-rounded or single-function.

Prerequisites:
  - version-2_integrated/FAPROTAX_report.txt
  - version-2_integrated/exported_table_clean/feature-table.tsv
  - version-2_integrated/exported_taxonomy/taxonomy.tsv
  - version-2_integrated/integrated_metadata.tsv

Run:
  python 26_ncycle_pgp_specificity.py

Output (version-2_integrated/):
  - NCycle_Completeness_Results.csv
  - NCycle_Completeness_Summary.csv
  - PGP_Specificity_Results.csv
  - PGP_Specificity_Summary.csv
  - Main_Fig_NCycle_PGP.png/.pdf
"""

import re
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import config

import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

OUTPUT_DIR = config.VERSION_DIR

# ── N-cycle steps (FAPROTAX function names) ──────────────────────────────────
NCYCLE_STEPS = [
    'nitrogen_fixation',
    'nitrification',
    'denitrification',
    'nitrate_reduction',
    'ureolysis',
]
NCYCLE_LABELS = {
    'nitrogen_fixation': 'N₂ Fixation',
    'nitrification':     'Nitrification',
    'denitrification':   'Denitrification',
    'nitrate_reduction': 'Nitrate Reduction',
    'ureolysis':         'Ureolysis',
}

# ── PGP genera by functional category ────────────────────────────────────────
PGP_CATEGORIES = {
    'N_Fixation':        ['Rhizobium', 'Bradyrhizobium', 'Mesorhizobium',
                          'Sinorhizobium', 'Azorhizobium', 'Frankia',
                          'Azotobacter', 'Herbaspirillum', 'Gluconacetobacter'],
    'P_Solubilization':  ['Bacillus', 'Pseudomonas', 'Burkholderia',
                          'Aspergillus', 'Trichoderma', 'Rhizobium',
                          'Enterobacter', 'Serratia'],
    'IAA_Phytohormone':  ['Pseudomonas', 'Bacillus', 'Azospirillum',
                          'Rhizobium', 'Streptomyces', 'Enterobacter',
                          'Herbaspirillum', 'Sphingomonas'],
    'Biocontrol':        ['Bacillus', 'Pseudomonas', 'Streptomyces',
                          'Lysobacter', 'Paenibacillus', 'Serratia',
                          'Trichoderma'],
    'ACC_Deaminase':     ['Pseudomonas', 'Bacillus', 'Burkholderia',
                          'Enterobacter', 'Herbaspirillum', 'Rhizobium',
                          'Streptomyces'],
}

COLORS = {'Space_Flight': '#3498db', 'Terrestrial_Soil': '#2ecc71'}
GROUP_LABELS = {'Space_Flight': 'Space Flight', 'Terrestrial_Soil': 'Terrestrial Soil'}


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def parse_faprotax_report(report_path):
    func2genera = {}
    current_func = None
    with open(report_path) as fh:
        for line in fh:
            line = line.rstrip()
            m = re.match(r'^# (\w+) \(\d+ records\):', line)
            if m:
                current_func = m.group(1)
                func2genera[current_func] = set()
                continue
            if current_func and line.startswith('    '):
                for part in line.split(';'):
                    part = part.strip()
                    if part.startswith('g__'):
                        name = part.split('__', 1)[-1].strip()
                        if name and name not in ('uncultured', 'metagenome',
                                                 'unidentified', ''):
                            func2genera[current_func].add(name)
    return func2genera


def aggregate_to_genus(feature_table_path, taxonomy_path):
    table = pd.read_csv(feature_table_path, sep='\t', skiprows=1, index_col=0)
    tax   = pd.read_csv(taxonomy_path, sep='\t', index_col=0)

    def get_genus(s):
        if pd.isna(s):
            return 'Unassigned'
        for part in reversed(s.split(';')):
            part = part.strip()
            if part.startswith('g__'):
                name = part.split('__', 1)[-1].strip()
                return name if name else 'Unassigned'
        return 'Unassigned'

    tax['Genus'] = tax['Taxon'].apply(get_genus)
    shared = table.index.intersection(tax.index)
    table  = table.loc[shared]
    table['Genus'] = tax.loc[shared, 'Genus']
    genus_table = table.groupby('Genus').sum()
    genus_table = genus_table.drop(index='Unassigned', errors='ignore')
    return genus_table.div(genus_table.sum(axis=0), axis=1) * 100


def mwu_compare(sp_vals, tr_vals):
    stat, p = mannwhitneyu(sp_vals, tr_vals, alternative='two-sided')
    return stat, p


def sig_label(p):
    if p < 0.001:
        return '***'
    elif p < 0.01:
        return '**'
    elif p < 0.05:
        return '*'
    return 'n.s.'


# ─────────────────────────────────────────────────────────────────────────────
# 1. N-Cycle Completeness
# ─────────────────────────────────────────────────────────────────────────────

def compute_ncycle_completeness(genus_rel, func2genera, meta):
    """
    Per sample: fraction of N-cycle steps with ≥1 genus present (0 to 1).
    Also returns per-step presence (0/1) for each sample.
    """
    rows = {}
    for sample in genus_rel.columns:
        step_scores = {}
        for step in NCYCLE_STEPS:
            assigned = func2genera.get(step, set())
            present  = assigned.intersection(genus_rel.index)
            if len(present) == 0:
                step_scores[step] = 0
            else:
                has_step = int((genus_rel.loc[list(present), sample] > 0).any())
                step_scores[step] = has_step
        rows[sample] = step_scores
    df = pd.DataFrame(rows).T  # samples × steps
    df['Completeness'] = df[NCYCLE_STEPS].mean(axis=1)
    df['Study_Group']  = meta.loc[df.index.intersection(meta.index), 'Study_Group']
    return df


# ─────────────────────────────────────────────────────────────────────────────
# 2. PGP Functional Specificity Index
# ─────────────────────────────────────────────────────────────────────────────

def compute_pgp_specificity(genus_rel, meta):
    """
    For each sample:
      - For each PGP category, sum relative abundance of PGP genera in category
        that are present in the genus table.
      - Category abundance = sum over present genera.
      - Shannon H' = -sum(p * log(p)) over 5 categories (p = proportion, ignore 0)
      - Specificity index = H'
    Returns DataFrame: index=samples, columns=[5 categories, Total_PGP_Abund, Specificity_H]
    """
    rows = {}
    for sample in genus_rel.columns:
        cat_abund = {}
        for cat, genera_list in PGP_CATEGORIES.items():
            present = genus_rel.index.intersection(genera_list)
            cat_abund[cat] = float(genus_rel.loc[present, sample].sum()) if len(present) else 0.0
        total = sum(cat_abund.values())
        if total > 0:
            props = np.array([v / total for v in cat_abund.values()])
            h = -np.nansum(props * np.log(props + 1e-12))  # Shannon
        else:
            h = 0.0
        row = dict(cat_abund)
        row['Total_PGP_Abund'] = total
        row['Specificity_H']   = h
        rows[sample] = row
    df = pd.DataFrame(rows).T
    df['Study_Group'] = meta.loc[df.index.intersection(meta.index), 'Study_Group']
    return df


# ─────────────────────────────────────────────────────────────────────────────
# 3. Plot
# ─────────────────────────────────────────────────────────────────────────────

def plot_results(ncycle_df, pgp_df, meta):
    sp_idx = meta.index[meta['Study_Group'] == 'Space_Flight']
    tr_idx = meta.index[meta['Study_Group'] == 'Terrestrial_Soil']

    fig, axes = plt.subplots(2, 2, figsize=(14, 10), constrained_layout=True)

    # ── A: N-cycle completeness boxplot ─────────────────────────────────────
    ax = axes[0, 0]
    nc_df = ncycle_df.dropna(subset=['Study_Group'])
    sp_comp = nc_df.loc[nc_df.index.intersection(sp_idx), 'Completeness']
    tr_comp = nc_df.loc[nc_df.index.intersection(tr_idx), 'Completeness']
    _, p_nc = mwu_compare(sp_comp, tr_comp)

    plot_data = pd.DataFrame({
        'Completeness': pd.concat([sp_comp, tr_comp]),
        'Group': (['Space Flight'] * len(sp_comp) + ['Terrestrial Soil'] * len(tr_comp))
    })
    for i, (grp, col) in enumerate([('Space Flight', '#3498db'), ('Terrestrial Soil', '#2ecc71')]):
        vals = plot_data[plot_data['Group'] == grp]['Completeness']
        bp = ax.boxplot(vals, positions=[i], widths=0.4,
                        patch_artist=True,
                        boxprops=dict(facecolor=col, alpha=0.7),
                        medianprops=dict(color='black', linewidth=2),
                        whiskerprops=dict(linewidth=1.5),
                        capprops=dict(linewidth=1.5),
                        flierprops=dict(marker='o', markersize=3, alpha=0.4))
        ax.scatter([i + np.random.uniform(-0.15, 0.15) for _ in vals],
                   vals, alpha=0.3, s=15, color='#555555')

    y_top = max(sp_comp.max(), tr_comp.max()) + 0.05
    ax.plot([0, 1], [y_top, y_top], color='black', lw=1)
    ax.text(0.5, y_top + 0.01, sig_label(p_nc), ha='center', fontsize=12, fontweight='bold')
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['Space Flight', 'Terrestrial Soil'], fontsize=11)
    ax.set_ylabel('N-Cycle Completeness Score (0–1)', fontsize=11)
    ax.set_title('N-Cycle Completeness\n(fraction of 5 steps with ≥1 genus)', fontsize=11, fontweight='bold')
    ax.set_ylim(0, 1.2)
    ax.grid(axis='y', alpha=0.3)
    ax.text(0.98, 0.02, f'p = {p_nc:.4f}', transform=ax.transAxes,
            ha='right', va='bottom', fontsize=9, color='gray')

    # ── B: N-cycle step presence heatmap ────────────────────────────────────
    ax = axes[0, 1]
    step_means = pd.DataFrame({
        'Space Flight':     nc_df.loc[nc_df.index.intersection(sp_idx), NCYCLE_STEPS].mean(),
        'Terrestrial Soil': nc_df.loc[nc_df.index.intersection(tr_idx), NCYCLE_STEPS].mean(),
    }, index=NCYCLE_STEPS).rename(index=NCYCLE_LABELS)

    import seaborn as sns
    import matplotlib.colors as mcolors

    # Use RdYlGn colormap and create heatmap
    cmap = plt.cm.get_cmap('RdYlGn')
    norm = mcolors.Normalize(vmin=0, vmax=1)

    sns.heatmap(step_means, ax=ax, annot=False,
                cmap='RdYlGn', vmin=0, vmax=1,
                linewidths=0.5, linecolor='white',
                cbar_kws={'label': 'Fraction of Samples with Step Present'})

    # Manually add text annotations (all black)
    for i, row_label in enumerate(step_means.index):
        for j, col_label in enumerate(step_means.columns):
            val = step_means.iloc[i, j]
            ax.text(j + 0.5, i + 0.5, f'{val:.2f}',
                   ha='center', va='center', color='black',
                   fontsize=10, fontweight='bold')

    ax.set_title('N-Cycle Step Presence\n(fraction of samples per group)', fontsize=11, fontweight='bold')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0, fontsize=11)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=10)

    # ── C: PGP specificity boxplot ───────────────────────────────────────────
    ax = axes[1, 0]
    pgp_filt = pgp_df.dropna(subset=['Study_Group'])
    sp_pgp = pgp_filt.loc[pgp_filt.index.intersection(sp_idx), 'Specificity_H']
    tr_pgp = pgp_filt.loc[pgp_filt.index.intersection(tr_idx), 'Specificity_H']
    _, p_pgp = mwu_compare(sp_pgp, tr_pgp)

    for i, (grp, vals, col) in enumerate([
        ('Space Flight', sp_pgp, '#3498db'),
        ('Terrestrial Soil', tr_pgp, '#2ecc71')
    ]):
        ax.boxplot(vals, positions=[i], widths=0.4,
                   patch_artist=True,
                   boxprops=dict(facecolor=col, alpha=0.7),
                   medianprops=dict(color='black', linewidth=2),
                   whiskerprops=dict(linewidth=1.5),
                   capprops=dict(linewidth=1.5),
                   flierprops=dict(marker='o', markersize=3, alpha=0.4))
        ax.scatter([i + np.random.uniform(-0.15, 0.15) for _ in vals],
                   vals.values, alpha=0.3, s=15, color='gray')

    y_top2 = max(sp_pgp.max(), tr_pgp.max()) + 0.1
    ax.plot([0, 1], [y_top2, y_top2], color='black', lw=1)
    ax.text(0.5, y_top2 + 0.02, sig_label(p_pgp), ha='center', fontsize=12, fontweight='bold')
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['Space Flight', 'Terrestrial Soil'], fontsize=11)
    ax.set_ylabel("PGP Functional Specificity (Shannon H')", fontsize=11)
    ax.set_title("PGP Functional Specificity Index\n(H' across 5 PGP categories)", fontsize=11, fontweight='bold')
    ax.set_ylim(bottom=0)
    ax.axhline(np.log(5), color='gray', lw=1, linestyle='--', alpha=0.6,
               label=f"Max H'=ln(5)={np.log(5):.2f}")
    ax.legend(fontsize=9)
    ax.grid(axis='y', alpha=0.3)
    ax.text(0.98, 0.02, f'p = {p_pgp:.4f}', transform=ax.transAxes,
            ha='right', va='bottom', fontsize=9, color='gray')

    # ── D: PGP category stacked bar ──────────────────────────────────────────
    ax = axes[1, 1]
    cats = list(PGP_CATEGORIES.keys())
    cat_labels = ['N Fixation', 'P Solubilization', 'IAA/Phytohormone',
                  'Biocontrol', 'ACC-Deaminase']
    cat_colors = ['#2196F3', '#4CAF50', '#FF9800', '#F44336', '#9C27B0']

    sp_means_cat = [pgp_filt.loc[pgp_filt.index.intersection(sp_idx), c].mean()
                    for c in cats]
    tr_means_cat = [pgp_filt.loc[pgp_filt.index.intersection(tr_idx), c].mean()
                    for c in cats]

    x    = np.array([0, 1])
    bottoms = np.zeros(2)
    for j, (cat, col, lbl) in enumerate(zip(cats, cat_colors, cat_labels)):
        vals_bar = np.array([sp_means_cat[j], tr_means_cat[j]])
        ax.bar(x, vals_bar, bottom=bottoms, color=col, alpha=0.8,
               label=lbl, width=0.5)
        bottoms += vals_bar

    ax.set_xticks([0, 1])
    ax.set_xticklabels(['Space Flight', 'Terrestrial Soil'], fontsize=11)
    ax.set_ylabel('Mean Relative Abundance (%)', fontsize=11)
    ax.set_title('PGP Category Composition\n(mean relative abundance per group)',
                 fontsize=11, fontweight='bold')
    ax.legend(fontsize=9, bbox_to_anchor=(1.01, 1), loc='upper left')
    ax.grid(axis='y', alpha=0.3)

    fig.suptitle(
        'N-Cycle Completeness and PGP Functional Specificity\nSpace Flight vs. Terrestrial Soil',
        fontsize=13, fontweight='bold'
    )

    for ext in ['png', 'pdf']:
        out = OUTPUT_DIR / f'Fig2B_NCycle_PGP.{ext}'
        fig.savefig(out, dpi=300 if ext == 'png' else None, bbox_inches='tight')
        print(f'  Saved: {out.name}')

    cb_dir = config.CURRENT_BIOLOGY_DIR
    cb_dir.mkdir(exist_ok=True)
    fig.savefig(cb_dir / 'Fig1c_NCycle_Heatmap.png', dpi=300, bbox_inches='tight')
    fig.savefig(cb_dir / 'Fig1c_NCycle_Heatmap.pdf', bbox_inches='tight')
    print(f'  Saved: {cb_dir / "Fig1c_NCycle_Heatmap.png"}')

    plt.close()


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    print('=' * 60)
    print('N-Cycle Completeness + PGP Functional Specificity')
    print('=' * 60)

    print('\n[1] Parsing FAPROTAX report...')
    func2genera = parse_faprotax_report(OUTPUT_DIR / 'FAPROTAX_report.txt')
    for step in NCYCLE_STEPS:
        n = len(func2genera.get(step, set()))
        print(f'  {step}: {n} genera')

    print('\n[2] Building genus-level table...')
    genus_rel = aggregate_to_genus(config.FEATURE_TABLE_CLEAN, config.TAXONOMY_FILE)
    print(f'  {genus_rel.shape[0]} genera × {genus_rel.shape[1]} samples')

    print('\n[3] Loading metadata...')
    meta = pd.read_csv(config.INTEGRATED_METADATA, sep='\t', index_col=0)
    meta = meta[meta['Study_Group'].isin(['Space_Flight', 'Terrestrial_Soil'])]
    # Root-zone filter: restrict Space_Flight to Plant Roots + Rhizosphere
    ROOT_ZONE_TYPES = {'Plant Roots', 'Rhizosphere'}
    meta = meta[
        (meta['Study_Group'] == 'Terrestrial_Soil') |
        ((meta['Study_Group'] == 'Space_Flight') & meta['sample_type'].isin(ROOT_ZONE_TYPES))
    ]
    genus_rel = genus_rel[genus_rel.columns.intersection(meta.index)]
    n_sp = (meta['Study_Group'] == 'Space_Flight').sum()
    n_tr = (meta['Study_Group'] == 'Terrestrial_Soil').sum()
    print(f'  {len(meta)} samples (Space={n_sp}, Terrestrial={n_tr})')

    print('\n[4] Computing N-cycle completeness...')
    ncycle_df = compute_ncycle_completeness(genus_rel, func2genera, meta)
    ncycle_df = ncycle_df.loc[ncycle_df.index.intersection(meta.index)]

    sp_idx = meta.index[meta['Study_Group'] == 'Space_Flight']
    tr_idx = meta.index[meta['Study_Group'] == 'Terrestrial_Soil']
    sp_comp = ncycle_df.loc[ncycle_df.index.intersection(sp_idx), 'Completeness']
    tr_comp = ncycle_df.loc[ncycle_df.index.intersection(tr_idx), 'Completeness']
    stat_nc, p_nc = mwu_compare(sp_comp, tr_comp)
    print(f'  Space Completeness:       {sp_comp.mean():.3f} ± {sp_comp.std():.3f}')
    print(f'  Terrestrial Completeness: {tr_comp.mean():.3f} ± {tr_comp.std():.3f}')
    print(f'  MWU p = {p_nc:.6f} ({sig_label(p_nc)})')

    print('\n  N-cycle step presence (fraction of samples):')
    for step in NCYCLE_STEPS:
        sp_pres = ncycle_df.loc[ncycle_df.index.intersection(sp_idx), step].mean()
        tr_pres = ncycle_df.loc[ncycle_df.index.intersection(tr_idx), step].mean()
        print(f'  {NCYCLE_LABELS[step]:25s}: Space={sp_pres:.2f}, Terrestrial={tr_pres:.2f}')

    print('\n[5] Computing PGP functional specificity...')
    pgp_df = compute_pgp_specificity(genus_rel, meta)
    pgp_df = pgp_df.loc[pgp_df.index.intersection(meta.index)]

    sp_pgp = pgp_df.loc[pgp_df.index.intersection(sp_idx), 'Specificity_H']
    tr_pgp = pgp_df.loc[pgp_df.index.intersection(tr_idx), 'Specificity_H']
    stat_pg, p_pg = mwu_compare(sp_pgp, tr_pgp)
    print(f"  Space H':       {sp_pgp.mean():.3f} ± {sp_pgp.std():.3f} (max={np.log(5):.2f})")
    print(f"  Terrestrial H': {tr_pgp.mean():.3f} ± {tr_pgp.std():.3f}")
    print(f"  MWU p = {p_pg:.6f} ({sig_label(p_pg)})")

    print('\n  PGP category abundances (mean %):')
    for cat in PGP_CATEGORIES:
        sp_m = pgp_df.loc[pgp_df.index.intersection(sp_idx), cat].mean()
        tr_m = pgp_df.loc[pgp_df.index.intersection(tr_idx), cat].mean()
        print(f'  {cat:25s}: Space={sp_m:.2f}%, Terrestrial={tr_m:.2f}%')

    print('\n[6] Saving results...')
    ncycle_out = ncycle_df.copy()
    ncycle_out.to_csv(OUTPUT_DIR / 'Fig2B_NCycle_Completeness_Results.csv')
    print('  Saved: Fig2B_NCycle_Completeness_Results.csv')

    nc_summary = pd.DataFrame({
        'Group':             ['Space_Flight', 'Terrestrial_Soil'],
        'Completeness_Mean': [round(sp_comp.mean(), 4), round(tr_comp.mean(), 4)],
        'Completeness_SD':   [round(sp_comp.std(), 4),  round(tr_comp.std(), 4)],
        'MWU_p':             [round(p_nc, 6), None],
        'sig':               [sig_label(p_nc), None],
    })
    nc_summary.to_csv(OUTPUT_DIR / 'Fig2B_NCycle_Completeness_Summary.csv', index=False)
    print('  Saved: Fig2B_NCycle_Completeness_Summary.csv')

    pgp_df.to_csv(OUTPUT_DIR / 'Fig2B_PGP_Specificity_Results.csv')
    print('  Saved: Fig2B_PGP_Specificity_Results.csv')

    pgp_summary = pd.DataFrame({
        'Group':         ['Space_Flight', 'Terrestrial_Soil'],
        'H_Mean':        [round(sp_pgp.mean(), 4), round(tr_pgp.mean(), 4)],
        'H_SD':          [round(sp_pgp.std(), 4),  round(tr_pgp.std(), 4)],
        'MWU_p':         [round(p_pg, 6), None],
        'sig':           [sig_label(p_pg), None],
    })
    pgp_summary.to_csv(OUTPUT_DIR / 'Fig2B_PGP_Specificity_Summary.csv', index=False)
    print('  Saved: Fig2B_PGP_Specificity_Summary.csv')

    print('\n[7] Plotting...')
    plot_results(ncycle_df, pgp_df, meta)

    print('\nDone.')
