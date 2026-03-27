"""
Functional Redundancy Analysis
================================
Quantifies how many distinct genera contribute to each FAPROTAX functional
group per sample, comparing Space Flight vs. Terrestrial Soil.

Ecological rationale:
  Functional redundancy = number of taxa that can perform the same function.
  Low redundancy → a single taxon loss can abolish an entire ecosystem service.
  This provides a structural (mechanistic) explanation for why functional
  collapse (observed in FAPROTAX Log2FC) is irreversible under spaceflight.

Method:
  1. Parse FAPROTAX_report.txt → function → set of genera assigned
  2. Build genus-level relative abundance table (all samples combined)
  3. For each sample × function: count genera with abundance > 0
  4. Compare Space Flight vs. Terrestrial Soil (Mann-Whitney U)

Prerequisites:
  - version-2_integrated/FAPROTAX_report.txt
  - version-2_integrated/exported_table_clean/feature-table.tsv
  - version-2_integrated/exported_taxonomy/taxonomy.tsv
  - version-2_integrated/integrated_metadata.tsv

Run:
  python 22_functional_redundancy.py

Output (version-2_integrated/):
  - Functional_Redundancy_Results.csv   : per-sample redundancy per function
  - Functional_Redundancy_Summary.csv   : group means, SD, p-values
  - Main_Fig_FunctionalRedundancy.png/.pdf
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
import seaborn as sns

OUTPUT_DIR = config.VERSION_DIR

# Functions to analyse (biologically informative subset)
FOCUS_FUNCTIONS = [
    'nitrification',
    'nitrogen_fixation',
    'denitrification',
    'nitrate_reduction',
    'fermentation',
    'aerobic_chemoheterotrophy',
    'ureolysis',
    'cellulolysis',
    'methylotrophy',
    'sulfate_respiration',
]

FUNC_LABELS = {
    'nitrification':           'Nitrification',
    'nitrogen_fixation':       'N₂ Fixation',
    'denitrification':         'Denitrification',
    'nitrate_reduction':       'Nitrate Reduction',
    'fermentation':            'Fermentation',
    'aerobic_chemoheterotrophy': 'Aerobic Chemoheterotrophy',
    'ureolysis':               'Ureolysis',
    'cellulolysis':            'Cellulolysis',
    'methylotrophy':           'Methylotrophy',
    'sulfate_respiration':     'Sulfate Respiration',
}

COLORS = {'Space_Flight': '#3498db', 'Terrestrial_Soil': '#2ecc71'}


# ─────────────────────────────────────────────────────────────────────────────
# 1. Parse FAPROTAX report → function → set of genera
# ─────────────────────────────────────────────────────────────────────────────

def parse_faprotax_report(report_path):
    """Return dict: function_name → set of genus names from lineage strings."""
    func2genera = {}
    current_func = None

    with open(report_path) as fh:
        for line in fh:
            line = line.rstrip()
            # Section header: "# nitrification (17 records):"
            m = re.match(r'^# (\w+) \(\d+ records\):', line)
            if m:
                current_func = m.group(1)
                func2genera[current_func] = set()
                continue
            # Taxon line (4-space indent): lineage string
            if current_func and line.startswith('    '):
                genus = None
                for part in line.split(';'):
                    part = part.strip()
                    if part.startswith('g__'):
                        name = part.split('__', 1)[-1].strip()
                        # skip empty / uninformative genus-level assignments
                        if name and name not in ('uncultured', 'metagenome',
                                                  'unidentified', ''):
                            genus = name
                if genus:
                    func2genera[current_func].add(genus)

    return func2genera


# ─────────────────────────────────────────────────────────────────────────────
# 2. Build genus-level feature table (combined, all samples)
# ─────────────────────────────────────────────────────────────────────────────

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
    # Relative abundance (%) per sample
    return genus_table.div(genus_table.sum(axis=0), axis=1) * 100


# ─────────────────────────────────────────────────────────────────────────────
# 3. Compute functional redundancy
# ─────────────────────────────────────────────────────────────────────────────

def compute_redundancy(genus_rel, func2genera, functions):
    """
    Returns DataFrame: index=samples, columns=functions
    Value = number of genera contributing to each function (abundance > 0).
    """
    records = {}
    for func in functions:
        assigned_genera = func2genera.get(func, set())
        present = genus_rel.index.intersection(assigned_genera)
        if len(present) == 0:
            records[func] = pd.Series(0, index=genus_rel.columns)
        else:
            # Count genera with abundance > 0 per sample
            records[func] = (genus_rel.loc[present] > 0).sum(axis=0)

    return pd.DataFrame(records)   # (samples × functions)


# ─────────────────────────────────────────────────────────────────────────────
# 4. Statistical comparison
# ─────────────────────────────────────────────────────────────────────────────

def compare_groups(redund_df, meta, functions):
    sp_idx = meta.index[meta['Study_Group'] == 'Space_Flight']
    tr_idx = meta.index[meta['Study_Group'] == 'Terrestrial_Soil']

    rows = []
    p_vals = []
    for func in functions:
        sp_vals = redund_df.loc[redund_df.index.intersection(sp_idx), func].dropna()
        tr_vals = redund_df.loc[redund_df.index.intersection(tr_idx), func].dropna()
        if len(sp_vals) == 0 or len(tr_vals) == 0:
            continue
        stat, p = mannwhitneyu(sp_vals, tr_vals, alternative='two-sided')
        rows.append({
            'Function':    func,
            'Space_Mean':  round(sp_vals.mean(), 4),
            'Space_SD':    round(sp_vals.std(), 4),
            'Terr_Mean':   round(tr_vals.mean(), 4),
            'Terr_SD':     round(tr_vals.std(), 4),
            'MWU_stat':    round(stat, 1),
            'p_raw':       round(p, 6),
        })
        p_vals.append(p)

    df = pd.DataFrame(rows)
    if len(p_vals) > 0:
        _, q_vals, _, _ = multipletests(p_vals, method='fdr_bh')
        df['q_BH'] = np.round(q_vals, 6)
    df['sig'] = df['q_BH'].apply(
        lambda q: '***' if q < 0.001 else ('**' if q < 0.01
                   else ('*' if q < 0.05 else 'n.s.')))
    return df.sort_values('p_raw')


# ─────────────────────────────────────────────────────────────────────────────
# 5. Plot
# ─────────────────────────────────────────────────────────────────────────────

def plot_redundancy(redund_df, meta, summary_df, functions):
    sp_idx = meta.index[meta['Study_Group'] == 'Space_Flight']
    tr_idx = meta.index[meta['Study_Group'] == 'Terrestrial_Soil']

    # ── Panel A: grouped bar chart (mean ± SE) for focus functions ──────────
    func_labels_ordered = [FUNC_LABELS[f] for f in functions
                           if f in FUNC_LABELS]
    func_keys_ordered   = [f for f in functions if f in FUNC_LABELS]

    sp_means = [redund_df.loc[redund_df.index.intersection(sp_idx), f].mean()
                for f in func_keys_ordered]
    tr_means = [redund_df.loc[redund_df.index.intersection(tr_idx), f].mean()
                for f in func_keys_ordered]
    sp_se    = [redund_df.loc[redund_df.index.intersection(sp_idx), f].sem()
                for f in func_keys_ordered]
    tr_se    = [redund_df.loc[redund_df.index.intersection(tr_idx), f].sem()
                for f in func_keys_ordered]

    x     = np.arange(len(func_keys_ordered))
    width = 0.35

    fig, axes = plt.subplots(1, 2, figsize=(16, 6), constrained_layout=True)

    # Bar chart
    ax = axes[0]
    bars_sp = ax.bar(x - width/2, sp_means, width, yerr=sp_se,
                     color=COLORS['Space_Flight'], alpha=0.8,
                     label='Space Flight', capsize=4, error_kw={'elinewidth': 1.2})
    bars_tr = ax.bar(x + width/2, tr_means, width, yerr=tr_se,
                     color=COLORS['Terrestrial_Soil'], alpha=0.8,
                     label='Terrestrial Soil', capsize=4, error_kw={'elinewidth': 1.2})

    # Significance stars from summary_df
    for i, func in enumerate(func_keys_ordered):
        row = summary_df[summary_df['Function'] == func]
        if len(row) == 0:
            continue
        sig = row.iloc[0]['sig']
        if sig != 'n.s.':
            y_top = max(sp_means[i] + sp_se[i], tr_means[i] + tr_se[i]) + 0.15
            ax.text(x[i], y_top, sig, ha='center', va='bottom', fontsize=11,
                    fontweight='bold', color='black')

    ax.set_xticks(x)
    ax.set_xticklabels(func_labels_ordered, rotation=35, ha='right', fontsize=10)
    ax.set_ylabel('Number of Contributing Genera\n(mean ± SE)', fontsize=11)
    ax.set_title('Functional Redundancy per Group', fontsize=12, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(axis='y', alpha=0.3)
    ax.set_ylim(bottom=0)

    # ── Panel B: heatmap of mean redundancy ─────────────────────────────────
    ax = axes[1]
    hm_data = pd.DataFrame({
        'Space Flight':    [redund_df.loc[redund_df.index.intersection(sp_idx), f].mean()
                            for f in func_keys_ordered],
        'Terrestrial Soil': [redund_df.loc[redund_df.index.intersection(tr_idx), f].mean()
                              for f in func_keys_ordered],
    }, index=func_labels_ordered)

    sns.heatmap(hm_data, ax=ax, annot=True, fmt='.1f',
                cmap='YlOrRd', linewidths=0.5, linecolor='white',
                cbar_kws={'label': 'Mean # Contributing Genera'})
    ax.set_title('Functional Redundancy Heatmap\n(mean contributing genera per sample)',
                 fontsize=12, fontweight='bold')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0, fontsize=11)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=10)

    fig.suptitle(
        'Functional Redundancy: Space Flight vs. Terrestrial Soil\n'
        '(Number of distinct genera per FAPROTAX functional group)',
        fontsize=13
    )

    for ext in ['png', 'pdf']:
        out = OUTPUT_DIR / f'Fig2C_FunctionalRedundancy.{ext}'
        fig.savefig(out, dpi=300 if ext == 'png' else None, bbox_inches='tight')
        print(f'  Saved: {out.name}')

    cb_dir = config.CURRENT_BIOLOGY_DIR
    cb_dir.mkdir(exist_ok=True)
    fig.savefig(cb_dir / 'Fig1d_FunctionalRedundancy.png', dpi=300, bbox_inches='tight')
    fig.savefig(cb_dir / 'Fig1d_FunctionalRedundancy.pdf', bbox_inches='tight')
    print(f'  Saved: {cb_dir / "Fig1d_FunctionalRedundancy.png"}')

    plt.close()


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    print('=' * 60)
    print('Functional Redundancy Analysis')
    print('=' * 60)

    print('\n[1] Parsing FAPROTAX report...')
    report_path = OUTPUT_DIR / 'FAPROTAX_report.txt'
    func2genera = parse_faprotax_report(report_path)
    print(f'  Functions parsed: {len(func2genera)}')
    for f in FOCUS_FUNCTIONS:
        n = len(func2genera.get(f, set()))
        print(f'  {f}: {n} genera assigned')

    print('\n[2] Building genus-level feature table...')
    genus_rel = aggregate_to_genus(config.FEATURE_TABLE_CLEAN, config.TAXONOMY_FILE)
    print(f'  {genus_rel.shape[0]} genera × {genus_rel.shape[1]} samples')

    print('\n[3] Loading metadata and filtering...')
    meta = pd.read_csv(config.INTEGRATED_METADATA, sep='\t', index_col=0)
    meta = meta[meta['Study_Group'].isin(['Space_Flight', 'Terrestrial_Soil'])]
    genus_rel = genus_rel[genus_rel.columns.intersection(meta.index)]
    print(f'  {len(meta)} samples retained '
          f'(Space={sum(meta.Study_Group=="Space_Flight")}, '
          f'Terrestrial={sum(meta.Study_Group=="Terrestrial_Soil")})')

    print('\n[4] Computing functional redundancy...')
    redund_df = compute_redundancy(genus_rel, func2genera, FOCUS_FUNCTIONS)
    redund_df = redund_df.loc[redund_df.index.intersection(meta.index)]
    print(redund_df.describe().round(2))

    print('\n[5] Statistical comparison...')
    summary_df = compare_groups(redund_df, meta, FOCUS_FUNCTIONS)
    print(summary_df.to_string(index=False))

    print('\n[6] Saving results...')
    # Per-sample results with group labels
    out_df = redund_df.copy()
    out_df.insert(0, 'Study_Group', meta.loc[out_df.index, 'Study_Group'])
    out_df.to_csv(OUTPUT_DIR / 'Fig2C_Functional_Redundancy_Results.csv')
    summary_df.to_csv(OUTPUT_DIR / 'Fig2C_Functional_Redundancy_Summary.csv', index=False)
    print('  Saved: Fig2C_Functional_Redundancy_Results.csv')
    print('  Saved: Fig2C_Functional_Redundancy_Summary.csv')

    # Additional: total genera with ANY function assignment per sample
    all_assigned = set().union(*func2genera.values())
    present_in_table = genus_rel.index.intersection(all_assigned)
    n_functional_genera = (genus_rel.loc[present_in_table] > 0).sum(axis=0)
    # Compare
    sp_func = n_functional_genera[meta.index[meta.Study_Group == 'Space_Flight']]
    tr_func = n_functional_genera[meta.index[meta.Study_Group == 'Terrestrial_Soil']]
    stat, p = mannwhitneyu(sp_func, tr_func, alternative='two-sided')
    print(f'\n  Total functional genera per sample:')
    print(f'    Space       : {sp_func.mean():.1f} ± {sp_func.std():.1f}')
    print(f'    Terrestrial : {tr_func.mean():.1f} ± {tr_func.std():.1f}')
    print(f'    MWU p = {p:.4f}')

    print('\n[7] Plotting...')
    plot_redundancy(redund_df, meta, summary_df, FOCUS_FUNCTIONS)

    print('\nDone.')
