"""
FAPROTAX Functional Trait Analysis
===================================
Maps ASVs to experimentally-validated functional traits (nitrogen fixation,
phosphate solubilization, plant growth promotion, etc.) using the FAPROTAX
database (Louca et al. 2016, Science).

Unlike PICRUSt2, FAPROTAX is based on curated literature-validated functional
assignments, making results more defensible as mechanistic evidence.

Input:
  - version-2_integrated/exported_table_clean/feature-table.tsv  (ASV counts)
  - version-2_integrated/exported_taxonomy/taxonomy.tsv           (SILVA taxonomy)
  - version-2_integrated/integrated_metadata.tsv                  (group labels)

Output (version-2_integrated/):
  - FAPROTAX_functional_table.tsv    : functions × samples matrix
  - FAPROTAX_group_comparison.csv    : mean relative abundance per group per function
  - FAPROTAX_key_functions.csv       : plant-relevant functions only (fold-change ranked)
  - Main_Fig_FAPROTAX.png / .pdf     : bar plot of key plant-relevant functions

Key functions tracked:
  nitrogen fixation, nitrification, denitrification,
  phosphate solubilization, plant growth promotion (PGPR),
  chemoheterotrophy, aerobic chemoheterotrophy,
  fermentation, sulfate reduction, methanogenesis
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import config

import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

FAPROTAX_DIR   = Path('/home/laugh/pepper/FAPROTAX_1.2.12')
COLLAPSE_SCRIPT = FAPROTAX_DIR / 'collapse_table.py'
FAPROTAX_DB    = FAPROTAX_DIR / 'FAPROTAX.txt'
OUTPUT_DIR     = config.VERSION_DIR

# Plant-relevant functional groups to highlight
PLANT_RELEVANT = [
    'nitrogen_fixation',
    'nitrification',
    'nitrate_reduction',
    'denitrification',
    'phosphorus_cycling',
    'plant_pathogen',
    'human_pathogens_all',
    'chemoheterotrophy',
    'aerobic_chemoheterotrophy',
    'fermentation',
    'sulfate_respiration',
    'dark_sulfide_oxidation',
    'aromatic_compound_degradation',
    'methylotrophy',
]

# N cycle functions for detailed subplot
N_CYCLE_FUNCTIONS = [
    'nitrification',
    'aerobic_ammonia_oxidation',
    'aerobic_nitrite_oxidation',
    'denitrification',
    'nitrogen_fixation',
    'nitrate_reduction',
]

# Control functions (stable reference) for contrast
CONTROL_FUNCTIONS = [
    'aerobic_chemoheterotrophy',
    'fermentation',
    'aromatic_compound_degradation',
    'cellulolysis',
    'chitinolysis',
]


# ─────────────────────────────────────────────
# 1. Prepare BIOM-style input for collapse_table.py
# ─────────────────────────────────────────────

def prepare_input_table(meta):
    """
    Create FAPROTAX-compatible TSV where row names = full SILVA taxonomy strings.
    ASVs sharing the same taxonomy string are summed (collapsed to taxon level).
    collapse_table.py then maps taxonomy strings to functional groups.
    """
    print('  Loading feature table and taxonomy...')
    table = pd.read_csv(config.FEATURE_TABLE_CLEAN, sep='\t', skiprows=1, index_col=0)

    # Filter to Space_Flight and Terrestrial_Soil only
    keep = meta.index[meta['Study_Group'].isin(['Space_Flight', 'Terrestrial_Soil'])].tolist()
    table = table[[c for c in table.columns if c in keep]]

    tax = pd.read_csv(config.TAXONOMY_FILE, sep='\t', index_col=0)
    shared = table.index.intersection(tax.index)
    table  = table.loc[shared]

    # Build full taxonomy string in SILVA format: k__; p__; c__; o__; f__; g__; s__
    def clean_taxon(t):
        if pd.isna(t):
            return 'Unassigned'
        return '; '.join(p.strip() for p in t.split(';'))

    tax['taxon_clean'] = tax['Taxon'].apply(clean_taxon)

    # Assign taxonomy to each ASV row and collapse (sum) by taxonomy string
    table = table.copy()
    table['taxon'] = tax.loc[shared, 'taxon_clean'].values
    collapsed = table.groupby('taxon').sum()   # row index = taxonomy string

    out_path = OUTPUT_DIR / 'faprotax_input.tsv'
    # FAPROTAX expects: first row = header, first column = taxonomy strings.
    # Use plain name (no #) so FAPROTAX reads first row as column headers.
    collapsed.index.name = 'OTU_ID'
    collapsed.to_csv(out_path, sep='\t')

    print(f'    Input table: {collapsed.shape[0]} unique taxa × {collapsed.shape[1]} samples')
    return out_path, collapsed.columns.tolist()


# ─────────────────────────────────────────────
# 2. Run collapse_table.py
# ─────────────────────────────────────────────

def run_faprotax(input_path):
    func_table_path    = OUTPUT_DIR / 'FAPROTAX_functional_table.tsv'
    report_path        = OUTPUT_DIR / 'FAPROTAX_report.txt'
    collapsed_otus_path = OUTPUT_DIR / 'FAPROTAX_collapsed_otus.tsv'

    cmd = [
        sys.executable, str(COLLAPSE_SCRIPT),
        '-i', str(input_path),
        '-o', str(func_table_path),
        '-g', str(FAPROTAX_DB),
        '-r', str(report_path),
        '--row_names_are_in_column', 'OTU_ID',
        '-f',   # force overwrite
        '-v',
    ]
    print('  Running FAPROTAX collapse_table.py ...')
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print('STDERR:', result.stderr[-2000:])
        raise RuntimeError('collapse_table.py failed')
    print(result.stdout[-500:] if result.stdout else '  (no stdout)')
    return func_table_path


# ─────────────────────────────────────────────
# 3. Compare groups
# ─────────────────────────────────────────────

def compare_groups(func_table_path, meta):
    print('  Computing group means...')
    func = pd.read_csv(func_table_path, sep='\t', index_col=0, comment='#')

    # Relative abundance (% within each sample)
    func_rel = func.div(func.sum(axis=0), axis=1) * 100

    sp_samples  = meta.index[meta['Study_Group'] == 'Space_Flight'].tolist()
    tr_samples  = meta.index[meta['Study_Group'] == 'Terrestrial_Soil'].tolist()

    sp_cols = [c for c in func_rel.columns if c in sp_samples]
    tr_cols = [c for c in func_rel.columns if c in tr_samples]

    result = pd.DataFrame({
        'Space_Flight_mean':      func_rel[sp_cols].mean(axis=1),
        'Terrestrial_Soil_mean':  func_rel[tr_cols].mean(axis=1),
    })
    epsilon = 1e-6
    result['Log2FC'] = np.log2(
        (result['Space_Flight_mean'] + epsilon) /
        (result['Terrestrial_Soil_mean'] + epsilon)
    )
    result['n_Space_samples']  = len(sp_cols)
    result['n_Terr_samples']   = len(tr_cols)

    out_path = OUTPUT_DIR / 'FAPROTAX_group_comparison.csv'
    result.sort_values('Log2FC', ascending=False).to_csv(out_path)
    print(f'    Saved: {out_path.name}  ({len(result)} functions)')
    return result


def extract_key_functions(result):
    present = [f for f in PLANT_RELEVANT if f in result.index]
    key = result.loc[present].sort_values('Log2FC', ascending=False)

    out_path = OUTPUT_DIR / 'FAPROTAX_key_functions.csv'
    key.to_csv(out_path)
    print(f'    Key functions found: {len(present)} / {len(PLANT_RELEVANT)} requested')
    print(key[['Space_Flight_mean', 'Terrestrial_Soil_mean', 'Log2FC']].to_string())
    return key


# ─────────────────────────────────────────────
# 4. Plot
# ─────────────────────────────────────────────

def plot_key_functions(key_df):
    if key_df.empty:
        print('  No key functions to plot.')
        return

    fig, ax = plt.subplots(figsize=(10, max(5, len(key_df) * 0.55)))

    labels = [f.replace('_', ' ').title() for f in key_df.index]
    x = np.arange(len(labels))
    width = 0.38

    sp_vals = key_df['Space_Flight_mean'].values
    tr_vals = key_df['Terrestrial_Soil_mean'].values

    ax.barh(x + width/2, sp_vals, width, color='#2980b9', label='Space Flight')
    ax.barh(x - width/2, tr_vals, width, color='#27ae60', label='Terrestrial Soil')

    ax.set_yticks(x)
    ax.set_yticklabels(labels, fontsize=10)
    ax.set_xlabel('Mean Relative Abundance (%)', fontsize=11)
    ax.set_title(
        'Functional Trait Comparison (FAPROTAX)\nSpace Flight vs. Terrestrial Soil',
        fontsize=13, fontweight='bold'
    )
    ax.legend(fontsize=10)
    ax.invert_yaxis()
    plt.tight_layout()

    png_path = OUTPUT_DIR / 'Fig2A_FAPROTAX_barplot.png'
    pdf_path = OUTPUT_DIR / 'Fig2A_FAPROTAX_barplot.pdf'
    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    plt.savefig(pdf_path, bbox_inches='tight')
    plt.close()
    print(f'    Saved: {png_path.name}')
    print(f'    Saved: {pdf_path.name}')


def plot_ncycle_detail(result):
    """Two-panel figure: N cycle functions (left) + control functions (right)."""
    ncycle  = [f for f in N_CYCLE_FUNCTIONS  if f in result.index]
    control = [f for f in CONTROL_FUNCTIONS  if f in result.index]
    if not ncycle:
        print('  No N cycle functions found — skipping detail plot.')
        return

    fig, axes = plt.subplots(1, 2, figsize=(14, max(5, max(len(ncycle), len(control)) * 0.65)))

    for ax, funcs, title in zip(
        axes,
        [ncycle, control],
        ['Nitrogen Cycling Functions', 'Reference Functions'],
    ):
        if not funcs:
            ax.set_visible(False)
            continue
        df   = result.loc[funcs]
        labels = [f.replace('_', ' ').title() for f in funcs]
        x    = np.arange(len(labels))
        w    = 0.38

        ax.barh(x + w/2, df['Space_Flight_mean'].values,  w, color='#2980b9', label='Space Flight')
        ax.barh(x - w/2, df['Terrestrial_Soil_mean'].values, w, color='#27ae60', label='Terrestrial Soil')
        ax.set_yticks(x)
        ax.set_yticklabels(labels, fontsize=10)
        ax.set_xlabel('Mean Relative Abundance (%)', fontsize=10)
        ax.set_title(title, fontsize=11, fontweight='bold')
        ax.invert_yaxis()
        ax.legend(fontsize=9)

    fig.suptitle('FAPROTAX Functional Trait Comparison\nSpace Flight vs. Terrestrial Soil',
                 fontsize=13, fontweight='bold', y=1.01)
    plt.tight_layout()

    png_path = OUTPUT_DIR / 'Fig2A_FAPROTAX_Ncycle.png'
    pdf_path = OUTPUT_DIR / 'Fig2A_FAPROTAX_Ncycle.pdf'
    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    plt.savefig(pdf_path, bbox_inches='tight')
    plt.close()
    print(f'    Saved: {png_path.name}')
    print(f'    Saved: {pdf_path.name}')


def plot_log2fc(result):
    """Horizontal log2FC barplot for N cycle + control functions."""
    funcs = [f for f in N_CYCLE_FUNCTIONS + CONTROL_FUNCTIONS if f in result.index]
    if not funcs:
        print('  No functions for log2FC plot.')
        return

    df = result.loc[funcs, 'Log2FC'].sort_values()
    labels = [f.replace('_', ' ').title() for f in df.index]
    colors = ['#c0392b' if v < 0 else '#2980b9' for v in df.values]

    # Mark N cycle vs control with linestyle
    is_ncycle = [f in N_CYCLE_FUNCTIONS for f in df.index]

    fig, ax = plt.subplots(figsize=(9, max(5, len(df) * 0.55)))
    bars = ax.barh(range(len(df)), df.values, color=colors, edgecolor='white', height=0.65)

    # Bold N cycle labels
    ax.set_yticks(range(len(df)))
    ticklabels = ax.set_yticklabels(labels, fontsize=10)
    for lbl, nc in zip(ticklabels, is_ncycle):
        if nc:
            lbl.set_fontweight('bold')

    ax.axvline(0, color='black', linewidth=0.8)
    ax.set_xlabel('Log₂ Fold Change (Space Flight / Terrestrial Soil)', fontsize=11)
    ax.set_title(
        'FAPROTAX Functional Trait Log₂FC\nSpace Flight vs. Terrestrial Soil',
        fontsize=13, fontweight='bold'
    )

    # Annotate values
    for i, v in enumerate(df.values):
        ax.text(v + (0.15 if v >= 0 else -0.15), i, f'{v:+.1f}',
                va='center', ha='left' if v >= 0 else 'right', fontsize=8.5)

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#2980b9', label='Space-enriched'),
        Patch(facecolor='#c0392b', label='Space-depleted'),
    ]
    ax.legend(handles=legend_elements, fontsize=9, loc='lower right')

    # Annotation: N cycle = bold
    ax.text(0.01, 0.99, '★ Bold = N cycle function', transform=ax.transAxes,
            fontsize=8, va='top', color='#555555')

    plt.tight_layout()

    png_path = OUTPUT_DIR / 'Fig2A_FAPROTAX_Log2FC.png'
    pdf_path = OUTPUT_DIR / 'Fig2A_FAPROTAX_Log2FC.pdf'
    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    plt.savefig(pdf_path, bbox_inches='tight')

    cb_dir = config.CURRENT_BIOLOGY_DIR
    cb_dir.mkdir(exist_ok=True)
    plt.savefig(cb_dir / 'Fig1b_FAPROTAX_Log2FC.png', dpi=300, bbox_inches='tight')
    plt.savefig(cb_dir / 'Fig1b_FAPROTAX_Log2FC.pdf', bbox_inches='tight')

    plt.close()
    print(f'    Saved: {png_path.name}')
    print(f'    Saved: {pdf_path.name}')
    print(f'    Saved: {cb_dir / "Fig1b_FAPROTAX_Log2FC.png"}')


# ─────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────

if __name__ == '__main__':
    print('=' * 60)
    print('FAPROTAX Functional Trait Analysis')
    print('=' * 60)

    meta = pd.read_csv(config.INTEGRATED_METADATA, sep='\t', index_col=0)

    print('\n[1] Preparing input table...')
    input_path, sample_ids = prepare_input_table(meta)

    print('\n[2] Running FAPROTAX...')
    func_table_path = run_faprotax(input_path)

    print('\n[3] Comparing Space vs. Terrestrial...')
    result = compare_groups(func_table_path, meta)

    print('\n[4] Extracting plant-relevant functions...')
    key_df = extract_key_functions(result)

    print('\n[5] Plotting...')
    plot_key_functions(key_df)
    plot_ncycle_detail(result)
    plot_log2fc(result)

    print('\nDone.')
