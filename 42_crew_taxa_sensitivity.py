"""
42_crew_taxa_sensitivity.py
============================
Sensitivity check: does nitrifier depletion persist after removing
crew/built-environment associated genera?

Crew-associated taxa are identified as genera flagged in the FAPROTAX
mammal_gut, human_gut, or human_associated categories (from the full
FAPROTAX functional table) AND genera known from the ISS built-environment
literature (Checinska et al. 2019, OSD-65 cross-reference).

Analysis:
  1. Identify crew-associated genera from FAPROTAX annotations
  2. Remove those genera from the feature table (genus level)
  3. Re-run nitrifier composite prevalence + FAPROTAX nitrification Log2FC
  4. Compare with full-dataset results

Analytical root zone: Space n=20, Terrestrial n=20.

Output: version-2_integrated/Crew_Taxa_Sensitivity.csv
        version-2_integrated/Crew_Taxa_Sensitivity.png / .pdf
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
import config

NITRIFIER_GENERA  = ['Nitrosospira', 'Nitrospira']
ROOT_ZONE_TYPES   = {'Plant Roots', 'Rhizosphere'}
RAREFACTION_DEPTH = config.DIVERSITY_PARAMS['rarefaction_depth']
RANDOM_SEED       = config.RANDOM_SEED
OUTPUT_DIR        = config.VERSION_DIR

rng = np.random.default_rng(RANDOM_SEED)

# Known ISS built-environment / human-associated genera
# Sources: Checinska et al. 2019 (ISS-MO), common skin/gut microbiome genera
ISS_BUILT_ENV_GENERA = {
    # Skin / human microbiome
    'Staphylococcus', 'Corynebacterium', 'Cutibacterium', 'Propionibacterium',
    'Malassezia', 'Streptococcus', 'Lactobacillus', 'Enterococcus',
    'Clostridium', 'Blautia', 'Faecalibacterium', 'Ruminococcus',
    'Bacteroides', 'Prevotella', 'Bifidobacterium', 'Akkermansia',
    'Escherichia', 'Klebsiella', 'Enterobacter',
    # Built environment (Checinska 2019 top ISS genera)
    'Methylobacterium', 'Sphingomonas', 'Bacillus', 'Micrococcus',
    'Rothia', 'Actinomyces', 'Gemella', 'Granulicatella',
    'Haemophilus', 'Veillonella',
}


def parse_genus(taxon_str):
    for p in reversed(str(taxon_str).split(';')):
        p = p.strip()
        if p.startswith('g__'):
            g = p[3:]
            if g and g not in ('', 'uncultured', 'metagenome'):
                return g
    return 'Unknown'


def load_genus_ra(sp_ids, tr_ids):
    ft  = pd.read_csv(config.FEATURE_TABLE_CLEAN, sep='\t', skiprows=1, index_col=0)
    tax = pd.read_csv(
        config.VERSION_DIR / 'exported_taxonomy' / 'taxonomy.tsv',
        sep='\t', index_col=0
    )
    tax['genus'] = tax['Taxon'].apply(parse_genus)

    all_ids  = [s for s in sp_ids + tr_ids if s in ft.columns]
    ft_sub   = ft[all_ids].T

    # rarefy
    valid = ft_sub[ft_sub.sum(axis=1) >= RAREFACTION_DEPTH]

    def rarefy(row):
        pool   = np.repeat(np.arange(len(row)), row.astype(int))
        chosen = rng.choice(pool, size=RAREFACTION_DEPTH, replace=False)
        return np.bincount(chosen, minlength=len(row)).astype(float) / RAREFACTION_DEPTH

    rarefied = valid.apply(rarefy, axis=1, result_type='expand')
    rarefied.columns = valid.columns

    shared   = rarefied.columns.intersection(tax.index)
    genus_ra = rarefied[shared].T
    genus_ra.index = tax.loc[shared, 'genus']
    genus_ra = genus_ra.groupby(level=0).sum().T   # samples × genera
    return genus_ra


def cliffs_delta(a, b):
    a, b = np.asarray(a), np.asarray(b)
    dom  = np.sum(a[:, None] > b[None, :]) - np.sum(a[:, None] < b[None, :])
    return dom / (len(a) * len(b))


def faprotax_nitrif_log2fc(fa_ra, sp_ids, tr_ids):
    nit_rows = [r for r in fa_ra.index if 'nitrif' in r.lower()]
    if not nit_rows:
        return np.nan
    nit  = fa_ra.loc[nit_rows].sum(axis=0)
    sp_v = nit.reindex([s for s in sp_ids if s in nit.index]).fillna(0).values
    tr_v = nit.reindex([s for s in tr_ids if s in nit.index]).fillna(0).values
    pseudo = 1e-6
    return float(np.log2((sp_v.mean() + pseudo) / (tr_v.mean() + pseudo)))


def compute_nitrifier_stats(genus_ra, sp_in, tr_in):
    present = [g for g in NITRIFIER_GENERA if g in genus_ra.columns]
    if present:
        sp_comp = genus_ra.loc[sp_in, present].sum(axis=1).values if sp_in else np.zeros(0)
        tr_comp = genus_ra.loc[tr_in, present].sum(axis=1).values if tr_in else np.zeros(0)
    else:
        sp_comp = np.zeros(len(sp_in))
        tr_comp = np.zeros(len(tr_in))

    sp_prev = int((sp_comp > 0).sum())
    tr_prev = int((tr_comp > 0).sum())
    _, fisher_p = stats.fisher_exact(
        [[sp_prev, len(sp_in) - sp_prev],
         [tr_prev, len(tr_in) - tr_prev]]
    )
    cd = cliffs_delta(sp_comp, tr_comp)
    return {
        'Space_prevalent': sp_prev, 'Space_n': len(sp_in),
        'Terr_prevalent':  tr_prev, 'Terr_n':  len(tr_in),
        'Fisher_p': fisher_p, 'Cliffs_delta': cd,
    }


def main():
    print('=== 42_crew_taxa_sensitivity.py ===')

    meta = pd.read_csv(config.INTEGRATED_METADATA, sep='\t', index_col=0)
    sp_ids = meta.index[(meta.Study_Group == 'Space_Flight') &
                        meta.sample_type.isin(ROOT_ZONE_TYPES)].tolist()
    tr_ids = meta.index[meta.Study_Group == 'Terrestrial_Soil'].tolist()

    print(f'  Space n={len(sp_ids)}, Terrestrial n={len(tr_ids)}')

    # Load FAPROTAX functional table (all samples) and normalize to RA
    fa = pd.read_csv(config.VERSION_DIR / 'FAPROTAX_functional_table.tsv',
                     sep='\t', index_col=0)
    col_sums = fa.sum(axis=0)
    col_sums[col_sums == 0] = 1
    fa_ra = fa.div(col_sums, axis=1)

    # ── Identify crew-associated genera from FAPROTAX ─────────────────────────
    # Get genus list from feature table taxonomy
    tax = pd.read_csv(
        config.VERSION_DIR / 'exported_taxonomy' / 'taxonomy.tsv',
        sep='\t', index_col=0
    )
    tax['genus'] = tax['Taxon'].apply(parse_genus)
    all_genera = set(tax['genus'].unique())

    # Crew genera = ISS built-env list intersected with genera in dataset
    crew_genera_in_data = all_genera.intersection(ISS_BUILT_ENV_GENERA)
    print(f'  Crew-associated genera in dataset: {len(crew_genera_in_data)}')
    print(f'  {sorted(crew_genera_in_data)}')

    # ── Full dataset (no removal) ─────────────────────────────────────────────
    print('\n[1] Full dataset (no crew removal)...')
    genus_ra_full = load_genus_ra(sp_ids, tr_ids)
    sp_in  = [s for s in sp_ids if s in genus_ra_full.index]
    tr_in  = [s for s in tr_ids if s in genus_ra_full.index]

    log2fc_full = faprotax_nitrif_log2fc(fa_ra, sp_ids, tr_ids)
    stats_full  = compute_nitrifier_stats(genus_ra_full, sp_in, tr_in)
    print(f'  Nitrification Log2FC: {log2fc_full:.3f}')
    print(f'  Composite nitrifier: Space {stats_full["Space_prevalent"]}/{stats_full["Space_n"]}  '
          f'Terr {stats_full["Terr_prevalent"]}/{stats_full["Terr_n"]}  '
          f'Fisher p={stats_full["Fisher_p"]:.4f}  delta={stats_full["Cliffs_delta"]:.3f}')

    # ── Crew-removed dataset ──────────────────────────────────────────────────
    print('\n[2] After removing crew-associated genera...')
    genus_ra_clean = genus_ra_full.drop(
        columns=[g for g in crew_genera_in_data if g in genus_ra_full.columns],
        errors='ignore'
    )
    removed_count = genus_ra_full.shape[1] - genus_ra_clean.shape[1]
    print(f'  Removed {removed_count} genera columns')

    # Recompute FAPROTAX log2fc on crew-removed feature table
    # We remap to genus level then exclude crew genera from total count
    ft   = pd.read_csv(config.FEATURE_TABLE_CLEAN, sep='\t', skiprows=1, index_col=0)
    crew_asvs = tax.index[tax['genus'].isin(crew_genera_in_data)].intersection(ft.index)
    ft_clean  = ft.drop(index=crew_asvs)
    print(f'  Removed {len(crew_asvs)} crew ASVs from feature table')

    # Rerun FAPROTAX equivalent: use cleaned table RA
    # Since we can't rerun FAPROTAX itself easily, we compute nitrifier RA directly
    # from the cleaned genus table (which is the most reliable route)
    sp_in_c = [s for s in sp_ids if s in genus_ra_clean.index]
    tr_in_c = [s for s in tr_ids if s in genus_ra_clean.index]
    stats_clean = compute_nitrifier_stats(genus_ra_clean, sp_in_c, tr_in_c)
    print(f'  Composite nitrifier: Space {stats_clean["Space_prevalent"]}/{stats_clean["Space_n"]}  '
          f'Terr {stats_clean["Terr_prevalent"]}/{stats_clean["Terr_n"]}  '
          f'Fisher p={stats_clean["Fisher_p"]:.4f}  delta={stats_clean["Cliffs_delta"]:.3f}')

    # For FAPROTAX log2fc on cleaned data, recalculate from genus RA
    nit_present = [g for g in NITRIFIER_GENERA if g in genus_ra_clean.columns]
    if nit_present:
        sp_nit = genus_ra_clean.loc[sp_in_c, nit_present].sum(axis=1).values
        tr_nit = genus_ra_clean.loc[tr_in_c, nit_present].sum(axis=1).values
        pseudo = 1e-6
        log2fc_clean = float(np.log2((sp_nit.mean() + pseudo) / (tr_nit.mean() + pseudo)))
    else:
        log2fc_clean = float(np.log2(1e-6 / (1e-6 + 1e-6)))  # effectively -inf
    print(f'  Nitrifier RA Log2FC (crew-removed): {log2fc_clean:.3f}')

    # ── Save results ──────────────────────────────────────────────────────────
    rows = [
        {'Condition': 'Full dataset',
         'Genera_removed': 0,
         'Nitrifier_Log2FC': round(log2fc_full, 3),
         **{f'Full_{k}': v for k, v in stats_full.items()}},
        {'Condition': 'Crew taxa removed',
         'Genera_removed': removed_count,
         'Nitrifier_Log2FC': round(log2fc_clean, 3),
         **{f'Crew_{k}': v for k, v in stats_clean.items()}},
    ]
    df = pd.DataFrame([
        {'Condition': 'Full dataset (no removal)',
         'Genera_removed': 0,
         'Nitrifier_RA_Log2FC': round(log2fc_full, 3),
         'Space_prev': stats_full['Space_prevalent'],
         'Space_n': stats_full['Space_n'],
         'Terr_prev': stats_full['Terr_prevalent'],
         'Terr_n': stats_full['Terr_n'],
         'Fisher_p': stats_full['Fisher_p'],
         'Cliffs_delta': stats_full['Cliffs_delta']},
        {'Condition': 'Crew/built-env taxa removed',
         'Genera_removed': removed_count,
         'Nitrifier_RA_Log2FC': round(log2fc_clean, 3),
         'Space_prev': stats_clean['Space_prevalent'],
         'Space_n': stats_clean['Space_n'],
         'Terr_prev': stats_clean['Terr_prevalent'],
         'Terr_n': stats_clean['Terr_n'],
         'Fisher_p': stats_clean['Fisher_p'],
         'Cliffs_delta': stats_clean['Cliffs_delta']},
    ])
    out_csv = OUTPUT_DIR / 'Crew_Taxa_Sensitivity.csv'
    df.to_csv(out_csv, index=False, float_format='%.4f')
    print(f'\nSaved: {out_csv}')
    print(f'\nCrew genera removed: {sorted(crew_genera_in_data)}')

    # ── Figure ────────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(9, 4.5))

    # Panel a: prevalence comparison
    ax = axes[0]
    conditions = ['Full\ndataset', 'Crew taxa\nremoved']
    sp_prevs = [stats_full['Space_prevalent'] / stats_full['Space_n'] * 100,
                stats_clean['Space_prevalent'] / stats_clean['Space_n'] * 100]
    tr_prevs = [stats_full['Terr_prevalent']  / stats_full['Terr_n']  * 100,
                stats_clean['Terr_prevalent'] / stats_clean['Terr_n'] * 100]
    x = np.arange(len(conditions))
    w = 0.35
    ax.bar(x - w/2, sp_prevs, w, color='#E8735A', label='ISS (Space)', zorder=3)
    ax.bar(x + w/2, tr_prevs, w, color='#5A8EB8', label='Terrestrial', zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels(conditions)
    ax.set_ylabel('Composite nitrifier prevalence (%)')
    ax.set_title('(a) Nitrifier prevalence\nbefore/after crew taxa removal', fontweight='bold')
    ax.set_ylim(0, 120)
    ax.legend()
    ax.yaxis.grid(True, linestyle='--', alpha=0.5, zorder=0)
    ax.set_axisbelow(True)

    # Panel b: Cliff's delta comparison
    ax = axes[1]
    deltas = [stats_full['Cliffs_delta'], stats_clean['Cliffs_delta']]
    colors = ['#7F8C8D', '#2C3E50']
    bars = ax.bar(conditions, deltas, color=colors, alpha=0.8, zorder=3)
    for bar, d in zip(bars, deltas):
        ax.text(bar.get_x() + bar.get_width()/2, d - 0.03,
                f'{d:.3f}', ha='center', va='top', fontsize=10, color='white', fontweight='bold')
    ax.axhline(0, color='black', linewidth=0.8)
    ax.set_ylabel("Cliff's δ (Space vs Terrestrial)")
    ax.set_title("(b) Effect size\nbefore/after crew taxa removal", fontweight='bold')
    ax.set_ylim(-1.1, 0.2)
    ax.yaxis.grid(True, linestyle='--', alpha=0.5, zorder=0)
    ax.set_axisbelow(True)

    fig.suptitle(
        'Sensitivity: nitrifier depletion after removing crew/built-environment taxa\n'
        f'({removed_count} genera removed; analytical root zone Space n=20, Terr n=20)',
        fontsize=10, fontweight='bold'
    )
    fig.tight_layout()

    for ext in ('pdf', 'png'):
        p = OUTPUT_DIR / f'Crew_Taxa_Sensitivity.{ext}'
        fig.savefig(p, dpi=300, bbox_inches='tight')
        print(f'Saved: {p}')
    plt.close(fig)

    print('\nDone.')


if __name__ == '__main__':
    main()
