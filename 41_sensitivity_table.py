"""
41_sensitivity_table.py
=======================
Sample-definition sensitivity analysis for key results.

Tests whether the main findings hold across four sample filters:
  A) Plant Roots only          (Space n=16, Terr n=20)
  B) Rhizosphere only          (Space n=4,  Terr n=20)
  C) Analytical root zone      (Space n=20, Terr n=20)  ← main analysis filter
  D) PERMANOVA root-zone       (Space n=36, Terr n=20)  ← includes Arcillite

For each filter, reports:
  - Nitrification Log2FC (FAPROTAX, from pre-computed per-sample values)
  - Composite nitrifier RA: mean Space vs Terr, MWU q, Cliff's delta
  - Faith PD: mean Space vs Terr, MWU q
  - Nitrifier composite prevalence: Space % vs Terr %

Output: version-2_integrated/Sensitivity_SampleDefinition.csv
        version-2_integrated/Sensitivity_SampleDefinition.png / .pdf
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
import config

NITRIFIER_GENERA  = ['Nitrosospira', 'Nitrospira']   # only genera with ASVs in this dataset
RAREFACTION_DEPTH = config.DIVERSITY_PARAMS['rarefaction_depth']
RANDOM_SEED       = config.RANDOM_SEED
OUTPUT_DIR        = config.VERSION_DIR

rng = np.random.default_rng(RANDOM_SEED)


# ── Helpers ──────────────────────────────────────────────────────────────────

def cliffs_delta(a, b):
    a, b = np.asarray(a), np.asarray(b)
    dom  = np.sum(a[:, None] > b[None, :]) - np.sum(a[:, None] < b[None, :])
    return dom / (len(a) * len(b))


def mwu_q_single(a, b):
    """Single-comparison MWU p (no multi-test correction here; reported as-is)."""
    _, p = stats.mannwhitneyu(a, b, alternative='two-sided')
    return p


def log2fc(sp, tr, pseudo=1e-6):
    return np.log2((np.mean(sp) + pseudo) / (np.mean(tr) + pseudo))


# ── Load raw data once ───────────────────────────────────────────────────────

def load_base():
    ft   = pd.read_csv(config.FEATURE_TABLE_CLEAN, sep='\t', skiprows=1, index_col=0)
    tax  = pd.read_csv(
        config.VERSION_DIR / 'exported_taxonomy' / 'taxonomy.tsv',
        sep='\t', index_col=0
    )
    meta = pd.read_csv(config.INTEGRATED_METADATA, sep='\t', index_col=0)
    faith = pd.read_csv(
        config.VERSION_DIR / 'exported_diversity' / 'alpha-diversity.tsv',
        sep='\t', index_col=0
    )
    # keep only faith_pd column
    if 'faith_pd' in faith.columns:
        faith = faith[['faith_pd']]
    return ft, tax, meta, faith


def rarefy_and_genus(ft, meta, keep_ids):
    """Rarefy samples in keep_ids and return genus-level RA table."""
    ft_sub = ft[[c for c in keep_ids if c in ft.columns]]
    ft_t   = ft_sub.T
    valid  = ft_t[ft_t.sum(axis=1) >= RAREFACTION_DEPTH]

    def rarefy_row(row):
        pool    = np.repeat(np.arange(len(row)), row.astype(int))
        chosen  = rng.choice(pool, size=RAREFACTION_DEPTH, replace=False)
        return np.bincount(chosen, minlength=len(row)).astype(float) / RAREFACTION_DEPTH

    rarefied = valid.apply(rarefy_row, axis=1, result_type='expand')
    rarefied.columns = valid.columns
    return rarefied


def genus_table_for_samples(ft, meta, sp_ids, tr_ids):
    from pathlib import Path
    tax = pd.read_csv(
        config.VERSION_DIR / 'exported_taxonomy' / 'taxonomy.tsv',
        sep='\t', index_col=0
    )

    all_ids  = [s for s in sp_ids + tr_ids if s in ft.columns]
    rarefied = rarefy_and_genus(ft, meta, all_ids)

    def parse_genus(t):
        for p in reversed(str(t).split(';')):
            p = p.strip()
            if p.startswith('g__'):
                g = p[3:]
                if g and g not in ('', 'uncultured', 'metagenome'):
                    return g
        return 'Unknown'

    tax['genus'] = tax['Taxon'].apply(parse_genus)
    shared = rarefied.columns.intersection(tax.index)
    genus_ra = rarefied[shared].T
    genus_ra.index = tax.loc[shared, 'genus']
    genus_ra = genus_ra.groupby(level=0).sum().T   # samples × genera

    sp_in = [s for s in sp_ids if s in genus_ra.index]
    tr_in = [s for s in tr_ids if s in genus_ra.index]
    return genus_ra, sp_in, tr_in


def faprotax_log2fc_for_filter(faprotax_path, sp_ids, tr_ids):
    """
    Compute nitrification Log2FC on per-sample RELATIVE ABUNDANCE from FAPROTAX
    functional table (functions × samples, raw counts).
    Normalizes each sample to relative abundance before computing group means.
    """
    try:
        fa = pd.read_csv(faprotax_path, sep='\t', index_col=0)
    except FileNotFoundError:
        return np.nan

    nit_rows = [r for r in fa.index if 'nitrif' in r.lower()]
    if not nit_rows:
        return np.nan

    # Normalize to relative abundance per sample (column sum = total functional counts)
    col_sums = fa.sum(axis=0)
    col_sums[col_sums == 0] = 1
    fa_ra = fa.div(col_sums, axis=1)

    nit = fa_ra.loc[nit_rows].sum(axis=0)
    sp_vals = nit.reindex([s for s in sp_ids if s in nit.index]).dropna().values
    tr_vals = nit.reindex([s for s in tr_ids if s in nit.index]).dropna().values
    if len(sp_vals) == 0 or len(tr_vals) == 0:
        return np.nan
    return log2fc(sp_vals, tr_vals)


# ── Per-filter analysis ───────────────────────────────────────────────────────

def analyze_filter(ft, meta, faith, filter_name, sp_ids, tr_ids):
    sp_ids = [s for s in sp_ids if s in ft.columns]
    tr_ids = [s for s in tr_ids if s in ft.columns]
    sp_n, tr_n = len(sp_ids), len(tr_ids)

    # Faith PD
    faith_sp = faith.reindex(sp_ids).dropna().iloc[:, 0].values
    faith_tr = faith.reindex(tr_ids).dropna().iloc[:, 0].values
    faith_p  = mwu_q_single(faith_sp, faith_tr)
    faith_d  = cliffs_delta(faith_sp, faith_tr)

    # Nitrifier composite RA
    genus_ra, sp_in, tr_in = genus_table_for_samples(ft, meta, sp_ids, tr_ids)
    present_gen = [g for g in NITRIFIER_GENERA if g in genus_ra.columns]
    if present_gen:
        nit_sp = genus_ra.loc[sp_in, present_gen].sum(axis=1).values if sp_in else np.zeros(0)
        nit_tr = genus_ra.loc[tr_in, present_gen].sum(axis=1).values if tr_in else np.zeros(0)
    else:
        nit_sp = np.zeros(sp_n)
        nit_tr = np.zeros(tr_n)

    nit_p = mwu_q_single(nit_sp, nit_tr) if (nit_sp.sum() + nit_tr.sum()) > 0 else np.nan
    nit_d = cliffs_delta(nit_sp, nit_tr)
    nit_prev_sp = f"{int((nit_sp > 0).sum())}/{sp_n} ({100*(nit_sp>0).mean():.0f}%)"
    nit_prev_tr = f"{int((nit_tr > 0).sum())}/{tr_n} ({100*(nit_tr>0).mean():.0f}%)"

    # FAPROTAX Log2FC — use the pre-computed per-sample functional table
    fa_path = config.VERSION_DIR / 'FAPROTAX_functional_table.tsv'
    nit_log2fc = faprotax_log2fc_for_filter(fa_path, sp_ids, tr_ids)

    return {
        'Filter': filter_name,
        'Space_n': sp_n,
        'Terr_n': tr_n,
        'Nitrification_Log2FC': round(nit_log2fc, 2) if not np.isnan(nit_log2fc) else 'NA',
        'Faith_PD_Space_mean': round(faith_sp.mean(), 1) if len(faith_sp) else 'NA',
        'Faith_PD_Terr_mean':  round(faith_tr.mean(), 1) if len(faith_tr) else 'NA',
        'Faith_PD_MWU_p': round(faith_p, 4),
        'Faith_PD_Cliffs_delta': round(faith_d, 3),
        'Nitrifier_composite_Space': nit_prev_sp,
        'Nitrifier_composite_Terr':  nit_prev_tr,
        'Nitrifier_MWU_p': round(nit_p, 4) if not np.isnan(nit_p) else 'NA',
        'Nitrifier_Cliffs_delta': round(nit_d, 3),
    }


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    print('=== 41_sensitivity_table.py ===')

    ft, tax, meta, faith = load_base()

    sp_plant_roots  = meta.index[(meta.Study_Group == 'Space_Flight') &
                                  (meta.sample_type == 'Plant Roots')].tolist()
    sp_rhizosphere  = meta.index[(meta.Study_Group == 'Space_Flight') &
                                  (meta.sample_type == 'Rhizosphere')].tolist()
    sp_analytical   = sp_plant_roots + sp_rhizosphere          # n=20
    sp_permanova_rz = meta.index[(meta.Study_Group == 'Space_Flight') &
                                  meta.sample_type.isin(
                                      {'Plant Roots', 'Rhizosphere', 'Arcillite'}
                                  )].tolist()                  # n=36
    tr_all = meta.index[meta.Study_Group == 'Terrestrial_Soil'].tolist()

    filters = [
        ('A: Plant Roots only',        sp_plant_roots,  tr_all),
        ('B: Rhizosphere only',         sp_rhizosphere,  tr_all),
        ('C: Analytical root zone\n(n=20, main analysis)', sp_analytical, tr_all),
        ('D: PERMANOVA root-zone\n(n=36, incl. Arcillite)', sp_permanova_rz, tr_all),
    ]

    rows = []
    for name, sp_ids, tr_ids in filters:
        label = name.replace('\n', ' ')
        print(f'\n  Computing filter: {label}  (Space n={len(sp_ids)}, Terr n={len(tr_ids)})')
        row = analyze_filter(ft, meta, faith, label, sp_ids, tr_ids)
        rows.append(row)
        print(f"    Nitrification Log2FC : {row['Nitrification_Log2FC']}")
        print(f"    Nitrifier composite  : Space {row['Nitrifier_composite_Space']}  "
              f"Terr {row['Nitrifier_composite_Terr']}")
        print(f"    Faith PD MWU p       : {row['Faith_PD_MWU_p']}")

    df = pd.DataFrame(rows)
    out_csv = OUTPUT_DIR / 'Sensitivity_SampleDefinition.csv'
    df.to_csv(out_csv, index=False)
    print(f'\nSaved: {out_csv}')

    # ── Figure: compact table-style heatmap ──────────────────────────────────
    metrics = ['Nitrification_Log2FC', 'Faith_PD_Cliffs_delta', 'Nitrifier_Cliffs_delta']
    labels  = ['Nitrification\nLog₂FC', 'Faith PD\nCliff\'s δ', 'Nitrifier composite\nCliff\'s δ']
    filter_labels = [
        'A: Plant Roots only\n(Space n=16)',
        'B: Rhizosphere only\n(Space n=4)',
        'C: Analytical root zone\n(Space n=20) ★',
        'D: PERMANOVA root-zone\n(Space n=36)',
    ]

    vals = np.zeros((len(rows), len(metrics)))
    for i, row in enumerate(rows):
        for j, m in enumerate(metrics):
            v = row[m]
            try:
                vals[i, j] = float(v)
            except (ValueError, TypeError):
                vals[i, j] = np.nan

    fig, ax = plt.subplots(figsize=(9, 4.5))
    cmap = plt.cm.RdBu_r
    im = ax.imshow(vals, cmap=cmap, aspect='auto', vmin=-15, vmax=5)

    ax.set_xticks(range(len(metrics)))
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(filter_labels, fontsize=9)

    for i in range(len(rows)):
        for j in range(len(metrics)):
            v = vals[i, j]
            txt = f'{v:.2f}' if not np.isnan(v) else 'NA'
            ax.text(j, i, txt, ha='center', va='center', fontsize=9,
                    color='white' if abs(v) > 7 else 'black')

    plt.colorbar(im, ax=ax, shrink=0.8, label='Value (Log₂FC or Cliff\'s δ)')
    ax.set_title(
        'Sensitivity analysis: key statistics across sample-definition filters\n'
        '★ = main analytical filter used throughout manuscript',
        fontsize=10, fontweight='bold'
    )
    fig.tight_layout()

    for ext in ('pdf', 'png'):
        p = OUTPUT_DIR / f'Sensitivity_SampleDefinition.{ext}'
        fig.savefig(p, dpi=300, bbox_inches='tight')
        print(f'Saved: {p}')
    plt.close(fig)

    print('\nDone.')


if __name__ == '__main__':
    main()
