"""
40_nitrifier_direct_analysis.py
================================
Direct taxon-level analysis of canonical nitrifier genera in ISS vs Terrestrial
pepper root-zone samples.

Outputs
-------
version-2_integrated/
  Fig_Nitrifier_Direct.pdf / .png   — 4-panel figure
  Nitrifier_Direct_Results.csv      — per-genus + composite stats with effect sizes

Analysis
--------
For Nitrosomonas, Nitrosospira, Nitrospira, Nitrobacter (and composite burden):
  - Prevalence per group + Fisher's exact test
  - Relative abundance (rarefied, genus-level) per sample
  - Mann-Whitney U + Cliff's delta (effect size) + bootstrap 95% CI on mean RA
  - Composite (sum of 4 genera): same stats

Sample filter: analytical root zone only
  Space_Flight: Plant Roots (n=16) + Rhizosphere (n=4) = n=20
  Terrestrial_Soil: all (n=20)
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
import config

# ── Parameters ──────────────────────────────────────────────────────────────
NITRIFIER_GENERA   = ['Nitrosomonas', 'Nitrosospira', 'Nitrospira', 'Nitrobacter']
ROOT_ZONE_TYPES    = {'Plant Roots', 'Rhizosphere'}
RAREFACTION_DEPTH  = config.DIVERSITY_PARAMS['rarefaction_depth']   # 1000
N_BOOTSTRAP        = 2000
RANDOM_SEED        = config.RANDOM_SEED
OUTPUT_DIR         = config.VERSION_DIR

rng = np.random.default_rng(RANDOM_SEED)


# ── Helpers ──────────────────────────────────────────────────────────────────

def cliffs_delta(a, b):
    """Cliff's delta effect size; range [-1, 1]; positive = a > b."""
    a, b = np.asarray(a), np.asarray(b)
    dominance = np.sum(a[:, None] > b[None, :]) - np.sum(a[:, None] < b[None, :])
    return dominance / (len(a) * len(b))


def bootstrap_mean_ci(x, n=N_BOOTSTRAP, rng=rng, ci=0.95):
    """Percentile bootstrap CI of the mean."""
    x = np.asarray(x, dtype=float)
    means = np.array([rng.choice(x, size=len(x), replace=True).mean()
                      for _ in range(n)])
    lo = np.percentile(means, 100 * (1 - ci) / 2)
    hi = np.percentile(means, 100 * (1 + ci) / 2)
    return lo, hi


def prevalence_fisher(present_sp, total_sp, present_terr, total_terr):
    """Fisher's exact test on a 2×2 prevalence table."""
    table = [[present_sp,  total_sp   - present_sp],
             [present_terr, total_terr - present_terr]]
    _, p = stats.fisher_exact(table, alternative='two-sided')
    return p


# ── Load data ────────────────────────────────────────────────────────────────

def load_data():
    # Feature table (raw counts, samples as columns)
    ft = pd.read_csv(config.FEATURE_TABLE_CLEAN, sep='\t', skiprows=1, index_col=0)

    # Taxonomy
    tax = pd.read_csv(
        config.VERSION_DIR / 'exported_taxonomy' / 'taxonomy.tsv',
        sep='\t', index_col=0
    )

    # Metadata
    meta = pd.read_csv(config.INTEGRATED_METADATA, sep='\t', index_col=0)

    # ── Analytical root zone filter ──────────────────────────────────────────
    space_ids = meta.index[
        (meta['Study_Group'] == 'Space_Flight') &
        meta['sample_type'].isin(ROOT_ZONE_TYPES)
    ].tolist()
    terr_ids  = meta.index[meta['Study_Group'] == 'Terrestrial_Soil'].tolist()
    keep_ids  = [s for s in (space_ids + terr_ids) if s in ft.columns]
    ft   = ft[keep_ids]
    meta = meta.loc[[s for s in keep_ids if s in meta.index]]

    # ── Rarefy ───────────────────────────────────────────────────────────────
    depth = RAREFACTION_DEPTH
    ft_t  = ft.T  # samples × ASVs
    valid = ft_t[ft_t.sum(axis=1) >= depth]

    def rarefy_row(row):
        pool = np.repeat(np.arange(len(row)), row.astype(int))
        chosen = rng.choice(pool, size=depth, replace=False)
        counts = np.bincount(chosen, minlength=len(row)).astype(float)
        return counts / depth   # relative abundance

    rarefied = valid.apply(rarefy_row, axis=1, result_type='expand')
    rarefied.columns = valid.columns

    # ── Genus-level collapse ──────────────────────────────────────────────────
    # Parse genus from SILVA taxonomy string
    def parse_genus(taxon_str):
        parts = [p.strip() for p in str(taxon_str).split(';')]
        # SILVA format: last meaningful level often g__GenusName
        for p in reversed(parts):
            if p.startswith('g__'):
                genus = p[3:]
                if genus and genus not in ('', 'uncultured', 'metagenome'):
                    return genus
        return 'Unknown'

    tax['genus'] = tax['Taxon'].apply(parse_genus)

    # Keep only ASVs present in rarefied table
    shared_asvs = rarefied.columns.intersection(tax.index)
    tax_sub = tax.loc[shared_asvs]

    # Group by genus, sum RA across ASVs
    genus_ra = rarefied[shared_asvs].T  # ASVs × samples
    genus_ra.index = tax_sub['genus']
    genus_ra = genus_ra.groupby(level=0).sum().T  # samples × genera

    # Attach group label
    meta_sub = meta.loc[genus_ra.index.intersection(meta.index)]
    genus_ra  = genus_ra.loc[meta_sub.index]

    space_mask = meta_sub['Study_Group'] == 'Space_Flight'
    terr_mask  = meta_sub['Study_Group'] == 'Terrestrial_Soil'

    return genus_ra, meta_sub, space_mask, terr_mask


# ── Per-genus statistics ──────────────────────────────────────────────────────

def compute_stats(genus_ra, space_mask, terr_mask):
    rows = []

    genera_to_test = NITRIFIER_GENERA + ['COMPOSITE']

    for genus in genera_to_test:
        if genus == 'COMPOSITE':
            # Sum of all 4 genera present in data
            present_gen = [g for g in NITRIFIER_GENERA if g in genus_ra.columns]
            sp_vals  = genus_ra.loc[space_mask, present_gen].sum(axis=1).values
            tr_vals  = genus_ra.loc[terr_mask,  present_gen].sum(axis=1).values
            label = 'Composite (4 genera)'
        else:
            if genus not in genus_ra.columns:
                sp_vals = np.zeros(space_mask.sum())
                tr_vals = np.zeros(terr_mask.sum())
            else:
                sp_vals = genus_ra.loc[space_mask, genus].values
                tr_vals = genus_ra.loc[terr_mask,  genus].values
            label = genus

        # Prevalence
        sp_prev = (sp_vals > 0).sum()
        tr_prev = (tr_vals > 0).sum()
        sp_n    = len(sp_vals)
        tr_n    = len(tr_vals)
        fisher_p = prevalence_fisher(sp_prev, sp_n, tr_prev, tr_n)

        # MWU
        stat, mwu_p = stats.mannwhitneyu(sp_vals, tr_vals, alternative='two-sided')

        # Cliff's delta (Space vs Terrestrial)
        delta = cliffs_delta(sp_vals, tr_vals)

        # Bootstrap CI on group means
        sp_ci = bootstrap_mean_ci(sp_vals)
        tr_ci = bootstrap_mean_ci(tr_vals)

        rows.append({
            'Genus': label,
            'Space_n': sp_n, 'Terr_n': tr_n,
            'Space_prevalent': sp_prev, 'Terr_prevalent': tr_prev,
            'Space_prev_pct': round(100 * sp_prev / sp_n, 1),
            'Terr_prev_pct':  round(100 * tr_prev / tr_n, 1),
            'Fisher_p': fisher_p,
            'Space_mean_RA': sp_vals.mean(),
            'Terr_mean_RA':  tr_vals.mean(),
            'Space_CI_lo': sp_ci[0], 'Space_CI_hi': sp_ci[1],
            'Terr_CI_lo':  tr_ci[0], 'Terr_CI_hi':  tr_ci[1],
            'MWU_p': mwu_p,
            'Cliffs_delta': delta,
            'Space_vals': sp_vals,
            'Terr_vals':  tr_vals,
        })

    return rows


def apply_fdr(rows):
    from scipy.stats import false_discovery_control
    mwu_ps = np.array([r['MWU_p'] for r in rows])
    fisher_ps = np.array([r['Fisher_p'] for r in rows])
    try:
        mwu_q    = false_discovery_control(mwu_ps, method='bh')
        fisher_q = false_discovery_control(fisher_ps, method='bh')
    except AttributeError:
        # scipy < 1.11 fallback
        from statsmodels.stats.multitest import multipletests
        _, mwu_q,    _, _ = multipletests(mwu_ps,    method='fdr_bh')
        _, fisher_q, _, _ = multipletests(fisher_ps, method='fdr_bh')

    for r, mq, fq in zip(rows, mwu_q, fisher_q):
        r['MWU_q']    = mq
        r['Fisher_q'] = fq
    return rows


def sig_stars(q):
    if q < 0.001: return '***'
    if q < 0.01:  return '**'
    if q < 0.05:  return '*'
    return 'ns'


# ── Figure ───────────────────────────────────────────────────────────────────

def make_figure(rows):
    genera_order = NITRIFIER_GENERA + ['COMPOSITE']
    label_map = {g: g for g in NITRIFIER_GENERA}
    label_map['COMPOSITE'] = 'Composite\n(4 genera)'

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # ── Panel A: Prevalence bar chart ────────────────────────────────────────
    ax = axes[0]
    x    = np.arange(len(rows))
    w    = 0.35
    sp_prev = [r['Space_prev_pct'] for r in rows]
    tr_prev = [r['Terr_prev_pct']  for r in rows]
    fisher_qs = [r['Fisher_q'] for r in rows]

    bars_sp = ax.bar(x - w/2, sp_prev, w, color='#E8735A', label='ISS (Space Flight)', zorder=3)
    bars_tr = ax.bar(x + w/2, tr_prev, w, color='#5A8EB8', label='Terrestrial',        zorder=3)

    for i, (sp, tr, q) in enumerate(zip(sp_prev, tr_prev, fisher_qs)):
        ymax = max(sp, tr) + 3
        stars = sig_stars(q)
        ax.text(x[i], ymax, stars, ha='center', va='bottom', fontsize=9)

    ax.set_xticks(x)
    ax.set_xticklabels([label_map.get(r['Genus'], r['Genus']) for r in rows],
                       fontsize=9)
    ax.set_ylabel('Prevalence (%)', fontsize=10)
    ax.set_title('(a) Prevalence of nitrifier genera\n(Fisher\'s exact, BH-corrected)',
                 fontsize=10, fontweight='bold')
    ax.set_ylim(0, 120)
    ax.legend(fontsize=8)
    ax.yaxis.grid(True, linestyle='--', alpha=0.5, zorder=0)
    ax.set_axisbelow(True)

    # ── Panel B: Mean RA + bootstrap CI + individual points ─────────────────
    ax = axes[1]
    sp_means = [r['Space_mean_RA'] * 100 for r in rows]
    tr_means = [r['Terr_mean_RA']  * 100 for r in rows]
    sp_lo    = [r['Space_CI_lo']   * 100 for r in rows]
    sp_hi    = [r['Space_CI_hi']   * 100 for r in rows]
    tr_lo    = [r['Terr_CI_lo']    * 100 for r in rows]
    tr_hi    = [r['Terr_CI_hi']    * 100 for r in rows]
    mwu_qs   = [r['MWU_q']         for r in rows]
    deltas   = [r['Cliffs_delta']  for r in rows]

    for i, row in enumerate(rows):
        sp_vals_pct = row['Space_vals'] * 100
        tr_vals_pct = row['Terr_vals']  * 100
        ax.scatter(np.full(len(sp_vals_pct), x[i] - w/2),
                   sp_vals_pct, color='#E8735A', alpha=0.4, s=18, zorder=4)
        ax.scatter(np.full(len(tr_vals_pct), x[i] + w/2),
                   tr_vals_pct, color='#5A8EB8', alpha=0.4, s=18, zorder=4)

    ax.bar(x - w/2, sp_means, w, color='#E8735A', alpha=0.6, zorder=3)
    ax.bar(x + w/2, tr_means, w, color='#5A8EB8', alpha=0.6, zorder=3)

    for i in range(len(rows)):
        ax.errorbar(x[i] - w/2, sp_means[i],
                    yerr=[[sp_means[i] - sp_lo[i]], [sp_hi[i] - sp_means[i]]],
                    fmt='none', color='#c0392b', linewidth=1.5, capsize=3, zorder=5)
        ax.errorbar(x[i] + w/2, tr_means[i],
                    yerr=[[tr_means[i] - tr_lo[i]], [tr_hi[i] - tr_means[i]]],
                    fmt='none', color='#1a5276', linewidth=1.5, capsize=3, zorder=5)

    for i, (q, d) in enumerate(zip(mwu_qs, deltas)):
        ymax = max(sp_means[i] + (sp_hi[i] - sp_means[i]),
                   tr_means[i] + (tr_hi[i] - tr_means[i])) + 0.005
        stars = sig_stars(q)
        ax.text(x[i], ymax, f'{stars}\nδ={d:.2f}', ha='center', va='bottom',
                fontsize=7, linespacing=1.4)

    ax.set_xticks(x)
    ax.set_xticklabels([label_map.get(r['Genus'], r['Genus']) for r in rows],
                       fontsize=9)
    ax.set_ylabel('Relative abundance (%)', fontsize=10)
    ax.set_title("(b) Mean relative abundance ± 95% bootstrap CI\n"
                 "(MWU, BH-corrected; Cliff's δ shown)",
                 fontsize=10, fontweight='bold')
    ax.yaxis.grid(True, linestyle='--', alpha=0.5, zorder=0)
    ax.set_axisbelow(True)

    fig.suptitle(
        'Direct nitrifier taxon analysis: ISS vs Terrestrial\n'
        'Analytical root zone: Space n=20 (Plant Roots + Rhizosphere), Terrestrial n=20',
        fontsize=11, fontweight='bold', y=1.01
    )
    fig.tight_layout()
    return fig


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    print('=== 40_nitrifier_direct_analysis.py ===')
    print(f'Rarefaction depth : {RAREFACTION_DEPTH}')
    print(f'Bootstrap iterations: {N_BOOTSTRAP}')

    print('\n[1] Loading data and rarefying...')
    genus_ra, meta_sub, space_mask, terr_mask = load_data()
    print(f'    Space n={space_mask.sum()}, Terrestrial n={terr_mask.sum()}')

    # Report which nitrifier genera have any reads
    for g in NITRIFIER_GENERA:
        if g in genus_ra.columns:
            sp_det = (genus_ra.loc[space_mask, g] > 0).sum()
            tr_det = (genus_ra.loc[terr_mask,  g] > 0).sum()
            print(f'    {g}: Space {sp_det}/{space_mask.sum()} detected, '
                  f'Terr {tr_det}/{terr_mask.sum()} detected')
        else:
            print(f'    {g}: NOT FOUND in rarefied genus table')

    print('\n[2] Computing statistics...')
    rows = compute_stats(genus_ra, space_mask, terr_mask)
    rows = apply_fdr(rows)

    # Print summary table
    print('\n  Genus                   | Space prev% | Terr prev% | Fisher q '
          '| MWU q   | Cliff δ')
    print('  ' + '-' * 75)
    for r in rows:
        print(f"  {r['Genus']:<24} | "
              f"{r['Space_prev_pct']:>11.1f} | "
              f"{r['Terr_prev_pct']:>10.1f} | "
              f"{r['Fisher_q']:>8.4f} | "
              f"{r['MWU_q']:>7.4f} | "
              f"{r['Cliffs_delta']:>7.3f}")

    print('\n[3] Saving CSV...')
    csv_rows = []
    for r in rows:
        csv_rows.append({k: v for k, v in r.items()
                         if k not in ('Space_vals', 'Terr_vals')})
    df_out = pd.DataFrame(csv_rows)
    out_csv = OUTPUT_DIR / 'Nitrifier_Direct_Results.csv'
    df_out.to_csv(out_csv, index=False, float_format='%.6f')
    print(f'    Saved: {out_csv}')

    print('\n[4] Generating figure...')
    fig = make_figure(rows)
    for ext in ('pdf', 'png'):
        p = OUTPUT_DIR / f'Fig_Nitrifier_Direct.{ext}'
        fig.savefig(p, dpi=300, bbox_inches='tight')
        print(f'    Saved: {p}')
    plt.close(fig)

    print('\nDone.')


if __name__ == '__main__':
    main()
