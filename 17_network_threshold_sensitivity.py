"""
Network Threshold Sensitivity Analysis
=======================================
Reconstructs Space Flight and Terrestrial Soil co-occurrence networks at
three Spearman correlation thresholds (|r| > 0.3, 0.4, 0.5) and evaluates
whether key findings (hub genus identity, fragility index, edge count) are
robust across threshold choices.

Key metrics evaluated per threshold:
  - Edge count (network density)
  - Top hub genus by degree
  - Rhodanobacter degree (known keystone taxon)
  - Fragility index = R_random / R_targeted (from Script 20)
  - R² consistency with PERMANOVA result

Ecological rationale:
  Network construction depends on the correlation threshold. If the main
  conclusion ("terrestrial network has a denser hub structure that collapses
  under targeted removal") only holds at |r|=0.4, it is threshold-dependent.
  Showing consistency across 0.3–0.5 strengthens confidence.

Prerequisites:
  - version-2_integrated/exported_table_space/feature-table.tsv
  - version-2_integrated/exported_table_terrestrial/feature-table.tsv
  - version-2_integrated/exported_taxonomy/taxonomy.tsv

Run:
  python 27_network_threshold_sensitivity.py

Output (version-2_integrated/):
  - NetworkSensitivity_Results.csv
  - Main_Fig_NetworkSensitivity.png/.pdf
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import config

import pandas as pd
import numpy as np
import networkx as nx
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec

OUTPUT_DIR    = config.VERSION_DIR
P_THRESHOLD   = config.NETWORK_PARAMS['p_threshold']
TOP_N         = 15
N_RANDOM_ITER = 500   # fewer iterations for speed
RANDOM_SEED   = config.RANDOM_SEED

R_THRESHOLDS  = [0.3, 0.4, 0.5]

COLORS  = {'Space_Flight': '#3498db', 'Terrestrial_Soil': '#2ecc71'}
LABELS  = {'Space_Flight': 'Space Flight', 'Terrestrial_Soil': 'Terrestrial Soil'}
MARKERS = {0.3: 'o', 0.4: 's', 0.5: '^'}


# ─────────────────────────────────────────────────────────────────────────────
# Data helpers (same as Scripts 10, 20)
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
                name = part.split('__')[-1].strip()
                return name if name else 'Unassigned'
        return 'Unassigned'

    tax['Genus'] = tax['Taxon'].apply(get_genus)
    shared = table.index.intersection(tax.index)
    table  = table.loc[shared]
    table['Genus'] = tax.loc[shared, 'Genus']
    genus_table = table.groupby('Genus').sum()
    genus_table = genus_table.drop(index='Unassigned', errors='ignore')
    return genus_table.div(genus_table.sum(axis=0), axis=1) * 100


def select_top_shared_genera(sp_rel, tr_rel, n=TOP_N):
    shared = sp_rel.index.intersection(tr_rel.index)
    combined_mean = (sp_rel.loc[shared].mean(axis=1) +
                     tr_rel.loc[shared].mean(axis=1)) / 2
    return combined_mean.sort_values(ascending=False).head(n).index


def build_network(genus_rel, top_genera, r_threshold):
    df   = genus_rel.loc[top_genera]
    n    = len(top_genera)
    data = df.values.T
    corr, pvals = spearmanr(data)
    if n == 2:
        corr  = np.array([[1, corr], [corr, 1]])
        pvals = np.array([[0, pvals], [pvals, 0]])
    G = nx.Graph()
    G.add_nodes_from(top_genera)
    for i in range(n):
        for j in range(i + 1, n):
            r, p = corr[i, j], pvals[i, j]
            if abs(r) > r_threshold and p < P_THRESHOLD and not np.isnan(r):
                G.add_edge(top_genera[i], top_genera[j], weight=r)
    return G


# ─────────────────────────────────────────────────────────────────────────────
# Attack simulations (same logic as Script 20)
# ─────────────────────────────────────────────────────────────────────────────

def targeted_attack(G):
    N     = len(G.nodes())
    G_c   = G.copy()
    sizes = [1.0]
    while G_c.number_of_nodes() > 0:
        node = max(G_c.degree(), key=lambda x: (x[1], x[0]))[0]
        G_c.remove_node(node)
        n_left = G_c.number_of_nodes()
        if n_left > 0:
            sizes.append(len(max(nx.connected_components(G_c), key=len)) / N)
        else:
            sizes.append(0.0)
    return float(np.mean(sizes))


def random_attack(G, n_iter=N_RANDOM_ITER, seed=RANDOM_SEED):
    N     = len(G.nodes())
    nodes = list(G.nodes())
    rng   = np.random.default_rng(seed)
    all_sizes = np.zeros((n_iter, N + 1))
    all_sizes[:, 0] = 1.0
    for it in range(n_iter):
        G_c   = G.copy()
        order = rng.permutation(nodes)
        for i, node in enumerate(order):
            G_c.remove_node(node)
            n_left = G_c.number_of_nodes()
            all_sizes[it, i + 1] = (
                len(max(nx.connected_components(G_c), key=len)) / N
                if n_left > 0 else 0.0
            )
    return float(np.mean(all_sizes))


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    print('=' * 60)
    print('Network Threshold Sensitivity Analysis')
    print('=' * 60)

    print('\n[1] Loading feature tables...')
    sp_rel  = aggregate_to_genus(config.FEATURE_TABLE_SPACE, config.TAXONOMY_FILE)
    tr_rel  = aggregate_to_genus(config.FEATURE_TABLE_TERR,  config.TAXONOMY_FILE)
    # Root-zone filter: restrict Space_Flight to Plant Roots + Rhizosphere
    ROOT_ZONE_TYPES = {'Plant Roots', 'Rhizosphere'}
    meta_all = pd.read_csv(config.INTEGRATED_METADATA, sep='\t', index_col=0)
    root_zone_ids = meta_all.index[
        (meta_all['Study_Group'] == 'Space_Flight') &
        meta_all['sample_type'].isin(ROOT_ZONE_TYPES)
    ].tolist()
    sp_rel = sp_rel[[c for c in sp_rel.columns if c in root_zone_ids]]
    print(f'  Space genus table (root-zone filter): {sp_rel.shape[0]} genera × {sp_rel.shape[1]} samples')
    top_gen = select_top_shared_genera(sp_rel, tr_rel)
    print(f'  Top {TOP_N} shared genera selected')
    print(f'  Genera: {list(top_gen)}')

    rows = []
    for r_thr in R_THRESHOLDS:
        print(f'\n[2] Building networks at |r| > {r_thr}...')
        G_sp = build_network(sp_rel, top_gen, r_thr)
        G_tr = build_network(tr_rel, top_gen, r_thr)

        n_edges_sp = G_sp.number_of_edges()
        n_edges_tr = G_tr.number_of_edges()
        print(f'  Space:       {G_sp.number_of_nodes()} nodes, {n_edges_sp} edges')
        print(f'  Terrestrial: {G_tr.number_of_nodes()} nodes, {n_edges_tr} edges')

        # Hub genus (highest degree)
        def top_hub(G):
            if G.number_of_edges() == 0:
                return 'none', 0
            node, deg = max(G.degree(), key=lambda x: x[1])
            return node, deg

        hub_sp, hub_sp_deg = top_hub(G_sp)
        hub_tr, hub_tr_deg = top_hub(G_tr)

        # Rhodanobacter degree
        rhod_sp = G_sp.degree('Rhodanobacter') if 'Rhodanobacter' in G_sp else 0
        rhod_tr = G_tr.degree('Rhodanobacter') if 'Rhodanobacter' in G_tr else 0

        # Robustness
        print(f'  Running attack simulations (r={r_thr})...')
        R_sp_t = targeted_attack(G_sp) if G_sp.number_of_edges() > 0 else np.nan
        R_tr_t = targeted_attack(G_tr) if G_tr.number_of_edges() > 0 else np.nan
        R_sp_r = random_attack(G_sp)   if G_sp.number_of_edges() > 0 else np.nan
        R_tr_r = random_attack(G_tr)   if G_tr.number_of_edges() > 0 else np.nan
        frag_sp = R_sp_r / R_sp_t if (R_sp_t and not np.isnan(R_sp_t) and R_sp_t > 0) else np.nan
        frag_tr = R_tr_r / R_tr_t if (R_tr_t and not np.isnan(R_tr_t) and R_tr_t > 0) else np.nan

        print(f'  Space R_targ={R_sp_t:.3f}, R_rand={R_sp_r:.3f}, Fragility={frag_sp:.3f}')
        print(f'  Terr  R_targ={R_tr_t:.3f}, R_rand={R_tr_r:.3f}, Fragility={frag_tr:.3f}')

        for group, G, n_edges, hub, hub_deg, rhod_deg, R_t, R_r, frag in [
            ('Space_Flight',     G_sp, n_edges_sp, hub_sp, hub_sp_deg, rhod_sp,
             R_sp_t, R_sp_r, frag_sp),
            ('Terrestrial_Soil', G_tr, n_edges_tr, hub_tr, hub_tr_deg, rhod_tr,
             R_tr_t, R_tr_r, frag_tr),
        ]:
            rows.append({
                'r_threshold':   r_thr,
                'Group':         group,
                'n_edges':       n_edges,
                'Top_Hub':       hub,
                'Hub_Degree':    hub_deg,
                'Rhodanobacter_Degree': rhod_deg,
                'R_targeted':    round(R_t, 4) if not np.isnan(R_t) else None,
                'R_random':      round(R_r, 4) if not np.isnan(R_r) else None,
                'Fragility':     round(frag, 4) if not np.isnan(frag) else None,
            })

    df = pd.DataFrame(rows)
    print('\n--- Results ---')
    print(df.to_string(index=False))

    print('\n[3] Saving results...')
    df.to_csv(OUTPUT_DIR / 'SuppS5_NetworkSensitivity_Results.csv', index=False)
    print('  Saved: SuppS5_NetworkSensitivity_Results.csv')
    df.to_csv(OUTPUT_DIR / 'SuppS2b_ThresholdSensitivity_Results.csv', index=False)
    print('  Saved: SuppS2b_ThresholdSensitivity_Results.csv')

    # ─────────────────────────────────────────────────────────────────────────
    # Plot
    # ─────────────────────────────────────────────────────────────────────────
    print('\n[4] Plotting...')

    fig = plt.figure(figsize=(14, 9))
    gs  = GridSpec(2, 3, figure=fig, hspace=0.45, wspace=0.35)

    ax_edges   = fig.add_subplot(gs[0, 0])
    ax_hub     = fig.add_subplot(gs[0, 1])
    ax_rhod    = fig.add_subplot(gs[0, 2])
    ax_rtarg   = fig.add_subplot(gs[1, 0])
    ax_rrand   = fig.add_subplot(gs[1, 1])
    ax_frag    = fig.add_subplot(gs[1, 2])

    metric_axs = [
        (ax_edges, 'n_edges',          'Edge Count',              False),
        (ax_hub,   'Hub_Degree',        'Top Hub Degree',          False),
        (ax_rhod,  'Rhodanobacter_Degree', 'Rhodanobacter Degree', False),
        (ax_rtarg, 'R_targeted',        'Robustness R (Targeted)', False),
        (ax_rrand, 'R_random',          'Robustness R (Random)',   False),
        (ax_frag,  'Fragility',         'Fragility Index',         True),
    ]

    for ax, col, ylabel, add_hline in metric_axs:
        for group in ['Space_Flight', 'Terrestrial_Soil']:
            sub = df[df['Group'] == group]
            ax.plot(sub['r_threshold'], sub[col],
                    color=COLORS[group], marker='o', lw=2, ms=7,
                    label=LABELS[group])
        ax.set_xlabel('|r| Threshold', fontsize=10)
        ax.set_ylabel(ylabel, fontsize=10)
        ax.set_title(ylabel, fontsize=10, fontweight='bold')
        ax.set_xticks(R_THRESHOLDS)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        if add_hline:
            ax.axhline(1.0, color='black', lw=1, linestyle='--', alpha=0.5,
                       label='Fragility=1')

    fig.suptitle(
        f'Network Threshold Sensitivity Analysis\n'
        f'Top {TOP_N} shared genera | p < {P_THRESHOLD} | |r| ∈ {R_THRESHOLDS}',
        fontsize=12, fontweight='bold'
    )

    for ext in ['png', 'pdf']:
        out = OUTPUT_DIR / f'SuppS5_NetworkSensitivity.{ext}'
        fig.savefig(out, dpi=300 if ext == 'png' else None, bbox_inches='tight')
        print(f'  Saved: {out.name}')

    # Current Biology supplementary figure output
    cb_dir = config.CURRENT_BIOLOGY_DIR
    cb_dir.mkdir(exist_ok=True)
    fig.savefig(cb_dir / 'SuppS2b_ThresholdSensitivity.png', dpi=300, bbox_inches='tight')
    fig.savefig(cb_dir / 'SuppS2b_ThresholdSensitivity.pdf', bbox_inches='tight')
    print(f'  Saved Current Biology: {cb_dir / "SuppS2b_ThresholdSensitivity.png"}')
    plt.close()

    print('\nDone.')
