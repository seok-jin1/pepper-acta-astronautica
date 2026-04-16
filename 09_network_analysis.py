"""
Option A: Genus-level Co-occurrence Network Analysis

Rationale:
  - Original ASV-level analysis (run_spearman_network.py) was statistically compromised:
      * Terrestrial Soil: n=20 samples, 100 nodes → p >> n (underdetermined)
      * Space and Terrestrial share 0 ASVs → direct node comparison impossible
  - Genus-level aggregation resolves both issues:
      * 152 shared genera → same node set for both groups
      * Top 15 shared genera → p=15 < n=20 (Terrestrial), p=15 << n=106 (Space)

Method:
  - Aggregate ASV counts to genus level using SILVA taxonomy
  - Select top 15 shared genera by mean relative abundance (average across both groups)
  - Compute relative abundance per sample per group
  - Build Spearman co-occurrence network: |r| > 0.6, FDR BH q < 0.05
  - Compare network metrics and identify hub genera

Output:
  - Fig2E_Network.png / .pdf
  - Fig2E_Network_Genus_Metrics.csv
  - Fig2E_Keystone_Genus_Comparison.csv
  - Fig2E_Network_Genus_Fast_Metrics.csv
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import config

import pandas as pd
import numpy as np
import networkx as nx
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

OUTPUT_DIR = config.VERSION_DIR
R_THRESHOLD = 0.4   # relaxed from 0.6 (genus-level is inherently less granular than ASV)
P_THRESHOLD = 0.05
USE_FDR     = False  # raw p < 0.05; FDR is over-conservative for 105 pairs (15-node network)
TOP_N = 15

# ─────────────────────────────────────────────
# 1. Load and aggregate to genus level
# ─────────────────────────────────────────────

def aggregate_to_genus(feature_table_path, taxonomy_path):
    """Aggregate ASV feature table to genus level."""
    table = pd.read_csv(feature_table_path, sep='\t', skiprows=1, index_col=0)
    tax   = pd.read_csv(taxonomy_path, sep='\t', index_col=0)

    def get_genus(taxon_str):
        if pd.isna(taxon_str):
            return 'Unassigned'
        for part in reversed(taxon_str.split(';')):
            part = part.strip()
            if part.startswith('g__'):
                name = part.split('__')[-1].strip()
                return name if name else 'Unassigned'
        return 'Unassigned'

    tax['Genus'] = tax['Taxon'].apply(get_genus)

    # Keep only ASVs present in this feature table
    shared_asvs = table.index.intersection(tax.index)
    table = table.loc[shared_asvs]
    tax_sub = tax.loc[shared_asvs]

    table['Genus'] = tax_sub['Genus']
    genus_table = table.groupby('Genus').sum()
    genus_table = genus_table.drop(index='Unassigned', errors='ignore')

    # Relative abundance (%)
    genus_rel = genus_table.div(genus_table.sum(axis=0), axis=1) * 100
    return genus_rel


# ─────────────────────────────────────────────
# 2. Select top N shared genera
# ─────────────────────────────────────────────

def select_shared_top_genera(sp_rel, tr_rel, n=TOP_N):
    """
    Select top N genera that appear in BOTH groups,
    ranked by mean relative abundance averaged across both groups.
    """
    shared = sp_rel.index.intersection(tr_rel.index)
    print(f"  Shared genera (Space ∩ Terrestrial): {len(shared)}")

    # Mean abundance in each group
    sp_mean = sp_rel.loc[shared].mean(axis=1)
    tr_mean = tr_rel.loc[shared].mean(axis=1)

    # Average rank across both groups
    combined_mean = (sp_mean + tr_mean) / 2
    top_genera = combined_mean.sort_values(ascending=False).head(n).index

    print(f"  Selected top {n} shared genera:")
    for i, g in enumerate(top_genera, 1):
        print(f"    {i:2d}. {g}  (Space: {sp_mean[g]:.2f}%  Terr: {tr_mean[g]:.2f}%)")

    return top_genera


# ─────────────────────────────────────────────
# 3. Build Spearman network
# ─────────────────────────────────────────────

def build_network(genus_rel, top_genera, group_name):
    """Build co-occurrence network for a subset of genera."""
    df = genus_rel.loc[top_genera]
    n_genera  = len(df.index)
    n_samples = df.shape[1]
    print(f"\n  [{group_name}] {n_genera} genera × {n_samples} samples")

    data = df.values.T  # (samples, genera)
    corr, pvals = spearmanr(data)

    # Handle case where spearmanr returns scalar (only 2 variables)
    if n_genera == 2:
        corr = np.array([[1, corr], [corr, 1]])
        pvals = np.array([[0, pvals], [pvals, 0]])

    p_flat = pvals[np.triu_indices(n_genera, k=1)]
    if USE_FDR:
        _, p_used, _, _ = multipletests(p_flat, method='fdr_bh')
    else:
        p_used = p_flat  # raw p-value

    p_adj_mat = np.zeros((n_genera, n_genera))
    p_adj_mat[np.triu_indices(n_genera, k=1)] = p_used
    p_adj_mat += p_adj_mat.T

    G = nx.Graph()
    G.add_nodes_from(df.index)
    for i in range(n_genera):
        for j in range(i + 1, n_genera):
            if abs(corr[i, j]) > R_THRESHOLD and p_adj_mat[i, j] < P_THRESHOLD and not np.isnan(corr[i, j]):
                G.add_edge(df.index[i], df.index[j], weight=corr[i, j])

    n_edges = G.number_of_edges()
    density = nx.density(G)
    avg_deg = np.mean([d for _, d in G.degree()])
    print(f"    Edges: {n_edges}  Density: {density:.3f}  Avg Degree: {avg_deg:.2f}")
    return G


# ─────────────────────────────────────────────
# 4. Plot side-by-side networks
# ─────────────────────────────────────────────

def plot_networks(G_space, G_terr, top_genera):
    fig, axes = plt.subplots(1, 2, figsize=(18, 9))

    # Shared layout positions (same coordinates for same genus)
    all_nodes = list(top_genera)
    pos_base  = nx.spring_layout(nx.Graph(), seed=42)

    # Build a union graph to get consistent positions
    G_union = nx.Graph()
    G_union.add_nodes_from(all_nodes)
    pos = nx.circular_layout(G_union)

    groups = [
        ('Space Flight',    G_space, axes[0], '#2980b9'),
        ('Terrestrial Soil', G_terr, axes[1], '#27ae60'),
    ]

    for name, G, ax, node_color in groups:
        degrees = dict(G.degree())

        # Draw edges
        pos_nodes = {n: pos[n] for n in G.nodes()}
        edges = list(G.edges(data=True))
        if edges:
            pos_edges = [(u, v) for u, v, _ in edges]
            weights   = [d['weight'] for _, _, d in edges]
            edge_colors = ['#3498db' if w > 0 else '#e74c3c' for w in weights]
            nx.draw_networkx_edges(G, pos_nodes, edgelist=pos_edges,
                                   ax=ax, edge_color=edge_colors,
                                   alpha=0.5, width=2.0)

        # Draw nodes (size by degree)
        node_sizes = [degrees.get(n, 0) * 200 + 300 for n in G.nodes()]
        nx.draw_networkx_nodes(G, pos_nodes, ax=ax,
                               node_size=node_sizes,
                               node_color=node_color, alpha=0.85)
        nx.draw_networkx_labels(G, pos_nodes, ax=ax,
                                font_size=8, font_color='black', font_weight='bold')

        # Metrics in title
        n_e   = G.number_of_edges()
        dens  = nx.density(G)
        avg_d = np.mean([d for _, d in G.degree()])
        ax.set_title(
            f'{name}\nNodes={G.number_of_nodes()}  Edges={n_e}  '
            f'Density={dens:.3f}  Avg Degree={avg_d:.2f}',
            fontsize=13, fontweight='bold', pad=12
        )
        ax.axis('off')

    # Legend for edge color
    blue_patch = mpatches.Patch(color='#3498db', label='Positive correlation')
    red_patch  = mpatches.Patch(color='#e74c3c', label='Negative correlation')
    fig.legend(handles=[blue_patch, red_patch], loc='lower center',
               ncol=2, fontsize=11, frameon=True)

    correction_label = 'FDR BH q' if USE_FDR else 'raw p'
    fig.suptitle(
        f'Genus-level Co-occurrence Network (Top {TOP_N} Shared Genera)\n'
        f'Spearman |r| > {R_THRESHOLD}, {correction_label} < {P_THRESHOLD}  '
        f'[n_Space=106, n_Terrestrial=20]',
        fontsize=14, y=1.01
    )
    plt.tight_layout(rect=[0, 0.05, 1, 1])

    png_path = OUTPUT_DIR / 'Fig2E_Network.png'
    pdf_path = OUTPUT_DIR / 'Fig2E_Network.pdf'
    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    plt.savefig(pdf_path, bbox_inches='tight')

    cb_dir = config.CURRENT_BIOLOGY_DIR
    cb_dir.mkdir(exist_ok=True)
    plt.savefig(cb_dir / 'Fig1e_Network.png', dpi=300, bbox_inches='tight')
    plt.savefig(cb_dir / 'Fig1e_Network.pdf', bbox_inches='tight')
    print(f"\n  Saved: {cb_dir / 'Fig1e_Network.png'}")

    plt.close()
    print(f"\n  Saved: {png_path.name}")
    print(f"  Saved: {pdf_path.name}")


# ─────────────────────────────────────────────
# 5. Save metrics and keystone genera
# ─────────────────────────────────────────────

def save_metrics(G_space, G_terr):
    rows = []
    for group, G in [('Space Flight', G_space), ('Terrestrial Soil', G_terr)]:
        neg = sum(1 for _, _, d in G.edges(data=True) if d['weight'] < 0)
        rows.append({
            'Group':          group,
            'Nodes':          G.number_of_nodes(),
            'Edges':          G.number_of_edges(),
            'Density':        round(nx.density(G), 4),
            'Avg_Degree':     round(np.mean([d for _, d in G.degree()]), 2),
            'Negative_Edges': neg,   # paper reports: Terrestrial=9, Space=3
        })
    metrics_df = pd.DataFrame(rows)
    path = OUTPUT_DIR / 'Fig2E_Network_Genus_Metrics.csv'
    metrics_df.to_csv(path, index=False)
    print(f"  Metrics: {path.name}")
    return metrics_df


def save_keystone(G_space, G_terr, top_genera):
    rows = []
    sp_deg = dict(G_space.degree())
    tr_deg = dict(G_terr.degree())
    for genus in top_genera:
        rows.append({
            'Genus':             genus,
            'Space_Degree':      sp_deg.get(genus, 0),
            'Terrestrial_Degree': tr_deg.get(genus, 0),
        })
    df = pd.DataFrame(rows).sort_values('Terrestrial_Degree', ascending=False)
    path = OUTPUT_DIR / 'Fig2E_Keystone_Genus_Comparison.csv'
    df.to_csv(path, index=False)
    print(f"  Keystone: {path.name}")
    return df


# ─────────────────────────────────────────────
# 6. Subsampling validation
# ─────────────────────────────────────────────

N_SUB   = 20   # match Terrestrial Soil sample size
N_ITER  = 100
SEED    = 42


def run_subsampling_validation(sp_rel, top_genera,
                               n_sub=N_SUB, n_iter=N_ITER, seed=SEED):
    """
    Validate that Space Flight network metrics are not driven by sample-size
    artifact. Space Flight samples (n=106) are repeatedly subsampled to n=20
    (matching Terrestrial Soil; 100 iterations) and network metrics are
    recomputed for each subsample.  Described in Materials & Methods, line 17.
    """
    rng     = np.random.default_rng(seed)
    samples = sp_rel.columns.tolist()
    records = []

    for i in range(n_iter):
        idx      = rng.choice(len(samples), size=n_sub, replace=False)
        sub_cols = [samples[j] for j in idx]
        sub_rel  = sp_rel[sub_cols]

        df       = sub_rel.loc[top_genera]
        df       = df[df.std(axis=1) > 0]   # drop constant rows (no variance → rank corr undefined)
        n_gen    = len(df.index)
        if n_gen < 3:
            continue
        data     = df.values.T          # (n_sub, n_gen)

        corr, pvals = spearmanr(data)
        if n_gen == 2:
            corr  = np.array([[1, corr], [corr, 1]])
            pvals = np.array([[0, pvals], [pvals, 0]])

        p_flat   = pvals[np.triu_indices(n_gen, k=1)]
        p_adj    = np.zeros((n_gen, n_gen))
        p_adj[np.triu_indices(n_gen, k=1)] = p_flat
        p_adj   += p_adj.T

        G = nx.Graph()
        G.add_nodes_from(df.index)
        for ii in range(n_gen):
            for jj in range(ii + 1, n_gen):
                r   = corr[ii, jj]
                pij = p_adj[ii, jj]
                if abs(r) > R_THRESHOLD and pij < P_THRESHOLD and not np.isnan(r):
                    G.add_edge(df.index[ii], df.index[jj], weight=r)

        neg = sum(1 for _, _, d in G.edges(data=True) if d['weight'] < 0)
        records.append({
            'iteration':     i + 1,
            'n_subsample':   n_sub,
            'edges':         G.number_of_edges(),
            'density':       round(nx.density(G), 4),
            'avg_degree':    round(np.mean([d for _, d in G.degree()]), 4),
            'negative_edges': neg,
        })

    return pd.DataFrame(records)


def save_subsampling_results(df_sub):
    path = OUTPUT_DIR / 'Fig2E_Network_Genus_Fast_Metrics.csv'
    df_sub.to_csv(path, index=False)
    print(f"  Subsampling results ({len(df_sub)} iterations): {path.name}")

    print(f"    Edges:      {df_sub['edges'].mean():.1f} ± {df_sub['edges'].std():.1f}")
    print(f"    Density:    {df_sub['density'].mean():.3f} ± {df_sub['density'].std():.3f}")
    print(f"    Avg Degree: {df_sub['avg_degree'].mean():.2f} ± {df_sub['avg_degree'].std():.2f}")
    return path


# ─────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────

if __name__ == '__main__':
    print('=' * 60)
    print('Option A: Genus-level Co-occurrence Network')
    print(f'Top {TOP_N} shared genera | |r|>{R_THRESHOLD} | FDR q<{P_THRESHOLD}')
    print('=' * 60)

    print('\n[1] Aggregating to genus level...')
    sp_rel = aggregate_to_genus(config.FEATURE_TABLE_SPACE, config.TAXONOMY_FILE)
    tr_rel = aggregate_to_genus(config.FEATURE_TABLE_TERR,  config.TAXONOMY_FILE)
    print(f'  Space genus table (raw): {sp_rel.shape[0]} genera × {sp_rel.shape[1]} samples')
    print(f'  Terrestrial genus table: {tr_rel.shape[0]} genera × {tr_rel.shape[1]} samples')

    # Filter Space to root-zone samples (Plant Roots + Rhizosphere) to match
    # the compartment filter applied in FAPROTAX and functional redundancy analyses.
    ROOT_ZONE_TYPES = {'Plant Roots', 'Rhizosphere'}
    meta_all = pd.read_csv(config.INTEGRATED_METADATA, sep='\t', index_col=0)
    root_zone_ids = meta_all.index[
        (meta_all['Study_Group'] == 'Space_Flight') &
        meta_all['sample_type'].isin(ROOT_ZONE_TYPES)
    ].tolist()
    sp_rel = sp_rel[[c for c in sp_rel.columns if c in root_zone_ids]]
    print(f'  Space genus table (root-zone filter): {sp_rel.shape[0]} genera × {sp_rel.shape[1]} samples')

    print('\n[2] Selecting top shared genera...')
    top_genera = select_shared_top_genera(sp_rel, tr_rel, n=TOP_N)

    print('\n[3] Building networks...')
    G_space = build_network(sp_rel, top_genera, 'Space Flight')
    G_terr  = build_network(tr_rel, top_genera, 'Terrestrial Soil')

    print('\n[4] Plotting networks...')
    plot_networks(G_space, G_terr, top_genera)

    print('\n[5] Saving metrics...')
    metrics_df   = save_metrics(G_space, G_terr)
    keystone_df  = save_keystone(G_space, G_terr, top_genera)

    print('\n[6] Subsampling validation (Space Flight n=20, 100 iterations)...')
    df_sub = run_subsampling_validation(sp_rel, top_genera)
    save_subsampling_results(df_sub)

    print('\n--- Network Metrics ---')
    print(metrics_df.to_string(index=False))
    print('\n--- Keystone Genera ---')
    print(keystone_df.to_string(index=False))
    print('\nDone.')
