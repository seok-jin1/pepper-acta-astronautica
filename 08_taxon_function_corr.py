import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

# Import centralized configuration (Peer Review Phase 2)
import config

def load_and_aggregate_taxa():
    """
    Aggregate ASV abundance to genus level.

    Returns:
        pd.DataFrame: Genus-level relative abundance (%)
    """
    print("Aggregating taxa to Genus level...")

    # Load using config paths
    table = pd.read_csv(config.FEATURE_TABLE_CLEAN, sep='\t', skiprows=1, index_col=0)
    tax = pd.read_csv(config.TAXONOMY_FILE, sep='\t', index_col=0)
    
    # Extract Genus name
    def get_genus(taxon_str):
        if pd.isna(taxon_str): return 'Unassigned'
        parts = taxon_str.split(';')
        for part in reversed(parts):
            if 'g__' in part and len(part.split('__')) > 1:
                name = part.split('__')[1].strip()
                return name if name != "" else "Unassigned"
        return 'Unassigned'
    
    tax['Genus'] = tax['Taxon'].apply(get_genus)
    
    # Merge and Group
    table_tax = table.join(tax[['Genus']], how='inner')
    genus_table = table_tax.groupby('Genus').sum()
    
    # Relative Abundance (%)
    genus_rel = genus_table.div(genus_table.sum(axis=0), axis=1) * 100
    return genus_rel

def load_pathways():
    """
    Load PICRUSt2 pathway predictions.

    Returns:
        pd.DataFrame: Pathway relative abundance (%)
    """
    print("Loading functional pathways...")

    # Load using config path (Phase 2)
    df = pd.read_csv(config.PICRUST2_PATHWAYS, sep='\t', index_col=0, compression='gzip')

    # Relative Abundance (%)
    df_rel = df.div(df.sum(axis=0), axis=1) * 100
    return df_rel

def perform_correlation(taxa_df, path_df):
    print("Calculating Spearman correlations with FDR correction (Benjamini-Hochberg)...")
    common_samples = taxa_df.columns.intersection(path_df.columns)
    taxa_sub = taxa_df[common_samples]
    path_sub = path_df[common_samples]

    # Top 10 Genera
    top_genera = taxa_sub.mean(axis=1).sort_values(ascending=False).head(12).index
    top_genera = [g for g in top_genera if g != 'Unassigned'][:10]

    # Pathway selection
    keywords = ['siderophore', 'iron', 'stress', 'oxidative', 'detoxification', 'superoxide']
    interest_paths_set = set()
    for kw in keywords:
        matches = path_sub.index[path_sub.index.str.contains(kw, case=False)]
        interest_paths_set.update(matches.tolist())
    # Deterministic ordering (sorted alphabetically) so results are reproducible
    interest_paths = sorted(interest_paths_set)[:15]

    # Supplement with top-abundance pathways until exactly 15 are selected.
    # This guarantees the paper's 150-pair (10 genera × 15 pathways) correlation matrix
    # and ensures consistent results when keyword matches are fewer than 15.
    if len(interest_paths) < 15:
        top_paths = path_sub.mean(axis=1).sort_values(ascending=False).index.tolist()
        for p in top_paths:
            if p not in interest_paths_set:
                interest_paths.append(p)
            if len(interest_paths) == 15:
                break

    if not interest_paths:
        print("Warning: No pathways found. Check that PICRUSt2 output is not empty.")
        return pd.DataFrame()

    # Step 1: 모든 pairwise 검정 수행
    results = []
    for genus in top_genera:
        for path in interest_paths:
            r, p = spearmanr(taxa_sub.loc[genus], path_sub.loc[path])
            results.append({'genus': genus, 'pathway': path, 'r': r, 'p': p})
    results_df = pd.DataFrame(results)

    # Step 2: BH FDR 보정
    reject, q_values, _, _ = multipletests(results_df['p'], alpha=0.05, method='fdr_bh')
    results_df['q'] = q_values
    n_tests = len(results_df)
    print(f"Total tests: {n_tests} | Significant (uncorrected p<0.05): {(results_df['p']<0.05).sum()} | Significant (FDR q<0.05): {reject.sum()}")

    # Step 3: 유의하지 않으면 0으로 마스킹
    corr_matrix = pd.DataFrame(0.0, index=top_genera, columns=interest_paths)
    for _, row in results_df[reject].iterrows():
        corr_matrix.loc[row['genus'], row['pathway']] = row['r']

    return corr_matrix

def plot_heatmap(corr_df):
    print("Generating Heatmap (Main Fig 6)...")
    plt.figure(figsize=(14, 10))
    # 가독성을 위해 이름 단축
    corr_df.columns = [c[:50] + '...' if len(c) > 50 else c for c in corr_df.columns]
    
    sns.heatmap(corr_df, annot=True, cmap='RdBu_r', center=0, fmt='.2f',
                cbar_kws={'label': 'Spearman Rho (FDR q<0.05)'})
    plt.title('Taxon-Function Correlation: Driving Forces of Spaceflight Adaptation\n(Benjamini-Hochberg FDR corrected, q < 0.05)',
              fontsize=14, fontweight='bold')
    plt.xlabel('Metabolic Pathways (PICRUSt2)', fontsize=12)
    plt.ylabel('Microbial Genera', fontsize=12)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()

    # Save using config path (Phase 2)
    fig_path = config.get_figure_path("Main", "2D", "Taxon_Function_Correlation")
    plt.savefig(fig_path, dpi=300)
    manuscript_dir = config.BASE_DIR / "manuscript_figures"
    manuscript_dir.mkdir(exist_ok=True)
    plt.savefig(manuscript_dir / "Fig2D_Taxon_Function_Correlation.pdf")
    print(f"✅ Figure saved to: {fig_path}")

    # Current Biology supplementary figure output
    cb_dir = config.CURRENT_BIOLOGY_DIR
    cb_dir.mkdir(exist_ok=True)
    plt.savefig(cb_dir / 'SuppS2a_TaxonFunction.png', dpi=300)
    plt.savefig(cb_dir / 'SuppS2a_TaxonFunction.pdf')
    print(f"✅ Saved Current Biology: {cb_dir / 'SuppS2a_TaxonFunction.png'}")
    plt.close()

if __name__ == "__main__":
    try:
        taxa_rel = load_and_aggregate_taxa()
        path_rel = load_pathways()
        corr_df = perform_correlation(taxa_rel, path_rel)
        plot_heatmap(corr_df)
        print("Analysis completed successfully. Main Fig 6 generated.")
    except Exception as e:
        print(f"Error: {e}")
