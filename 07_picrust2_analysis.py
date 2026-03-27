import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import gzip

# Import centralized configuration (Peer Review Phase 2)
import config

def load_data():
    """
    Load PICRUSt2 pathway predictions and metadata.

    Returns:
        tuple: (pathway_abundance_df, metadata_df)
    """
    # PICRUSt2 pathway abundance (using config paths - Phase 2)
    df = pd.read_csv(config.PICRUST2_PATHWAYS, sep='\t', index_col=0, compression='gzip')

    # Metadata
    meta = pd.read_csv(config.INTEGRATED_METADATA, sep='\t', index_col=0)

    return df, meta

def find_differential_pathways(df, meta):
    """
    Identify differentially abundant metabolic pathways between environments.

    Calculates Log2 Fold Change (Space/Terrestrial) for all MetaCyc pathways.

    Args:
        df: Pathway abundance matrix
        meta: Sample metadata with Study_Group column

    Returns:
        pd.DataFrame: Pathways with fold change statistics, sorted by Log2FC
    """
    print("Analyzing differential pathways...")
    # 공통 샘플만 유지
    common_samples = df.columns.intersection(meta.index)
    df = df[common_samples]
    meta = meta.loc[common_samples]
    
    # Relative abundance (%)
    df_rel = df.div(df.sum(axis=0), axis=1) * 100
    
    # Group별 평균 계산
    groups = meta['Study_Group'].unique()
    group_means = {}
    for g in groups:
        samples = meta[meta['Study_Group'] == g].index
        group_means[g] = df_rel[samples].mean(axis=1)
    
    means_df = pd.DataFrame(group_means)
    
    # Space vs Terrestrial 차이 계산 (Log2 Fold Change 유사)
    # 0 분모 방지를 위해 아주 작은 값 추가
    epsilon = 1e-6
    means_df['Diff'] = means_df['Space_Flight'] - means_df['Terrestrial_Soil']
    means_df['Log2FC'] = np.log2((means_df['Space_Flight'] + epsilon) / (means_df['Terrestrial_Soil'] + epsilon))
    
    return means_df.sort_values('Log2FC', ascending=False)

def plot_top_pathways(diff_df):
    """
    Generate bar plot of top differentially abundant pathways.

    Generates: Main Figure 4 (Metabolic Pathway Comparison)

    Args:
        diff_df: DataFrame with Log2FC values (output of find_differential_pathways)
    """
    print("Plotting top 20 upregulated pathways in Space...")

    # Get top pathways from config parameter
    n_top = config.PICRUST2_PARAMS.get("top_pathways", 30) // 2  # Split between up/down

    top_space = diff_df.head(n_top)
    top_terr = diff_df.tail(n_top)

    combined = pd.concat([top_space, top_terr])

    plt.figure(figsize=(12, 10))
    sns.barplot(x='Log2FC', y=combined.index, data=combined, palette='vlag')
    plt.title('Differentially Abundant Pathways (Space vs Terrestrial)')
    plt.xlabel('Log2 Fold Change (Space/Terrestrial)')
    plt.tight_layout()

    # Save using config path (Phase 2)
    fig_path = config.get_figure_path("Main", "1C", "Functional_Pathways")
    plt.savefig(fig_path, dpi=300)
    manuscript_dir = config.BASE_DIR / "manuscript_figures"
    manuscript_dir.mkdir(exist_ok=True)
    plt.savefig(manuscript_dir / "Fig1C_Functional_Pathways.pdf")
    print(f"✅ Figure saved to: {fig_path}")

    # Current Biology supplementary figure output
    cb_dir = config.CURRENT_BIOLOGY_DIR
    cb_dir.mkdir(exist_ok=True)
    plt.savefig(cb_dir / 'SuppS1b_PICRUSt2_Pathways.png', dpi=300)
    plt.savefig(cb_dir / 'SuppS1b_PICRUSt2_Pathways.pdf')
    print(f"✅ Saved Current Biology: {cb_dir / 'SuppS1b_PICRUSt2_Pathways.png'}")
    plt.close()

if __name__ == "__main__":
    try:
        df, meta = load_data()
        diff_df = find_differential_pathways(df, meta)
        plot_top_pathways(diff_df)
        
        # 특정 키워드 검색 (Siderophore, Stress)
        print("\n--- Key Pathway Search Results ---")
        keywords = ['siderophore', 'iron', 'stress', 'oxidative', 'antibiotic']
        for kw in keywords:
            matches = diff_df[diff_df.index.str.contains(kw, case=False)]
            if not matches.empty:
                print(f"\n[{kw.upper()}] 관련 경로:")
                print(matches[['Space_Flight', 'Terrestrial_Soil', 'Log2FC']].head(5))
                
        print("\nFunctional analysis completed.")
    except Exception as e:
        print(f"Analysis pending: Results file not ready yet. ({e})")
