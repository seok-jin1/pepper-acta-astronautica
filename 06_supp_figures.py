import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# 스타일 설정
sns.set(style="whitegrid")
plt.rcParams.update({'font.size': 12})

# --- S.Fig 1: DADA2 Denoising Stats ---
def plot_denoising_stats():
    print("Generating S.Fig 1 (Denoising Stats)...")
    try:
        stats = pd.read_csv('version-2_integrated/exported_denoising_stats/stats.tsv', sep='	', skiprows=[1], index_col=0)
        
        # 비율 계산 (최종 생존율)
        stats['Retention Rate (%)'] = (stats['non-chimeric'] / stats['input']) * 100
        stats = stats.sort_values('Retention Rate (%)')
        
        plt.figure(figsize=(15, 8))
        
        # Stacked Bar Plot을 위한 데이터 재구조화 (절대 수치 비교가 아닌 단계별 감소 추이)
        # 여기서는 Retention Rate만 간단히 시각화
        colors = plt.cm.viridis(np.linspace(0, 1, len(stats)))
        
        bars = plt.bar(stats.index, stats['Retention Rate (%)'], color='cadetblue', edgecolor='none')
        
        plt.title('S.Fig 1: Data Retention Rate After Denoising (DADA2)', fontsize=16)
        plt.xlabel('Samples', fontsize=12)
        plt.ylabel('Percentage of Reads Retained (%)', fontsize=12)
        plt.ylim(0, 100)
        plt.xticks([]) # 샘플 이름이 너무 많아 겹치므로 숨김 (또는 일부만 표시)
        plt.axhline(y=stats['Retention Rate (%)'].mean(), color='r', linestyle='--', label=f'Mean: {stats["Retention Rate (%)"].mean():.1f}%')
        plt.legend()
        
        plt.tight_layout()
        plt.savefig('version-2_integrated/Denoising_Stats.png')
        plt.close()
        print("Done.")
    except Exception as e:
        print(f"Failed to generate S.Fig 1: {e}")

# --- S.Fig 3: Plant DNA Contamination ---
def plot_plant_contamination():
    print("Generating S.Fig 3 (Plant DNA Contamination)...")
    try:
        # Load Taxonomy
        tax = pd.read_csv('version-2_integrated/exported_taxonomy/taxonomy.tsv', sep='	', index_col=0)
        
        # Load Raw Table
        table = pd.read_csv('version-2_integrated/exported_table_raw/feature-table.tsv', sep='	', skiprows=1, index_col=0)
        
        # Identify Contaminants
        mito_ids = tax[tax['Taxon'].str.contains('Mitochondria', case=False)].index
        chloro_ids = tax[tax['Taxon'].str.contains('Chloroplast', case=False)].index
        
        # Calculate Counts
        total_counts = table.sum()
        mito_counts = table.loc[table.index.intersection(mito_ids)].sum()
        chloro_counts = table.loc[table.index.intersection(chloro_ids)].sum()
        
        # Dataframe for plotting
        contam_df = pd.DataFrame({
            'Mitochondria': mito_counts,
            'Chloroplast': chloro_counts,
            'Microbial': total_counts - mito_counts - chloro_counts
        })
        
        # Normalize to 100%
        contam_df_pct = contam_df.div(contam_df.sum(axis=1), axis=0) * 100
        contam_df_pct = contam_df_pct.sort_values('Microbial') # 미생물 비율 낮은 순으로 정렬
        
        # Plotting
        plt.figure(figsize=(16, 8))
        contam_df_pct.plot(kind='bar', stacked=True, color=['#e74c3c', '#2ecc71', '#3498db'], width=0.8, figsize=(16,8))
        
        plt.title('S.Fig 3: Proportion of Plant DNA Contamination per Sample', fontsize=16)
        plt.ylabel('Proportion of Reads (%)')
        plt.xlabel('Samples')
        plt.xticks([]) # 샘플 이름 숨김
        plt.legend(title='Category', loc='upper left', bbox_to_anchor=(1, 1))
        
        plt.tight_layout()
        plt.savefig('version-2_integrated/Plant_DNA_Contamination.png')
        plt.close()
        print("Done.")

    except Exception as e:
        print(f"Failed to generate S.Fig 3: {e}")

# --- S.Fig 4: Taxa Bar Plot (Phylum Level) ---
def plot_taxa_bar():
    print("Generating S.Fig 4 (Taxa Bar Plot)...")
    try:
        # Load Taxonomy & Clean Table
        # Filter to Space_Flight and Terrestrial_Soil only (exclude Ground_Seed)
        meta = pd.read_csv('version-2_integrated/integrated_metadata.tsv', sep='	', index_col=0)
        keep_samples = meta.index[meta['Study_Group'].isin(['Space_Flight', 'Terrestrial_Soil'])].tolist()
        tax = pd.read_csv('version-2_integrated/exported_taxonomy/taxonomy.tsv', sep='	', index_col=0)
        table = pd.read_csv('version-2_integrated/exported_table_clean/feature-table.tsv', sep='	', skiprows=1, index_col=0)
        table = table[[c for c in table.columns if c in keep_samples]]
        
        # Map Feature IDs to Phylum
        def get_phylum(taxon_str):
            if pd.isna(taxon_str): return 'Unassigned'
            parts = taxon_str.split(';')
            for part in parts:
                if part.strip().startswith('p__'):
                    return part.strip().split('__')[1] if len(part.strip().split('__')) > 1 else 'Unassigned'
            return 'Unassigned'

        tax['Phylum'] = tax['Taxon'].apply(get_phylum)
        
        # Aggregate by Phylum
        table_tax = table.join(tax[['Phylum']], how='inner')
        phylum_counts = table_tax.groupby('Phylum').sum()
        
        # Filter Top 10 Phyla + Others
        top_n = 10
        total_abundance = phylum_counts.sum(axis=1).sort_values(ascending=False)
        top_phyla = total_abundance.head(top_n).index
        
        phylum_counts_top = phylum_counts.loc[top_phyla]
        others = phylum_counts.loc[~phylum_counts.index.isin(top_phyla)].sum()
        phylum_counts_top.loc['Others'] = others
        
        # Transpose for plotting (Samples as X, Phyla as stacked bars)
        phylum_pct = phylum_counts_top.T.div(phylum_counts_top.T.sum(axis=1), axis=0) * 100
        
        # Sort within each group by Proteobacteria proportion, then concatenate
        # (Space Flight first, then Terrestrial Soil) — matches paper Fig 1b layout
        if 'Proteobacteria' in phylum_pct.columns:
            space_ids_plot = [s for s in phylum_pct.index
                              if s in meta.index and meta.loc[s, 'Study_Group'] == 'Space_Flight']
            terr_ids_plot  = [s for s in phylum_pct.index
                              if s in meta.index and meta.loc[s, 'Study_Group'] == 'Terrestrial_Soil']
            sp_sorted   = phylum_pct.loc[space_ids_plot].sort_values('Proteobacteria', ascending=False)
            terr_sorted = phylum_pct.loc[terr_ids_plot].sort_values('Proteobacteria', ascending=False)
            phylum_pct  = pd.concat([sp_sorted, terr_sorted])
            
        # Plotting
        plt.figure(figsize=(16, 8))
        phylum_pct.plot(kind='bar', stacked=True, cmap='tab20', width=0.9, figsize=(16,8))
        
        plt.title('S.Fig 4: Relative Abundance of Top 10 Phyla', fontsize=16)
        plt.ylabel('Relative Abundance (%)')
        plt.xlabel('Samples')
        plt.xticks([])
        plt.legend(title='Phylum', bbox_to_anchor=(1.01, 1), loc='upper left')
        
        plt.tight_layout()
        plt.savefig('version-2_integrated/Fig1B_Taxa_Barplot.png')
        plt.savefig('version-2_integrated/Fig1B_Taxa_Barplot.pdf')
        import os; os.makedirs('manuscript_figures', exist_ok=True)
        plt.savefig('manuscript_figures/Fig1B_Taxa_Barplot.pdf')
        plt.close()
        print("Done.")
        
    except Exception as e:
        print(f"Failed to generate S.Fig 4: {e}")

if __name__ == "__main__":
    plot_denoising_stats()
    plot_plant_contamination()
    plot_taxa_bar()
