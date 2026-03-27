#!/usr/bin/env bash
# =============================================================================
# QIIME2 Amplicon Processing Pipeline
# ISS Chile Pepper Rhizosphere Microbiome Study
#
# Usage:
#   conda activate qiime2-amplicon
#   bash 05_qiime2_pipeline.sh
#
# Prerequisites:
#   - FASTQ files downloaded (00_download_data.py)
#   - manifest.tsv generated (02_generate_manifest.py)
#   - SILVA 138 classifier downloaded (see README.md)
#
# Output directory: version-2_integrated/
# =============================================================================

set -euo pipefail

WORKDIR="version-2_integrated"
CLASSIFIER="classifiers/silva-138-99-classifier.qza"  # Trained TaxonomicClassifier (see README.md Step 3)

mkdir -p "${WORKDIR}"

echo "====================================================="
echo " Step 1: Import paired-end reads"
echo "====================================================="

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path "${WORKDIR}/manifest.tsv" \
  --output-path "${WORKDIR}/demux.qza" \
  --input-format PairedEndFastqManifestPhred33V2

qiime demux summarize \
  --i-data "${WORKDIR}/demux.qza" \
  --o-visualization "${WORKDIR}/demux-summary.qzv"

echo "====================================================="
echo " Step 2: DADA2 Denoising"
echo "====================================================="
# Trim lengths determined from demux quality plots
# V3-V4 region: forward ~240bp, reverse ~200bp

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs "${WORKDIR}/demux.qza" \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 240 \
  --p-trunc-len-r 200 \
  --p-n-threads 0 \
  --o-table "${WORKDIR}/table.qza" \
  --o-representative-sequences "${WORKDIR}/rep-seqs.qza" \
  --o-denoising-stats "${WORKDIR}/denoising-stats.qza"

qiime metadata tabulate \
  --m-input-file "${WORKDIR}/denoising-stats.qza" \
  --o-visualization "${WORKDIR}/denoising-stats.qzv"

echo "====================================================="
echo " Step 3: Taxonomic Classification (SILVA 138)"
echo "====================================================="

qiime feature-classifier classify-sklearn \
  --i-classifier "${CLASSIFIER}" \
  --i-reads "${WORKDIR}/rep-seqs.qza" \
  --p-n-jobs -1 \
  --o-classification "${WORKDIR}/taxonomy.qza"

echo "====================================================="
echo " Step 3b: Export Raw Table (needed by 06_supp_figures.py)"
echo "====================================================="
# Export the pre-filtering table for plant DNA contamination analysis (Supp Fig 3)

qiime tools export \
  --input-path "${WORKDIR}/table.qza" \
  --output-path "${WORKDIR}/exported_table_raw"
biom convert \
  -i "${WORKDIR}/exported_table_raw/feature-table.biom" \
  -o "${WORKDIR}/exported_table_raw/feature-table.tsv" \
  --to-tsv

echo "====================================================="
echo " Step 4: Filter Host-Derived Reads"
echo "====================================================="
# Remove mitochondria and chloroplast sequences

qiime taxa filter-table \
  --i-table "${WORKDIR}/table.qza" \
  --i-taxonomy "${WORKDIR}/taxonomy.qza" \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table "${WORKDIR}/table-no-mitochondria-no-chloroplast.qza"

qiime taxa filter-seqs \
  --i-sequences "${WORKDIR}/rep-seqs.qza" \
  --i-taxonomy "${WORKDIR}/taxonomy.qza" \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-sequences "${WORKDIR}/rep-seqs-no-mitochondria-no-chloroplast.qza"

echo "====================================================="
echo " Step 5: Phylogenetic Tree"
echo "====================================================="

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences "${WORKDIR}/rep-seqs-no-mitochondria-no-chloroplast.qza" \
  --o-alignment "${WORKDIR}/aligned-rep-seqs.qza" \
  --o-masked-alignment "${WORKDIR}/masked-aligned-rep-seqs.qza" \
  --o-tree "${WORKDIR}/unrooted-tree.qza" \
  --o-rooted-tree "${WORKDIR}/rooted-tree.qza"

echo "====================================================="
echo " Step 6: Feature Table Summary"
echo "====================================================="

qiime feature-table summarize \
  --i-table "${WORKDIR}/table-no-mitochondria-no-chloroplast.qza" \
  --m-sample-metadata-file "${WORKDIR}/integrated_metadata.tsv" \
  --o-visualization "${WORKDIR}/table-no-mit-no-chlo-summary.qzv"

# Export the QZV to get sample-frequency-detail.html (needed by 06_plot_depth.py)
mkdir -p "${WORKDIR}/table_summary_stats"
qiime tools export \
  --input-path "${WORKDIR}/table-no-mit-no-chlo-summary.qzv" \
  --output-path "${WORKDIR}/table_summary_stats"

echo "====================================================="
echo " Step 7: Core Diversity Metrics (rarefaction depth: 1000)"
echo "====================================================="
# Rarefaction depth = 1000 reads (retains ~86% of samples)
# Determined from Step 06_plot_depth.py

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny "${WORKDIR}/rooted-tree.qza" \
  --i-table "${WORKDIR}/table-no-mitochondria-no-chloroplast.qza" \
  --p-sampling-depth 1000 \
  --m-metadata-file "${WORKDIR}/integrated_metadata.tsv" \
  --output-dir "${WORKDIR}/core-metrics-results"

echo "====================================================="
echo " Step 8: Export Tables for Downstream Python Analysis"
echo "====================================================="

# Full combined table (ISS + Terrestrial)
qiime tools export \
  --input-path "${WORKDIR}/table-no-mitochondria-no-chloroplast.qza" \
  --output-path "${WORKDIR}/exported_table_clean"
biom convert \
  -i "${WORKDIR}/exported_table_clean/feature-table.biom" \
  -o "${WORKDIR}/exported_table_clean/feature-table.tsv" \
  --to-tsv

# Taxonomy
qiime tools export \
  --input-path "${WORKDIR}/taxonomy.qza" \
  --output-path "${WORKDIR}/exported_taxonomy"

# Rep seqs (for dominant ASV identification)
qiime tools export \
  --input-path "${WORKDIR}/rep-seqs-no-mitochondria-no-chloroplast.qza" \
  --output-path "${WORKDIR}/exported_sequences"

# Denoising stats
qiime tools export \
  --input-path "${WORKDIR}/denoising-stats.qza" \
  --output-path "${WORKDIR}/exported_denoising_stats"

# Diversity metrics — export all vectors needed by downstream scripts
# Shannon entropy (needed by 06_supp_figures.py indirectly via metadata)
qiime tools export \
  --input-path "${WORKDIR}/core-metrics-results/shannon_vector.qza" \
  --output-path "${WORKDIR}/exported_diversity"

# Faith's PD (needed by 10_supp_faith_pd.py and 11_supp_temporal.py)
qiime tools export \
  --input-path "${WORKDIR}/core-metrics-results/faith_pd_vector.qza" \
  --output-path "${WORKDIR}/exported_diversity"
# Exports as: exported_diversity/alpha-diversity.tsv

# Bray-Curtis PCoA ordination (needed by 11_supp_temporal.py)
qiime tools export \
  --input-path "${WORKDIR}/core-metrics-results/bray_curtis_pcoa_results.qza" \
  --output-path "${WORKDIR}/exported_diversity"
# Exports as: exported_diversity/ordination.txt

echo "====================================================="
echo " Step 9: Filter Group-Specific Tables"
echo "====================================================="

# Space Flight only
qiime feature-table filter-samples \
  --i-table "${WORKDIR}/table-no-mitochondria-no-chloroplast.qza" \
  --m-metadata-file "${WORKDIR}/integrated_metadata.tsv" \
  --p-where "[Study_Group]='Space_Flight'" \
  --o-filtered-table "${WORKDIR}/table-space-flight.qza"

qiime tools export \
  --input-path "${WORKDIR}/table-space-flight.qza" \
  --output-path "${WORKDIR}/exported_table_space"
biom convert \
  -i "${WORKDIR}/exported_table_space/feature-table.biom" \
  -o "${WORKDIR}/exported_table_space/feature-table.tsv" \
  --to-tsv

# Terrestrial Soil only
qiime feature-table filter-samples \
  --i-table "${WORKDIR}/table-no-mitochondria-no-chloroplast.qza" \
  --m-metadata-file "${WORKDIR}/integrated_metadata.tsv" \
  --p-where "[Study_Group]='Terrestrial_Soil'" \
  --o-filtered-table "${WORKDIR}/table-terrestrial-soil.qza"

qiime tools export \
  --input-path "${WORKDIR}/table-terrestrial-soil.qza" \
  --output-path "${WORKDIR}/exported_table_terrestrial"
biom convert \
  -i "${WORKDIR}/exported_table_terrestrial/feature-table.biom" \
  -o "${WORKDIR}/exported_table_terrestrial/feature-table.tsv" \
  --to-tsv

echo "====================================================="
echo " Step 10: PICRUSt2 Functional Prediction"
echo "====================================================="
# Activate picrust2 environment before running:
#   conda activate picrust2
# Or run separately as:
#   bash -c "conda run -n picrust2 picrust2_pipeline.py ..."

echo ""
echo "NOTE: Run PICRUSt2 in the 'picrust2' conda environment:"
echo ""
echo "  conda activate picrust2"
echo "  picrust2_pipeline.py \\"
echo "    -s ${WORKDIR}/exported_sequences/dna-sequences.fasta \\"
echo "    -i ${WORKDIR}/exported_table_clean/feature-table.biom \\"
echo "    -o ${WORKDIR}/picrust2_out \\"
echo "    -p 8"
echo ""
echo "  conda deactivate"
echo ""
echo "After PICRUSt2 completes, proceed with Python analysis scripts."

echo "====================================================="
echo " QIIME2 Pipeline Complete!"
echo " Next: python 06_supp_figures.py"
echo "====================================================="
