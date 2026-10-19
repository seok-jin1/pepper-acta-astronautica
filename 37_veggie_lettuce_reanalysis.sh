#!/usr/bin/env bash
set -euo pipefail

# Script 37: QIIME2 reanalysis of Veggie lettuce root microbiome
# OSD-267 (VEG-01A, SpaceX-3) + OSD-269 (VEG-03A, SpaceX-8)
# Supports Acta Astronautica Prediction P3 (cross-mission generalization)

# Paths come from config.py. Run from the repository root: the FASTQ paths in
# manifest.tsv are written relative to $PWD.
cd "$(dirname "${BASH_SOURCE[0]}")"
D=$(python -c 'import config; print(config.VEGGIE_LETTUCE_DIR)')
CLASSIFIER=$(python -c 'import config; print(config.SILVA_V4_CLASSIFIER)')
EXPORT_DIR="$D/exports"
mkdir -p "$EXPORT_DIR"

echo "[37] Step 1: Import paired-end reads..."
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path "$D/manifest.tsv" \
  --output-path "$D/demux.qza" \
  --input-format PairedEndFastqManifestPhred33V2

echo "[37] Step 2: Summarize demux..."
qiime demux summarize \
  --i-data "$D/demux.qza" \
  --o-visualization "$D/demux-summary.qzv"

echo "[37] Step 3: DADA2 denoise (trimmed input — minimal truncation)..."
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs "$D/demux.qza" \
  --p-trim-left-f 0 --p-trim-left-r 0 \
  --p-trunc-len-f 230 --p-trunc-len-r 200 \
  --o-table "$D/table.qza" \
  --o-representative-sequences "$D/rep-seqs.qza" \
  --o-denoising-stats "$D/denoising-stats.qza" \
  --p-n-threads 8

echo "[37] Step 4: Feature table summary..."
qiime feature-table summarize \
  --i-table "$D/table.qza" \
  --m-sample-metadata-file "$D/metadata_harmonized.tsv" \
  --o-visualization "$D/table-summary.qzv"

echo "[37] Step 5: Taxonomic classification (SILVA-138, 515-806)..."
qiime feature-classifier classify-sklearn \
  --i-classifier "$CLASSIFIER" \
  --i-reads "$D/rep-seqs.qza" \
  --o-classification "$D/taxonomy.qza" \
  --p-n-jobs 4

echo "[37] Step 6: Rarefy (sampling depth = 1000, matching pepper pipeline)..."
qiime feature-table rarefy \
  --i-table "$D/table.qza" \
  --p-sampling-depth 1000 \
  --o-rarefied-table "$D/table-rarefied.qza"

echo "[37] Step 7: Export..."
qiime tools export --input-path "$D/table-rarefied.qza" --output-path "$EXPORT_DIR/"
qiime tools export --input-path "$D/taxonomy.qza" --output-path "$EXPORT_DIR/"
biom convert -i "$EXPORT_DIR/feature-table.biom" -o "$EXPORT_DIR/feature-table.tsv" --to-tsv

echo "[37] Done. Exports in $EXPORT_DIR/"
echo "  feature-table.tsv: $(wc -l < "$EXPORT_DIR/feature-table.tsv") lines"
echo "  taxonomy.tsv:      $(wc -l < "$EXPORT_DIR/taxonomy.tsv") lines"
