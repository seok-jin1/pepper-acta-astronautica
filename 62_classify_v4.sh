#!/usr/bin/env bash
# 62_classify_v4.sh — SILVA 138 V4 (515F/806R, 240 nt) naive Bayes classifier + ASV classification
#
# Usage:  bash 62_classify_v4.sh train      # build the classifier once (~30-60 min, ~10 GB RAM)
#         bash 62_classify_v4.sh classify   # classify $V4_WORK/dada2/asv.fasta
#
# Inputs  : $V4_REF/silva-138-99-seqs.qza, $V4_REF/silva-138-99-tax.qza  (QIIME 2 data resources)
# Outputs : $V4_REF/silva-138-99-v4-240-nb-classifier.qza
#           $V4_WORK/taxonomy/taxonomy.tsv
set -euo pipefail
V4_REF=${V4_REF:-$HOME/pepper-v4/ref}
V4_WORK=${V4_WORK:-$HOME/pepper-v4/work}
QIIME=${QIIME:-$HOME/miniconda3/envs/qiime2-amplicon/bin/qiime}
CLF="$V4_REF/silva-138-99-v4-240-nb-classifier.qza"

case "${1:-}" in
  train)
    "$QIIME" feature-classifier extract-reads \
      --i-sequences "$V4_REF/silva-138-99-seqs.qza" \
      --p-f-primer GTGYCAGCMGCCGCGGTAA --p-r-primer GGACTACNVGGGTWTCTAAT \
      --p-trunc-len 240 --p-min-length 100 --p-n-jobs 4 \
      --o-reads "$V4_REF/silva-138-99-v4-240-seqs.qza"
    "$QIIME" feature-classifier fit-classifier-naive-bayes \
      --i-reference-reads "$V4_REF/silva-138-99-v4-240-seqs.qza" \
      --i-reference-taxonomy "$V4_REF/silva-138-99-tax.qza" \
      --o-classifier "$CLF"
    ;;
  classify)
    mkdir -p "$V4_WORK/taxonomy"
    "$QIIME" tools import --type 'FeatureData[Sequence]' \
      --input-path "$V4_WORK/dada2/asv.fasta" --output-path "$V4_WORK/taxonomy/asv-seqs.qza"
    "$QIIME" feature-classifier classify-sklearn --p-n-jobs 4 \
      --i-classifier "$CLF" --i-reads "$V4_WORK/taxonomy/asv-seqs.qza" \
      --o-classification "$V4_WORK/taxonomy/taxonomy.qza"
    "$QIIME" tools export --input-path "$V4_WORK/taxonomy/taxonomy.qza" --output-path "$V4_WORK/taxonomy"
    ;;
  *) echo "usage: $0 train|classify" >&2; exit 2 ;;
esac
