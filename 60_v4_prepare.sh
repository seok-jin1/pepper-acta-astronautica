#!/usr/bin/env bash
# 60_v4_prepare.sh — bring every dataset onto the same 16S V4 window before denoising
#
# All datasets are cut to the 240 nt that immediately follow the 515F primer
# (E. coli ~534-773), so ASVs from different studies are directly comparable.
#   OSD-772 (ISS pepper) and OSD-267/269 (Veggie lettuce): 2x251 V4 reads whose R1 already starts
#     right after 515F (primers removed by the provider) -> R1 is used as is (truncated in 61).
#   PRJNA1145089 (field pepper): V3-V4 reads (341F-805R, primers removed). cutadapt locates the
#     515F site inside each read and keeps only the downstream V4 part.
#
# Inputs : $V4_RAW/{OSD-772,PRJNA1145089}/*.fastq.gz, lettuce raw folders (see LETTUCE_RAW)
# Outputs: $V4_WORK/input/<study>/<sample>.fastq.gz  (one single-end file per sample)
set -euo pipefail
V4_RAW=${V4_RAW:-$HOME/pepper-v4/raw}
V4_WORK=${V4_WORK:-$HOME/pepper-v4/work}
LETTUCE_RAW=${LETTUCE_RAW:-$HOME/pepper-aa-work/version-2_integrated/external/veggie_lettuce}
CUTADAPT=${CUTADAPT:-$HOME/miniconda3/envs/qiime2-amplicon/bin/cutadapt}
PRIMER_515F=GTGYCAGCMGCCGCGGTAA

mkdir -p "$V4_WORK"/input/{OSD-772,PRJNA1145089,OSD-267,OSD-269} "$V4_WORK/logs"

# OSD-772: GLDS-675_GAmplicon_<sample>_S<n>_L001_R1_raw.fastq.gz -> <sample>_S<n>_L001.fastq.gz
for f in "$V4_RAW"/OSD-772/*_R1_raw.fastq.gz; do
  s=$(basename "$f" _R1_raw.fastq.gz); s=${s#GLDS-675_GAmplicon_}
  ln -sf "$f" "$V4_WORK/input/OSD-772/$s.fastq.gz"
done

# Veggie lettuce (R1 only)
for study in OSD-267 OSD-269; do
  # provider raw reads only (primer-free R1); the GeneLab *_trimmed files are a different product
  for f in "$LETTUCE_RAW/$study"/raw/*_L001_R1_001.fastq.gz "$LETTUCE_RAW/$study"/raw/*_L001_R1_raw.fastq.gz; do
    [ -e "$f" ] || continue
    s=$(basename "$f"); s=${s%%_L001_R1_*}; s=${s#GLDS-*_Amplicon_}
    ln -sf "$f" "$V4_WORK/input/$study/$s.fastq.gz"
  done
done

# PRJNA1145089: keep the V4 part after the 515F site
: > "$V4_WORK/logs/cutadapt_PRJNA1145089.tsv"
for f in "$V4_RAW"/PRJNA1145089/*.fastq.gz; do
  s=$(basename "$f" .fastq.gz)
  "$CUTADAPT" -g "$PRIMER_515F" -e 0.1 --overlap 15 --discard-untrimmed -j 4 \
    -o "$V4_WORK/input/PRJNA1145089/$s.fastq.gz" "$f" --report=minimal \
    | tail -1 | sed "s/^/$s\t/" >> "$V4_WORK/logs/cutadapt_PRJNA1145089.tsv"
done

for d in "$V4_WORK"/input/*; do echo "$(basename "$d"): $(ls "$d" | wc -l) samples"; done
