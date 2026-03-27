"""
Generate QIIME2 PairedEndFastqManifestPhred33V2 manifest.

Required input files (create before running):
  osd_r1_list.txt   — absolute paths to OSD-772 R1 FASTQ files, one per line
  osd_r2_list.txt   — absolute paths to OSD-772 R2 FASTQ files, same order as R1
  external_list.txt — SRR accession IDs for PRJNA1145089, one per line
                      (expects split paired files: <acc>_1.fastq.gz / <acc>_2.fastq.gz
                       obtained via: fasterq-dump --split-files <accession>)

Output:
  version-2_integrated/manifest.tsv — 3-column PairedEnd QIIME2 manifest
"""

import os
import sys

EXPECTED_OSD_COUNT = 109  # 106 Space_Flight + 3 Ground_Control

manifest_path = "version-2_integrated/manifest.tsv"

# ── Validate OSD-772 R1/R2 lists before writing ──────────────────────────────
with open("osd_r1_list.txt") as f:
    r1_lines = [l.strip() for l in f if l.strip() and not l.startswith("#")]
with open("osd_r2_list.txt") as f:
    r2_lines = [l.strip() for l in f if l.strip() and not l.startswith("#")]

if len(r1_lines) != len(r2_lines):
    sys.exit(
        f"ERROR: osd_r1_list.txt has {len(r1_lines)} entries but "
        f"osd_r2_list.txt has {len(r2_lines)}. Lists must be the same length."
    )
if len(r1_lines) != EXPECTED_OSD_COUNT:
    print(
        f"WARNING: Expected {EXPECTED_OSD_COUNT} OSD-772 samples, "
        f"got {len(r1_lines)}. Check osd_r1_list.txt.example for the full list."
    )

pairing_errors = []
for i, (r1, r2) in enumerate(zip(r1_lines, r2_lines), 1):
    try:
        sid_r1 = os.path.basename(r1).split("GAmplicon_")[1].split("_R1_")[0]
        sid_r2 = os.path.basename(r2).split("GAmplicon_")[1].split("_R2_")[0]
        if sid_r1 != sid_r2:
            pairing_errors.append(f"  line {i}: R1 sample '{sid_r1}' != R2 sample '{sid_r2}'")
    except IndexError:
        pairing_errors.append(f"  line {i}: could not parse sample ID from '{r1}' or '{r2}'")

if pairing_errors:
    sys.exit("ERROR: R1/R2 pairing mismatch (files must be in the same sample order):\n" +
             "\n".join(pairing_errors))

print(f"  OSD-772: {len(r1_lines)} R1/R2 pairs validated OK")

# ── Write manifest ────────────────────────────────────────────────────────────
with open(manifest_path, "w") as f:
    # PairedEndFastqManifestPhred33V2 requires 3 columns
    f.write("sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n")

    # 1. OSD-772 (Space Flight & Ground Seed)
    # Filename example: GLDS-675_GAmplicon_Q1_Fruit1_S5_L001_R1_raw.fastq.gz
    for r1_path, r2_path in zip(r1_lines, r2_lines):
        filename = os.path.basename(r1_path)
        try:
            sample_id = filename.split("GAmplicon_")[1].split("_R1_")[0]
            sample_id = sample_id.replace(" ", "_")
            f.write(f"{sample_id}\t{r1_path}\t{r2_path}\n")
        except IndexError:
            print(f"Warning: Could not parse ID from {filename}")

    # 2. External Terrestrial Reference (PRJNA1145089)
    # Expects: external_data/PRJNA1145089/<accession>_1.fastq.gz (R1)
    #          external_data/PRJNA1145089/<accession>_2.fastq.gz (R2)
    with open("external_list.txt", "r") as ext_f:
        for line in ext_f:
            acc = line.strip()
            if not acc:
                continue
            base_dir = os.path.abspath(f"external_data/PRJNA1145089")
            r1 = os.path.join(base_dir, f"{acc}_1.fastq.gz")
            r2 = os.path.join(base_dir, f"{acc}_2.fastq.gz")
            sample_id = f"Terrestrial_{acc}"
            f.write(f"{sample_id}\t{r1}\t{r2}\n")

print(f"Manifest written to: {manifest_path}")
print("Verify R1/R2 file pairs exist before running 05_qiime2_pipeline.sh")