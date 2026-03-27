# ISS Chile Pepper Rhizosphere Microbiome

**Can rhizosphere microbiomes support plants in space? Microbiomic evidence for functional network collapse in ISS-grown *Capsicum annuum***

Seok-Jin Kang, Hongchul Shin
*Current Biology* Correspondence (2026, in preparation) | Previously: *Plant Physiology Letters*
Department of Biotechnology, Korea University, Seoul, Republic of Korea

---

## Overview

This repository contains the complete bioinformatics pipeline used in the study. We perform a multi-level comparative analysis of the *Capsicum annuum* rhizosphere microbiome grown aboard the International Space Station (ISS) versus terrestrial soil controls, integrating:

- **16S rRNA amplicon sequencing** (QIIME2 + DADA2, V3–V4)
- **FAPROTAX functional trait annotation** (curated, v1.2.12; Louca et al. 2016)
- **Predicted functional pathways** (PICRUSt2 / MetaCyc)
- **Community assembly null models** (βNTI/RCbray; picante)
- **Genus-level co-occurrence network analysis** (Spearman |r| > 0.4, p < 0.05)
- **Network threshold sensitivity simulations**
- **Functional redundancy and N-cycle completeness quantification**

**Key finding:** Spaceflight restructures the rhizosphere toward a stochastically assembled, generalist-dominated community in which nitrogen-cycling capacity has structurally collapsed — nitrification abundance is >1000× lower in Space (FAPROTAX Log₂FC = −10.68, q < 0.001), nitrification functional redundancy drops from 4.80 to 0.20 contributing genera, and the keystone hub *Rhodanobacter* is absent from the Space co-occurrence network at all tested correlation thresholds.

---

## Datasets

| Dataset | Description | n (analysis) | Access |
|---------|-------------|--------------|--------|
| **OSD-772** | ISS-grown *C. annuum* rhizosphere (NASA VEGGIE) | 106 (Space Flight) | [NASA OSDR](https://osdr.nasa.gov/bio/repo/data/studies/OSD-772) |
| **PRJNA1145089** | Terrestrial soil-grown *C. annuum* rhizosphere | 20 (Terrestrial Soil) | [NCBI SRA](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1145089) |

> **Note on OSD-772 sample heterogeneity:** The 106 Space Flight samples span diverse sample types (rhizosphere, arcillite, wick, stem, fruit, leaf, seed, foam, swab, water). The root-zone filter (ROOT_ZONE n=36) and strict rhizosphere filter (STRICT n=4) yield the same direction of difference vs. Terrestrial (tested during revision).

> **Note on sample groups:** The integrated metadata contains three groups — `Space_Flight` (n=106), `Ground_Seed` (ISS experiment ground controls within OSD-772, n=3), and `Terrestrial_Soil` (n=20). The manuscript compares **Space_Flight vs. Terrestrial_Soil** only. `Ground_Seed` is excluded from all main analyses.

> **Raw data not included.** Download instructions are in [Step 0](#step-0-download-terrestrial-reference-data-prjna1145089) and [Step 1](#step-1-prepare-iss-data-osd-772) below.

---

## Repository Structure

```
.
├── README.md
├── CLAUDE.md                         # Claude Code project instructions
├── RESEARCH_NOTES.md                 # Post-rejection analysis plan and notes
├── config.py                         # Centralized paths and analysis parameters
├── .gitignore
│
├── envs/
│   ├── qiime2-amplicon.yml           # Core environment (QIIME2, Python, R scripts)
│   ├── picrust2.yml                  # Functional prediction (PICRUSt2)
│   └── pepper-network.yml            # R-based SpiecEasi network (optional)
│
├── osd772_metadata.tsv               # OSD-772 sample metadata (NASA OSDR, bundled)
│
├── 00_download_data.py               # Download PRJNA1145089 from NCBI ENA
├── 01_create_metadata.py             # Merge OSD-772 and terrestrial metadata
├── 02_generate_manifest.py           # Generate paired-end QIIME2 manifest
├── 03_sync_ids.py                    # Synchronize sample IDs across datasets
├── 04_fix_ids.py                     # Fix sample ID formatting issues
├── 05_qiime2_pipeline.sh             # Full QIIME2 amplicon processing pipeline
├── 06_supp_figures.py                # → Fig1B_Taxa_Barplot.png/.pdf  [Dropped from Current Biology]
├── 07_picrust2_analysis.py           # → SuppS1b_PICRUSt2_Pathways.png/.pdf  [Supp Fig S1b]
├── 08_taxon_function_corr.py         # → SuppS2a_TaxonFunction.png/.pdf  [Supp Fig S2a]
├── 09_network_analysis.py            # → Fig1e_Network.png/.pdf  [Fig 1e]
├── 10_supp_faith_pd.py               # → SuppS1a_FaithPD.png/.pdf  [Supp Fig S1a]
├── 11_supp_temporal.py               # → SuppS1d_Temporal.png/.pdf  [Supp Fig S1d]
├── 12_faprotax_analysis.py           # → Fig1b_FAPROTAX_Log2FC.png/.pdf  [Fig 1b]
├── 13_bnti_rcbray_analysis.R         # → SuppS2c_BNTI_RCbray.png/.pdf  [Supp Fig S2c]
├── 14_functional_redundancy.py       # → Fig1d_FunctionalRedundancy.png/.pdf  [Fig 1d]
├── 15_permanova.R                    # → SuppS1c_PERMANOVA.png/.pdf  [Supp Fig S1c]
├── 16_ncycle_pgp_specificity.py      # → Fig1c_NCycle_Heatmap.png/.pdf  [Fig 1c]
├── 17_network_threshold_sensitivity.py # → SuppS2b_ThresholdSensitivity.png/.pdf  [Supp Fig S2b]
├── rename_r_outputs.sh               # Copy R outputs to current_biology/ naming
│
├── current_biology/                   # Current Biology Correspondence outputs
│   ├── Fig1b_FAPROTAX_Log2FC.png/.pdf
│   ├── Fig1c_NCycle_Heatmap.png/.pdf
│   ├── Fig1d_FunctionalRedundancy.png/.pdf
│   ├── Fig1e_Network.png/.pdf
│   ├── SuppS1a_FaithPD.png/.pdf
│   ├── SuppS1b_PICRUSt2_Pathways.png/.pdf
│   ├── SuppS1c_PERMANOVA.png/.pdf
│   ├── SuppS1d_Temporal.png/.pdf
│   ├── SuppS2a_TaxonFunction.png/.pdf
│   ├── SuppS2b_ThresholdSensitivity.png/.pdf
│   └── SuppS2c_BNTI_RCbray.png/.pdf
│
└── version-2_integrated/             # All intermediate TSVs, CSVs, and figures
    # All figures now output to version-2_integrated/ (results/ is deprecated)
```

---

## Environment Setup

```bash
# Core analysis (QIIME2, Python scripts 06–17, R scripts via Rscript)
conda env create -f envs/qiime2-amplicon.yml
conda activate qiime2-amplicon

# Functional prediction (PICRUSt2 v2.5+)
conda env create -f envs/picrust2.yml

# R-based SpiecEasi network (optional, not used in main figures)
conda env create -f envs/pepper-network.yml
```

All Python scripts (06–17) and R scripts (13, 15) run in the `qiime2-amplicon` environment unless noted. Required R packages: `vegan`, `picante`, `ggplot2`, `dplyr`, `patchwork`.

---

## Full Pipeline — End-to-End Reproduction

> **All scripts must be run from the repository root** (the directory containing `config.py`):
> ```bash
> cd /path/to/pepper/
> ```
> Scripts use paths relative to this root (e.g., `version-2_integrated/`). Running from a different directory will cause file-not-found errors.

### Step 0: Download Terrestrial Reference Data (PRJNA1145089)

```bash
conda activate qiime2-amplicon
python 00_download_data.py
```

Downloads 16S amplicon FASTQ files from NCBI ENA to `external_data/PRJNA1145089/`. To obtain paired-end R1/R2 files required for QIIME2:

```bash
while read acc; do
  fasterq-dump --split-files --outdir external_data/PRJNA1145089 "$acc"
  gzip external_data/PRJNA1145089/${acc}_1.fastq external_data/PRJNA1145089/${acc}_2.fastq
done < external_list.txt
```

---

### Step 1: Prepare ISS Data (OSD-772)

Download raw paired-end FASTQ files from [NASA OSDR OSD-772](https://osdr.nasa.gov/bio/repo/data/studies/OSD-772) into `external_data/OSD-772/`:

```
external_data/OSD-772/GLDS-675_GAmplicon_<sample-id>_R1_raw.fastq.gz  (109 files)
external_data/OSD-772/GLDS-675_GAmplicon_<sample-id>_R2_raw.fastq.gz  (109 files)
```

`osd_r1_list.txt` and `osd_r2_list.txt` are pre-configured for this path and included in the repository.

---

### Step 2: Integrate Metadata

```bash
python 01_create_metadata.py
python 02_generate_manifest.py
python 03_sync_ids.py
python 04_fix_ids.py
```

Output: `version-2_integrated/integrated_metadata.tsv`
→ Three groups: `Space_Flight` (n=106), `Ground_Seed` (n=3), `Terrestrial_Soil` (n=20)

---

### Step 3: QIIME2 Amplicon Processing

```bash
conda activate qiime2-amplicon

# Download SILVA 138 V3-V4 reference (one-time):
mkdir -p classifiers
wget "https://data.qiime2.org/2024.5/common/silva-138-99-seqs.qza" \
  -O classifiers/silva-138-99-seqs.qza
wget "https://data.qiime2.org/2024.5/common/silva-138-99-tax.qza" \
  -O classifiers/silva-138-99-tax.qza

# Train Naive Bayes classifier (one-time, ~1 hr, ~16 GB RAM):
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads  classifiers/silva-138-99-seqs.qza \
  --i-reference-taxonomy classifiers/silva-138-99-tax.qza \
  --o-classifier classifiers/silva-138-99-classifier.qza

# Run full pipeline:
bash 05_qiime2_pipeline.sh
```

Pipeline steps (in order):
1. Import paired-end reads (PairedEndFastqManifestPhred33V2)
2. DADA2 denoising (trunc-len-f 240, trunc-len-r 200)
3. Export raw table → `exported_table_raw/`
4. Taxonomic classification (SILVA 138)
5. Filter mitochondria & chloroplast
6. Phylogenetic tree (MAFFT + FastTree)
7. Core diversity metrics (rarefaction depth: 1,000 reads)
8. Export tables as TSV (combined, Space-only, Terrestrial-only)
9. PICRUSt2 functional prediction (run separately in `picrust2` env)

---

### Steps 4–11: Figure Scripts (run in order)

```bash
python 06_supp_figures.py        # Dropped from Current Biology — phylum barplot
python 07_picrust2_analysis.py   # Supp Fig S1b — MetaCyc pathway Log₂FC
python 08_taxon_function_corr.py # Supp Fig S2a — taxon–function heatmap
python 09_network_analysis.py    # Fig 1e — co-occurrence network
python 10_supp_faith_pd.py       # Supp Fig S1a — Faith's PD

# Export rooted tree first (required for scripts 11 and 13):
qiime tools export \
  --input-path version-2_integrated/rooted-tree.qza \
  --output-path version-2_integrated/exported_tree

python 11_supp_temporal.py       # Supp Fig S1d — Q1–Q4 temporal

# Fig 1b: FAPROTAX functional trait annotation
# Requires FAPROTAX 1.2.12 at /home/laugh/pepper/FAPROTAX_1.2.12/
python 12_faprotax_analysis.py

# Supp Fig S2c: βNTI/RCbray community assembly null model (~10–30 min)
Rscript 13_bnti_rcbray_analysis.R

# Fig 1d: Functional redundancy
python 14_functional_redundancy.py

# Supp Fig S1c: PERMANOVA + subsampling validation (1000 iterations, ~5–10 min)
Rscript 15_permanova.R

# Fig 1c: N-cycle completeness + PGP specificity
python 16_ncycle_pgp_specificity.py

# Supp Fig S2b: Network threshold sensitivity (|r| = 0.3 / 0.4 / 0.5)
python 17_network_threshold_sensitivity.py
```

---

## Output Files

All files in `version-2_integrated/` unless otherwise noted. Filenames are prefixed with their figure number for easy identification.

### Main Figures

| Figure | Script | Figure output | Data output |
|--------|--------|---------------|-------------|
| **Fig 1b** | `12_faprotax_analysis.py` | `Fig1b_FAPROTAX_Log2FC.png/.pdf` | `FAPROTAX_functional_table.tsv`<br>`FAPROTAX_group_comparison.csv`<br>`FAPROTAX_key_functions.csv`<br>`FAPROTAX_report.txt` |
| **Fig 1c** | `16_ncycle_pgp_specificity.py` | `Fig1c_NCycle_Heatmap.png/.pdf` | `Fig1c_NCycle_Completeness_Results.csv`<br>`Fig1c_NCycle_Completeness_Summary.csv`<br>`Fig1c_PGP_Specificity_Results.csv`<br>`Fig1c_PGP_Specificity_Summary.csv` |
| **Fig 1d** | `14_functional_redundancy.py` | `Fig1d_FunctionalRedundancy.png/.pdf` | `Fig1d_Functional_Redundancy_Results.csv`<br>`Fig1d_Functional_Redundancy_Summary.csv` |
| **Fig 1e** | `09_network_analysis.py` | `Fig1e_Network.png/.pdf` | `Fig1e_Keystone_Genus_Comparison.csv`<br>`Fig1e_Network_Genus_Metrics.csv` |

### Supplementary Figures

| Supp Fig | Script | Figure output | Data output |
|----------|--------|---------------|-------------|
| **S1a** | `10_supp_faith_pd.py` | `SuppS1a_FaithPD.png/.pdf` | — |
| **S1b** | `07_picrust2_analysis.py` | `SuppS1b_PICRUSt2_Pathways.png/.pdf` | — |
| **S1c** | `15_permanova.R` | `SuppS1c_PERMANOVA.png/.pdf` | `SuppS1c_PERMANOVA_results.csv`<br>`SuppS1c_PERMANOVA_subsampling.csv` |
| **S1d** | `11_supp_temporal.py` | `SuppS1d_Temporal.png/.pdf` | — |
| **S2a** | `08_taxon_function_corr.py` | `SuppS2a_TaxonFunction.png/.pdf` | — |
| **S2b** | `17_network_threshold_sensitivity.py` | `SuppS2b_ThresholdSensitivity.png/.pdf` | `SuppS2b_ThresholdSensitivity_Results.csv` |
| **S2c** | `13_bnti_rcbray_analysis.R` | `SuppS2c_BNTI_RCbray.png/.pdf` | `SuppS2c_BNTI_RCbray_results.csv`<br>`SuppS2c_BNTI_RCbray_summary.csv` |

---

## Figure–Script Mapping (Current Biology Correspondence)

Scripts output individual panels to `current_biology/`. Panels are combined manually.

### Main Figure 1

| Panel | Script | Output (`current_biology/`) | Key result |
|-------|--------|---------------------------|------------|
| **1a** | — | Experimental design schematic | Conceptual diagram |
| **1b** | `12_faprotax_analysis.py` | `Fig1b_FAPROTAX_Log2FC.png/.pdf` | Nitrification Log₂FC = −10.68 |
| **1c** | `16_ncycle_pgp_specificity.py` | `Fig1c_NCycle_Heatmap.png/.pdf` | N-cycle completeness: Space 73% vs. Terr 100% |
| **1d** | `14_functional_redundancy.py` | `Fig1d_FunctionalRedundancy.png/.pdf` | Nitrification 0.20 vs. 4.80 genera |
| **1e** | `09_network_analysis.py` | `Fig1e_Network.png/.pdf` | Rhodanobacter hub collapse |

### Supplementary Figure S1

| Panel | Script | Output (`current_biology/`) | Key result |
|-------|--------|---------------------------|------------|
| **S1a** | `10_supp_faith_pd.py` | `SuppS1a_FaithPD.png/.pdf` | Faith's PD: Space 10.05 vs. Terr 24.98 |
| **S1b** | `07_picrust2_analysis.py` | `SuppS1b_PICRUSt2_Pathways.png/.pdf` | MetaCyc pathway Log₂FC |
| **S1c** | `15_permanova.R` + `rename_r_outputs.sh` | `SuppS1c_PERMANOVA.png/.pdf` | PERMANOVA R²=0.565, p=0.001 |
| **S1d** | `11_supp_temporal.py` | `SuppS1d_Temporal.png/.pdf` | Q1–Q4 temporal stability |

### Supplementary Figure S2

| Panel | Script | Output (`current_biology/`) | Key result |
|-------|--------|---------------------------|------------|
| **S2a** | `08_taxon_function_corr.py` | `SuppS2a_TaxonFunction.png/.pdf` | Taxon–function correlation heatmap |
| **S2b** | `17_network_threshold_sensitivity.py` | `SuppS2b_ThresholdSensitivity.png/.pdf` | Network metrics across |r| = 0.3/0.4/0.5 |
| **S2c** | `13_bnti_rcbray_analysis.R` + `rename_r_outputs.sh` | `SuppS2c_BNTI_RCbray.png/.pdf` | Terr 93.2% deterministic vs. Space <0.1% |

### Dropped from Current Biology version
- **Fig 1b** (phylum barplot, `06_supp_figures.py`) — removed
- **Fig 1c** (PICRUSt2, `07_picrust2_analysis.py`) — moved to Supp S1b

---

## Legacy Figure–Script Mapping (Plant Physiology Letters)

### Main Figures

| Figure | Panel | Script | Key result |
|--------|-------|--------|------------|
| **Fig 1a** | Overview schematic | — | Conceptual diagram |
| **Fig 1b** | Phylum barplot | `06_supp_figures.py` | Proteobacteria shift |
| **Fig 1c** | PICRUSt2 Log₂FC | `07_picrust2_analysis.py` | Nitrifier denitrification −11.5 |
| **Fig 1b** | FAPROTAX Log₂FC barplot | `12_faprotax_analysis.py` | Nitrification Log₂FC = −10.68 |
| **Fig 1c** | N-cycle step presence heatmap | `16_ncycle_pgp_specificity.py` | Space 73% vs. Terrestrial 100% |
| **Fig 1d** | Functional redundancy barplot | `14_functional_redundancy.py` | Nitrification 0.20 vs. 4.80 genera |
| **Supp S2a** | Taxon–function heatmap | `08_taxon_function_corr.py` | Rhodanobacter–nitrification correlation |
| **Fig 1e** | Co-occurrence network | `09_network_analysis.py` | Rhodanobacter hub collapse |

### Supplementary Figures

| Supp Fig | Script | Description |
|----------|--------|-------------|
| **S1** | `10_supp_faith_pd.py` | Faith's PD: Space 10.05 vs. Terrestrial 24.98 |
| **S2** | `15_permanova.R` | PERMANOVA bootstrap validation (1000 iterations) |
| **S3** | `11_supp_temporal.py` | Q1–Q4 temporal stability (Faith's PD, Shannon) |
| **S4** | `13_bnti_rcbray_analysis.R` | βNTI/RCbray distribution: Terrestrial 93.2% deterministic vs. Space <0.1% |
| **S5** | `17_network_threshold_sensitivity.py` | Network metrics across |r| = 0.3/0.4/0.5 |

---

## Key Results Summary

| Analysis | Space Flight | Terrestrial Soil | Statistic |
|----------|-------------|-----------------|-----------|
| PERMANOVA R² (ROOT_ZONE) | — | — | F=70.15, R²=0.565, p=0.001 |
| Shannon diversity | 2.66 ± — | 6.23 ± — | p < 0.001 |
| Nitrification (FAPROTAX Log₂FC) | — | — | −10.68, q < 0.001 |
| N-cycle completeness | 0.730 ± 0.193 | 1.000 ± 0.000 | p < 0.001 |
| Nitrification redundancy (genera) | 0.20 | 4.80 | q < 0.001 |
| Deterministic assembly (βNTI) | <0.1% | 93.2% | — |
| Rhodanobacter hub degree | 0 (all thresholds) | 6–8 | — |

---

## Reproducibility

- Core analysis parameters and paths are centralized in `config.py` (`RANDOM_SEED = 42`); scripts 00–06 and 10–11 use hardcoded relative paths
- Validate all required input files before running:
  ```bash
  python config.py
  ```
- MetaCyc pathway descriptions: `version-2_integrated/metacyc_pathways_info.txt.gz`

---

## Funding

This work was supported by the Google Cloud Research Credits program (award GCP477932969).

---

## Contact

Hongchul Shin (corresponding author)
Department of Biotechnology, College of Life Sciences and Biotechnology, Korea University, Seoul 02841, Republic of Korea
saekomi5@korea.ac.kr
