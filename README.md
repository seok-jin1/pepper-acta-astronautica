# ISS-grown *Capsicum annuum* rhizosphere: ammonia-oxidation bottleneck

Analysis code and intermediate outputs for

> Seok-Jin Kang, Hongchul Shin. **A structural ammonia-oxidation bottleneck in ISS-grown *Capsicum annuum*: implications for BLSS nitrogen-cycle design.** Manuscript, 2026.

Department of Biotechnology, College of Life Sciences and Biotechnology, Korea University, Seoul, Republic of Korea

---

## Overview

We re-analyse the root-zone 16S rRNA amplicon data of chile pepper grown on the International Space Station (NASA OSD-772) against a terrestrial pepper rhizosphere reference (PRJNA1145089), and ask whether the ammonia- and nitrite-oxidising taxa that a bioregenerative life-support system (BLSS) depends on persist in flight. The analysis layers are:

- 16S rRNA V4 amplicon processing (QIIME 2 2024.10, DADA2, SILVA 138), rarefied to 1,000 reads
- Functional prediction with FAPROTAX 1.2.12 and PICRUSt2 2.5.2, combined into a consensus table
- Direct ASV-level detection of the canonical nitrifier genera and a Nitrospira phylogeny
- Genome-level N-cycle gene evidence from NCBI RefSeq / PGAP annotations of 13 rhizosphere genera
- ISS interior-surface cross-reference (OSD-65) and FEAST source tracking
- WGCNA modules, plus exploratory community-assembly and co-occurrence analyses
- A reliability block diagram (RBD) and a closed-budget nitrogen mass-balance ODE

All functional analyses use the **analytical root zone**: Space_Flight Plant Roots (n = 16) + Rhizosphere (n = 4) versus Terrestrial rhizosphere (n = 20).

### Key results and where they come from

| Result | Value | File | Script |
|---|---|---|---|
| Canonical nitrifier genera detected (composite of 4 genera) | ISS 0/20 vs Terrestrial 19/20, Cliff's δ = −0.95 | `version-2_integrated/Nitrifier_Direct_Results.csv` | `40` |
| FAPROTAX nitrification Log₂FC (Space/Terrestrial) | −13.08 | `version-2_integrated/FAPROTAX_group_comparison.csv` | `12` |
| Nitrification contributing genera per sample | 0.10 vs 4.80 | `version-2_integrated/Fig2C_Functional_Redundancy_Summary.csv` | `14` |
| N-cycle completeness | 76% vs 100% | `version-2_integrated/Fig2B_NCycle_Completeness_Summary.csv` | `16` |
| Nitrifier signal across rarefaction depths 1k / 5k / 10k | 0% ISS vs 95–100% Terrestrial | `version-2_integrated/SuppTable_S1c_rarefaction_sensitivity.csv` | `44` |
| Genera carrying *amoA* among 13 queried RefSeq genera | 2 (*Nitrosomonas*, *Nitrosospira*) | `acta_astronautica/gtdb_ncycle_matrix.csv` | `36` |
| Nitrifier genera in the ISS surface microbiome (OSD-65) | none detected | `acta_astronautica/iss_env_ncycle_genera.csv` | `39` |
| Nitrospira ASVs | 19 (15 canonical NOB, 4 comammox) | `acta_astronautica/nitrospira_tree/nitrospira_asv_clade.csv` | `50`* |
| Crop generations to 50% N-cycle failure (p = 0.10 per taxon per generation) | 8 | `acta_astronautica/simulation/rbd_failure_thresholds.csv` | `48` |
| First-harvest N with nitrifier inoculation on day 5 / day 10 (healthy: 914 mg N L⁻¹) | 914 / 512 mg N L⁻¹ | `acta_astronautica/simulation/ode_scenario_summary.csv` | `49` |

\* `50_nitrospira_phylogeny_reconstructed.py` is a reconstruction; see [Notes](#notes).

---

## Datasets

| Dataset | Description | Samples | Access |
|---|---|---|---|
| **OSD-772** | ISS-grown *C. annuum*, Advanced Plant Habitat; Illumina MiSeq, V4 | 106 Space_Flight (20 in the analytical root zone) | [NASA OSDR](https://osdr.nasa.gov/bio/repo/data/studies/OSD-772) |
| **PRJNA1145089** | Terrestrial soil-grown *C. annuum* rhizosphere | 20 | [NCBI SRA](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1145089) |
| **OSD-65** | ISS-MO Microbial Tracking-1A/1B interior surfaces (GeneLab-processed tables) | 196 | [NASA OSDR](https://osdr.nasa.gov/bio/repo/data/studies/OSD-65) |
| **OSD-267, OSD-269** | Veggie lettuce roots (VEG-01A, VEG-03A); exploratory only, not used in the manuscript | 4 flight + 4 ground | [NASA OSDR](https://osdr.nasa.gov/bio/repo/data/studies/OSD-267) |

The integrated metadata also contains `Ground_Seed` (n = 3, OSD-772 ground controls), which is excluded from all comparisons. Raw FASTQ files are not included; see [Reproducing the analysis](#reproducing-the-analysis).

---

## Manuscript figures

Panels were assembled into the final figures by hand; the files below are the script outputs each figure is built from.

| Manuscript | Content | Repository output | Script(s) |
|---|---|---|---|
| **Fig. 1** | N-cycle depletion pattern: Faith's PD, Bray–Curtis PCoA, FAPROTAX Log₂FC, FAPROTAX–PICRUSt2 consensus | `figures/aa/Fig2_failure_signature.{pdf,png}` | `figures/aa/Fig2_failure_signature.py` (inputs from `05`, `12`, `15`, `16`, `35`) |
| **Fig. 2** | Cross-filter PERMANOVA robustness | `version-2_integrated/SuppS2_PERMANOVA.{pdf,png}` | `15_permanova.R` |
| **Fig. 3** | Q1–Q4 quadrant dynamics of the ISS rhizosphere | `version-2_integrated/SuppS1d_Temporal.{pdf,png}` | `11_supp_temporal.py` |
| **Fig. 4** | Direct ASV-level detection of nitrifier genera | `figures/aa/Fig_Nitrifier_Direct.{pdf,png}` | `40_nitrifier_direct_analysis.py` |
| **Fig. 5** | Phylogenetic placement of Nitrospira ASVs | `figures/aa/SuppS1e_nitrospira_tree.{pdf,png}` | `50_nitrospira_phylogeny_reconstructed.py`* |
| **Fig. 6** | Functional redundancy, RefSeq/PGAP N-cycle gene matrix, PGP-trait entropy | `figures/aa/Fig3_root_cause.{pdf,png}` | `figures/aa/Fig3_root_cause.py` (inputs from `14`, `16`, `36`) |
| **Fig. 7** | ISS interior-surface reference lacks nitrifier genera | `figures/aa/Fig4_system_consequence.{pdf,png}` | `figures/aa/Fig4_system_consequence.py` (input from `39`) |
| **Fig. 8** | WGCNA modules of the pepper root zone | `acta_astronautica/pepper_wgcna/SuppS6_WGCNA_combined.{pdf,png}` | `46_pepper_wgcna.R`, `47_pepper_wgcna_visualizations.py` |
| **Fig. 9** | Reliability block diagram and N mass-balance ODE | `acta_astronautica/simulation/rbd_diagram`, `rbd_reliability_curves`, `ode_trajectories` (`.pdf`/`.png`) | `48_rbd_simulation.py`, `49_nitrogen_ode_simulation.py` |
| **Fig. 10** | Evidence synthesis and BLSS design considerations | `figures/aa/Fig6_mechanism_synthesis.{pdf,png}` | `figures/aa/Fig6_mechanism_synthesis.py` (wording was edited by hand in the final figure) |
| **Fig. S1** | Exploratory community-assembly proxy | `figures/aa/SuppFig_S4_assembly.{pdf,png}` | `figures/aa/SuppFig_S4_assembly.py` (input from `13`) |
| **Fig. S2** | Exploratory co-occurrence network and threshold sensitivity | `figures/aa/SuppFig_S5_network.{pdf,png}` | `figures/aa/SuppFig_S5_network.py` (inputs from `09`, `17`) |

Output file names follow the scripts' internal numbering, not the manuscript's. Other files in `figures/aa/` (`Fig1_system_context`, `Fig2d_consensus`, `Fig5_generalization`, `Fig7_systems_panel`) are not used in the manuscript.

**Tables.** Table 1 (ODE parameters) is defined in `49_nitrogen_ode_simulation.py`. The supplementary tables are `version-2_integrated/SuppTable_S1*.csv`: S1a sample-definition sensitivity (`41`), S1b non-independence diagnostics (`43`), S1c rarefaction sensitivity (`44`), S1d effect sizes (`45`); S1e–S1i are copies of the Nitrospira clade calls and the ODE outputs in `acta_astronautica/`.

---

## Repository layout

```
.
├── config.py                  # paths and parameters shared by all scripts
├── 00–05                      # data download, metadata, QIIME 2 processing
├── 06–17                      # core pepper analyses (diversity, FAPROTAX, PICRUSt2, assembly, networks)
├── 22–24                      # FEAST source tracking (+ analysis_expansion_utils.py)
├── 35–50                      # analyses added for this manuscript (see below)
├── figures/aa/                # main-figure scripts and outputs
├── acta_astronautica/         # outputs of 35–50 (consensus, RefSeq matrix, ISS env, WGCNA, simulations, tree, lettuce)
├── version-2_integrated/      # QIIME 2 exports, intermediate tables and figure outputs
├── current_biology/           # panel outputs from an earlier draft of the analysis (not used here)
├── envs/                      # conda environments
└── tests/                     # unit tests for 35, 36 and the lettuce harmonisation
```

---

## Reproducing the analysis

Run every script from the repository root (the directory containing `config.py`).

### Environments and external resources

```bash
conda env create -f envs/qiime2-amplicon.yml   # QIIME 2, Python and R analyses
conda env create -f envs/picrust2.yml          # PICRUSt2
conda env create -f envs/pepper-feast.yml      # FEAST (scripts 22–24)
conda env create -f envs/pepper-network.yml    # optional SpiecEasi networks
```

- **FAPROTAX 1.2.12**: download separately and point `FAPROTAX_DIR` at it (default: `../FAPROTAX_1.2.12`).
- **SILVA 138 classifiers** in `classifiers/`: `05_qiime2_pipeline.sh` expects `classifiers/silva-138-99-classifier.qza`; script 37 uses `SILVA_V4_CLASSIFIER` (default `classifiers/silva-138-99-515-806-nb-classifier.qza`).
- **NCBI Datasets CLI** for script 36; **MAFFT** and **FastTree 2** for script 50 (both ship with the QIIME 2 environment).
- R packages: `vegan`, `picante`, `WGCNA`, `ggplot2`, `dplyr`, `patchwork`.

`python config.py` checks that the required input files are in place.

### Order

1. **Data and processing**: `00_download_data.py` (PRJNA1145089), OSD-772 FASTQ into `external_data/OSD-772/`, then `01`–`04` (metadata and manifest), `05_qiime2_pipeline.sh` (QIIME 2) and PICRUSt2.
2. **Core analyses**: `07` PICRUSt2 pathways, `09` network, `10` Faith's PD, `11` quadrant dynamics, `12` FAPROTAX, `13` assembly null model, `14` functional redundancy, `15` PERMANOVA, `16` N-cycle completeness and PGP specificity, `17` network threshold sensitivity.
3. **Source tracking**: `22` → `23` (FEAST) → `24`.
4. **Manuscript analyses**:
   `35` PICRUSt2/FAPROTAX consensus ·
   `36` RefSeq/PGAP N-cycle gene matrix ·
   `39` ISS surface (OSD-65) cross-reference ·
   `40` direct nitrifier detection ·
   `41` sample-definition sensitivity ·
   `42` crew/built-environment taxa sensitivity ·
   `43` non-independence check ·
   `44` rarefaction sensitivity ·
   `45` effect sizes ·
   `46`–`47` WGCNA ·
   `48` RBD ·
   `49` N mass-balance ODE ·
   `50` Nitrospira phylogeny (`--compare` checks it against the saved outputs).
5. **Figures**: `python figures/aa/<script>.py`.
6. **Exploratory lettuce reanalysis** (not in the manuscript): `37_veggie_lettuce_reanalysis.sh`, then `38`–`38j`.

---

## Notes

- **Snapshot and later additions.** Tag `acta-astronautica-2026-04-23` is the analysis state on 23 April 2026. Commits after it add files that were used for the manuscript but sat outside that working copy: FEAST scripts `22`–`24` with `analysis_expansion_utils.py` and `envs/pepper-feast.yml`, scripts `40` and `42` with their outputs, the result tables (`*.csv`), and the path fixes in `config.py`, `37`, `38` and the lettuce manifest.
- **Reconstructed script.** The commands that built the Nitrospira tree were not saved. `50_nitrospira_phylogeny_reconstructed.py` re-implements them from the Methods and the kept intermediates; it reproduces the saved alignment, tree topology (Robinson–Foulds distance 0) and all 19 clade calls. Its reference set is `acta_astronautica/nitrospira_tree/reference_set.fasta`.
- **Lettuce reanalysis.** After rarefaction the lettuce comparison has 4 flight and 2 ground samples, and no nitrifier genus is detected in either group, so it cannot test the pepper finding. Results are in `acta_astronautica/lettuce/`; `pepper_lettuce_parallel_summary.csv` sets them beside pepper.
- **Not included**: manuscript sources, raw FASTQ, QIIME 2 artifacts (`*.qza`, `*.qzv`) and PICRUSt2 intermediates; the pipeline scripts regenerate them.

---

## Funding

This work was supported by the Google Cloud Research Credits program (award GCP477932969).

## Contact

Hongchul Shin (corresponding author), Department of Biotechnology, College of Life Sciences and Biotechnology, Korea University, Seoul 02841, Republic of Korea · saekomi5@korea.ac.kr
