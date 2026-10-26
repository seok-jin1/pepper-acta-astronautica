# Nitrifiers in the root zone of ISS-grown chile pepper

Analysis code for

> Seok-Jin Kang, Hongchul Shin. **Nitrifiers are nearly absent from the root zone of chile pepper grown on the International Space Station: implications for nitrogen management in bioregenerative life support.** Manuscript, 2026.

Department of Biotechnology, College of Life Sciences and Biotechnology, Korea University, Seoul, Republic of Korea

> **This branch (`reanalysis-v4`) holds the pipeline behind the current manuscript: scripts 60–67.**
> Scripts `00`–`50` on this branch belong to an earlier, superseded version of the study and are kept only for provenance. They are not used by the current manuscript; see [Superseded pipeline](#superseded-pipeline).

---

## What the analysis does

Public 16S rRNA gene amplicon datasets cover different regions of the gene, so their amplicon sequence variants (ASVs) start at different positions and cannot be compared directly. This pipeline trims every pepper and lettuce dataset to the same 240-nt window of the V4 region (*Escherichia coli* 534–773, immediately after the 515F primer), denoises each study separately, and merges the tables so that identical sequences become the same ASV. It then screens pre-specified nitrifier lineages and reports every absence with an explicit detection limit.

Two definitions are fixed in advance:

- **Canonical genera** — *Nitrosomonas*, *Nitrosospira*, *Nitrospira*. *Nitrobacter* is **not** assessable in this window: 3 of its 13 distinct V4 reference sequences in SILVA 138 are identical to *Bradyrhizobium*, *Afipia* or other Xanthobacteraceae references. Script `67` reports this check.
- **Nitrifier lineages** (broad, deliberately inclusive) — all Nitrosomonadaceae including uncultured groups (Ellin6067, MND1, IS-44, mle1-7, 966-1, DSSD61), *Nitrosococcus*, *Nitrospira*, *Nitrobacter*, *Nitrococcus*, *Ca.* Nitrotoga, *Nitrolancea*, and the archaeal class Nitrososphaeria.

A sensitivity set adds cultured AOA genera, and a refined set drops Group 1.1c and SILVA "*Nitrolancea*".

### Headline results

| Result | Value | Output file | Script |
|---|---|---|---|
| Canonical genera, whole ISS science carrier | 0/75 samples; 11.49 M reads; pooled RA ≤3.2 × 10⁻⁷ | `revision/table2.tsv` | `67` |
| Canonical genera, ISS root zone vs field | 0/20 (≤1.05 × 10⁻⁶) vs 20/20 (7.1 × 10⁻³) | `revision/table2.tsv` | `67` |
| Primary test, growth quadrants vs field cultivars | 0/4 vs 5/5, Fisher p = 0.008 | `revision/table3_tests.tsv` | `67` |
| Sensitivity, canonical genera plus cultured AOA | root zone 1/20; quadrants 1/4 vs 5/5, p = 0.048 | `revision/table3_tests.tsv` | `67` |
| Nitrifier lineages, ISS root zone | 27 reads in 4/20 samples; 7.7 × 10⁻⁶ | `revision/table2.tsv` | `67` |
| Field / ISS root zone, nitrifier lineages | 3,731-fold; cluster bootstrap lower bound 1,306-fold | `revision/cluster_bootstrap.tsv` | `67` |
| Arcillite vs root zone, same quadrants | higher in 4/4 quadrants; sign test p = 0.0625; pooled ratio 20.4 | `revision/within_mission_quadrants.tsv` | `67` |
| Veggie lettuce ground-control roots | 0/4 (two samples informative, >100,000 reads each) | `revision/table2.tsv` | `67` |
| ISS interior surfaces (OSD-65) | 0/34 libraries; ≤3.0 × 10⁻⁶ | `osd65/summary.tsv` | `65` |
| PH-04 fertiliser nitrogen | 1.54 g N per quadrant, 47% as ammonium | `revision/table4_nitrogen_budget.tsv` | `67` |

`revision/numbers.tsv` lists every value quoted in the manuscript next to the quantity it comes from.

---

## Datasets

| Dataset | Used for | Samples used | Access |
|---|---|---|---|
| **OSD-772** | ISS *C. annuum*, Advanced Plant Habitat (PH-04); the primary dataset | 75 science-carrier samples (20 root zone, 16 arcillite, 16 wick, 12 foam, 3 water, 8 swabs) plus 34 tissue and seed libraries screened only for nitrifier reads | [NASA OSDR](https://osdr.nasa.gov/bio/repo/data/studies/OSD-772) |
| **PRJNA1145089** | Field pepper rhizosphere; reference for an established soil community, not a ground control | 20 (5 cultivars × 4 plants) | [NCBI SRA](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1145089) |
| **OSD-267, OSD-269** | Veggie lettuce roots (VEG-01A, VEG-03A); hardware ground control | 4 flight + 4 ground-control root samples | [OSD-267](https://osdr.nasa.gov/bio/repo/data/studies/OSD-267), [OSD-269](https://osdr.nasa.gov/bio/repo/data/studies/OSD-269) |
| **OSD-65** | ISS interior surfaces; presence/absence context only, GeneLab-processed tables | 34 surface libraries (8 controls and 4 unlabelled libraries excluded) | [NASA OSDR](https://osdr.nasa.gov/bio/repo/data/studies/OSD-65) |

No new sequencing was performed.

---

## Getting the data

Raw reads are not stored in this repository. Download them into the layout the scripts expect:

```
$V4_RAW/
  OSD-772/        GLDS-675_GAmplicon_<sample>_S<n>_L001_R1_raw.fastq.gz   (R1 only is used)
  PRJNA1145089/   SRR*.fastq.gz                                          (20 merged-read runs)
$LETTUCE_RAW/
  OSD-267/raw/    GLDS-267_Amplicon_<sample>_S<n>_L001_R1_001.fastq.gz
  OSD-269/raw/    GLDS-269_Amplicon_<sample>_S<n>_L001_R1_raw.fastq.gz
$ISS_ENV_DIR/processed/
                  GLDS-65_GAmplicon_16S-counts.tsv
                  GLDS-65_GAmplicon_16S-taxonomy.tsv
$V4_REF/
                  silva-138-99-seqs.qza
                  silva-138-99-tax.qza
```

- **OSDR studies (OSD-772, OSD-267, OSD-269, OSD-65)** — download the raw amplicon files (or, for OSD-65, the processed count and taxonomy tables) from the study pages linked above, or through the [OSDR API](https://visualization.osdr.nasa.gov/biodata/api/). Use the `*_R1_raw` / `*_R1_001` files; the GeneLab `*_trimmed` products are a different processing stage and must not be mixed in.
- **PRJNA1145089** — fetch the 20 runs from ENA or SRA, for example with `00_download_data.py` (kept from the earlier pipeline) or any FASTQ fetcher.
- **SILVA 138 99%** — the QIIME 2 reference artifacts from the [QIIME 2 data resources](https://docs.qiime2.org/2024.10/data-resources/).

Set the environment variables before running anything; each has a default under `$HOME` (see the header of each script):

```bash
export V4_RAW=/path/to/raw            # raw FASTQ, per study
export V4_WORK=/path/to/work          # all intermediate and final outputs
export V4_REF=/path/to/silva          # SILVA artifacts and the trained classifier
export LETTUCE_RAW=/path/to/lettuce   # OSD-267 / OSD-269 raw folders
export ISS_ENV_DIR=/path/to/osd65     # OSD-65 processed tables
```

---

## Software

| Step | Environment | Key versions |
|---|---|---|
| `60`, `62` | `envs/qiime2-amplicon.yml` | QIIME 2 2024.10, cutadapt 4.9 |
| `61` | same environment (R) | R 4.3, DADA2 1.30 |
| `63`–`67` | any Python 3.12 with pandas and SciPy | Python 3.12.2, pandas 2.1.4, SciPy 1.16.0 |

Training the classifier in `62` needs roughly 10 GB of RAM and 30–60 minutes. The DADA2 step is the other long one; set `V4_THREADS` to the number of cores you want it to use.

---

## Running the pipeline

Run the scripts in order; each reads what the previous one wrote.

```bash
bash 60_v4_prepare.sh            # trim every dataset to the shared 240-nt V4 window
Rscript 61_dada2_v4.R            # denoise per study, merge, md5 ASV ids
bash 62_classify_v4.sh train     # build the SILVA 138 V4-240 classifier (once)
bash 62_classify_v4.sh classify  # classify the ASVs
python 63_v4_checkpoint.py       # sample groups, shared ASVs, first nitrifier screen
python 64_v4_detection_stats.py  # detection limits, rarefaction, cluster tests
python 65_v4_osd65.py            # ISS interior surfaces
python 66_v4_figures.py          # manuscript figures 1-4
python 67_v4_revision_stats.py   # definitions, sensitivity analyses, bootstrap,
                                 # nitrogen budget, and every manuscript table
```

### What each script writes

| Script | Writes to `$V4_WORK/` |
|---|---|
| `60_v4_prepare.sh` | `input/<study>/<sample>.fastq.gz`, `logs/cutadapt_PRJNA1145089.tsv` |
| `61_dada2_v4.R` | `dada2/{feature-table.tsv, asv.fasta, track.tsv, seqtab_*.rds}` |
| `62_classify_v4.sh` | `taxonomy/taxonomy.tsv`; classifier into `$V4_REF` |
| `63_v4_checkpoint.py` | `checkpoint/{sample_groups.tsv, nitrifier_by_sample.tsv, nitrifier_prevalence.tsv, shared_asvs.tsv}` |
| `64_v4_detection_stats.py` | `stats/{detection_limits, rarefied_prevalence, detection_power, cluster_tests, iss_compartments, iss_quadrant_pairs, iss_nitrifier_asvs, lettuce_controls}.tsv` |
| `65_v4_osd65.py` | `osd65/{libraries, nitrifier_by_library, summary, nitrifier_asvs}.tsv` |
| `66_v4_figures.py` | `figures/Fig1–Fig4 (.png, .pdf)`, `figures/table1_datasets.tsv` |
| `67_v4_revision_stats.py` | `revision/*.tsv`, `revision/nitrifier_asvs.fasta`, `revision/numbers.tsv` |

Outputs are regenerated by the scripts and are not tracked here.

---

## Where the manuscript's figures, tables and numbers come from

| Manuscript item | Output | Script |
|---|---|---|
| **Fig. 1** Datasets and the shared V4 window | `figures/Fig1_study_design_v4_window.{png,pdf}` | `66` |
| **Fig. 2** Nitrifier lineages by compartment | `figures/Fig2_nitrifiers_by_compartment.{png,pdf}` | `66` |
| **Fig. 3** Within the PH-04 carrier: quadrants and lineage composition | `figures/Fig3_within_mission.{png,pdf}` | `66` |
| **Fig. 4** Detection power in the ISS root zone | `figures/Fig4_detection_power.{png,pdf}` | `66` |
| **Table 1** Datasets | `revision/table1.tsv` | `67` |
| **Table 2** Detection with 95% intervals | `revision/table2.tsv` | `67` |
| **Table 3** Between-environment and within-mission tests | `revision/table3_tests.tsv` | `67` |
| **Table 4** Illustrative nitrogen quantities | `revision/table4_nitrogen_budget.tsv` | `67` |
| **Table 5** Design options | written by hand in the manuscript build script; the supporting citations are listed there | — |
| Every in-text number | `revision/numbers.tsv` | `67` |
| Supplementary per-sample table | `revision/supplementary_per_sample.tsv` | `67` |
| Supplementary nitrifier ASV sequences | `revision/nitrifier_asvs.fasta` | `67` |

Additional outputs that back statements in the text: `revision/cluster_bootstrap.tsv` (cluster-robust intervals and fold differences), `revision/within_mission_quadrants.tsv` (quadrant-level comparisons with the sample-level values reported as descriptive), `revision/power_v5.tsv` (detection power against both definitions), `revision/rarefaction_v5.tsv` (100 rarefactions at 20,000 and 1,000 reads), `revision/iss_nitrifier_nearest_field.tsv` (mismatches between each ISS nitrifier ASV and its nearest field ASV).

---

## Notes

- **Analysis plan.** The nitrifier lineages, the cluster units, the primary test and the within-mission hypothesis were fixed in a written plan before the harmonised data were generated. The same public datasets had been analysed earlier without harmonisation, so the plan was informed by that analysis; the manuscript states this and lists the deviations.
- **Denominator.** All relative abundances use reads assigned to Bacteria or Archaea after removal of mitochondria and chloroplasts.
- **Clusters.** ISS samples are grouped by growth quadrant (Q1–Q4) and field samples by cultivar, which is a conservative aggregation unit; the plot layout of the field study was not reported.
- **Detection limits are read-based.** Each OSD-772 library was amplified from two 1-ng aliquots of DNA that also contained plant DNA, so bounds below about 10⁻⁵ describe the sequencing rather than the sampled community. The manuscript says so.
- **No blanks.** OSD-772 contains no sequenced extraction or PCR blanks. Script `67` runs a cross-talk check instead, comparing the observed and expected counts of the dominant ISS nitrifier ASV in wick and swab libraries from the same sequencing runs.
- **Randomness.** Rarefaction and the cluster bootstrap in `67` are seeded (`SEED = 20261021`); `61` uses `set.seed(42)`.

## Superseded pipeline

Scripts `00`–`50` and the `version-2_integrated/`, `acta_astronautica/` and `figures/` outputs come from an earlier version of this study, which compared the two pepper datasets without harmonising the amplicon region and added FAPROTAX, PICRUSt2, a RefSeq gene panel, FEAST source tracking, WGCNA, a reliability block diagram and a nitrogen ODE. Internal review found that the two datasets shared no ASVs under that processing, which invalidated the cross-study comparisons; those layers were dropped rather than repaired. The files remain on `main` and under the tag `acta-astronautica-2026-04-23` for provenance. Do not mix their outputs with the results above.

## Data availability

OSD-772 ([10.26030/mjvk-v435](https://doi.org/10.26030/mjvk-v435)), OSD-267, OSD-269 and OSD-65 are in the NASA Open Science Data Repository; PRJNA1145089 is in the NCBI Sequence Read Archive. All are public; this repository adds no new sequence data.
