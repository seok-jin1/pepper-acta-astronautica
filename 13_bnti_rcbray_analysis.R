#!/usr/bin/env Rscript
# =============================================================================
# βNTI + RCbray Null Model Analysis
# ISS Chile Pepper Rhizosphere Microbiome Study
#
# Tests whether community assembly is driven by deterministic selection
# (|βNTI| > 2) vs. stochastic processes (RCbray).
#
# Interpretation (Chase et al. 2011 framework):
#   βNTI < -2  : Homogeneous selection (strong environmental filter)
#   βNTI > +2  : Variable selection
#   |βNTI| ≤ 2, RCbray > +0.95 : Dispersal limitation
#   |βNTI| ≤ 2, RCbray < -0.95 : Homogenizing dispersal
#   otherwise  : Undominated stochasticity
#
# Prerequisites:
#   - version-2_integrated/exported_table_clean/feature-table.tsv
#   - version-2_integrated/integrated_metadata.tsv
#   - version-2_integrated/exported_tree/tree.nwk
#     (export from QIIME2: see README Step 3b)
#
# Run in qiime2-amplicon environment:
#   conda activate qiime2-amplicon
#   Rscript 16_bnti_rcbray_analysis.R
#
# Output (version-2_integrated/):
#   - BNTI_RCbray_results.csv      : per-sample-pair βNTI and RCbray values
#   - BNTI_RCbray_summary.csv      : assembly process fractions per group
#   - BNTI_RCbray_plot.png/.pdf    : visualization
# =============================================================================

suppressPackageStartupMessages({
  library(picante)
  library(vegan)
  library(ggplot2)
  library(dplyr)
})

set.seed(42)

WORKDIR <- "version-2_integrated"
N_RAND  <- 999  # randomizations for βNTI and RCbray

cat("============================================================\n")
cat("βNTI + RCbray Null Model Analysis\n")
cat("============================================================\n\n")

# ---------------------------------------------------------------------------
# 1. Load data
# ---------------------------------------------------------------------------
cat("[1] Loading data...\n")

# Feature table (ASVs × samples)
otu_raw <- read.table(
  file.path(WORKDIR, "exported_table_clean/feature-table.tsv"),
  sep = "\t", header = TRUE, row.names = 1, skip = 1, check.names = FALSE, comment.char = ""
)
# picante expects samples as rows
otu <- t(otu_raw)
cat(sprintf("  Feature table: %d samples × %d ASVs\n", nrow(otu), ncol(otu)))

# Metadata
meta <- read.table(
  file.path(WORKDIR, "integrated_metadata.tsv"),
  sep = "\t", header = TRUE, row.names = 1
)
# Keep Space_Flight and Terrestrial_Soil only
meta <- meta[meta$Study_Group %in% c("Space_Flight", "Terrestrial_Soil"), , drop = FALSE]
otu  <- otu[rownames(otu) %in% rownames(meta), ]
meta <- meta[rownames(otu), , drop = FALSE]
cat(sprintf("  After filtering: %d samples (Space=%d, Terrestrial=%d)\n",
  nrow(otu),
  sum(meta$Study_Group == "Space_Flight"),
  sum(meta$Study_Group == "Terrestrial_Soil")))

# Phylogenetic tree
tree_path <- file.path(WORKDIR, "exported_tree/tree.nwk")
if (!file.exists(tree_path)) {
  cat("\nERROR: Phylogenetic tree not found at:", tree_path, "\n")
  cat("Export from QIIME2 first:\n")
  cat("  qiime tools export \\\n")
  cat("    --input-path version-2_integrated/rooted-tree.qza \\\n")
  cat("    --output-path version-2_integrated/exported_tree\n\n")
  quit(status = 1)
}
tree <- read.tree(tree_path)
cat(sprintf("  Tree loaded: %d tips\n", length(tree$tip.label)))

# Match tree to OTU table
otu_matched <- match.phylo.comm(tree, otu)
otu  <- otu_matched$comm
tree <- otu_matched$phy
cat(sprintf("  After tree matching: %d samples × %d ASVs\n", nrow(otu), ncol(otu)))

# ---------------------------------------------------------------------------
# 2. βNTI (Beta Nearest Taxon Index)
# ---------------------------------------------------------------------------
cat(sprintf("\n[2] Computing βNTI (%d randomizations, seed=42)...\n", N_RAND))
cat("    This may take 10–30 minutes depending on sample size.\n")

phydist <- cophenetic(tree)

bnti_result <- ses.mntd(
  samp    = otu,
  dis     = phydist,
  null.model = "taxa.labels",
  abundance.weighted = TRUE,
  runs    = N_RAND
)

# βNTI = (obs - mean_null) / sd_null × sqrt(2) — standard formula
bnti_values <- (bnti_result$mntd.obs - bnti_result$mntd.rand.mean) /
               bnti_result$mntd.rand.sd

cat("  Done.\n")

# ---------------------------------------------------------------------------
# 3. RCbray (Raup-Crick Bray-Curtis)
# ---------------------------------------------------------------------------
cat(sprintf("\n[3] Computing RCbray (%d randomizations)...\n", N_RAND))

# RCbray: fraction of null communities with Bray-Curtis distance < observed
bc_obs  <- as.matrix(vegdist(otu, method = "bray"))
rc_bray <- matrix(NA, nrow = nrow(otu), ncol = nrow(otu),
                  dimnames = list(rownames(otu), rownames(otu)))

for (i in 1:(nrow(otu) - 1)) {
  for (j in (i + 1):nrow(otu)) {
    pair_otu <- otu[c(i, j), , drop = FALSE]
    # Only use ASVs present in at least one of the two samples
    pair_otu <- pair_otu[, colSums(pair_otu) > 0, drop = FALSE]

    null_bc <- replicate(N_RAND, {
      null_comm <- pair_otu
      null_comm[1, ] <- sample(pair_otu[1, ])
      null_comm[2, ] <- sample(pair_otu[2, ])
      as.numeric(vegdist(null_comm, method = "bray"))
    })

    rc_bray[i, j] <- rc_bray[j, i] <-
      (sum(null_bc < bc_obs[i, j]) - sum(null_bc > bc_obs[i, j])) / N_RAND
  }
}
cat("  Done.\n")

# ---------------------------------------------------------------------------
# 4. Classify assembly processes
# ---------------------------------------------------------------------------
cat("\n[4] Classifying assembly processes...\n")

results <- data.frame(
  sample1    = character(),
  sample2    = character(),
  group1     = character(),
  group2     = character(),
  pair_type  = character(),
  bnti       = numeric(),
  rcbray     = numeric(),
  process    = character(),
  stringsAsFactors = FALSE
)

for (i in 1:(nrow(otu) - 1)) {
  for (j in (i + 1):nrow(otu)) {
    s1 <- rownames(otu)[i]; s2 <- rownames(otu)[j]
    g1 <- meta[s1, "Study_Group"]; g2 <- meta[s2, "Study_Group"]
    pair_type <- if (g1 == g2) g1 else "Between"

    b <- bnti_values[i]   # βNTI is per-sample, use average of pair
    b2 <- bnti_values[j]
    bnti_pair <- (b + b2) / 2
    rc <- rc_bray[i, j]

    process <- dplyr::case_when(
      bnti_pair < -2          ~ "Homogeneous selection",
      bnti_pair > +2          ~ "Variable selection",
      rc > +0.95              ~ "Dispersal limitation",
      rc < -0.95              ~ "Homogenizing dispersal",
      TRUE                    ~ "Undominated stochasticity"
    )

    results <- rbind(results, data.frame(
      sample1 = s1, sample2 = s2,
      group1 = g1, group2 = g2,
      pair_type = pair_type,
      bnti = bnti_pair, rcbray = rc,
      process = process,
      stringsAsFactors = FALSE
    ))
  }
}

# Summary by pair type
summary_df <- results %>%
  group_by(pair_type, process) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(pair_type) %>%
  mutate(fraction = n / sum(n)) %>%
  ungroup()

cat("\n--- Assembly Process Fractions ---\n")
print(as.data.frame(summary_df), row.names = FALSE)

# ---------------------------------------------------------------------------
# 5. Save results
# ---------------------------------------------------------------------------
cat("\n[5] Saving results...\n")

write.csv(results,    file.path(WORKDIR, "SuppS4_BNTI_RCbray_results.csv"),  row.names = FALSE)
write.csv(summary_df, file.path(WORKDIR, "SuppS4_BNTI_RCbray_summary.csv"), row.names = FALSE)
cat("  Saved: SuppS4_BNTI_RCbray_results.csv\n")
cat("  Saved: SuppS4_BNTI_RCbray_summary.csv\n")

# ---------------------------------------------------------------------------
# 6. Plot
# ---------------------------------------------------------------------------
cat("\n[6] Plotting...\n")

plot_data <- results[results$pair_type != "Between", ]
plot_data$pair_type <- factor(plot_data$pair_type,
  levels = c("Space_Flight", "Terrestrial_Soil"))

p <- ggplot(plot_data, aes(x = bnti, fill = pair_type)) +
  geom_histogram(binwidth = 0.5, alpha = 0.7, position = "identity", color = "white") +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_fill_manual(values = c("Space_Flight" = "#3498db", "Terrestrial_Soil" = "#2ecc71")) +
  labs(
    title = "βNTI Distribution: Space Flight vs. Terrestrial Soil",
    subtitle = sprintf("Dashed lines at βNTI = ±2 (deterministic selection threshold, %d randomizations)", N_RAND),
    x = "βNTI", y = "Number of sample pairs",
    fill = "Group"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "top")

ggsave(file.path(WORKDIR, "SuppS4_BNTI_RCbray.png"), p, width = 10, height = 6, dpi = 300)
ggsave(file.path(WORKDIR, "SuppS4_BNTI_RCbray.pdf"), p, width = 10, height = 6)
cat("  Saved: SuppS4_BNTI_RCbray.png\n")
cat("  Saved: SuppS4_BNTI_RCbray.pdf\n")

cat("\n============================================================\n")
cat("βNTI + RCbray analysis complete.\n")
cat("============================================================\n")
