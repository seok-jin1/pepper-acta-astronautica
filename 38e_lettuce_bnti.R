#!/usr/bin/env Rscript
# =============================================================================
# 38e_lettuce_bnti.R
# βNTI + RCbray null-model analysis for lettuce (OSD-267 + OSD-269)
#
# Adapts 13_bnti_rcbray_analysis.R to the lettuce DADA2-generated pipeline.
# FLAGGED AS EXPLORATORY — n=4 per condition is below the threshold for
# stable null-model inference; reported for parallelism with pepper only.
#
# Prerequisites:
#   - version-2_integrated/external/veggie_lettuce/exports/feature-table.tsv
#   - version-2_integrated/external/veggie_lettuce/exported_tree/tree.nwk
#   - version-2_integrated/external/veggie_lettuce/metadata_harmonized.tsv
#
# Run:
#   conda activate qiime2-amplicon
#   Rscript 38e_lettuce_bnti.R
#
# Output (acta_astronautica/lettuce/):
#   - lettuce_bnti_results.csv     : per-pair βNTI + RCbray
#   - lettuce_bnti_summary.csv     : assembly process fractions
# =============================================================================

suppressPackageStartupMessages({
  library(picante)
  library(vegan)
  library(dplyr)
})

set.seed(42)

WORKDIR  <- "version-2_integrated/external/veggie_lettuce"
OUT_DIR  <- "acta_astronautica/lettuce"
N_RAND   <- 99

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("============================================================\n")
cat("38e: βNTI + RCbray for lettuce (EXPLORATORY)\n")
cat("============================================================\n\n")

# ---------------------------------------------------------------------------
# 1. Load data
# ---------------------------------------------------------------------------
cat("[1] Loading data...\n")

ft_path <- file.path(WORKDIR, "exports/feature-table.tsv")
otu_raw <- read.table(ft_path, sep = "\t", header = TRUE, row.names = 1,
                      skip = 1, check.names = FALSE, comment.char = "")
otu <- t(otu_raw)  # samples × ASVs
cat(sprintf("  Feature table: %d samples × %d ASVs\n", nrow(otu), ncol(otu)))

meta <- read.table(file.path(WORKDIR, "metadata_harmonized.tsv"),
                   sep = "\t", header = TRUE, row.names = 1)
otu  <- otu[rownames(otu) %in% rownames(meta), , drop = FALSE]
meta <- meta[rownames(otu), , drop = FALSE]

# Tree
tree_path <- file.path(WORKDIR, "exported_tree/tree.nwk")
if (!file.exists(tree_path)) {
  stop(sprintf("Tree not found at %s. Generate via `qiime phylogeny align-to-tree-mafft-fasttree`.", tree_path))
}
tree <- read.tree(tree_path)

# Match tree tips to table ASVs
common <- intersect(colnames(otu), tree$tip.label)
cat(sprintf("  Shared ASVs between table and tree: %d\n", length(common)))
otu <- otu[, common, drop = FALSE]
tree <- drop.tip(tree, setdiff(tree$tip.label, common))

cat(sprintf("  Samples: Space_Flight=%d, Terrestrial_Soil=%d\n",
            sum(meta$Study_Group == "Space_Flight"),
            sum(meta$Study_Group == "Terrestrial_Soil")))

# Drop zero-count ASVs
keep_asvs <- colSums(otu) > 0
otu <- otu[, keep_asvs]
tree <- drop.tip(tree, setdiff(tree$tip.label, colnames(otu)))

# ---------------------------------------------------------------------------
# 2. βNTI via picante comdistnt
# ---------------------------------------------------------------------------
cat("\n[2] Computing βNTI (n_rand=", N_RAND, ")...\n", sep = "")

phydist <- cophenetic(tree)
bMNTD_obs <- as.matrix(comdistnt(otu, phydist, abundance.weighted = TRUE))

bMNTD_rand <- array(NA, c(nrow(otu), nrow(otu), N_RAND))
for (i in seq_len(N_RAND)) {
  shuffled_names <- sample(colnames(phydist))
  ph <- phydist
  rownames(ph) <- shuffled_names
  colnames(ph) <- shuffled_names
  bMNTD_rand[, , i] <- as.matrix(comdistnt(otu, ph, abundance.weighted = TRUE))
  if (i %% 20 == 0) cat(sprintf("  rand %d/%d\n", i, N_RAND))
}

bNTI <- matrix(NA, nrow(otu), nrow(otu))
for (i in 1:(nrow(otu) - 1)) {
  for (j in (i + 1):nrow(otu)) {
    m <- mean(bMNTD_rand[i, j, ])
    s <- sd(bMNTD_rand[i, j, ])
    if (s > 0) {
      bNTI[i, j] <- (bMNTD_obs[i, j] - m) / s
      bNTI[j, i] <- bNTI[i, j]
    }
  }
}
rownames(bNTI) <- rownames(otu)
colnames(bNTI) <- rownames(otu)

# ---------------------------------------------------------------------------
# 3. RCbray via vegan raupcrick
# ---------------------------------------------------------------------------
cat("\n[3] Computing RCbray...\n")
rc <- raupcrick(as.data.frame(otu > 0), null = "r1", nsimul = 199, chase = TRUE)
rc <- as.matrix(rc)
# Rescale to -1..+1 (vegan returns 0..1 probability → shift)
rc <- 2 * (rc - 0.5)
rownames(rc) <- rownames(otu)
colnames(rc) <- rownames(otu)

# ---------------------------------------------------------------------------
# 4. Classify assembly processes (Stegen et al. 2013)
# ---------------------------------------------------------------------------
cat("\n[4] Classifying assembly processes...\n")

rows <- list()
samples <- rownames(otu)
for (i in 1:(length(samples) - 1)) {
  for (j in (i + 1):length(samples)) {
    b <- bNTI[i, j]
    r <- rc[i, j]
    grp_i <- meta[samples[i], "Study_Group"]
    grp_j <- meta[samples[j], "Study_Group"]
    pair_type <- if (grp_i == grp_j) grp_i else "cross"

    if (is.na(b)) {
      process <- "undefined"
    } else if (b < -2) {
      process <- "homogeneous_selection"
    } else if (b > 2) {
      process <- "variable_selection"
    } else if (!is.na(r) && r > 0.95) {
      process <- "dispersal_limitation"
    } else if (!is.na(r) && r < -0.95) {
      process <- "homogenizing_dispersal"
    } else {
      process <- "undominated_stochastic"
    }

    rows[[length(rows) + 1]] <- data.frame(
      sample_1 = samples[i], sample_2 = samples[j],
      group_1 = grp_i, group_2 = grp_j,
      pair_type = pair_type,
      bNTI = round(b, 3),
      RCbray = round(r, 3),
      process = process,
      stringsAsFactors = FALSE
    )
  }
}
result <- do.call(rbind, rows)
write.csv(result, file.path(OUT_DIR, "lettuce_bnti_results.csv"), row.names = FALSE)
cat(sprintf("  Saved: %s/lettuce_bnti_results.csv (%d pairs)\n", OUT_DIR, nrow(result)))

# ---------------------------------------------------------------------------
# 5. Summary by pair type
# ---------------------------------------------------------------------------
cat("\n[5] Summary by pair type...\n")
summary_df <- result %>%
  group_by(pair_type, process) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(pair_type) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  arrange(pair_type, desc(pct))

print(summary_df, n = 30)
write.csv(summary_df, file.path(OUT_DIR, "lettuce_bnti_summary.csv"), row.names = FALSE)
cat(sprintf("  Saved: %s/lettuce_bnti_summary.csv\n", OUT_DIR))

cat("\n============================================================\n")
cat("βNTI/RCbray analysis complete (EXPLORATORY — n=4 per condition).\n")
cat("============================================================\n")
