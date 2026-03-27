#!/usr/bin/env Rscript
# =============================================================================
# PERMANOVA + Subsampling Validation
# ISS Chile Pepper Rhizosphere Microbiome Study
#
# Required analyses:
#   1. PERMANOVA (adonis2) – Bray-Curtis, full dataset (ALL) and ROOT_ZONE filter
#   2. Subsampling validation – subsample Space to n=20 (1000 iterations) to
#      verify that the result is not driven by sample size imbalance
#
# Ecological rationale:
#   Beta-diversity analysis (PCoA, Fig 1B) shows visual separation, but without
#   PERMANOVA there is no formal test of whether groups differ in composition.
#   PERMANOVA H0: no difference in community composition between groups.
#   Subsampling validates that Space n=106 vs. Terrestrial n=20 imbalance does
#   not inflate R² or deflate p.
#
# Prerequisites:
#   - version-2_integrated/exported_table_clean/feature-table.tsv
#   - version-2_integrated/integrated_metadata.tsv
#
# Run:
#   conda activate qiime2-amplicon
#   Rscript 25_permanova.R
#
# Output (version-2_integrated/):
#   - PERMANOVA_results.csv       : full PERMANOVA table (F, R², p for each filter)
#   - PERMANOVA_subsampling.csv   : subsampling distribution summary
#   - Main_Fig_PERMANOVA.png/.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(vegan)
  library(ggplot2)
  library(dplyr)
})

set.seed(42)
WORKDIR   <- "version-2_integrated"
N_PERM    <- 999
N_BOOT    <- 1000
N_SUB     <- 20     # subsample Space_Flight down to n=20

cat("============================================================\n")
cat("PERMANOVA + Subsampling Validation\n")
cat("============================================================\n\n")

# ---------------------------------------------------------------------------
# 1. Load data
# ---------------------------------------------------------------------------
cat("[1] Loading data...\n")

otu_raw <- read.table(
  file.path(WORKDIR, "exported_table_clean/feature-table.tsv"),
  sep = "\t", header = TRUE, row.names = 1,
  skip = 1, check.names = FALSE, comment.char = ""
)
otu <- t(otu_raw)  # samples × ASVs

meta <- read.table(
  file.path(WORKDIR, "integrated_metadata.tsv"),
  sep = "\t", header = TRUE, row.names = 1
)
meta <- meta[meta$Study_Group %in% c("Space_Flight", "Terrestrial_Soil"), , drop = FALSE]
otu  <- otu[rownames(otu) %in% rownames(meta), ]
meta <- meta[rownames(otu), , drop = FALSE]

otu_rel <- otu / rowSums(otu)
cat(sprintf("  ALL filter: %d samples (Space=%d, Terrestrial=%d)\n",
  nrow(otu), sum(meta$Study_Group == "Space_Flight"),
  sum(meta$Study_Group == "Terrestrial_Soil")))

# ROOT_ZONE filter: Rhizosphere + Plant_Roots + Arcillite
root_zone_types <- c("Rhizosphere", "Plant_Roots", "Arcillite", "Plant Roots")
meta_rz <- meta[meta$Study_Group == "Terrestrial_Soil" |
                (meta$Study_Group == "Space_Flight" &
                   meta$sample_type %in% root_zone_types), , drop = FALSE]
otu_rz  <- otu_rel[rownames(meta_rz), ]
cat(sprintf("  ROOT_ZONE filter: %d samples (Space=%d, Terrestrial=%d)\n",
  nrow(meta_rz),
  sum(meta_rz$Study_Group == "Space_Flight"),
  sum(meta_rz$Study_Group == "Terrestrial_Soil")))

# ---------------------------------------------------------------------------
# 2. PERMANOVA – ALL filter
# ---------------------------------------------------------------------------
cat("\n[2] PERMANOVA (ALL filter)...\n")

bc_all  <- vegdist(otu_rel, method = "bray")
perm_all <- adonis2(bc_all ~ Study_Group, data = meta, permutations = N_PERM)
cat("  adonis2 result:\n")
print(perm_all)

F_all  <- perm_all[1, "F"]
R2_all <- perm_all[1, "R2"]
p_all  <- perm_all[1, "Pr(>F)"]

# ---------------------------------------------------------------------------
# 3. PERMANOVA – ROOT_ZONE filter
# ---------------------------------------------------------------------------
cat("\n[3] PERMANOVA (ROOT_ZONE filter)...\n")

bc_rz   <- vegdist(otu_rz, method = "bray")
meta_rz_df <- as.data.frame(meta_rz)
perm_rz <- adonis2(bc_rz ~ Study_Group, data = meta_rz_df, permutations = N_PERM)
cat("  adonis2 result:\n")
print(perm_rz)

F_rz  <- perm_rz[1, "F"]
R2_rz <- perm_rz[1, "R2"]
p_rz  <- perm_rz[1, "Pr(>F)"]

# ---------------------------------------------------------------------------
# 4. Subsampling validation (ALL filter, Space downsampled to n=20)
# ---------------------------------------------------------------------------
cat(sprintf("\n[4] Subsampling validation (n_boot=%d, n_sub=%d)...\n", N_BOOT, N_SUB))

sp_ids <- rownames(meta)[meta$Study_Group == "Space_Flight"]
tr_ids <- rownames(meta)[meta$Study_Group == "Terrestrial_Soil"]

boot_F  <- numeric(N_BOOT)
boot_R2 <- numeric(N_BOOT)
boot_p  <- numeric(N_BOOT)

for (i in seq_len(N_BOOT)) {
  sub_sp  <- sample(sp_ids, N_SUB, replace = FALSE)
  sub_ids <- c(sub_sp, tr_ids)
  sub_rel <- otu_rel[sub_ids, ]
  sub_meta <- meta[sub_ids, , drop = FALSE]
  sub_bc   <- vegdist(sub_rel, method = "bray")
  sub_perm <- adonis2(sub_bc ~ Study_Group, data = sub_meta, permutations = 199)
  boot_F[i]  <- sub_perm[1, "F"]
  boot_R2[i] <- sub_perm[1, "R2"]
  boot_p[i]  <- sub_perm[1, "Pr(>F)"]
  if (i %% 100 == 0) cat(sprintf("  iter %d/%d\n", i, N_BOOT))
}

boot_F  <- boot_F[!is.na(boot_F)]
boot_R2 <- boot_R2[!is.na(boot_R2)]
boot_p  <- boot_p[!is.na(boot_p)]
cat(sprintf("  Subsampling F:  median=%.2f, 95%% CI=[%.2f, %.2f]\n",
  median(boot_F), quantile(boot_F, 0.025), quantile(boot_F, 0.975)))
cat(sprintf("  Subsampling R²: median=%.4f, 95%% CI=[%.4f, %.4f]\n",
  median(boot_R2), quantile(boot_R2, 0.025), quantile(boot_R2, 0.975)))
cat(sprintf("  Fraction significant (p<0.05): %.1f%%\n",
  100 * mean(boot_p < 0.05)))

# ---------------------------------------------------------------------------
# 5. Save results
# ---------------------------------------------------------------------------
cat("\n[5] Saving results...\n")

perm_df <- data.frame(
  Filter   = c("ALL", "ROOT_ZONE"),
  n_Space  = c(sum(meta$Study_Group == "Space_Flight"),
               sum(meta_rz$Study_Group == "Space_Flight")),
  n_Terr   = c(sum(meta$Study_Group == "Terrestrial_Soil"),
               sum(meta_rz$Study_Group == "Terrestrial_Soil")),
  F_stat   = round(c(F_all,  F_rz),  3),
  R2       = round(c(R2_all, R2_rz), 4),
  p_value  = c(p_all, p_rz),
  stringsAsFactors = FALSE
)
write.csv(perm_df, file.path(WORKDIR, "SuppS2_PERMANOVA_results.csv"), row.names = FALSE)
cat("  Saved: SuppS2_PERMANOVA_results.csv\n")

sub_df <- data.frame(
  Metric = c("F_median", "F_CI_lo", "F_CI_hi",
             "R2_median", "R2_CI_lo", "R2_CI_hi",
             "Frac_sig"),
  Value  = round(c(median(boot_F), quantile(boot_F, 0.025), quantile(boot_F, 0.975),
                   median(boot_R2), quantile(boot_R2, 0.025), quantile(boot_R2, 0.975),
                   mean(boot_p < 0.05)), 4)
)
write.csv(sub_df, file.path(WORKDIR, "SuppS2_PERMANOVA_subsampling.csv"), row.names = FALSE)
cat("  Saved: SuppS2_PERMANOVA_subsampling.csv\n")

# ---------------------------------------------------------------------------
# 6. Plot
# ---------------------------------------------------------------------------
cat("\n[6] Plotting...\n")

sig_label <- function(p) {
  if (p < 0.001) "***"
  else if (p < 0.01) "**"
  else if (p < 0.05) "*"
  else "n.s."
}

# Panel A: PERMANOVA R² bar chart
perm_plot_df <- perm_df %>%
  mutate(Filter_label = c(
    sprintf("ALL\n(Space n=%d, Terr n=%d)", n_Space[1], n_Terr[1]),
    sprintf("ROOT_ZONE\n(Space n=%d, Terr n=%d)", n_Space[2], n_Terr[2])
  ),
  sig = sapply(p_value, sig_label),
  label = sprintf("F=%.1f\nR²=%.3f\n%s", F_stat, R2, sig))

p1 <- ggplot(perm_plot_df, aes(x = Filter_label, y = R2)) +
  geom_col(fill = c("#3498db", "#2ecc71"), alpha = 0.8, width = 0.5) +
  geom_text(aes(label = label), vjust = -0.3, size = 3.5) +
  scale_y_continuous(limits = c(0, max(perm_plot_df$R2) * 1.5),
                     labels = scales::percent_format()) +
  labs(title = "PERMANOVA (Bray-Curtis)",
       subtitle = sprintf("%d permutations", N_PERM),
       x = "Sample Filter", y = "R² (Variance Explained)") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

# Panel B: Subsampling F distribution
sub_plot_df <- data.frame(F_val = boot_F)
p2 <- ggplot(sub_plot_df, aes(x = F_val)) +
  geom_histogram(bins = 40, fill = "#3498db", alpha = 0.7, color = "white") +
  geom_vline(xintercept = median(boot_F), color = "navy", linewidth = 1.2, linetype = "dashed") +
  geom_vline(xintercept = F_all, color = "red", linewidth = 1.2) +
  annotate("text", x = F_all * 1.02, y = Inf, label = sprintf("Full F=%.1f", F_all),
           color = "red", hjust = 0, vjust = 1.5, size = 3.5) +
  annotate("text", x = median(boot_F) * 1.02, y = Inf,
           label = sprintf("Subsample\nmedian=%.1f", median(boot_F)),
           color = "navy", hjust = 0, vjust = 1.5, size = 3.5) +
  labs(title = sprintf("Subsampling Validation\n(Space downsampled to n=%d, %d iterations)", N_SUB, N_BOOT),
       x = "PERMANOVA F statistic", y = "Count") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

# Combine
if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  fig <- p1 + p2 + plot_annotation(
    title = "PERMANOVA: Community Composition Differs Between Space Flight and Terrestrial Soil",
    subtitle = sprintf("adonis2 (Bray-Curtis), %d permutations; subsampling validates size-balanced result",
                       N_PERM),
    theme = theme(plot.title = element_text(face = "bold", size = 12))
  )
} else if (requireNamespace("gridExtra", quietly = TRUE)) {
  library(gridExtra)
  fig <- gridExtra::arrangeGrob(
    grobs = list(p1, p2),
    ncol = 2
  )
} else {
  fig <- p1
}

ggsave(file.path(WORKDIR, "SuppS2_PERMANOVA.png"), fig,
       width = 10, height = 5, dpi = 300)
ggsave(file.path(WORKDIR, "SuppS2_PERMANOVA.pdf"), fig,
       width = 10, height = 5)
cat("  Saved: Supp_S2_PERMANOVA.png/.pdf\n")

cat("\n--- Summary ---\n")
print(perm_df, row.names = FALSE)

cat("\n============================================================\n")
cat("PERMANOVA analysis complete.\n")
cat("============================================================\n")
