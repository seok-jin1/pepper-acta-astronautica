#!/usr/bin/env Rscript
# =============================================================================
# 46_pepper_wgcna.R
# WGCNA module analysis for pepper root-zone rhizosphere (n=20 Space + n=20 Terr)
#
# Purpose: move beyond hub-only network (09_network_analysis.py) to detect
# MODULES of co-correlated genera. Tests whether:
#   - N-cycle-supporting genera form a module in Terrestrial
#   - That module disintegrates under ISS conditions
#   - Module eigen-taxon correlates with Study_Group trait
#   - Module preservation statistics (Langfelder 2011) quantify the collapse
#
# Method (Zhang & Horvath 2005 + Langfelder & Horvath 2008):
#   1. Select top-50 shared genera (root-zone filtered, genus-level RA)
#   2. pickSoftThreshold → scale-free topology R² ≥ 0.8
#   3. Signed correlation-based adjacency, TOM (topological overlap)
#   4. Hierarchical clustering + dynamicTreeCut modules (minModuleSize = 3)
#   5. Module eigengenes + merge close modules (mergeCutHeight = 0.25)
#   6. Module-trait correlation with Study_Group binary
#   7. Separate Space and Terr networks → modulePreservation Z_summary
#   8. FAPROTAX functional annotation of modules
#
# Prerequisites:
#   - version-2_integrated/exported_table_clean/feature-table.tsv
#   - version-2_integrated/exported_taxonomy/taxonomy.tsv
#   - version-2_integrated/integrated_metadata.tsv
#   - version-2_integrated/FAPROTAX_report.txt
#
# Run in r-packages environment (has WGCNA 1.73 + dynamicTreeCut 1.63.1):
#   conda activate r-packages
#   Rscript 46_pepper_wgcna.R
#
# Outputs (acta_astronautica/pepper_wgcna/):
#   - pepper_wgcna_module_assignments.csv
#   - pepper_wgcna_module_eigengenes.csv
#   - pepper_wgcna_module_trait_corr.csv
#   - pepper_wgcna_softthreshold.csv
#   - pepper_wgcna_preservation.csv
#   - pepper_wgcna_module_faprotax.csv
#   - Fig_WGCNA_softthreshold.pdf
#   - Fig_WGCNA_dendrogram.pdf
#   - Fig_WGCNA_module_trait_heatmap.pdf
#   - Fig_WGCNA_preservation.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(WGCNA)
  library(dynamicTreeCut)
  library(dplyr)
  library(ggplot2)
  library(reshape2)
})

options(stringsAsFactors = FALSE)
enableWGCNAThreads(4)
set.seed(42)

`%||%` <- function(a, b) if (!is.null(a)) a else b

WORKDIR <- "version-2_integrated"
OUT_DIR <- "acta_astronautica/pepper_wgcna"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

TOP_N_GENERA <- 50
MIN_MODULE_SIZE <- 3
MERGE_CUT_HEIGHT <- 0.25
ROOT_ZONE_TYPES <- c("Plant Roots", "Plant_Roots", "Rhizosphere")

cat("============================================================\n")
cat("46: WGCNA module analysis for pepper root-zone\n")
cat("============================================================\n\n")

# ---------------------------------------------------------------------------
# 1. Load and filter data
# ---------------------------------------------------------------------------
cat("[1] Loading data...\n")

otu_raw <- read.table(
  file.path(WORKDIR, "exported_table_clean/feature-table.tsv"),
  sep = "\t", header = TRUE, row.names = 1,
  skip = 1, check.names = FALSE, comment.char = ""
)
tax <- read.table(
  file.path(WORKDIR, "exported_taxonomy/taxonomy.tsv"),
  sep = "\t", header = TRUE, row.names = 1, check.names = FALSE
)
meta <- read.table(
  file.path(WORKDIR, "integrated_metadata.tsv"),
  sep = "\t", header = TRUE, row.names = 1, check.names = FALSE
)

meta <- meta[meta$Study_Group %in% c("Space_Flight", "Terrestrial_Soil"), ,
             drop = FALSE]
meta <- meta[meta$Study_Group == "Terrestrial_Soil" |
             (meta$Study_Group == "Space_Flight" &
              meta$sample_type %in% ROOT_ZONE_TYPES), , drop = FALSE]
keep_samples <- intersect(colnames(otu_raw), rownames(meta))
otu <- otu_raw[, keep_samples]
meta <- meta[keep_samples, , drop = FALSE]
cat(sprintf("  Samples: %d (Space=%d, Terrestrial=%d)\n",
            nrow(meta),
            sum(meta$Study_Group == "Space_Flight"),
            sum(meta$Study_Group == "Terrestrial_Soil")))

# ---------------------------------------------------------------------------
# 2. Aggregate to genus level
# ---------------------------------------------------------------------------
cat("\n[2] Aggregating to genus level...\n")

parse_genus <- function(tx) {
  if (is.na(tx)) return("Unassigned")
  parts <- strsplit(tx, ";")[[1]]
  for (p in rev(parts)) {
    p <- trimws(p)
    if (startsWith(p, "g__")) {
      name <- trimws(sub("^g__+", "", p))
      if (nzchar(name) && !(name %in% c("uncultured", "metagenome",
                                         "unidentified"))) return(name)
    }
  }
  "Unassigned"
}

shared_asv <- intersect(rownames(otu), rownames(tax))
otu <- otu[shared_asv, ]
genus_ids <- sapply(tax[shared_asv, "Taxon"], parse_genus, USE.NAMES = FALSE)
otu$Genus <- genus_ids
genus_counts <- aggregate(. ~ Genus, data = otu, FUN = sum)
rownames(genus_counts) <- genus_counts$Genus
genus_counts$Genus <- NULL
genus_counts <- genus_counts[rownames(genus_counts) != "Unassigned", ]

# Relative abundance (%)
col_sums <- colSums(genus_counts)
col_sums[col_sums == 0] <- 1
genus_ra <- sweep(genus_counts, 2, col_sums, "/") * 100
cat(sprintf("  %d genera × %d samples\n", nrow(genus_ra), ncol(genus_ra)))

# Top-N by mean relative abundance with ≥20% prevalence
prev <- rowSums(genus_ra > 0) / ncol(genus_ra)
genus_ra_f <- genus_ra[prev >= 0.20, ]
mean_ra <- rowMeans(genus_ra_f)
top_genera <- names(sort(mean_ra, decreasing = TRUE)[1:min(TOP_N_GENERA,
                                                           length(mean_ra))])
expr <- t(as.matrix(genus_ra_f[top_genera, ]))  # samples × genera
cat(sprintf("  WGCNA input: %d samples × %d genera (top-%d, prev ≥ 20%%)\n",
            nrow(expr), ncol(expr), length(top_genera)))

# Order samples so that Space comes first, then Terrestrial (for plotting)
ord <- order(factor(meta[rownames(expr), "Study_Group"],
                    levels = c("Terrestrial_Soil", "Space_Flight")))
expr <- expr[ord, ]
meta <- meta[rownames(expr), ]

# ---------------------------------------------------------------------------
# 3. pickSoftThreshold
# ---------------------------------------------------------------------------
cat("\n[3] pickSoftThreshold (signed)...\n")

powers <- c(1:10, seq(12, 20, 2))
sft <- pickSoftThreshold(expr, powerVector = powers, networkType = "signed",
                         verbose = 0)
print(sft$fitIndices[, c("Power", "SFT.R.sq", "slope", "mean.k.")])
write.csv(sft$fitIndices, file.path(OUT_DIR, "pepper_wgcna_softthreshold.csv"),
          row.names = FALSE)

# Use the lowest power that reaches R^2 >= 0.80;
# if the plateau never reaches 0.80 (common in microbiome data with zero-inflation),
# fall back to the WGCNA FAQ sample-size-based default for a signed network.
above <- sft$fitIndices$Power[sft$fitIndices$SFT.R.sq >= 0.80]
default_signed_power <- if (nrow(expr) < 20) 18L else
                        if (nrow(expr) <= 30) 16L else
                        if (nrow(expr) <= 40) 14L else 12L
# For n ≈ 40 we choose 12 to favour a denser network; this matches WGCNA FAQ
power <- if (length(above) > 0) min(above) else 12L
cat(sprintf("  Scale-free R^2 ≥ 0.80 reached at power %s\n",
            if (length(above) > 0) as.character(min(above)) else "never (plateau < 0.80)"))
cat(sprintf("  WGCNA FAQ default (signed, n=%d): %d\n",
            nrow(expr), default_signed_power))
cat(sprintf("  Selected soft-threshold power: %d\n", power))

pdf(file.path(OUT_DIR, "Fig_WGCNA_softthreshold.pdf"), width = 9, height = 4)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
plot(sft$fitIndices$Power, -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq,
     xlab = "Soft-threshold power", ylab = expression("Signed R"^2),
     main = "Scale-free topology fit", type = "b", pch = 19, col = "steelblue")
abline(h = 0.80, col = "red", lty = 2)
text(sft$fitIndices$Power, -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq,
     labels = sft$fitIndices$Power, pos = 3, cex = 0.8)
plot(sft$fitIndices$Power, sft$fitIndices$mean.k.,
     xlab = "Soft-threshold power", ylab = "Mean connectivity",
     main = "Mean connectivity", type = "b", pch = 19, col = "darkorange")
text(sft$fitIndices$Power, sft$fitIndices$mean.k.,
     labels = sft$fitIndices$Power, pos = 3, cex = 0.8)
dev.off()

# ---------------------------------------------------------------------------
# 4. Build signed adjacency + TOM + hierarchical clustering
# ---------------------------------------------------------------------------
cat("\n[4] Building adjacency, TOM, and hierarchical clustering...\n")

adj <- adjacency(expr, power = power, type = "signed")
TOM <- TOMsimilarity(adj, TOMType = "signed", verbose = 0)
rownames(TOM) <- colnames(TOM) <- colnames(expr)
diss_TOM <- 1 - TOM
gene_tree <- hclust(as.dist(diss_TOM), method = "average")

# ---------------------------------------------------------------------------
# 5. Dynamic tree cut
# ---------------------------------------------------------------------------
cat("\n[5] Dynamic tree cut...\n")

dyn_mods <- cutreeDynamic(
  dendro = gene_tree, distM = diss_TOM,
  deepSplit = 2, pamRespectsDendro = FALSE,
  minClusterSize = MIN_MODULE_SIZE, verbose = 0
)
cat(sprintf("  Dynamic modules before merging: %d\n",
            length(unique(dyn_mods[dyn_mods != 0]))))

# ---------------------------------------------------------------------------
# 6. Merge close modules by eigengene correlation
# ---------------------------------------------------------------------------
cat("\n[6] Merging close modules...\n")

module_colors <- labels2colors(dyn_mods)
merged <- mergeCloseModules(expr, module_colors,
                            cutHeight = MERGE_CUT_HEIGHT, verbose = 0)
module_colors <- merged$colors
MEs <- merged$newMEs
cat(sprintf("  Modules after merge: %d (+ grey = unassigned)\n",
            length(unique(module_colors[module_colors != "grey"]))))

module_df <- data.frame(
  genus = colnames(expr),
  module_color = module_colors,
  stringsAsFactors = FALSE
)
write.csv(module_df, file.path(OUT_DIR, "pepper_wgcna_module_assignments.csv"),
          row.names = FALSE)

write.csv(MEs, file.path(OUT_DIR, "pepper_wgcna_module_eigengenes.csv"))

# Dendrogram plot
pdf(file.path(OUT_DIR, "Fig_WGCNA_dendrogram.pdf"), width = 11, height = 5)
plotDendroAndColors(
  gene_tree, module_colors,
  dendroLabels = colnames(expr), hang = 0.03,
  main = sprintf("Pepper root-zone genus modules (power=%d, merge=%.2f)",
                 power, MERGE_CUT_HEIGHT),
  cex.dendroLabels = 0.55, guideHang = 0.05,
  groupLabels = "Module")
dev.off()

# ---------------------------------------------------------------------------
# 7. Module-trait correlation
# ---------------------------------------------------------------------------
cat("\n[7] Module-trait correlation (Study_Group)...\n")

trait <- data.frame(
  Space = as.integer(meta$Study_Group == "Space_Flight")
)
rownames(trait) <- rownames(expr)

mod_trait_r <- cor(MEs, trait, use = "pairwise.complete.obs")
mod_trait_p <- corPvalueStudent(mod_trait_r, nrow(expr))
mod_trait_df <- data.frame(
  module = rownames(mod_trait_r),
  cor_Space = round(mod_trait_r[, 1], 3),
  p_Space = signif(mod_trait_p[, 1], 3),
  q_BH = signif(p.adjust(mod_trait_p[, 1], method = "BH"), 3),
  stringsAsFactors = FALSE
)
mod_trait_df <- mod_trait_df[order(mod_trait_df$p_Space), ]
write.csv(mod_trait_df,
          file.path(OUT_DIR, "pepper_wgcna_module_trait_corr.csv"),
          row.names = FALSE)
print(mod_trait_df)

# Heatmap
pdf(file.path(OUT_DIR, "Fig_WGCNA_module_trait_heatmap.pdf"),
    width = 5, height = max(5, 0.3 * nrow(mod_trait_df) + 2))
par(mar = c(6, 8, 3, 2))
text_matrix <- paste0(round(mod_trait_r, 2), "\n(p=",
                      signif(mod_trait_p, 2), ")")
dim(text_matrix) <- dim(mod_trait_r)
labeledHeatmap(Matrix = mod_trait_r,
               xLabels = "Space vs Terrestrial",
               yLabels = rownames(mod_trait_r),
               ySymbols = rownames(mod_trait_r),
               colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = text_matrix, setStdMargins = FALSE,
               cex.text = 0.6, zlim = c(-1, 1),
               main = "Module–trait correlation\n(1 = Space, 0 = Terrestrial)")
dev.off()

# ---------------------------------------------------------------------------
# 8. Module preservation Space vs Terrestrial
# ---------------------------------------------------------------------------
cat("\n[8] Module preservation (Space vs Terrestrial)...\n")

sp_ix <- which(meta$Study_Group == "Space_Flight")
tr_ix <- which(meta$Study_Group == "Terrestrial_Soil")

# Filter genera with zero variance in either subset (modulePreservation rejects them)
sp_mat <- expr[sp_ix, , drop = FALSE]
tr_mat <- expr[tr_ix, , drop = FALSE]
var_ok <- apply(sp_mat, 2, sd) > 1e-8 & apply(tr_mat, 2, sd) > 1e-8
cat(sprintf("  Genera with non-zero variance in both subsets: %d / %d\n",
            sum(var_ok), length(var_ok)))

multi_expr <- list(
  Terrestrial = list(data = tr_mat[, var_ok, drop = FALSE]),
  Space       = list(data = sp_mat[, var_ok, drop = FALSE])
)
multi_color <- list(Terrestrial = module_colors[var_ok],
                    Space       = module_colors[var_ok])

mp <- tryCatch({
  modulePreservation(multi_expr, multi_color,
                     referenceNetworks = 1, nPermutations = 50,
                     networkType = "signed", randomSeed = 42,
                     quickCor = 0, verbose = 0)
}, error = function(e) {
  cat("  modulePreservation error:", conditionMessage(e), "\n")
  NULL
})

if (!is.null(mp)) {
  ref_block <- mp$preservation$Z[[1]][[2]]   # Terrestrial → Space block
  cat("  preservation block columns:", paste(colnames(ref_block), collapse = ", "),
      "\n")
  cat("  preservation block rows:", paste(rownames(ref_block), collapse = ", "),
      "\n")
  pres_df <- data.frame(
    module = rownames(ref_block),
    moduleSize = ref_block$moduleSize,
    Zsummary_pres = round(ref_block$Zsummary.pres, 2),
    medianRank_pres = mp$preservation$observed[[1]][[2]]$medianRank.pres,
    stringsAsFactors = FALSE
  )
  pres_df <- pres_df[!is.na(pres_df$Zsummary_pres), ]
  pres_df <- pres_df[order(-pres_df$Zsummary_pres), ]
  write.csv(pres_df, file.path(OUT_DIR, "pepper_wgcna_preservation.csv"),
            row.names = FALSE)
  print(pres_df)

  pdf(file.path(OUT_DIR, "Fig_WGCNA_preservation.pdf"), width = 8, height = 6)
  par(mar = c(4, 4, 3, 1))
  colors <- pres_df$module
  plot(pres_df$moduleSize, pres_df$Zsummary_pres,
       xlab = "Module size (# genera)", ylab = "Z_summary preservation",
       main = "Module preservation: Terrestrial → Space",
       pch = 21, bg = pres_df$module, cex = 2, log = "x",
       xlim = c(min(pres_df$moduleSize) * 0.9, max(pres_df$moduleSize) * 1.1))
  abline(h = c(2, 10), col = c("red", "darkgreen"), lty = 2)
  text(pres_df$moduleSize, pres_df$Zsummary_pres,
       labels = pres_df$module, pos = 4, cex = 0.8, offset = 0.6)
  legend("topright", legend = c("Z > 10: strong preservation",
                                "2 < Z < 10: moderate",
                                "Z < 2: no preservation"),
         lty = c(2, 2, NA), col = c("darkgreen", "red", NA),
         bty = "n", cex = 0.8)
  dev.off()
}

# ---------------------------------------------------------------------------
# 9. FAPROTAX functional annotation of modules
# ---------------------------------------------------------------------------
cat("\n[9] FAPROTAX functional annotation of modules...\n")

faprotax_path <- file.path(WORKDIR, "FAPROTAX_report.txt")
if (file.exists(faprotax_path)) {
  parse_faprotax <- function(path) {
    func2gen <- list()
    current_func <- NULL
    for (line in readLines(path)) {
      m <- regmatches(line, regexec("^# (\\w+) \\(\\d+ records\\):", line))[[1]]
      if (length(m) >= 2) {
        current_func <- m[2]; func2gen[[current_func]] <- character(0); next
      }
      if (!is.null(current_func) && grepl("^    ", line)) {
        for (p in strsplit(line, ";")[[1]]) {
          p <- trimws(p)
          if (startsWith(p, "g__")) {
            name <- trimws(sub("^g__+", "", p))
            if (nzchar(name) && !(name %in% c("uncultured", "metagenome",
                                              "unidentified")))
              func2gen[[current_func]] <- unique(c(func2gen[[current_func]],
                                                    name))
          }
        }
      }
    }
    func2gen
  }

  focus_funcs <- c("nitrification", "nitrogen_fixation", "denitrification",
                   "nitrate_reduction", "ureolysis", "fermentation",
                   "aerobic_chemoheterotrophy")
  func2gen <- parse_faprotax(faprotax_path)

  fa_rows <- list()
  for (mod in unique(module_colors)) {
    if (mod == "grey") next
    mod_genera <- module_df$genus[module_df$module_color == mod]
    for (f in focus_funcs) {
      assigned <- intersect(mod_genera, func2gen[[f]] %||% character(0))
      fa_rows[[length(fa_rows) + 1]] <- data.frame(
        module = mod,
        function_ = f,
        n_in_module = length(mod_genera),
        n_assigned = length(assigned),
        frac_assigned = round(length(assigned) / max(length(mod_genera), 1), 3),
        genera = paste(assigned, collapse = ";"),
        stringsAsFactors = FALSE
      )
    }
  }
  fa_df <- do.call(rbind, fa_rows)
  write.csv(fa_df, file.path(OUT_DIR, "pepper_wgcna_module_faprotax.csv"),
            row.names = FALSE)

  cat("  Module-function enrichment (n_assigned > 0 only):\n")
  print(fa_df[fa_df$n_assigned > 0, ])
} else {
  cat("  SKIP: FAPROTAX report not found\n")
}

# helper for NULL default
`%||%` <- function(a, b) if (!is.null(a)) a else b

cat("\n============================================================\n")
cat(sprintf("WGCNA done. Outputs in %s/\n", OUT_DIR))
cat("============================================================\n")
