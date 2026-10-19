#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(utils)
})

outdir <- file.path("version-2_integrated", "source_tracking")
counts_path <- file.path(outdir, "genus_counts_space.tsv")
meta_path <- file.path(outdir, "primary_source_sink_metadata.tsv")
install_feast <- identical(Sys.getenv("INSTALL_FEAST", "0"), "1")

if (!file.exists(counts_path) || !file.exists(meta_path)) {
  stop("Required source-tracking inputs are missing. Run 22_build_source_tracking_inputs.py first.")
}

if (!requireNamespace("FEAST", quietly = TRUE)) {
  if (install_feast) {
    if (!requireNamespace("remotes", quietly = TRUE)) {
      install.packages("remotes", repos = "https://cloud.r-project.org")
    }
    remotes::install_github("cozygene/FEAST")
  }
}

if (!requireNamespace("FEAST", quietly = TRUE)) {
  stop(
    paste(
      "FEAST is not installed in the active R environment.",
      "Create the dedicated FEAST environment with:",
      "conda env create -f envs/pepper-feast.yml",
      "Then run:",
      "conda activate pepper-feast",
      "Rscript 23_run_feast.R",
      "If you still want in-place install in the active R environment, rerun with:",
      "INSTALL_FEAST=1 Rscript 23_run_feast.R",
      sep = "\n"
    )
  )
}

counts <- read.table(
  counts_path,
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  comment.char = ""
)

metadata <- read.table(
  meta_path,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

rownames(metadata) <- metadata$SampleID
metadata <- metadata[rownames(counts), , drop = FALSE]
feast_meta <- metadata[, c("Env", "SourceSink", "id")]
env_labels <- sort(unique(metadata$Env), decreasing = TRUE)

strip_env_suffix <- function(values, labels) {
  escaped <- vapply(labels, function(label) {
    gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", label)
  }, character(1))
  pattern <- paste0("_(", paste(escaped, collapse = "|"), ")$")
  sub(pattern, "", values, perl = TRUE)
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)

set.seed(42)
feast_result <- FEAST::FEAST(
  C = as.matrix(counts),
  metadata = feast_meta,
  different_sources_flag = 0,
  dir_path = paste0(outdir, "/"),
  outfile = "genus_primary"
)
invisible(feast_result)

setwd(old_wd)

feast_matrix_path <- file.path(outdir, "genus_primary_source_contributions_matrix.txt")
if (!file.exists(feast_matrix_path)) {
  stop("FEAST completed but did not write the expected source contribution matrix.")
}

feast_matrix <- read.table(
  feast_matrix_path,
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  comment.char = "",
  quote = "\""
)

source_ids <- colnames(feast_matrix)
source_env_map <- setNames(rep(NA_character_, length(source_ids)), source_ids)
source_env_map["Unknown"] <- "unknown"
source_sample_ids <- strip_env_suffix(setdiff(source_ids, "Unknown"), env_labels)
source_env_map[setdiff(source_ids, "Unknown")] <- metadata[source_sample_ids, "Env"]
if (any(is.na(source_env_map))) {
  stop("Failed to map one or more FEAST source columns back to source environments.")
}

sink_ids <- strip_env_suffix(rownames(feast_matrix), env_labels)
feast_matrix_by_env <- sapply(unique(source_env_map), function(source_env) {
  rowSums(feast_matrix[, source_env_map == source_env, drop = FALSE])
})
feast_matrix_by_env <- as.data.frame(feast_matrix_by_env, check.names = FALSE)
rownames(feast_matrix_by_env) <- sink_ids

write.table(
  feast_matrix_by_env,
  file = file.path(outdir, "feast_contributions_matrix.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

sink_metadata <- metadata[metadata$SourceSink == "Sink", c("SampleID", "sample_type", "location")]
sink_metadata$SampleID <- rownames(sink_metadata)

long_df <- as.data.frame(as.table(as.matrix(feast_matrix_by_env)), stringsAsFactors = FALSE)
colnames(long_df) <- c("SampleID", "Source", "Contribution")
long_df <- merge(long_df, sink_metadata, by = "SampleID", all.x = TRUE)

write.table(
  long_df,
  file = file.path(outdir, "feast_source_contributions.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

summary_df <- aggregate(
  Contribution ~ sample_type + location + Source,
  data = long_df,
  FUN = mean
)

write.table(
  summary_df,
  file = file.path(outdir, "feast_source_contributions_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Saved:\n")
cat("  ", file.path(outdir, "feast_contributions_matrix.tsv"), "\n", sep = "")
cat("  ", file.path(outdir, "feast_source_contributions.tsv"), "\n", sep = "")
cat("  ", file.path(outdir, "feast_source_contributions_summary.tsv"), "\n", sep = "")

# ---------------------------------------------------------------------------
# Stability check: re-run FEAST with two additional seeds
# ---------------------------------------------------------------------------
cat("\nRunning stability check (seeds 43, 44)...\n")
stability_results <- list()
stability_results[["seed_42"]] <- colMeans(feast_matrix_by_env)

for (extra_seed in c(43, 44)) {
  setwd(old_wd)
  set.seed(extra_seed)
  feast_extra <- FEAST::FEAST(
    C = as.matrix(counts),
    metadata = feast_meta,
    different_sources_flag = 0,
    dir_path = paste0(outdir, "/"),
    outfile = paste0("genus_stability_", extra_seed)
  )
  invisible(feast_extra)
  setwd(old_wd)

  extra_path <- file.path(outdir, paste0("genus_stability_", extra_seed, "_source_contributions_matrix.txt"))
  if (file.exists(extra_path)) {
    extra_mat <- read.table(extra_path, sep = "\t", header = TRUE, row.names = 1,
                            check.names = FALSE, comment.char = "", quote = "\"")
    extra_env <- sapply(unique(source_env_map), function(env) {
      rowSums(extra_mat[, source_env_map == env, drop = FALSE])
    })
    extra_env <- as.data.frame(extra_env, check.names = FALSE)
    stability_results[[paste0("seed_", extra_seed)]] <- colMeans(extra_env)
    file.remove(extra_path)
  }
}

stability_df <- do.call(rbind, lapply(names(stability_results), function(nm) {
  data.frame(run_seed = nm, t(stability_results[[nm]]), check.names = FALSE)
}))
write.table(stability_df, file = file.path(outdir, "feast_stability_check.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: feast_stability_check.tsv\n")
