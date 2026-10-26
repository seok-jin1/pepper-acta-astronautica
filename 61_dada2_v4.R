# 61_dada2_v4.R — denoise each study separately on the shared 240-nt V4 window, then merge
#
# Each study gets its own error model and chimera removal (different runs / platforms);
# the per-study tables are merged afterwards, so identical sequences become the same ASV.
#
# Inputs : $V4_WORK/input/<study>/*.fastq.gz   (from 60_v4_prepare.sh)
# Outputs: $V4_WORK/dada2/seqtab_<study>.rds, seqtab_merged.rds
#          $V4_WORK/dada2/feature-table.tsv   (ASV id = md5 of sequence; rows = ASVs, columns = samples)
#          $V4_WORK/dada2/asv.fasta, $V4_WORK/dada2/track.tsv
suppressMessages({ library(dada2); library(digest) })

work <- Sys.getenv("V4_WORK", file.path(Sys.getenv("HOME"), "pepper-v4", "work"))
threads <- as.integer(Sys.getenv("V4_THREADS", "6"))
trunc_len <- 240
studies <- c("OSD-772", "PRJNA1145089", "OSD-267", "OSD-269")
out_dir <- file.path(work, "dada2")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(42)

tabs <- list(); tracks <- list()
for (study in studies) {
  fq <- sort(list.files(file.path(work, "input", study), pattern = "fastq.gz$", full.names = TRUE))
  if (length(fq) == 0) { message("skip ", study, " (no input)"); next }
  ids <- sub("\\.fastq\\.gz$", "", basename(fq))
  filt <- file.path(work, "filtered", study, basename(fq))
  dir.create(dirname(filt[1]), showWarnings = FALSE, recursive = TRUE)
  ft <- filterAndTrim(fq, filt, truncLen = trunc_len, maxN = 0, maxEE = 2, truncQ = 2,
                      rm.phix = TRUE, compress = TRUE, multithread = threads)
  keep <- ft[, "reads.out"] > 0
  err <- learnErrors(filt[keep], nbases = 1e8, multithread = threads, randomize = TRUE)
  dd <- dada(filt[keep], err = err, pool = "pseudo", multithread = threads)
  st <- makeSequenceTable(dd)
  rownames(st) <- ids[keep]
  st_nochim <- removeBimeraDenovo(st, method = "consensus", multithread = threads)
  saveRDS(st_nochim, file.path(out_dir, paste0("seqtab_", study, ".rds")))
  getN <- function(x) sum(getUniques(x))
  denoised <- if (sum(keep) == 1) getN(dd) else sapply(dd, getN)
  tr <- data.frame(study = study, sample = ids, input = ft[, "reads.in"], filtered = ft[, "reads.out"],
                   denoised = 0, nonchim = 0)
  tr$denoised[keep] <- denoised
  tr$nonchim[keep] <- rowSums(st_nochim)
  tabs[[study]] <- st_nochim; tracks[[study]] <- tr
  message(sprintf("%s: %d samples, %d ASVs, %.1f%% reads kept", study, sum(keep), ncol(st_nochim),
                  100 * sum(tr$nonchim) / sum(tr$input)))
}

merged <- mergeSequenceTables(tables = tabs, repeats = "error")
saveRDS(merged, file.path(out_dir, "seqtab_merged.rds"))
seqs <- colnames(merged)
asv_id <- vapply(seqs, digest, character(1), algo = "md5", serialize = FALSE)
writeLines(paste0(">", asv_id, "\n", seqs), file.path(out_dir, "asv.fasta"))
ft_out <- data.frame(`#OTU ID` = asv_id, t(merged), check.names = FALSE)
write.table(ft_out, file.path(out_dir, "feature-table.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, tracks), file.path(out_dir, "track.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("merged: %d samples x %d ASVs", nrow(merged), ncol(merged)))
