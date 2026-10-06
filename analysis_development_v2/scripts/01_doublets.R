suppressPackageStartupMessages({library(Matrix); library(scDblFinder); library(SingleCellExperiment); library(BiocParallel)})
source("analysis_development_v2/scripts/config.R")
d <- readRDS(datapath("01_qc_input.rds"))
calls <- list()
for (id in unique(d$metadata$sample)) {
  out <- file.path(ROOT, "data", paste0("doublets_", id, ".rds"))
  if (file.exists(out)) { calls[[id]] <- readRDS(out); next }
  ix <- which(d$metadata$sample == id)
  logmsg("Doublet classification:", id, "barcodes:", length(ix))
  set.seed(SEED + match(id, sort(unique(d$metadata$sample))))
  sce <- SingleCellExperiment(list(counts = d$counts[, ix, drop = FALSE]))
  sce <- scDblFinder(sce, dbr = .08, dbr.sd = 1, nfeatures = 2000, dims = 20,
    BPPARAM = SerialParam(), verbose = TRUE)
  z <- data.frame(cell_id = colnames(sce), sample = id,
    doublet_score = sce$scDblFinder.score, doublet_call = as.character(sce$scDblFinder.class))
  stopifnot(identical(z$cell_id, colnames(d$counts)[ix]), all(is.finite(z$doublet_score)))
  saveRDS(z, out); calls[[id]] <- z
  logmsg(id, "predicted doublets:", sum(z$doublet_call == "doublet"))
  rm(sce); invisible(gc())
}
calls <- do.call(rbind, calls)
save_table(calls, "doublet_calls.csv")
sm <- do.call(rbind, lapply(split(calls, calls$sample), function(z)
  data.frame(sample = z$sample[1], qc_barcodes = nrow(z), predicted_doublets = sum(z$doublet_call == "doublet"),
    predicted_doublet_pct = 100 * mean(z$doublet_call == "doublet"), singlets = sum(z$doublet_call == "singlet"),
    dbr_prior = .08, dbr_sd = 1)))
save_table(sm, "doublet_summary.csv")
print(sm, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(ROOT, "logs", "sessionInfo_doublets.txt"))
