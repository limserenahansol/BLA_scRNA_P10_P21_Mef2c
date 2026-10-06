suppressPackageStartupMessages({library(Matrix); library(scDblFinder); library(SingleCellExperiment); library(BiocParallel)})
source("analysis_development_v2/scripts/config.R")
d <- readRDS(datapath("01_qc_input.rds"))
id <- "P10_2"; ix <- which(d$metadata$sample == id)
set.seed(SEED + match(id, sort(unique(d$metadata$sample))))
w <- character()
sce <- withCallingHandlers(scDblFinder(SingleCellExperiment(list(counts = d$counts[, ix, drop = FALSE])),
  dbr = .08, dbr.sd = 1, nfeatures = 2000, dims = 20, BPPARAM = SerialParam(), verbose = TRUE),
  warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
prior <- readRDS(datapath(paste0("doublets_", id, ".rds")))
save_table(data.frame(warning = names(table(w)), occurrences = as.integer(table(w))), "doublet_warning_audit.csv")
save_table(data.frame(sample = id, identical_calls = identical(as.character(sce$scDblFinder.class), prior$doublet_call),
  max_score_difference = max(abs(sce$scDblFinder.score - prior$doublet_score)), warnings_captured = length(w)), "doublet_repeat_validation.csv")
stopifnot(identical(as.character(sce$scDblFinder.class), prior$doublet_call))
