# Run from the repository root. No package installation or original-input writes.
args <- commandArgs(trailingOnly = TRUE)
first <- if (length(args)) as.integer(args[1]) else 0L
stopifnot(is.finite(first), first >= 0L, first <= 5L)
required <- c("Seurat", "Matrix", "harmony", "uwot", "RANN", "digest", "scDblFinder", "SingleCellExperiment", "BiocParallel", "future")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install these dependencies first: ", paste(missing, collapse = ", "))
source("analysis_development_v2/scripts/config.R")
if (first <= 0 && any(!file.exists(c(NEW_INPUT, OLD_INPUT, PRIOR_INPUT)))) stop("Set BLA_NEW_RDS, BLA_OLD_RDS, BLA_PRIOR_RDS to your downloaded input files.")
if (!file.exists(file.path(ROOT, "resources", "ref_subclass_mean_log2_markers.csv"))) stop("Unzip the reference resource bundle into analysis_development_v2/resources/ first.")
rscript <- file.path(R.home("bin"), "Rscript.exe")
if (!file.exists(rscript)) rscript <- file.path(R.home("bin"), "Rscript")
steps <- c("00_inventory_qc.R", "01_doublets.R", "02_annotate.R", "03_embeddings.R", "04_expression.R")
for (i in seq_along(steps)) if (i - 1 >= first) {
  log <- file.path(ROOT, "logs", paste0(steps[i], ".log"))
  cat("Running", steps[i], "→", log, "\n")
  code <- system2(rscript, shQuote(file.path("analysis_development_v2", "scripts", steps[i])), stdout = log, stderr = log)
  if (code != 0) stop("Step failed: ", steps[i], ". Read ", log, "; no downstream step was run.")
}
python <- Sys.getenv("BLA_PYTHON", "python")
code <- system2(python, shQuote("analysis_development_v2/scripts/05_figures.py"),
  stdout = file.path(ROOT, "logs", "05_figures.log"), stderr = file.path(ROOT, "logs", "05_figures.log"))
if (code != 0) stop("Figure generation failed; inspect logs/05_figures.log.")
cat("Finished. See analysis_development_v2/results/. No KO/WT inference or cell-level significance tests were performed.\n")
