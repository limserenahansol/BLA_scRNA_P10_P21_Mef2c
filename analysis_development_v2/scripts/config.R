options(stringsAsFactors = FALSE, width = 160)
Sys.setenv(OMP_NUM_THREADS = "4", OPENBLAS_NUM_THREADS = "4")
if (requireNamespace("future", quietly = TRUE)) future::plan("sequential")
ROOT <- normalizePath(Sys.getenv("BLA_V2_ROOT", "analysis_development_v2"), mustWork = TRUE)
NEW_INPUT <- Sys.getenv("BLA_NEW_RDS", "K:/scRNA_BLA_phd/hansol2_combined_seurat.rds")
OLD_INPUT <- Sys.getenv("BLA_OLD_RDS", "K:/scRNA_BLA_phd/hansol/seu_harm_qc.rds")
PRIOR_INPUT <- Sys.getenv("BLA_PRIOR_RDS", "K:/scRNA_BLA_phd/P10_P21_amygdala_2026-10/B06_P10_P21_annotated.rds")
SEED <- 1062026L
for (d in c("data", "results", "results/tables", "results/figures", "logs", "resources", "release"))
  dir.create(file.path(ROOT, d), recursive = TRUE, showWarnings = FALSE)
tabpath <- function(x) file.path(ROOT, "results", "tables", x)
datapath <- function(x) file.path(ROOT, "data", x)
logmsg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), paste(..., collapse = " "), "\n")
save_table <- function(d, name) write.csv(d, tabpath(name), row.names = FALSE, na = "")
candidate_genes <- c("Mef2c", "Mef2a", "Mef2d", "Arc", "Homer1", "Nr4a1", "Bdnf", "Pcdh10", "Npas4", "Fos", "Egr1",
  "Ncam1", "Sema6d", "Cadm2", "Pcdh7", "Grip1", "C1qa", "C1qb", "C1qc", "C3", "C4b", "Itgam", "Trem2", "Tyrobp",
  "Cd68", "Cx3cr1", "P2ry12", "Csf1r", "Mertk", "Megf10", "Axl", "Gas6", "Sema3e", "Plxnd1", "Nrp1", "Nrp2",
  "Kirrel3", "Sdk2", "Sema5b", "Cdh8", "Cdh9", "Cdh13", "Cdh18", "C1ql3", "Slit1", "Slit2", "Slit3", "Robo1", "Robo2",
  "Ntng1", "Ntng2", "Lrrc4c", "Lrrc4", "Adgrb3", "Igsf21", "Nrxn2", "Cntn4", "Cntn6", "Ptprg", "Plxna1", "Plxna3")
age_levels <- c("E18", "P0", "P10", "P21")
