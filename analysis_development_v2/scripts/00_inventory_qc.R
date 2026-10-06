suppressPackageStartupMessages({library(Seurat); library(Matrix); library(digest)})
source("analysis_development_v2/scripts/config.R")
set.seed(SEED)
logmsg("Reading new and original old inputs")
new <- readRDS(NEW_INPUT)
old <- readRDS(OLD_INPUT)
old <- subset(old, cells = colnames(old)[as.character(old$time) %in% c("P10", "P21")])
prior <- readRDS(PRIOR_INPUT)
n <- LayerData(new, assay = "RNA", layer = "counts")
o <- LayerData(old, assay = "RNA", layer = "counts")
shared <- intersect(rownames(n), rownames(o))
save_table(data.frame(gene = rownames(n), present_in_old = rownames(n) %in% shared), "gene_coverage.csv")
infiles <- c(new = NEW_INPUT, old = OLD_INPUT, prior = PRIOR_INPUT)
save_table(data.frame(input = names(infiles), path = unname(infiles),
  bytes = file.info(infiles)$size, sha256 = vapply(infiles, digest, "", file = TRUE, algo = "sha256")), "input_manifest.csv")

checks <- lapply(list(new = n, old = o), function(m) c(finite = all(is.finite(m@x)), nonnegative = all(m@x >= 0),
  integer = all(m@x == floor(m@x)), unique_genes = !anyDuplicated(rownames(m)), unique_cells = !anyDuplicated(colnames(m))))
stopifnot(all(unlist(checks)), identical(colnames(n), rownames(new[[]])), identical(colnames(o), rownames(old[[]])))
save_table(data.frame(input = rep(names(checks), each = 5), check = rep(names(checks[[1]]), 2),
  passed = unlist(checks, use.names = FALSE)), "input_integrity.csv")

nm <- data.frame(sample = as.character(new$sample), age = as.character(new$age),
  run = as.character(new$run), index = as.character(new$index), source = "new_2026",
  original_cell_id = colnames(n), nCount_full = Matrix::colSums(n), nFeature_full = Matrix::colSums(n > 0),
  prior_cell_type = NA_character_, row.names = paste0("new_", colnames(n)))
om <- data.frame(sample = as.character(old$orig.ident), age = as.character(old$time),
  run = NA_character_, index = NA_character_, source = "old_2023",
  original_cell_id = colnames(o), nCount_full = Matrix::colSums(o), nFeature_full = Matrix::colSums(o > 0),
  prior_cell_type = as.character(prior$cell_type[match(colnames(o), colnames(prior))]), row.names = paste0("old_", colnames(o)))
stopifnot(all(nm$nCount_full == new$nCount_RNA), all(nm$nFeature_full == new$nFeature_RNA),
  all(om$nCount_full == old$nCount_RNA), all(om$nFeature_full == old$nFeature_RNA))
logmsg("Shared genes:", length(shared))
n <- n[shared, ]; o <- o[shared, ]

# Comparing same cell-barcode suffixes tests resequencing, rather than merely
# comparing different sample prefixes. Random same-library pairs are a baseline.
strip_id <- function(x) sub("^.*_", "", x)
dup <- list()
for (ns in unique(nm$sample)) for (os in unique(om$sample)) {
  ni <- which(nm$sample == ns); oi <- which(om$sample == os)
  nb <- strip_id(colnames(n)[ni]); ob <- strip_id(colnames(o)[oi]); b <- intersect(nb, ob)
  if (!length(b)) next
  a <- n[, ni[match(b, nb)], drop = FALSE]; z <- o[, oi[match(b, ob)], drop = FALSE]
  cs <- Matrix::colSums(a * z) / sqrt(Matrix::colSums(a * a) * Matrix::colSums(z * z))
  perm <- sample(seq_len(ncol(z)))
  base <- Matrix::colSums(a * z[, perm, drop = FALSE]) / sqrt(Matrix::colSums(a * a) * Matrix::colSums(z[, perm, drop = FALSE]^2))
  dup[[length(dup) + 1L]] <- data.frame(new_sample = ns, old_sample = os, matched_barcodes = length(b),
    median_cosine = median(cs), q95_cosine = unname(quantile(cs, .95)), n_cosine_over097 = sum(cs > .97),
    random_pair_median = median(base), suspected_resequencing = sum(cs > .97) >= 20)
}
dup <- do.call(rbind, dup)
save_table(dup, "same_barcode_resequencing_check.csv")
if (any(dup$suspected_resequencing)) stop("Possible resequenced library: inspect same_barcode_resequencing_check.csv before pooling.")
colnames(n) <- rownames(nm); colnames(o) <- rownames(om)
counts <- cbind(n, o); md <- rbind(nm, om)
md$nCount_shared <- Matrix::colSums(counts)
md$nFeature_shared <- Matrix::colSums(counts > 0)
mt <- grep("^mt-", shared)
md$percent_mt <- 100 * Matrix::colSums(counts[mt, , drop = FALSE]) / md$nCount_shared
md$qc_primary <- md$nFeature_shared >= 500 & md$nCount_shared >= 1000 & md$percent_mt <= 20
md$qc_loose <- md$nFeature_shared >= 300 & md$nCount_shared >= 500 & md$percent_mt <= 25
md$qc_strict_mt <- md$qc_primary & md$percent_mt <= 15
summary <- do.call(rbind, lapply(split(seq_len(nrow(md)), md$sample), function(i) {
  v <- md[i, ]
  data.frame(sample = v$sample[1], age = v$age[1], source = v$source[1], run = v$run[1], index = v$index[1],
    input_barcodes = nrow(v), qc_loose = sum(v$qc_loose), qc_primary = sum(v$qc_primary), strict_mt15 = sum(v$qc_strict_mt),
    median_umi_before = median(v$nCount_shared), median_genes_before = median(v$nFeature_shared),
    median_mt_before = median(v$percent_mt), median_umi_after = median(v$nCount_shared[v$qc_primary]),
    median_genes_after = median(v$nFeature_shared[v$qc_primary]), median_mt_after = median(v$percent_mt[v$qc_primary]))
}))
save_table(summary, "sample_qc_summary.csv")
save_table(cbind(cell_id = rownames(md), md), "all_input_barcodes_qc.csv")
saveRDS(list(counts = counts[, md$qc_primary, drop = FALSE], metadata = md[md$qc_primary, ], shared_genes = shared),
  datapath("01_qc_input.rds"))
writeLines(capture.output(sessionInfo()), file.path(ROOT, "logs", "sessionInfo_inventory.txt"))
print(summary, row.names = FALSE)
logmsg("Inventory and QC checkpoint complete")
