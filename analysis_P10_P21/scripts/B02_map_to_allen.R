# Correlation mapping of P10/P21 BLA cells to the Allen WMB amygdala reference
# (CTXsp + sAMY ROIs, adult): subclass, then supertype within the subclass.
suppressPackageStartupMessages({ library(Seurat); library(dplyr); library(Matrix) })
ROOT <- Sys.getenv("BLA_ROOT", "K:/scRNA_BLA_phd/P10_P21_amygdala_2026-10")
seu <- readRDS(file.path(ROOT, "B01_P10_P21_processed.rds"))
X <- LayerData(seu, "data")

corr_map <- function(ref, X) {
  g <- intersect(rownames(ref), rownames(X))
  R <- as.matrix(ref[g, ]); R <- scale(R)                     # gene-wise z over columns? -> no: centre per group
  R <- apply(as.matrix(ref[g, ]), 2, function(v) (v - mean(v)) / sd(v))
  out <- matrix(0, ncol(X), ncol(R), dimnames = list(colnames(X), colnames(R)))
  for (s in seq(1, ncol(X), 4000)) {
    idx <- s:min(s + 3999, ncol(X))
    C <- as.matrix(X[g, idx]); C <- apply(C, 2, function(v) (v - mean(v)) / (sd(v) + 1e-9))
    out[idx, ] <- crossprod(C, R) / (length(g) - 1)
  }
  list(r = out, n_genes = length(g))
}
ref_sub <- read.csv(file.path(ROOT, "ref/ref_subclass_mean_log2_markers.csv"), row.names = 1, check.names = FALSE)
ref_sup <- read.csv(file.path(ROOT, "ref/ref_supertype_mean_log2_markers.csv"), row.names = 1, check.names = FALSE)
st2sc <- read.csv(file.path(ROOT, "ref/supertype_to_subclass.csv"), check.names = FALSE)

m1 <- corr_map(ref_sub, X); cat("subclass mapping on", m1$n_genes, "genes\n")
top2 <- t(apply(m1$r, 1, function(v) { o <- order(v, decreasing = TRUE)[1:2]; c(o, v[o]) }))
seu$allen_subclass <- colnames(m1$r)[top2[, 1]]
seu$allen_subclass_r <- top2[, 3]
seu$allen_subclass_margin <- top2[, 3] - top2[, 4]
seu$allen_class <- st2sc$class[match(seu$allen_subclass, st2sc$subclass)]
seu$allen_nt <- st2sc$neurotransmitter[match(seu$allen_subclass, st2sc$subclass)]

m2 <- corr_map(ref_sup, X); cat("supertype mapping on", m2$n_genes, "genes\n")
sup_sc <- st2sc$subclass[match(colnames(m2$r), st2sc$supertype)]
seu$allen_supertype <- vapply(seq_len(ncol(seu)), function(i) {
  ok <- which(sup_sc == seu$allen_subclass[i]); if (!length(ok)) return(NA_character_)
  colnames(m2$r)[ok[which.max(m2$r[i, ok])]] }, "")

# Cluster-level: correlate cluster mean expression with the reference
cl_mean <- sapply(split(seq_len(ncol(seu)), seu$cl), function(ix) Matrix::rowMeans(X[, ix, drop = FALSE]))
mc <- corr_map(ref_sub, Matrix(cl_mean, sparse = TRUE))$r
cl_tab <- data.frame(cl = rownames(mc),
  n = as.vector(table(seu$cl)[rownames(mc)]),
  pct_P10 = round(100 * tapply(seu$age == "P10", seu$cl, mean)[rownames(mc)], 1),
  top1 = colnames(mc)[apply(mc, 1, which.max)],
  r1 = round(apply(mc, 1, max), 3),
  top2 = colnames(mc)[apply(mc, 1, function(v) order(v, decreasing = TRUE)[2])],
  r2 = round(apply(mc, 1, function(v) sort(v, decreasing = TRUE)[2]), 3))
vote <- seu@meta.data %>% count(cl, allen_subclass) %>% group_by(cl) %>%
  mutate(frac = n / sum(n)) %>% slice_max(frac, n = 1, with_ties = FALSE)
cl_tab$cell_vote <- vote$allen_subclass[match(cl_tab$cl, vote$cl)]
cl_tab$cell_vote_frac <- round(vote$frac[match(cl_tab$cl, vote$cl)], 2)
cl_tab <- cl_tab[order(as.integer(cl_tab$cl)), ]
write.csv(cl_tab, file.path(ROOT, "tables/B02_cluster_to_allen_subclass.csv"), row.names = FALSE)
print(cl_tab, row.names = FALSE)
saveRDS(list(cell_sub = m1$r, cl_sub = mc), file.path(ROOT, "B02_allen_corr.rds"))
saveRDS(seu, file.path(ROOT, "B02_P10_P21_mapped.rds"))
