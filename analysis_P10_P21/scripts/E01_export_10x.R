# Export the annotated P10 + P21 object as 10x-style files (any tool can read them).
suppressPackageStartupMessages({ library(Seurat); library(Matrix) })
ROOT <- Sys.getenv("BLA_ROOT", "K:/scRNA_BLA_phd/P10_P21_amygdala_2026-10")
s <- readRDS(file.path(ROOT, "B06_P10_P21_annotated.rds"))
d <- file.path(ROOT, "data", "BLA_P10_P21_counts_10x"); dir.create(d, showWarnings = FALSE, recursive = TRUE)
m <- LayerData(s, "counts")
writeMM(m, file.path(d, "matrix.mtx"))
write.table(data.frame(rownames(m), rownames(m), "Gene Expression"), file.path(d, "features.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
writeLines(colnames(m), file.path(d, "barcodes.tsv"))
md <- s@meta.data
md$umap_1 <- Embeddings(s, "umap")[, 1]; md$umap_2 <- Embeddings(s, "umap")[, 2]
write.csv(cbind(barcode = colnames(s), md), file.path(d, "metadata.csv"), row.names = FALSE)
for (f in c("matrix.mtx", "features.tsv", "barcodes.tsv")) R.utils::gzip(file.path(d, f), overwrite = TRUE)
cat("exported", nrow(m), "genes x", ncol(m), "cells\n")
