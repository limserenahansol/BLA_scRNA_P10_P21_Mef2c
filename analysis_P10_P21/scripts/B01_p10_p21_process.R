# P10 + P21 BLA cells from the existing QC'd object (Tobias, 2023: mito <20%,
# hb <20%, DoubletFinder singlets) -> Seurat v5, Harmony on sample, finer clusters.
suppressPackageStartupMessages({ library(Seurat); library(harmony); library(dplyr) })
set.seed(1)
ROOT <- Sys.getenv("BLA_ROOT", "K:/scRNA_BLA_phd/P10_P21_amygdala_2026-10")
old <- readRDS(Sys.getenv("BLA_INPUT_RDS", file.path(ROOT, "data", "seu_harm_qc.rds")))
old <- subset(old, subset = time %in% c("P10", "P21"))

counts <- GetAssayData(old, assay = "RNA", layer = "counts")
md <- old@meta.data[, c("orig.ident", "time", "nCount_RNA", "nFeature_RNA",
                        "percent_mito", "percent_ribo", "S.Score", "G2M.Score", "Phase")]
md$old_cluster_res0.4 <- old$RNA_snn_res.0.4
names(md)[names(md) == "orig.ident"] <- "sample"
names(md)[names(md) == "time"] <- "age"
rm(old); invisible(gc())

seu <- CreateSeuratObject(counts, meta.data = md)
seu$sample <- factor(seu$sample, levels = c("P10s1", "P21s1", "P21s2"))
seu$age <- factor(seu$age, levels = c("P10", "P21"))
cat("cells per sample:\n"); print(table(seu$sample))

seu[["RNA"]] <- split(seu[["RNA"]], f = seu$sample)
seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, nfeatures = 3000, verbose = FALSE)
seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, npcs = 50, verbose = FALSE)
seu <- IntegrateLayers(seu, method = HarmonyIntegration, orig.reduction = "pca",
                       new.reduction = "harmony", verbose = FALSE)
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:30, verbose = FALSE)
for (r in c(0.5, 1.0, 2.0)) seu <- FindClusters(seu, resolution = r, verbose = FALSE)
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:30, verbose = FALSE)
seu <- JoinLayers(seu)
seu$cl <- seu$RNA_snn_res.1
cat("clusters res1:", nlevels(seu$cl), "\n")
saveRDS(seu, file.path(ROOT, "B01_P10_P21_processed.rds"))
cat("percent_mito summary:\n"); print(summary(seu$percent_mito))
