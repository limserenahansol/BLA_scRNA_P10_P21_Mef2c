suppressPackageStartupMessages({ library(Seurat); library(harmony); library(dplyr) })
set.seed(1)
ROOT <- Sys.getenv("BLA_ROOT", "K:/scRNA_BLA_phd/P10_P21_amygdala_2026-10")
seu <- readRDS(file.path(ROOT, "B02_P10_P21_mapped.rds"))
glut_cl <- c(1, 2, 8, 14, 15, 16, 22, 26)
g <- subset(seu, subset = cl %in% glut_cl)
g[["RNA"]] <- split(g[["RNA"]], f = g$sample)
g <- NormalizeData(g, verbose = FALSE) |> FindVariableFeatures(nfeatures = 2000, verbose = FALSE) |>
  ScaleData(verbose = FALSE) |> RunPCA(npcs = 30, verbose = FALSE)
g <- IntegrateLayers(g, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", verbose = FALSE)
g <- FindNeighbors(g, reduction = "harmony", dims = 1:20, verbose = FALSE)
for (r in c(0.4, 0.8)) g <- FindClusters(g, resolution = r, verbose = FALSE)
g <- RunUMAP(g, reduction = "harmony", dims = 1:20, verbose = FALSE)
g <- JoinLayers(g); g$gcl <- g$RNA_snn_res.0.8; Idents(g) <- "gcl"
saveRDS(g, file.path(ROOT, "B05_gaba.rds"))
genes <- c("Gad2","Slc17a7","Slc17a6","Foxp2","Tshz1","Meis2","Pbx3","Chst9","Six3","Sp9","Cyp26b1","Prkcd","Sst","Crh","Tac2","Nts","Penk","Pdyn","Tac1","Drd1","Drd2","Adora2a","Gpr88","Isl1","Ebf1","Lhx6","Lhx8","Nr2e1","Sox6","Nfib","Pvalb","Vip","Cck","Npy","Lamp5","Sncg","Calb2","Reln","Nos1","Chodl","Htr3a","Adarb2","Prox1","Kcnc2","Esr1","Otp","Mef2c","Pdgfra","Inpp5d","Slc1a3")
genes <- intersect(genes, rownames(g)); X <- LayerData(g, "data")[genes, ]
pc <- sapply(levels(g$gcl), function(k) round(100 * Matrix::rowMeans(X[, g$gcl == k] > 0)))
mk <- FindAllMarkers(g, only.pos = TRUE, max.cells.per.ident = 300, logfc.threshold = 0.4, min.pct = 0.2, verbose = FALSE) |> filter(p_val_adj < 0.01)
write.csv(mk, file.path(ROOT, "tables/B05_gaba_subcluster_markers.csv"), row.names = FALSE)
top <- mk |> mutate(s = avg_log2FC * (pct.1 - pct.2)) |> group_by(cluster) |> slice_max(s, n = 12) |> summarise(top = paste(gene, collapse = " "))
info <- g@meta.data |> group_by(gcl) |> summarise(n = n(), P10 = sum(age == "P10"), P21 = sum(age == "P21"), nF = median(nFeature_RNA),
  from_cl = paste(names(sort(table(cl), decreasing = TRUE))[1:2], collapse = "/"),
  allen_sup = names(sort(table(allen_supertype), decreasing = TRUE))[1],
  allen_sup_frac = round(max(table(allen_supertype)) / n(), 2))
write.csv(cbind(info, t(pc)), file.path(ROOT, "tables/B05_gaba_subcluster_summary.csv"), row.names = FALSE)
print(as.data.frame(info)); print(t(pc)); print(as.data.frame(top))
