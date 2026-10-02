suppressPackageStartupMessages({ library(Seurat); library(harmony); library(dplyr) })
set.seed(1)
ROOT <- Sys.getenv("BLA_ROOT", "K:/scRNA_BLA_phd/P10_P21_amygdala_2026-10")
seu <- readRDS(file.path(ROOT, "B02_P10_P21_mapped.rds"))
glut_cl <- c(5, 6, 7, 9, 11, 12, 17, 23, 33)
g <- subset(seu, subset = cl %in% glut_cl)
g[["RNA"]] <- split(g[["RNA"]], f = g$sample)
g <- NormalizeData(g, verbose = FALSE) |> FindVariableFeatures(nfeatures = 2000, verbose = FALSE) |>
  ScaleData(verbose = FALSE) |> RunPCA(npcs = 30, verbose = FALSE)
g <- IntegrateLayers(g, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", verbose = FALSE)
g <- FindNeighbors(g, reduction = "harmony", dims = 1:20, verbose = FALSE)
for (r in c(0.4, 0.8)) g <- FindClusters(g, resolution = r, verbose = FALSE)
g <- RunUMAP(g, reduction = "harmony", dims = 1:20, verbose = FALSE)
g <- JoinLayers(g); g$gcl <- g$RNA_snn_res.0.8; Idents(g) <- "gcl"
saveRDS(g, file.path(ROOT, "B04_glut.rds"))
genes <- c("Slc17a7","Slc17a6","Tbr1","Neurod6","Rspo2","Ppp1r1b","Fezf2","Lypd1","Etv1","Cck","Lmo3","Otof","Tshz2","Satb1","Satb2","Nr4a2","Bcl11b",
  "Lhx9","Lhx2","Nr2f2","Sema5a","Sema3e","Cdh9","Cdh8","Slit3","C1ql3","Trpc5","Npnt","Grp","Car3","Tfap2d","Sox11","Dcx","Penk","Htr2c","Nos1","Zfp804b","Mef2c","Gad2","Meis2","Foxp2")
genes <- intersect(genes, rownames(g)); X <- LayerData(g, "data")[genes, ]
pc <- sapply(levels(g$gcl), function(k) round(100 * Matrix::rowMeans(X[, g$gcl == k] > 0)))
mk <- FindAllMarkers(g, only.pos = TRUE, max.cells.per.ident = 300, logfc.threshold = 0.4, min.pct = 0.2, verbose = FALSE) |> filter(p_val_adj < 0.01)
write.csv(mk, file.path(ROOT, "tables/B04_glut_subcluster_markers.csv"), row.names = FALSE)
top <- mk |> mutate(s = avg_log2FC * (pct.1 - pct.2)) |> group_by(cluster) |> slice_max(s, n = 12) |> summarise(top = paste(gene, collapse = " "))
info <- g@meta.data |> group_by(gcl) |> summarise(n = n(), P10 = sum(age == "P10"), P21 = sum(age == "P21"), nF = median(nFeature_RNA),
  from_cl = paste(names(sort(table(cl), decreasing = TRUE))[1:2], collapse = "/"),
  allen_sup = names(sort(table(allen_supertype), decreasing = TRUE))[1],
  allen_sup_frac = round(max(table(allen_supertype)) / n(), 2))
write.csv(cbind(info, t(pc)), file.path(ROOT, "tables/B04_glut_subcluster_summary.csv"), row.names = FALSE)
print(as.data.frame(info)); print(t(pc)); print(as.data.frame(top))
