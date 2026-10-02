suppressPackageStartupMessages({ library(Seurat); library(dplyr) })
ROOT <- Sys.getenv("BLA_ROOT", "K:/scRNA_BLA_phd/P10_P21_amygdala_2026-10")
seu <- readRDS(file.path(ROOT, "B02_P10_P21_mapped.rds")); Idents(seu) <- "cl"
genes <- c("Snap25","Slc17a7","Slc17a6","Gad1","Gad2","Tbr1","Neurod6","Satb1","Rspo2","Ppp1r1b","Fezf2","Lypd1","Etv1","Cck","Lmo3","Otof",
 "Satb2","Nr4a2","Car3","Cux2","Lhx9","Nr2f2","Prkcd","Sst","Crh","Tac2","Nts","Penk","Pdyn","Tac1","Six3","Gpr88","Drd1","Drd2","Adora2a",
 "Foxp2","Tshz1","Meis2","Pbx3","Chst9","Lhx6","Lhx8","Sox6","Nfib","Pvalb","Vip","Npy","Lamp5","Sncg","Calb2","Reln","Nos1","Chodl","Htr3a",
 "Dcx","Sox11","Neurod1","Eomes","Aqp4","Gja1","Gfap","Olig2","Pdgfra","Bmp4","Enpp6","Mog","Plp1","Cx3cr1","P2ry12","Mrc1","Cldn5","Vtn","Pdgfrb","Dcn","Foxj1","Ttr","Mki67","Mef2c")
genes <- intersect(genes, rownames(seu))
X <- LayerData(seu, "data")[genes, ]
av <- sapply(levels(seu$cl), function(k) Matrix::rowMeans(expm1(X[, seu$cl == k])))
pc <- sapply(levels(seu$cl), function(k) Matrix::rowMeans(X[, seu$cl == k] > 0))
write.csv(round(pc * 100), file.path(ROOT, "tables/B03_cluster_pct_curated_markers.csv"))
mk <- FindAllMarkers(seu, only.pos = TRUE, max.cells.per.ident = 300, logfc.threshold = 0.5, min.pct = 0.25, verbose = FALSE)
mk <- mk %>% filter(p_val_adj < 0.01) %>% group_by(cluster) %>% arrange(desc(avg_log2FC), .by_group = TRUE)
write.csv(mk, file.path(ROOT, "tables/B03_cluster_markers_res1.csv"), row.names = FALSE)
top <- mk %>% mutate(score = avg_log2FC * (pct.1 - pct.2)) %>% slice_max(score, n = 12) %>% summarise(top = paste(gene, collapse = " "))
write.csv(top, file.path(ROOT, "tables/B03_cluster_top12.csv"), row.names = FALSE)
