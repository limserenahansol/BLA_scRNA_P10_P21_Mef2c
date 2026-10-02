#!/usr/bin/env Rscript
# =============================================================================
# E18 / P0 scRNA-seq - first-pass Seurat preprocessing
#
# Input : Cell Ranger outs per sample (filtered_feature_bc_matrix.h5)
# Output: QC'd, doublet-filtered, normalised, clustered Seurat object, UMAPs,
#         neuron vs non-neuron call, broad cell class, marker tables, QC tables.
# This is a hand-off object for the collaborator, not a final analysis:
# thresholds are deliberately permissive and every removed cell is logged.
#
# Run (Windows):
#   Rscript 05_seurat_preprocess.R
#   Rscript 05_seurat_preprocess.R --cr_dir=K:/.../results/cellranger --out_dir=K:/.../results/seurat
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(Matrix)
  library(dplyr); library(ggplot2); library(patchwork)
  library(SingleCellExperiment); library(scDblFinder); library(harmony)
})
options(future.globals.maxSize = 8 * 1024^3)
set.seed(1)

# ---- Paths & parameters ------------------------------------------------------
get_script_dir <- function() {
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a))) else getwd()
}
PROJ_DIR <- dirname(get_script_dir())
args <- commandArgs(TRUE)
arg <- function(key, default) {
  hit <- grep(paste0("^--", key, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", key, "="), "", hit) else default
}
CR_DIR      <- arg("cr_dir",  file.path(PROJ_DIR, "results", "cellranger"))
OUT_DIR     <- arg("out_dir", file.path(PROJ_DIR, "results", "seurat"))
SAMPLES_CSV <- arg("samples", file.path(PROJ_DIR, "samples.csv"))

QC <- list(
  min_features = as.numeric(arg("min_features", 500)),   # genes per cell
  min_counts   = as.numeric(arg("min_counts",   1000)),  # UMIs per cell
  max_mt       = as.numeric(arg("max_mt",       10)),    # % mitochondrial (use 5 for nuclei)
  max_hb       = 5,                                       # % haemoglobin
  mad_upper    = 5                                        # drop cells > median + 5 MAD log10(UMI) (per sample)
)
N_HVG <- 3000; N_PCS <- 30; RESOLUTION <- 0.5
REMOVE_DOUBLETS <- TRUE

dir.create(file.path(OUT_DIR, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "tables"),  recursive = TRUE, showWarnings = FALSE)
fig <- function(name) file.path(OUT_DIR, "figures", name)
tab <- function(name) file.path(OUT_DIR, "tables", name)
save_plot <- function(p, name, w, h) {
  ggsave(fig(paste0(name, ".png")), p, width = w, height = h, dpi = 200, bg = "white")
  ggsave(fig(paste0(name, ".pdf")), p, width = w, height = h)
}
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

# ---- Sample sheet -------------------------------------------------------------
sheet <- read.csv(SAMPLES_CSV, stringsAsFactors = FALSE, strip.white = TRUE) %>%
  distinct(sample_id, stage)
msg("samples: ", paste(sheet$sample_id, collapse = ", "))

read_sample <- function(sid) {
  h5  <- file.path(CR_DIR, sid, "filtered_feature_bc_matrix.h5")
  if (!file.exists(h5)) h5 <- file.path(CR_DIR, sid, "outs", "filtered_feature_bc_matrix.h5")
  mtx <- sub("\\.h5$", "", h5)
  m <- if (file.exists(h5)) Read10X_h5(h5) else if (dir.exists(mtx)) Read10X(mtx) else NULL
  if (is.list(m)) m <- m[["Gene Expression"]]
  m
}

# ---- 1. Load + per-cell QC metrics -------------------------------------------
objs <- list()
for (i in seq_len(nrow(sheet))) {
  sid <- sheet$sample_id[i]
  m <- read_sample(sid)
  if (is.null(m)) { msg("WARNING: no Cell Ranger matrix for ", sid, " - skipped"); next }
  so <- CreateSeuratObject(m, project = sid, min.cells = 3, min.features = 200)
  so$sample <- sid; so$stage <- sheet$stage[i]
  so[["percent_mt"]]   <- PercentageFeatureSet(so, pattern = "^mt-")
  so[["percent_ribo"]] <- PercentageFeatureSet(so, pattern = "^Rp[sl]")
  so[["percent_hb"]]   <- PercentageFeatureSet(so, pattern = "^Hb[ab]-")
  objs[[sid]] <- so
  msg(sid, ": ", ncol(so), " cells loaded")
}
stopifnot(length(objs) > 0)

qc_feats <- c("nFeature_RNA", "nCount_RNA", "percent_mt", "percent_hb")
qc_violin <- function(o, title) {
  (VlnPlot(o, features = qc_feats, group.by = "sample", pt.size = 0, ncol = 4) &
     theme(axis.title.x = element_blank())) + plot_annotation(title = title)
}

# ---- 2. Cell filtering + doublets, per sample ---------------------------------
qc_log <- list()
for (sid in names(objs)) {
  so <- objs[[sid]]; md <- so@meta.data
  lc <- log10(md$nCount_RNA)
  upper <- median(lc) + QC$mad_upper * mad(lc)
  fail <- data.frame(
    low_genes = md$nFeature_RNA < QC$min_features,
    low_umi   = md$nCount_RNA   < QC$min_counts,
    high_umi  = lc > upper,
    high_mt   = md$percent_mt   > QC$max_mt,
    high_hb   = md$percent_hb   > QC$max_hb)
  so$qc_pass <- rowSums(fail) == 0

  # scDblFinder on cells that pass QC (doublet rate set from cell number)
  so_q <- subset(so, cells = colnames(so)[so$qc_pass])
  sce <- scDblFinder(SingleCellExperiment(list(counts = LayerData(so_q, "counts"))),
                     verbose = FALSE)
  so$doublet_class <- NA_character_
  so$doublet_score <- NA_real_
  so$doublet_class[colnames(so_q)] <- as.character(sce$scDblFinder.class)
  so$doublet_score[colnames(so_q)] <- sce$scDblFinder.score
  keep <- so$qc_pass & (!REMOVE_DOUBLETS | so$doublet_class %in% "singlet")

  qc_log[[sid]] <- data.frame(
    sample = sid, stage = so$stage[1], cells_cellranger = ncol(so),
    fail_low_genes = sum(fail$low_genes), fail_low_umi = sum(fail$low_umi),
    fail_high_umi = sum(fail$high_umi), fail_high_mt = sum(fail$high_mt),
    fail_high_hb = sum(fail$high_hb), pass_qc = sum(so$qc_pass),
    doublets = sum(so$doublet_class %in% "doublet", na.rm = TRUE),
    cells_kept = sum(keep), pct_kept = round(100 * mean(keep), 1),
    median_umi_kept = median(so$nCount_RNA[keep]),
    median_genes_kept = median(so$nFeature_RNA[keep]),
    median_pct_mt_kept = round(median(so$percent_mt[keep]), 2),
    high_umi_cutoff = round(10^upper))
  objs[[sid]] <- so
  objs[[paste0(sid, "__kept")]] <- subset(so, cells = colnames(so)[keep])
  msg(sid, ": ", sum(keep), "/", ncol(so), " cells kept (",
      qc_log[[sid]]$doublets, " doublets)")
}
qc_tab <- bind_rows(qc_log); rownames(qc_tab) <- NULL
write.csv(qc_tab, tab("qc_cell_filtering_per_sample.csv"), row.names = FALSE)
print(qc_tab[, c("sample", "cells_cellranger", "pass_qc", "doublets", "cells_kept",
                 "median_umi_kept", "median_genes_kept")])

raw_ids  <- intersect(sheet$sample_id, names(objs))
kept_ids <- paste0(raw_ids, "__kept")
raw_all  <- merge(objs[[raw_ids[1]]], objs[raw_ids[-1]], add.cell.ids = raw_ids)
raw_all$sample <- factor(raw_all$sample, levels = raw_ids)
save_plot(qc_violin(raw_all, "QC before filtering (Cell Ranger cells)"), "01_qc_violin_before", 14, 4.5)
p_sc <- FeatureScatter(raw_all, "nCount_RNA", "nFeature_RNA", group.by = "qc_pass",
                       split.by = "sample", pt.size = 0.2, shuffle = TRUE) +
  scale_x_log10() + scale_y_log10()
save_plot(p_sc, "01_qc_scatter_umi_genes", 4 * length(kept_ids), 4)
rm(raw_all)

# ---- 3. Merge + normalise + PCA ------------------------------------------------
seu <- merge(objs[[kept_ids[1]]], objs[kept_ids[-1]], add.cell.ids = sub("__kept$", "", kept_ids))
rm(objs); invisible(gc())
seu$sample <- factor(seu$sample, levels = sheet$sample_id)
stage_order <- c("E18", "P0", "P10", "P21", "P56")
seu$stage  <- factor(seu$stage, levels = c(intersect(stage_order, unique(seu$stage)),
                                            setdiff(unique(seu$stage), stage_order)))
save_plot(qc_violin(seu, "QC after filtering"), "02_qc_violin_after", 14, 4.5)

msg("normalise / HVG / PCA on ", ncol(seu), " cells")
seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, nfeatures = N_HVG, verbose = FALSE)
seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, npcs = 50, verbose = FALSE)
save_plot(ElbowPlot(seu, ndims = 50), "03_pca_elbow", 6, 4)

# Unintegrated UMAP: shows raw sample/stage structure (batch check)
seu <- RunUMAP(seu, reduction = "pca", dims = 1:N_PCS, reduction.name = "umap.unintegrated",
               verbose = FALSE)

# ---- 4. Harmony across samples, cluster ----------------------------------------
msg("harmony integration across samples")
seu <- IntegrateLayers(seu, method = HarmonyIntegration, orig.reduction = "pca",
                       new.reduction = "harmony", verbose = FALSE)
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:N_PCS, verbose = FALSE)
seu <- FindClusters(seu, resolution = RESOLUTION, verbose = FALSE)
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:N_PCS, reduction.name = "umap",
               verbose = FALSE)
seu <- JoinLayers(seu)
seu$cluster <- seu$seurat_clusters

# ---- 5. Neuron vs non-neuron + broad class ------------------------------------
# E18/P0: many neurons are immature (low Rbfox3), so pan-neuronal genes that are
# on from neuronal birth (Stmn2, Tubb3, Elavl3/4, Dcx, Gap43) carry the score.
markers <- list(
  Neuron          = c("Snap25", "Syt1", "Stmn2", "Tubb3", "Elavl3", "Elavl4", "Rbfox3",
                      "Map2", "Gap43", "Dcx", "Ina", "Nefl", "Celf4"),
  Progenitor      = c("Sox2", "Pax6", "Nes", "Hes5", "Hes1", "Fabp7", "Vim", "Notch1"),
  Cycling         = c("Mki67", "Top2a", "Cenpf", "Birc5", "Ccnb1"),
  Astrocyte       = c("Aldh1l1", "Aqp4", "Gja1", "Slc1a3", "Gfap", "Aldoc"),
  OPC_Oligo       = c("Olig1", "Olig2", "Pdgfra", "Sox10", "Cspg4"),
  Microglia_Macro = c("Cx3cr1", "P2ry12", "C1qa", "C1qb", "Csf1r", "Aif1", "Tmem119"),
  Endothelial     = c("Cldn5", "Flt1", "Pecam1", "Kdr", "Esam"),
  Mural           = c("Pdgfrb", "Rgs5", "Vtn", "Kcnj8", "Acta2"),
  VLMC_Fibro      = c("Col1a1", "Col1a2", "Dcn", "Lum", "Igf2"),
  Ependymal_CP    = c("Foxj1", "Ttr", "Kcnj13", "Folr1"),
  Erythrocyte     = c("Hba-a1", "Hba-a2", "Hbb-bs", "Hbb-bt", "Alas2"))
neuron_type_markers <- list(
  Glutamatergic = c("Slc17a7", "Slc17a6", "Neurod6", "Neurod2", "Tbr1", "Satb2"),
  GABAergic     = c("Gad1", "Gad2", "Slc32a1", "Dlx6os1", "Sp9"))

present <- function(g) intersect(g, rownames(seu))
markers <- lapply(markers, present); markers <- markers[lengths(markers) >= 2]
neuron_type_markers <- lapply(neuron_type_markers, present)
write.csv(data.frame(class = rep(names(markers), lengths(markers)), gene = unlist(markers)),
          tab("marker_genes_used.csv"), row.names = FALSE)

seu <- AddModuleScore(seu, features = markers, name = "score_", seed = 1)
score_cols <- paste0("score_", seq_along(markers))
names(seu@meta.data)[match(score_cols, names(seu@meta.data))] <- paste0("score_", names(markers))
score_cols <- paste0("score_", names(markers))

# Class is assigned per cluster (mean module score), which is far more stable
# than per-cell argmax on sparse data. Per-cell scores stay in meta.data.
cl_scores <- seu@meta.data %>% group_by(cluster) %>%
  summarise(across(all_of(score_cols), mean), n_cells = n(), .groups = "drop")
cl_scores$broad_class <- sub("^score_", "", score_cols[max.col(as.matrix(cl_scores[, score_cols]))])
top2 <- apply(as.matrix(cl_scores[, score_cols]), 1, function(x) sort(x, decreasing = TRUE)[1:2])
cl_scores$class_margin <- round(top2[1, ] - top2[2, ], 3)
seu$broad_class <- factor(cl_scores$broad_class[match(seu$cluster, cl_scores$cluster)],
                          levels = names(markers))
seu$neuron_vs_non <- factor(ifelse(seu$broad_class == "Neuron", "Neuron", "Non-neuron"),
                            levels = c("Neuron", "Non-neuron"))

if (all(lengths(neuron_type_markers) >= 2)) {
  seu <- AddModuleScore(seu, features = neuron_type_markers, name = "ntype_", seed = 1)
  nt <- seu@meta.data %>% filter(neuron_vs_non == "Neuron") %>% group_by(cluster) %>%
    summarise(glut = mean(ntype_1), gaba = mean(ntype_2), .groups = "drop") %>%
    mutate(neuron_type = ifelse(glut > gaba, "Glutamatergic", "GABAergic"))
  seu$neuron_type <- ifelse(seu$neuron_vs_non == "Neuron",
                            nt$neuron_type[match(seu$cluster, nt$cluster)], "Non-neuron")
  names(seu@meta.data)[names(seu@meta.data) %in% c("ntype_1", "ntype_2")] <-
    c("score_Glutamatergic", "score_GABAergic")
  cl_scores$neuron_type <- nt$neuron_type[match(cl_scores$cluster, nt$cluster)]
}
write.csv(cl_scores, tab("cluster_class_scores.csv"), row.names = FALSE)

# ---- 6. Figures ----------------------------------------------------------------
pal_class <- c(Neuron = "#C0392B", Progenitor = "#2E86C1", Cycling = "#85C1E9",
               Astrocyte = "#27AE60", OPC_Oligo = "#F39C12", Microglia_Macro = "#8E44AD",
               Endothelial = "#16A085", Mural = "#A04000", VLMC_Fibro = "#7F8C8D",
               Ependymal_CP = "#D4AC0D", Erythrocyte = "#E91E63")
pal_nn <- c(Neuron = "#C0392B", `Non-neuron` = "#5D6D7E")
pt <- if (ncol(seu) > 50000) 0.05 else 0.2
um <- function(group, red = "umap", cols = NULL, label = FALSE, title = group)
  DimPlot(seu, reduction = red, group.by = group, cols = cols, pt.size = pt,
          label = label, repel = TRUE, shuffle = TRUE, raster = FALSE) +
  ggtitle(title) + coord_equal() + NoAxes()

save_plot(um("neuron_vs_non", cols = pal_nn, title = "Neuron vs non-neuron") |
          um("broad_class", cols = pal_class, title = "Broad cell class"),
          "04_umap_neuron_vs_nonneuron", 14, 6)
save_plot(um("cluster", label = TRUE, title = "Clusters (Harmony, res 0.5)") |
          um("stage", title = "Stage") | um("sample", title = "Sample"),
          "05_umap_cluster_stage_sample", 20, 6)
save_plot(DimPlot(seu, reduction = "umap", group.by = "neuron_vs_non", split.by = "sample",
                  cols = pal_nn, pt.size = pt, raster = FALSE) + NoAxes(),
          "06_umap_split_by_sample", 4 * nlevels(seu$sample), 4.5)
save_plot(um("sample", red = "umap.unintegrated", title = "Sample - before Harmony") |
          um("broad_class", red = "umap.unintegrated", cols = pal_class, title = "Class - before Harmony"),
          "07_umap_unintegrated_batch_check", 14, 6)
if ("neuron_type" %in% names(seu@meta.data))
  save_plot(um("neuron_type", cols = c(Glutamatergic = "#C0392B", GABAergic = "#2874A6",
                                       `Non-neuron` = "grey85"), title = "Neuron type (cluster-level)"),
            "08_umap_neuron_type", 7, 6)

key_genes <- present(c("Snap25", "Stmn2", "Tubb3", "Rbfox3", "Slc17a6", "Gad2",
                       "Sox2", "Mki67", "Aqp4", "Olig2", "Cx3cr1", "Cldn5"))
save_plot(FeaturePlot(seu, key_genes, ncol = 4, pt.size = 0.1, order = TRUE, raster = FALSE) &
            NoAxes() & NoLegend(), "09_featureplot_markers", 16, 12)
save_plot(FeaturePlot(seu, score_cols, ncol = 4, pt.size = 0.1, raster = FALSE) & NoAxes(),
          "10_featureplot_class_scores", 16, 4 * ceiling(length(score_cols) / 4))

dot_genes <- unique(unlist(lapply(markers, head, 4)))
Idents(seu) <- "cluster"
p_dot <- DotPlot(seu, features = dot_genes, group.by = "cluster") + RotatedAxis() +
  labs(x = NULL, y = "cluster") + theme(axis.text.x = element_text(size = 8))
save_plot(p_dot, "11_dotplot_markers_by_cluster", 18, 0.3 * nlevels(seu$cluster) + 3)

comp <- seu@meta.data %>% dplyr::count(sample, stage, broad_class) %>%
  group_by(sample) %>% mutate(frac = n / sum(n)) %>% ungroup()
write.csv(comp, tab("composition_per_sample.csv"), row.names = FALSE)
p_bar <- ggplot(comp, aes(sample, frac, fill = broad_class)) + geom_col(width = 0.8) +
  scale_fill_manual(values = pal_class) + scale_y_continuous(labels = scales::percent) +
  facet_grid(~stage, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = "fraction of cells", fill = NULL) + theme_classic()
save_plot(p_bar, "12_composition_per_sample", 8, 5)

# ---- 7. Cluster markers ---------------------------------------------------------
msg("cluster markers (downsampled to 300 cells/cluster)")
mk <- FindAllMarkers(subset(seu, downsample = 300), only.pos = TRUE, logfc.threshold = 0.5,
                     min.pct = 0.2, verbose = FALSE)
mk <- mk %>% filter(p_val_adj < 0.05) %>% arrange(cluster, desc(avg_log2FC))
write.csv(mk, tab("cluster_markers_all.csv"), row.names = FALSE)
top10 <- mk %>% group_by(cluster) %>% slice_head(n = 10) %>%
  summarise(top10 = paste(gene, collapse = ", "), .groups = "drop") %>%
  left_join(cl_scores %>% select(any_of(c("cluster", "n_cells", "broad_class", "class_margin", "neuron_type"))),
            by = "cluster")
write.csv(top10, tab("cluster_top10_markers_with_class.csv"), row.names = FALSE)

# ---- 8. Save ----------------------------------------------------------------------
write.csv(cbind(cell = colnames(seu), seu@meta.data), tab("cell_metadata.csv"), row.names = FALSE)
summ <- seu@meta.data %>% dplyr::count(sample, stage, neuron_vs_non) %>%
  tidyr::pivot_wider(names_from = neuron_vs_non, values_from = n, values_fill = 0)
write.csv(summ, tab("neuron_vs_nonneuron_counts.csv"), row.names = FALSE)
print(as.data.frame(summ))

# scale.data is dropped to keep the file small; counts + log-normalised data,
# PCA, Harmony and both UMAPs are kept.
seu_save <- DietSeurat(seu, layers = c("counts", "data"),
                       dimreducs = c("pca", "harmony", "umap", "umap.unintegrated"))
seu_save@misc$preprocessing <- list(
  date = as.character(Sys.Date()), qc = QC, n_hvg = N_HVG, n_pcs = N_PCS,
  resolution = RESOLUTION, doublets = "scDblFinder per sample, removed",
  integration = "Harmony on sample (reduction 'harmony'); 'umap.unintegrated' = before",
  class_call = "cluster-level argmax of AddModuleScore over marker sets (tables/marker_genes_used.csv)",
  seurat = as.character(packageVersion("Seurat")))
VariableFeatures(seu_save) <- VariableFeatures(seu)
saveRDS(seu_save, file.path(OUT_DIR, "E18_P0_seurat_preprocessed.rds"))
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
msg("done -> ", OUT_DIR)
