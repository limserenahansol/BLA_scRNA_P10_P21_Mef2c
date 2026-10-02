# Figures 1-4: cell types at P10, marker validation, BLA glutamatergic subtypes,
# composition P10 vs P21.
suppressPackageStartupMessages({ library(Seurat); library(dplyr); library(ggplot2); library(patchwork); library(tidyr) })
ROOT <- Sys.getenv("BLA_ROOT", "K:/scRNA_BLA_phd/P10_P21_amygdala_2026-10")
seu <- readRDS(file.path(ROOT, "B06_P10_P21_annotated.rds"))
FIG <- function(n) file.path(ROOT, "figures", n)
sv <- function(p, n, w, h) { ggsave(FIG(paste0(n, ".png")), p, width = w, height = h, dpi = 220, bg = "white")
                             ggsave(FIG(paste0(n, ".pdf")), p, width = w, height = h) }
cts <- levels(seu$cell_type)
pal <- setNames(c(colorRampPalette(c("#7B1E1E", "#E4572E", "#F3A712"))(10),
                  colorRampPalette(c("#1B3A6B", "#2E86AB", "#7FC8F8"))(15),
                  colorRampPalette(c("#2D6A4F", "#95D5B2", "#B7B7A4"))(10)), cts)
write.csv(data.frame(cell_type = cts, color = pal), file.path(ROOT, "tables/palette.csv"), row.names = FALSE)
p10 <- subset(seu, subset = age == "P10")
lab_theme <- theme(plot.title = element_text(size = 13, face = "bold"))

# ---- F1: UMAP ------------------------------------------------------------------
um <- function(o, title, leg = TRUE) {
  o$ct_num <- factor(paste(match(o$cell_type, cts), o$cell_type), levels = paste(seq_along(cts), cts))
  p <- DimPlot(o, group.by = "ct_num", cols = setNames(pal, paste(seq_along(cts), cts)), pt.size = 0.25, raster = FALSE, shuffle = TRUE) +
    ggtitle(title) + coord_equal() + NoAxes() + lab_theme
  if (leg) p + guides(color = guide_legend(ncol = 1, override.aes = list(size = 3))) +
    theme(legend.text = element_text(size = 8)) else p + NoLegend()
}
lab_df <- p10@meta.data %>% mutate(u1 = Embeddings(p10, "umap")[, 1], u2 = Embeddings(p10, "umap")[, 2]) %>%
  group_by(cell_type) %>% summarise(u1 = median(u1), u2 = median(u2), n = n()) %>% filter(n >= 40)
lab_df$short <- match(lab_df$cell_type, cts)
p1 <- um(p10, sprintf("P10 BLA - %s cells, %d cell types", format(ncol(p10), big.mark = ","), nlevels(droplevels(p10$cell_type))), leg = TRUE) +
  ggrepel::geom_text_repel(data = lab_df, aes(u1, u2, label = short), inherit.aes = FALSE, size = 3.4, fontface = "bold",
                           max.overlaps = 40, min.segment.length = 0, segment.size = 0.2)
sv(p1, "F1_P10_umap_cell_types", 13, 8.5)
sv(DimPlot(seu, group.by = "cell_type", split.by = "age", cols = pal, pt.size = 0.2, raster = FALSE) +
     NoAxes() + NoLegend() + ggtitle("Same map, P10 vs P21"), "F1b_umap_split_P10_P21", 12, 6)

# ---- F2: marker validation dot plot (P10) -----------------------------------------
markers <- list(
  Neuron = c("Snap25", "Slc17a7", "Slc17a6", "Gad2"),
  `BLA/BMA Glut` = c("Rspo2", "Etv1", "Lypd1", "Tshz2", "Satb1", "Cdh8", "Fgf10", "Otof", "Esr1", "Satb2", "Nr4a2", "Igfbpl1"),
  ITC = c("Foxp2", "Tshz1", "Pbx3", "Meis2", "Chst9"),
  `CeA/striatal-like` = c("Prkcd", "Calcrl", "Penk", "Tac1", "Drd1", "Drd2", "Adora2a", "Ebf1"),
  Interneuron = c("Sst", "Pvalb", "Maf", "Vip", "Cck", "Cnr1", "Reln", "Lamp5", "Npy"),
  Glia = c("Slc1a3", "Hexb", "Pdgfra", "Top2a", "Bcas1", "Plp1", "Flt1", "Vtn", "Tmem212"))
markers <- lapply(markers, intersect, rownames(seu))
Idents(p10) <- "cell_type"; Idents(p10) <- factor(Idents(p10), levels = rev(levels(droplevels(p10$cell_type))))
p2 <- DotPlot(p10, features = markers, cols = c("grey92", "#B2182B"), dot.scale = 4.5, cluster.idents = FALSE) +
  RotatedAxis() + labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_text(size = 7.5), axis.text.y = element_text(size = 8.5),
        strip.text = element_text(size = 8, face = "bold"), legend.position = "bottom") +
  ggtitle("Known amygdala markers (P10)") + lab_theme
sv(p2, "F2_P10_marker_dotplot", 17, 10)

# Allen mapping agreement: for each cell type, top Allen subclasses of its cells
al <- p10@meta.data %>% count(cell_type, allen_subclass) %>% group_by(cell_type) %>%
  mutate(frac = n / sum(n)) %>% slice_max(frac, n = 2, with_ties = FALSE) %>%
  summarise(allen_top = paste0(allen_subclass, " (", round(100 * frac), "%)", collapse = "; "))
write.csv(al, file.path(ROOT, "tables/B07_celltype_vs_allen_subclass_P10.csv"), row.names = FALSE)

# ---- F3: BLA glutamatergic sub-populations --------------------------------------
gl <- subset(seu, subset = class == "Glutamatergic" & !is.na(glut_umap1))
gl[["gumap"]] <- CreateDimReducObject(as.matrix(gl@meta.data[, c("glut_umap1", "glut_umap2")]) |>
                                        `colnames<-`(c("gumap_1", "gumap_2")), key = "gumap_", assay = "RNA")
gl$cell_type <- droplevels(gl$cell_type)
pa <- DimPlot(gl, reduction = "gumap", group.by = "cell_type", cols = pal, pt.size = 0.5, label = FALSE, raster = FALSE) +
  coord_equal() + NoAxes() + ggtitle("Glutamatergic neurons (P10 + P21)") + lab_theme +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3))) + theme(legend.text = element_text(size = 8))
pb <- DimPlot(gl, reduction = "gumap", group.by = "age", cols = c(P10 = "#E4572E", P21 = "#2E86AB"), pt.size = 0.4,
              shuffle = TRUE, raster = FALSE) + coord_equal() + NoAxes() + ggtitle("Age") + lab_theme
gmk <- intersect(c("Slc17a7", "Slc17a6", "Rspo2", "Etv1", "Cdh9", "Sema3e", "Tshz2", "Satb1", "Lypd1", "Cdh8", "Tshz3",
                   "Fgf10", "Nrp1", "Otof", "Trhr", "Col23a1", "Vgll3", "Prr16", "Zfp804b", "Scn5a", "Esr1", "Reln",
                   "Igfbpl1", "Epha3", "Satb2", "Nr4a2", "Ppp1r1b", "Fezf2", "Mef2c"), rownames(gl))
Idents(gl) <- "cell_type"
pc <- DotPlot(gl, features = gmk, cols = c("grey92", "#B2182B"), dot.scale = 5) + RotatedAxis() +
  labs(x = NULL, y = NULL) + theme(axis.text.x = element_text(size = 8.5, face = "italic"), axis.text.y = element_text(size = 9))
sv((pa | pb) / pc + plot_layout(heights = c(1.3, 1)), "F3_BLA_glut_subtypes", 15, 11)
pf <- FeaturePlot(gl, reduction = "gumap", features = c("Rspo2", "Etv1", "Tshz2", "Otof", "Fgf10", "Esr1", "Slc17a6", "Ppp1r1b", "Mef2c"),
                  ncol = 3, pt.size = 0.3, order = TRUE, raster = FALSE) & NoAxes() & coord_equal()
sv(pf, "F3b_BLA_glut_featureplots", 12, 12)
sup <- gl@meta.data %>% count(cell_type, allen_supertype) %>% group_by(cell_type) %>% mutate(frac = n / sum(n)) %>%
  slice_max(frac, n = 3, with_ties = FALSE)
write.csv(sup, file.path(ROOT, "tables/B07_glut_subtype_vs_allen_supertype.csv"), row.names = FALSE)

# ---- F4: composition P10 vs P21 ----------------------------------------------------
comp <- seu@meta.data %>% count(sample, age, class, cell_type) %>% group_by(sample) %>%
  mutate(frac_all = n / sum(n)) %>% group_by(sample, class) %>% mutate(frac_class = n / sum(n)) %>% ungroup()
write.csv(comp, file.path(ROOT, "tables/B07_composition_per_sample.csv"), row.names = FALSE)
p4a <- ggplot(comp %>% distinct(sample, age, class, .keep_all = FALSE) %>%
                left_join(comp %>% group_by(sample, class) %>% summarise(f = sum(frac_all), .groups = "drop"), by = c("sample", "class")),
              aes(sample, f, fill = class)) + geom_col(width = 0.75) +
  scale_fill_manual(values = c(Glutamatergic = "#C0392B", GABAergic = "#2874A6", `Non-neuronal` = "#7F8C8D")) +
  scale_y_continuous(labels = scales::percent) + labs(x = NULL, y = "fraction of all cells", fill = NULL) + theme_classic()
neur <- comp %>% filter(class != "Non-neuronal") %>% group_by(sample) %>% mutate(f = n / sum(n)) %>% ungroup()
fc <- neur %>% group_by(cell_type, age) %>% summarise(f = mean(f), .groups = "drop") %>%
  pivot_wider(names_from = age, values_from = f, values_fill = 0) %>%
  mutate(log2_P21_vs_P10 = log2((P21 + 1e-3) / (P10 + 1e-3)))
write.csv(fc, file.path(ROOT, "tables/B07_neuron_composition_P10_vs_P21.csv"), row.names = FALSE)
p4b <- ggplot(neur, aes(cell_type, f, color = age)) + geom_point(size = 2.2, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c(P10 = "#E4572E", P21 = "#2E86AB")) + scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "fraction of neurons", color = NULL, subtitle = "dot = sample (P10 n=1, P21 n=2)") +
  theme_classic() + theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8))
sv(p4a + p4b + plot_layout(widths = c(1, 4)), "F4_composition_P10_vs_P21", 16, 6)
cat("figures 1-4 done\n")
