# Figures 5-8: P10 vs P21 DE per cell type; transcription factors (Mef2 family,
# cell-type TFs); axon-guidance / adhesion cues expressed by BLA cells during
# the P5-P10 window of ectopic Mef2c-cKO innervation; Mef2c co-varying genes.
suppressPackageStartupMessages({ library(Seurat); library(dplyr); library(ggplot2); library(patchwork); library(tidyr); library(Matrix) })
set.seed(1)
ROOT <- Sys.getenv("BLA_ROOT", "K:/scRNA_BLA_phd/P10_P21_amygdala_2026-10")
seu <- readRDS(file.path(ROOT, "B06_P10_P21_annotated.rds"))
pal <- deframe <- NULL
pal <- with(read.csv(file.path(ROOT, "tables/palette.csv")), setNames(color, cell_type))
FIG <- function(n) file.path(ROOT, "figures", n); TAB <- function(n) file.path(ROOT, "tables", n)
sv <- function(p, n, w, h) { ggsave(FIG(paste0(n, ".png")), p, width = w, height = h, dpi = 220, bg = "white")
                             ggsave(FIG(paste0(n, ".pdf")), p, width = w, height = h) }
cts <- levels(seu$cell_type)
counts <- LayerData(seu, "counts"); X <- LayerData(seu, "data")

# pseudo-bulk CPM per sample x cell type
grp <- interaction(seu$cell_type, seu$sample, sep = "|", drop = TRUE)
G <- sparseMatrix(i = seq_len(ncol(seu)), j = as.integer(grp), x = 1, dims = c(ncol(seu), nlevels(grp)))
pb <- as.matrix(counts %*% G); colnames(pb) <- levels(grp)
cpm <- t(t(pb) / colSums(pb)) * 1e6
ncell <- table(grp)

# per-type x age: % expressing and mean log-normalised
type_age <- interaction(seu$cell_type, seu$age, sep = "|", drop = TRUE)
GA <- sparseMatrix(i = seq_len(ncol(seu)), j = as.integer(type_age), x = 1, dims = c(ncol(seu), nlevels(type_age)))
nA <- colSums(GA)
B <- X; B@x[] <- 1
pctA <- as.matrix(B %*% GA) %*% diag(1 / nA); colnames(pctA) <- levels(type_age)
avgA <- as.matrix(X %*% GA) %*% diag(1 / nA); colnames(avgA) <- levels(type_age)
long_stat <- function(genes) {
  genes <- intersect(genes, rownames(seu))
  data.frame(gene = rep(genes, ncol(pctA)), key = rep(colnames(pctA), each = length(genes)),
             pct = as.vector(pctA[genes, ]), avg = as.vector(avgA[genes, ])) %>%
    separate(key, c("cell_type", "age"), sep = "\\|") %>%
    mutate(cell_type = factor(cell_type, levels = rev(cts)), gene = factor(gene, levels = genes))
}
dotp <- function(d, title, split_age = TRUE) {
  d <- d %>% group_by(gene) %>% mutate(z = (avg - mean(avg)) / (sd(avg) + 1e-9)) %>% ungroup()
  p <- ggplot(d, aes(gene, cell_type, size = 100 * pct, color = z)) + geom_point() +
    scale_size_area(max_size = 4.5, name = "% cells") + scale_color_gradient2(low = "#2166AC", mid = "grey90", high = "#B2182B", name = "z") +
    labs(x = NULL, y = NULL, title = title) + theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, face = "italic"), panel.grid = element_line(linewidth = 0.2))
  if (split_age) p + facet_grid(~age) else p
}

# ---- F5: P10 vs P21 DE per cell type ------------------------------------------------
de_all <- list()
for (ct in cts) {
  c10 <- colnames(seu)[seu$cell_type == ct & seu$age == "P10"]; c21 <- colnames(seu)[seu$cell_type == ct & seu$age == "P21"]
  if (length(c10) < 30 || length(c21) < 30) next
  m <- FindMarkers(seu, ident.1 = c10, ident.2 = c21, logfc.threshold = 0.5, min.pct = 0.2,
                   max.cells.per.ident = 400, verbose = FALSE)
  m$gene <- rownames(m); m$cell_type <- ct
  k10 <- paste0(ct, "|P10s1"); k1 <- paste0(ct, "|P21s1"); k2 <- paste0(ct, "|P21s2")
  if (all(c(k10, k1, k2) %in% colnames(cpm))) {
    m$lfc_vs_P21s1 <- log2((cpm[m$gene, k10] + 10) / (cpm[m$gene, k1] + 10))
    m$lfc_vs_P21s2 <- log2((cpm[m$gene, k10] + 10) / (cpm[m$gene, k2] + 10))
    m$consistent <- m$p_val_adj < 0.01 & sign(m$lfc_vs_P21s1) == sign(m$avg_log2FC) &
      sign(m$lfc_vs_P21s2) == sign(m$avg_log2FC) & abs(m$lfc_vs_P21s1) > 0.4 & abs(m$lfc_vs_P21s2) > 0.4
  } else m$consistent <- m$p_val_adj < 0.01
  de_all[[ct]] <- m
}
de <- bind_rows(de_all) %>% mutate(direction = ifelse(avg_log2FC > 0, "higher at P10", "higher at P21"))
write.csv(de, TAB("B08_DE_P10_vs_P21_all_celltypes.csv"), row.names = FALSE)
des <- de %>% filter(consistent)
write.csv(des, TAB("B08_DE_P10_vs_P21_consistent.csv"), row.names = FALSE)
cnt <- des %>% count(cell_type, direction) %>% mutate(cell_type = factor(cell_type, levels = rev(cts)))
p5a <- ggplot(cnt, aes(cell_type, ifelse(direction == "higher at P10", n, -n), fill = direction)) + geom_col() + coord_flip() +
  scale_fill_manual(values = c(`higher at P10` = "#E4572E", `higher at P21` = "#2E86AB")) +
  labs(x = NULL, y = "number of DE genes (consistent vs both P21 replicates)", fill = NULL) + theme_classic(base_size = 9) +
  theme(legend.position = "top")
shared <- des %>% count(gene, direction) %>% arrange(desc(n)) %>% group_by(direction) %>% slice_head(n = 15)
write.csv(shared, TAB("B08_DE_shared_across_celltypes.csv"), row.names = FALSE)
neur_ct <- cts[1:25]
hm <- long_stat(shared$gene) %>% filter(cell_type %in% neur_ct) %>% select(gene, cell_type, age, avg) %>%
  pivot_wider(names_from = age, values_from = avg) %>% mutate(lfc = (P10 - P21) / log(2))
p5b <- ggplot(hm, aes(gene, cell_type, fill = pmax(pmin(lfc, 3), -3))) + geom_tile() +
  scale_fill_gradient2(low = "#2E86AB", mid = "white", high = "#E4572E", name = "log2 P10/P21\n(mean expr)") +
  labs(x = NULL, y = NULL, title = "Genes changing P10 -> P21 in many neuron types") + theme_minimal(base_size = 8.5) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, face = "italic"))
sv(p5a + p5b + plot_layout(widths = c(1, 2)), "F5_DE_P10_vs_P21", 17, 8)

# ---- F6: transcription factors ---------------------------------------------------------
tf <- read.csv(file.path(ROOT, "resources", "mouse_TF_list_1321.csv"))
tfs <- intersect(tf$gene, rownames(seu))
mef2 <- c("Mef2a", "Mef2b", "Mef2c", "Mef2d")
sv(dotp(long_stat(c(mef2, "Satb2", "Bcl11b", "Fezf2", "Cux2", "Tbr1", "Foxp2", "Foxp1", "Tshz1", "Meis2")),
        "Mef2 family and neuronal identity TFs, P10 vs P21"), "F6a_Mef2_family_dotplot", 10, 9.5)
p10 <- subset(seu, subset = age == "P10"); Idents(p10) <- "cell_type"
tfm <- FindAllMarkers(p10, features = tfs, only.pos = TRUE, logfc.threshold = 0.4, min.pct = 0.25,
                      max.cells.per.ident = 200, verbose = FALSE) %>% filter(p_val_adj < 0.01)
tfm$family <- tf$family[match(tfm$gene, tf$gene)]
write.csv(tfm, TAB("B08_TF_markers_per_celltype_P10.csv"), row.names = FALSE)
top_tf <- tfm %>% mutate(s = avg_log2FC * (pct.1 - pct.2)) %>% group_by(cluster) %>% slice_max(s, n = 3) %>% ungroup()
gl <- unique(top_tf$gene[order(match(top_tf$cluster, cts))])
z <- avgA[gl, grep("\\|P10$", colnames(avgA))]; colnames(z) <- sub("\\|P10$", "", colnames(z))
z <- t(scale(t(z))); zd <- as.data.frame(as.table(z)); names(zd) <- c("gene", "cell_type", "z")
zd$cell_type <- factor(zd$cell_type, levels = rev(cts)); zd$gene <- factor(zd$gene, levels = gl)
p6b <- ggplot(zd, aes(gene, cell_type, fill = pmin(z, 3))) + geom_tile() +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", name = "z") +
  labs(x = NULL, y = NULL, title = "Top cell-type-specific TFs at P10 (from 1,321 curated mouse TFs)") +
  theme_minimal(base_size = 8) + theme(axis.text.x = element_text(angle = 70, hjust = 1, face = "italic", size = 6.5))
sv(p6b, "F6b_celltype_TFs_P10", 18, 8)

# ---- F7: guidance / adhesion cues in BLA cells ------------------------------------------
cues <- list(
  `Semaphorin` = c("Sema3a", "Sema3c", "Sema3d", "Sema3e", "Sema3f", "Sema5a", "Sema5b", "Sema6a", "Sema6d", "Sema7a"),
  `Ephrin` = c("Efna1", "Efna2", "Efna3", "Efna5", "Efnb1", "Efnb2", "Efnb3"),
  `Slit/Netrin` = c("Slit1", "Slit2", "Slit3", "Ntn1", "Ntn4", "Ntng1", "Ntng2"),
  `Secreted other` = c("Cxcl12", "Reln", "Wnt5a", "Wnt7b", "Fgf10", "Bdnf", "Ntf3", "Nrg1", "Nrg3", "Cbln1", "Cbln2", "Cbln4", "C1ql2", "C1ql3"),
  `Adhesion / synaptic organiser` = c("Flrt2", "Flrt3", "Tenm2", "Tenm3", "Tenm4", "Lrrtm4", "Nlgn1", "Cdh6", "Cdh8", "Cdh9", "Cdh10", "Cdh11",
                                      "Cdh12", "Cdh13", "Cdh18", "Pcdh7", "Pcdh9", "Pcdh10", "Pcdh17", "Kirrel3", "Sdk1", "Sdk2",
                                      "Cntn4", "Cntn5", "Cntn6", "Nectin3", "Igsf21"))
cue_df <- data.frame(gene = unlist(cues), family = rep(names(cues), lengths(cues))) %>% filter(gene %in% rownames(seu))
receptors <- c(Sema3a = "Nrp1/Plxna1-4", Sema3c = "Nrp1-2/Plxnd1", Sema3d = "Nrp1/Plxnd1", Sema3e = "Plxnd1", Sema3f = "Nrp2/Plxna3",
  Sema5a = "Plxna1/a3", Sema5b = "Plxna1/a3", Sema6a = "Plxna2/a4", Sema6d = "Plxna1", Sema7a = "Plxnc1/Itgb1",
  Efna1 = "EphA", Efna2 = "EphA", Efna3 = "EphA", Efna5 = "EphA (Epha3-7)", Efnb1 = "EphB", Efnb2 = "EphB/Epha4", Efnb3 = "EphB/Epha4",
  Slit1 = "Robo1-3", Slit2 = "Robo1-3", Slit3 = "Robo1-3", Ntn1 = "Dcc/Unc5a-d/Neo1", Ntn4 = "Unc5/Neo1", Ntng1 = "Lrrc4c (NGL-1)", Ntng2 = "Lrrc4 (NGL-2)",
  Cxcl12 = "Cxcr4/Ackr3", Reln = "Vldlr/Lrp8", Wnt5a = "Ryk/Fzd", Wnt7b = "Fzd", Fgf10 = "Fgfr2", Bdnf = "Ntrk2", Ntf3 = "Ntrk3/Ntrk2",
  Nrg1 = "Erbb4", Nrg3 = "Erbb4", Cbln1 = "Nrxn + Grid2", Cbln2 = "Nrxn + Grid1", Cbln4 = "Nrxn + Dcc/Neo1", C1ql2 = "Adgrb3", C1ql3 = "Adgrb3",
  Flrt2 = "Unc5/Adgrl", Flrt3 = "Unc5/Adgrl", Tenm2 = "Adgrl1-3", Tenm3 = "Tenm3/Adgrl", Tenm4 = "Adgrl", Lrrtm4 = "Nrxn/Ptprs",
  Nlgn1 = "Nrxn1-3", Kirrel3 = "Kirrel3 (homophilic)", Sdk1 = "Sdk1 (homophilic)", Sdk2 = "Sdk2 (homophilic)", Nectin3 = "Nectin1",
  Igsf21 = "Nrxn2", Cntn4 = "Ptprg/App", Cntn5 = "Ptprg", Cntn6 = "Ptprg/Chl1")
cue_df$partner_on_axon <- ifelse(cue_df$gene %in% names(receptors), receptors[cue_df$gene],
                                 ifelse(grepl("^Cdh|^Pcdh", cue_df$gene), "homophilic (same cadherin/protocadherin)", ""))
targets <- cts[1:25]
cs <- long_stat(cue_df$gene) %>% filter(cell_type %in% targets) %>% left_join(cue_df, by = "gene")
cs$family <- factor(cs$family, levels = names(cues))
# developmental change per cell type and enrichment in a type vs other neurons, at P10
wide <- cs %>% select(gene, cell_type, age, pct, avg) %>% pivot_wider(names_from = age, values_from = c(pct, avg))
wide <- wide %>% group_by(gene) %>% mutate(pct_other_max_P10 = sapply(seq_along(pct_P10), function(i) max(pct_P10[-i])),
                                           spec_pp_P10 = 100 * (pct_P10 - pct_other_max_P10),
                                           log2_P10_vs_P21 = (avg_P10 - avg_P21) / log(2)) %>% ungroup() %>%
  left_join(cue_df, by = "gene")
write.csv(wide, TAB("B08_guidance_cues_by_celltype_P10_P21.csv"), row.names = FALSE)
p7 <- ggplot(cs %>% filter(age == "P10") %>% group_by(gene) %>% mutate(z = (avg - mean(avg)) / (sd(avg) + 1e-9)),
             aes(gene, cell_type, size = 100 * pct, color = pmin(z, 3))) + geom_point() +
  facet_grid(~family, scales = "free_x", space = "free_x") +
  scale_size_area(max_size = 4, name = "% cells") + scale_color_gradient2(low = "#2166AC", mid = "grey90", high = "#B2182B", name = "z") +
  labs(x = NULL, y = NULL, title = "Axon guidance / adhesion cues expressed by amygdala neurons at P10") +
  theme_bw(base_size = 8.5) + theme(axis.text.x = element_text(angle = 70, hjust = 1, face = "italic"),
                                    strip.text = element_text(size = 7.5, face = "bold"), panel.grid = element_line(linewidth = 0.15))
sv(p7, "F7a_guidance_cues_P10", 19, 8.5)
focus <- c("ITC Foxp2/Tshz1 (a)", "ITC Foxp2/Kcnh5 (b)", cts[1:6])
cand <- wide %>% filter(cell_type %in% focus, pct_P10 >= 0.3) %>%
  arrange(desc(spec_pp_P10)) %>% group_by(cell_type) %>% slice_head(n = 6) %>% ungroup() %>%
  select(cell_type, gene, family, partner_on_axon, pct_P10, pct_P21, spec_pp_P10, log2_P10_vs_P21)
write.csv(cand, TAB("B08_candidate_cues_ITC_BLAglut.csv"), row.names = FALSE)
dyn <- wide %>% filter(cell_type %in% focus) %>% mutate(cell_type = factor(cell_type, levels = rev(focus)))
p7b <- ggplot(dyn, aes(gene, cell_type, fill = pmax(pmin(log2_P10_vs_P21, 2), -2))) + geom_tile(color = "white") +
  facet_grid(~family, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(low = "#2E86AB", mid = "white", high = "#E4572E", name = "log2\nP10/P21") +
  labs(x = NULL, y = NULL, title = "Which cues are higher at P10 (ectopic-innervation window) than at P21") +
  theme_minimal(base_size = 8.5) + theme(axis.text.x = element_text(angle = 70, hjust = 1, face = "italic"),
                                         strip.text = element_text(size = 7.5, face = "bold"))
sv(p7b, "F7b_guidance_cues_P10_vs_P21", 19, 4.5)

# ---- F8: Mef2c ----------------------------------------------------------------------------
mc <- long_stat("Mef2c") %>% mutate(cell_type = factor(cell_type, levels = cts))
p8a <- ggplot(mc, aes(cell_type, 100 * pct, fill = age)) + geom_col(position = position_dodge(0.8), width = 0.75) +
  scale_fill_manual(values = c(P10 = "#E4572E", P21 = "#2E86AB")) +
  labs(x = NULL, y = "% cells with Mef2c", fill = NULL, title = "Mef2c in amygdala cell types") +
  theme_classic(base_size = 9) + theme(axis.text.x = element_text(angle = 60, hjust = 1))
# genes co-varying with Mef2c inside BLA glutamatergic neurons at P10 (depth regressed out)
bg <- colnames(seu)[seu$cell_type %in% cts[1:6] & seu$age == "P10"]
Xb <- X[, bg]; keepg <- rownames(Xb)[Matrix::rowMeans(Xb > 0) >= 0.15]
Xm <- as.matrix(Xb[keepg, ]); dep <- log10(seu$nCount_RNA[match(bg, colnames(seu))])
ct_f <- factor(seu$cell_type[match(bg, colnames(seu))])
D <- model.matrix(~ dep + ct_f)
res <- Xm - t(D %*% solve(crossprod(D), crossprod(D, t(Xm))))
r <- cor(t(res), res["Mef2c", ])[, 1]
cor_df <- data.frame(gene = names(r), r = r, pct = Matrix::rowMeans(Xb[keepg, ] > 0)) %>% filter(gene != "Mef2c") %>% arrange(desc(r))
cor_df$is_TF <- cor_df$gene %in% tfs; cor_df$is_cue <- cor_df$gene %in% cue_df$gene
write.csv(cor_df, TAB("B08_Mef2c_covarying_genes_BLAglut_P10.csv"), row.names = FALSE)
top <- bind_rows(head(cor_df, 25) %>% mutate(side = "positive"), tail(cor_df, 10) %>% mutate(side = "negative"))
p8b <- ggplot(top, aes(reorder(gene, r), r, fill = ifelse(is_cue, "guidance/adhesion", ifelse(is_TF, "TF", "other")))) +
  geom_col() + coord_flip() + scale_fill_manual(values = c(`guidance/adhesion` = "#E4572E", TF = "#6C3483", other = "grey60"), name = NULL) +
  labs(x = NULL, y = "correlation with Mef2c (residual, depth + subtype removed)",
       title = "Genes co-varying with Mef2c in P10 BLA glutamatergic neurons") + theme_classic(base_size = 9)
sv(p8a / p8b + plot_layout(heights = c(1, 1.6)), "F8_Mef2c", 12, 12)
cat("F5-F8 done; consistent DE genes:", nrow(des), "\n")
print(cand, n = 60)
