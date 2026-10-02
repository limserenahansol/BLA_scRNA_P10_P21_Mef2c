# Scenario 2 (refinement / pruning) and scenario 3 (cue - receptor) read-outs on
# the amygdala (target) side, P10 vs P21.
suppressPackageStartupMessages({ library(Seurat); library(dplyr); library(ggplot2); library(patchwork); library(tidyr); library(Matrix) })
ROOT <- Sys.getenv("BLA_ROOT", "K:/scRNA_BLA_phd/P10_P21_amygdala_2026-10")
seu <- readRDS(file.path(ROOT, "B06_P10_P21_annotated.rds"))
cts <- levels(seu$cell_type); X <- LayerData(seu, "data")
FIG <- function(n) file.path(ROOT, "figures", n); TAB <- function(n) file.path(ROOT, "tables", n)
sv <- function(p, n, w, h) { ggsave(FIG(paste0(n, ".png")), p, width = w, height = h, dpi = 220, bg = "white")
                             ggsave(FIG(paste0(n, ".pdf")), p, width = w, height = h) }

stat <- function(genes, types) {
  genes <- intersect(genes, rownames(seu))
  cells <- seu$cell_type %in% types
  key <- interaction(droplevels(seu$cell_type[cells]), seu$age[cells], seu$sample[cells], sep = "|", drop = TRUE)
  Xs <- X[genes, cells, drop = FALSE]
  res <- lapply(levels(key), function(k) {
    ix <- which(key == k)
    data.frame(gene = genes, key = k, pct = Matrix::rowMeans(Xs[, ix, drop = FALSE] > 0),
               avg = Matrix::rowMeans(Xs[, ix, drop = FALSE]), n = length(ix))
  })
  bind_rows(res) %>% separate(key, c("cell_type", "age", "sample"), sep = "\\|") %>%
    mutate(gene = factor(gene, levels = genes), cell_type = factor(cell_type, levels = rev(types)))
}
agg_age <- function(d) d %>% group_by(gene, cell_type, age) %>%
  summarise(pct = weighted.mean(pct, n), avg = weighted.mean(avg, n), n = sum(n), .groups = "drop")
dot_age <- function(d, title, sub = NULL) {
  d <- agg_age(d) %>% group_by(gene) %>% mutate(z = (avg - mean(avg)) / (sd(avg) + 1e-9)) %>% ungroup()
  ggplot(d, aes(gene, cell_type, size = 100 * pct, color = pmin(pmax(z, -2.5), 2.5))) + geom_point() + facet_grid(~age) +
    scale_size_area(max_size = 5, name = "% cells") +
    scale_color_gradient2(low = "#2166AC", mid = "grey90", high = "#B2182B", name = "z") +
    labs(x = NULL, y = NULL, title = title, subtitle = sub) + theme_bw(base_size = 9.5) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, face = "italic"), panel.grid = element_line(linewidth = 0.2),
          plot.title = element_text(face = "bold"))
}
# replicate-aware change: P10s1 vs each P21 sample, mean log-expression
change <- function(d) d %>% select(gene, cell_type, sample, avg, pct) %>%
  pivot_wider(names_from = sample, values_from = c(avg, pct)) %>%
  mutate(lfc_vs_P21s1 = (avg_P10s1 - avg_P21s1) / log(2), lfc_vs_P21s2 = (avg_P10s1 - avg_P21s2) / log(2),
         consistent = sign(lfc_vs_P21s1) == sign(lfc_vs_P21s2) & pmin(abs(lfc_vs_P21s1), abs(lfc_vs_P21s2)) > 0.3)

# ---- Scenario 2: refinement machinery ------------------------------------------------
glia <- c("Microglia", "Astrocyte")
prune_genes <- c("C1qa", "C1qb", "C1qc", "C3", "C4b", "Itgam", "Trem2", "Tyrobp", "Cd68", "Cx3cr1", "P2ry12",
                 "Csf1r", "Mertk", "Megf10", "Axl", "Gas6", "Mef2c", "Mef2a")
d_pr <- stat(prune_genes, glia)
mef2_targets <- c("Mef2c", "Arc", "Homer1", "Nr4a1", "Bdnf", "Pcdh10", "Npas4", "Fos", "Egr1", "Ncam1", "Sema6d", "Cadm2", "Pcdh7", "Grip1")
neur_focus <- c(cts[1:6], "ITC Foxp2/Tshz1 (a)", "ITC Foxp2/Kcnh5 (b)")
d_m2 <- stat(mef2_targets, neur_focus)
write.csv(bind_rows(change(d_pr), change(d_m2)), TAB("B10_scenario2_genes_P10_vs_P21.csv"), row.names = FALSE)
p2a <- dot_age(d_pr, "a  Pruning machinery in BLA glia",
               "complement (C1q, C3, C4b), microglial receptors (Itgam/CR3, Trem2, Cx3cr1) and astrocyte phagocytosis (Megf10, Mertk)")
p2b <- dot_age(d_m2, "b  Target neurons: MEF2 activity-dependent genes and the P10-high adhesion programme",
               "Arc / Homer1 / Nr4a1 / Bdnf / Pcdh10 = MEF2-regulated synapse-elimination genes; Ncam1 / Sema6d / Cadm2 / Pcdh7 = shared P10-high genes (F5)")
# per-sample bars for the key glial genes so the single P10 animal is visible
key <- d_pr %>% filter(gene %in% c("C1qa", "C1qb", "C4b", "C3", "Trem2", "Megf10", "Mertk"))
p2c <- ggplot(key, aes(sample, 100 * pct, fill = age)) + geom_col(width = 0.7) + facet_grid(cell_type ~ gene) +
  scale_fill_manual(values = c(P10 = "#E4572E", P21 = "#2E86AB")) +
  labs(x = NULL, y = "% cells", fill = NULL, title = "c  Per sample (P10 n=1, P21 n=2)") + theme_bw(base_size = 8.5) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), plot.title = element_text(face = "bold"))
sv((p2a | p2c) / p2b + plot_layout(heights = c(1, 1.25)), "F10_scenario2_pruning_window", 17, 11)

# ---- Scenario 3: cues on the BLA side ---------------------------------------------------
cue_genes <- c("Sema3e", "Plxnd1", "Nrp1", "Kirrel3", "Sdk2", "Sema6d", "Sema5b", "Cdh8", "Cdh9", "Cdh13", "Cdh18", "C1ql3", "Slit2", "Ntng1")
types3 <- c(cts[1:10], cts[11:16])
d3 <- stat(cue_genes, types3)
ch3 <- change(d3); write.csv(ch3, TAB("B10_scenario3_cues_P10_vs_P21.csv"), row.names = FALSE)
p3a <- dot_age(d3, "a  BLA-side cues and local receptors, P10 vs P21",
               "Sema3e (Rspo2 BLA) signals through Plxnd1, and Nrp1 decides repulsion vs attraction; Kirrel3 / Sdk2 / cadherins are homophilic")
hm <- ch3 %>% filter(cell_type %in% neur_focus) %>% mutate(lfc = (lfc_vs_P21s1 + lfc_vs_P21s2) / 2,
                                                          cell_type = factor(cell_type, levels = rev(neur_focus)))
p3b <- ggplot(hm, aes(gene, cell_type, fill = pmax(pmin(lfc, 1.5), -1.5))) + geom_tile(color = "white") +
  geom_text(aes(label = ifelse(consistent, "*", "")), size = 4) +
  scale_fill_gradient2(low = "#2E86AB", mid = "white", high = "#E4572E", name = "log2\nP10 / P21") +
  labs(x = NULL, y = NULL, title = "b  Which cues peak at P10 (ectopic-innervation window)?",
       subtitle = "* = same direction and |log2FC| > 0.3 against BOTH P21 samples") +
  theme_minimal(base_size = 9.5) + theme(axis.text.x = element_text(angle = 60, hjust = 1, face = "italic"), plot.title = element_text(face = "bold"))
sv(p3a / p3b + plot_layout(heights = c(1.5, 1)), "F11_scenario3_BLA_cues", 13, 12)
print(hm %>% filter(consistent) %>% arrange(desc(abs(lfc))) %>% select(cell_type, gene, lfc, pct_P10s1, pct_P21s1, pct_P21s2) %>% as.data.frame() %>% head(40))
print(agg_age(d_pr) %>% filter(gene %in% c("C1qa", "C4b", "C3", "Trem2", "Megf10", "Mertk")) %>% as.data.frame())
