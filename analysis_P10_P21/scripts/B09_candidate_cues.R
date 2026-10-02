# Rank guidance/adhesion cues that mark ITC and BLA glutamatergic sub-populations
# at P10, against the MEAN of the other amygdala neuron types (the max-based
# score in B08 is uninformative for long, near-ubiquitous neuronal genes).
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(tidyr) })
ROOT <- Sys.getenv("BLA_ROOT", "K:/scRNA_BLA_phd/P10_P21_amygdala_2026-10")
w <- read.csv(file.path(ROOT, "tables/B08_guidance_cues_by_celltype_P10_P21.csv"), check.names = FALSE)
cts <- read.csv(file.path(ROOT, "tables/palette.csv"))$cell_type
focus <- c(cts[1:6], "ITC Foxp2/Tshz1 (a)", "ITC Foxp2/Kcnh5 (b)")
w <- w %>% group_by(gene) %>%
  mutate(other_avg = (sum(avg_P10) - avg_P10) / (n() - 1), other_pct = (sum(pct_P10) - pct_P10) / (n() - 1),
         enrich_log2 = log2((expm1(avg_P10) + 0.05) / (expm1(other_avg) + 0.05)),
         pct_minus_other = 100 * (pct_P10 - other_pct)) %>% ungroup()
cand <- w %>% filter(cell_type %in% focus, pct_P10 >= 0.3, enrich_log2 >= 0.7, pct_minus_other >= 10) %>%
  arrange(cell_type, desc(enrich_log2)) %>% group_by(cell_type) %>% slice_head(n = 6) %>% ungroup() %>%
  mutate(trend = case_when(log2_P10_vs_P21 > 0.5 ~ "higher at P10", log2_P10_vs_P21 < -0.5 ~ "higher at P21", TRUE ~ "stable"))
write.csv(cand %>% select(cell_type, gene, family, partner_on_axon, pct_P10, pct_P21, pct_minus_other, enrich_log2, log2_P10_vs_P21, trend),
          file.path(ROOT, "tables/B09_candidate_cues_ranked.csv"), row.names = FALSE)
genes <- unique(cand$gene)
d <- w %>% filter(cell_type %in% focus, gene %in% genes) %>%
  mutate(cell_type = factor(cell_type, levels = rev(focus)), gene = factor(gene, levels = genes))
p <- ggplot(d, aes(gene, cell_type, size = 100 * pct_P10, color = pmax(pmin(enrich_log2, 3), -1))) + geom_point() +
  geom_point(data = d %>% semi_join(cand, by = c("gene", "cell_type")), shape = 21, color = "black", stroke = 0.6, fill = NA) +
  scale_size_area(max_size = 6, name = "% cells (P10)") +
  scale_color_gradient2(low = "#2166AC", mid = "grey90", high = "#B2182B", midpoint = 0, name = "log2 enrichment\nvs other neurons") +
  labs(x = NULL, y = NULL, title = "Guidance / adhesion cues that mark ITC and BLA glutamatergic sub-populations at P10",
       subtitle = "circled = passes >=30% of cells, >=+10 pp and >=1.6x over the mean of the other amygdala neuron types") +
  theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 60, hjust = 1, face = "italic"), panel.grid = element_line(linewidth = 0.2))
ggsave(file.path(ROOT, "figures/F7c_candidate_cues_ITC_BLAglut.png"), p, width = 14, height = 5.5, dpi = 220, bg = "white")
ggsave(file.path(ROOT, "figures/F7c_candidate_cues_ITC_BLAglut.pdf"), p, width = 14, height = 5.5)
print(as.data.frame(cand %>% select(cell_type, gene, partner_on_axon, pct_P10, pct_minus_other, enrich_log2, trend)))
