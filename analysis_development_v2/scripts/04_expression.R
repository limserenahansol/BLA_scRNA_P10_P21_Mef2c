suppressPackageStartupMessages({library(Seurat); library(Matrix)})
source("analysis_development_v2/scripts/config.R")
logmsg("Sample-level expression and selection sensitivity")
s <- readRDS(datapath("03_all_qc_annotated.rds"))
C <- LayerData(s, layer = "counts"); X <- LayerData(s, layer = "data"); md <- s[[]]
stopifnot(identical(rownames(md), colnames(C)), identical(rownames(C), rownames(X)))
tf <- read.csv(file.path(ROOT, "resources", "mouse_TF_list_1321.csv"))
marker <- read.csv(tabpath("annotation_marker_sets.csv"))
extra <- c("Foxp2", "Tshz1", "Pbx3", "Meis2", "Rspo2", "Etv1", "Satb1", "Satb2", "Tbr1", "Lhx6", "Bcl11b",
  "Nr2f1", "Nr2f2", "Neurod1", "Neurod2", "Sox11", "Dcx", "Igfbpl1", "Snap25", "Syt1", "Elavl3", "Stmn2")
genes <- intersect(unique(c(candidate_genes, tf$gene, marker$gene, extra)), rownames(C))
save_table(data.frame(gene = unique(c(candidate_genes, tf$gene, marker$gene, extra)),
  measured_shared = unique(c(candidate_genes, tf$gene, marker$gene, extra)) %in% rownames(C)), "expression_gene_coverage.csv")
selected <- C[genes, , drop = FALSE]; detected <- selected > 0

# Hypergeometric expectation, not stochastic subsampling. Cells with fewer than
# 1,000 shared-gene UMIs have no standardized-detection estimate.
N <- Matrix::colSums(C)
D <- selected
nc <- rep(N, diff(D@p)); kc <- D@x
D@x <- -expm1(lchoose(pmax(nc - kc, 0), 1000) - lchoose(nc, 1000))
D@x[nc - kc < 1000 & nc >= 1000] <- 1
D@x[nc < 1000] <- 0
stopifnot(all(is.finite(D@x)), all(D@x >= 0 & D@x <= 1), max(abs(Matrix::colSums(C) - md$nCount_shared)) == 0)
selections <- list(primary_singlets = md$doublet_call == "singlet",
  strict_mt15_singlets = md$doublet_call == "singlet" & md$percent_mt <= 15,
  include_putative_doublets = rep(TRUE, nrow(md)))

summaries <- list(); group_info <- list(); bulk <- list(); composition <- list()
for (sel in names(selections)) {
  use <- selections[[sel]]
  for (id in unique(md$sample)) {
    ii <- which(use & md$sample == id)
    if (!length(ii)) next
    z <- as.data.frame(table(factor(md$broad_class[ii], levels = sort(unique(md$broad_class)))))
    names(z) <- c("population", "n_cells")
    z$sample <- id; z$age <- md$age[ii[1]]; z$source <- md$source[ii[1]]
    z$selection <- sel; z$denominator <- length(ii); z$percent <- 100 * z$n_cells / length(ii)
    composition[[length(composition) + 1L]] <- z
    groups <- split(ii, md$broad_class[ii])
    for (p in c("ITC-like GABA (provisional)", "Rspo2-positive glut")) {
      jj <- ii[md$population[ii] == p]
      if (length(jj)) groups[[p]] <- jj
    }
    for (pop in names(groups)) {
      ix <- groups[[pop]]; ng <- length(ix); nstd <- sum(N[ix] >= 1000)
      total <- sum(N[ix]); sums <- Matrix::rowSums(C[, ix, drop = FALSE])
      cpm <- sums / total * 1e6
      meta <- data.frame(selection = sel, sample = id, age = md$age[ix[1]], source = md$source[ix[1]], population = pop,
        n_cells = ng, n_depth1000_eligible = nstd, total_umi = total, adequate_cells30 = ng >= 30,
        grouping_note = if (pop == "Rspo2-positive glut") "Selected on Rspo2 detection; Rspo2 expression is an anchor, not an independent discovery" else
          if (grepl("ITC-like", pop)) "Marker/reference-supported candidate; spatial ITC identity unverified" else "Broad class")
      key <- paste(sel, id, pop, sep = "|")
      if (sel == "primary_singlets") { bulk[[key]] <- sums; group_info[[key]] <- meta }
      zg <- data.frame(gene = genes, counts_sum = sums[genes], cpm = cpm[genes], log2_cpm05 = log2(cpm[genes] + .5),
        detection_pct = 100 * Matrix::rowMeans(detected[, ix, drop = FALSE]),
        detection_depth1000_pct = if (nstd) 100 * Matrix::rowSums(D[, ix, drop = FALSE]) / nstd else NA_real_,
        mean_lognorm = Matrix::rowMeans(X[genes, ix, drop = FALSE]),
        annotation_anchor = genes %in% unique(c(marker$gene, "Foxp2", "Tshz1", "Pbx3", "Meis2", "Rspo2")))
      summaries[[length(summaries) + 1L]] <- cbind(meta[rep(1, nrow(zg)), ], zg, row.names = NULL)
    }
  }
  logmsg("Finished selection:", sel)
}
e <- do.call(rbind, summaries); co <- do.call(rbind, composition)
eg <- gzfile(tabpath("sample_population_expression.csv.gz"), "wt")
write.csv(e, eg, row.names = FALSE, na = ""); close(eg)
save_table(e[e$gene %in% candidate_genes, ], "candidate_expression_per_library.csv")
tg <- gzfile(tabpath("TF_expression_per_library.csv.gz"), "wt")
write.csv(e[e$gene %in% tf$gene, ], tg, row.names = FALSE, na = ""); close(tg)
save_table(co, "population_composition_per_library.csv")
gi <- do.call(rbind, group_info); save_table(gi, "pseudobulk_group_metadata.csv")
b <- do.call(cbind, bulk)
stopifnot(identical(colnames(b), rownames(gi)), all(b >= 0), all(colSums(b) == gi$total_umi))
saveRDS(list(counts = b, groups = gi, unit = "library; animal independence unknown", genes = rownames(C)), datapath("sample_population_pseudobulk.rds"))
con <- gzfile(tabpath("all_gene_pseudobulk_counts.csv.gz"), "wt")
write.csv(data.frame(gene = rownames(b), b, check.names = FALSE), con, row.names = FALSE); close(con)

# P10 overlaps sources; never hide this by merging the libraries into one dot.
# P21 has only two old-source libraries. All effects are descriptive.
effects <- list()
for (sel in names(selections)) for (pop in unique(e$population)) {
  zz <- e[e$selection == sel & e$population == pop & e$adequate_cells30 & e$gene %in% candidate_genes, ]
  for (g in unique(zz$gene)) {
    z <- zz[zz$gene == g, ]; p21 <- z[z$age == "P21", ]
    if (nrow(p21) < 2) next
    for (src in c("old_2023", "new_2026")) {
      p10 <- z[z$age == "P10" & z$source == src, ]
      if (!nrow(p10)) next
      delta <- outer(log2(p10$cpm + .5), log2(p21$cpm + .5), "-")
      effects[[length(effects) + 1L]] <- data.frame(selection = sel, population = pop, gene = g, P10_source = src,
        n_P10_libraries = nrow(p10), n_P21_libraries = nrow(p21),
        log2_CPM_ratio_P10_over_P21 = log2((mean(p10$cpm) + .5) / (mean(p21$cpm) + .5)),
        pairwise_log2_ratio_min = min(delta), pairwise_log2_ratio_max = max(delta),
        pairwise_direction_consistent = all(delta > 0) || all(delta < 0),
        depth1000_detection_difference_pp = mean(p10$detection_depth1000_pct) - mean(p21$detection_depth1000_pct),
        P10_CPM_mean = mean(p10$cpm), P21_CPM_mean = mean(p21$cpm),
        annotation_anchor = any(z$annotation_anchor), interpretation = "descriptive; no KO contrast; animal n unknown")
    }
  }
}
fx <- do.call(rbind, effects); save_table(fx, "candidate_P10_P21_descriptive_effects.csv")
robust <- list()
for (pop in unique(fx$population)) for (g in candidate_genes) {
  z <- fx[fx$population == pop & fx$gene == g, ]; if (!nrow(z)) next
  need <- expand.grid(selection = names(selections), P10_source = c("old_2023", "new_2026"))
  available <- nrow(z) == nrow(need) && all(paste(need$selection, need$P10_source) %in% paste(z$selection, z$P10_source))
  same <- available && (all(z$log2_CPM_ratio_P10_over_P21 > 0) || all(z$log2_CPM_ratio_P10_over_P21 < 0))
  robust[[length(robust) + 1L]] <- data.frame(population = pop, gene = g, all_sources_and_selections_available = available,
    same_effect_direction_across_sources_and_selections = same,
    all_P10_P21_library_pair_directions_agree = available && all(z$pairwise_direction_consistent) && same,
    min_log2_ratio = min(z$log2_CPM_ratio_P10_over_P21), max_log2_ratio = max(z$log2_CPM_ratio_P10_over_P21),
    annotation_anchor = any(z$annotation_anchor))
}
save_table(do.call(rbind, robust), "candidate_source_QC_sensitivity.csv")
save_table(data.frame(check = c("pseudobulk_group_UMI_matches", "depth1000_probability_bounds", "RNA_counts_match_metadata"), passed = TRUE), "expression_validation.csv")
writeLines(capture.output(sessionInfo()), file.path(ROOT, "logs", "sessionInfo_expression.txt"))
logmsg("Expression complete:", nrow(e), "gene/sample/population/selection summaries")
