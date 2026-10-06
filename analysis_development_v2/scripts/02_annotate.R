suppressPackageStartupMessages({library(Seurat); library(Matrix)})
source("analysis_development_v2/scripts/config.R")
set.seed(SEED)
logmsg("Reading QC checkpoint for annotation")
d <- readRDS(datapath("01_qc_input.rds"))
s <- CreateSeuratObject(d$counts, meta.data = d$metadata, project = "BLA_development_v2")
s <- NormalizeData(s, scale.factor = 10000, verbose = FALSE)
X <- LayerData(s, layer = "data")
rm(d); invisible(gc())

# These marker sets are excluded from atlas training and its held-out gene split.
markers <- list(
  Glutamatergic = c("Slc17a7", "Slc17a6", "Neurod2", "Neurod6", "Tbr1"),
  GABAergic = c("Gad1", "Gad2", "Slc32a1", "Dlx1", "Dlx2"),
  Astrocyte = c("Aldh1l1", "Slc1a3", "Aqp4", "Gja1", "Aldoc"),
  OPC = c("Pdgfra", "Cspg4", "Olig1", "Olig2", "Sox10"),
  Oligodendrocyte = c("Mbp", "Mog", "Plp1", "Mag", "Opalin"),
  Microglia = c("C1qa", "C1qb", "Csf1r", "Cx3cr1", "P2ry12", "Aif1"),
  Endothelial = c("Cldn5", "Pecam1", "Flt1", "Kdr", "Esam"),
  Pericyte = c("Pdgfrb", "Rgs5", "Abcc9", "Kcnj8", "Vtn"),
  Fibroblast = c("Dcn", "Col1a1", "Col1a2", "Lum"),
  Ependymal = c("Foxj1", "Tmem212", "Ccdc153", "Pifo"),
  Choroid = c("Ttr", "Folr1", "Kcnj13"),
  Progenitor = c("Sox2", "Hes1", "Hes5", "Nes", "Fabp7"),
  Erythroid = c("Hba-a1", "Hba-a2", "Hbb-bs", "Hbb-bt", "Alas2"))
markers <- lapply(markers, intersect, rownames(X))
save_table(do.call(rbind, lapply(names(markers), function(k) data.frame(broad_class = k, gene = markers[[k]]))), "annotation_marker_sets.csv")
scores <- sapply(markers, function(g) Matrix::colMeans(X[g, , drop = FALSE]))
detected <- sapply(markers, function(g) Matrix::colSums(X[g, , drop = FALSE] > 0))
top <- max.col(scores, ties.method = "first")
top_score <- scores[cbind(seq_len(nrow(scores)), top)]
second <- apply(scores, 1, function(v) sort(v, decreasing = TRUE)[2])
s$marker_class <- colnames(scores)[top]
s$marker_score <- top_score
s$marker_margin <- top_score - second
s$marker_n_detected <- detected[cbind(seq_len(nrow(scores)), top)]

ref <- read.csv(file.path(ROOT, "resources", "ref_subclass_mean_log2_markers.csv"), row.names = 1, check.names = FALSE)
excluded <- unique(c(candidate_genes, unlist(markers), "Foxp2", "Tshz1", "Pbx3", "Meis2", "Rspo2", "Etv1"))
genes <- setdiff(intersect(rownames(ref), rownames(X)), excluded)
genes <- genes[!grepl("^mt-|^Rpl|^Rps", genes)]
genes <- sample(sort(genes))
ga <- genes[seq_along(genes) %% 2 == 1]; gb <- setdiff(genes, ga)
save_table(data.frame(gene = genes, split = ifelse(genes %in% ga, "training", "held_out")), "reference_gene_split.csv")

cor_map <- function(reference, query, genes, block = 2500L) {
  R <- as.matrix(reference[genes, , drop = FALSE])
  R <- sweep(R, 2, colMeans(R), "-")
  R <- sweep(R, 2, sqrt(colSums(R^2)), "/")
  out <- matrix(0, ncol(query), ncol(R), dimnames = list(colnames(query), colnames(R)))
  for (start in seq.int(1L, ncol(query), by = block)) {
    ix <- start:min(start + block - 1L, ncol(query))
    Q <- as.matrix(query[genes, ix, drop = FALSE])
    Q <- sweep(Q, 2, colMeans(Q), "-")
    Q <- sweep(Q, 2, pmax(sqrt(colSums(Q^2)), 1e-12), "/")
    out[ix, ] <- crossprod(Q, R)
  }
  out
}
best <- function(corr) {
  k <- max.col(corr, ties.method = "first")
  data.frame(label = colnames(corr)[k], r = corr[cbind(seq_len(nrow(corr)), k)],
    margin = corr[cbind(seq_len(nrow(corr)), k)] - apply(corr, 1, function(v) sort(v, decreasing = TRUE)[2]))
}
broad_of_atlas <- function(v) {
  out <- rep("Other", length(v))
  out[grepl("Glut", v)] <- "Glutamatergic"; out[grepl("Gaba", v)] <- "GABAergic"
  out[grepl("IMN", v)] <- "Immature_neuron"
  keys <- c(Astrocyte = "Astro", OPC = "OPC", Oligodendrocyte = "Oligo", Microglia = "Microglia|BAM",
    Endothelial = "Endo", Pericyte = "Peri|SMC", Fibroblast = "VLMC", Ependymal = "Ependymal", Choroid = "Choroid")
  for (k in names(keys)) out[grepl(keys[k], v)] <- k
  out
}
logmsg("Atlas training split:", length(ga), "genes; held-out:", length(gb))
a <- best(cor_map(ref, X, ga)); b <- best(cor_map(ref, X, gb))
s$atlas_subclass <- a$label; s$atlas_r <- a$r; s$atlas_margin <- a$margin
s$atlas_subclass_heldout <- b$label
s$atlas_class <- broad_of_atlas(a$label); s$atlas_class_heldout <- broad_of_atlas(b$label)
s$atlas_split_agree <- s$atlas_class == s$atlas_class_heldout
s$marker_atlas_agree <- s$marker_class == s$atlas_class
strong_marker <- s$marker_n_detected >= 2 & s$marker_score >= .3 & s$marker_margin >= .15
supported <- s$atlas_split_agree & s$marker_atlas_agree & s$marker_n_detected >= 2
s$broad_class <- ifelse(supported | strong_marker, s$marker_class, "Ambiguous")
s$annotation_support <- ifelse(supported, "marker_and_reference", ifelse(strong_marker, "marker_only", "ambiguous"))
pan <- intersect(c("Snap25", "Syt1", "Elavl3", "Stmn2"), rownames(X))
young <- intersect(c("Dcx", "Neurod1", "Sox11", "Igfbpl1"), rownames(X))
s$immature_neuron_score <- Matrix::colMeans(X[young, , drop = FALSE])
im <- s$broad_class == "Ambiguous" & Matrix::colSums(X[pan, , drop = FALSE] > 0) >= 2 &
  Matrix::colSums(X[young, , drop = FALSE] > 0) >= 2
s$broad_class[im] <- "Immature_neuron"; s$annotation_support[im] <- "immature_markers"

# Previous fine labels are predictions, not a claim that an adult identity exists
# at E18/P0. Candidate genes and independent marker genes are excluded here too.
prior <- readRDS(PRIOR_INPUT)
P <- LayerData(prior, layer = "counts")[rownames(X), , drop = FALSE]
P <- log1p(t(t(P) / Matrix::colSums(P)) * 10000)
labels <- as.character(prior$cell_type)
labels[grepl("^ITC", labels)] <- "ITC-like Foxp2/Tshz1"
centroid <- function(mat, label) sapply(split(seq_len(ncol(mat)), label), function(i) Matrix::rowMeans(mat[, i, drop = FALSE]))
cent <- centroid(P, labels)
pa <- best(cor_map(cent, X, ga)); pb <- best(cor_map(cent, X, gb))
s$previous_label_prediction <- pa$label; s$previous_label_r <- pa$r; s$previous_label_margin <- pa$margin
s$previous_label_heldout <- pb$label; s$fine_split_agree <- pa$label == pb$label
save_table(data.frame(gene = rownames(cent), cent), "previous_label_centroids.csv")

# Validation of the atlas-only baseline showed poor agreement, especially in
# young/low-depth cells. Do not tune labels to a UMAP. Use previous P10/P21
# centroids as a second reference, still excluding the independent markers.
previous_broad <- function(v) {
  z <- rep("Other", length(v))
  z[grepl("Glut", v)] <- "Glutamatergic"
  z[grepl("^IN |^IN\\.|ITC|GABA|SPN", v)] <- "GABAergic"
  for (k in c("Astrocyte", "OPC", "Oligodendrocyte", "Microglia", "Endothelial", "Ependymal", "Pericyte"))
    z[grepl(k, v)] <- k
  z[grepl("COP|NFOL", v)] <- "Oligodendrocyte"
  z[grepl("progenitor", v, ignore.case = TRUE)] <- "Progenitor"
  z
}
s$previous_broad_prediction <- previous_broad(pa$label)
s$previous_broad_heldout <- previous_broad(pb$label)
s$previous_broad_split_agree <- s$previous_broad_prediction == s$previous_broad_heldout
s$marker_previous_agree <- s$marker_class == s$previous_broad_prediction
prev_supported <- s$previous_broad_split_agree & s$marker_previous_agree & s$marker_n_detected >= 2
s$broad_class[prev_supported] <- s$marker_class[prev_supported]
s$annotation_support[prev_supported] <- "marker_and_previous_reference"
C <- LayerData(s, layer = "counts")
pan_n <- Matrix::colSums(C[pan, , drop = FALSE] > 0)
glut_tx <- Matrix::colSums(C[intersect(c("Slc17a7", "Slc17a6"), rownames(C)), , drop = FALSE])
gaba_tx <- Matrix::colSums(C[intersect(c("Gad1", "Gad2", "Slc32a1"), rownames(C)), , drop = FALSE])
for (k in c("Glutamatergic", "GABAergic")) {
  tx <- if (k == "Glutamatergic") glut_tx else gaba_tx
  other_tx <- if (k == "Glutamatergic") gaba_tx else glut_tx
  keep <- s$previous_broad_split_agree & s$previous_broad_prediction == k & pan_n >= 2 & tx >= 2 & tx > 2 * other_tx
  s$broad_class[keep] <- k
  s$annotation_support[keep] <- "transmitter_pan_neuron_and_previous_reference"
}
s$glut_transmitter_umi <- glut_tx; s$gaba_transmitter_umi <- gaba_tx; s$pan_neuron_n_detected <- pan_n
rm(C)

s$population <- s$broad_class
itc_anchor <- Matrix::colSums(X[intersect(c("Foxp2", "Tshz1", "Pbx3", "Meis2"), rownames(X)), , drop = FALSE] > 0) >= 2
itc <- s$broad_class == "GABAergic" & itc_anchor &
  ((s$fine_split_agree & grepl("^ITC", s$previous_label_prediction)) |
   (grepl("STR-PAL Chst9", s$atlas_subclass) & grepl("STR-PAL Chst9", s$atlas_subclass_heldout)))
s$population[itc] <- "ITC-like GABA (provisional)"
rs <- s$broad_class == "Glutamatergic" & as.vector(X["Rspo2", ] > 0)
s$population[rs] <- "Rspo2-positive glut"
s$annotation_age_status <- ifelse(s$age %in% c("E18", "P0"), "developmental_identity_tentative", "marker_supported_broad_identity")

validation <- do.call(rbind, lapply(split(seq_len(ncol(s)), s$sample), function(i) data.frame(sample = s$sample[i[1]],
  age = s$age[i[1]], n = length(i), split_broad_agreement = mean(s$atlas_split_agree[i]),
  independent_marker_atlas_agreement = mean(s$marker_atlas_agree[i]), ambiguous_pct = 100 * mean(s$broad_class[i] == "Ambiguous"),
  marker_only_pct = 100 * mean(s$annotation_support[i] == "marker_only"), fine_split_agreement = mean(s$fine_split_agree[i]))))
validation$previous_broad_split_agreement <- sapply(split(seq_len(ncol(s)), s$sample), function(i) mean(s$previous_broad_split_agree[i]))
validation$independent_marker_previous_agreement <- sapply(split(seq_len(ncol(s)), s$sample), function(i) mean(s$marker_previous_agree[i]))
save_table(validation, "annotation_validation_per_sample.csv")
loo <- list()
for (id in c("P10s1", "P21s1", "P21s2")) {
  train <- as.character(prior$sample) != id
  cc <- centroid(P[, train, drop = FALSE], labels[train])
  test <- which(as.character(prior$sample) == id)
  pr <- best(cor_map(cc, P[, test, drop = FALSE], ga))
  loo[[id]] <- data.frame(sample = id, n = length(test), previous_label_agreement = mean(pr$label == labels[test]),
    evaluation = "leave-one-old-library-out; inherited labels, not independent ground truth")
}
save_table(do.call(rbind, loo), "leave_one_old_library_out.csv")
save_table(data.frame(cell_id = colnames(s), sample = s$sample, scores), "cell_marker_scores.csv")
stopifnot(inherits(s, "Seurat"), ncol(s) == nrow(s[[]]))
saveRDS(s, datapath("02_annotated_qc.rds"), compress = FALSE)
check <- readRDS(datapath("02_annotated_qc.rds"))
stopifnot(inherits(check, "Seurat"), identical(colnames(check), colnames(s)), identical(check[[]], s[[]]))
rm(check); invisible(gc())
saveRDS(list(markers = markers, training_genes = ga, heldout_genes = gb, reference = ref), datapath("annotation_reference.rds"))
print(validation, row.names = FALSE)
print(table(s$sample, s$broad_class))
writeLines(capture.output(sessionInfo()), file.path(ROOT, "logs", "sessionInfo_annotation.txt"))
logmsg("Annotation checkpoint complete")
