suppressPackageStartupMessages({library(Seurat); library(Matrix); library(harmony); library(uwot); library(RANN); library(digest)})
source("analysis_development_v2/scripts/config.R")
set.seed(SEED)
s <- readRDS(datapath("02_annotated_qc.rds"))
db <- read.csv(tabpath("doublet_calls.csv"))
ix <- match(colnames(s), db$cell_id); stopifnot(!anyNA(ix))
s$doublet_call <- db$doublet_call[ix]; s$doublet_score <- db$doublet_score[ix]
stopifnot(inherits(s, "Seurat"), nrow(db) == ncol(s))
saveRDS(s, datapath("03_all_qc_annotated.rds"), compress = FALSE)
s <- subset(s, cells = colnames(s)[s$doublet_call == "singlet"])
s$age <- factor(s$age, levels = age_levels)
logmsg("Embedding predicted singlets:", ncol(s))
rna_hash_before <- digest(LayerData(s, layer = "counts"), algo = "xxhash64")
data_hash_before <- digest(LayerData(s, layer = "data"), algo = "xxhash64")
s <- FindVariableFeatures(s, nfeatures = 3000, verbose = FALSE)
vg <- setdiff(VariableFeatures(s), candidate_genes)
vg <- vg[!grepl("^mt-|^Rpl|^Rps", vg)]
VariableFeatures(s) <- vg
save_table(data.frame(gene = vg), "embedding_genes.csv")
s <- ScaleData(s, features = vg, verbose = FALSE)
s <- RunPCA(s, features = vg, npcs = 30, seed.use = SEED, verbose = FALSE)
pc <- Embeddings(s, "pca")
cache <- datapath("BLA_E18_P0_P10_P21_seurat_v2.rds")
reused <- FALSE
if (file.exists(cache)) {
  old_map <- readRDS(cache)
  if (identical(colnames(old_map), colnames(s)) &&
      identical(rna_hash_before, digest(LayerData(old_map, layer = "counts"), algo = "xxhash64")) &&
      identical(rownames(Loadings(old_map, "pca")), vg)) {
    for (nm in c("pca", "harmony_source", "umap_raw", "umap_source")) s[[nm]] <- old_map[[nm]]
    for (nm in names(old_map@graphs)) s@graphs[[nm]] <- old_map@graphs[[nm]]
    s$seurat_clusters <- old_map$seurat_clusters
    reused <- TRUE
    logmsg("Reused identical RNA embeddings; refreshed annotations and independent diagnostics")
  }
  rm(old_map); invisible(gc())
}
if (!reused) {
logmsg("Source-level Harmony; no age or per-library regression")
s <- harmony::RunHarmony(s, group.by.vars = "source", reduction.use = "pca", dims.use = 1:30,
  reduction.save = "harmony_source", theta = 1, lambda = 1, max_iter = 20, ncores = 4, verbose = TRUE)
hc <- Embeddings(s, "harmony_source")
for (name in c("raw", "source")) {
  logmsg("UMAP:", name)
  set.seed(SEED)
  inp <- if (name == "raw") pc else hc
  e <- uwot::umap(inp[, 1:30], n_neighbors = 30, min_dist = .3, metric = "cosine", n_epochs = 200,
    n_threads = 4, n_sgd_threads = 1, verbose = TRUE)
  rownames(e) <- colnames(s); colnames(e) <- paste0("UMAP", 1:2)
  s[[paste0("umap_", name)]] <- CreateDimReducObject(e, key = paste0("UMAP", toupper(name), "_"), assay = "RNA")
}
s <- FindNeighbors(s, reduction = "harmony_source", dims = 1:30, graph.name = c("source_nn", "source_snn"), verbose = FALSE)
s <- FindClusters(s, graph.name = "source_snn", resolution = .6, random.seed = SEED, verbose = FALSE)
}

# Independent checks use marker/reference-supported classes, not the graph
# clusters that were optimised. The P10 mixing test balances old/new sources.
set.seed(SEED)
classes <- s$broad_class
valid <- which(classes != "Ambiguous")
eval_ix <- if (length(valid) > 10000) sample(valid, 10000) else valid
diagnostics <- list()
for (nm in c("pca", "harmony_source")) {
  v <- Embeddings(s, nm)[eval_ix, 1:30]
  nn <- RANN::nn2(v, k = 16)$nn.idx[, -1, drop = FALSE]
  retention <- mean(matrix(classes[eval_ix][nn], nrow(nn)) == classes[eval_ix])
  diagnostics[[length(diagnostics) + 1]] <- data.frame(reduction = nm, test = "marker_class_neighbour_retention", n = length(eval_ix), value = retention)
}
p10 <- which(s$age == "P10" & classes != "Ambiguous")
balanced <- unlist(lapply(split(p10, classes[p10]), function(i) {
  a <- i[s$source[i] == "old_2023"]; b <- i[s$source[i] == "new_2026"]
  n <- min(length(a), length(b), 1500L)
  if (n < 50) return(integer())
  c(sample(a, n), sample(b, n))
}))
if (length(balanced) > 100) for (nm in c("pca", "harmony_source")) {
  nn <- RANN::nn2(Embeddings(s, nm)[balanced, 1:30], k = 16)$nn.idx[, -1, drop = FALSE]
  other_source <- mean(matrix(s$source[balanced][nn], nrow(nn)) != s$source[balanced])
  retention <- mean(matrix(classes[balanced][nn], nrow(nn)) == classes[balanced])
  diagnostics[[length(diagnostics) + 1]] <- data.frame(reduction = nm, test = "P10_balanced_other_source_neighbours", n = length(balanced), value = other_source)
  diagnostics[[length(diagnostics) + 1]] <- data.frame(reduction = nm, test = "P10_balanced_class_retention", n = length(balanced), value = retention)
}
diagnostics <- do.call(rbind, diagnostics)
save_table(diagnostics, "embedding_validation.csv")
baseline <- diagnostics$value[diagnostics$reduction == "pca" & diagnostics$test == "marker_class_neighbour_retention"]
corrected <- diagnostics$value[diagnostics$reduction == "harmony_source" & diagnostics$test == "marker_class_neighbour_retention"]
s@misc$primary_map <- if (corrected < baseline - .05) "umap_raw" else "umap_source"
s@misc$primary_map_reason <- "Source map used only if marker-class neighbour retention falls by no more than 5 percentage points versus PCA. Neither map establishes developmental biology."
stopifnot(identical(rna_hash_before, digest(LayerData(s, layer = "counts"), algo = "xxhash64")),
  identical(data_hash_before, digest(LayerData(s, layer = "data"), algo = "xxhash64")))
save_table(data.frame(check = c("counts_unchanged_by_embedding", "normalised_RNA_unchanged_by_embedding"), passed = TRUE), "embedding_rna_integrity.csv")
coord <- cbind(cell_id = colnames(s), s[[]],
  umap_raw1 = Embeddings(s, "umap_raw")[, 1], umap_raw2 = Embeddings(s, "umap_raw")[, 2],
  umap_source1 = Embeddings(s, "umap_source")[, 1], umap_source2 = Embeddings(s, "umap_source")[, 2])
save_table(coord, "cell_metadata.csv")
saveRDS(s, datapath("BLA_E18_P0_P10_P21_seurat_v2.rds"), compress = FALSE)
check <- readRDS(datapath("BLA_E18_P0_P10_P21_seurat_v2.rds"))
stopifnot(inherits(check, "Seurat"), identical(check[[]], s[[]]),
  identical(rna_hash_before, digest(LayerData(check, layer = "counts"), algo = "xxhash64")),
  identical(data_hash_before, digest(LayerData(check, layer = "data"), algo = "xxhash64")))
rm(check); invisible(gc())
writeLines(capture.output(sessionInfo()), file.path(ROOT, "logs", "sessionInfo_embeddings.txt"))
print(diagnostics, row.names = FALSE)
logmsg("Embedding complete; primary map:", s@misc$primary_map)
