# Final cell-type labels for P10 + P21 BLA cells: res-1 clusters for glia,
# glutamatergic (B04) and GABAergic (B05) sub-clusters for neurons. Doublets
# (two-lineage co-expression) and low-quality clusters are removed and logged.
suppressPackageStartupMessages({ library(Seurat); library(dplyr) })
ROOT <- Sys.getenv("BLA_ROOT", "K:/scRNA_BLA_phd/P10_P21_amygdala_2026-10")
seu  <- readRDS(file.path(ROOT, "B02_P10_P21_mapped.rds"))
glut <- readRDS(file.path(ROOT, "B04_glut.rds"))
gaba <- readRDS(file.path(ROOT, "B05_gaba.rds"))

nn_map <- c(`0` = "Astrocyte", `13` = "Astrocyte", `3` = "Microglia", `4` = "OPC", `29` = "OPC",
            `28` = "OPC cycling", `19` = "COP/NFOL", `18` = "Oligodendrocyte", `20` = "Endothelial",
            `27` = "Pericyte/mural", `25` = "Ependymal", `30` = "Cycling progenitor")
glut_map <- c(`3` = "BLA Glut Rspo2/Etv1", `0` = "BLA Glut Tshz2/Satb1", `4` = "BLA Glut Cdh8/Tshz3",
              `5` = "BLA Glut Fgf10/Nrp1", `6` = "BLA Glut Otof/Trhr", `7` = "BLA Glut Vgll3/Prr16",
              `2` = "BMA Glut Slc17a6/Zfp804b", `10` = "PA/BMAp Glut Esr1/Reln",
              `13` = "Immature Glut Igfbpl1/Epha3", `14` = "CLA/EPd-like Glut Satb2/Nr4a2")
gaba_map <- c(`0` = "ITC Foxp2/Tshz1 (a)", `2` = "ITC Foxp2/Kcnh5 (b)",
              `1` = "CeA/MeA GABA Meis2/Tshz2", `7` = "CeA GABA Prkcd/Calcrl",
              `4` = "D2 SPN-like Adora2a", `9` = "D1 SPN-like Ebf1/Tac1",
              `3` = "IN Sst", `5` = "IN Maf/Calb1 (Pvalb-lineage)", `11` = "IN Moxd1/Gna14 (MGE)",
              `15` = "IN Pvalb chandelier", `8` = "IN Vip", `12` = "IN Reln/Npas1 (CGE)",
              `13` = "IN Cck/Cnr1 (CGE)", `10` = "IN Lamp5", `17` = "Unresolved Nts/Pappa2")

ct <- setNames(rep(NA_character_, ncol(seu)), colnames(seu))
why <- setNames(rep("", ncol(seu)), colnames(seu))
x <- as.character(seu$cl); ct[x %in% names(nn_map)] <- nn_map[x[x %in% names(nn_map)]]
why[x %in% c("10", "21")] <- "low quality (few genes / mito)"
why[x %in% c("24", "31", "32")] <- "doublet (neuron or OPC + microglia/OPC genes)"
g <- as.character(glut$gcl); ok <- g %in% names(glut_map)
ct[colnames(glut)[ok]] <- glut_map[g[ok]]
why[colnames(glut)[!ok]] <- ifelse(g[!ok] == "8", "low quality glut", "doublet (glut + ITC/glia genes)")
g <- as.character(gaba$gcl); ok <- g %in% names(gaba_map)
ct[colnames(gaba)[ok]] <- gaba_map[g[ok]]
why[colnames(gaba)[!ok]] <- "doublet (GABA + astro/OPC/microglia genes)"

removed <- data.frame(cell = names(ct)[is.na(ct)], reason = why[is.na(ct)], sample = seu$sample[is.na(ct)])
write.csv(removed %>% count(reason, sample), file.path(ROOT, "tables/B06_removed_cells.csv"), row.names = FALSE)
seu$cell_type <- ct
seu <- subset(seu, cells = names(ct)[!is.na(ct)])

order_ct <- c(glut_map[c("3","0","4","5","6","7","2","10","13","14")],
              gaba_map[c("0","2","1","7","4","9","3","5","11","15","8","12","13","10","17")],
              unique(nn_map))
seu$cell_type <- factor(seu$cell_type, levels = order_ct)
seu$class <- factor(ifelse(seu$cell_type %in% glut_map, "Glutamatergic",
                    ifelse(seu$cell_type %in% gaba_map, "GABAergic", "Non-neuronal")),
                    levels = c("Glutamatergic", "GABAergic", "Non-neuronal"))
seu$glut_umap1 <- NA; seu$glut_umap2 <- NA
e <- Embeddings(glut, "umap"); keep <- intersect(rownames(e), colnames(seu))
seu$glut_umap1[match(keep, colnames(seu))] <- e[keep, 1]; seu$glut_umap2[match(keep, colnames(seu))] <- e[keep, 2]
saveRDS(seu, file.path(ROOT, "B06_P10_P21_annotated.rds"))
tab <- as.data.frame.matrix(table(seu$cell_type, seu$sample))
write.csv(tab, file.path(ROOT, "tables/B06_cell_type_counts_per_sample.csv"))
print(tab); cat("removed:", nrow(removed), "of", length(ct), "\n")
