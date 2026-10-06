# E18 / P0 / P10 / P21 developmental update — 6 October 2026

This is the new exploratory analysis. The previous P10/P21 pipeline and data-v1
remain available; their results should not be mixed with v2 thresholds/labels.
Read [METHODS.md](METHODS.md) before interpreting cell fractions or wiring hypotheses.

## What is included?

- R input audit, shared-gene QC, per-library putative doublet calls, conservative
  marker/reference-supported labels, raw/source UMAPs and independent diagnostics.
- Library-level candidate-gene and TF expression; old/new P10 shown separately;
  P10/P21 descriptive CPM ratios with mito and doublet sensitivity.
- Figures as PNG/PDF, numerical CSV/CSV.GZ tables, source JSON and a new
  18-slide collaborator PowerPoint in `results/`.
- A processed Seurat object and the original added counts-only input, published
  as public [data-v2 release downloads](https://github.com/limserenahansol/BLA_scRNA_P10_P21_Mef2c/releases/tag/data-v2-2026-10-06).
  Large data are excluded from git. See [SEURAT_DOWNLOADS.md](SEURAT_DOWNLOADS.md)
  for file descriptions, loading examples and checksum verification.

The new input contains 156,727 barcodes/32,285 genes, with E18_1/E18_2,
P0_1/P0_2/P0_3 and P10_1/P10_2. Add the old P10s1/P21s1/P21s2 libraries:
173,431 input barcodes → 90,605 QC-pass → 76,446 putative singlets; 25,029 shared
genes. These are **library counts, not confirmed animal replicate counts**.
There is no labelled KO/WT comparison or connectivity measurement.

## Start with results, without rerunning

1. Open `results/BLA_development_E18_P0_P10_P21_2026-10-06_final.pptx` (or its companion PDF).
2. Read `results/RESULTS_SUMMARY.md` and the QC/annotation/source diagnostics.
3. For a candidate gene, filter `candidate_expression_per_library.csv` on `gene`,
   `population`, `sample`, and `selection`. Never pool individual cells into
   independent replicates. `annotation_anchor=TRUE` flags genes used in labeling.
4. `candidate_P10_P21_descriptive_effects.csv` contains true raw-pseudobulk CPM
   ratios and every library-pair direction/range. `candidate_source_QC_sensitivity.csv`
   indicates consistency, **not statistical significance or causation**.
5. `TF_expression_per_library.csv.gz` contains the supplied mouse TF screen.
   Gray heatmap entries indicate missing/insufficient-cell groups, not zero expression.

## Data and required reference resources

Download the previous source and annotated objects from
[data-v1](https://github.com/limserenahansol/BLA_scRNA_P10_P21_Mef2c/releases/tag/data-v1).
Download the new source `hansol2_combined_seurat.rds` and processed object
`BLA_E18_P0_P10_P21_seurat_v2.rds` from
[data-v2-2026-10-06](https://github.com/limserenahansol/BLA_scRNA_P10_P21_Mef2c/releases/tag/data-v2-2026-10-06).
Resource bundle contains the adult Allen centroid tables reused from v1 and
the supplied mouse TF list. Unzip its files directly into `resources/`.
Do not distribute unpublished data by making a private repository public.

```r
library(Seurat)
obj <- readRDS("BLA_E18_P0_P10_P21_seurat_v2.rds")
table(obj$sample, obj$broad_class)
table(obj$age, obj$population)
DimPlot(obj, reduction = "umap_source", group.by = "broad_class", raster = TRUE)
DimPlot(obj, reduction = "umap_raw", group.by = "source", raster = TRUE)
FeaturePlot(obj, reduction = "umap_source", features = "Mef2c", raster = TRUE)
counts <- LayerData(obj, assay = "RNA", layer = "counts")
# Annotation uncertainty, library/source and doublet calls are in obj[[]].
# Gene expression/DE must use RNA, not Harmony coordinates.
```

The original new object retains 7,256 additional genes; inspect those there,
not as old/new fold changes. RDS files may be gzip-compressed externally or
uncompressed; R `readRDS()` recognizes both. Keep original input files unchanged.

## Rerun the analysis

Use R >=4.5 and Seurat v5. Exact installed versions are in `results/validation/`.
Required R packages: Seurat, Matrix, harmony, uwot, RANN, digest, scDblFinder,
SingleCellExperiment, BiocParallel, future. Install scDblFinder/SCE/BiocParallel
through Bioconductor; do not silently upgrade an existing analysis environment.
Python plotting requires numpy, pandas, matplotlib. No Cell Ranger is needed
to analyze this RDS; it does not replace Cell Ranger/empty-droplet validation.
The existing raw-data Cell Ranger pipeline still uses **mm10-2020-A**, as requested.
The new RDS does not store its reference build: shared gene symbols alone do not
prove genome/annotation compatibility, so obtain the processing metadata.

Run from the repository root. Example in R:

```r
Sys.setenv(
  BLA_NEW_RDS = "D:/shared/hansol2_combined_seurat.rds",
  BLA_OLD_RDS = "D:/shared/seu_harm_qc.rds",
  BLA_PRIOR_RDS = "D:/shared/B06_P10_P21_annotated.rds",
  BLA_V2_ROOT = "analysis_development_v2",
  BLA_PYTHON = "python"
)
source("analysis_development_v2/scripts/run_pipeline.R")
```

Or set the same environment variables in your terminal and run:

```text
Rscript analysis_development_v2/scripts/run_pipeline.R
```

To restart at step 04 after valid checkpoints exist:

```text
Rscript analysis_development_v2/scripts/run_pipeline.R 4
```

Each step runs in a fresh R process and stops on failure. Logs/checkpoints are
separate from originals. Checkpoints use uncompressed RDS after a compressed
checkpoint failed CRC/reload validation. Budget tens of GB disk and >=64 GB RAM
for this large dataset; peak requirements depend on platform. Per-library
doublet results can resume from cached files. Never reuse caches after changing
input data/QC/parameters: start in a fresh output directory and save provenance.

The checks use sample-level pseudobulk, excluded markers, a held-out gene split,
out-of-library consistency, before/after embedding diagnostics and selection
sensitivity. They deliberately retain failed baselines and ambiguous outcomes.

To prepare local collaborator data after step 05: run `06_export_data.R`, then
`package_outputs.py`, then `07_validate_release.R`. These commands do not upload
anything. The prepared `release/private_data/` folder includes the slim Seurat
object, original new source, all-QC counts, pseudobulk and per-cell audit tables.
Only the two Seurat files listed above have been published from that folder.
The remaining local audit/checkpoint exports are not uploaded by these commands.

## Rebuild the PowerPoint

`scripts/build_presentation.mjs` uses the bundled Codex `@oai/artifact-tool`
runtime, with editable tables/charts and embedded scientific figures. Set
`BLA_PRESENTATION_SKILL` and `BLA_RUNTIME_PYTHON` if your runtime paths differ.
The finalizer verifies package/layout/native chart data and exports slide previews.
For a repeat build, set `BLA_PPTX_OUTPUT` to a new filename; validated decks and
their receipts are intentionally not overwritten.
This special presentation runtime is **not needed** for the R/Python analysis.
The delivered PPTX can be edited in PowerPoint without that runtime. Dense UMAP
and heatmap images have their full-resolution PDF/PNG and numerical sources.

## Interpretation boundaries

P0 low-count retention, very large new P10 barcode sets, high P10_2 putative
doublet calls, poor atlas-only annotation agreement, and age/source confounding
prevent definitive tissue composition or developmental differential-expression
claims. Request Cell Ranger web summaries/unfiltered matrices and animal/sex/
genotype/reference metadata before confirmatory analysis. Adult cortical atlas
context does not identify the actual developing cortical projection neurons.
Expression does not establish axon pruning, cue response, receptor complexes,
TF binding, direct MEF2 targets or the cause of a KO phenotype.
