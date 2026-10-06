# Developmental update, 2026-10-06

Objective: add E18/P0 and two P10 libraries to the previous P10/P21 amygdala
analysis, describe populations and candidate MEF2/pruning/guidance genes, and
deliver a reproducible Seurat object, results and collaborator presentation.

## Design and units

The experimental summaries use libraries/samples, not individual cells as
replicates. Animal IDs, sex, genotype and Cell Ranger reference/version are not
stored in the new object. Distinct libraries are not automatically distinct
animals. These data do not contain a labelled Mef2c KO/WT comparison.
P21 occurs only in the old dataset. P10 overlaps the datasets and permits an
explicit old/new comparison. Age and source remain partly confounded.

## Approaches considered

1. Merge and analyse uncorrected RNA: transparent, preserves age, but batch and
   depth can separate cells. This is the expression baseline and source of all
   expression summaries.
2. Harmony across every library: convenient for visualisation, but library and
   age are confounded. Rejected as the primary developmental analysis.
3. Shared-gene RNA, independent marker/atlas annotation, and a source-level
   Harmony embedding: chosen. The correction is for a secondary map only.
   Both uncorrected and corrected maps and their P10-only diagnostics are saved.

## Predeclared parameters and validation

- Common gene symbols only for cross-dataset analysis. Retain the original
  new object, including its extra genes, as a separate release input.
- Main QC: >=500 detected shared genes, >=1,000 shared-gene UMIs, <=20% mito.
  Sensitivity: >=300 genes, >=500 UMIs, <=25% mito; and strict <=15% mito.
  QC metrics are recalculated from counts. No arbitrary upper-gene cutoff.
- Test same-barcode old/new libraries for highly similar count vectors.
  Record overlap and cosine similarity. Flag possible resequenced libraries;
  do not infer biological independence from barcode naming.
- scDblFinder per capture/library, seed 1062026. Unknown loading prevents a
  reliable prior rate. Use dbr=0.08 with dbr.sd=1 (wide prior, classifier-driven
  threshold), retain scores/calls and report include/exclude sensitivity.
- LogNormalize scale 10,000, 3,000 variable genes, PCA 30 components, UMAP seed
  1062026. Exclude mitochondrial/ribosomal genes and the candidate hypothesis
  genes from the embedding/annotation feature set where applicable.
- Adult Allen subclass and previous-label centroids guide annotation. A fixed
  separate gene split and curated markers test agreement. Adult labels in E18
  and P0 are explicitly tentative. An ambiguous label is a permitted outcome.
- Evaluate counts/metadata agreement and no duplicate IDs; QC before/after;
  marker-versus-reference agreement and held-out genes; old/new P10 cell-type
  concordance; RNA versus Harmony neighbourhood diagnostics; and sensitivity
  of candidate effects to mito and doublet selection.
- Expression uses sample-by-population raw pseudobulk counts (CPM), with
  detection percentages and analytical expected detection at 1,000 UMIs.
  Report sample values, ranges and effect directions. No cell-level p-values.
- P10/P21 changes are log2 ratios of pseudobulk CPM (+0.5 CPM), not differences
  of mean log-normalised values mislabeled as fold changes. Show old P10 and
  new P10 separately and whether each agrees with both P21 samples.
- Neither gene expression nor adult cortex co-expression demonstrates pruning,
  attraction/repulsion, axonal connectivity, direct MEF2 targets or KO effects.

Original input files and the previous release remain unchanged. Consequential
differences from v1 are documented in METHODS.md and the new release notes.
