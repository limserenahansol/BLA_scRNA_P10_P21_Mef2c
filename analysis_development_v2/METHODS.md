# Developmental update: methods and limitations

## Inputs and design

The new `hansol2_combined_seurat.rds` is a counts-only Seurat v5 object:
32,285 genes, 156,727 barcodes, seven libraries. It contains E18_1/E18_2,
P0_1/P0_2/P0_3 and P10_1/P10_2. It does **not** contain a labelled Mef2c KO/WT
comparison, animal IDs, spatial coordinates, Cell Ranger web summaries, FASTQs,
or a documented reference build. A filename is not proof of genotype or anatomy.

From the previous object, retain P10s1/P21s1/P21s2 (16,704 cells before new QC).
The old P56 libraries are not added to the developmental contrast. Use the
25,029 common gene symbols. The 7,256 new-only genes are present in the original
new input but are **not** compared with fictitious zero counts in the old data.
Distinct capture libraries are not proven independent animals. All summaries
are library-level descriptive observations; no biological p-values are reported.

## Integrity, QC and putative doublets

Counts must be finite, nonnegative integers with unique cell/gene names and
exactly agree with stored per-cell UMI/gene metadata. Recalculate QC on the
shared genes. Main QC: >=500 genes, >=1,000 UMIs, <=20% mitochondrial UMIs.
Also report loose >=300 genes/>=500 UMIs/<=25% mito retention, and repeat
expression at <=15% mito. The loose population is **not** reannotated; it is a
retention diagnostic, not a validated alternative population census.

The initial source objects' complete gzip streams pass an independent Python
CRC check. A compressed generated checkpoint failed both an independent CRC
check and Seurat reload. It was regenerated with `compress=FALSE`; important
checkpoints are reloaded and metadata/count hashes checked. This is an output
I/O repair, not a change to original counts. The precise platform cause is not
established. Input R decompression warnings are retained in the audit logs.

Same-barcode suffixes across each old/new library pair are compared using count
cosines and randomized pairs. No pair has highly concordant (>0.97 cosine)
profiles suggestive of resequencing. This does **not** establish animal independence.

Run scDblFinder per capture/library, nfeatures=2,000, dims=20, seed 1062026 +
sample position. Because loading is unknown, dbr=0.08, dbr.sd=1 permits
classifier-driven thresholds; 8% is not a forced removal rate. Keep every score
and call. Primary maps use putative singlets, and expression is repeated including
all QC-pass cells. P10_2 has 32.1% classifier-labelled doublets: a substantial
selection warning, not a measured true doublet rate. A same-seed independent
P10_2 rerun yields identical calls/scores. Warnings are xgboost argument-deprecation
warnings (captured in `doublet_warning_audit.csv`), not hidden numerical errors.

This RDS alone does not provide the unfiltered droplet/background data needed
to validate empty-droplet calling or ambient-RNA correction. Low-count P0 and
the very large new P10 barcode collections need Cell Ranger reports and raw
matrices before a definitive cell census. Passing simple QC is not proof of a cell.

## Annotation and its failed baseline

LogNormalize RNA to 10,000 UMIs. Adult Allen amygdala subclass centroids are a
baseline and context, **not** ground truth for developmental cells. Randomly split
2,914 reference genes into 1,457 training/1,457 held-out genes. Exclude independent
curated broad markers and the hypothesis genes from both splits. Preserve
training and held-out predictions, correlations, margins and agreement.

Atlas-only per-cell correlation gave poor agreement with independent markers.
This failed baseline is reported, not described as successful label transfer.
Add previous P10/P21 fine-label centroids as a second reference, with the same
excluded genes and separate gene splits. These inherit v1 label uncertainty.
Leave one old library out to check recovery of inherited labels; this is an
out-of-library consistency check, **not** independent biological truth.

Broad assignments require >=2 marker genes and either agreement with both atlas
gene splits, agreement with both previous-reference broad gene splits, or a
marker mean log expression >=0.3 and margin >=0.15. A transmitter-supported
neuronal assignment additionally permits >=2 pan-neuronal genes, >=2 UMIs of
Slc17a7/Slc17a6 (glut) or Gad1/Gad2/Slc32a1 (GABA), a >2-fold excess over the
opposite transmitter count, and matching previous-reference broad gene splits.
Immature neurons require >=2 pan-neuronal and >=2 Dcx/Neurod1/Sox11/Igfbpl1 genes.
All other cells remain ambiguous. Inconsistent cells are not forced into a type.
Marker-only labels remain provisional; developmental fine identities are tentative.

An ITC-like candidate requires GABA identity, >=2 Foxp2/Tshz1/Pbx3/Meis2 anchors,
and agreement of previous ITC-like predictions or adult STR-PAL Chst9 predictions.
No spatial ITC identity is established. Rspo2-positive glut cells are a nested
subset of glutamatergic cells selected on observed Rspo2; its expression in this
subset is an annotation anchor, not an independent enrichment discovery. These
nested subsets are never added to broad composition denominators.

## Embedding: useful visualization, limited correction

3,000 variable genes, excluding mito/ribosomal and candidate genes; scale only
those features; PCA=30. Harmony uses `source` (old/new), theta=1, lambda=1,
max_iter=20. Do not regress out age or every library. UMAP: cosine, neighbors=30,
min_dist=0.3, epochs=200, seed=1062026, threads=4, one SGD thread. Preserve raw
and source-corrected UMAPs; graph clusters at resolution=0.6 are exploratory.

Compare broad-class neighbor retention before/after and old/new neighbor mixing
in a balanced P10-only subset. Use the corrected map only if retention falls by
no more than five percentage points. Meeting this gate does **not** prove batch
correction sufficient: residual mixing/poor label purity are separately reported.
RNA counts and normalized RNA must remain hash-identical after correction.
Neither UMAP separation nor overlap establishes maturation, connectivity or cell types.

## Expression, TFs, and sensitivity

Sum uncorrected raw RNA counts within each library/population and compute CPM
against all shared-gene UMIs. Require >=30 cells for displayed developmental
effects. Detection percentage is also reported. Standardized detection is the
analytical hypergeometric expectation at exactly 1,000 UMIs, computed only for
cells with >=1,000 shared-gene UMIs; save its eligible denominator. This controls
one depth effect but not gene-specific capture, reference, dissociation or ambient RNA.

Average library CPM with equal library weights. P10/P21 effects are
`log2((mean_P10_CPM + 0.5)/(mean_P21_CPM + 0.5))`, separately for old/new P10.
Report every P10/P21 library-pair direction and the effect range. No cell-level
p-values or replicate-independent confidence intervals are claimed. Stronger
descriptive candidates retain direction across sources, mito thresholds,
doublet inclusion, and all library pairs. This is a sensitivity criterion, not
a causal or statistically validated target ranking.

P21 is only in the old source. Old/new mito medians and depth differ sharply.
Any new-P10 versus P21 difference is partly confounded with source. Missing
young types, marker dropout, low-cell groups and selection bias remain visible.
TF expression is screened using the supplied 1,321 mouse TF list; expression
does not demonstrate TF activity, binding, direct MEF2 targets or KO effects.

## Corrections to interpretation of v1

- Mean-log-expression differences divided by log(2) are not gene fold changes.
  Use the raw-count library-level CPM ratios defined above.
- Removal of whole marker-mixed clusters does not prove those cells were doublets.
  V2 retains all source cells through QC and stores per-library putative calls.
- MEF2-dependent synapse regulation in prior studies is not direct evidence for
  cortical-axon pruning in BLA. MEF2C effects can be circuit/context dependent.
- Nrp1/Plxnd1 transcript detection or ratios do not establish receptor complexes
  or an attraction/repulsion switch; dropout does not establish protein absence.
- Kirrel3 expression suggests a candidate recognition mechanism, not proven
  cortex-to-ITC targeting. BLA samples do not reveal the projecting cortical
  neurons' genotype, axon abundance, synapses or functional receptor activity.

## Sources

- [Allen whole adult mouse-brain atlas, Yao et al., 2023](https://www.nature.com/articles/s41586-023-06812-z).
- [scDblFinder official vignette](https://www.bioconductor.org/packages/release/bioc/vignettes/scDblFinder/inst/doc/scDblFinder.html).
- [Seurat reference mapping documentation](https://satijalab.org/seurat/articles/multimodal_reference_mapping).
- [R serialization documentation](https://stat.ethz.ch/R-manual/R-devel/library/base/help/saveRDS.html).
- [Flavell et al., 2006: MEF2 and activity-dependent synapse regulation](https://pubmed.ncbi.nlm.nih.gov/16484497/).
- [Tsai et al., 2012: Pcdh10 and MEF2-dependent synapse elimination in experimental neurons](https://pubmed.ncbi.nlm.nih.gov/23260144/); not a BLA axon-pruning measurement.
- [Harrington et al., 2016: circuit/context-dependent MEF2C phenotypes](https://pmc.ncbi.nlm.nih.gov/articles/PMC5094851/).
- [Chauvet et al., 2007: neuropilin/PlexinD1 context in Sema3E responses](https://doi.org/10.1016/j.neuron.2007.10.019).
- [Martin et al., 2015: Kirrel3 and target-specific hippocampal synapses](https://elifesciences.org/articles/09395).

Version-specific session information is copied into `results/validation/`.
Original inputs and data-v1 are preserved.
