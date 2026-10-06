# Seurat downloads: developmental v2

Public release: [data-v2-2026-10-06](https://github.com/limserenahansol/BLA_scRNA_P10_P21_Mef2c/releases/tag/data-v2-2026-10-06).
These are GitHub Release assets, not files stored in the git history.

## Which file should I use?

- **`BLA_E18_P0_P10_P21_seurat_v2.rds` (696,108,887 bytes):** use this to
  inspect the expanded analysis. Seurat v5, 76,446 putative singlets and
  25,029 shared genes across 10 libraries. Contains RNA counts, normalized RNA,
  per-cell annotations/QC/source metadata and dimensional reductions.
  Stage totals: E18 6,936; P0 15,199; P10 47,757; P21 6,554.
  Flagged doublets are excluded; ambiguous annotations are retained explicitly.
- **`hansol2_combined_seurat.rds` (372,752,891 bytes):** byte-identical original
  added counts-only Seurat object: 156,727 barcodes and 32,285 genes from seven
  E18/P0/P10 libraries. It does not include P21 or the new filtered annotations/
  UMAPs. This is a count-matrix input, **not Illumina BCL/FASTQ data**.
  It preserves 7,256 genes absent from the shared-gene comparative object.
- **`SEURAT_SHA256.txt`:** SHA-256 checksums of both RDS downloads.

The 3 GB uncompressed diagnostic copy and other intermediate exports are not
needed to inspect the main object and are not included in this release.
Older P10/P21/P56 inputs and legacy subcluster objects remain in
[data-v1](https://github.com/limserenahansol/BLA_scRNA_P10_P21_Mef2c/releases/tag/data-v1).

## Load and inspect in R

Use R 4.5 and Seurat 5 (validated with R 4.5.1, Seurat 5.3.0,
SeuratObject 5.1.0). Download the main RDS into your working directory.
No Cell Ranger run is needed to open these count-based objects.

```r
library(Seurat)
obj <- readRDS("BLA_E18_P0_P10_P21_seurat_v2.rds")
stopifnot(inherits(obj, "Seurat"), all(dim(obj) == c(25029L, 76446L)))
table(obj$age)
table(obj$sample, obj$broad_class)
table(obj$age, obj$population)
names(obj@reductions)
DimPlot(obj, reduction = "umap_source", group.by = "broad_class", raster = TRUE)
DimPlot(obj, reduction = "umap_raw", group.by = "source", raster = TRUE)
FeaturePlot(obj, reduction = "umap_source", features = "Mef2c", raster = TRUE)
counts <- LayerData(obj, assay = "RNA", layer = "counts")
metadata <- obj[[]]
```

R recognizes the gzip-compressed RDS directly; do not manually unzip it.
The sparse counts, normalized layer and metadata are already in the main object.
Loading uses more RAM than the downloaded file size: allow several GB for
inspection and substantially more for reprocessing (the pipeline recommends
at least 64 GB RAM). Do not densify the whole count matrix.

## Verify the downloaded files

Compare local hashes with `SEURAT_SHA256.txt`:

```powershell
Get-FileHash -Algorithm SHA256 -LiteralPath './BLA_E18_P0_P10_P21_seurat_v2.rds'
Get-FileHash -Algorithm SHA256 -LiteralPath './hansol2_combined_seurat.rds'
```

On Linux, with both files and the checksum manifest in the same directory:

```sh
sha256sum -c SEURAT_SHA256.txt
```

The processed export was previously roundtrip-tested with `readRDS()` against
its uncompressed copy: RNA counts, normalized data, metadata and every embedding
were unchanged. See `scripts/07_validate_release.R` for those checks.

## Reanalysis and interpretation

See [README.md](README.md) for the complete R pipeline and reference bundles,
and [METHODS.md](METHODS.md) for thresholds and validation. For rerunning all
steps, download the older source and annotated object from data-v1 too.
Preserve original inputs and keep outputs separate.

The biological replicate/animal IDs and new-input reference build are unknown.
Ten libraries do not establish ten independent animals. P21 is available only
from the older source, so age/source effects are confounded. Fine ITC/BLA
identities remain provisional; expression nominates experiments, not axonal
connections, direct MEF2 targets or a knockout mechanism. There is no labelled
KO/WT comparison in these files. Use RNA expression for gene analysis, not
Harmony coordinates; compare library-level summaries rather than treating
cells as independent animal replicates.

These downloads are publicly accessible. The repository's contact and requested
restrictions on redistribution/publication remain in its main README.
