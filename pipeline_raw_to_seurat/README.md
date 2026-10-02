# E18 / P0 scRNA-seq — raw BCL → Cell Ranger → first-pass Seurat

10x Chromium 3' v3 libraries sequenced on a NextSeq 500 (NB500982), delivered as
raw Illumina run folders (BCL, not FASTQ).

| Run | Date | Flowcell | Libraries (SI-GA single index) |
|---|---|---|---|
| P165 | 2019-12-19 | H2GVNBGXF (high output) | P0_1 (E6), P0_2 (E9) |
| P173 | — | — | P10_1 (E10), P10_2 (E11) |
| P188 | 2020-02-12 | H5FVJBGXF (high output) | P0_3 (F2), E18_1 (F3), E18_2 (F4) |
| P189 | 2020-02-27 | H2NV5AFX2 (mid output)  | top-up of P188 libraries (checked by step 02) |

Read structure (all runs): R1 28 bp (barcode+UMI), I1 8 bp, R2 130 bp (cDNA).
`samples.csv` holds the sample → run → index mapping. A sample listed on two
runs gets its FASTQs pooled in `cellranger count`.

## Pipeline

| Step | Script | Where | What it does |
|---|---|---|---|
| 0 | `pipeline/00_install_tools.sh` | Linux/WSL | bcl2fastq2 (micromamba), Cell Ranger, mouse reference mm10-2020-A (matches the earlier BLA data) |
| 1 | `pipeline/01_untar_runs.sh` | Linux/WSL | extracts run tarballs (skips focus images) |
| 2 | `pipeline/02_detect_indices.sh` | Linux/WSL | converts 2 tiles, counts I1 and checks which SI-GA wells are really in each run vs `samples.csv` |
| 3 | `pipeline/03_demultiplex.sh` | Linux/WSL | BCL → FASTQ (4 oligos per SI-GA well → one sample), only indices present in that run |
| 4 | `pipeline/04_cellranger_count.sh` | Linux/WSL | `cellranger count` (introns included, no BAM) and copies shareable outputs to `results/cellranger/` |
| 5 | `pipeline/05_seurat_preprocess.R` | R (Windows or Linux) | QC, doublets, normalisation, Harmony, clusters, UMAP, neuron vs non-neuron |

All paths and resources are set in `pipeline/config.sh`. Run steps 1–4 with
`bash pipeline/run_all.sh 2>&1 | tee run_all.log`. Every step skips finished work.

Cell Ranger needs Linux. It does not run on native Windows. On this PC, use WSL2
(one-time setup from an **admin** PowerShell, then reboot):

```
wsl --install -d Ubuntu-24.04
```

By default WSL can use only half the RAM. Cell Ranger wants ≥64 GB, so create
`C:\Users\<you>\.wslconfig` with `[wsl2]` / `memory=110GB` / `processors=30`, then
run `wsl --shutdown`. Stanford Sherlock (`module load cellranger`) works the same
way. Set `TAR_DIR`, `WORK_DIR`, and `REF_DIR` in `config.sh`.

## Step 05 — what the Seurat object contains

`results/seurat/E18_P0_seurat_preprocessed.rds` (Seurat v5, layers `counts` + `data`)

* **Cell QC (per sample, logged in `tables/qc_cell_filtering_per_sample.csv`)**:
  ≥500 genes, ≥1000 UMIs, ≤10 % mito, ≤5 % haemoglobin, and log10 UMI no more
  than 5 MAD above the sample median. Doublets are called with scDblFinder per
  sample and removed. The thresholds are deliberately permissive. Change them
  with `--min_features=`, `--min_counts=`, or `--max_mt=` (use 5 for nuclei).
* LogNormalize → 3000 HVGs → PCA (50) → **Harmony on sample** → 30 PCs → Louvain
  (res 0.5) → UMAP. Reductions: `pca`, `harmony`, `umap` (Harmony),
  `umap.unintegrated` (before Harmony, for batch checks).
* **meta.data**: `sample`, `stage`, QC metrics, `doublet_score`, `cluster`,
  `broad_class`, `neuron_vs_non`, `neuron_type` (Glutamatergic/GABAergic), and
  per-cell module scores `score_*`.
* **Class calls** come from the cluster-level argmax of AddModuleScore over
  marker sets (`tables/marker_genes_used.csv`, `tables/cluster_class_scores.csv`).
  The neuron set uses genes that turn on early (Stmn2, Tubb3, Elavl3/4, Dcx,
  Gap43, Snap25), because E18/P0 neurons are largely immature (low Rbfox3). These
  are first-pass labels and are not curated annotation.

Figures (`results/seurat/figures/`, PNG + PDF): QC violins/scatter,
04 neuron vs non-neuron UMAP, 05 clusters/stage/sample, 06 split by sample,
07 pre-Harmony batch check, 08 neuron type, 09–10 marker and score feature plots,
11 marker dot plot, 12 composition per sample.

Tables (`results/seurat/tables/`): QC per sample, cluster markers (all + top 10
with class), composition, neuron/non-neuron counts, cell metadata.

## Hand-off to collaborator

Send `results/` (cellranger `*.h5`, `web_summary.html`, `metrics_summary.csv`
per sample, `all_samples_metrics_summary.csv`, and the `seurat/` folder) plus
this README and `samples.csv`. To start from counts instead of the Seurat object:

```r
m <- Seurat::Read10X_h5("results/cellranger/P0_1/filtered_feature_bc_matrix.h5")
```
