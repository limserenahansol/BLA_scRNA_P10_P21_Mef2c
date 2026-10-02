# BLA scRNA-seq across development (E18, P0, P10, P21): cell types, TFs and cues for the Mef2c-cKO question

Hansol Lim (Stanford). 10x Chromium 3' v3 single-cell RNA-seq of mouse basolateral amygdala.

This repository contains:

1. **Pipeline code**
   - Raw Illumina BCL → FASTQ → Cell Ranger → first-pass Seurat (`pipeline_raw_to_seurat/`)
   - The P10 / P21 analysis (`analysis_P10_P21/`)
2. **Data** to check or re-analyse, as GitHub Release assets (see [Data](#data); the files are too large for git)
3. **Slides**: `docs/BLA_P10_P21_celltypes_TF_cues_Mef2c_scenarios.pptx` (18 slides) and all tables in one workbook, `docs/BLA_P10_P21_tables_for_collaborator.xlsx`
4. **Results**: every figure (`results/figures/*.png`) and table (`results/tables/*.csv`)

> **Status, Oct 2026.** P10 (n = 1) and P21 (n = 2) are analysed. E18 (n = 2), P0 (n = 3) and two more P10 animals are raw BCL that still has to be run through `pipeline_raw_to_seurat/`.

---

## The biological question

In Mef2c conditional KO mice, L2/3 callosal axons stop targeting the contralateral S1/S2. Ectopic axons appear around the BLA and FoxP2+ intercalated cells (ITC) from about P5, and they increase and branch by P8–P10. The neurons they come from may not be in SSc.

Plan: look at the amygdala first, then go back to cortex.

1. Which amygdala cell populations exist at P10, and how do they differ from P21?
2. Which transcription factors, Mef2c included, mark each population?
3. Which guidance or adhesion cues do BLA/ITC cells express in the P5–P10 window?
4. Which cortical neurons carry the matching partner receptors?

## Main results (P10 vs P21)

| | Finding | Figure |
|---|---|---|
| Cell types | 35 types in 13,625 cells: 10 BLA/BMA glutamatergic, 2 ITC (Foxp2+ Tshz1+ Pbx3+), CeA/MeA GABA, 9 interneuron types, glia | F1–F3 |
| BLA glut | Rspo2/Etv1 group (92% maps to Allen 014 LA-BLA-BMA-PA Glut), plus Tshz2/Satb1, Cdh8/Tshz3, Fgf10/Nrp1, Otof/Trhr, Vgll3/Prr16 | F3 |
| Mef2c | >90% of cells in Tshz2/Satb1, Cdh8/Tshz3 and Otof/Trhr BLA neurons, Pvalb chandelier cells and microglia; 50–60% in ITC and Rspo2 BLA | F6a, F8 |
| Scenario 2 (pruning) | P10 BLA is pre-pruning: astrocyte Mertk 30% vs 75%, Megf10 15% vs ~30%, microglial C1q 8% vs 17–37%. Target neurons express a P10-high adhesion programme (Ncam1, Sema6d, Cadm2, Pcdh7) | F10 |
| Scenario 3 (cues) | ITC: Kirrel3, Sdk2, Sema6d, Sema5b, Cdh18. Rspo2 BLA: C1ql3, Cdh8, Cdh9, Sema3e, Slit2 | F7c, F11 |
| Cortex side (adult Allen) | Cdh9 is in 61–64% of L2/3 IT cells in BLA-input areas (TEa-PERI-ECT, PL-ILA-ORB), 46% in AI, but 9% in SSp. SSp L4/5 IT is mostly Plxnd1+ Nrp1− (repulsion-type) | F12 |

**Caveats.**
- P10 is one animal, so P10 vs P21 is descriptive.
- The Allen reference and the cortex data are adult.
- Correlations with Mef2c are weak (|r| < 0.12).

All of these results are hypotheses to test, for example by HCR at P5–P10 or in the cKO.

---

## Data

Download from the **Releases** page of this repository (tag `data-v1`). Put the files in `data/`.

| File | What it is | Size |
|---|---|---|
| `BLA_P10_P21_annotated_seurat.rds` | **Main object.** Seurat v5, P10 + P21, 13,625 cells, final `cell_type` labels, Allen mapping, UMAP and Harmony | ~190 MB |
| `BLA_P10_P21_counts_10x.zip` | The same cells as 10x files (`matrix.mtx.gz`, `features.tsv.gz`, `barcodes.tsv.gz`) plus `metadata.csv` with labels and UMAP coordinates. Works in any tool | ~70 MB |
| `BLA_P10_P21_P56_QC_seurat_2023.rds` | The input object: P10, P21 and P56 after the original QC (T. Straub 2023; Cell Ranger 7, mm10-2020-A, DoubletFinder, Harmony) | ~610 MB |
| `BLA_glut_subclusters.rds`, `BLA_gaba_subclusters.rds` | Re-clustered glutamatergic and GABAergic neurons, with their own UMAPs | 85 / 60 MB |
| `allen_amygdala_reference.zip` | Allen WMB amygdala centroids (ROIs CTXsp + sAMY, 84 subclasses / 233 supertypes) and per-cell cortex values used for F12 | ~250 MB |
| `figures_pdf.zip` | Vector PDF of every figure | small |

Load the main object:

```r
library(Seurat)
seu <- readRDS("data/BLA_P10_P21_annotated_seurat.rds")
table(seu$cell_type, seu$age)
DimPlot(subset(seu, age == "P10"), group.by = "cell_type", label = TRUE) + NoLegend()
DotPlot(seu, features = c("Mef2c", "Foxp2", "Rspo2", "Kirrel3"), group.by = "cell_type")
```

Or in Python:

```python
import scanpy as sc, pandas as pd
ad = sc.read_10x_mtx("data/BLA_P10_P21_counts_10x")
ad.obs = pd.read_csv("data/BLA_P10_P21_counts_10x/metadata.csv", index_col=0).loc[ad.obs_names]
```

Key metadata columns:
- `sample`: P10s1, P21s1, P21s2
- `age`: P10 or P21
- `cell_type`: 35 labels
- `class`: Glutamatergic, GABAergic or Non-neuronal
- `allen_subclass`, `allen_supertype`, `allen_subclass_r`: correlation mapping to the Allen reference
- `glut_umap1/2`: UMAP coordinates for the glutamatergic-only map
- `cl`: whole-data cluster, resolution 1

---

## Pipeline 1 — raw BCL → Seurat (`pipeline_raw_to_seurat/`)

Use this for the E18 / P0 / new-P10 runs, which arrive as NextSeq run folders (`P165_*.tar.gz`, …). Details are in `pipeline_raw_to_seurat/README.md`.

| Step | Script | Does |
|---|---|---|
| 0 | `00_install_tools.sh` | bcl2fastq2 (micromamba), Cell Ranger, mouse reference **mm10-2020-A** (the same reference as the older data) |
| 1 | `01_untar_runs.sh` | extract the run tarballs |
| 2 | `02_detect_indices.sh` + `detect_10x_index.py` | convert 2 tiles and check which SI-GA wells are really in each run against `samples.csv` |
| 3 | `03_demultiplex.sh` | BCL → FASTQ (4 oligos per SI-GA well) |
| 4 | `04_cellranger_count.sh` | `cellranger count`. A library sequenced on two runs, e.g. P188 + P189, is pooled |
| 5 | `05_seurat_preprocess.R` | QC, scDblFinder, Harmony, clusters, neuron vs non-neuron UMAP |

- Steps 1–4 need Linux: WSL2, a cluster such as `server/run_on_server.sbatch`, or Sherlock.
- Step 5 runs in R on any OS.
- Sample → run → index mapping is in `samples.csv`, e.g. P0_1 = SI-GA-E6, P0_2 = E9, P10_1 = E10, P10_2 = E11, P0_3 = F2, E18_1 = F3, E18_2 = F4.

```bash
cd pipeline_raw_to_seurat/pipeline
bash 00_install_tools.sh /path/to/cellranger-x.y.z.tar.gz
bash run_all.sh 2>&1 | tee run_all.log
```

## Pipeline 2 — P10 / P21 analysis (`analysis_P10_P21/scripts/`)

Set the working folder once. It must contain `data/`, `ref/`, `resources/`, `tables/` and `figures/`.

```bash
export BLA_ROOT=/path/to/workdir            # R and Python scripts read this
export ABC_ATLAS_CACHE=/path/to/abc_atlas_cache   # only for the A01 / C02 Allen steps
```

| Script | Does | Main output |
|---|---|---|
| `A01_allen_amygdala_reference.py STR` / `CTXsp`, `A02_merge_reference.py` | Allen WMB amygdala reference centroids. Skip these if you unzip `allen_amygdala_reference.zip` into `ref/` | `ref/ref_*` |
| `B01_p10_p21_process.R` | P10 + P21 from the input object, then Seurat v5, Harmony and clusters | `B01_*.rds` |
| `B02_map_to_allen.R` | correlation mapping to Allen subclass and supertype | `B02_*.rds`, cluster table |
| `B03_cluster_markers.R` | curated amygdala markers and top markers per cluster | `tables/B03_*` |
| `B04_glut_subcluster.R`, `B05_gaba_subcluster.R` | neuron re-clustering | `B04/B05_*.rds` |
| `B06_annotate.R` | final labels; removes doublets and low-quality cells (logged) | `B06_P10_P21_annotated.rds` |
| `B07_fig_celltypes.R` | F1–F4: UMAP, markers, glut subtypes, composition | figures |
| `B08_fig_dev_tf_cues.R` | F5–F8: DE P10 vs P21, TFs, guidance cues, Mef2c | figures, tables |
| `B09_candidate_cues.R` | F7c: cue ranking for ITC and BLA glut | `tables/B09_*` |
| `B10_scenarios_BLA.R` | F10–F11: pruning window, cue dynamics | figures |
| `C01`–`C03_*.py` | F9, F12: adult Allen cortex partners by subclass and area | figures |
| `D01_build_deck.py`, `E01_export_10x.R` | slides; 10x export | |

Requirements:
- R ≥ 4.3: Seurat 5, harmony, dplyr, ggplot2, patchwork, tidyr, ggrepel, R.utils
- Python ≥ 3.10: pandas, numpy, anndata, scipy, matplotlib, pyarrow, python-pptx

## Methods in brief

- **QC.** The input object had mito < 20%, hb < 20% and DoubletFinder singlets. The P10/P21 re-analysis also removed clusters that co-express two lineages (e.g. glutamatergic + ITC, or neuron + microglia/OPC) and low-gene clusters: 3,079 of 16,704 cells (`results/tables/B06_removed_cells.csv`).
- **Clustering.** Seurat v5: LogNormalize → 3,000 HVG → PCA → Harmony on sample → Louvain → UMAP.
- **Labels.** Marker genes, cross-checked by Pearson correlation of each cell and each cluster with Allen WMB adult amygdala centroids (ROIs CTXsp + sAMY, about 199k cells, 2,983 marker genes).
- **P10 vs P21 DE.** Wilcoxon per cell type, kept only when the pseudo-bulk direction holds against both P21 replicates.
- **TFs.** 1,321 curated mouse TFs (`analysis_P10_P21/resources/mouse_TF_list_1321.csv`).

## Citation / contact

Unpublished data: please do not redistribute or publish without asking. Hansol Lim, hsollim@stanford.edu.

References used in the slides:
- Yao et al. 2023, Nature (Allen WMB atlas)
- Kim et al. 2016, Nat Neurosci (Rspo2 / Ppp1r1b BLA)
- Flavell et al. 2006, Science; Pfeiffer et al. 2010, Neuron; Tsai et al. 2012, Neuron (MEF2 and synapse elimination)
- Chauvet et al. 2007, Neuron (Sema3e–Plxnd1/Nrp1)
- Martin et al. 2015, eLife (Kirrel3)
- Alcamo et al. 2008, Neuron; Britanova et al. 2008, Neuron (Satb2)
