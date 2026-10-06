# Developmental analysis v2 — 2026-10-06

Adds seven E18/P0/P10 count libraries to the previous three P10/P21 libraries.
173,431 input barcodes; 90,605 QC-pass; 76,446 putative singlets; 25,029 shared
genes. This is an exploratory developmental resource, not a KO/WT comparison.

Included: 18-slide PowerPoint and native PowerPoint-exported PDF, full numerical
summary tables, PNG/PDF figures, reproducible R/Python pipeline, documentation,
reference resources, session information and validation reports.

Validation: original-stream CRC and count/metadata checks; separate marker/gene
reference checks; raw/corrected embedding diagnostics; source/mito/doublet effect
sensitivity; independently reconstructed pseudobulk CPM; analytical versus Monte
Carlo detection standardization; round-trip compressed Seurat verification;
native editable chart workbooks and PowerPoint opening with no layout warnings.
Original files and data-v1 remain unchanged.

Important limitations: many ambiguous cells, low-count P0, 32.1% putative P10_2
doublets, poor atlas-only annotation agreement, incomplete source mixing,
age/source confounding and unknown animal independence. Gene expression does
not establish pruning, cue responses, direct MEF2 targets or the KO phenotype.
Core examples: Pcdh10 direction is consistent; ITC-like Kirrel3 direction reverses
between old/new P10. See the full per-library values and failed sensitivity checks.

Privacy: the existing repository is public. This release includes code, slides,
resources and **library-level summary results only**. New original/processed
cell-level scRNA data are prepared locally and are not included without an
explicit public-sharing decision. The 696 MB collaborator Seurat object and
supporting audit files are ready for direct/private sharing.

Download the results/pipeline ZIP and reference-resource ZIP below. Their
SHA-256 checksums are supplied. Start with `analysis_development_v2/README.md`.
