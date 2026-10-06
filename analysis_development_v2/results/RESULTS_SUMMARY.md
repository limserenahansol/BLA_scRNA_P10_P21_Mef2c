# Developmental v2: exploratory findings

## Conclusion

The new E18/P0/P10 counts-only object is now included alongside previous
P10/P21 libraries. Broad marker-supported populations and TF/candidate expression
are available, but these data cannot identify the cause of a MEF2C-KO wiring
phenotype. Cell census and fine developmental identities are not definitive.

## Data and checks

10 libraries across E18/P0/P10/P21, 173,431 input barcodes, 90,605 shared-gene
QC-pass barcodes, **76,446 putative singlets**, 25,029 common genes.
No labelled KO/WT contrast or stored animal IDs. Original inputs are preserved.
Complete original gzip streams and count/metadata integrity pass independent checks.
Large generated checkpoints are uncompressed and round-trip validated after a
compressed checkpoint failed CRC/reload. Count and normalized-RNA hashes are
unchanged after embedding correction and after slimming the collaborator object.

| Library | QC-pass | Putative singlets | Important caveat |
|---|---:|---:|---|
| E18_1 | 3,779 | 3,332 | Adult subtype identities tentative |
| E18_2 | 4,018 | 3,604 | Adult subtype identities tentative |
| P0_1 | 7,784 | 6,943 | Low-count input; only 35.8% retained |
| P0_2 | 4,888 | 4,296 | Only 8.0% of input barcodes pass QC |
| P0_3 | 4,604 | 3,960 | Better count depth than other P0 libraries |
| P10_1 | 33,548 | 31,017 | Very large capture; many ambiguous cells |
| P10_2 | 17,171 | 11,654 | 32.1% classifier-labelled doublets |
| P10s1 | 6,546 | 5,086 | Very low mitochondrial capture; many ambiguous cells |
| P21s1 | 2,824 | 2,287 | Previous source only |
| P21s2 | 5,443 | 4,267 | Previous source only |

Repeated P10_2 doublet classification yields identical calls and scores; captured
warnings concern deprecated xgboost arguments. Calls remain putative because
loading/cell-calling/background information is unavailable.

## Populations and annotation uncertainty

Marker-supported candidates include glutamatergic and GABAergic neurons,
immature neurons, astrocytes, OPCs/oligodendrocytes, microglia, vascular and other
non-neuronal cells. ITC-like GABA and Rspo2-positive glut populations are nested,
provisional subsets; they are not added to tissue composition denominators.

The adult-atlas-only baseline has poor independent marker agreement (3.8–58.3%
across QC-pass libraries). The previous developmental-reference broad gene
splits agree more often (59.1–97.2%), but inherit v1 annotation uncertainty.
Leave-one-old-library-out fine-label recovery is 84.6–90.5% against inherited
labels, not independent ground truth. Among primary putative singlets,
45.7–47.3% of new P10 and 70.3% of old P10 remain ambiguous. No forced subtype
claim is made for these cells. E18 annotation is particularly developmental.

Source Harmony increases balanced P10 old/new neighbor mixing from 3.4% to
10.7%, **still far from sufficient mixing**. Independent broad-class neighbor
retention changes from 55.9% to 55.7%. Passing the predeclared five-percentage-point
loss gate makes the corrected map an available visualization, not validated
integration. Raw and corrected maps are both supplied.

## Core candidate findings, including failures

Values below are `log2((P10 library-mean CPM+0.5)/(P21 library-mean CPM+0.5))`.
Old/new P10 are kept separate. No cell-level p-values, animal-level significance,
KO effect, or direct transcriptional regulation is inferred.

| Gene / population | Old P10 vs P21 | New P10 vs P21 | Source/QC interpretation |
|---|---:|---:|---|
| Mef2c, glutamatergic | −0.07 | +0.46 | Not direction-consistent |
| Mef2c, GABAergic | +0.19 | +0.53 | Positive in all source/QC/library-pair checks; modest effect |
| Pcdh10, glutamatergic | +0.64 | +1.31 | Positive in all direction checks; candidate for follow-up |
| Kirrel3, ITC-like | +0.43 | −0.41 | Reverses across source; not a stable developmental increase |
| Sema3e, glutamatergic | +0.52 | +0.04 | New P10/P21 library pairs are not direction-consistent |
| Plxnd1, ITC-like | +2.37 | +2.79 | Direction-consistent; inspect sparse counts/identity/ambient RNA |
| Nrp1, glutamatergic | +0.57 | −0.01 | Not direction-consistent |

Pcdh10 also retains a positive direction in GABAergic and provisional ITC-like
groups. Direction consistency is a robustness diagnostic, **not** statistical
significance or proof of pruning. Marker-dependent selection, sparse genes,
ambient RNA and broad-class uncertainty can change an apparent effect.
All candidates and failures are in the per-library and sensitivity tables.

The P10 TF figure includes prespecified lineage/developmental genes and the
supplied mouse TF screen. It describes transcription, not TF activity. Adult
cortical receptor-expression context from v1 is separately labelled adult
reference, not the developing cortical projection neurons in this experiment.

## What the collaborator should test next

First establish matched-age KO/WT axon/bouton/synapse changes normalized to
cortical labeling and injection. Then spatially verify target cells and identify
the actual cortical projection neurons. Only then test guidance/matching versus
branch/synapse retention using genotype-specific expression/protein measurements
and a focused perturbation/rescue. These are proposed tests, not observed results.

Before confirmatory scRNA analysis, obtain animal/sex/genotype/capture/reference
metadata and Cell Ranger reports/unfiltered droplet matrices. The present object
does not support validated empty-droplet or ambient-RNA correction.

See [METHODS.md](../METHODS.md) for full definitions and primary references,
[README.md](../README.md) for data usage, and the 18-slide PowerPoint in this folder.
