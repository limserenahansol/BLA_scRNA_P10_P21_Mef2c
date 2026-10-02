"""Allen WMB-10Xv3 amygdala reference: per-subclass and per-supertype mean
log2 expression and fraction expressing, for cells dissected from ROI CTXsp
(LA/BLA/BMA/PA/CLA/EP) and sAMY (CEA/MEA/AAA/IA/BA).  Adult (P56) reference.

usage: python A01_allen_amygdala_reference.py STR|CTXsp
"""
import os
import sys, numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
from pathlib import Path

CACHE = Path(os.environ.get("ABC_ATLAS_CACHE", r"C:\Users\hsollim\Downloads\abc_atlas_cache"))
OUT = Path(os.environ.get("BLA_ROOT", r"K:\scRNA_BLA_phd\P10_P21_amygdala_2026-10")) / "ref"
part = sys.argv[1]
roi = {"STR": "sAMY", "CTXsp": "CTXsp"}[part]
h5 = CACHE / f"expression_matrices/WMB-10Xv3/20230630/WMB-10Xv3-{part}-log2.h5ad"

meta = pd.read_csv(CACHE / "metadata/WMB-10X/20231215/views/cell_metadata_with_cluster_annotation.csv",
                   usecols=["cell_label", "feature_matrix_label", "region_of_interest_acronym",
                            "class", "subclass", "supertype"])
meta = meta[(meta.feature_matrix_label == f"WMB-10Xv3-{part}") & (meta.region_of_interest_acronym == roi)]
meta = meta.set_index("cell_label")

a = ad.read_h5ad(h5, backed="r")
sym = a.var["gene_symbol"].astype(str).values
idx = np.where(a.obs_names.isin(meta.index))[0]
print(part, roi, "cells in ROI:", len(idx), "genes:", a.n_vars, flush=True)

levels = {"subclass": meta["subclass"], "supertype": meta["supertype"]}
groups = {lv: pd.Categorical(meta.loc[a.obs_names[idx], lv].values) for lv in levels}
sums = {lv: np.zeros((len(g.categories), a.n_vars)) for lv, g in groups.items()}
nnz = {lv: np.zeros((len(g.categories), a.n_vars)) for lv, g in groups.items()}

step = 20000
for s in range(0, len(idx), step):
    rows = idx[s:s + step]
    X = a.X[rows] if not sp.issparse(a.X) else a[rows].X
    X = sp.csr_matrix(X)
    B = X.copy(); B.data = np.ones_like(B.data)
    for lv, g in groups.items():
        codes = g.codes[s:s + step]
        G = sp.csr_matrix((np.ones(len(codes)), (codes, np.arange(len(codes)))),
                          shape=(len(g.categories), len(codes)))
        sums[lv] += (G @ X).toarray()
        nnz[lv] += (G @ B).toarray()
    print(f"  {min(s + step, len(idx))}/{len(idx)}", flush=True)

for lv, g in groups.items():
    n = pd.Series(g).value_counts().reindex(g.categories).values.astype(float)
    mean = pd.DataFrame(sums[lv] / n[:, None], index=g.categories, columns=sym).T
    pct = pd.DataFrame(nnz[lv] / n[:, None], index=g.categories, columns=sym).T
    mean = mean.groupby(level=0).max(); pct = pct.groupby(level=0).max()
    mean.to_parquet(OUT / f"allen_{part}_{lv}_mean_log2.parquet")
    pct.to_parquet(OUT / f"allen_{part}_{lv}_pct.parquet")
    pd.DataFrame({lv: g.categories, "n_cells": n, "roi": roi}).to_csv(OUT / f"allen_{part}_{lv}_n.csv", index=False)
print("done", part)
