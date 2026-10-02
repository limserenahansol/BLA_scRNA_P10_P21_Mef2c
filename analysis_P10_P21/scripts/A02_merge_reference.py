"""Merge CTXsp + sAMY Allen centroids, keep groups with >=30 cells, pick
marker genes (one-vs-rest within the reference), write CSVs for R mapping."""
import os
import numpy as np, pandas as pd
from pathlib import Path
R = Path(os.environ.get("BLA_ROOT", r"K:\scRNA_BLA_phd\P10_P21_amygdala_2026-10")) / "ref"
for lv, min_n in [("subclass", 30), ("supertype", 30)]:
    means, pcts, ns = [], [], []
    for part in ["CTXsp", "STR"]:
        m = pd.read_parquet(R / f"allen_{part}_{lv}_mean_log2.parquet")
        p = pd.read_parquet(R / f"allen_{part}_{lv}_pct.parquet")
        n = pd.read_csv(R / f"allen_{part}_{lv}_n.csv").set_index(lv)["n_cells"]
        means.append(m * n.values); pcts.append(p * n.values); ns.append(n)
    genes = means[0].index.intersection(means[1].index)
    N = pd.concat(ns, axis=1).fillna(0).sum(axis=1)
    S = sum(x.reindex(index=genes, columns=N.index).fillna(0) for x in means)
    P = sum(x.reindex(index=genes, columns=N.index).fillna(0) for x in pcts)
    keep = N[N >= min_n].index
    M = (S[keep] / N[keep]); P = (P[keep] / N[keep])
    roi = pd.concat([pd.read_csv(R / f"allen_{p}_{lv}_n.csv").assign(part=p) for p in ["CTXsp", "STR"]])
    roi = roi.pivot_table(index=lv, columns="roi", values="n_cells", aggfunc="sum").reindex(keep).fillna(0)
    # markers: per group top 40 genes by (mean - max of others), detected in >=40% of the group
    mk = set()
    for g in keep:
        other = M.drop(columns=g).max(axis=1)
        sc = (M[g] - other)[P[g] >= 0.4]
        mk |= set(sc.sort_values(ascending=False).head(40).index)
    mk = sorted(mk)
    M.to_csv(R / f"ref_{lv}_mean_log2_all_genes.csv.gz")
    P.to_csv(R / f"ref_{lv}_pct_all_genes.csv.gz")
    M.loc[mk].to_csv(R / f"ref_{lv}_mean_log2_markers.csv")
    pd.DataFrame({"n_cells": N[keep]}).join(roi).to_csv(R / f"ref_{lv}_ncells.csv")
    print(lv, "groups:", len(keep), "marker genes:", len(mk))
