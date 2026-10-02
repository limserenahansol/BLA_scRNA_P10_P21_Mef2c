"""Per-cell values of candidate genes in adult Allen cortical IT/ET neurons,
import os
from the 10Xv2 Isocortex shards (the only chemistry covering every area,
including SSp and TEa-PERI-ECT). Output: one parquet with cells x genes + ROI."""
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
from pathlib import Path
CACHE = Path(os.environ.get("ABC_ATLAS_CACHE", r"C:\Users\hsollim\Downloads\abc_atlas_cache"))
OUT = Path(os.environ.get("BLA_ROOT", r"K:\scRNA_BLA_phd\P10_P21_amygdala_2026-10")) / "ref" / "cortex_IT_ET_candidate_genes_v2.parquet"
GENES = ["Mef2c", "Mef2a", "Mef2d", "Plxnd1", "Nrp1", "Nrp2", "Plxna1", "Plxna4", "Kirrel3", "Sdk2", "Cdh8", "Cdh9",
         "Cdh13", "Cdh18", "Cdh6", "Robo1", "Robo2", "Lrrc4c", "Pcdh10", "Arc", "Homer1", "Nr4a1", "Bdnf",
         "Satb2", "Cux2", "Rorb", "Bcl11b", "Fezf2", "Lmo4", "Snap25"]
SUB = ["007 L2/3 IT CTX Glut", "006 L4/5 IT CTX Glut", "005 L5 IT CTX Glut", "004 L6 IT CTX Glut", "022 L5 ET CTX Glut",
       "003 L5/6 IT TPE-ENT Glut"]
meta = pd.read_csv(CACHE / "metadata/WMB-10X/20231215/views/cell_metadata_with_cluster_annotation.csv",
                   usecols=["cell_label", "feature_matrix_label", "region_of_interest_acronym", "subclass", "supertype", "donor_label"],
                   low_memory=False)
meta = meta[meta.subclass.isin(SUB) & meta.feature_matrix_label.str.startswith("WMB-10Xv2-Isocortex")].set_index("cell_label")
parts = []
for k in range(1, 5):
    a = ad.read_h5ad(CACHE / f"expression_matrices/WMB-10Xv2/20230630/WMB-10Xv2-Isocortex-{k}-log2.h5ad", backed="r")
    sym = pd.Index(a.var["gene_symbol"].astype(str))
    gi = [int(np.where(sym == g)[0][0]) for g in GENES if g in set(sym)]
    gn = [sym[i] for i in gi]
    rows = np.where(a.obs_names.isin(meta.index))[0]
    print(f"shard {k}: {len(rows)} cells", flush=True)
    for s in range(0, len(rows), 30000):
        r = rows[s:s + 30000]
        X = a[r].X
        X = X[:, gi].toarray() if sp.issparse(X) else np.asarray(X)[:, gi]
        df = pd.DataFrame(X.astype(np.float32), index=a.obs_names[r], columns=gn)
        parts.append(df.join(meta[["region_of_interest_acronym", "subclass", "supertype", "donor_label"]]))
        print(f"  {min(s + 30000, len(rows))}/{len(rows)}", flush=True)
pd.concat(parts).to_parquet(OUT)
print("done", OUT)
