"""Adult Allen WMB Isocortex (all 10Xv2 + 10Xv3 shards, precomputed per-subclass
sums by the GPCR project): % cells and mean log2 for Mef2c and the axon-side
partners of the BLA/ITC candidate cues, across cortical glutamatergic subclasses."""
import numpy as np, pandas as pd, glob, os
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = os.environ.get("BLA_ROOT", r"K:\scRNA_BLA_phd\P10_P21_amygdala_2026-10")
S = N = Z = None
for f in glob.glob(os.path.join(os.environ.get("ALLEN_GENOMEWIDE_DIR", r"C:\Users\hsollim\Desktop\cursor\Genelist_analysis_WMB\v3\outputs\allen_genomewide"), "Isocortex*.npz")):
    z = np.load(f, allow_pickle=True)
    lab = list(z["subclass_labels"]); genes = list(z["genes"])
    s = pd.DataFrame(z["subclass_sum"], index=lab, columns=genes)
    nz = pd.DataFrame(z["subclass_nnz"], index=lab, columns=genes)
    n = pd.Series(z["subclass_n"], index=lab)
    S = s if S is None else S.add(s, fill_value=0); Z = nz if Z is None else Z.add(nz, fill_value=0)
    N = n if N is None else N.add(n, fill_value=0)
glut = [l for l in N.index if l.endswith("Glut") and N[l] >= 200]
order = ["007 L2/3 IT CTX Glut", "006 L4/5 IT CTX Glut", "005 L5 IT CTX Glut", "004 L6 IT CTX Glut",
         "022 L5 ET CTX Glut", "032 L5 NP CTX Glut", "030 L6 CT CTX Glut", "029 L6b CTX Glut",
         "003 L5/6 IT TPE-ENT Glut", "002 IT EP-CLA Glut", "001 CLA-EPd-CTX Car3 Glut"]
glut = [g for g in order if g in glut] + [g for g in glut if g not in order]
genes = {"Mef2c": "TF (KO gene)",
         "Robo1": "Slit1-3", "Robo2": "Slit1-3", "Lrrc4c": "Ntng1", "Lrrc4": "Ntng2", "Adgrb3": "C1ql3",
         "Nrxn2": "Igsf21", "Ptprg": "Cntn4/5/6", "Chl1": "Cntn6", "Plxnd1": "Sema3e", "Nrp1": "Sema3a",
         "Plxna1": "Sema5a/5b/6d", "Plxna3": "Sema5a/5b", "Plxna4": "Sema6a",
         "Cdh8": "Cdh8 (homophilic)", "Cdh9": "Cdh9 (homophilic)", "Cdh13": "Cdh13 (homophilic)",
         "Cdh18": "Cdh18 (homophilic)", "Cdh6": "Cdh6 (homophilic)", "Kirrel3": "Kirrel3 (homophilic)",
         "Sdk2": "Sdk2 (homophilic)", "Pcdh7": "Pcdh7 (homophilic)"}
g = [x for x in genes if x in S.columns]
pct = (Z.loc[glut, g].T / N[glut]).T * 100
mean = (S.loc[glut, g].T / N[glut]).T
out = pct.round(1).T; out.insert(0, "binds_BLA_cue", [genes[x] for x in g])
out.to_csv(os.path.join(ROOT, "tables", "C01_cortex_glut_receptor_pct_adult_Allen.csv"))
mean.T.round(3).to_csv(os.path.join(ROOT, "tables", "C01_cortex_glut_receptor_meanlog2_adult_Allen.csv"))
fig, ax = plt.subplots(figsize=(15, 6.2))
yy, xx = np.meshgrid(range(len(glut)), range(len(g)), indexing="ij")
z = (mean - mean.mean()) / (mean.std() + 1e-9)
sc = ax.scatter(xx.ravel(), yy.ravel(), s=(pct.values.ravel() / 100) * 260, c=z.values.ravel(), cmap="RdBu_r", vmin=-2, vmax=2, edgecolors="none")
ax.set_xticks(range(len(g))); ax.set_xticklabels([f"{x}\n({genes[x]})" for x in g], rotation=70, ha="right", fontsize=8, style="italic")
ax.set_yticks(range(len(glut))); ax.set_yticklabels([f"{l}  (n={int(N[l]):,})" for l in glut], fontsize=9); ax.invert_yaxis()
ax.axvline(0.5, color="k", lw=0.8)
for p in (25, 50, 100): ax.scatter([], [], s=p / 100 * 260, c="grey", label=f"{p}%")
ax.legend(title="% cells", loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False)
plt.colorbar(sc, ax=ax, fraction=0.02, pad=0.12, label="z (mean log2 across subclasses)")
ax.set_title("Adult cortex (Allen WMB Isocortex): Mef2c and axon-side partners of BLA/ITC cues, by glutamatergic subclass", fontsize=11, loc="left")
ax.grid(alpha=0.2); plt.tight_layout()
for ext in ("png", "pdf"): plt.savefig(os.path.join(ROOT, "figures", f"F9_cortex_receptors_adult_Allen.{ext}"), dpi=220, bbox_inches="tight")
print(out.to_string())
