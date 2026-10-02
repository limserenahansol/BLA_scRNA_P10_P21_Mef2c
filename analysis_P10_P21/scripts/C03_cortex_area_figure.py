"""Scenario 3, cortex side: by cortical area (Allen ROI) and IT/ET subclass, the
fraction of adult neurons carrying the partners of BLA/ITC cues, and how those
genes co-vary with Mef2c inside L2/3 and L4/5 IT neurons."""
import os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = os.environ.get("BLA_ROOT", r"K:\scRNA_BLA_phd\P10_P21_amygdala_2026-10")
d = pd.read_parquet(os.path.join(ROOT, "ref", "cortex_IT_ET_candidate_genes_v2.parquet"))
d = d[d.Snap25 > 0]
AREAS = ["SSp", "SS-GU-VISC", "MOp", "MO-FRP", "AUD", "VIS", "VIS-PTLp", "RSP", "ACA", "PL-ILA-ORB", "AI", "TEa-PERI-ECT"]
BLA_INPUT = {"AI", "TEa-PERI-ECT", "PL-ILA-ORB", "ACA"}   # areas with known direct BLA input (rodent tracing literature)
SUBS = ["007 L2/3 IT CTX Glut", "006 L4/5 IT CTX Glut", "005 L5 IT CTX Glut", "004 L6 IT CTX Glut", "022 L5 ET CTX Glut"]
d["Plxnd1+Nrp1+"] = ((d.Plxnd1 > 0) & (d.Nrp1 > 0)).astype(float)
d["Plxnd1+Nrp1-"] = ((d.Plxnd1 > 0) & (d.Nrp1 == 0)).astype(float)
show = ["Mef2c", "Plxnd1", "Plxnd1+Nrp1+", "Plxnd1+Nrp1-", "Nrp1", "Kirrel3", "Sdk2", "Cdh8", "Cdh9", "Cdh18", "Robo1", "Robo2", "Lrrc4c"]
rows = []
for sub in SUBS[:2]:
    for a in AREAS:
        x = d[(d.subclass == sub) & (d.region_of_interest_acronym == a)]
        if len(x) < 100: continue
        r = {"subclass": sub, "area": a, "n": len(x)}
        for g in show:
            r[g] = 100 * (x[g] > 0).mean() if g in x and "+" not in g else 100 * x[g].mean()
        rows.append(r)
T = pd.DataFrame(rows); T.to_csv(os.path.join(ROOT, "tables", "C03_cortex_area_partner_pct.csv"), index=False)

fig, axes = plt.subplots(1, 3, figsize=(20, 7.2), gridspec_kw={"width_ratios": [1, 1, 0.8]})
for ax, sub in zip(axes[:2], SUBS[:2]):
    t = T[T.subclass == sub].set_index("area").reindex([a for a in AREAS if a in set(T[T.subclass == sub].area)])
    M = t[show].values
    Z = (M - M.mean(0)) / (M.std(0) + 1e-9)
    im = ax.imshow(Z, cmap="RdBu_r", vmin=-2, vmax=2, aspect="auto")
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            ax.text(j, i, f"{M[i, j]:.0f}", ha="center", va="center", fontsize=7.5)
    ax.set_xticks(range(len(show))); ax.set_xticklabels(show, rotation=60, ha="right", fontsize=9, style="italic")
    ax.set_yticks(range(len(t))); ax.set_yticklabels([f"{a}{'  *' if a in BLA_INPUT else ''}  (n={int(n):,})" for a, n in zip(t.index, t.n)], fontsize=9)
    ax.set_title(f"{sub}: % cells (colour = z across areas)", fontsize=10.5, loc="left")
plt.colorbar(im, ax=axes[1], fraction=0.03, pad=0.02, label="z across areas")

# co-variation with Mef2c inside L2/3 + L4/5 IT, area and supertype regressed out
x = d[d.subclass.isin(SUBS[:2])].copy()
genes = ["Plxnd1", "Nrp1", "Nrp2", "Plxna1", "Plxna4", "Kirrel3", "Sdk2", "Cdh8", "Cdh9", "Cdh13", "Cdh18", "Cdh6", "Robo1", "Robo2", "Lrrc4c",
         "Pcdh10", "Arc", "Homer1", "Nr4a1", "Bdnf", "Satb2", "Cux2", "Rorb", "Lmo4"]
D = pd.get_dummies(x[["supertype", "region_of_interest_acronym", "donor_label"]].astype(str), drop_first=True).astype(float)
D["snap25"] = x.Snap25.values; D["c"] = 1.0
A = D.values
def resid(v):
    b, *_ = np.linalg.lstsq(A, v, rcond=None); return v - A @ b
rm = resid(x.Mef2c.values.astype(float))
cor = {g: np.corrcoef(rm, resid(x[g].values.astype(float)))[0, 1] for g in genes}
c = pd.Series(cor).sort_values()
c.to_csv(os.path.join(ROOT, "tables", "C03_Mef2c_covariation_L23_L45_IT_adult.csv"), header=["r_residual"])
axes[2].barh(c.index, c.values, color=["#B2182B" if v > 0 else "#2166AC" for v in c.values])
axes[2].axvline(0, color="k", lw=0.8); axes[2].tick_params(axis="y", labelsize=9)
axes[2].set_title(f"Co-variation with Mef2c in adult L2/3 + L4/5 IT\n(n={len(x):,}; supertype, area, donor, depth removed)", fontsize=10.5, loc="left")
axes[2].set_xlabel("residual correlation r")
fig.suptitle("Cortex side of scenario 3 (adult Allen WMB 10Xv2). * = area with known direct input to BLA", fontsize=12, x=0.01, ha="left")
plt.tight_layout()
for ext in ("png", "pdf"): plt.savefig(os.path.join(ROOT, "figures", f"F12_scenario3_cortex_by_area.{ext}"), dpi=200, bbox_inches="tight")
print(T.round(1).to_string(index=False)); print(c.round(3).to_string())
