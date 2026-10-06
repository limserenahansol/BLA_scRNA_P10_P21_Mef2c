"""Readable descriptive figures; all numbers come from saved library-level tables."""
from pathlib import Path
import json, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D

ROOT = Path(os.environ.get("BLA_V2_ROOT", "analysis_development_v2")).resolve()
T = ROOT / "results/tables"; F = ROOT / "results/figures"; F.mkdir(parents=True, exist_ok=True)
ORDER = ["E18_1", "E18_2", "P0_1", "P0_2", "P0_3", "P10_1", "P10_2", "P10s1", "P21s1", "P21s2"]
AGES = ["E18", "P0", "P10", "P21"]
COLORS = dict(zip(["Glutamatergic", "GABAergic", "Immature_neuron", "Astrocyte", "Progenitor", "OPC", "Oligodendrocyte", "Microglia", "Endothelial", "Pericyte", "Fibroblast", "Ependymal", "Choroid", "Erythroid", "Ambiguous"],
                 ["#D97706", "#5367B4", "#A259A6", "#00A896", "#81B29A", "#70A0D5", "#35638C", "#A05A2C", "#D44C74", "#A6814C", "#807366", "#6A9D4E", "#B8A840", "#B33030", "#C3C7CC"]))
plt.rcParams.update({"font.family": "Arial", "font.size": 12, "axes.titlesize": 16, "axes.labelsize": 13,
                     "pdf.fonttype": 42, "ps.fonttype": 42, "axes.spines.top": False, "axes.spines.right": False})
def read(name): return pd.read_csv(T / name, low_memory=False)
def save(fig, name):
    fig.savefig(F / f"{name}.png", dpi=160, bbox_inches="tight", facecolor="white")
    fig.savefig(F / f"{name}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)

q = read("sample_qc_summary.csv").set_index("sample").reindex(ORDER)
db = read("doublet_summary.csv").set_index("sample").reindex(ORDER)
fig, ax = plt.subplots(1, 3, figsize=(17, 5))
x = np.arange(len(ORDER))
for j, (col, label, color) in enumerate([(q.input_barcodes, "Input barcodes", "#A8B0B9"), (q.qc_primary, "QC-pass", "#3B82A6"), (db.singlets, "Putative singlets", "#087E8B")]):
    ax[0].bar(x + (j-1)*.25, col, width=.24, label=label, color=color)
ax[0].set_yscale("log"); ax[0].set_ylabel("Barcodes (log scale)"); ax[0].set_title("Retention, not tissue abundance")
ax[0].legend(fontsize=10)
ax[1].bar(x-.17, q.median_mt_before, width=.34, color="#A8B0B9", label="Before QC")
ax[1].bar(x+.17, q.median_mt_after, width=.34, color="#3B82A6", label="After QC")
ax[1].set_ylabel("Median mitochondrial UMIs (%)"); ax[1].set_title("Strong old/new technical difference"); ax[1].legend(fontsize=10)
ax[2].bar(x, db.predicted_doublet_pct, color=["#C0504D" if v>25 else "#5367B4" for v in db.predicted_doublet_pct])
ax[2].set_ylabel("Classifier-labelled doublets (%)"); ax[2].set_title("Calls are putative; loading unknown")
for a in ax: a.set_xticks(x, ORDER, rotation=55, ha="right"); a.set_xlabel("Library")
fig.tight_layout(); save(fig, "V2_F01_QC")

md = read("cell_metadata.csv")
classes = [k for k in COLORS if k in md.broad_class.unique()]
lims = {}
for nm in ["raw", "source"]:
    xx = md[f"umap_{nm}1"]; yy = md[f"umap_{nm}2"]
    lims[nm] = ((xx.min()-.4, xx.max()+.4), (yy.min()-.4, yy.max()+.4))
def umap(ax, data, nm, field, palette, title, legends=False):
    for k, color in palette.items():
        v = data[data[field] == k]
        if len(v): ax.scatter(v[f"umap_{nm}1"], v[f"umap_{nm}2"], s=1.8, c=color, alpha=.48, rasterized=True, linewidths=0)
    ax.set(xlim=lims[nm][0], ylim=lims[nm][1], xlabel="UMAP 1 (unitless)", ylabel="UMAP 2 (unitless)", title=title)
    if legends:
        handles=[Line2D([], [], color=c, marker="o", linestyle="", markersize=5, label=k.replace("_", " ")) for k,c in palette.items() if k in data[field].values]
        ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1,1), fontsize=10, frameon=False)
fig, ax = plt.subplots(1,2,figsize=(16,6))
umap(ax[0], md, "raw", "broad_class", COLORS, "Uncorrected RNA baseline")
umap(ax[1], md, "source", "broad_class", COLORS, "Source-corrected map (visualization only)", True)
fig.suptitle("E18–P21 broad populations; ambiguous cells retained", fontsize=18)
fig.tight_layout(); save(fig, "V2_F02_UMAP_baseline_corrected")
age_palette = dict(zip(AGES,["#7356A6", "#E0A535", "#087E8B", "#C65050"]))
fig, ax = plt.subplots(1,2,figsize=(15,6))
umap(ax[0], md, "source", "age", age_palette, "Stage on the same map", True)
umap(ax[1], md, "source", "source", {"new_2026":"#087E8B", "old_2023":"#C65050"}, "Source remains partly confounded with age", True)
fig.tight_layout(); save(fig, "V2_F03_UMAP_age_source")
fig, ax = plt.subplots(2,2,figsize=(12,11))
for a, age in zip(ax.flat, AGES):
    v=md[md.age==age]; umap(a,v,"source","broad_class",COLORS,f"{age}: {len(v):,} putative singlets")
handles=[Line2D([],[],color=COLORS[k],marker="o",linestyle="",markersize=5,label=k.replace("_"," ")) for k in classes]
fig.legend(handles=handles,loc="lower center",ncol=5,fontsize=10,frameon=False)
fig.tight_layout(rect=(0,.09,1,1)); save(fig,"V2_F04_UMAP_stages")

co = read("population_composition_per_library.csv")
pv=co[co.selection=="primary_singlets"].pivot(index="sample",columns="population",values="percent").reindex(ORDER).fillna(0)
fig, ax = plt.subplots(figsize=(14,5.5)); bottom=np.zeros(len(ORDER))
for k in classes:
    if k in pv: ax.bar(x,pv[k],bottom=bottom,color=COLORS[k],label=k.replace("_"," ")); bottom+=pv[k].values
ax.set(ylim=(0,100),ylabel="Fraction of retained barcodes (%)",xlabel="Library",title="Captured composition varies strongly by library and QC")
ax.set_xticks(x,ORDER,rotation=45,ha="right"); ax.legend(loc="upper left",bbox_to_anchor=(1,1),fontsize=10,frameon=False)
fig.tight_layout(); save(fig,"V2_F05_population_composition")

e=read("sample_population_expression.csv.gz")
p=e[(e.selection=="primary_singlets") & e.adequate_cells30]
groups=["Glutamatergic", "GABAergic", "ITC-like GABA (provisional)", "Rspo2-positive glut", "Astrocyte", "Microglia"]
def heatmap(genes, name, title, ages=AGES, populations=groups):
    cols=[(pop,age) for pop in populations for age in ages]
    av=p[p.gene.isin(genes)].groupby(["gene","population","age"],observed=True).cpm.mean()
    a=np.array([[np.log2(av.get((g,pop,age),np.nan)+.5) for pop,age in cols] for g in genes])
    fig,ax=plt.subplots(figsize=(max(10,len(cols)*.52), max(5,len(genes)*.33)))
    cm=plt.get_cmap("viridis").copy(); cm.set_bad("#E4E4E4")
    im=ax.imshow(a,aspect="auto",cmap=cm,vmin=-1,vmax=11)
    ax.set_yticks(range(len(genes)),genes); ax.set_xticks(range(len(cols)),[age for _,age in cols])
    for i,pop in enumerate(populations):
        ax.text(i*len(ages)+(len(ages)-1)/2,-1.5,pop.replace(" (provisional)","*").replace("Glutamatergic","Glut").replace("GABAergic","GABA"),ha="center",fontsize=10)
        if i: ax.axvline(i*len(ages)-.5,color="white",lw=2)
    ax.set_title(title,pad=42); ax.set_xlabel("Stage; equal library weight within each group (gray = insufficient cells)")
    fig.colorbar(im,ax=ax,pad=.015,label="log2(mean library CPM + 0.5)")
    fig.tight_layout(); save(fig,name)
heatmap(["Sema3e","Plxnd1","Nrp1","Nrp2","Kirrel3","Sdk2","Cdh8","Cdh9","Cdh13","Slit1","Slit2","Robo1","Robo2","Ntng1","Ntng2","Lrrc4c"],"V2_F06_guidance","Guidance/adhesion candidate expression — not receptor activity")
heatmap(["Mef2c","Mef2a","Mef2d","Arc","Homer1","Pcdh10","Npas4","Bdnf","C1qa","C1qb","C3","Cx3cr1","P2ry12","Mertk","Megf10","Trem2","Tyrobp"],"V2_F07_pruning","Activity/synapse-remodelling and glial machinery — not observed pruning")

fig,axs=plt.subplots(2,3,figsize=(15,8),sharey=True)
for ax,pop in zip(axs.flat,groups):
    v=p[(p.gene=="Mef2c") & (p.population==pop)]
    for src, col, mk in [("new_2026","#087E8B","o"),("old_2023","#C65050","s")]:
        z=v[v.source==src]
        for j,age in enumerate(AGES):
            vals=z[z.age==age].cpm.values
            if not len(vals): continue
            off=np.linspace(-.10,.10,len(vals)) if len(vals)>1 else np.array([0])
            ax.scatter(j+off,np.log2(vals+.5),color=col,marker=mk,s=50,label=src if j==0 else None)
    ax.set_xticks(range(4),AGES); ax.set_title(pop.replace(" (provisional)","*"),fontsize=13)
    ax.set_xlabel("Stage"); ax.set_ylabel("Mef2c log2(CPM + 0.5)")
fig.suptitle("Mef2c: every point is a library, not a cell replicate",fontsize=18)
fig.legend(handles=[Line2D([],[],color="#087E8B",marker="o",linestyle="",label="New source"),Line2D([],[],color="#C65050",marker="s",linestyle="",label="Old source")],loc="lower center",ncol=2)
fig.tight_layout(rect=(0,.04,1,.95)); save(fig,"V2_F08_Mef2c_per_library")

tfs=["Mef2c","Mef2a","Mef2d","Satb1","Satb2","Tbr1","Bcl11b","Nr2f1","Nr2f2","Neurod2","Sox11","Lhx6","Foxp2","Etv1","Sox9","Sox10","Spi1","Irf8"]
v=p[(p.age=="P10") & p.gene.isin(tfs)]
av=v.groupby(["gene","population"],observed=True)[["cpm","detection_depth1000_pct"]].mean()
fig,ax=plt.subplots(figsize=(12,7))
for j,pop in enumerate(groups):
    for i,g in enumerate(tfs):
        if (g,pop) in av.index:
            r=av.loc[(g,pop)]; ax.scatter(j,i,s=12+2.4*r.detection_depth1000_pct,c=np.log2(r.cpm+.5),vmin=-1,vmax=11,cmap="viridis",edgecolors="none")
ax.set_xticks(range(len(groups)),[s.replace(" (provisional)","*").replace("Glutamatergic","Glut").replace("GABAergic","GABA") for s in groups],rotation=25,ha="right")
ax.set_yticks(range(len(tfs)),tfs); ax.invert_yaxis(); ax.set_title("P10 transcription factors: equal library weight; no KO contrast")
fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(-1,11),cmap="viridis"),ax=ax,label="log2(mean library CPM + 0.5)")
ax.legend(handles=[ax.scatter([],[],s=12+2.4*k,color="#777",label=f"{k}%") for k in [10,50,90]],title="Expected detection\nat 1,000 UMIs",loc="upper left",bbox_to_anchor=(1.24,1),frameon=False)
fig.tight_layout(); save(fig,"V2_F09_P10_TFs")

# A shorter, prespecified slide view. The full 18-gene plot and all TFs remain.
short_tfs=["Mef2c","Mef2a","Mef2d","Satb1","Satb2","Tbr1","Bcl11b","Sox11","Lhx6","Foxp2","Etv1","Spi1"]
fig,ax=plt.subplots(figsize=(15,4.7))
for j,pop in enumerate(groups):
    for i,g in enumerate(short_tfs):
        if (g,pop) in av.index:
            r=av.loc[(g,pop)]; ax.scatter(j,i,s=15+3*r.detection_depth1000_pct,c=np.log2(r.cpm+.5),vmin=-1,vmax=11,cmap="viridis",edgecolors="none")
ax.set_xticks(range(len(groups)),["Glut","GABA","ITC-like*","Rspo2+ Glut","Astrocyte","Microglia"],fontsize=16)
ax.set_yticks(range(len(short_tfs)),short_tfs,fontsize=16); ax.invert_yaxis()
ax.set_title("P10 TF expression: prespecified subset; full screen in tables",fontsize=19,pad=14)
cb=fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(-1,11),cmap="viridis"),ax=ax,pad=.03)
cb.set_label("log2(CPM + 0.5)",fontsize=16); cb.ax.tick_params(labelsize=14)
ax.legend(handles=[ax.scatter([],[],s=15+3*k,color="#777",label=f"{k}%") for k in [10,50,90]],title="Expected detection\nat 1,000 UMIs",loc="upper left",bbox_to_anchor=(1.23,1),fontsize=15,title_fontsize=15,frameon=False)
fig.tight_layout(); save(fig,"V2_F09b_P10_TFs_slide")

fx=read("candidate_P10_P21_descriptive_effects.csv")
rob=read("candidate_source_QC_sensitivity.csv")
fig,axs=plt.subplots(1,3,figsize=(16,5.5))
for ax,pop in zip(axs,["Glutamatergic","GABAergic","ITC-like GABA (provisional)"]):
    z=fx[(fx.selection=="primary_singlets") & (fx.population==pop)]
    zz=z.pivot(index="gene",columns="P10_source",values="log2_CPM_ratio_P10_over_P21").dropna()
    if not {"old_2023","new_2026"}.issubset(zz.columns): continue
    rr=rob[rob.population==pop].set_index("gene").reindex(zz.index)
    flag=rr.all_P10_P21_library_pair_directions_agree.fillna(False)
    ax.scatter(zz.old_2023,zz.new_2026,c=np.where(flag,"#087E8B","#BEC4C9"),s=38)
    for g in ["Mef2c","Pcdh10","Sema3e","Kirrel3","Plxnd1","Nrp1","Sdk2"]:
        if g in zz.index: ax.annotate(g,zz.loc[g,["old_2023","new_2026"]],xytext=(3,4),textcoords="offset points",fontsize=10)
    ax.axhline(0,color="#666",lw=.7); ax.axvline(0,color="#666",lw=.7)
    lim=max(1,np.nanmax(np.abs(zz[["old_2023","new_2026"]].values)))+.3
    ax.set(xlim=(-lim,lim),ylim=(-lim,lim),xlabel="Old P10 / P21 log2 CPM ratio",ylabel="New P10 / P21 log2 CPM ratio",title=pop.replace(" (provisional)","*")); ax.plot([-lim,lim],[-lim,lim],"--",color="#AAA",lw=.8)
fig.suptitle("Developmental direction must survive source and QC checks",fontsize=18)
fig.legend(handles=[Line2D([],[],color="#087E8B",marker="o",linestyle="",label="All source / mito / doublet / library-pair directions agree"),Line2D([],[],color="#BEC4C9",marker="o",linestyle="",label="Mixed direction or incomplete support")],loc="lower center",ncol=2,fontsize=10)
fig.tight_layout(rect=(0,.05,1,.94)); save(fig,"V2_F10_source_sensitivity")

v=read("annotation_validation_per_sample.csv").set_index("sample").reindex(ORDER)
fig,ax=plt.subplots(figsize=(13,5))
for col,label,color in [("independent_marker_atlas_agreement","Adult reference vs marker (baseline)","#B9BEC4"),("previous_broad_split_agreement","Previous-reference independent gene splits","#5367B4"),("independent_marker_previous_agreement","Previous reference vs excluded markers","#087E8B")]:
    ax.plot(x,100*v[col],"o-",color=color,label=label)
ax.set(ylim=(0,100),ylabel="Agreement (%)",xlabel="Library",title="Reference checks: disagreement is retained, not hidden")
ax.set_xticks(x,ORDER,rotation=45,ha="right"); ax.legend(fontsize=10)
fig.tight_layout(); save(fig,"V2_F11_annotation_validation")

# Machine-readable chart data for editable PowerPoint tables/charts.
report={"input_barcodes":int(q.input_barcodes.sum()),"qc_barcodes":int(q.qc_primary.sum()),"putative_singlets":int(db.singlets.sum()),"common_genes":25029,
 "libraries":len(ORDER),"stages":AGES,"new_libraries":7,"old_libraries":3,
 "qc":q.reset_index().replace({np.nan:None}).to_dict("records"),"doublets":db.reset_index().to_dict("records"),
 "composition":co[co.selection=="primary_singlets"].to_dict("records"),
 "annotation":v.reset_index().replace({np.nan:None}).to_dict("records"),"embedding":read("embedding_validation.csv").to_dict("records"),
 "core_effects":fx[(fx.selection=="primary_singlets") & fx.gene.isin(["Mef2c","Pcdh10","Kirrel3","Sema3e","Plxnd1","Nrp1"])].replace({np.nan:None}).to_dict("records"),
 "sensitivity":rob.replace({np.nan:None}).to_dict("records"),
 "mef2c":p[p.gene=="Mef2c"].replace({np.nan:None}).to_dict("records"),
 "robust_candidates":rob[rob.all_P10_P21_library_pair_directions_agree & ~rob.annotation_anchor].to_dict("records")}
(T/"presentation_data.json").write_text(json.dumps(report,indent=2,allow_nan=False),encoding="utf-8")
print(json.dumps({k:report[k] for k in ["input_barcodes","qc_barcodes","putative_singlets","common_genes","libraries"]}))
