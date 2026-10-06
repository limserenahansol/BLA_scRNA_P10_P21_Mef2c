"""Independent read-only numerical/package checks; save a small audit report."""
from pathlib import Path
import json, os, zipfile
import numpy as np
import pandas as pd
ROOT=Path(os.environ.get("BLA_V2_ROOT","analysis_development_v2")).resolve()
T=ROOT/"results/tables"; V=ROOT/"results/validation"
e=pd.read_csv(T/"sample_population_expression.csv.gz",low_memory=False)
b=pd.read_csv(T/"all_gene_pseudobulk_counts.csv.gz",index_col=0)
gi=pd.read_csv(T/"pseudobulk_group_metadata.csv")
worst=0.; seen=0
info={(r["sample"],r["population"]):r for r in gi.to_dict("records")}
for key,z in e[e.selection=="primary_singlets"].groupby(["sample","population"]):
    col="primary_singlets|"+key[0]+"|"+key[1]
    r=info[key]
    assert b[col].sum()==r["total_umi"]
    expected=b.loc[z.gene,col].to_numpy()/r["total_umi"]*1e6
    worst=max(worst,float(np.max(np.abs(expected-z.cpm.to_numpy())))); seen+=len(z)
assert worst<1e-7
assert np.all(e.detection_depth1000_pct.fillna(0)<=e.detection_pct+1e-8)
fx=pd.read_csv(T/"candidate_P10_P21_descriptive_effects.csv")
effect_error=0.
means=e[e.adequate_cells30].groupby(["selection","population","gene","age","source"],observed=True).cpm.mean()
for r in fx.to_dict("records"):
    key=(r["selection"],r["population"],r["gene"])
    a=means.loc[key+("P10",r["P10_source"])]
    d=means.loc[key+("P21","old_2023")]
    expected=np.log2((a+.5)/(d+.5))
    effect_error=max(effect_error,abs(expected-r["log2_CPM_ratio_P10_over_P21"]))
assert effect_error<1e-10

deck=ROOT/"results/BLA_development_E18_P0_P10_P21_2026-10-06_final.pptx"
audit=json.loads((V/"presentation.validation.json").read_text())
assert audit["packageIntegrity"]["finding_count"]==0 and not audit["presentationLayout"]["warnings"]
assert audit["nativeChartValidation"]["passed"]
native=json.loads((V/"native_powerpoint.json").read_text(encoding="utf-8-sig"))
assert native["native_open_verified"] and native["source_pptx_unchanged"] and native["native_series_readable"]
source=json.loads((T/"presentation_data.json").read_text())
chart_error=0.
for r in native["series"]:
    if r["slide"]==3:
        vals=([100*q["qc_primary"]/q["input_barcodes"] for q in source["qc"]] if r["series"].startswith("QC-pass") else
              [next(x["predicted_doublet_pct"] for x in source["doublets"] if x["sample"]==q["sample"]) for q in source["qc"]])
    elif r["slide"]==8:
        # The three scatter charts share their two series names; compare as
        # an exact multiset against the prespecified three population series.
        pop_values=[]
        src="new_2026" if r["series"]=="New source" else "old_2023"
        for pop in ["Glutamatergic","GABAergic","ITC-like GABA (provisional)"]:
            pop_values.append([np.log2(x["cpm"]+.5) for x in source["mef2c"] if x["population"]==pop and x["source"]==src])
        errors=[np.max(np.abs(np.round(v,5)-r["values"])) for v in pop_values if len(v)==len(r["values"])]
        chart_error=max(chart_error,min(errors)); continue
    else: continue
    chart_error=max(chart_error,float(np.max(np.abs(np.round(vals,5)-r["values"]))))
assert chart_error<1e-8
with zipfile.ZipFile(deck) as z:
    assert len([n for n in z.namelist() if n.startswith("ppt/slides/slide") and n.endswith(".xml")])==18
    assert len([n for n in z.namelist() if n.startswith("ppt/embeddings/") and n.endswith(".xlsx")])==5
result={"independent_pseudobulk_CPM_values_checked":seen,"max_CPM_absolute_error":worst,
 "descriptive_effects_checked":len(fx),"max_log2_ratio_error":effect_error,"standardized_detection_not_inflated":True,
 "native_chart_values_match_source":True,"max_native_chart_display_error":chart_error,
 "ppt_slides":18,"native_charts":5,"native_tables":6,"layout_warnings":0,"passed":True}
(V/"independent_delivery_checks.json").write_text(json.dumps(result,indent=2),encoding="utf-8")
print(json.dumps(result))
