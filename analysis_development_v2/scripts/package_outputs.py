"""Package generated outputs; no input modification and no GitHub upload."""
from pathlib import Path
import gzip, hashlib, json, os, shutil, zipfile
ROOT=Path(os.environ.get("BLA_V2_ROOT","analysis_development_v2")).resolve()
release=ROOT/"release"; release.mkdir(exist_ok=True)
private=release/"private_data"
raw=private/"BLA_E18_P0_P10_P21_seurat_v2.uncompressed.rds"
target=private/"BLA_E18_P0_P10_P21_seurat_v2.rds"
if raw.exists():
    if not target.exists() or target.stat().st_mtime < raw.stat().st_mtime:
        with raw.open("rb") as source, target.open("wb") as dest:
            with gzip.GzipFile(filename="",mode="wb",fileobj=dest,compresslevel=3,mtime=0) as zipped:
                shutil.copyfileobj(source,zipped,8*1024*1024)
    # Independent compressed-stream and byte-level comparison to valid raw RDS.
    ha=hashlib.sha256(); hb=hashlib.sha256()
    with raw.open("rb") as f:
        for b in iter(lambda:f.read(8*1024*1024),b""): ha.update(b)
    with gzip.open(target,"rb") as f:
        for b in iter(lambda:f.read(8*1024*1024),b""): hb.update(b)
    if ha.digest()!=hb.digest(): raise RuntimeError("Compressed export differs from validated uncompressed RDS")
    print(json.dumps({"export":str(target),"bytes":target.stat().st_size,"compressed_stream_valid":True,"uncompressed_bytes_identical":True}),flush=True)

resources=[p for p in (ROOT/"resources").iterdir() if p.is_file()]
with zipfile.ZipFile(release/"BLA_v2_reference_resources.zip","w",zipfile.ZIP_DEFLATED,compresslevel=5) as z:
    for p in sorted(resources): z.write(p,p.name)
public_results=[p for p in (ROOT/"results").rglob("*") if p.is_file()]
cell_level={"all_input_barcodes_qc.csv","doublet_calls.csv","cell_metadata.csv","cell_marker_scores.csv","previous_label_centroids.csv"}
public_results=[p for p in public_results if p.name not in cell_level and not p.name.startswith("all_gene_pseudobulk")]
with zipfile.ZipFile(release/"BLA_v2_results_and_pipeline.zip","w",zipfile.ZIP_DEFLATED,compresslevel=5) as z:
    for p in sorted(public_results): z.write(p,str(p.relative_to(ROOT)))
    for p in sorted((ROOT/"scripts").iterdir()):
        if p.is_file() and p.suffix in {".R",".py",".mjs",".ps1"}: z.write(p,"scripts/"+p.name)
    for name in ["README.md","METHODS.md","ANALYSIS_PLAN.md","RELEASE_NOTES.md"]: z.write(ROOT/name,name)
assets=[release/"BLA_v2_reference_resources.zip",release/"BLA_v2_results_and_pipeline.zip",ROOT/"results/BLA_development_E18_P0_P10_P21_2026-10-06_final.pptx",ROOT/"results/BLA_development_E18_P0_P10_P21_2026-10-06_final.pdf"]
if private.exists(): assets+=sorted(p for p in private.iterdir() if p.is_file() and ".uncompressed." not in p.name and p.name!="CHECKSUMS.sha256")
manifest=[]
for p in assets:
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(8*1024*1024),b""): h.update(b)
    manifest.append({"file":p.name,"bytes":p.stat().st_size,"sha256":h.hexdigest(),"public_release_eligible":p.parent!=private})
(release/"manifest.json").write_text(json.dumps(manifest,indent=2),encoding="utf-8")
(release/"CHECKSUMS.sha256").write_text("".join(f"{x['sha256']}  {x['file']}\n" for x in manifest if x["public_release_eligible"]),encoding="utf-8")
if private.exists():
    (private/"CHECKSUMS.sha256").write_text("".join(f"{x['sha256']}  {x['file']}\n" for x in manifest if not x["public_release_eligible"]),encoding="utf-8")
print(json.dumps({"public_assets":[x for x in manifest if x["public_release_eligible"]],"private_assets_prepared":sum(not x["public_release_eligible"] for x in manifest)}),flush=True)
