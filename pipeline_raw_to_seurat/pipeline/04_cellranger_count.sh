#!/usr/bin/env bash
# cellranger count per sample. A sample sequenced on two runs (e.g. P188 high
# output + P189 mid-output top-up) gets both FASTQ folders in one --fastqs list,
# so its reads are pooled into one library.
# Intronic reads are counted (Cell Ranger >= 7 default) - right for both cells
# and nuclei. No BAM is written (saves ~30 GB per sample).
set -euo pipefail
source "$(dirname "$0")/config.sh"
command -v "$CELLRANGER" >/dev/null || { echo "cellranger not on PATH"; exit 1; }
[[ -d "$REF_DIR" ]] || { echo "reference not found: $REF_DIR"; exit 1; }
mkdir -p "$WORK_DIR/cellranger"; cd "$WORK_DIR/cellranger"

SAMPLES=$(tail -n +2 "$SAMPLES_CSV" | cut -d, -f1 | sort -u)
for sid in $SAMPLES; do
  if [[ -f "$sid/outs/filtered_feature_bc_matrix.h5" ]]; then echo "[count] $sid done already"; continue; fi
  dirs=$(find "$WORK_DIR/fastq" -name "${sid}_S*_R1_001.fastq.gz" -printf '%h\n' | sort -u | paste -sd, -)
  [[ -n "$dirs" ]] || { echo "[count] WARNING: no FASTQs for $sid - skipped"; continue; }
  echo "[count] $sid  <- $dirs"
  rm -rf "$sid"
  "$CELLRANGER" count --id="$sid" --sample="$sid" --fastqs="$dirs" \
    --transcriptome="$REF_DIR" --chemistry=auto --create-bam=false \
    --localcores="$THREADS" --localmem="$MEM_GB" > "$sid.log" 2>&1 \
    || { echo "[count] $sid failed; tail of log:"; tail -30 "$sid.log"; exit 1; }
done

# ---- Collect the small, shareable outputs ----------------------------------
dest="$OUT_DIR/cellranger"; mkdir -p "$dest"
python3 - "$WORK_DIR/cellranger" "$dest" $SAMPLES <<'EOF'
import csv, shutil, sys, os
src, dest, samples = sys.argv[1], sys.argv[2], sys.argv[3:]
keep = ["filtered_feature_bc_matrix.h5", "raw_feature_bc_matrix.h5",
        "metrics_summary.csv", "web_summary.html", "molecule_info.h5"]
rows = []
for s in samples:
    o = os.path.join(src, s, "outs")
    if not os.path.isdir(o):
        continue
    d = os.path.join(dest, s); os.makedirs(d, exist_ok=True)
    for f in keep:
        if os.path.exists(os.path.join(o, f)):
            shutil.copy2(os.path.join(o, f), d)
    for sub in ("filtered_feature_bc_matrix", "raw_feature_bc_matrix"):
        if os.path.isdir(os.path.join(o, sub)) and not os.path.isdir(os.path.join(d, sub)):
            shutil.copytree(os.path.join(o, sub), os.path.join(d, sub))
    with open(os.path.join(o, "metrics_summary.csv")) as fh:
        m = next(csv.DictReader(fh)); m = {"sample": s, **m}; rows.append(m)
if rows:
    with open(os.path.join(dest, "all_samples_metrics_summary.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    for m in rows:
        print(f"[count] {m['sample']:<8} cells={m.get('Estimated Number of Cells','?'):>8}  "
              f"reads/cell={m.get('Mean Reads per Cell','?'):>8}  genes/cell={m.get('Median Genes per Cell','?'):>6}  "
              f"saturation={m.get('Sequencing Saturation','?')}")
EOF
echo "[count] shareable outputs -> $dest"
