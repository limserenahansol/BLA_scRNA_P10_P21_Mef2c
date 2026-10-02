#!/usr/bin/env bash
# BCL -> FASTQ per run, using the SI-GA wells in samples.csv (each well = 4
# oligos, all written to the same Sample_ID). Indices that step 02 found absent
# from a run are skipped. Output: $WORK_DIR/fastq/<run>/<sample>_S*_L00*_{R1,R2,I1}_001.fastq.gz
set -euo pipefail
source "$(dirname "$0")/config.sh"
TOOL=$(pick_demux_tool)
mkdir -p "$WORK_DIR/fastq"

oligos() {  # oligos SI-GA-E6 -> 4 lines
  awk -F, -v w="$1" '$1==w {print $2; print $3; print $4; print $5}' "$SI_TABLE"
}

for run in $(runs_in_sheet); do
  out="$WORK_DIR/fastq/$run"
  if [[ -f "$out/.demux_done" ]]; then echo "[demux] $run done already"; continue; fi
  rf=$(run_folder "$run"); [[ -n "$rf" ]] || { echo "[demux] $run not extracted"; exit 1; }
  det="$WORK_DIR/index_check/detected_$run.csv"
  [[ -f "$det" ]] || echo "[demux] WARNING: no $det (step 02 skipped) - using all samples.csv rows"

  rows=()
  while IFS=, read -r sid r stage idx; do
    idx=$(echo "$idx" | tr -d '\r'); [[ "$r" == "$run" ]] || continue
    if [[ -f "$det" ]] && ! grep -q "^$run,$idx,.*,present$" "$det"; then
      echo "[demux] $run: skipping $sid ($idx) - not present in this run"; continue
    fi
    [[ -n "$(oligos "$idx")" ]] || { echo "[demux] ERROR: $idx not in $SI_TABLE"; exit 1; }
    rows+=("$sid,$idx")
  done < <(tail -n +2 "$SAMPLES_CSV")
  [[ ${#rows[@]} -gt 0 ]] || { echo "[demux] $run: no samples to demultiplex"; continue; }

  rm -rf "$out"; mkdir -p "$out"
  sheet="$out/SampleSheet_$run.csv"
  if [[ "$TOOL" == bcl2fastq ]]; then
    { echo "[Data]"; echo "Sample_ID,Sample_Name,index"
      for x in "${rows[@]}"; do sid=${x%%,*}; for o in $(oligos "${x##*,}"); do echo "$sid,$sid,$o"; done; done
    } > "$sheet"
    echo "[demux] $run: bcl2fastq, ${#rows[@]} samples"
    bcl2fastq --runfolder-dir "$rf" --output-dir "$out" --sample-sheet "$sheet" \
      --use-bases-mask "$BASES_MASK_BCL2FASTQ" --create-fastq-for-index-reads \
      --minimum-trimmed-read-length 8 --mask-short-adapter-reads 8 \
      --ignore-missing-positions --ignore-missing-controls --ignore-missing-filter --ignore-missing-bcls \
      -r 6 -p "$THREADS" -w 6 > "$out/bcl2fastq.log" 2>&1 \
      || { echo "[demux] bcl2fastq failed; tail of log:"; tail -20 "$out/bcl2fastq.log"; exit 1; }
  else
    { echo "[Header]"; echo "FileFormatVersion,2"
      echo "[BCLConvert_Settings]"; echo "OverrideCycles,$OVERRIDE_CYCLES_BCLCONVERT"
      echo "CreateFastqForIndexReads,1"; echo "MinimumTrimmedReadLength,8"; echo "MaskShortReads,8"
      echo "[BCLConvert_Data]"; echo "Sample_ID,Index"
      for x in "${rows[@]}"; do sid=${x%%,*}; for o in $(oligos "${x##*,}"); do echo "$sid,$o"; done; done
    } > "$sheet"
    echo "[demux] $run: bcl-convert, ${#rows[@]} samples"
    bcl-convert --bcl-input-directory "$rf" --output-directory "$out" --sample-sheet "$sheet" \
      --force > "$out/bcl-convert.log" 2>&1 \
      || { echo "[demux] bcl-convert failed; tail of log:"; tail -20 "$out/bcl-convert.log"; exit 1; }
  fi

  # Per-sample read counts from the converter's own report
  stats=$(find "$out" -path '*Stats/Stats.json' -o -path '*Reports/Demultiplex_Stats.csv' | head -1)
  [[ -n "$stats" ]] && cp "$stats" "$out/demux_stats_$run.${stats##*.}"
  if [[ "$stats" == *.json ]]; then
    python3 - "$stats" <<'EOF'
import json, sys, collections
d = json.load(open(sys.argv[1])); c = collections.Counter(); tot = 0
for lane in d["ConversionResults"]:
    tot += lane["TotalClustersPF"]
    for s in lane["DemuxResults"]:
        c[s["SampleId"]] += s["NumberReads"]
for s, n in c.items():
    print(f"[demux] {s}: {n:,} reads ({n / tot:.1%} of PF clusters)")
print(f"[demux] undetermined: {(tot - sum(c.values())) / tot:.1%}")
EOF
  elif [[ -n "$stats" ]]; then
    column -s, -t "$stats" | head -20
  fi
  touch "$out/.demux_done"
  if [[ "$CLEAN_BCL" == 1 ]]; then rm -rf "$rf"; echo "[demux] removed $rf"; fi
done
