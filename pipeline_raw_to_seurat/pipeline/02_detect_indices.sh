#!/usr/bin/env bash
# Quick check (~1-2 min per run): convert 2 tiles of lane 1 without
# demultiplexing, count the I1 index read, and compare with samples.csv.
# Writes $WORK_DIR/index_check/detected_<run>.csv, which step 03 uses to skip
# indices that are not actually in a run.
set -euo pipefail
source "$(dirname "$0")/config.sh"
TOOL=$(pick_demux_tool)
D="$WORK_DIR/index_check"; mkdir -p "$D"
TILES='s_1_1110[12]'

for run in $(runs_in_sheet); do
  rf=$(run_folder "$run"); [[ -n "$rf" ]] || { echo "[detect] $run not extracted"; exit 1; }
  out="$D/$run"; rm -rf "$out"; mkdir -p "$out"
  echo "[detect] $run ($TOOL, tiles $TILES)"
  if [[ "$TOOL" == bcl2fastq ]]; then
    printf '[Data]\nSample_ID,Sample_Name\nALL,ALL\n' > "$out/SampleSheet.csv"
    bcl2fastq --runfolder-dir "$rf" --output-dir "$out/fastq" \
      --sample-sheet "$out/SampleSheet.csv" --use-bases-mask "$BASES_MASK_BCL2FASTQ" \
      --tiles "$TILES" --create-fastq-for-index-reads \
      --minimum-trimmed-read-length 8 --mask-short-adapter-reads 8 \
      --ignore-missing-positions --ignore-missing-controls --ignore-missing-filter --ignore-missing-bcls \
      -r 4 -p "$THREADS" -w 4 > "$out/bcl2fastq.log" 2>&1
  else
    cat > "$out/SampleSheet.csv" <<EOF
[Header]
FileFormatVersion,2
[BCLConvert_Settings]
OverrideCycles,$OVERRIDE_CYCLES_BCLCONVERT
CreateFastqForIndexReads,1
[BCLConvert_Data]
Sample_ID
ALL
EOF
    bcl-convert --bcl-input-directory "$rf" --output-directory "$out/fastq" \
      --sample-sheet "$out/SampleSheet.csv" --tiles "$TILES" --force > "$out/bcl-convert.log" 2>&1
  fi
  i1=$(find "$out/fastq" -name 'ALL_*_I1_001.fastq.gz' | head -1)
  [[ -n "$i1" ]] || { echo "[detect] no I1 FASTQ produced; see $out/*.log"; exit 1; }
  python3 "$PIPE_DIR/detect_10x_index.py" --fastq "$i1" --table "$SI_TABLE" \
    --samples "$SAMPLES_CSV" --run "$run" --out "$D/detected_$run.csv" \
    --min-fraction "$MIN_INDEX_FRACTION"
  rm -rf "$out/fastq"
done
echo "[detect] tables in $D/detected_*.csv"
