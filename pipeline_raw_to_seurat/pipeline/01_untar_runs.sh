#!/usr/bin/env bash
# Extract the three NextSeq run tarballs into $WORK_DIR/bcl.
# Focus/thumbnail images are skipped (not needed for demultiplexing).
set -euo pipefail
source "$(dirname "$0")/config.sh"
mkdir -p "$WORK_DIR/bcl"

for run in $(runs_in_sheet); do
  if [[ -n "$(run_folder "$run")" && -f "$(run_folder "$run")/RTAComplete.txt" && -f "$(run_folder "$run")/.untar_done" ]]; then
    echo "[untar] $run already extracted"; continue
  fi
  tarball=$(ls "$TAR_DIR"/${run}_*.tar.gz "$TAR_DIR"/${run}_*.tgz "$TAR_DIR"/${run}*.tar 2>/dev/null | head -1 || true)
  if [[ -z "$tarball" ]]; then
    echo "[untar] ERROR: no tarball for $run in $TAR_DIR (expected ${run}_*.tar.gz)"; exit 1
  fi
  echo "[untar] $run <- $tarball ($(du -h "$tarball" | cut -f1))"
  # tar exits non-zero on a truncated download, so no separate gzip -t pass
  tar -xf "$tarball" -C "$WORK_DIR/bcl" \
      --exclude='*/Images/*' --exclude='*/Thumbnail_Images/*' \
    || { echo "[untar] ERROR: $tarball is truncated/corrupt (download not finished?)"; exit 1; }
  rf=$(run_folder "$run")
  [[ -f "$rf/RTAComplete.txt" ]] || { echo "[untar] ERROR: $rf has no RTAComplete.txt"; exit 1; }
  touch "$rf/.untar_done"
  if [[ -f "$rf/SampleSheet.csv" ]]; then
    echo "[untar] note: $run ships its own SampleSheet.csv (ignored; samples.csv is used):"
    sed -n '/\[Data\]/,$p' "$rf/SampleSheet.csv" | head -20
  fi
done
echo "[untar] done"; du -sh "$WORK_DIR"/bcl/*
