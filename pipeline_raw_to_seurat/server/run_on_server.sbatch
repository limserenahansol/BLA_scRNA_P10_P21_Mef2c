#!/usr/bin/env bash
# Run the whole pipeline on a SLURM cluster (e.g. Stanford Sherlock): download
# the three run tarballs straight to the cluster, then BCL -> FASTQ ->
# cellranger -> Seurat. Only results/ (a few GB) needs to come back.
#
#   1. Copy this project folder to the cluster:
#        scp -r E18_P0_2026 <sunetid>@login.sherlock.stanford.edu:~/
#   2. Put the three download links (one per line) in server/urls.txt
#   3. sbatch ~/E18_P0_2026/server/run_on_server.sbatch
#
#SBATCH --job-name=e18p0_scrna
#SBATCH --partition=normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=2-00:00:00
#SBATCH --output=%x_%j.log
set -euo pipefail

PROJ="${PROJ:-$HOME/E18_P0_2026}"
export WORK_DIR="${SCRATCH:-$HOME}/scRNA_E18_P0_work"
export TAR_DIR="$WORK_DIR/tars"
export OUT_DIR="$PROJ/results"
export THREADS="${SLURM_CPUS_PER_TASK:-16}"
export MEM_GB=$(( ${SLURM_MEM_PER_NODE:-98304} / 1024 - 6 ))
export PATH="$HOME/bin:$PATH"
mkdir -p "$TAR_DIR"

# ---- Tools: cluster modules if present, else the user-space installer --------
module load biology 2>/dev/null || true
module load cellranger 2>/dev/null || true
module load bcl2fastq 2>/dev/null || module load bcl2fastq2 2>/dev/null || true
module load R 2>/dev/null || true
if ! command -v cellranger >/dev/null || ! command -v bcl2fastq >/dev/null || [[ ! -d "${REF_DIR:-$HOME/refs/refdata-gex-mm10-2020-A}" ]]; then
  CR_TAR=$(ls "$HOME"/cellranger-*.tar* 2>/dev/null | head -1 || true)
  bash "$PROJ/pipeline/00_install_tools.sh" ${CR_TAR:+"$CR_TAR"}
fi

# ---- Download (resumable; re-submitting continues partial files) -------------
URLS="$PROJ/server/urls.txt"
if [[ -s "$URLS" ]]; then
  while read -r url; do
    [[ -z "$url" || "$url" == \#* ]] && continue
    echo "[download] $url"
    ( cd "$TAR_DIR" && wget -c --content-disposition --tries=20 --timeout=60 -q --show-progress "$url" )
  done < "$URLS"
fi
ls -lh "$TAR_DIR"

# ---- Pipeline -----------------------------------------------------------------
bash "$PROJ/pipeline/run_all.sh"
echo "[done] copy back: $OUT_DIR"
