#!/usr/bin/env bash
# Steps 01-04 (Linux / WSL). Each step skips work that is already done, so this
# can simply be re-run after a failure. Step 05 (Seurat) runs in R afterwards.
#   bash run_all.sh 2>&1 | tee run_all.log
set -euo pipefail
P="$(cd "$(dirname "$0")" && pwd)"
export PATH="$HOME/bin:$PATH"
bash "$P/01_untar_runs.sh"
bash "$P/02_detect_indices.sh"
bash "$P/03_demultiplex.sh"
bash "$P/04_cellranger_count.sh"
source "$P/config.sh"
if command -v Rscript >/dev/null 2>&1; then
  Rscript "$P/05_seurat_preprocess.R" --cr_dir="$OUT_DIR/cellranger" --out_dir="$OUT_DIR/seurat"
else
  echo "Cell Ranger done. Now run step 05 in R (Windows):"
  echo "  Rscript $P/05_seurat_preprocess.R"
fi
