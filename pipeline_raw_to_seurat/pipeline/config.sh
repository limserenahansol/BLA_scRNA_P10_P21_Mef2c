# ===== E18 / P0 scRNA-seq (10x 3' v3, NextSeq 500) — pipeline config =====
# Source this file; every pipeline step reads it. Edit paths here only.
#
# Runs (from RunInfo.xml): R1=28, I1=8 (single index, SI-GA plate), R2=130
#   P165_191219_NB500982_0166_AH2GVNBGXF   high output, P0_1 + P0_2
#   P188_200212_NB500982_0189_AH5FVJBGXF   high output, P0_3 + E18_1 + E18_2
#   P189_200227_NB500982_0195_AH2NV5AFX2   mid output (top-up), same libraries as P188?
#                                          -> 02_detect_indices.sh verifies this

PIPE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJ_DIR="$(dirname "$PIPE_DIR")"

# Folder that holds the downloaded P165_*.tar.gz / P188_*.tar.gz / P189_*.tar.gz
# (WSL sees Windows drives as /mnt/<letter>/...)
TAR_DIR="${TAR_DIR:-/mnt/c/Users/hsollim/Downloads}"

# Scratch space: extracted BCL run folders + FASTQs + cellranger runs.
# Needs ~3x the tarball size. On WSL, a Linux-native path ($HOME/...) is much
# faster than /mnt/<drive> for cellranger; outputs get copied back in step 05.
WORK_DIR="${WORK_DIR:-$HOME/scRNA_E18_P0_work}"

# Final, shareable results (cellranger outs + Seurat object + figures)
OUT_DIR="${OUT_DIR:-$PROJ_DIR/results}"

SAMPLES_CSV="$PROJ_DIR/samples.csv"
SI_TABLE="$PROJ_DIR/resources/SI-GA_single_index.csv"

# Tools: cellranger is downloaded from 10x (free form); bcl2fastq2 or bcl-convert
# from Illumina/conda. Leave blank to use whatever is on PATH.
CELLRANGER="${CELLRANGER:-cellranger}"
DEMUX_TOOL="${DEMUX_TOOL:-auto}"          # auto | bcl2fastq | bcl-convert

# Reference: mm10-2020-A, same as the earlier P10/P21/P56 BLA data, so the
# datasets can be merged gene-for-gene. (Current 10x mouse = GRCm39-2024-A.)
REF_DIR="${REF_DIR:-$HOME/refs/refdata-gex-mm10-2020-A}"
REF_URL="https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz"

# Resources (this PC: 32 cores / 128 GB; leave headroom for Windows)
THREADS="${THREADS:-28}"
MEM_GB="${MEM_GB:-100}"

# Read structure
BASES_MASK_BCL2FASTQ="Y28,I8,Y130"
OVERRIDE_CYCLES_BCLCONVERT="Y28;I8;Y130"

# An index must carry at least this fraction of a run's reads to be demultiplexed
MIN_INDEX_FRACTION=0.005

# Delete extracted BCL folders after FASTQs are made (tarballs are never touched)
CLEAN_BCL="${CLEAN_BCL:-0}"

mkdir -p "$WORK_DIR" "$OUT_DIR"

pick_demux_tool() {
  if [[ "$DEMUX_TOOL" != "auto" ]]; then echo "$DEMUX_TOOL"; return; fi
  if command -v bcl2fastq >/dev/null 2>&1; then echo bcl2fastq
  elif command -v bcl-convert >/dev/null 2>&1; then echo bcl-convert
  else echo "ERROR: neither bcl2fastq nor bcl-convert on PATH (run 00_install_tools.sh)" >&2; return 1
  fi
}

run_folder() {  # run_folder P165 -> full path of extracted run folder
  ls -d "$WORK_DIR/bcl/$1"_*/ 2>/dev/null | head -1 | sed 's:/$::'
}

runs_in_sheet() {
  tail -n +2 "$SAMPLES_CSV" | cut -d, -f2 | tr -d '\r' | sort -u
}
