#!/usr/bin/env bash
# Install bcl2fastq2 (via micromamba), Cell Ranger, and the 10x mouse reference
# into $HOME. Linux / WSL Ubuntu only. Safe to re-run.
#
#   bash 00_install_tools.sh /path/to/cellranger-9.x.x.tar.gz
#
# Cell Ranger tarball: https://www.10xgenomics.com/support/software/cell-ranger/downloads
# (10x asks for name/email before giving the download link.)
set -euo pipefail
source "$(dirname "$0")/config.sh"

CR_TAR="${1:-}"
mkdir -p "$HOME/bin" "$HOME/refs"

# ---- bcl2fastq2 -------------------------------------------------------------
if ! command -v bcl2fastq >/dev/null 2>&1; then
  if [[ ! -x "$HOME/bin/micromamba" ]]; then
    echo "[install] micromamba"
    curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj -C "$HOME" bin/micromamba
  fi
  export MAMBA_ROOT_PREFIX="$HOME/micromamba"
  echo "[install] bcl2fastq2 -> env 'demux'"
  "$HOME/bin/micromamba" create -y -n demux -c bih-cubi -c conda-forge bcl2fastq2 python=3.11 \
    || "$HOME/bin/micromamba" create -y -n demux -c dranew -c conda-forge bcl2fastq python=3.11
  ln -sf "$MAMBA_ROOT_PREFIX/envs/demux/bin/bcl2fastq" "$HOME/bin/bcl2fastq"
fi

# ---- Cell Ranger ------------------------------------------------------------
if ! command -v cellranger >/dev/null 2>&1; then
  if [[ -z "$CR_TAR" ]]; then
    echo "[install] cellranger not found. Download the tarball from 10x and re-run:"
    echo "          bash $0 /path/to/cellranger-x.y.z.tar.gz"; exit 1
  fi
  tar -xf "$CR_TAR" -C "$HOME"   # .tar or .tar.gz
  CR_BIN=$(ls -d "$HOME"/cellranger-*/bin 2>/dev/null | sort -V | tail -1)
  ln -sf "$CR_BIN/cellranger" "$HOME/bin/cellranger"
fi

# ---- Reference --------------------------------------------------------------
if [[ ! -d "$REF_DIR" ]]; then
  echo "[install] reference -> $REF_DIR (~11 GB download)"
  curl -L -o "$HOME/refs/ref.tar.gz" "$REF_URL"
  tar -xzf "$HOME/refs/ref.tar.gz" -C "$(dirname "$REF_DIR")"
  rm -f "$HOME/refs/ref.tar.gz"
fi

grep -q 'HOME/bin' "$HOME/.bashrc" || echo 'export PATH="$HOME/bin:$PATH"' >> "$HOME/.bashrc"
export PATH="$HOME/bin:$PATH"
echo "[install] done:"
bcl2fastq --version 2>&1 | head -2 || true
cellranger --version || true
ls "$REF_DIR"
