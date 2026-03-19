#!/usr/bin/env bash
set -euo pipefail

# Simple installer for Allele_auto
# Usage: ./install.sh [--no-path] [--env-name NAME]

ENV_NAME="allele_env"
ADD_PATH=1

while [[ $# -gt 0 ]]; do
  case "$1" in
    --no-path) ADD_PATH=0; shift ;;
    --env-name) ENV_NAME="$2"; shift 2 ;;
    -n) ENV_NAME="$2"; shift 2 ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

HERE=$(cd "$(dirname "$0")" && pwd)

echo "Installer running in: $HERE"

if ! command -v conda >/dev/null 2>&1; then
  echo "ERROR: conda not found. Please install Miniconda/Anaconda first." >&2
  exit 1
fi

echo "Initializing conda shell hook..."
eval "$(conda shell.bash hook)"

# pick environment.yml: prefer root, fallback to bin/
ENV_FILE="$HERE/environment.yml"
if [ ! -f "$ENV_FILE" ]; then
  ENV_FILE="$HERE/bin/environment.yml"
fi

if [ ! -f "$ENV_FILE" ]; then
  echo "ERROR: environment.yml not found in repo root or bin/." >&2
  exit 1
fi

echo "Creating/updating conda environment '$ENV_NAME' from $ENV_FILE ..."
if conda env list | awk '{print $1}' | grep -q "^${ENV_NAME}$"; then
  conda env update -f "$ENV_FILE" -n "$ENV_NAME"
else
  conda env create -f "$ENV_FILE" -n "$ENV_NAME"
fi

echo "Marking main script executable..."
if [ -f "$HERE/bin/allele_identification.sh" ]; then
  chmod +x "$HERE/bin/allele_identification.sh"
fi

if [ $ADD_PATH -eq 1 ]; then
  SHELL_RC="$HOME/.bash_profile"
  if [ -f "$HOME/.bashrc" ] && [ -z "$(grep -F "$HERE/bin" "$HOME/.bash_profile" 2>/dev/null || true)" ]; then
    SHELL_RC="$HOME/.bashrc"
  fi
  echo "Adding $HERE/bin to PATH in $SHELL_RC (if not already present)..."
  grep -F "$HERE/bin" "$SHELL_RC" >/dev/null 2>&1 || echo "export PATH=\"$HERE/bin:\$PATH\"" >> "$SHELL_RC"
  echo "To use immediately: source $SHELL_RC"
else
  echo "Skipping PATH modification (run scripts via $HERE/bin/)"
fi

echo "Installation complete. Activate environment with: conda activate $ENV_NAME"
echo "Example run: bash $HERE/bin/allele_identification.sh -p 00_data/chrpairs.txt -a SA -b SB"
