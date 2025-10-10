#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"   # tests
DATA_DIR="$ROOT/data"
FRAG_GZ="$DATA_DIR/toy_fragments.tsv.gz"
CPEAKS_GZ="$DATA_DIR/cpeaks_hg38.bed.gz"
OUTDIR="$ROOT/output"

echo "ROOT: $ROOT"
echo "FRAG: $FRAG_GZ"
echo "CPEAKS: $CPEAKS_GZ"
echo "OUTDIR: $OUTDIR"

# cleanup + prepare
rm -rf "$OUTDIR"
mkdir -p "$OUTDIR"

ln -sf "$CPEAKS_GZ" ./cpeaks_hg38.bed.gz

# Run with normal mode to reduce memory usage
python main.py --fragment_path "$FRAG_GZ" --output "$OUTDIR" --mode normal

# Basic assertions: default output_name is 'cell-cpeaks' -> check that file
MTX_FILE="$OUTDIR/cell-cpeaks.mtx"
if [[ ! -f "$MTX_FILE" ]]; then
  echo "Smoke test failed: expected output file not found: $MTX_FILE" >&2
  ls -R "$OUTDIR" || true
  exit 2
fi
if [[ ! -s "$MTX_FILE" ]]; then
  echo "Smoke test failed: output file is empty: $MTX_FILE" >&2
  exit 3
fi

echo "Smoke test passed. Output: $MTX_FILE"
