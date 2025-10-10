#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
FRAG_GZ="$ROOT/data/toy_fragments.tsv.gz"
CPEAKS_GZ="$ROOT/data/cpeaks_hg38.bed.gz"
OUTDIR="$ROOT/output"

echo "ROOT: $ROOT"
echo "FRAG: $FRAG_GZ"
echo "CPEAKS: $CPEAKS_GZ"
echo "OUTDIR: $OUTDIR"

# Ensure output dir
rm -rf "$OUTDIR"
mkdir -p "$OUTDIR"

# Place cpeaks file at repo root because main.py expects cpeaks_hg38.bed.gz at repo root
cp "$CPEAKS_GZ" "$ROOT/../cpeaks_hg38.bed.gz" 2>/dev/null || cp "$CPEAKS_GZ" "./cpeaks_hg38.bed.gz"

# Run the main script (adjust python path if needed)
python main.py --fragment_path "$FRAG_GZ" --output "$OUTDIR"

# Basic assertions: expect output mtx file (main.py writes <output_name>.mtx inside savepath)
# default output_name is 'cell-cpeaks' -> file: $OUTDIR/cell-cpeaks.mtx
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
