#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DATA="$ROOT/tests/data/toy_fragments.tsv"
OUTDIR="$ROOT/tests/output"
WEIGHTS="$ROOT/tests/weights/model.pth"

mkdir -p "$OUTDIR"
mkdir -p "$(dirname "$WEIGHTS")"

# 如果需要权重且尚未提交到仓库，可通过环境变量 WEIGHTS_URL 下载（CI 中我们会传）
if [[ -n "${WEIGHTS_URL:-}" && ! -f "$WEIGHTS" ]]; then
  echo "Downloading weights from WEIGHTS_URL..."
  curl -L -o "$WEIGHTS" "$WEIGHTS_URL"
fi

# --------- 这里替换为你项目中实际调用 map2cpeak 的命令 ---------
# 举例（请按实际 CLI/参数替换）：
map2cpeak --fragments "$DATA" --out "$OUTDIR/peaks.bed" --weights "$WEIGHTS"
# ----------------------------------------------------------------

# 简单断言：输出文件存在并且非空
if [[ ! -s "$OUTDIR/peaks.bed" ]]; then
  echo "Smoke test failed: output file missing or empty" >&2
  ls -l "$OUTDIR" || true
  exit 2
fi

echo "Smoke test passed: $OUTDIR/peaks.bed"
