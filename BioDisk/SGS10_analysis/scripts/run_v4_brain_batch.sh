#!/bin/bash
# ============================================================
# run_v4_brain_batch.sh
# ============================================================
# 批量跑 v4 二次注释 - 脑样本
#
# 已验证通过 (不重跑):
#   - brain3  (2026-04-12 02:07, v4.1-c 最终版)
#
# 本脚本跑剩下 4 个:
#   - con_brain  (正常对照)
#   - bla_brain  (空白对照)
#   - nao56      (PDX)
#   - brain8     (PDX)
#
# 用法:
#   cd ~/SGS10_analysis
#   bash scripts/run_v4_brain_batch.sh
#
# 输出 (每个样本):
#   <sample>/mouse/<sample>_manual_annotated_v4.rds
#   <sample>/<sample>_v4_annotation_report.md
#   <sample>/<sample>_v4_umap.pdf
#   scripts/logs/v4_<sample>.log
# ============================================================

set -e

PROJECT_ROOT="/mnt/BioDisk/wangmingyuan/SGS10_analysis"
SCRIPT_DIR="${PROJECT_ROOT}/scripts"
LOG_DIR="${SCRIPT_DIR}/logs"

mkdir -p "${LOG_DIR}"

cd "${PROJECT_ROOT}"

SAMPLES=(con_brain bla_brain nao56 brain8)

echo "=========================================="
echo "SGS10 v4 Brain Batch Run"
echo "Start time: $(date '+%Y-%m-%d %H:%M:%S')"
echo "Samples:    ${SAMPLES[*]}"
echo "=========================================="
echo ""

for sample in "${SAMPLES[@]}"; do
  echo ""
  echo "=========================================="
  echo "[$(date '+%H:%M:%S')] Running: $sample"
  echo "=========================================="

  log_file="${LOG_DIR}/v4_${sample}.log"

  start_ts=$(date +%s)
  if Rscript scripts/secondary_annotation_v4.R "$sample" 2>&1 | tee "$log_file"; then
    end_ts=$(date +%s)
    elapsed=$((end_ts - start_ts))
    echo ""
    echo "[$(date '+%H:%M:%S')] ✓ $sample DONE (${elapsed}s)"
    grep -E "(Cells:|Markers used|confidence:|Total flipped|rescore_agree|subtype_rescore)" "$log_file" \
      | head -20 || true
  else
    echo ""
    echo "[$(date '+%H:%M:%S')] ✗ $sample FAILED — see $log_file"
    exit 1
  fi
done

echo ""
echo "=========================================="
echo "All brain samples complete"
echo "End time: $(date '+%Y-%m-%d %H:%M:%S')"
echo "=========================================="
echo ""
echo "Reports:"
for sample in "${SAMPLES[@]}"; do
  report="${PROJECT_ROOT}/${sample}/${sample}_v4_annotation_report.md"
  if [ -f "$report" ]; then
    echo "  ✓ $report"
  else
    echo "  ✗ missing: $report"
  fi
done
