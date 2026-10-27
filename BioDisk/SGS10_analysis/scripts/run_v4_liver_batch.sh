#!/bin/bash
# ============================================================
# run_v4_liver_batch.sh
# ============================================================
# 批量跑 v4 二次注释 - 肝脏样本
#
# 已验证通过 (不重跑):
#   - con_liver  (2026-04-12 00:31)
#   - liver3     (2026-04-12 00:36)
#
# 本脚本跑剩下 3 个:
#   - bla_liver  (纯鼠空白对照)
#   - liver8     (PDX, 有 mouse + human)
#   - gan56      (PDX, 有 mouse + human)
#
# 用法:
#   cd ~/SGS10_analysis
#   bash scripts/run_v4_liver_batch.sh
#
# 输出 (每个样本):
#   <sample>/mouse/<sample>_mouse_annotated_v4.rds
#   <sample>/human/<sample>_human_annotated_v4.rds  (PDX only)
#   <sample>/<sample>_v4_annotation_report.md
#   <sample>/<sample>_v4_umap.pdf
#   scripts/logs/v4_<sample>.log
# ============================================================

set -e  # 任何一个样本失败立即停止

PROJECT_ROOT="/mnt/BioDisk/wangmingyuan/SGS10_analysis"
SCRIPT_DIR="${PROJECT_ROOT}/scripts"
LOG_DIR="${SCRIPT_DIR}/logs"

mkdir -p "${LOG_DIR}"

cd "${PROJECT_ROOT}"

SAMPLES=(bla_liver liver8 gan56)

echo "=========================================="
echo "SGS10 v4 Liver Batch Run"
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
    # 提取关键数字
    grep -E "(Cells:|Markers used|confidence:|Total flipped|rescore_agree)" "$log_file" \
      | head -20 || true
  else
    echo ""
    echo "[$(date '+%H:%M:%S')] ✗ $sample FAILED — see $log_file"
    exit 1
  fi
done

echo ""
echo "=========================================="
echo "All liver samples complete"
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
