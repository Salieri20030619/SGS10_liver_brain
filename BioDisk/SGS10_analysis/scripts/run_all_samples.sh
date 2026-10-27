#!/bin/bash
# ==================== 主调度脚本 ====================
# 文件名: ~/SGS10_analysis/scripts/run_all_samples.sh
# 用法:
#   bash run_all_samples.sh                # 跑全部 10 个样本
#   bash run_all_samples.sh liver          # 只跑 5 个肝样本
#   bash run_all_samples.sh brain          # 只跑 5 个脑样本
#   bash run_all_samples.sh single gan56   # 只跑指定单样本

set -e
set -u

SCRIPT_DIR="/mnt/BioDisk/wangmingyuan/SGS10_analysis/scripts"
LOG_DIR="/mnt/BioDisk/wangmingyuan/SGS10_analysis/logs"
N_CORES=8       # 每个样本使用的核心数,根据机器调整

mkdir -p "$LOG_DIR"

LIVER_SAMPLES=(gan56 liver3 liver8 bla_liver con_liver)
BRAIN_SAMPLES=(nao56 brain3 brain8 bla_brain con_brain)

run_liver() {
    local sample=$1
    local log="$LOG_DIR/${sample}_liver.log"
    echo "[$(date '+%F %T')] >>> 开始肝样本: $sample"
    if Rscript "$SCRIPT_DIR/liver_pdx_analysis.R" "$sample" "$N_CORES" \
        > "$log" 2>&1; then
        echo "[$(date '+%F %T')] ✓ 肝样本完成: $sample (日志: $log)"
    else
        echo "[$(date '+%F %T')] ✗ 肝样本失败: $sample (日志: $log)"
    fi
}

run_brain() {
    local sample=$1
    local log="$LOG_DIR/${sample}_brain.log"
    echo "[$(date '+%F %T')] >>> 开始脑样本: $sample"
    if Rscript "$SCRIPT_DIR/brain_analysis.R" "$sample" "$N_CORES" \
        > "$log" 2>&1; then
        echo "[$(date '+%F %T')] ✓ 脑样本完成: $sample (日志: $log)"
    else
        echo "[$(date '+%F %T')] ✗ 脑样本失败: $sample (日志: $log)"
    fi
}

MODE="${1:-all}"

case "$MODE" in
    all)
        for s in "${LIVER_SAMPLES[@]}"; do run_liver "$s"; done
        for s in "${BRAIN_SAMPLES[@]}"; do run_brain "$s"; done
        ;;
    liver)
        for s in "${LIVER_SAMPLES[@]}"; do run_liver "$s"; done
        ;;
    brain)
        for s in "${BRAIN_SAMPLES[@]}"; do run_brain "$s"; done
        ;;
    single)
        if [ $# -lt 2 ]; then
            echo "用法: bash run_all_samples.sh single <sample_name>"; exit 1
        fi
        s="$2"
        if [[ " ${LIVER_SAMPLES[*]} " =~ " $s " ]]; then
            run_liver "$s"
        elif [[ " ${BRAIN_SAMPLES[*]} " =~ " $s " ]]; then
            run_brain "$s"
        else
            echo "未知样本: $s"; exit 1
        fi
        ;;
    *)
        echo "用法: bash run_all_samples.sh [all|liver|brain|single <sample>]"
        exit 1
        ;;
esac

echo
echo "========================================"
echo "全部任务结束: $(date '+%F %T')"
echo "日志目录: $LOG_DIR"
echo "========================================"
