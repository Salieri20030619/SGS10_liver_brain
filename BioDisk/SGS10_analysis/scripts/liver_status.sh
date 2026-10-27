#!/bin/bash

# 快速查看所有肝脏样本状态

ANALYSIS_DIR="/mnt/BioDisk/wangmingyuan/SGS10_analysis"
SAMPLES=("gan56" "liver3" "liver8" "bla_liver" "con_liver")
LOG_DIR="${ANALYSIS_DIR}/logs/liver"

echo "=========================================="
echo "肝脏样本分析状态"
echo "时间: $(date '+%Y-%m-%d %H:%M:%S')"
echo "=========================================="
printf "%-15s %-12s %-15s %-20s\n" "样本" "Screen状态" "完成状态" "manual_annotated大小"
echo "--------------------------------------------------------"

for sample in "${SAMPLES[@]}"; do
    # Screen状态
    if screen -ls 2>/dev/null | grep -q "liver_${sample}"; then
        screen_status="运行中"
    else
        screen_status="未运行"
    fi
    
    # 完成状态
    manual_file="${ANALYSIS_DIR}/${sample}/${sample}_manual_annotated.rds"
    if [ -f "${manual_file}" ]; then
        file_size=$(du -h ${manual_file} 2>/dev/null | cut -f1)
        complete_status="已完成 (${file_size})"
    else
        complete_status="未完成"
    fi
    
    printf "%-15s %-12s %-15s %-20s\n" "${sample}" "${screen_status}" "${complete_status}" ""
done

echo "=========================================="
echo "Screen会话列表:"
screen -ls 2>/dev/null | grep liver || echo "  无活跃会话"
echo "=========================================="

# 显示最近日志
if [ -d "${LOG_DIR}" ]; then
    echo "最近更新的日志:"
    ls -lt ${LOG_DIR}/*.log 2>/dev/null | head -5 | awk '{print "  " $9 " (" $5 ")"}'
fi
