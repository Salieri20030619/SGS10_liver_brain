#!/bin/bash

# 一键运行所有肝脏样本（后台运行）

ANALYSIS_DIR="/mnt/BioDisk/wangmingyuan/SGS10_analysis"
SCRIPT_DIR="${ANALYSIS_DIR}/scripts"
LOG_DIR="${ANALYSIS_DIR}/logs/liver"
SAMPLES=("gan56" "liver3" "liver8" "bla_liver" "con_liver")
CORES=25

# 创建日志目录
mkdir -p ${LOG_DIR}

echo "=========================================="
echo "一键运行所有肝脏样本"
echo "=========================================="

# 先停止可能正在运行的
for sample in "${SAMPLES[@]}"; do
    screen -S "liver_${sample}" -X quit 2>/dev/null
done

echo "开始运行 ${#SAMPLES[@]} 个样本..."
echo ""

for sample in "${SAMPLES[@]}"; do
    log_file="${LOG_DIR}/${sample}_$(date +%Y%m%d_%H%M%S).log"
    
    echo "启动 ${sample}..."
    
    cd ${SCRIPT_DIR}
    
    screen -dmS "liver_${sample}" bash -c "
        echo '========================================' > ${log_file}
        echo '肝脏样本分析: ${sample}' >> ${log_file}
        echo '开始时间: $(date)' >> ${log_file}
        echo '========================================' >> ${log_file}
        Rscript liver_samples_analysis_complete.R ${sample} ${CORES} >> ${log_file} 2>&1
        echo '完成时间: $(date)' >> ${log_file}
    "
    
    echo "  ✓ ${sample} 已启动 (screen: liver_${sample}, 日志: ${log_file})"
    sleep 2
done

echo ""
echo "=========================================="
echo "所有样本已启动！"
echo ""
echo "查看状态: bash liver_status.sh"
echo "查看所有screen: screen -ls | grep liver"
echo "查看日志目录: ls -lh ${LOG_DIR}/"
echo "=========================================="
