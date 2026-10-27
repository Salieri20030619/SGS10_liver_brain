#!/bin/bash

# 清理所有肝脏样本的中间文件

ANALYSIS_DIR="/mnt/BioDisk/wangmingyuan/SGS10_analysis"
SAMPLES=("gan56" "liver3" "liver8" "bla_liver" "con_liver")
LOG_DIR="${ANALYSIS_DIR}/logs/liver"

echo "=========================================="
echo "清理所有肝脏样本的中间文件"
echo "=========================================="

read -p "确认要清理所有样本并准备重新运行吗? (y/n): " confirm

if [ "$confirm" != "y" ]; then
    echo "取消操作"
    exit 0
fi

echo "停止所有运行中的样本..."
for sample in "${SAMPLES[@]}"; do
    screen -S "liver_${sample}" -X quit 2>/dev/null && echo "  停止 liver_${sample}"
done

echo ""
echo "清理中间文件..."
for sample in "${SAMPLES[@]}"; do
    sample_dir="${ANALYSIS_DIR}/${sample}"
    if [ -d "${sample_dir}" ]; then
        echo "清理 ${sample}..."
        cd ${sample_dir}
        
        # 删除中间结果，保留原始数据
        rm -f ${sample}_manual_annotated.rds
        rm -f ${sample}_complete.rds
        rm -f ${sample}_checkpoint.txt
        rm -f ${sample}_manual_proportions.csv
        rm -f ${sample}_consensus_stats.csv
        rm -f ${sample}_newly_annotated_cells.csv
        rm -f ${sample}_threshold_quality_check.pdf
        rm -f ${sample}_annotation_comparison.pdf
        rm -f ${sample}_subtype_heatmap.pdf
        rm -f ${sample}_complete_report.txt
        
        echo "  ✓ ${sample} 清理完成"
    else
        echo "  ✗ ${sample} 目录不存在"
    fi
done

# 清理日志
rm -rf ${LOG_DIR}
mkdir -p ${LOG_DIR}

echo ""
echo "✅ 清理完成！"
echo "运行 bash batch_run_all_liver_samples.sh 重新开始分析"
