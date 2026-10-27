#!/bin/bash

# 清理所有样本的中间文件，准备重新运行

ANALYSIS_DIR="/mnt/BioDisk/wangmingyuan/SGS10_analysis"
SAMPLES=("nao56" "brain3" "brain8" "bla_brain" "con_brain")

echo "=========================================="
echo "清理所有样本的中间文件"
echo "=========================================="

read -p "确认要清理所有样本并准备重新运行吗? (y/n): " confirm

if [ "$confirm" != "y" ]; then
    echo "取消操作"
    exit 0
fi

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
        rm -f ${sample}_microglia_*
        rm -f ${sample}_neuron_*
        rm -f ${sample}_astro_*
        
        echo "  ✓ ${sample} 清理完成"
    else
        echo "  ✗ ${sample} 目录不存在"
    fi
done

# 清理日志
rm -rf ${ANALYSIS_DIR}/logs
mkdir -p ${ANALYSIS_DIR}/logs

echo ""
echo "清理完成！"
echo "运行 bash batch_run_all_brain_samples.sh 重新开始分析"
