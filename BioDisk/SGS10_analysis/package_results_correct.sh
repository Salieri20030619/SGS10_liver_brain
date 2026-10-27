#!/bin/bash

# 正确的打包脚本
ANALYSIS_DIR="/mnt/BioDisk/wangmingyuan/SGS10_analysis"
OUTPUT_DIR="${ANALYSIS_DIR}/results_for_download"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

# 创建输出目录
mkdir -p ${OUTPUT_DIR}

echo "=========================================="
echo "打包脑样本分析结果（排除RDS大文件）"
echo "工作目录: ${ANALYSIS_DIR}"
echo "输出目录: ${OUTPUT_DIR}"
echo "时间: $(date)"
echo "=========================================="

SAMPLES=("nao56" "brain3" "brain8" "bla_brain" "con_brain")

for sample in "${SAMPLES[@]}"; do
    sample_dir="${ANALYSIS_DIR}/${sample}"
    
    if [ ! -d "${sample_dir}" ]; then
        echo "⚠️ ${sample} 目录不存在，跳过"
        continue
    fi
    
    echo ""
    echo "处理 ${sample}..."
    
    # 创建样本子目录
    sample_output="${OUTPUT_DIR}/${sample}"
    mkdir -p ${sample_output}
    
    # 复制所有PDF文件
    cp ${sample_dir}/${sample}_*.pdf ${sample_output}/ 2>/dev/null
    pdf_count=$(ls ${sample_output}/${sample}_*.pdf 2>/dev/null | wc -l)
    echo "  ✓ 复制 ${pdf_count} 个PDF文件"
    
    # 复制所有CSV文件
    cp ${sample_dir}/${sample}_*.csv ${sample_output}/ 2>/dev/null
    csv_count=$(ls ${sample_output}/${sample}_*.csv 2>/dev/null | wc -l)
    echo "  ✓ 复制 ${csv_count} 个CSV文件"
    
    # 复制报告文件
    cp ${sample_dir}/${sample}_complete_report.txt ${sample_output}/ 2>/dev/null
    cp ${sample_dir}/${sample}_checkpoint.txt ${sample_output}/ 2>/dev/null
    
    # 复制精细亚群RDS（这些比较小）
    for subtype in microglia neuron astro; do
        if [ -f "${sample_dir}/${sample}_${subtype}_subset.rds" ]; then
            size=$(du -h ${sample_dir}/${sample}_${subtype}_subset.rds | cut -f1)
            cp ${sample_dir}/${sample}_${subtype}_subset.rds ${sample_output}/ 2>/dev/null
            echo "  ✓ 复制 ${subtype}_subset.rds (${size})"
        fi
    done
done

echo ""
echo "=========================================="
echo "打包完成！"
echo "输出目录: ${OUTPUT_DIR}"
echo "=========================================="

# 创建汇总报告
SUMMARY_FILE="${OUTPUT_DIR}/README_${TIMESTAMP}.txt"
cat > ${SUMMARY_FILE} << 'EOF'
脑样本分析结果汇总
==================

本目录包含所有5个脑样本的分析结果（排除大的Seurat对象RDS文件）

各样本文件夹内容说明：
----------------------
- *_annotation_comparison.pdf : SingleR与手动注释对比图
- *_subtype_heatmap.pdf : 亚型标记基因热图
- *_threshold_quality_check.pdf : 阈值调整质量检查图
- *_manual_proportions.csv : 细胞类型比例统计
- *_consensus_stats.csv : 共识注释统计
- *_newly_annotated_cells.csv : 新增注释细胞列表
- *_complete_report.txt : 完整分析报告
- *_microglia_states.pdf : 小胶质细胞状态图
- *_neuron_subtypes.pdf : 神经元亚型图
- *_astro_regions.pdf : 星形胶质细胞区域图
- *_microglia_subset.rds : 小胶质细胞亚群对象
- *_neuron_subset.rds : 神经元亚群对象
- *_astro_subset.rds : 星形胶质细胞亚群对象

完整的Seurat对象（*_manual_annotated.rds）大小约1.1-1.3GB，未包含在此包中。
如需完整对象，请单独下载。

分析时间: $(date)
