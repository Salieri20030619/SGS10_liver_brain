#!/bin/bash

# 打包所有肝脏样本的分析结果（排除RDS文件）
# 用法: bash package_liver_results.sh

ANALYSIS_DIR="/mnt/BioDisk/wangmingyuan/SGS10_analysis"
OUTPUT_DIR="${ANALYSIS_DIR}/liver_results_download"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

# 样本列表
SAMPLES=("gan56" "liver3" "liver8" "bla_liver" "con_liver")

# 创建输出目录
mkdir -p ${OUTPUT_DIR}

echo "=========================================="
echo "打包肝脏样本分析结果（排除RDS文件）"
echo "时间: $(date)"
echo "输出目录: ${OUTPUT_DIR}"
echo "=========================================="

# 创建汇总README
README_FILE="${OUTPUT_DIR}/README_${TIMESTAMP}.txt"
cat > ${README_FILE} << 'EOF'
肝脏样本单细胞分析结果汇总
==========================

分析时间: 2026-03-26
样本列表: gan56, liver3, liver8, bla_liver, con_liver

各样本文件夹内容说明:
--------------------
PDF文件（可视化结果）:
  - annotation_comparison.pdf : SingleR与手动注释对比图
  - subtype_heatmap.pdf : 亚型标记基因热图
  - threshold_quality_check.pdf : 阈值调整质量检查
  - clusters.pdf : UMAP聚类图
  - singler_celltypes.pdf : SingleR自动注释图
  - QC.pdf : 质控图
  - variable_features.pdf : 高变基因图
  - elbow.pdf : PCA肘部图

CSV文件（统计数据）:
  - manual_proportions.csv : 最终细胞类型比例（最重要）
  - singler_proportions.csv : SingleR注释比例
  - newly_annotated_cells.csv : 新增注释细胞列表
  - consensus_stats.csv : 共识注释统计

TXT文件（报告）:
  - complete_report.txt : 完整分析报告（最重要）
  - checkpoint.txt : 分析步骤检查点

注意: 完整的Seurat对象（*_manual_annotated.rds，约200-800MB）未包含在此包中。
如需完整对象进行下游分析，请单独下载。

