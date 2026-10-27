#!/bin/bash

# 设置固定路径
BASE_DIR="/mnt/BioDisk/wangmingyuan/SGS10_analysis"
OUTPUT_DIR="${BASE_DIR}/download_results"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

# 清理旧目录
rm -rf ${OUTPUT_DIR}
mkdir -p ${OUTPUT_DIR}

echo "=========================================="
echo "打包脑样本分析结果"
echo "基础目录: ${BASE_DIR}"
echo "输出目录: ${OUTPUT_DIR}"
echo "时间: $(date)"
echo "=========================================="

# 样本列表
SAMPLES=("nao56" "brain3" "brain8" "bla_brain" "con_brain")

# 创建汇总README
README_FILE="${OUTPUT_DIR}/README_${TIMESTAMP}.txt"
cat > ${README_FILE} << 'EOF'
脑样本单细胞分析结果汇总
======================

分析时间: 2026-03-25
样本列表: nao56, brain3, brain8, bla_brain, con_brain

各样本文件夹内容说明:
--------------------
PDF文件（可视化结果）:
  - annotation_comparison.pdf : SingleR与手动注释对比图
  - subtype_heatmap.pdf : 亚型标记基因热图
  - threshold_quality_check.pdf : 阈值调整质量检查
  - clusters.pdf : UMAP聚类图
  - singler_celltypes.pdf : SingleR自动注释图
  - microglia_states.pdf : 小胶质细胞状态图
  - neuron_subtypes.pdf : 神经元亚型图
  - astro_regions.pdf : 星形胶质细胞区域图
  - QC.pdf : 质控图
  - variable_features.pdf : 高变基因图
  - elbow.pdf : PCA肘部图

CSV文件（统计数据）:
  - manual_proportions.csv : 最终细胞类型比例（最重要）
  - consensus_stats.csv : 共识注释统计
  - newly_annotated_cells.csv : 新增注释细胞列表
  - singler_proportions.csv : SingleR注释比例
  - microglia_state_proportions.csv : 小胶质细胞状态比例
  - neuron_proportions.csv : 神经元亚型比例
  - astro_proportions.csv : 星形胶质细胞区域比例

TXT文件（报告）:
  - complete_report.txt : 完整分析报告（最重要）
  - checkpoint.txt : 分析步骤检查点

注意: 完整的Seurat对象（*_manual_annotated.rds，约1.1-1.3GB）未包含在此包中。

