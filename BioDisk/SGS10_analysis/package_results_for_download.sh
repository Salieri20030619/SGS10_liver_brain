#!/bin/bash

# ==================== 打包肝脑分析结果 ====================
# 文件名: package_results_for_download.sh
# 用途: 打包所有分析结果（排除RDS大文件）到D盘下载目录

echo "=========================================="
echo "打包肝脑分析结果"
echo "=========================================="

# 结果目录
RESULT_DIR="liver_brain_complete_results"
DOWNLOAD_DIR="liver_brain_download"

# 创建下载目录
mkdir -p ${DOWNLOAD_DIR}
mkdir -p ${DOWNLOAD_DIR}/plots
mkdir -p ${DOWNLOAD_DIR}/deg_tables
mkdir -p ${DOWNLOAD_DIR}/enrichment
mkdir -p ${DOWNLOAD_DIR}/crosstalk
mkdir -p ${DOWNLOAD_DIR}/threshold_summary

echo ""
echo "1. 复制PDF文件（火山图、热图）..."
cp ${RESULT_DIR}/plots/volcano/*.pdf ${DOWNLOAD_DIR}/plots/ 2>/dev/null
cp ${RESULT_DIR}/plots/heatmap/*.pdf ${DOWNLOAD_DIR}/plots/ 2>/dev/null
echo "   ✓ 复制火山图和热图"

echo ""
echo "2. 复制CSV文件（差异基因表）..."
cp ${RESULT_DIR}/deg_tables/*.csv ${DOWNLOAD_DIR}/deg_tables/ 2>/dev/null
echo "   ✓ 复制差异基因表"

echo ""
echo "3. 复制富集分析结果..."
cp ${RESULT_DIR}/liver/*.csv ${DOWNLOAD_DIR}/enrichment/ 2>/dev/null
cp ${RESULT_DIR}/brain/*.csv ${DOWNLOAD_DIR}/enrichment/ 2>/dev/null
cp ${RESULT_DIR}/crosstalk/*.csv ${DOWNLOAD_DIR}/enrichment/ 2>/dev/null
echo "   ✓ 复制富集结果"

echo ""
echo "4. 复制肝脑交互分析结果..."
cp ${RESULT_DIR}/crosstalk/*.pdf ${DOWNLOAD_DIR}/crosstalk/ 2>/dev/null
cp ${RESULT_DIR}/crosstalk/*.csv ${DOWNLOAD_DIR}/crosstalk/ 2>/dev/null
echo "   ✓ 复制肝脑交互结果"

echo ""
echo "5. 复制细胞数量统计..."
cp ${RESULT_DIR}/threshold_summary/*.csv ${DOWNLOAD_DIR}/threshold_summary/ 2>/dev/null
echo "   ✓ 复制细胞数量统计"

echo ""
echo "6. 复制分析报告..."
cp ${RESULT_DIR}/analysis_report.txt ${DOWNLOAD_DIR}/ 2>/dev/null
echo "   ✓ 复制分析报告"

echo ""
echo "7. 创建汇总README..."
cat > ${DOWNLOAD_DIR}/README.txt << 'EOF'
========================================
肝脑肿瘤完整分析结果
========================================

分析时间: 2026-03-31
方法: limma差异分析 + 逐级宽松阈值 + KEGG/GO富集

【分析参数】
- 阈值层级: Very_Strict (|log2FC|>1, adj.P<0.05) 到 Very_Relaxed (|log2FC|>0.4, P<0.05)
- 自动选择第一个产生≥10个DEG的阈值
- 最小细胞数: 10

【分析的细胞类型】
1. All_Cells - 全部细胞
2. Hepatocytes - 肝细胞（最主要细胞）
3. Endothelial - 内皮细胞（最主要亚型）
4. Kupffer - 枯否细胞（最主要亚型）

【比较类型】
- 肝脏内部: 肿瘤 vs 空白、肿瘤 vs 对照、空白 vs 对照
- 脑内部: 肿瘤对应脑 vs 脑空白、肿瘤对应脑 vs 脑对照、脑空白 vs 脑对照
- 跨器官: gan56 vs nao56, liver3 vs brain3, liver8 vs brain8

【输出文件说明】

1. plots/ - 可视化图
   - *_volcano.pdf : 火山图（标注top10基因）
   - *_heatmap.pdf : 热图（top30差异基因）

2. deg_tables/ - 差异基因表
   - *_all_genes.csv : 所有基因的差异分析结果
   - 包含: logFC, P.Value, adj.P.Val, sig (Up/Down/NS)

3. enrichment/ - 富集分析结果
   - *_KEGG.csv : KEGG通路富集结果
   - *_GO_BP.csv : GO生物学过程富集结果

4. crosstalk/ - 肝脑交互分析
   - crosstalk_shared_pathways.csv : 三个跨器官比较的共同通路
   - crosstalk_pathways_bubble.pdf : 共同通路气泡图

5. threshold_summary/ - 细胞数量统计
   - cell_counts.csv : 各样本各细胞类型细胞数

6. analysis_report.txt - 完整分析报告

【使用建议】
1. 先查看 analysis_report.txt 了解整体结果
2. 查看 plots/ 目录下的火山图和热图
3. 查看 enrichment/ 目录下的富集结果
4. 查看 crosstalk/ 目录下的肝脑交互结果

========================================
