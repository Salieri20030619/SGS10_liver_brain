#!/bin/bash
# 文件名: ~/SGS10_analysis/scripts/package_liver_results.sh
# 用法: bash package_liver_results.sh

cd /mnt/BioDisk/wangmingyuan/SGS10_analysis

# 创建打包目录
PACKAGE_DIR="liver_samples_results_$(date +%Y%m%d_%H%M%S)"
mkdir -p $PACKAGE_DIR

echo "开始打包肝脏样本结果到: $PACKAGE_DIR"

# 复制每个样本的结果
for sample in gan56 liver3 liver8 bla_liver con_liver; do
    echo "处理样本: $sample"
    
    mkdir -p $PACKAGE_DIR/$sample
    
    # 复制主要结果文件
    cp $sample/${sample}_complete.rds $PACKAGE_DIR/$sample/ 2>/dev/null
    cp $sample/${sample}_celltype_proportions.csv $PACKAGE_DIR/$sample/ 2>/dev/null
    cp $sample/${sample}_immune_subtype_proportions.csv $PACKAGE_DIR/$sample/ 2>/dev/null
    cp $sample/${sample}_analysis_report.txt $PACKAGE_DIR/$sample/ 2>/dev/null
    cp $sample/*.pdf $PACKAGE_DIR/$sample/ 2>/dev/null
    cp $sample/${sample}_all_markers.csv $PACKAGE_DIR/$sample/ 2>/dev/null
    cp $sample/${sample}_top10_markers.csv $PACKAGE_DIR/$sample/ 2>/dev/null
done

# 创建汇总报告
cat > $PACKAGE_DIR/README.md << 'EOF'
# 肝脏样本单细胞转录组分析结果

## 分析说明
- 使用Seurat v5进行单细胞分析
- 人类基因转换为小鼠同源基因后注释
- 包含主要细胞类型和免疫细胞亚群分析

## 样本信息
- gan56: 肿瘤样本（人鼠混合）
- liver3: 肝脏样本3（人鼠混合）
- liver8: 肝脏样本8（人鼠混合）
- bla_liver: 空白对照肝脏
- con_liver: 对照肝脏

## 文件说明
每个样本文件夹包含：
- *_complete.rds: Seurat对象，可在R中继续分析
- *_celltype_proportions.csv: 细胞类型比例表
- *_immune_subtype_proportions.csv: 免疫细胞亚群比例表
- *_analysis_report.txt: 详细分析报告
- *.pdf: 所有可视化图表

## 软件要求
- R >= 4.4
- Seurat >= 5.4.0
- ggplot2, dplyr, SingleR等包
EOF

# 打包
tar -czf ${PACKAGE_DIR}.tar.gz $PACKAGE_DIR/

# 显示结果
echo ""
echo "========================================="
echo "打包完成！"
echo "打包文件: $(pwd)/${PACKAGE_DIR}.tar.gz"
echo "文件大小: $(du -h ${PACKAGE_DIR}.tar.gz | cut -f1)"
echo "========================================="
echo ""
echo "下载命令（在本地Windows PowerShell中执行）:"
echo "scp username@server_ip:$(pwd)/${PACKAGE_DIR}.tar.gz D:/"
