#!/bin/bash

# 批量运行所有肝脏样本（从头开始）
SAMPLES=("gan56" "liver3" "liver8" "bla_liver" "con_liver")
CORES=25

echo "=========================================="
echo "批量运行所有肝脏样本（从头开始）"
echo "样本列表: ${SAMPLES[@]}"
echo "=========================================="

for sample in "${SAMPLES[@]}"; do
    echo ""
    echo "=========================================="
    echo "处理 ${sample}..."
    echo "=========================================="
    
    cd /mnt/BioDisk/wangmingyuan/SGS10_analysis/${sample}
    
    # 删除所有现有文件（保留原始数据raw.rds）
    echo "清理旧文件..."
    rm -f ${sample}_qc.rds ${sample}_norm.rds ${sample}_pca.rds ${sample}_cluster.rds ${sample}_umap.rds
    rm -f ${sample}_singler.rds ${sample}_manual_annotated.rds ${sample}_complete.rds
    rm -f ${sample}_*.pdf ${sample}_*.csv ${sample}_*.txt
    rm -f ${sample}_checkpoint.txt
    
    echo "开始分析..."
    Rscript /mnt/BioDisk/wangmingyuan/SGS10_analysis/scripts/liver_complete_analysis_clean.R ${sample} ${CORES}
    
    cd ..
done

echo ""
echo "=========================================="
echo "所有样本分析完成！"
echo "=========================================="
