#!/bin/bash

# 快速打包脚本 - 只打包PDF和CSV，完全排除RDS
OUTPUT_DIR="/mnt/BioDisk/wangmingyuan/SGS10_analysis/download_results"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

mkdir -p ${OUTPUT_DIR}

echo "=========================================="
echo "快速打包PDF和CSV文件（排除所有RDS）"
echo "=========================================="

SAMPLES=("nao56" "brain3" "brain8" "bla_brain" "con_brain")

# 收集所有需要打包的文件
ALL_FILES=""

for sample in "${SAMPLES[@]}"; do
    sample_dir="/mnt/BioDisk/wangmingyuan/SGS10_analysis/${sample}"
    
    if [ -d "${sample_dir}" ]; then
        echo "收集 ${sample} 的文件..."
        
        # 创建样本子目录
        mkdir -p ${OUTPUT_DIR}/${sample}
        
        # 复制PDF文件
        cp ${sample_dir}/${sample}_*.pdf ${OUTPUT_DIR}/${sample}/ 2>/dev/null
        pdf_count=$(ls ${OUTPUT_DIR}/${sample}/*.pdf 2>/dev/null | wc -l)
        
        # 复制CSV文件
        cp ${sample_dir}/${sample}_*.csv ${OUTPUT_DIR}/${sample}/ 2>/dev/null
        csv_count=$(ls ${OUTPUT_DIR}/${sample}/*.csv 2>/dev/null | wc -l)
        
        # 复制报告文件
        cp ${sample_dir}/${sample}_*report.txt ${OUTPUT_DIR}/${sample}/ 2>/dev/null
        cp ${sample_dir}/${sample}_checkpoint.txt ${OUTPUT_DIR}/${sample}/ 2>/dev/null
        
        echo "  ✓ ${sample}: ${pdf_count} PDF, ${csv_count} CSV"
    fi
done

echo ""
echo "=========================================="
echo "创建压缩包..."
echo "=========================================="

cd ${OUTPUT_DIR}

# 方法1：直接打包整个目录（会自动包含所有文件）
ZIP_FILE="brain_results_${TIMESTAMP}.zip"
zip -r ${ZIP_FILE} ./* > /dev/null 2>&1

if [ -f "${ZIP_FILE}" ]; then
    ZIP_SIZE=$(du -h ${ZIP_FILE} | cut -f1)
    echo "✓ 压缩包创建成功: ${ZIP_FILE}"
    echo "  大小: ${ZIP_SIZE}"
    echo ""
    echo "文件列表:"
    unzip -l ${ZIP_FILE} | head -20
else
    echo "⚠️ 压缩失败，但文件已保存在: ${OUTPUT_DIR}"
fi

echo ""
echo "=========================================="
echo "下载方式:"
echo "=========================================="
echo "1. 下载整个目录:"
echo "   scp -r wangmingyuan@server:${OUTPUT_DIR} ./"
echo ""
echo "2. 下载压缩包:"
echo "   scp wangmingyuan@server:${OUTPUT_DIR}/${ZIP_FILE} ./"
echo ""
echo "3. 查看文件:"
echo "   cd ${OUTPUT_DIR} && ls -lh"
echo "=========================================="

# 显示目录内容
ls -lh ${OUTPUT_DIR}/
