#!/bin/bash
# ==================== 下载所有脑样本结果到本地 ====================
# 用法: 在本地电脑运行: scp -P 8022 wangmingyuan@g2.cnsknowall.com:/mnt/BioDisk/wangmingyuan/SGS10_analysis/scripts/download_brain_results.sh ./
#       bash download_brain_results.sh

echo "=========================================="
echo "下载脑样本分析结果"
echo "时间: $(date)"
echo "=========================================="
echo ""

# 创建本地目录
DOWNLOAD_DIR="$HOME/Downloads/brain_samples_$(date +%Y%m%d)"
mkdir -p "$DOWNLOAD_DIR"

samples=("nao56" "brain3" "brain8" "bla_brain" "con_brain")

for sample in "${samples[@]}"; do
    echo "下载 $sample ..."
    
    # 下载PDF文件
    scp -P 8022 wangmingyuan@g2.cnsknowall.com:"/mnt/BioDisk/wangmingyuan/SGS10_analysis/$sample/${sample}_*.pdf" "$DOWNLOAD_DIR/" 2>/dev/null
    
    # 下载CSV文件
    scp -P 8022 wangmingyuan@g2.cnsknowall.com:"/mnt/BioDisk/wangmingyuan/SGS10_analysis/$sample/${sample}_*.csv" "$DOWNLOAD_DIR/" 2>/dev/null
    
    # 下载RDS文件（可选，注释掉如果不需要）
    # scp -P 8022 wangmingyuan@g2.cnsknowall.com:"/mnt/BioDisk/wangmingyuan/SGS10_analysis/$sample/${sample}_*.rds" "$DOWNLOAD_DIR/" 2>/dev/null
done

echo ""
echo "下载完成！文件保存在: $DOWNLOAD_DIR"
ls -lh "$DOWNLOAD_DIR"
