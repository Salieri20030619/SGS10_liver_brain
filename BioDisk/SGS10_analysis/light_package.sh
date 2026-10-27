#!/bin/bash
cd ~/SGS10_analysis

echo "========================================="
echo "创建轻量级打包文件（不含 RDS）"
echo "========================================="

# 检查目录是否存在
if [ ! -d "pseudotime_multicurve_optimized" ]; then
    echo "错误：目录不存在"
    exit 1
fi

# 显示原始大小
echo "原始目录大小:"
du -sh pseudotime_multicurve_optimized/

# 统计 RDS 文件
echo ""
echo "RDS 文件信息:"
find pseudotime_multicurve_optimized -name "*.rds" -exec ls -lh {} \; | awk '{print $5, $9}'
RDS_SIZE=$(find pseudotime_multicurve_optimized -name "*.rds" -exec du -b {} + | awk '{sum+=$1} END {print sum/1024/1024 " MB"}')
echo "RDS 文件总大小: $RDS_SIZE"

# 创建轻量级打包
echo ""
echo "创建轻量级打包文件..."
find pseudotime_multicurve_optimized -type f \
  \( -name "*.pdf" -o -name "*.png" -o -name "*.csv" -o -name "*.txt" \) \
  | tar -czf pseudotime_multicurve_optimized_light.tar.gz -T -

# 显示打包结果
echo ""
echo "打包完成！"
ls -lh pseudotime_multicurve_optimized_light.tar.gz

echo ""
echo "文件统计:"
tar -tzf pseudotime_multicurve_optimized_light.tar.gz | wc -l
echo "个文件"

echo ""
echo "========================================="
echo "下载命令:"
echo "scp -P 8022 wangmingyuan@g2.cnsknowall.com:~/SGS10_analysis/pseudotime_multicurve_optimized_light.tar.gz D:/"
echo "========================================="
