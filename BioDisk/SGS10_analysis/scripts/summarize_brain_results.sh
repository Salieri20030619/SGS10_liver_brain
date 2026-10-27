#!/bin/bash
# ==================== 汇总脑样本分析结果 ====================

samples=("nao56" "brain3" "brain8" "bla_brain" "con_brain")

echo "=========================================="
echo "脑样本分析结果汇总"
echo "时间: $(date)"
echo "=========================================="
echo ""

# 统计各样本完成情况
echo "1. 完成状态:"
echo "------------------------------------------"
printf "%-10s | %-12s | %-10s | %-15s\n" "样本" "细胞数" "完成时间" "文件大小"
echo "------------------------------------------"

for sample in "${samples[@]}"; do
    complete_file="/mnt/BioDisk/wangmingyuan/SGS10_analysis/$sample/${sample}_complete.rds"
    if [ -f "$complete_file" ]; then
        # 获取文件修改时间
        mod_time=$(stat -c "%y" "$complete_file" | cut -d'.' -f1)
        size=$(du -h "$complete_file" | cut -f1)
        
        # 获取细胞数（从qc文件）
        qc_file="/mnt/BioDisk/wangmingyuan/SGS10_analysis/$sample/${sample}_qc.rds"
        if [ -f "$qc_file" ]; then
            cells=$(Rscript -e "
                obj <- readRDS('$qc_file')
                cat(ncol(obj))
            " 2>/dev/null)
        else
            cells="-"
        fi
        
        printf "%-10s | %-12s | %-10s | %-15s\n" "$sample" "$cells" "$mod_time" "$size"
    else
        printf "%-10s | %-12s | %-10s | %-15s\n" "$sample" "-" "未完成" "-"
    fi
done

echo ""

# 统计各样本生成的文件
echo "2. 生成的文件统计:"
echo "------------------------------------------"
printf "%-10s | %-8s | %-8s | %-8s | %-10s\n" "样本" "PDF" "CSV" "RDS" "总大小"
echo "------------------------------------------"

for sample in "${samples[@]}"; do
    sample_dir="/mnt/BioDisk/wangmingyuan/SGS10_analysis/$sample"
    if [ -d "$sample_dir" ]; then
        pdf_count=$(ls "$sample_dir"/*.pdf 2>/dev/null | wc -l)
        csv_count=$(ls "$sample_dir"/*.csv 2>/dev/null | wc -l)
        rds_count=$(ls "$sample_dir"/*.rds 2>/dev/null | wc -l)
        total_size=$(du -sh "$sample_dir" 2>/dev/null | cut -f1)
        printf "%-10s | %-8s | %-8s | %-8s | %-10s\n" "$sample" "$pdf_count" "$csv_count" "$rds_count" "$total_size"
    fi
done

echo "------------------------------------------"
echo ""
echo "3. 细胞类型比例 (top 3):"
echo "------------------------------------------"

for sample in "${samples[@]}"; do
    prop_file="/mnt/BioDisk/wangmingyuan/SGS10_analysis/$sample/${sample}_celltype_proportions.csv"
    if [ -f "$prop_file" ]; then
        echo "$sample:"
        cat "$prop_file" | head -4 | column -t -s ','
        echo ""
    fi
done
