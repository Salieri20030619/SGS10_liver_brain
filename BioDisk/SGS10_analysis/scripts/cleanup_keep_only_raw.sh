#!/bin/bash
# 文件名: ~/SGS10_analysis/scripts/cleanup_keep_only_raw.sh
# 用法: bash cleanup_keep_only_raw.sh [--dry-run]

BASE_DIR="/mnt/BioDisk/wangmingyuan/SGS10_analysis"
DRY_RUN=false

if [ "$1" == "--dry-run" ]; then
    DRY_RUN=true
    echo "🔍 试运行模式: 只显示将要删除的文件"
    echo ""
fi

echo "========================================="
echo "彻底清理 - 只保留原始数据文件"
echo "时间: $(date)"
echo "========================================="
echo ""

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# 脑样本和肝脏样本
BRAIN_SAMPLES=("nao56" "brain3" "brain8" "con_brain" "bla_brain")
LIVER_SAMPLES=("gan56" "liver3" "liver8" "con_liver" "bla_liver")
ALL_SAMPLES=("${BRAIN_SAMPLES[@]}" "${LIVER_SAMPLES[@]}")

# 统计变量
total_deleted=0
total_size=0
total_kept=0

# 删除文件函数
delete_file() {
    local file=$1
    local file_type=$2
    if [ -f "$file" ]; then
        size=$(stat -c%s "$file" 2>/dev/null || stat -f%z "$file" 2>/dev/null)
        size_mb=$((size / 1024 / 1024))
        
        if [ "$DRY_RUN" = true ]; then
            echo "  [DRY-RUN] 将删除 $file_type: $file ($size_mb MB)"
        else
            rm -f "$file"
            echo "  ✓ 删除 $file_type: $file ($size_mb MB)"
        fi
        return $size_mb
    fi
    return 0
}

# 清理单个样本目录
cleanup_sample() {
    local sample=$1
    local sample_dir="$BASE_DIR/$sample"
    
    if [ ! -d "$sample_dir" ]; then
        return
    fi
    
    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "清理样本: $sample"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    sample_deleted=0
    sample_size=0
    
    # 1. 删除所有分析结果文件
    echo ""
    echo "📁 删除分析结果文件..."
    
    # 删除所有RDS文件（除了原始数据）
    for rds in "$sample_dir"/*.rds; do
        if [ -f "$rds" ]; then
            # 检查是否是原始数据文件
            if [[ "$rds" == *"_raw.rds" ]]; then
                echo "  🔒 保留原始数据: $(basename "$rds")"
                ((total_kept++))
            else
                size_mb=$(delete_file "$rds" "RDS")
                if [ $? -gt 0 ]; then
                    ((sample_deleted++))
                    ((sample_size += size_mb))
                fi
            fi
        fi
    done
    
    # 2. 删除所有PDF文件
    echo ""
    echo "📊 删除PDF文件..."
    for pdf in "$sample_dir"/*.pdf; do
        if [ -f "$pdf" ]; then
            size_mb=$(delete_file "$pdf" "PDF")
            if [ $? -gt 0 ]; then
                ((sample_deleted++))
                ((sample_size += size_mb))
            fi
        fi
    done
    
    # 3. 删除所有CSV文件
    echo ""
    echo "📈 删除CSV文件..."
    for csv in "$sample_dir"/*.csv; do
        if [ -f "$csv" ]; then
            size_mb=$(delete_file "$csv" "CSV")
            if [ $? -gt 0 ]; then
                ((sample_deleted++))
                ((sample_size += size_mb))
            fi
        fi
    done
    
    # 4. 删除所有TXT报告文件
    echo ""
    echo "📄 删除报告文件..."
    for txt in "$sample_dir"/*.txt; do
        if [ -f "$txt" ]; then
            size_mb=$(delete_file "$txt" "TXT")
            if [ $? -gt 0 ]; then
                ((sample_deleted++))
                ((sample_size += size_mb))
            fi
        fi
    done
    
    # 5. 删除检查点文件
    echo ""
    echo "📌 删除检查点文件..."
    checkpoint="$sample_dir/${sample}_checkpoint.txt"
    if [ -f "$checkpoint" ]; then
        size_mb=$(delete_file "$checkpoint" "检查点")
        if [ $? -gt 0 ]; then
            ((sample_deleted++))
            ((sample_size += size_mb))
        fi
    fi
    
    # 6. 删除其他临时文件
    echo ""
    echo "🗑️ 删除临时文件..."
    temp_patterns=("*.tmp" "*.temp" "*.log" "*.out" "*.err")
    for pattern in "${temp_patterns[@]}"; do
        for temp in "$sample_dir"/$pattern; do
            if [ -f "$temp" ]; then
                size_mb=$(delete_file "$temp" "临时文件")
                if [ $? -gt 0 ]; then
                    ((sample_deleted++))
                    ((sample_size += size_mb))
                fi
            fi
        done
    done
    
    # 7. 删除子目录（亚群分析结果）
    echo ""
    echo "📁 删除子目录..."
    subdirs=("microglia_subset" "neuron_subset" "astro_subset" 
             "kupffer_subset" "hepatocyte_subset" "stellate_subset")
    for subdir in "${subdirs[@]}"; do
        if [ -d "$sample_dir/$subdir" ]; then
            if [ "$DRY_RUN" = true ]; then
                echo "  [DRY-RUN] 将删除目录: $sample_dir/$subdir"
            else
                rm -rf "$sample_dir/$subdir"
                echo "  ✓ 删除目录: $sample_dir/$subdir"
            fi
            ((sample_deleted++))
        fi
    done
    
    # 统计
    echo ""
    if [ "$DRY_RUN" = false ]; then
        echo "  📊 $sample: 已删除 $sample_deleted 个文件 (约 $sample_size MB)"
    else
        echo "  📊 $sample: 将删除 $sample_deleted 个文件 (约 $sample_size MB)"
    fi
    
    total_deleted=$((total_deleted + sample_deleted))
    total_size=$((total_size + sample_size))
}

# 清理日志目录
cleanup_logs() {
    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "清理日志目录"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    log_dir="$BASE_DIR/logs"
    if [ -d "$log_dir" ]; then
        log_count=0
        log_size=0
        
        for log in "$log_dir"/*; do
            if [ -f "$log" ]; then
                size_mb=$(delete_file "$log" "日志")
                if [ $? -gt 0 ]; then
                    ((log_count++))
                    ((log_size += size_mb))
                fi
            fi
        done
        
        # 删除空目录
        if [ "$DRY_RUN" = false ]; then
            rmdir "$log_dir" 2>/dev/null
        fi
        
        echo ""
        if [ "$DRY_RUN" = true ]; then
            echo "  📊 日志目录: 将删除 $log_count 个日志文件 (约 $log_size MB)"
        else
            echo "  📊 日志目录: 已删除 $log_count 个日志文件 (约 $log_size MB)"
        fi
        
        total_deleted=$((total_deleted + log_count))
        total_size=$((total_size + log_size))
    fi
}

# 清理结果目录
cleanup_results() {
    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "清理结果目录"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    results_dirs=("differential_results" "pseudotime_results" "results")
    
    for dir_name in "${results_dirs[@]}"; do
        dir_path="$BASE_DIR/$dir_name"
        if [ -d "$dir_path" ]; then
            if [ "$DRY_RUN" = true ]; then
                echo "  [DRY-RUN] 将删除目录: $dir_path"
            else
                rm -rf "$dir_path"
                echo "  ✓ 删除目录: $dir_path"
            fi
            ((total_deleted++))
        fi
    done
}

# 清理脚本目录中的临时文件
cleanup_scripts() {
    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "清理脚本目录临时文件"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    script_dir="$BASE_DIR/scripts"
    if [ -d "$script_dir" ]; then
        # 删除日志文件
        for log in "$script_dir"/*.log; do
            if [ -f "$log" ]; then
                size_mb=$(delete_file "$log" "脚本日志")
                if [ $? -gt 0 ]; then
                    ((total_deleted++))
                    ((total_size += size_mb))
                fi
            fi
        done
        
        # 删除临时文件
        for tmp in "$script_dir"/*.tmp "$script_dir"/*.temp; do
            if [ -f "$tmp" ]; then
                size_mb=$(delete_file "$tmp" "脚本临时文件")
                if [ $? -gt 0 ]; then
                    ((total_deleted++))
                    ((total_size += size_mb))
                fi
            fi
        done
    fi
}

# 显示磁盘空间
show_disk_usage() {
    echo ""
    echo "💾 磁盘使用情况:"
    df -h "$BASE_DIR" | awk 'NR==2 {printf "  当前: %s / %s (%s)\n", $3, $2, $5}'
}

# 显示原始数据统计
show_raw_data_stats() {
    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "保留的原始数据文件"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    for sample in "${ALL_SAMPLES[@]}"; do
        sample_dir="$BASE_DIR/$sample"
        raw_file="$sample_dir/${sample}_raw.rds"
        if [ -f "$raw_file" ]; then
            size=$(ls -lh "$raw_file" | awk '{print $5}')
            echo "  ✓ $sample: ${sample}_raw.rds ($size)"
        fi
    done
}

# 主清理流程
main() {
    # 显示清理前的磁盘使用
    echo "清理前:"
    show_disk_usage
    echo ""
    
    # 清理每个样本
    for sample in "${ALL_SAMPLES[@]}"; do
        cleanup_sample "$sample"
    done
    
    # 清理日志和结果目录
    cleanup_logs
    cleanup_results
    cleanup_scripts
    
    # 显示清理后的统计
    echo ""
    echo "========================================="
    if [ "$DRY_RUN" = true ]; then
        echo "🔍 试运行完成！"
        echo "将删除 $total_deleted 个文件/目录 (约 $total_size MB)"
        echo ""
        echo "使用以下命令执行实际清理:"
        echo "  bash cleanup_keep_only_raw.sh"
    else
        echo "✅ 清理完成！"
        echo "已删除 $total_deleted 个文件/目录 (约 $total_size MB)"
    fi
    echo "========================================="
    
    # 显示清理后的磁盘使用
    if [ "$DRY_RUN" = false ]; then
        echo ""
        echo "清理后:"
        show_disk_usage
        show_raw_data_stats
    fi
}

main
