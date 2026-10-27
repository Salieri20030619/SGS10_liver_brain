#!/bin/bash
# ==================== 运行所有5个脑样本 ====================
# 文件名: ~/SGS10_analysis/scripts/run_all_brain_samples.sh
# 用法: ./run_all_brain_samples.sh

# 设置颜色输出
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# 设置时间戳
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_DIR="/mnt/BioDisk/wangmingyuan/SGS10_analysis/logs"
mkdir -p "$LOG_DIR"

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}    开始运行所有5个脑样本              ${NC}"
echo -e "${BLUE}    时间: $(date)                      ${NC}"
echo -e "${BLUE}========================================${NC}\n"

# 检查核心数
TOTAL_CORES=$(nproc)
SAMPLE_CORES=$((TOTAL_CORES / 10))
echo -e "${YELLOW}系统总核心数:${NC} $TOTAL_CORES"
echo -e "${YELLOW}每个样本使用:${NC} $SAMPLE_CORES 核心"
echo -e "${YELLOW}总内存使用预估:${NC} 约150GB\n"

# 定义样本列表
samples=(
    "nao56"
    "brain3"
    "brain8"
    "bla_brain"
    "con_brain"
)

# 清理函数
cleanup() {
    echo -e "\n${YELLOW}清理旧的screen会话...${NC}"
    for sample in "${samples[@]}"; do
        screen -S "$sample" -X quit 2>/dev/null
        echo "  已停止: $sample"
    done
}

# 停止已有的screen会话
cleanup

# 删除旧的检查点文件（可选，如果需要完全重新开始）
read -p "是否删除所有旧的检查点文件？(y/n): " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo -e "${YELLOW}删除旧的检查点文件...${NC}"
    for sample in "${samples[@]}"; do
        rm -f "/mnt/BioDisk/wangmingyuan/SGS10_analysis/$sample/${sample}_checkpoint.txt" 2>/dev/null
        echo "  已删除: $sample 的检查点"
    done
fi

# 启动所有样本
echo -e "\n${GREEN}启动所有5个脑样本...${NC}\n"

for sample in "${samples[@]}"; do
    echo -e "${BLUE}启动 ${sample}...${NC}"
    
    # 创建日志文件
    LOG_FILE="${LOG_DIR}/${sample}_${TIMESTAMP}.log"
    
    # 在screen中运行
    screen -dmS "$sample" bash -c "
        echo '========================================' > '$LOG_FILE'
        echo '开始运行 $sample 于 $(date)' >> '$LOG_FILE'
        echo '========================================' >> '$LOG_FILE'
        echo '' >> '$LOG_FILE'
        
        cd /mnt/BioDisk/wangmingyuan/SGS10_analysis
        
        # 激活环境并运行
        source ~/miniconda3/bin/activate r_4.4
        Rscript scripts/brain_samples_analysis_complete.R $sample $SAMPLE_CORES >> '$LOG_FILE' 2>&1
        
        echo '' >> '$LOG_FILE'
        echo '========================================' >> '$LOG_FILE'
        echo '完成 $sample 于 $(date)' >> '$LOG_FILE'
        echo '========================================' >> '$LOG_FILE'
    "
    
    # 检查是否启动成功
    if [ $? -eq 0 ]; then
        echo -e "${GREEN}  ✓ $sample 已启动${NC} (日志: ${LOG_FILE})"
    else
        echo -e "${RED}  ✗ $sample 启动失败${NC}"
    fi
    
    # 等待2秒避免I/O冲突
    sleep 2
done

echo -e "\n${GREEN}所有样本已启动！${NC}"
echo -e "${YELLOW}查看运行状态:${NC} screen -ls"
echo -e "${YELLOW}查看日志目录:${NC} ls -lh $LOG_DIR"
echo -e "${YELLOW}实时查看所有日志:${NC} tail -f $LOG_DIR/*.log"

# ==================== 生成监控脚本 ====================

cat > ~/SGS10_analysis/scripts/monitor_brain_samples.sh << 'EOF'
#!/bin/bash
# ==================== 监控所有脑样本进度 ====================

# 颜色设置
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

samples=("nao56" "brain3" "brain8" "bla_brain" "con_brain")

while true; do
    clear
    echo -e "${BLUE}========================================${NC}"
    echo -e "${BLUE}      脑样本运行监控 $(date)${NC}"
    echo -e "${BLUE}========================================${NC}\n"
    
    # 显示screen会话
    echo -e "${YELLOW}Screen会话状态:${NC}"
    screen -ls | grep -E "(nao56|brain3|brain8|bla_brain|con_brain)" || echo "  无运行中"
    echo ""
    
    # 显示系统资源
    echo -e "${YELLOW}系统资源:${NC}"
    free -h | head -2
    echo ""
    
    # 显示各样本进度
    echo -e "${YELLOW}样本进度:${NC}"
    echo "------------------------------------------------------------------------"
    printf "%-12s | %-40s | %-10s\n" "样本" "当前步骤" "文件大小"
    echo "------------------------------------------------------------------------"
    
    for sample in "${samples[@]}"; do
        checkpoint_file="/mnt/BioDisk/wangmingyuan/SGS10_analysis/$sample/${sample}_checkpoint.txt"
        complete_file="/mnt/BioDisk/wangmingyuan/SGS10_analysis/$sample/${sample}_complete.rds"
        
        if [ -f "$complete_file" ]; then
            # 已完成
            size=$(du -h "$complete_file" | cut -f1)
            printf "${GREEN}%-12s | %-40s | %-10s${NC}\n" "$sample" "✅ 已完成" "$size"
        elif [ -f "$checkpoint_file" ]; then
            # 进行中
            steps=$(cat "$checkpoint_file" | tr '\n' ' ')
            steps=${steps:0:40}
            if [ -f "/mnt/BioDisk/wangmingyuan/SGS10_analysis/$sample/${sample}_complete.rds" ]; then
                size=$(du -h "/mnt/BioDisk/wangmingyuan/SGS10_analysis/$sample/${sample}_complete.rds" 2>/dev/null | cut -f1)
                size=${size:-"-"}
            else
                size="-"
            fi
            printf "${YELLOW}%-12s | %-40s | %-10s${NC}\n" "$sample" "$steps" "$size"
        else
            # 未开始
            printf "${RED}%-12s | %-40s | %-10s${NC}\n" "$sample" "⏳ 未开始" "-"
        fi
    done
    echo "------------------------------------------------------------------------"
    
    # 显示最新日志
    echo -e "\n${YELLOW}最新日志片段:${NC}"
    LOG_DIR="/mnt/BioDisk/wangmingyuan/SGS10_analysis/logs"
    if [ -d "$LOG_DIR" ]; then
        latest_log=$(ls -t "$LOG_DIR"/*.log 2>/dev/null | head -1)
        if [ -f "$latest_log" ]; then
            echo "  最近日志: $(basename "$latest_log")"
            echo "  ................................"
            tail -n 5 "$latest_log" 2>/dev/null | sed 's/^/  /'
        else
            echo "  暂无日志"
        fi
    fi
    
    sleep 30
done
EOF

chmod +x ~/SGS10_analysis/scripts/monitor_brain_samples.sh

# ==================== 生成结果汇总脚本 ====================

cat > ~/SGS10_analysis/scripts/summarize_brain_results.sh << 'EOF'
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
EOF

chmod +x ~/SGS10_analysis/scripts/summarize_brain_results.sh

# ==================== 生成下载脚本 ====================

cat > ~/SGS10_analysis/scripts/download_brain_results.sh << 'EOF'
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
EOF

chmod +x ~/SGS10_analysis/scripts/download_brain_results.sh

echo -e "\n${GREEN}辅助脚本已生成:${NC}"
echo "  1. 监控脚本: monitor_brain_samples.sh"
echo "  2. 汇总脚本: summarize_brain_results.sh"
echo "  3. 下载脚本: download_brain_results.sh"
echo -e "\n${YELLOW}开始监控进度:${NC} ./monitor_brain_samples.sh"
echo -e "${YELLOW}查看汇总结果:${NC} ./summarize_brain_results.sh"
