#!/bin/bash
# ==================== 监控所有肝脏样本进度 ====================
# 文件名: ~/SGS10_analysis/scripts/monitor_liver_samples.sh

# 颜色设置
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
PURPLE='\033[0;35m'
NC='\033[0m'

samples=("gan56" "liver3" "liver8" "bla_liver" "con_liver")

steps=(
    "data_loaded:📥 加载数据"
    "qc_done:🔍 质控"
    "norm_done:📊 标准化"
    "pca_done:📉 PCA降维"
    "cluster_done:🔬 聚类"
    "umap_done:🗺️ UMAP"
    "singler_done:🏷️ SingleR注释"
    "marker_done:🧬 Marker基因"
    "kupffer_done:🦠 枯否细胞亚群"
    "hepatocyte_done:🧪 肝细胞亚群"
    "stellate_done:⭐ 星状细胞亚群"
)

while true; do
    clear
    echo -e "${BLUE}╔════════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${BLUE}║          肝脏样本完整进度监控 $(date '+%Y-%m-%d %H:%M:%S')          ║${NC}"
    echo -e "${BLUE}╚════════════════════════════════════════════════════════════════╝${NC}\n"
    
    echo -e "${YELLOW}系统资源:${NC}"
    free -h | head -2 | sed 's/^/  /'
    echo ""
    
    printf "${CYAN}%-12s${NC}" "样本"
    for step_info in "${steps[@]}"; do
        step_icon=$(echo $step_info | cut -d':' -f2)
        printf " ${CYAN}%s${NC}" "$step_icon"
    done
    printf " ${CYAN}⏱️ 时间${NC} ${CYAN}📊 细胞${NC}\n"
    printf "%0.s-" {1..130}
    echo ""
    
    for sample in "${samples[@]}"; do
        if screen -ls 2>/dev/null | grep -q "$sample"; then
            status_color=$YELLOW
            status="▶"
        else
            complete_file="/mnt/BioDisk/wangmingyuan/SGS10_analysis/$sample/${sample}_complete.rds"
            if [ -f "$complete_file" ]; then
                status_color=$GREEN
                status="✓"
            else
                status_color=$RED
                status="✗"
            fi
        fi
        printf "${status_color}%-12s${NC}" "$sample($status)"
        
        checkpoint_file="/mnt/BioDisk/wangmingyuan/SGS10_analysis/$sample/${sample}_checkpoint.txt"
        for step_info in "${steps[@]}"; do
            step_code=$(echo $step_info | cut -d':' -f1)
            if [ -f "$checkpoint_file" ] && grep -q "$step_code" "$checkpoint_file"; then
                printf " ${GREEN}✓${NC}   "
            else
                printf " ${RED}·${NC}   "
            fi
        done
        
        log_file=$(ls -t /mnt/BioDisk/wangmingyuan/SGS10_analysis/logs/${sample}_*.log 2>/dev/null | head -1)
        if [ -f "$log_file" ]; then
            start_time=$(grep -m 1 "开始处理" "$log_file" 2>/dev/null | grep -oE '[0-9]{4}-[0-9]{2}-[0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2}')
            if [ -n "$start_time" ]; then
                start_sec=$(date -d "$start_time" +%s 2>/dev/null)
                now_sec=$(date +%s)
                elapsed=$((now_sec - start_sec))
                printf " ${CYAN}%02d:%02d${NC}" $((elapsed/60)) $((elapsed%60))
            else
                printf " ${CYAN}--:--${NC}"
            fi
        else
            printf " ${CYAN}--:--${NC}"
        fi
        
        qc_file="/mnt/BioDisk/wangmingyuan/SGS10_analysis/$sample/${sample}_qc.rds"
        if [ -f "$qc_file" ]; then
            cells=$(Rscript -e "obj <- readRDS('$qc_file'); cat(ncol(obj))" 2>/dev/null)
            printf " ${PURPLE}%6s${NC}" "$cells"
        else
            printf " ${PURPLE}     -${NC}"
        fi
        echo ""
    done
    
    printf "%0.s-" {1..130}
    echo -e "\n"
    
    echo -e "${YELLOW}最新日志:${NC}"
    tail -n 5 /mnt/BioDisk/wangmingyuan/SGS10_analysis/logs/*.log 2>/dev/null | tail -5 || echo "  暂无日志"
    
    sleep 20
done
