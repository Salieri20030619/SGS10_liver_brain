#!/bin/bash
# ==================== 监控所有脑样本完整进度 ====================

# 颜色设置
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m'

samples=("nao56" "brain3" "brain8" "bla_brain" "con_brain")

# 定义所有步骤
steps=(
    "data_loaded:📥 加载数据"
    "qc_done:🔍 质控"
    "norm_done:📊 标准化"
    "pca_done:📉 PCA降维"
    "cluster_done:🔬 聚类"
    "umap_done:🗺️ UMAP"
    "singler_done:🏷️ SingleR注释"
    "marker_done:🧬 Marker基因"
    "microglia_done:🦠 小胶质细胞亚群"
    "neuron_done:🧠 神经元亚群"
    "astro_done:⭐ 星形胶质细胞亚群"
)

while true; do
    clear
    echo -e "${BLUE}========================================${NC}"
    echo -e "${BLUE}      脑样本完整进度监控 $(date)${NC}"
    echo -e "${BLUE}========================================${NC}\n"
    
    # 显示系统资源
    echo -e "${YELLOW}系统资源:${NC}"
    free -h | head -2 | sed 's/^/  /'
    echo ""
    
    # 显示各样本进度表格
    printf "${CYAN}%-12s${NC}" "样本"
    for step_info in "${steps[@]}"; do
        step_name=$(echo $step_info | cut -d':' -f1)
        step_icon=$(echo $step_info | cut -d':' -f2)
        printf " ${CYAN}%s${NC}" "$step_icon"
    done
    printf " ${CYAN}⏱️ 时间${NC}\n"
    
    # 分隔线
    printf "%0.s-" {1..120}
    echo ""
    
    for sample in "${samples[@]}"; do
        # 检查screen状态
        if screen -ls 2>/dev/null | grep -q "$sample"; then
            status_color=$YELLOW
            status="▶"
        else
            # 检查是否全部完成
            complete_file="/mnt/BioDisk/wangmingyuan/SGS10_analysis/$sample/${sample}_complete.rds"
            if [ -f "$complete_file" ]; then
                status_color=$GREEN
                status="✓"
            else
                status_color=$RED
                status="✗"
            fi
        fi
        
        # 样本名
        printf "${status_color}%-12s${NC}" "$sample($status)"
        
        # 检查每个步骤
        checkpoint_file="/mnt/BioDisk/wangmingyuan/SGS10_analysis/$sample/${sample}_checkpoint.txt"
        
        for step_info in "${steps[@]}"; do
            step_code=$(echo $step_info | cut -d':' -f1)
            
            if [ -f "$checkpoint_file" ] && grep -q "$step_code" "$checkpoint_file"; then
                # 步骤已完成
                printf " ${GREEN}✓${NC}   "
            else
                # 检查是否有对应的输出文件（备用判断）
                case $step_code in
                    "neuron_done")
                        if [ -f "/mnt/BioDisk/wangmingyuan/SGS10_analysis/$sample/${sample}_neuron_subset.rds" ]; then
                            printf " ${GREEN}✓${NC}   "
                        else
                            printf " ${RED}·${NC}   "
                        fi
                        ;;
                    "astro_done")
                        if [ -f "/mnt/BioDisk/wangmingyuan/SGS10_analysis/$sample/${sample}_astro_subset.rds" ]; then
                            printf " ${GREEN}✓${NC}   "
                        else
                            printf " ${RED}·${NC}   "
                        fi
                        ;;
                    "microglia_done")
                        if [ -f "/mnt/BioDisk/wangmingyuan/SGS10_analysis/$sample/${sample}_microglia_subset.rds" ]; then
                            printf " ${GREEN}✓${NC}   "
                        else
                            printf " ${RED}·${NC}   "
                        fi
                        ;;
                    "marker_done")
                        if [ -f "/mnt/BioDisk/wangmingyuan/SGS10_analysis/$sample/${sample}_all_markers.csv" ]; then
                            printf " ${GREEN}✓${NC}   "
                        else
                            printf " ${RED}·${NC}   "
                        fi
                        ;;
                    "complete_done")
                        if [ -f "/mnt/BioDisk/wangmingyuan/SGS10_analysis/$sample/${sample}_complete.rds" ]; then
                            printf " ${GREEN}✓${NC}   "
                        else
                            printf " ${RED}·${NC}   "
                        fi
                        ;;
                    *)
                        printf " ${RED}·${NC}   "
                        ;;
                esac
            fi
        done
        
        # 显示运行时间
        log_file="/mnt/BioDisk/wangmingyuan/SGS10_analysis/logs/${sample}.log"
        if [ -f "$log_file" ]; then
            # 获取日志中最早的时间戳
            start_time=$(grep -m 1 "开始处理" "$log_file" 2>/dev/null | grep -oE '[0-9]{4}-[0-9]{2}-[0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2}')
            if [ -n "$start_time" ]; then
                start_sec=$(date -d "$start_time" +%s 2>/dev/null)
                now_sec=$(date +%s)
                if [ -n "$start_sec" ]; then
                    elapsed=$((now_sec - start_sec))
                    elapsed_min=$((elapsed / 60))
                    elapsed_sec=$((elapsed % 60))
                    printf " ${CYAN}%02d:%02d${NC}" $elapsed_min $elapsed_sec
                else
                    printf " ${CYAN}--:--${NC}"
                fi
            else
                printf " ${CYAN}--:--${NC}"
            fi
        else
            printf " ${CYAN}--:--${NC}"
        fi
        
        echo ""
    done
    
    # 分隔线
    printf "%0.s-" {1..120}
    echo ""
    
    # 图例
    echo -e "\n${YELLOW}图例:${NC}"
    echo "  ${GREEN}✓${NC} 已完成  ${RED}·${NC} 未完成  ${YELLOW}▶${NC} 运行中  ${GREEN}✓${NC} 已完成"
    echo ""
    
    # 显示当前最新日志
    echo -e "${YELLOW}最新运行日志:${NC}"
    for sample in "${samples[@]}"; do
        log_file="/mnt/BioDisk/wangmingyuan/SGS10_analysis/logs/${sample}.log"
        if [ -f "$log_file" ]; then
            # 获取最后几行中包含关键信息的行
            last_line=$(tail -n 3 "$log_file" 2>/dev/null | grep -E "步骤|完成|神经元|星形胶质|小胶质" | tail -1)
            if [ -n "$last_line" ]; then
                echo -e "  ${BLUE}${sample}:${NC} ${last_line:0:80}"
            fi
        fi
    done
    
    sleep 20
done
