#!/bin/bash

# 批量运行所有脑样本的脚本
# 用法: bash batch_run_all_brain_samples.sh

# 设置工作目录
ANALYSIS_DIR="/mnt/BioDisk/wangmingyuan/SGS10_analysis"
SCRIPT_DIR="${ANALYSIS_DIR}/scripts"
SCRIPT_NAME="brain_samples_analysis_optimized_final.R"
LOG_DIR="${ANALYSIS_DIR}/logs"

# 创建日志目录
mkdir -p ${LOG_DIR}

# 定义要运行的样本列表
SAMPLES=("nao56" "brain3" "brain8" "bla_brain" "con_brain")

# 核心数（每个样本使用的核心数）
CORES=25

# 是否运行精细注释
RUN_FINE_ANNOTATION="TRUE"

# 颜色输出
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo "=========================================="
echo "批量运行脑样本分析"
echo "样本数量: ${#SAMPLES[@]}"
echo "每个样本核心数: ${CORES}"
echo "精细注释: ${RUN_FINE_ANNOTATION}"
echo "日志目录: ${LOG_DIR}"
echo "=========================================="

# 清理函数
cleanup_sample() {
    local sample=$1
    local sample_dir="${ANALYSIS_DIR}/${sample}"
    
    echo -e "${YELLOW}清理 ${sample} 的中间文件...${NC}"
    
    if [ -d "${sample_dir}" ]; then
        cd ${sample_dir}
        
        # 删除中间RDS文件（保留raw, qc, norm, pca, cluster, umap, singler）
        # 但删除manual_annotated和complete以便重新生成
        rm -f ${sample}_manual_annotated.rds
        rm -f ${sample}_complete.rds
        rm -f ${sample}_checkpoint.txt
        rm -f ${sample}_manual_proportions.csv
        rm -f ${sample}_consensus_stats.csv
        rm -f ${sample}_newly_annotated_cells.csv
        rm -f ${sample}_threshold_quality_check.pdf
        rm -f ${sample}_annotation_comparison.pdf
        rm -f ${sample}_subtype_heatmap.pdf
        
        # 删除精细注释结果（如果需要完全重新运行）
        rm -f ${sample}_microglia_subset.rds
        rm -f ${sample}_neuron_subset.rds
        rm -f ${sample}_astro_subset.rds
        rm -f ${sample}_microglia_elbow.pdf
        rm -f ${sample}_microglia_states.pdf
        rm -f ${sample}_neuron_subtypes.pdf
        rm -f ${sample}_astro_regions.pdf
        rm -f ${sample}_microglia_state_markers.csv
        rm -f ${sample}_microglia_state_proportions.csv
        rm -f ${sample}_neuron_proportions.csv
        rm -f ${sample}_astro_proportions.csv
        
        echo -e "${GREEN}✓ ${sample} 清理完成${NC}"
    else
        echo -e "${RED}✗ ${sample} 目录不存在${NC}"
    fi
}

# 运行单个样本的函数
run_sample() {
    local sample=$1
    local log_file="${LOG_DIR}/${sample}_$(date +%Y%m%d_%H%M%S).log"
    
    echo -e "${GREEN}[$(date '+%Y-%m-%d %H:%M:%S')] 开始运行 ${sample}${NC}"
    
    cd ${SCRIPT_DIR}
    
    # 使用screen运行
    screen -dmS "${sample}_analysis" bash -c "
        echo '开始分析 ${sample}' > ${log_file}
        echo '时间: $(date)' >> ${log_file}
        echo '========================================' >> ${log_file}
        Rscript ${SCRIPT_NAME} ${sample} ${CORES} ${RUN_FINE_ANNOTATION} >> ${log_file} 2>&1
        echo '========================================' >> ${log_file}
        echo '分析完成时间: $(date)' >> ${log_file}
        echo '✓ ${sample} 分析完成' >> ${log_file}
    "
    
    echo -e "${GREEN}[$(date '+%Y-%m-%d %H:%M:%S')] ${sample} 已在screen会话中启动，会话名: ${sample}_analysis${NC}"
    echo "   查看日志: tail -f ${log_file}"
    echo "   连接会话: screen -r ${sample}_analysis"
}

# 主程序
echo ""
echo "选择操作:"
echo "1) 清理所有样本的中间文件并重新运行"
echo "2) 仅重新运行（不清理）"
echo "3) 仅清理中间文件（不运行）"
echo "4) 查看运行状态"
echo "5) 停止所有运行"
read -p "请输入选择 (1-5): " choice

case $choice in
    1)
        echo -e "${YELLOW}清理所有样本的中间文件...${NC}"
        for sample in "${SAMPLES[@]}"; do
            cleanup_sample ${sample}
        done
        
        echo -e "${YELLOW}开始运行所有样本...${NC}"
        for sample in "${SAMPLES[@]}"; do
            run_sample ${sample}
            sleep 2  # 避免同时启动过多进程
        done
        
        echo -e "\n${GREEN}所有样本已在screen会话中启动！${NC}"
        echo "查看运行状态: bash ${SCRIPT_DIR}/monitor_samples.sh"
        ;;
    2)
        echo -e "${YELLOW}开始运行所有样本（不清理）...${NC}"
        for sample in "${SAMPLES[@]}"; do
            run_sample ${sample}
            sleep 2
        done
        ;;
    3)
        echo -e "${YELLOW}清理所有样本的中间文件...${NC}"
        for sample in "${SAMPLES[@]}"; do
            cleanup_sample ${sample}
        done
        echo -e "${GREEN}清理完成！${NC}"
        ;;
    4)
        echo -e "${YELLOW}查看运行状态...${NC}"
        screen -ls | grep -E "nao56|brain3|brain8|bla_brain|con_brain"
        ;;
    5)
        echo -e "${RED}停止所有运行中的样本...${NC}"
        for sample in "${SAMPLES[@]}"; do
            screen -S "${sample}_analysis" -X quit 2>/dev/null && echo "停止 ${sample}" || echo "${sample} 未运行"
        done
        ;;
    *)
        echo -e "${RED}无效选择${NC}"
        exit 1
        ;;
esac

echo ""
echo "完成时间: $(date)"
