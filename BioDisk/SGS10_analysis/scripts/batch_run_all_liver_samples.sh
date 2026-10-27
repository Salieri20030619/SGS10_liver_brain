#!/bin/bash

# 批量运行所有肝脏样本的脚本
# 用法: bash batch_run_all_liver_samples.sh

# 设置工作目录
ANALYSIS_DIR="/mnt/BioDisk/wangmingyuan/SGS10_analysis"
SCRIPT_DIR="${ANALYSIS_DIR}/scripts"
SCRIPT_NAME="liver_samples_analysis_complete.R"
LOG_DIR="${ANALYSIS_DIR}/logs/liver"

# 创建日志目录
mkdir -p ${LOG_DIR}

# 定义要运行的样本列表
SAMPLES=("gan56" "liver3" "liver8" "bla_liver" "con_liver")

# 核心数（每个样本使用的核心数）
CORES=25

# 颜色输出
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo "=========================================="
echo "批量运行肝脏样本分析"
echo "样本数量: ${#SAMPLES[@]}"
echo "每个样本核心数: ${CORES}"
echo "日志目录: ${LOG_DIR}"
echo "=========================================="

# 清理函数
cleanup_sample() {
    local sample=$1
    local sample_dir="${ANALYSIS_DIR}/${sample}"
    
    echo -e "${YELLOW}清理 ${sample} 的中间文件...${NC}"
    
    if [ -d "${sample_dir}" ]; then
        cd ${sample_dir}
        
        # 删除中间结果，保留原始数据
        rm -f ${sample}_manual_annotated.rds
        rm -f ${sample}_complete.rds
        rm -f ${sample}_checkpoint.txt
        rm -f ${sample}_manual_proportions.csv
        rm -f ${sample}_consensus_stats.csv
        rm -f ${sample}_newly_annotated_cells.csv
        rm -f ${sample}_threshold_quality_check.pdf
        rm -f ${sample}_annotation_comparison.pdf
        rm -f ${sample}_subtype_heatmap.pdf
        rm -f ${sample}_complete_report.txt
        
        echo -e "${GREEN}  ✓ ${sample} 清理完成${NC}"
    else
        echo -e "${RED}  ✗ ${sample} 目录不存在${NC}"
    fi
}

# 检查样本状态
check_sample_status() {
    local sample=$1
    local sample_dir="${ANALYSIS_DIR}/${sample}"
    local manual_file="${sample_dir}/${sample}_manual_annotated.rds"
    
    if [ -f "${manual_file}" ]; then
        echo -e "${GREEN}已完成${NC}"
    else
        echo -e "${YELLOW}未完成${NC}"
    fi
}

# 运行单个样本的函数
run_sample() {
    local sample=$1
    local log_file="${LOG_DIR}/${sample}_$(date +%Y%m%d_%H%M%S).log"
    
    echo -e "${GREEN}[$(date '+%Y-%m-%d %H:%M:%S')] 开始运行 ${sample}${NC}"
    
    cd ${SCRIPT_DIR}
    
    # 使用screen运行
    screen -dmS "liver_${sample}" bash -c "
        echo '========================================' > ${log_file}
        echo '肝脏样本分析: ${sample}' >> ${log_file}
        echo '开始时间: $(date)' >> ${log_file}
        echo '========================================' >> ${log_file}
        Rscript ${SCRIPT_NAME} ${sample} ${CORES} >> ${log_file} 2>&1
        echo '========================================' >> ${log_file}
        echo '完成时间: $(date)' >> ${log_file}
        echo '✓ ${sample} 分析完成' >> ${log_file}
        
        # 检查是否成功
        if [ -f \"${ANALYSIS_DIR}/${sample}/${sample}_manual_annotated.rds\" ]; then
            echo '✅ 成功生成 manual_annotated.rds' >> ${log_file}
        else
            echo '❌ 未生成 manual_annotated.rds' >> ${log_file}
        fi
    "
    
    echo -e "${GREEN}  ✓ ${sample} 已在screen会话中启动，会话名: liver_${sample}${NC}"
    echo "    查看日志: tail -f ${log_file}"
    echo "    连接会话: screen -r liver_${sample}"
}

# 显示所有样本状态
show_status() {
    echo ""
    echo -e "${BLUE}=========================================="
    echo "当前样本状态"
    echo "==========================================${NC}"
    printf "%-15s %-15s\n" "样本" "状态"
    echo "----------------------------------------"
    for sample in "${SAMPLES[@]}"; do
        status=$(check_sample_status ${sample})
        printf "%-15s %-15s\n" "${sample}" "${status}"
    done
    echo ""
}

# 停止所有运行中的样本
stop_all() {
    echo -e "${RED}停止所有运行中的样本...${NC}"
    for sample in "${SAMPLES[@]}"; do
        screen -S "liver_${sample}" -X quit 2>/dev/null && echo "  停止 liver_${sample}" || echo "  liver_${sample} 未运行"
    done
}

# 主菜单
while true; do
    clear
    echo "=========================================="
    echo "    肝脏样本批量分析管理脚本"
    echo "=========================================="
    echo ""
    show_status
    echo ""
    echo "请选择操作:"
    echo "  1) 清理所有样本并重新运行"
    echo "  2) 仅重新运行（不清理）"
    echo "  3) 仅清理中间文件（不运行）"
    echo "  4) 查看运行状态"
    echo "  5) 查看实时日志"
    echo "  6) 停止所有运行"
    echo "  7) 退出"
    echo ""
    read -p "请输入选择 (1-7): " choice
    
    case $choice in
        1)
            echo -e "${YELLOW}清理所有样本的中间文件...${NC}"
            for sample in "${SAMPLES[@]}"; do
                cleanup_sample ${sample}
            done
            
            echo -e "${YELLOW}开始运行所有样本...${NC}"
            for sample in "${SAMPLES[@]}"; do
                run_sample ${sample}
                sleep 2
            done
            
            echo -e "\n${GREEN}所有样本已在screen会话中启动！${NC}"
            echo "查看状态: screen -ls | grep liver"
            echo "查看日志: tail -f ${LOG_DIR}/liver*.log"
            echo ""
            read -p "按回车继续..." dummy
            ;;
        2)
            echo -e "${YELLOW}开始运行所有样本（不清理）...${NC}"
            for sample in "${SAMPLES[@]}"; do
                run_sample ${sample}
                sleep 2
            done
            echo ""
            read -p "按回车继续..." dummy
            ;;
        3)
            echo -e "${YELLOW}清理所有样本的中间文件...${NC}"
            for sample in "${SAMPLES[@]}"; do
                cleanup_sample ${sample}
            done
            echo -e "${GREEN}清理完成！${NC}"
            echo ""
            read -p "按回车继续..." dummy
            ;;
        4)
            echo ""
            echo -e "${BLUE}Screen会话列表:${NC}"
            screen -ls | grep liver || echo "  无活跃会话"
            echo ""
            echo -e "${BLUE}最近日志文件:${NC}"
            ls -lt ${LOG_DIR}/*.log 2>/dev/null | head -5 | awk '{print "  " $9 " (" $5 ")"}'
            echo ""
            read -p "按回车继续..." dummy
            ;;
        5)
            echo ""
            echo "选择要查看日志的样本:"
            select sample in "${SAMPLES[@]}" "退出"; do
                if [ "$sample" = "退出" ]; then
                    break
                elif [ -n "$sample" ]; then
                    latest_log=$(ls -t ${LOG_DIR}/${sample}_*.log 2>/dev/null | head -1)
                    if [ -f "${latest_log}" ]; then
                        tail -f "${latest_log}"
                    else
                        echo "未找到 ${sample} 的日志文件"
                    fi
                    break
                fi
            done
            ;;
        6)
            stop_all
            echo ""
            read -p "按回车继续..." dummy
            ;;
        7)
            echo -e "${GREEN}退出${NC}"
            exit 0
            ;;
        *)
            echo -e "${RED}无效选择${NC}"
            sleep 1
            ;;
    esac
done
