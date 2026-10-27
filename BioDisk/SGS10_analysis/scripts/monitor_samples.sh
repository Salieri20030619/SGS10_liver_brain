#!/bin/bash

# 监控所有样本运行状态的脚本
# 用法: bash monitor_samples.sh [interval]

INTERVAL=${1:-30}  # 默认每30秒刷新一次

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# 样本列表
SAMPLES=("nao56" "brain3" "brain8" "bla_brain" "con_brain")

# 分析目录
ANALYSIS_DIR="/mnt/BioDisk/wangmingyuan/SGS10_analysis"
LOG_DIR="${ANALYSIS_DIR}/logs"

# 创建日志目录
mkdir -p ${LOG_DIR}

# 检查样本状态
check_sample() {
    local sample=$1
    local sample_dir="${ANALYSIS_DIR}/${sample}"
    local checkpoint_file="${sample_dir}/${sample}_checkpoint.txt"
    local manual_file="${sample_dir}/${sample}_manual_annotated.rds"
    local log_file="${LOG_DIR}/${sample}_*.log"
    
    # 检查screen会话
    if screen -ls | grep -q "${sample}_analysis"; then
        status="${GREEN}运行中${NC}"
        status_code="running"
    else
        status="${YELLOW}未运行${NC}"
        status_code="stopped"
    fi
    
    # 检查完成状态
    if [ -f "${manual_file}" ]; then
        completed="${GREEN}已完成${NC}"
    else
        completed="${RED}未完成${NC}"
    fi
    
    # 检查进度（从日志获取）
    progress="N/A"
    latest_log=$(ls -t ${LOG_DIR}/${sample}_*.log 2>/dev/null | head -1)
    if [ -f "${latest_log}" ]; then
        if grep -q "步骤8: 手动精细注释" ${latest_log} 2>/dev/null; then
            progress="步骤8"
        elif grep -q "步骤7: SingleR自动注释" ${latest_log} 2>/dev/null; then
            progress="步骤7"
        elif grep -q "步骤6: UMAP降维" ${latest_log} 2>/dev/null; then
            progress="步骤6"
        elif grep -q "步骤5: 聚类" ${latest_log} 2>/dev/null; then
            progress="步骤5"
        elif grep -q "分析完成" ${latest_log} 2>/dev/null; then
            progress="${GREEN}完成${NC}"
        fi
    fi
    
    # 检查文件大小
    manual_size="N/A"
    if [ -f "${manual_file}" ]; then
        manual_size=$(du -h ${manual_file} 2>/dev/null | cut -f1)
    fi
    
    printf "  %-15s | %-10s | %-10s | %-15s | %-10s\n" \
        "${sample}" "${status}" "${completed}" "${progress}" "${manual_size}"
}

# 获取CPU和内存使用情况
get_system_stats() {
    echo -e "${BLUE}系统资源使用情况:${NC}"
    echo "  CPU使用率: $(top -bn1 | grep "Cpu(s)" | awk '{print $2}' | cut -d'%' -f1)%"
    echo "  内存使用: $(free -h | awk 'NR==2{print $3 "/" $2}')"
    echo "  磁盘使用: $(df -h ${ANALYSIS_DIR} | awk 'NR==2{print $5 " (" $3 "/" $2 ")"}')"
    echo ""
}

# 显示活跃的screen会话
show_screen_sessions() {
    echo -e "${BLUE}活跃的Screen会话:${NC}"
    screen -ls | grep -E "nao56|brain3|brain8|bla_brain|con_brain" || echo "  无活跃会话"
    echo ""
}

# 显示最近的日志
show_recent_logs() {
    echo -e "${BLUE}最近的日志文件:${NC}"
    ls -lt ${LOG_DIR}/*.log 2>/dev/null | head -5 | awk '{print "  " $9 " (" $5 ")"}'
    echo ""
}

# 监控主循环
monitor() {
    while true; do
        clear
        echo "=========================================="
        echo "        脑样本分析监控面板"
        echo "=========================================="
        echo "时间: $(date '+%Y-%m-%d %H:%M:%S')"
        echo "刷新间隔: ${INTERVAL}秒"
        echo "=========================================="
        echo ""
        
        show_screen_sessions
        get_system_stats
        
        echo -e "${BLUE}样本运行状态:${NC}"
        printf "  %-15s | %-10s | %-10s | %-15s | %-10s\n" "样本" "Screen状态" "完成状态" "当前步骤" "文件大小"
        echo "  ---------------------------------------------------------------------------------"
        
        for sample in "${SAMPLES[@]}"; do
            check_sample ${sample}
        done
        
        echo ""
        show_recent_logs
        
        echo "=========================================="
        echo "命令:"
        echo "  q - 退出监控"
        echo "  l - 查看最新日志"
        echo "  c - 清除已完成样本的日志"
        echo "  s - 查看screen会话详情"
        echo "=========================================="
        
        read -t ${INTERVAL} -n 1 key
        if [ "$key" = "q" ]; then
            echo -e "\n退出监控"
            break
        elif [ "$key" = "l" ]; then
            echo -e "\n选择样本查看日志:"
            select sample in "${SAMPLES[@]}" "取消"; do
                if [ "$sample" = "取消" ]; then
                    break
                elif [ -n "$sample" ]; then
                    latest_log=$(ls -t ${LOG_DIR}/${sample}_*.log 2>/dev/null | head -1)
                    if [ -f "${latest_log}" ]; then
                        less +G "${latest_log}"
                    else
                        echo "未找到 ${sample} 的日志"
                    fi
                    break
                fi
            done
        elif [ "$key" = "c" ]; then
            echo -e "\n清理已完成样本的日志..."
            for sample in "${SAMPLES[@]}"; do
                manual_file="${ANALYSIS_DIR}/${sample}/${sample}_manual_annotated.rds"
                if [ -f "${manual_file}" ]; then
                    rm -f ${LOG_DIR}/${sample}_*.log
                    echo "  已清理 ${sample} 的日志"
                fi
            done
        elif [ "$key" = "s" ]; then
            echo -e "\nScreen会话详情:"
            screen -ls | grep -E "nao56|brain3|brain8|bla_brain|con_brain"
            echo ""
            read -p "按回车继续..." dummy
        fi
    done
}

# 运行监控
monitor
