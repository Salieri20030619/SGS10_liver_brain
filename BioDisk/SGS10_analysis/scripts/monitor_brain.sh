#!/bin/bash
# 文件名: ~/SGS10_analysis/scripts/monitor_brain.sh

brain_samples=("nao56" "brain3" "brain8" "bla_brain" "con_brain")

while true; do
    clear
    echo "========== 5个脑样本运行状态 =========="
    echo "时间: $(date)"
    echo ""
    
    # 显示screen会话
    echo "Screen会话:"
    screen -ls | grep -E "(nao56|brain3|brain8|bla_brain|con_brain)" || echo "  无运行中"
    echo ""
    
    # 显示资源使用
    echo "系统资源:"
    free -h | head -2
    echo ""
    echo "CPU使用:"
    top -bn1 | grep "Cpu(s)" | awk '{print "  " $0}'
    echo ""
    
    # 显示各样本进度
    echo "样本进度:"
    for sample in "${brain_samples[@]}"; do
        if [ -f "/mnt/BioDisk/wangmingyuan/SGS10_analysis/$sample/${sample}_checkpoint.txt" ]; then
            steps=$(cat "/mnt/BioDisk/wangmingyuan/SGS10_analysis/$sample/${sample}_checkpoint.txt" | tr '\n' ' ')
            echo "  $sample: $steps"
            
            # 显示当前步骤的文件大小
            if [[ $steps == *"data_loaded"* ]] && [ -f "/mnt/BioDisk/wangmingyuan/SGS10_analysis/$sample/${sample}_raw.rds" ]; then
                size=$(du -h "/mnt/BioDisk/wangmingyuan/SGS10_analysis/$sample/${sample}_raw.rds" | cut -f1)
                echo "    raw: $size"
            fi
            if [[ $steps == *"qc_done"* ]] && [ -f "/mnt/BioDisk/wangmingyuan/SGS10_analysis/$sample/${sample}_qc.rds" ]; then
                size=$(du -h "/mnt/BioDisk/wangmingyuan/SGS10_analysis/$sample/${sample}_qc.rds" | cut -f1)
                echo "    qc: $size"
            fi
        else
            echo "  $sample: 未开始"
        fi
    done
    
    # 显示最新日志
    echo ""
    echo "最新日志:"
    tail -n 5 ~/SGS10_analysis/logs/*.log 2>/dev/null | head -10 || echo "  无日志"
    
    sleep 30
done
