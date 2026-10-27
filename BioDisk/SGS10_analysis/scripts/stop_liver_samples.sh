#!/bin/bash
# ==================== 停止所有肝脏样本 ====================
# 文件名: ~/SGS10_analysis/scripts/stop_liver_samples.sh

samples=("gan56" "liver3" "liver8" "bla_liver" "con_liver")

echo "正在停止所有肝脏样本..."

for sample in "${samples[@]}"; do
    if screen -ls 2>/dev/null | grep -q "$sample"; then
        screen -S "$sample" -X quit
        echo "  ✓ $sample 已停止"
    else
        echo "  ✗ $sample 未运行"
    fi
done

echo "所有样本已停止！"
