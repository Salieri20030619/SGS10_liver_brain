#!/bin/bash
# 文件名: ~/SGS10_analysis/scripts/run_brain_samples.sh

# 设置核心数（每个样本使用10%）
export TOTAL_CORES=$(nproc)
export SAMPLE_CORES=$((TOTAL_CORES / 10))
echo "总核心数: $TOTAL_CORES"
echo "每个脑样本使用: $SAMPLE_CORES 核心"
echo ""

# 5个脑样本列表
brain_samples=(
    "nao56"
    "brain3"
    "brain8"
    "bla_brain"
    "con_brain"
)

# 创建日志目录
mkdir -p ~/SGS10_analysis/logs

# 启动所有脑样本
echo "开始启动所有5个脑样本..."
for sample in "${brain_samples[@]}"; do
    echo "启动 $sample ..."
    screen -dmS "$sample" bash -c "
        cd /mnt/BioDisk/wangmingyuan/SGS10_analysis
        source ~/miniconda3/bin/activate r_4.4
        Rscript ~/SGS10_analysis/scripts/brain_samples_analysis.R $sample $SAMPLE_CORES > logs/${sample}.log 2>&1
    "
    sleep 2  # 避免同时启动造成I/O冲击
done

echo ""
echo "所有5个脑样本已启动！"
echo "查看状态: screen -ls"
echo "查看日志: tail -f ~/SGS10_analysis/logs/*.log"
echo ""
echo "各样本screen会话名:"
for sample in "${brain_samples[@]}"; do
    echo "  $sample"
done
