#!/bin/bash
# ==================== 运行所有5个肝脏样本 ====================
# 文件名: ~/SGS10_analysis/scripts/run_liver_samples.sh

# 颜色设置
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_DIR="/mnt/BioDisk/wangmingyuan/SGS10_analysis/logs"
mkdir -p "$LOG_DIR"

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}    开始运行所有5个肝脏样本           ${NC}"
echo -e "${BLUE}    时间: $(date)                      ${NC}"
echo -e "${BLUE}========================================${NC}\n"

TOTAL_CORES=$(nproc)
SAMPLE_CORES=$((TOTAL_CORES / 10))
echo -e "${YELLOW}系统总核心数:${NC} $TOTAL_CORES"
echo -e "${YELLOW}每个样本使用:${NC} $SAMPLE_CORES 核心\n"

samples=("gan56" "liver3" "liver8" "bla_liver" "con_liver")

# 停止已有的会话
for sample in "${samples[@]}"; do
    screen -S "$sample" -X quit 2>/dev/null
done

# 删除检查点
read -p "是否删除所有旧的检查点文件？(y/n): " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo -e "${YELLOW}删除旧的检查点文件...${NC}"
    for sample in "${samples[@]}"; do
        rm -f "/mnt/BioDisk/wangmingyuan/SGS10_analysis/$sample/${sample}_checkpoint.txt" 2>/dev/null
    done
fi

# 启动样本
echo -e "\n${GREEN}启动所有5个肝脏样本...${NC}\n"
for sample in "${samples[@]}"; do
    echo -e "${BLUE}启动 ${sample}...${NC}"
    LOG_FILE="${LOG_DIR}/${sample}_${TIMESTAMP}.log"
    
    screen -dmS "$sample" bash -c "
        echo '========================================' > '$LOG_FILE'
        echo '开始运行 $sample 于 $(date)' >> '$LOG_FILE'
        cd /mnt/BioDisk/wangmingyuan/SGS10_analysis
        source ~/miniconda3/bin/activate r_4.4
        Rscript scripts/liver_samples_analysis.R $sample $SAMPLE_CORES TRUE >> '$LOG_FILE' 2>&1
        echo '完成 $sample 于 $(date)' >> '$LOG_FILE'
    "
    
    if [ $? -eq 0 ]; then
        echo -e "${GREEN}  ✓ $sample 已启动${NC}"
    else
        echo -e "${RED}  ✗ $sample 启动失败${NC}"
    fi
    sleep 2
done

echo -e "\n${GREEN}所有样本已启动！${NC}"
echo -e "${YELLOW}查看状态:${NC} screen -ls"
echo -e "${YELLOW}查看日志:${NC} tail -f $LOG_DIR/*.log"
