#!/bin/bash
# ==================== 一键复制所有肝脏样本数据并运行分析 ====================
# 文件名: ~/SGS10_analysis/scripts/setup_and_run_liver.sh
# 用法: bash setup_and_run_liver.sh

# 颜色设置
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}  肝脏样本数据复制与运行脚本${NC}"
echo -e "${BLUE}  时间: $(date)${NC}"
echo -e "${BLUE}========================================${NC}\n"

# ==================== 1. 创建目录 ====================
echo -e "${YELLOW}1. 创建数据目录...${NC}"
mkdir -p /mnt/BioDisk/wangmingyuan/SGS10_analysis/data/liver_samples/
mkdir -p /mnt/BioDisk/wangmingyuan/SGS10_analysis/logs
echo -e "${GREEN}  ✓ 目录创建完成${NC}\n"

# ==================== 2. 复制所有肝脏样本数据 ====================
echo -e "${YELLOW}2. 复制肝脏样本数据...${NC}"

# 样本列表
declare -A samples=(
    ["liver3"]="25090559_Liver_3_expression"
    ["liver8"]="25090559_Liver_8_expression"
    ["bla_liver"]="25090559_Bla_liver_1_expression"
    ["con_liver"]="25090559_Con_liver_1_expression"
)

for sample in "${!samples[@]}"; do
    dir_name="${samples[$sample]}"
    src_path="/mnt/BioData/bio/SGS10例/SGS25090559-8例样本数据/data/matrix/${dir_name}"
    dst_path="/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/liver_samples/${dir_name}"
    
    echo -e "  复制 ${sample}..."
    if [ -d "$dst_path" ]; then
        echo -e "    ${YELLOW}目录已存在，跳过${NC}"
    else
        cp -r "$src_path" "$dst_path" 2>/dev/null
        if [ $? -eq 0 ]; then
            echo -e "    ${GREEN}✓ 复制完成${NC}"
        else
            echo -e "    ${RED}✗ 复制失败${NC}"
        fi
    fi
done

# gan56 已经在正确位置
echo -e "  gan56 ${GREEN}✓ 已存在${NC}\n"

# ==================== 3. 解压zip文件 ====================
echo -e "${YELLOW}3. 解压数据文件...${NC}"
cd /mnt/BioDisk/wangmingyuan/SGS10_analysis/data/liver_samples/

for sample in liver3 liver8 bla_liver con_liver; do
    case $sample in
        liver3) dir="25090559_Liver_3_expression"; zip="Liver_3_filtered_feature_bc_matrix.zip" ;;
        liver8) dir="25090559_Liver_8_expression"; zip="Liver_8_filtered_feature_bc_matrix.zip" ;;
        bla_liver) dir="25090559_Bla_liver_1_expression"; zip="Bla_liver_1_filtered_feature_bc_matrix.zip" ;;
        con_liver) dir="25090559_Con_liver_1_expression"; zip="Con_liver_1_filtered_feature_bc_matrix.zip" ;;
    esac
    
    echo -e "  解压 ${sample}..."
    if [ -d "${dir}/${zip%.zip}" ]; then
        echo -e "    ${YELLOW}已解压，跳过${NC}"
    else
        cd "$dir"
        unzip -q "$zip" 2>/dev/null
        if [ $? -eq 0 ]; then
            echo -e "    ${GREEN}✓ 解压完成${NC}"
        else
            echo -e "    ${RED}✗ 解压失败${NC}"
        fi
        cd ..
    fi
done
cd ~

# ==================== 4. 验证数据 ====================
echo -e "\n${YELLOW}4. 验证数据完整性...${NC}"
Rscript -e "
library(Seurat)
cat('验证肝脏样本数据:\n')
samples <- c('liver3', 'liver8', 'bla_liver', 'con_liver')
for(sample in samples) {
  dir <- switch(sample,
    liver3 = '/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/liver_samples/25090559_Liver_3_expression/Liver_3_filtered_feature_bc_matrix',
    liver8 = '/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/liver_samples/25090559_Liver_8_expression/Liver_8_filtered_feature_bc_matrix',
    bla_liver = '/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/liver_samples/25090559_Bla_liver_1_expression/Bla_liver_1_filtered_feature_bc_matrix',
    con_liver = '/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/liver_samples/25090559_Con_liver_1_expression/Con_liver_1_filtered_feature_bc_matrix'
  )
  if(dir.exists(dir)) {
    data <- Read10X(dir)
    cat(sprintf('  %s: ✓ %d 基因, %d 细胞\n', sample, nrow(data), ncol(data)))
  } else {
    cat(sprintf('  %s: ✗ 目录不存在\n', sample))
  }
}
"

# ==================== 5. 更新脚本中的数据路径 ====================
echo -e "\n${YELLOW}5. 更新分析脚本中的数据路径...${NC}"
cd ~/SGS10_analysis/scripts

# 备份原脚本
cp liver_samples_analysis.R liver_samples_analysis.R.bak

# 更新数据路径部分
sed -i '/data_paths <- list/,/)/c\
data_paths <- list(\
  gan56 = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/gan56/gan_56_filtered_feature_bc_matrix",\
  liver3 = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/liver_samples/25090559_Liver_3_expression/Liver_3_filtered_feature_bc_matrix",\
  liver8 = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/liver_samples/25090559_Liver_8_expression/Liver_8_filtered_feature_bc_matrix",\
  bla_liver = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/liver_samples/25090559_Bla_liver_1_expression/Bla_liver_1_filtered_feature_bc_matrix",\
  con_liver = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/liver_samples/25090559_Con_liver_1_expression/Con_liver_1_filtered_feature_bc_matrix"\
)' liver_samples_analysis.R

echo -e "${GREEN}  ✓ 脚本已更新${NC}\n"

# ==================== 6. 创建停止脚本 ====================
echo -e "${YELLOW}6. 创建停止脚本...${NC}"
cat > ~/SGS10_analysis/scripts/stop_liver_samples.sh << 'EOF'
#!/bin/bash
# 停止所有肝脏样本的screen会话

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

samples=("gan56" "liver3" "liver8" "bla_liver" "con_liver")

echo -e "${YELLOW}正在停止所有肝脏样本...${NC}"

for sample in "${samples[@]}"; do
    if screen -ls 2>/dev/null | grep -q "$sample"; then
        screen -S "$sample" -X quit
        echo -e "  ${GREEN}✓ $sample 已停止${NC}"
    else
        echo -e "  ${RED}✗ $sample 未运行${NC}"
    fi
done

echo -e "${GREEN}所有样本已停止！${NC}"
EOF
chmod +x ~/SGS10_analysis/scripts/stop_liver_samples.sh
echo -e "${GREEN}  ✓ 停止脚本已创建${NC}\n"

# ==================== 7. 创建运行脚本（如果不存在） ====================
if [ ! -f ~/SGS10_analysis/scripts/run_liver_samples.sh ]; then
    echo -e "${YELLOW}7. 创建运行脚本...${NC}"
    cat > ~/SGS10_analysis/scripts/run_liver_samples.sh << 'EOF'
#!/bin/bash
# ==================== 运行所有5个肝脏样本 ====================

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
EOF
    chmod +x ~/SGS10_analysis/scripts/run_liver_samples.sh
    echo -e "${GREEN}  ✓ 运行脚本已创建${NC}\n"
fi

# ==================== 8. 创建监控脚本（如果不存在） ====================
if [ ! -f ~/SGS10_analysis/scripts/monitor_liver_samples.sh ]; then
    echo -e "${YELLOW}8. 创建监控脚本...${NC}"
    cat > ~/SGS10_analysis/scripts/monitor_liver_samples.sh << 'EOF'
#!/bin/bash
# ==================== 监控所有肝脏样本进度 ====================

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
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
    echo -e "${BLUE}========================================${NC}"
    echo -e "${BLUE}      肝脏样本完整进度监控 $(date)${NC}"
    echo -e "${BLUE}========================================${NC}\n"
    
    echo -e "${YELLOW}系统资源:${NC}"
    free -h | head -2 | sed 's/^/  /'
    echo ""
    
    printf "${CYAN}%-12s${NC}" "样本"
    for step_info in "${steps[@]}"; do
        step_icon=$(echo $step_info | cut -d':' -f2)
        printf " ${CYAN}%s${NC}" "$step_icon"
    done
    printf " ${CYAN}⏱️ 时间${NC}\n"
    printf "%0.s-" {1..120}
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
        
        # 运行时间
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
        echo ""
    done
    
    printf "%0.s-" {1..120}
    echo -e "\n\n${YELLOW}图例:${NC} ${GREEN}✓${NC} 已完成  ${RED}·${NC} 未完成  ${YELLOW}▶${NC} 运行中"
    echo -e "\n${YELLOW}最新日志:${NC}"
    tail -n 5 /mnt/BioDisk/wangmingyuan/SGS10_analysis/logs/*.log 2>/dev/null | tail -5 || echo "  暂无日志"
    sleep 20
done
EOF
    chmod +x ~/SGS10_analysis/scripts/monitor_liver_samples.sh
    echo -e "${GREEN}  ✓ 监控脚本已创建${NC}\n"
fi

# ==================== 9. 清理旧检查点 ====================
echo -e "${YELLOW}9. 清理旧的检查点文件...${NC}"
rm -f /mnt/BioDisk/wangmingyuan/SGS10_analysis/*/*_checkpoint.txt
echo -e "${GREEN}  ✓ 检查点已清理${NC}\n"

# ==================== 10. 启动分析 ====================
echo -e "${YELLOW}10. 准备启动分析...${NC}"
read -p "是否现在启动所有5个肝脏样本分析？(y/n): " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    cd ~/SGS10_analysis/scripts
    ./run_liver_samples.sh
    echo -e "\n${GREEN}启动完成！使用以下命令监控进度：${NC}"
    echo "  cd ~/SGS10_analysis/scripts"
    echo "  ./monitor_liver_samples.sh"
else
    echo -e "\n${YELLOW}您可以稍后手动启动：${NC}"
    echo "  cd ~/SGS10_analysis/scripts"
    echo "  ./run_liver_samples.sh"
    echo "  ./monitor_liver_samples.sh"
fi

echo -e "\n${BLUE}========================================${NC}"
echo -e "${GREEN}✅ 设置完成！${NC}"
echo -e "${BLUE}========================================${NC}"
