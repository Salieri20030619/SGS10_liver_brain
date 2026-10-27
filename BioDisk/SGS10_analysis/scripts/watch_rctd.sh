#!/bin/bash
# ============================================================
# watch_rctd.sh - RCTD batch 进度监视脚本
# ============================================================
# 用法:
#   一次性查看:  bash watch_rctd.sh
#   持续刷新:    bash watch_rctd.sh -w         (每 60s 刷新)
#   自定义间隔:  bash watch_rctd.sh -w 30      (每 30s 刷新)
#   只看摘要:    bash watch_rctd.sh -s         (简洁版, 单屏输出)
# ============================================================

PROJECT_ROOT="$HOME/SGS10_analysis"
LOG_FILE="$PROJECT_ROOT/scripts/logs/rctd_batch.log"
OUT_DIR="$PROJECT_ROOT/spatial_results/rctd"

# 10 个 spatial 样本
SAMPLES=(
  "Liver_3" "Liver_8" "gan56" "Bla_liver_1" "Con_liver_1"
  "Brain_3" "Brain_8" "nao56" "Bla_brain_1" "Con_brain_1"
)
TOTAL=${#SAMPLES[@]}

# 颜色
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
BOLD='\033[1m'
RESET='\033[0m'

# ============================================================
# 单次状态输出
# ============================================================
show_status() {
  clear
  echo -e "${BOLD}${CYAN}========================================${RESET}"
  echo -e "${BOLD}${CYAN}   RCTD Batch 监视  [$(date +%H:%M:%S)]${RESET}"
  echo -e "${BOLD}${CYAN}========================================${RESET}"
  
  # --- 1. 服务器负载 ---
  echo ""
  echo -e "${BOLD}📊 服务器负载${RESET}"
  uptime | awk -F'load average:' '{print "  " $2}'
  free -h | awk 'NR==2 {printf "  内存: 用 %s / 总 %s, 可用 %s\n", $3, $2, $7}'
  
  # --- 2. RCTD 进程 ---
  echo ""
  echo -e "${BOLD}🏃 RCTD 进程${RESET}"
  
  # 所有 R 进程 (Rscript 主 + workRSOCK worker)
  R_PROCS=$(ps -u $USER -o pid,pcpu,stat,time,cmd --no-headers --sort=-pcpu 2>/dev/null | \
            grep -E "R --|Rscript" | grep -v grep)
  
  if [ -z "$R_PROCS" ]; then
    echo -e "  ${RED}✗ R 进程不存在 (已结束或未启动)${RESET}"
  else
    # 显示 top 5 R 进程
    echo "$R_PROCS" | head -5 | awk '{
      cmd_type = "other"
      for (i=5; i<=NF; i++) {
        if ($i ~ /workRSOCK/) { cmd_type = "worker"; break }
        if ($i ~ /rctd_batch/) { cmd_type = "main"; break }
      }
      printf "  PID=%-8s CPU=%-6s%% STAT=%-5s TIME=%-12s TYPE=%s\n", $1, $2, $3, $4, cmd_type
    }'
    
    # 总 CPU
    total_cpu=$(echo "$R_PROCS" | awk '{sum+=$2} END{printf "%.0f", sum}')
    n_proc=$(echo "$R_PROCS" | wc -l)
    cores=$(awk "BEGIN{printf \"%.0f\", $total_cpu/100}")
    echo -e "  ${CYAN}总计: $n_proc 进程, 总 CPU=${total_cpu}% (= $cores 核等价)${RESET}"
  fi
  
  # Screen 状态
  SCREEN_STATUS=$(screen -ls 2>/dev/null | grep rctd_batch)
  if [ -n "$SCREEN_STATUS" ]; then
    echo -e "  ${GREEN}Screen:${RESET}"
    echo "$SCREEN_STATUS" | awk '{print "    " $0}'
  fi
  
  # --- 3. Reference 构建状态 ---
  echo ""
  echo -e "${BOLD}📚 Reference${RESET}"
  if [ -f "$OUT_DIR/ref_brain_22ct.rds" ]; then
    size=$(du -h "$OUT_DIR/ref_brain_22ct.rds" | awk '{print $1}')
    t=$(stat -c %y "$OUT_DIR/ref_brain_22ct.rds" | cut -d. -f1)
    echo -e "  ${GREEN}✓ Brain ref${RESET}  ($size, saved $t)"
  else
    echo -e "  ${YELLOW}⏳ Brain ref${RESET}  (还在构建...)"
  fi
  
  if [ -f "$OUT_DIR/ref_liver_33ct.rds" ]; then
    size=$(du -h "$OUT_DIR/ref_liver_33ct.rds" | awk '{print $1}')
    t=$(stat -c %y "$OUT_DIR/ref_liver_33ct.rds" | cut -d. -f1)
    echo -e "  ${GREEN}✓ Liver ref${RESET}  ($size, saved $t)"
  else
    echo -e "  ${YELLOW}⏳ Liver ref${RESET}  (还在构建...)"
  fi
  
  # --- 4. 样本完成情况 ---
  echo ""
  echo -e "${BOLD}🔬 10 样本进度${RESET}"
  completed=0
  running=""
  
  for s in "${SAMPLES[@]}"; do
    if [ -f "$OUT_DIR/per_sample/${s}_rctd_results.rds" ]; then
      size=$(du -h "$OUT_DIR/per_sample/${s}_rctd_results.rds" 2>/dev/null | awk '{print $1}')
      t=$(stat -c %y "$OUT_DIR/per_sample/${s}_rctd_results.rds" 2>/dev/null | cut -d. -f1 | cut -d' ' -f2)
      echo -e "  ${GREEN}✓${RESET} $s  ($size, done $t)"
      completed=$((completed + 1))
    else
      # 判断是不是正在跑 (log 里有没有这个 sample 的"RCTD on")
      if [ -f "$LOG_FILE" ]; then
        latest_sample=$(grep "=== RCTD on" "$LOG_FILE" 2>/dev/null | tail -1 | sed 's/.*RCTD on //; s/ (.*//')
        if [ "$latest_sample" = "$s" ]; then
          echo -e "  ${YELLOW}⏳${RESET} $s  ${CYAN}(正在跑)${RESET}"
          running="$s"
        else
          echo -e "  ${RED}○${RESET} $s  (等待)"
        fi
      else
        echo -e "  ${RED}○${RESET} $s  (等待)"
      fi
    fi
  done
  
  echo ""
  echo -e "  ${BOLD}进度: $completed / $TOTAL${RESET}  $(printf '[%-10s]' $(printf '=%.0s' $(seq 1 $((completed * 10 / TOTAL)))))"
  
  # --- 5. 预估剩余时间 ---
  if [ $completed -ge 1 ] && [ $completed -lt $TOTAL ]; then
    echo ""
    # 看日志第一个 "✓ done" 和最后一个 "✓ done" 时间差, 算每样本平均时间
    first_done=$(grep "✓ done" "$LOG_FILE" 2>/dev/null | head -1 | grep -oP '\[\d+:\d+:\d+\]' | tr -d '[]')
    last_done=$(grep "✓ done" "$LOG_FILE" 2>/dev/null | tail -1 | grep -oP '\[\d+:\d+:\d+\]' | tr -d '[]')
    
    if [ -n "$first_done" ] && [ -n "$last_done" ] && [ "$first_done" != "$last_done" ] && [ $completed -ge 2 ]; then
      # 计算秒数差
      first_s=$(echo $first_done | awk -F: '{print $1*3600+$2*60+$3}')
      last_s=$(echo $last_done | awk -F: '{print $1*3600+$2*60+$3}')
      diff_s=$((last_s - first_s))
      avg_s=$((diff_s / (completed - 1)))
      
      remaining=$((TOTAL - completed))
      est_s=$((avg_s * remaining))
      est_h=$((est_s / 3600))
      est_m=$(((est_s % 3600) / 60))
      
      echo -e "  ${BLUE}⏱  平均每样本: $((avg_s / 60)) min${RESET}"
      echo -e "  ${BLUE}⏱  预估剩余: ${est_h}h ${est_m}min${RESET}"
    fi
  fi
  
  # --- 6. 最近 log ---
  echo ""
  echo -e "${BOLD}📜 Log 最后 15 行${RESET}"
  if [ -f "$LOG_FILE" ]; then
    tail -15 "$LOG_FILE" | awk '{print "  " $0}'
  else
    echo -e "  ${RED}(log 文件不存在: $LOG_FILE)${RESET}"
  fi
  
  # --- 7. Error 检测 ---
  echo ""
  if [ -f "$LOG_FILE" ]; then
    error_count=$(grep -c "错误\|Error\|FAILED\|✗" "$LOG_FILE" 2>/dev/null)
    if [ "$error_count" -gt 0 ]; then
      echo -e "${BOLD}${RED}⚠️  发现 $error_count 条错误/警告${RESET}"
      echo "  最近 5 条:"
      grep -E "错误|Error|FAILED|✗" "$LOG_FILE" | tail -5 | awk '{print "    " $0}'
    else
      echo -e "${GREEN}✓ 无错误${RESET}"
    fi
  fi
  
  echo ""
}

# ============================================================
# 简洁一行版 (适合快速查)
# ============================================================
show_summary() {
  completed=0
  for s in "${SAMPLES[@]}"; do
    [ -f "$OUT_DIR/per_sample/${s}_rctd_results.rds" ] && completed=$((completed + 1))
  done
  
  RCTD_PROC=$(ps aux | grep "Rscript.*rctd_batch" | grep -v grep | wc -l)
  brain_ref=$([ -f "$OUT_DIR/ref_brain_22ct.rds" ] && echo "✓" || echo "⏳")
  liver_ref=$([ -f "$OUT_DIR/ref_liver_33ct.rds" ] && echo "✓" || echo "⏳")
  load=$(uptime | awk -F'load average:' '{print $2}' | awk -F, '{print $1}' | tr -d ' ')
  
  state="🏃跑中"
  [ "$RCTD_PROC" -eq 0 ] && state="💤结束"
  
  echo "[$(date +%H:%M:%S)] $state | Brain$brain_ref Liver$liver_ref | 样本 $completed/10 | load=$load"
  
  # 最后一行 log
  if [ -f "$LOG_FILE" ]; then
    tail -1 "$LOG_FILE" | awk '{print "  └─ " $0}'
  fi
}

# ============================================================
# 主逻辑
# ============================================================
case "$1" in
  -w|--watch)
    interval=${2:-60}
    echo "监视中 (每 ${interval}s 刷新, Ctrl+C 退出)..."
    sleep 2
    while true; do
      show_status
      sleep $interval
    done
    ;;
  -s|--summary)
    show_summary
    ;;
  -h|--help)
    head -15 "$0" | tail -11
    ;;
  *)
    show_status
    ;;
esac
