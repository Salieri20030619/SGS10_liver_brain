#!/bin/bash
# =============================================================================
# 启动 RCTD 100% 重跑 (后台 nohup, 12-20 小时)
# 
# 改动 vs 原 rctd_batch_10samples.R:
#   - SUBSAMPLE_RATE: 0.2 → 1.0 (不抽样, 全部 spots)
#   - MAX_CORES: 8 → 4 (系统 load 244, 不要跟 gromacs 抢太凶)
#   - OUT_DIR: spatial_results/rctd → spatial_results/rctd_full (避免覆盖)
#
# 后台运行 + 失败重试 + nice 19 (低优先级)
# =============================================================================

PR="/mnt/BioDisk/wangmingyuan/SGS10_analysis"
cd "$PR"

# ---------- 1. 创建重跑脚本 ----------
SCRIPT_OUT="$PR/scripts/rctd_batch_10samples_FULL.R"

# 复制原脚本, 改关键参数
sed -e 's|SUBSAMPLE_RATE <- 0.2|SUBSAMPLE_RATE <- 1.0|' \
    -e 's|MAX_CORES <- 8|MAX_CORES <- 4|' \
    -e 's|OUT_DIR <- file.path(PROJECT_ROOT, "spatial_results", "rctd")|OUT_DIR <- file.path(PROJECT_ROOT, "spatial_results", "rctd_full")|' \
    "$PR/scripts/rctd_batch_10samples.R" > "$SCRIPT_OUT"

echo "Created: $SCRIPT_OUT"
echo ""
echo "===== 改动校验 ====="
echo "--- SUBSAMPLE_RATE: ---"
grep "SUBSAMPLE_RATE <-" "$SCRIPT_OUT" | head -3
echo "--- MAX_CORES: ---"
grep "MAX_CORES <-" "$SCRIPT_OUT" | head -3
echo "--- OUT_DIR: ---"
grep 'OUT_DIR <- file.path' "$SCRIPT_OUT" | head -3

echo ""
echo "===== 启动准备 ====="

# ---------- 2. 创建输出目录 ----------
mkdir -p "$PR/spatial_results/rctd_full/per_sample"
mkdir -p "$PR/spatial_results/rctd_full/enriched"
mkdir -p "$PR/spatial_results/rctd_full/summary"
mkdir -p "$PR/scripts/logs"

# ---------- 3. 复用现有 reference (软链接, 不重跑 reference 构建) ----------
# Reference 与抽样无关, 直接复用 20% 版的
if [ -f "$PR/spatial_results/rctd/ref_brain_22ct.rds" ]; then
    ln -sf "$PR/spatial_results/rctd/ref_brain_22ct.rds" "$PR/spatial_results/rctd_full/ref_brain_22ct.rds"
    echo "✓ Linked brain reference"
fi
if [ -f "$PR/spatial_results/rctd/ref_liver_33ct.rds" ]; then
    ln -sf "$PR/spatial_results/rctd/ref_liver_33ct.rds" "$PR/spatial_results/rctd_full/ref_liver_33ct.rds"
    echo "✓ Linked liver reference"
fi

# ---------- 4. 显示启动命令 (不直接启动, 要你确认) ----------
LOG_FILE="$PR/scripts/logs/rctd_full_$(date +%Y%m%d_%H%M%S).log"

echo ""
echo "==============================================================================="
echo "启动命令 (检查后用 nohup 后台跑):"
echo "==============================================================================="
cat <<EOF

cd $PR
nohup nice -n 19 Rscript scripts/rctd_batch_10samples_FULL.R \\
    > $LOG_FILE 2>&1 &
echo "PID=\$!" > $PR/scripts/logs/rctd_full_pid.txt
echo "RCTD 100% 重跑已启动, PID=\$!"
echo "查看日志: tail -f $LOG_FILE"
echo "查看进度: grep '^\[' $LOG_FILE | tail -20"
echo "查看 PID: cat $PR/scripts/logs/rctd_full_pid.txt"
echo "杀掉:    kill \$(cat $PR/scripts/logs/rctd_full_pid.txt | grep -oP '\\d+')"

EOF

echo "==============================================================================="
echo "如果决定启动, 直接复制上面的 cd ... nohup ... 命令运行"
echo "==============================================================================="
echo ""
echo "启动前请最后核对:"
echo "  - 当前 load average 是否过高 (>250 不建议启动): $(uptime | awk -F'load average:' '{print $2}')"
echo "  - 磁盘空间足够 (>50GB free): $(df -h $PR | tail -1)"
echo "  - 系统时间充裕 (12-20 小时不会重启): $(uptime)"
