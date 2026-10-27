#!/bin/bash
# 三任务全面诊断: Phase 3 v3 + RCTD 100% + v8 全量
PR="/mnt/BioDisk/wangmingyuan/SGS10_analysis"
cd "$PR"

NOW=$(date +%H:%M:%S)
echo "==============================================================================="
echo "三任务进度诊断 @ $NOW"
echo "==============================================================================="

# ---------- 0. 系统资源 ----------
echo ""
echo "===== 0. 系统资源 ====="
echo "Load: $(uptime | awk -F'load average:' '{print $2}')"
free -h | head -2

# ---------- 1. R 进程 (按 CPU 时间排序) ----------
echo ""
echo "===== 1. R 进程 (Rscript 主进程) ====="
ps aux | grep wangmingyuan | grep "Rscript" | grep -v grep | grep -v RSOCK | \
    awk '{printf "PID=%s 启动=%s CPU时长=%s CPU%%=%s%% RAM%%=%s%% CMD=%s\n", \
                  $2, $9, $10, $3, $4, $NF}' | head -10

n_r=$(ps aux | grep wangmingyuan | grep Rscript | grep -v grep | grep -v RSOCK | wc -l)
echo ""
echo "总主 R 进程: $n_r (不算 RSOCK worker)"

# ---------- 2. Phase 3 v3 marker ----------
echo ""
echo "==============================================================================="
echo "===== [1/3] Phase 3 v3 微环境 marker 基因 ====="
echo "==============================================================================="

P3_DIR="$PR/D1_v2_RCTD_composition/microenv_markers"
echo ""
echo "--- 输出文件 (按修改时间排序) ---"
ls -lat "$P3_DIR/" 2>&1 | head -15

echo ""
echo "--- 关键 csv 是否生成 ---"
for f in liver_microenv_markers_all.csv \
         liver_microenv_markers_top20.csv \
         liver_M05_tumor_specific_markers.csv \
         brain_microenv_markers_all.csv \
         brain_microenv_markers_top20.csv \
         brain_microenv_in_tumor_markers.csv; do
    if [ -f "$P3_DIR/$f" ]; then
        rows=$(($(wc -l < "$P3_DIR/$f") - 1))
        echo "  ✓ $f ($rows rows)"
    else
        echo "  ✗ $f (未生成)"
    fi
done

echo ""
echo "--- v3_2 log 末尾 30 行 ---"
P3_LOG="$PR/scripts/logs/microenv_marker_v3_2.log"
if [ -f "$P3_LOG" ]; then
    echo "Log: $P3_LOG ($(stat -c '%y' $P3_LOG | cut -d. -f1))"
    echo ""
    tail -n 30 "$P3_LOG"
else
    echo "v3_2 log 不存在!"
    ls -la "$PR/scripts/logs/"*marker*.log 2>&1
fi

# ---------- 3. RCTD 100% ----------
echo ""
echo "==============================================================================="
echo "===== [2/3] RCTD 100% 重跑 ====="
echo "==============================================================================="

RCTD_FULL="$PR/spatial_results/rctd_full"
if [ -d "$RCTD_FULL" ]; then
    n_done=$(ls "$RCTD_FULL/per_sample/"*_rctd_results.rds 2>/dev/null | wc -l)
    echo "已完成 sample: $n_done / 10"
    
    if [ "$n_done" -gt 0 ]; then
        echo ""
        echo "--- 已完成 sample 列表 ---"
        ls "$RCTD_FULL/per_sample/"*_rctd_results.rds 2>/dev/null | \
            xargs -I {} basename {} _rctd_results.rds | sed 's/^/  ✓ /'
    fi
    
    echo ""
    echo "--- per_sample 目录 ---"
    ls -lat "$RCTD_FULL/per_sample/" 2>&1 | head -8
fi

echo ""
echo "--- RCTD log 末尾 25 行 ---"
RCTD_LOG=$(ls -t "$PR/scripts/logs/rctd_full"*.log 2>/dev/null | head -1)
if [ -n "$RCTD_LOG" ]; then
    echo "Log: $RCTD_LOG"
    echo "Size: $(ls -la $RCTD_LOG | awk '{print $5}') bytes"
    echo "Last modified: $(stat -c '%y' $RCTD_LOG | cut -d. -f1)"
    echo ""
    tail -n 25 "$RCTD_LOG"
else
    echo "(no rctd_full log found)"
fi

# ---------- 4. v8 全量 ----------
echo ""
echo "==============================================================================="
echo "===== [3/3] v8 客户友好版全量 ====="
echo "==============================================================================="

V8_DIR="$PR/taskB_he_overlay_v8"
if [ -d "$V8_DIR" ]; then
    n_done_v8=$(ls "$V8_DIR/"*/E_hex_microenv_on_HE.png 2>/dev/null | wc -l)
    echo "已完成 sample (panel E): $n_done_v8 / 10"
    
    if [ "$n_done_v8" -gt 0 ]; then
        echo ""
        echo "--- 已完成 sample 目录 ---"
        for d in "$V8_DIR"/*/; do
            sname=$(basename "$d")
            n_png=$(ls "$d"*.png 2>/dev/null | wc -l)
            echo "  $sname ($n_png panels)"
        done
    fi
    
    echo ""
    echo "--- run_full log 末尾 15 行 ---"
    if [ -f "$V8_DIR/run_full.log" ]; then
        tail -n 15 "$V8_DIR/run_full.log"
    elif [ -f "$V8_DIR/test_Liver_3.log" ]; then
        echo "(只有测试 log, 全量未启动)"
        tail -n 5 "$V8_DIR/test_Liver_3.log"
    else
        echo "(无 log)"
    fi
else
    echo "v8 目录不存在"
fi

echo ""
echo "==============================================================================="
echo "完成"
echo "==============================================================================="
