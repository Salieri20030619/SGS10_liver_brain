#!/bin/bash
echo "============================================================"
echo "   RCTD 完整状态检查  [$(date +%H:%M:%S)]"
echo "============================================================"

echo ""
echo "【1. 概览】"
bash ~/SGS10_analysis/scripts/watch_rctd.sh -s

echo ""
echo "【2. R 进程】"
ps -u wangmingyuan -o pid,pcpu,stat,time,cmd --sort=-pcpu | grep -E "R --|Rscript" | head -10 | awk '{
  cmd_type = "unknown"
  for (i=5; i<=NF; i++) {
    if ($i ~ /workRSOCK/) { cmd_type = "worker"; break }
    if ($i ~ /rctd_batch/) { cmd_type = "main"; break }
  }
  printf "  PID=%-8s CPU=%-6s STAT=%-5s TIME=%-12s TYPE=%s\n", $1, $2, $3, $4, cmd_type
}'

echo ""
echo "【3. R 进程总 CPU】"
ps -u wangmingyuan -o pcpu,cmd --no-headers | grep -E "R --|Rscript" | awk '
  {sum+=$1; n++}
  END {printf "  进程数=%d, 总CPU=%.0f%% (= %.0f 核等价)\n", n, sum, sum/100}'

echo ""
echo "【4. 负载 / 内存】"
uptime | awk -F'load average:' '{print "  load:" $2}'
free -h | awk 'NR==2 {printf "  内存: 用 %s / 可用 %s\n", $3, $7}'

echo ""
echo "【5. 输出文件】"
echo "  Reference:"
ls -lh ~/SGS10_analysis/spatial_results/rctd/ref_*.rds 2>/dev/null | awk '{print "    " $5 "  " $NF}'
echo "  Per-sample:"
ls -lh ~/SGS10_analysis/spatial_results/rctd/per_sample/*_rctd_results.rds 2>/dev/null | awk '{print "    " $5 "  " $NF}' || echo "    (还没有)"

echo ""
echo "【6. Log 关键标记】"
grep -E "=== RCTD on|✓ done|Building|saving" ~/SGS10_analysis/scripts/logs/rctd_batch.log 2>/dev/null | tail -10 | awk '{print "  " $0}'

echo ""
echo "【7. Log 最新 15 行】"
tail -15 ~/SGS10_analysis/scripts/logs/rctd_batch.log 2>/dev/null | awk '{print "  " $0}'

echo ""
echo "============================================================"
