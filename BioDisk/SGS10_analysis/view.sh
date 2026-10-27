#!/bin/bash
# 查看 Phase 2 微环境关键细胞分析的全部结果
PR="/mnt/BioDisk/wangmingyuan/SGS10_analysis"
P2_DIR="$PR/D1_v2_RCTD_composition/microenv_key_celltypes"

echo "==============================================================================="
echo "Phase 2 微环境关键细胞结果查看"
echo "==============================================================================="

# ---------- 输出文件清单 ----------
echo ""
echo "===== 输出文件 ====="
ls -la "$P2_DIR/"

# ---------- 1. summary.md ----------
echo ""
echo "===== 1. summary.md (高层总结) ====="
cat "$P2_DIR/summary.md"

# ---------- 2. Liver top5 ----------
echo ""
echo "===== 2. LIVER 微环境 top 5 关键细胞 ====="
echo ""
cat "$P2_DIR/liver_microenv_celltypes_top5.csv" | column -t -s ',' | head -45

# ---------- 3. Brain top5 ----------
echo ""
echo "===== 3. BRAIN 微环境 top 5 关键细胞 ====="
echo ""
cat "$P2_DIR/brain_microenv_celltypes_top5.csv" | column -t -s ',' | head -40

# ---------- 4. Liver tumor vs baseline 差异 ----------
echo ""
echo "===== 4. LIVER tumor vs baseline 差异 (top 3 per microenv by abs log2FC) ====="
echo ""
echo "(看 tumor 状态下哪些细胞在每微环境中变化最大)"
python3 -c "
import csv
with open('$P2_DIR/liver_microenv_tumor_vs_baseline.csv') as f:
    rows = list(csv.DictReader(f))
# 按 microenv 分组, 取 abs log2FC top 3
from collections import defaultdict
g = defaultdict(list)
for r in rows:
    g[r['microenv_id']].append(r)
for me_id in sorted(g.keys()):
    rs = g[me_id]
    label = rs[0]['microenv_label']
    print(f'\n--- {me_id} {label} ---')
    rs_sorted = sorted(rs, key=lambda x: abs(float(x.get('log2FC_tumor_vs_baseline', 0))), reverse=True)[:5]
    for r in rs_sorted:
        log2fc = float(r['log2FC_tumor_vs_baseline'])
        sign = '↑ tumor' if log2fc > 0 else '↓ tumor'
        print(f'  {r[\"celltype\"]:25s}  baseline={float(r[\"baseline\"]):.3f}  tumor={float(r[\"tumor\"]):.3f}  log2FC={log2fc:+.2f}  {sign}')
" 2>&1

# ---------- 5. Brain tumor vs baseline ----------
echo ""
echo "===== 5. BRAIN tumor vs baseline 差异 ====="
python3 -c "
import csv
with open('$P2_DIR/brain_microenv_tumor_vs_baseline.csv') as f:
    rows = list(csv.DictReader(f))
from collections import defaultdict
g = defaultdict(list)
for r in rows:
    g[r['microenv_id']].append(r)
for me_id in sorted(g.keys()):
    rs = g[me_id]
    label = rs[0]['microenv_label']
    print(f'\n--- {me_id} {label} ---')
    rs_sorted = sorted(rs, key=lambda x: abs(float(x.get('log2FC_tumor_vs_baseline', 0))), reverse=True)[:5]
    for r in rs_sorted:
        log2fc = float(r['log2FC_tumor_vs_baseline'])
        sign = '↑ tumor' if log2fc > 0 else '↓ tumor'
        print(f'  {r[\"celltype\"]:25s}  baseline={float(r[\"baseline\"]):.3f}  tumor={float(r[\"tumor\"]):.3f}  log2FC={log2fc:+.2f}  {sign}')
" 2>&1

echo ""
echo "==============================================================================="
echo "完成"
echo "==============================================================================="
