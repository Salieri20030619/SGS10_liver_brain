#!/bin/bash
# ============================================================
# verify_three_axes.sh
# 交叉验证你提出的 3 条肝脑机制轴在我的 Stage 1 CSV 里的数字
# ============================================================
# 用法:
#   cd ~/SGS10_analysis
#   bash scripts/verify_three_axes.sh
# ============================================================

cd /mnt/BioDisk/wangmingyuan/SGS10_analysis

OUT=/tmp/verify_three_axes_output.txt
rm -f $OUT
exec > >(tee -a $OUT) 2>&1

echo "============================================================"
echo "验证 3 条肝脑机制轴 vs Stage 1 CSV"
echo "$(date '+%Y-%m-%d %H:%M:%S')"
echo "============================================================"
echo ""

# ============================================================
# 轴 A: Efferocytosis 小胶质细胞吞噬功能受损
# ============================================================
echo "########## 轴 A: Efferocytosis ##########"
echo "声称: 脑 Microglia_Homeostatic tumor_vs_control 里"
echo "  Elmo1 FC=-0.70, padj=6.1e-19"
echo "  Dock4 FC=-0.67, padj=2.5e-10"
echo "  Tlr7  FC=-0.81"
echo "  Ifitm10 FC=-1.26"
echo ""

CSV_A="diff_stage1_results/deg_v4_refined/brain/brain_tumor_vs_control_Microglia_Homeostatic.csv"
if [ -f "$CSV_A" ]; then
  echo "## [1.1] 查: Elmo1 / Dock4 / Tlr7 / Ifitm10 / Il6st 在 $(basename $CSV_A)"
  head -1 "$CSV_A"
  grep -E "^\"(Elmo1|Dock4|Tlr7|Ifitm10|Il6st|Mertk|Gas6|Axl|Rac1)\"," "$CSV_A" || \
    echo "(none of these genes found)"

  echo ""
  echo "## [1.2] 该 CSV 的 Top 30 行 (按 p_val_adj 排序, 已是默认)"
  head -31 "$CSV_A"

  echo ""
  echo "## [1.3] 统计该 CSV 中 Up 和 Down 基因数"
  awk -F',' 'NR>1 {gsub(/"/, "", $NF); print $NF}' "$CSV_A" | sort | uniq -c

  echo ""
  echo "## [1.4] 统计 padj<0.05 的 down 基因数 (严格阈值, 看是否接近 '17 个下调 1 个上调')"
  awk -F',' 'NR>1 {
    padj=$6; gsub(/"/, "", padj);
    fc=$2; gsub(/"/, "", fc);
    sig=$NF; gsub(/"/, "", sig);
    if (padj+0 < 0.05 && fc+0 < -0.5) n_down++;
    if (padj+0 < 0.05 && fc+0 > 0.5)  n_up++;
  } END {printf "  padj<0.05 & |FC|>0.5: Up=%d, Down=%d\n", n_up+0, n_down+0}' "$CSV_A"
else
  echo "  ⚠ FILE NOT FOUND: $CSV_A"
fi
echo ""
echo ""

# ============================================================
# 轴 B: Thyroid hormone signaling
# ============================================================
echo "########## 轴 B: Thyroid hormone ##########"
echo "声称:"
echo "  shared KEGG Thyroid hormone signaling pathway"
echo "    Liver p.adjust=1.82e-07, Brain p.adjust=2.54e-03"
echo "    combined_score=9.33"
echo "  脑 Endothelial 里 Thra FC=+1.11"
echo ""

CSV_B="diff_stage1_results/deg_v4_refined/brain/brain_tumor_vs_control_Endothelial.csv"
if [ -f "$CSV_B" ]; then
  echo "## [2.1] 查: Thra / Thrb / Dio1 / Dio2 / Dio3 在 $(basename $CSV_B)"
  head -1 "$CSV_B"
  grep -E "^\"(Thra|Thrb|Dio1|Dio2|Dio3|Thrsp|Thrap3)\"," "$CSV_B" || \
    echo "(none of these genes found)"
else
  echo "  ⚠ FILE NOT FOUND: $CSV_B"
fi

echo ""
echo "## [2.2] shared KEGG 中查 thyroid/hormone"
SHARED="diff_stage1_results/viz/07_kegg_venn/shared_pathways.csv"
if [ -f "$SHARED" ]; then
  echo "  -- header --"
  head -1 "$SHARED"
  echo "  -- 匹配 thyroid/hormone --"
  grep -iE "thyroid|hormone|endocrine" "$SHARED" || echo "(no match)"

  echo ""
  echo "## [2.3] shared_pathways.csv 的完整 rank 列表 (全部 44 条)"
  cat "$SHARED"
else
  echo "  ⚠ FILE NOT FOUND: $SHARED"
fi
echo ""
echo ""

# ============================================================
# 轴 C: BBB 代偿强化 / 脑内皮免疫激活
# ============================================================
echo "########## 轴 C: BBB 代偿 + 脑内皮免疫激活 ##########"
echo "声称:"
echo "  脑 Endothelial tumor_vs_control:"
echo "    Cldn5   FC=+1.53, padj=8.49e-10"
echo "    Ctnnb1  FC=+0.91, padj=9.51e-03"
echo "    Cdh5    FC=+0.69, padj=4.62e-05"
echo "    Flt1    FC=+0.36"
echo "    Cxcl12  FC=+0.87"
echo "    H2-D1 / H2-K1 上调"
echo "  shared KEGG Adherens junction:"
echo "    Liver p=6.99e-08, Brain p=9.24e-03"
echo ""

if [ -f "$CSV_B" ]; then
  echo "## [3.1] 查: Cldn5 / Ctnnb1 / Cdh5 / Flt1 / Cxcl12 / H2-D1 / H2-K1 在 $(basename $CSV_B)"
  grep -E "^\"(Cldn5|Ctnnb1|Cdh5|Flt1|Cxcl12|H2-D1|H2-K1|H2-Q7|Ocln|Tjp1|Tjp2)\"," "$CSV_B" || \
    echo "(none found)"

  echo ""
  echo "## [3.2] 脑 Endothelial 的 Top 20 上调 DEG"
  awk -F',' 'NR>1 {
    fc=$2; gsub(/"/, "", fc);
    sig=$NF; gsub(/"/, "", sig);
    if (sig == "Up") print
  }' "$CSV_B" | head -20

  echo ""
  echo "## [3.3] 统计 Up vs Down"
  awk -F',' 'NR>1 {gsub(/"/, "", $NF); print $NF}' "$CSV_B" | sort | uniq -c
fi
echo ""

echo "## [3.4] shared KEGG 中查 Adherens / junction / cadherin"
if [ -f "$SHARED" ]; then
  grep -iE "adheren|junction|cadherin|tight" "$SHARED" || echo "(no match)"
fi
echo ""
echo ""

# ============================================================
# 附加: Microglia_Homeostatic 里 Il6st 的精确值 (核对之前的 hit)
# ============================================================
echo "########## 附加: 核对 Il6st 原 hit ##########"
if [ -f "$CSV_A" ]; then
  echo "## [4.1] Il6st in Microglia_Homeostatic"
  grep "^\"Il6st\"," "$CSV_A"
fi
echo ""

# ============================================================
# 完成
# ============================================================
echo "============================================================"
echo "验证完成"
echo "输出同时保存到: $OUT"
echo "============================================================"
