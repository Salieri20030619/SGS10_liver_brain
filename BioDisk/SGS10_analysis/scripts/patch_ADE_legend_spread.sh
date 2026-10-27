#!/usr/bin/env bash
# patch_ADE_legend_spread.sh —— A/D/E 图例横向平铺 (两边顶到头)
# 前提: 已运行过 patch_ADE_legend_bottom.sh
# 用法: bash scripts/patch_ADE_legend_spread.sh
set -e
F=/mnt/BioDisk/wangmingyuan/SGS10_analysis/scripts/he_overlay_RCTD100.R

cp -n "$F" "${F}.bak_spread" || true
if grep -q 'LEGEND_SPREAD' "$F"; then
  echo "已含 LEGEND_SPREAD, 跳过。重打请先: cp ${F}.bak_spread $F"
  exit 0
fi

# ---- Panel A: ncol 6→8, 加 byrow+title.position ----
sed -i 's/alpha = 1), ncol = 6))/alpha = 1), ncol = 8, byrow = TRUE, title.position = "top"))/' "$F"

# ---- Panel D: ncol 3→6, 加 byrow+title.position ----
sed -i 's/size = 3, alpha = 1), ncol = 3))/size = 3, alpha = 1), ncol = 6, byrow = TRUE, title.position = "top"))/' "$F"

# ---- Panel E: ncol 3→6, 加 byrow+title.position ----
sed -i 's/list(alpha = 1), ncol = 3))/list(alpha = 1), ncol = 6, byrow = TRUE, title.position = "top"))/' "$F"

# ---- A/D/E theme: 加 legend.justification 左对齐顶边 ----
sed -i 's/theme(legend.position = "bottom")  # LEGEND_BOTTOM_A/theme(legend.position = "bottom", legend.justification = c(0, 0.5))  # LEGEND_BOTTOM_A LEGEND_SPREAD/' "$F"
sed -i 's/theme(legend.position = "bottom")  # LEGEND_BOTTOM_DE/theme(legend.position = "bottom", legend.justification = c(0, 0.5))  # LEGEND_BOTTOM_DE LEGEND_SPREAD/' "$F"
# pE 那行末尾无注释, 用 $ 锚定
sed -i 's/theme(legend.position = "bottom")$/theme(legend.position = "bottom", legend.justification = c(0, 0.5))/' "$F"

echo "=== 补丁完成, 确认 ==="
grep -n 'byrow\|title.position\|justification\|LEGEND_SPREAD' "$F"
