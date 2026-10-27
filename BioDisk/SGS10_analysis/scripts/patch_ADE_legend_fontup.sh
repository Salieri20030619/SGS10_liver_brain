#!/usr/bin/env bash
# patch_ADE_legend_fontup.sh —— 图例字号+2, 加列间距, D/E改居中
# 前提: 已运行 patch_ADE_legend_spread.sh
# 用法: bash scripts/patch_ADE_legend_fontup.sh
set -e
F=/mnt/BioDisk/wangmingyuan/SGS10_analysis/scripts/he_overlay_RCTD100.R

cp -n "$F" "${F}.bak_fontup" || true
if grep -q 'LEGEND_FONTUP' "$F"; then
  echo "已含 LEGEND_FONTUP, 跳过。重打请先: cp ${F}.bak_fontup $F"
  exit 0
fi

# 1) plot_theme: legend.title 13->15, legend.text 11->13
sed -i 's/legend.title  = element_text(size = 13, face = "bold")/legend.title  = element_text(size = 15, face = "bold")/' "$F"
sed -i 's/legend.text   = element_text(size = 11)/legend.text   = element_text(size = 13)/' "$F"

# 2) pA: 加列间距, 保持左对齐
sed -i 's/legend.justification = c(0, 0.5))  # LEGEND_BOTTOM_A LEGEND_SPREAD/legend.justification = c(0, 0.5), legend.key.spacing.x = unit(10, "pt"))  # LEGEND_BOTTOM_A LEGEND_SPREAD LEGEND_FONTUP/' "$F"

# 3) pD: 改居中 + 加列间距
sed -i 's/legend.justification = c(0, 0.5))  # LEGEND_BOTTOM_DE LEGEND_SPREAD/legend.justification = "center", legend.key.spacing.x = unit(10, "pt"))  # LEGEND_BOTTOM_DE LEGEND_SPREAD LEGEND_FONTUP/' "$F"

# 4) pE: 改居中 + 加列间距 (行末无注释, 用$锚定)
sed -i 's/legend.justification = c(0, 0.5))$/legend.justification = "center", legend.key.spacing.x = unit(10, "pt"))/' "$F"

echo "=== 补丁完成, 确认 ==="
grep -n 'legend.text\|legend.title\|LEGEND_FONTUP\|key.spacing\|justification' "$F"
