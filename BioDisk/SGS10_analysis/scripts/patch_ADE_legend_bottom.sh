#!/usr/bin/env bash
# patch_ADE_legend_bottom.sh —— A/D/E 三张图图例全部移到底部横向展开
# 兼容"A 已单独打过补丁"的情况
# 用法: bash scripts/patch_ADE_legend_bottom.sh
set -e
F=/mnt/BioDisk/wangmingyuan/SGS10_analysis/scripts/he_overlay_RCTD100.R

cp -n "$F" "${F}.bak_ADE_leg" || true

# ---------- Panel A (若已含 LEGEND_BOTTOM_A 则跳过) ----------
if ! grep -q 'LEGEND_BOTTOM_A' "$F"; then
  sed -i 's/ncol = 2))/ncol = 6)) +\n    theme(legend.position = "bottom")  # LEGEND_BOTTOM_A/' "$F"
  sed -i '/A_top_CT_on_HE/{n; s/height = plot_height)/height = plot_height + 3)/;}' "$F"
  echo "  [A] 已打补丁"
else
  echo "  [A] 已含 LEGEND_BOTTOM_A, 跳过"
fi

# ---------- Panel D + E ----------
if ! grep -q 'LEGEND_BOTTOM_DE' "$F"; then
  # D: size=3,alpha=1 是其独有特征 (E只有alpha=1)
  sed -i 's/size = 3, alpha = 1), ncol = 1))/size = 3, alpha = 1), ncol = 3)) +\n    theme(legend.position = "bottom")  # LEGEND_BOTTOM_DE/' "$F"
  # E: list(alpha=1) 是其独有特征
  sed -i 's/list(alpha = 1), ncol = 1))/list(alpha = 1), ncol = 3)) +\n    theme(legend.position = "bottom")/' "$F"
  # D 保存高度 +2
  sed -i '/D_microenvironment_on_HE/{n; s/height = plot_height)/height = plot_height + 2)/;}' "$F"
  # E 保存高度 +2
  sed -i '/E_hex_microenv_on_HE/{n; s/height = plot_height)/height = plot_height + 2)/;}' "$F"
  echo "  [D+E] 已打补丁"
else
  echo "  [D+E] 已含 LEGEND_BOTTOM_DE, 跳过"
fi

echo "=== 确认插入 ==="
grep -n 'LEGEND_BOTTOM\|plot_height + [0-9]' "$F"
