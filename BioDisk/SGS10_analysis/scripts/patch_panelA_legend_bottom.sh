#!/usr/bin/env bash
# patch_panelA_legend_bottom.sh —— Panel A Top CT 图例移到底部横向展开
# 用法: bash scripts/patch_panelA_legend_bottom.sh
set -e
F=/mnt/BioDisk/wangmingyuan/SGS10_analysis/scripts/he_overlay_RCTD100.R

cp -n "$F" "${F}.bak_panelA_leg" || true
if grep -q 'LEGEND_BOTTOM_A' "$F"; then
  echo "已含 LEGEND_BOTTOM_A, 跳过。重打请先: cp ${F}.bak_panelA_leg $F"
  exit 0
fi

# 1) Panel A guides: ncol 2->6, 追加 theme(legend.position="bottom")
sed -i 's/ncol = 2))/ncol = 6)) +\n    theme(legend.position = "bottom")  # LEGEND_BOTTOM_A/' "$F"

# 2) Panel A save_fig_all 高度 +3 英寸 (容纳底部图例; 肝/脑均适用)
sed -i '/A_top_CT_on_HE/{n; s/height = plot_height)/height = plot_height + 3)/;}' "$F"

echo "=== 补丁完成, 确认 ==="
grep -n 'LEGEND_BOTTOM_A\|plot_height + 3' "$F"
