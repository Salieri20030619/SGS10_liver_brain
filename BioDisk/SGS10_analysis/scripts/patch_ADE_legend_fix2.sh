#!/usr/bin/env bash
# patch_ADE_legend_fix2.sh —— A溢出修复+全局加粗+D/E改1行居中
# 前提: 已运行 patch_ADE_legend_fontup.sh
# 用法: bash scripts/patch_ADE_legend_fix2.sh
set -e
F=/mnt/BioDisk/wangmingyuan/SGS10_analysis/scripts/he_overlay_RCTD100.R

cp -n "$F" "${F}.bak_fix2" || true
if grep -q 'LEGEND_FIX2' "$F"; then
  echo "已含 LEGEND_FIX2, 跳过。重打请先: cp ${F}.bak_fix2 $F"
  exit 0
fi

# 1) 全局: legend.text 加粗
sed -i 's/legend.text   = element_text(size = 13)/legend.text   = element_text(size = 13, face = "bold")/' "$F"

# 2) pA: 间距 10pt->5pt (防溢出)
sed -i 's/unit(10, "pt"))  # LEGEND_BOTTOM_A/unit(5, "pt"))  # LEGEND_BOTTOM_A LEGEND_FIX2/' "$F"

# 3) pD: ncol 6->7 (7条目1行, 消除孤行) + 间距 10->5pt
sed -i 's/size = 3, alpha = 1), ncol = 6, byrow/size = 3, alpha = 1), ncol = 7, byrow/' "$F"
sed -i 's/unit(10, "pt"))  # LEGEND_BOTTOM_DE/unit(5, "pt"))  # LEGEND_BOTTOM_DE/' "$F"

# 4) pE: ncol 6->7 + 间距 10->5pt (行末无注释, 用$锚定)
sed -i 's/list(alpha = 1), ncol = 6, byrow/list(alpha = 1), ncol = 7, byrow/' "$F"
sed -i 's/unit(10, "pt"))$/unit(5, "pt"))/' "$F"

echo "=== 补丁完成, 确认 ==="
grep -n 'face.*bold\|ncol = 7\|ncol = 8\|unit(5\|LEGEND_FIX2' "$F"
