#!/usr/bin/env bash
# patch_5sample_rotate180.sh —— 给 microenv_5sample_montage.R 加"脑样本转180度"
# 用法: bash scripts/patch_5sample_rotate180.sh
set -e
F=/mnt/BioDisk/wangmingyuan/SGS10_analysis/scripts/microenv_5sample_montage.R

cp -n "$F" "${F}.bak_rotate180" || true
if grep -q 'rotate180' "$F"; then
  echo "已含 rotate180, 跳过避免重复。重打请先: cp ${F}.bak_rotate180 $F"
  exit 0
fi

# 1) build_panel 函数签名 加 rotate180 参数 (默认FALSE)
sed -i 's/me_labels, me_title) {/me_labels, me_title, rotate180 = FALSE) {/' "$F"

# 2) merge 之后 插入坐标镜像块 (仅 rotate180=TRUE 时执行)
sed -i '/by = "barcode")/a\
  if (rotate180) {   # [client] 脑样本坐标绕各样本中心镜像 -> 整体转180度\
    df$x <- max(df$x) + min(df$x) - df$x\
    df$y <- max(df$y) + min(df$y) - df$y\
  }' "$F"

# 3) lapply 调用 传 rotate180 flag (brain=TRUE, liver=FALSE)
sed -i 's/me_title = me_title)/me_title = me_title, rotate180 = (organ == "brain"))/' "$F"

echo "=== 打补丁完成, 确认插入 ==="
grep -n 'rotate180' "$F"
