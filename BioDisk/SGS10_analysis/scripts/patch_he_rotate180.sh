#!/usr/bin/env bash
# patch_he_rotate180.sh —— 给 he_overlay_RCTD100.R 加"脑HE系列整体转180度"
# 用法: bash patch_he_rotate180.sh
set -e
F=/mnt/BioDisk/wangmingyuan/SGS10_analysis/scripts/he_overlay_RCTD100.R

cp -n "$F" "${F}.bak_rotate180" || true
if grep -q 'ROTATE180' "$F"; then
  echo "已含 ROTATE180, 跳过避免重复。重打请先: cp ${F}.bak_rotate180 $F"
  exit 0
fi

# 1) size 日志行后 -> 定义 ROTATE180 (仅脑)
sed -i '/log_msg(sprintf("  size: %d x %d", size_x, size_y))/a\
  ROTATE180 <- (cfg$tissue == "brain")   # [client] 脑 HE 系列整体转 180度\
  if (ROTATE180) log_msg("  ROTATE180: brain HE + points mirrored")' "$F"

# 2) 读 HE 后 -> 栅格转180(行/列反序)
sed -i '/he_img <- png::readPNG(cfg\$he)/a\
  if (ROTATE180) {\
    .d <- dim(he_img)\
    he_img <- if (length(.d) == 3) he_img[.d[1]:1, .d[2]:1, , drop = FALSE] else he_img[.d[1]:1, .d[2]:1]\
  }' "$F"

# 3) RCTD 日志行后 -> 点坐标绕中心镜像
sed -i '/log_msg(sprintf("  RCTD: %d spots, %d CT", nrow(coords), ncol(weights)))/a\
  if (ROTATE180) { coords$x <- size_x - coords$x; coords$y <- size_y - coords$y }  # [client] 点坐标镜像' "$F"

echo "=== 打补丁完成, 确认插入 ==="
grep -n 'ROTATE180\|.d <- dim(he_img)\|coords$x <- size_x - coords$x' "$F"
