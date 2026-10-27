#!/usr/bin/env bash
# patch_5sample_fix_bla.sh —— 纠正: Bla_brain_1 坐标本已正确, 不应转180度
# 在已打 patch_5sample_rotate180.sh 的文件上继续改
# 用法: bash scripts/patch_5sample_fix_bla.sh
set -e
F=/mnt/BioDisk/wangmingyuan/SGS10_analysis/scripts/microenv_5sample_montage.R

cp -n "$F" "${F}.bak_fix_bla" || true
if grep -q 'ROTATE_BRAIN_SAMPLES' "$F"; then
  echo "已含 ROTATE_BRAIN_SAMPLES, 跳过。重打请先: cp ${F}.bak_fix_bla $F"
  exit 0
fi

# 1) 在 load_coords 函数定义前插入 ROTATE_BRAIN_SAMPLES 向量 (排除 Bla_brain_1)
sed -i '/# ---- 取某样本坐标/i\
# Bla_brain_1 原始坐标朝向已正确, 其余4个脑样本需转180度\
ROTATE_BRAIN_SAMPLES <- c("Con_brain_1", "Brain_3", "Brain_8", "nao56")\
' "$F"

# 2) build_panel 签名: 去掉 rotate180 参数 (改为内部按样本名判断)
sed -i 's/, rotate180 = FALSE) {/) {/' "$F"

# 3) if(rotate180) -> if(sample %in% ROTATE_BRAIN_SAMPLES)
sed -i 's/if (rotate180) {   # \[client\]/if (sample %in% ROTATE_BRAIN_SAMPLES) {  # [client]/' "$F"

# 4) lapply 调用: 去掉 rotate180 参数 (内部自判断, 不再需要外部传入)
sed -i 's/, rotate180 = (organ == "brain")//' "$F"

echo "=== 纠正补丁完成, 确认 ==="
grep -n 'ROTATE_BRAIN_SAMPLES\|rotate180\|Bla_brain' "$F"
