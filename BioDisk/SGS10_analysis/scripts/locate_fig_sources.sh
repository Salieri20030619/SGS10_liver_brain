#!/bin/bash
# locate_fig_sources.sh
# 定位所有 27 张论文图的 (1) 绘图脚本 (2) 数据源
# 用法: bash locate_fig_sources.sh > /tmp/fig_sources.txt 2>&1

ROOT=/mnt/BioDisk/wangmingyuan
cd $ROOT

# 27 个图的搜索关键词 (取 PNG 文件名核心部分)
KEYWORDS=(
  "Liver_UMAP"
  "Brain_UMAP"
  "Liver_microenv_heatmap"
  "Brain_microenv_heatmap"
  "Liver_3_spatial_overlay"
  "Brain_3_spatial_overlay"
  "Liver_markers_heatmap"
  "Liver_3_hex_bin"
  "HSC_DEG_lollipop"
  "HSC_volcano"
  "HSC_spatial_5x5"
  "HSC_myeloid_dumbbell"
  "HSC_myeloid_hub_dumbbell"
  "silent_inflammation"
  "silent_inflammation_dotplot"
  "X2_brain_spatial"
  "TF_dual_mode"
  "Excit_L45_TF_lollipop"
  "Brain_3_hex_bin"
  "X1_TGFb_spatial"
  "X1_TGFb_bar"
  "thyroid_cascade"
  "efferocytosis_3layers"
  "Brain_CellChat_3groups"
  "Liver_CellChat_3groups"
  "species_separation"
  "DAM_cleaning"
  "Brain_markers_heatmap"
)

echo "============================================"
echo "STAGE 1: 查找每张图对应的 R/Py 绘图脚本"
echo "============================================"
for kw in "${KEYWORDS[@]}"; do
  echo ""
  echo "### [$kw] ###"
  # 找包含该关键词的 .R 或 .py 脚本 (排除 png/pdf/checkpoint)
  find . -type f \( -name "*.R" -o -name "*.py" -o -name "*.Rmd" \) \
    2>/dev/null | xargs grep -l -i "$kw" 2>/dev/null | head -5
done

echo ""
echo ""
echo "============================================"
echo "STAGE 2: 列出所有看起来是绘图脚本的文件"
echo "  (按修改时间倒序, 前 50)"
echo "============================================"
find . -type f \( -name "*.R" -o -name "*.py" \) \
  -newer /etc/hostname 2>/dev/null \
  | xargs ls -lt 2>/dev/null | head -50

echo ""
echo ""
echo "============================================"
echo "STAGE 3: 关键数据目录 (CSV / RDS 索引)"
echo "============================================"
echo ""
echo "--- 所有 RDS 文件 (前 40, 按大小倒序) ---"
find . -name "*.rds" 2>/dev/null | xargs ls -lS 2>/dev/null | head -40

echo ""
echo "--- 关键 CSV 文件 (含 DEG/lollipop/dumbbell/cascade/efferocyt/inflammation/cellchat 关键词) ---"
find . -name "*.csv" 2>/dev/null \
  | grep -iE "DEG|lollipop|dumbbell|cascade|efferocyt|inflammation|cellchat|TF_|microenv|silent|paired_corr" \
  | head -30

echo ""
echo "--- taskB 目录结构 (深度 3) ---"
find taskB -maxdepth 3 -type d 2>/dev/null | head -40

echo ""
echo "--- SGS10_analysis 目录结构 (深度 3) ---"
find SGS10_analysis -maxdepth 3 -type d 2>/dev/null | head -40

echo ""
echo "============================================"
echo "STAGE 4: 'scripts' 'paper_figs' 这类目录"
echo "============================================"
find . -maxdepth 4 -type d \( -name "*script*" -o -name "*paper*" -o -name "*figs*" -o -name "*plot*" \) 2>/dev/null | head -20

echo ""
echo "DONE."
