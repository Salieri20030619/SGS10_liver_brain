#!/usr/bin/env bash
# check_restyle_v1.sh
# 自检 figures/restyle_v1/ : 完整文件树 + 缺图检查 + 多余文件标记
set -u
cd /mnt/BioDisk/wangmingyuan/SGS10_analysis/figures/restyle_v1 || { echo "目录不存在"; exit 1; }

FMT="pdf tiff png"
FLAT="liver_umap_grouped_v3_9 brain_umap_grouped_v3_9 \
liver_microenv_5sample brain_microenv_5sample \
liver_microenv_stacked brain_microenv_stacked \
liver_microenv_celltype_heatmap brain_microenv_celltype_heatmap \
liver_microenv_celltype_dotplot brain_microenv_celltype_dotplot \
liver_microenv_markers_heatmap brain_microenv_markers_heatmap \
Fig2a_L1_HSC_myofib Fig2b_L4_immune_vascular \
liver_3groups_circle brain_3groups_circle \
Fig4_liver_celltype_map Fig4_brain_celltype_map"
SUB="A_top_CT_on_HE D_microenvironment_on_HE E_hex_microenv_on_HE HE_only"
SUBDIRS="Liver_3 Brain_3"

echo "########## 1) 完整文件树 ##########"
find . -type f \( -name '*.pdf' -o -name '*.tiff' -o -name '*.png' \) | sort
echo
echo "文件总数: $(find . -type f \( -name '*.pdf' -o -name '*.tiff' -o -name '*.png' \) | wc -l)"

echo
echo "########## 2) 缺图检查 (有 ✗ 即缺该格式) ##########"
exp=""
chk(){ local p="$1"; for f in $FMT; do
  if [ -f "$p.$f" ]; then exp="$exp ./$p.$f"; else echo "  ✗ 缺: $p.$f"; fi
done; }
for b in $FLAT; do chk "$b"; done
for d in $SUBDIRS; do for b in $SUB; do chk "$d/$b"; done; done
ll=$(find . -type f \( -iname '*lollipop*' -o -iname '*HSC*DEG*' \) | sort)
echo "  -- img13 lollipop 候选: ${ll:-未找到!}"
exp="$exp $ll"

echo
echo "########## 3) 不在清单内的文件 (可能多余/旧版残留, 需你判断) ##########"
find . -type f \( -name '*.pdf' -o -name '*.tiff' -o -name '*.png' \) | sort | while read -r f; do
  echo "$exp" | tr ' ' '\n' | grep -qxF "$f" || echo "  ⚠ $f"
done
echo "  (典型残留: B_*/C_* 旧面板, Fig2c*/Fig3*, liver/brain_cellchat_chord)"

echo
echo "########## 4) img23/24 (在 cellchat 目录, 不在 restyle_v1) ##########"
ls -1 ../../taskB/cellchat/pathway_plots/brain/ 2>/dev/null || echo "  目录不存在?"
