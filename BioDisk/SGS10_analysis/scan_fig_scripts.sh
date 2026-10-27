#!/usr/bin/env bash
# ============================================================
# scan_fig_scripts.sh
# 用途：在工作目录下按"图标题特征字符串"定位生成各张图的脚本
# 用法：bash scan_fig_scripts.sh [搜索根目录]
#   默认搜索根目录 = /mnt/BioDisk/wangmingyuan
# 输出：每个锚点命中的 文件:行号:内容
# 说明：覆盖 .R .r .py .ipynb .sh .Rmd 等常见脚本；ipynb 按纯文本 grep
# ============================================================
set -u
ROOT="${1:-/mnt/BioDisk/wangmingyuan}"
echo "搜索根目录: $ROOT"
echo "============================================================"

# 锚点：图序号 | 描述 | 用于 grep 的特征字符串（固定字符串模式 -F）
declare -a ANCHORS=(
"img1|肝UMAP图谱|Cell Atlas"
"img1_2|肝UMAP细胞数|136711 cells"
"img2|脑UMAP细胞数|258990 cells"
"img5_6|空间TopCT|HE + Top CT"
"img5_6b|空间TopCT|argmax"
"img7_8_18|空间微环境|Microenvironment"
"img9|微环境占比柱|Proportion"
"img10_19|微环境x细胞热图|Microenv x Cell Type"
"img10_19b|微环境x细胞热图|mean weight"
"img11_20|Top5气泡图|Top 5 cell types per microenv"
"img13|HSC棒棒糖DEG|myofibroblast transition"
"img13b|HSC棒棒糖DEG|Top 10 DEG"
"img14|空间基因网格L1|HSC myofibroblast activation"
"img16|空间基因网格L4|immune-vascular hub"
"img17|脑六边形微环境|HE + Hex"
"img15_22|CellChat圈图|Number of interactions"
"img15_22b|CellChat圈图|Interaction weight"
"cellchat|CellChat通用|netVisual"
"doheatmap|Seurat热图|DoHeatmap"
"zone|分区注释|host_dom"
)

# 文件类型过滤
INCLUDES=(--include=*.R --include=*.r --include=*.py --include=*.ipynb \
          --include=*.sh --include=*.Rmd --include=*.qmd --include=*.txt)

for entry in "${ANCHORS[@]}"; do
  IFS='|' read -r tag desc pat <<< "$entry"
  echo ""
  echo ">>> [$tag] $desc  —— 锚点: \"$pat\""
  grep -rInF "${INCLUDES[@]}" -- "$pat" "$ROOT" 2>/dev/null \
    | grep -v '/\.git/' \
    | head -20
done

echo ""
echo "============================================================"
echo "扫描完成。把上面命中的脚本路径整理出来发我，我逐个改。"
echo "若某锚点 0 命中，可能字符串被拼接/变量化，告诉我该图序号我换策略。"
