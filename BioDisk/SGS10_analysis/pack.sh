#!/bin/bash
# ============================================================================
# SGS10 - 打包 D4_GO_output_v2 完整结果
# 用法: bash pack_D4_GO.sh
# 产出: D4_GO_output_v2_<date>.tar.gz (放在项目根)
# ============================================================================

set -eu
PROJECT_ROOT="/mnt/BioDisk/wangmingyuan/SGS10_analysis"
SRC_DIR="$PROJECT_ROOT/D4_GO_output_v2"
DATE_TAG=$(date +%Y%m%d_%H%M%S)
PKG_NAME="D4_GO_output_v2_${DATE_TAG}.tar.gz"
PKG_PATH="$PROJECT_ROOT/$PKG_NAME"

cd "$PROJECT_ROOT"

# ---------- 1. 检查源目录 ----------
echo "==============================================================================="
echo "SGS10 D4_GO_output_v2 打包脚本"
echo "时间: $(date)"
echo "==============================================================================="
if [ ! -d "$SRC_DIR" ]; then
    echo "[ERR] $SRC_DIR 不存在"
    exit 1
fi
echo ""
echo "源目录: $SRC_DIR"
du -sh "$SRC_DIR" 2>&1
echo ""
echo "源目录结构 (深度 2):"
find "$SRC_DIR" -maxdepth 2 -type d | sort
echo ""
echo "文件统计:"
echo "  enrichment/ csv 数: $(find $SRC_DIR/enrichment -name '*.csv' 2>/dev/null | wc -l)"
echo "  figures/    pdf 数: $(find $SRC_DIR/figures -name '*.pdf' 2>/dev/null | wc -l)"
echo "  cross/      文件数: $(find $SRC_DIR/cross -type f 2>/dev/null | wc -l)"
echo "  logs/       log 数: $(find $SRC_DIR/logs -type f 2>/dev/null | wc -l)"
echo "  task_summary.csv: $([ -f $SRC_DIR/task_summary.csv ] && echo OK || echo MISS)"

# ---------- 2. 同时收一份 D3 必要的输入文件作为 reference (方便溯源) ----------
echo ""
echo "顺手附带 D3 关键输入 (filtered_DEG + signatures, 用于溯源):"
TMP_REF=$(mktemp -d)
mkdir -p "$TMP_REF/D3_reference/filtered_DEG"
mkdir -p "$TMP_REF/D3_reference/signatures"
cp -p "$PROJECT_ROOT/_D3_FINAL/signature_output/filtered_DEG/"*.csv "$TMP_REF/D3_reference/filtered_DEG/" 2>/dev/null || true
cp -p "$PROJECT_ROOT/_D3_FINAL/signature_output/signatures/"*.txt "$TMP_REF/D3_reference/signatures/" 2>/dev/null || true
echo "  filtered_DEG: $(ls $TMP_REF/D3_reference/filtered_DEG/ 2>/dev/null | wc -l) csv"
echo "  signatures:   $(ls $TMP_REF/D3_reference/signatures/ 2>/dev/null | wc -l) txt"

# ---------- 3. 写一份 README ----------
cat > "$TMP_REF/D3_reference/README.txt" <<'EOF'
SGS10 D3 富集分析输入文件 (作 D4_GO_output_v2 溯源用)
==================================================

filtered_DEG/   D3 step2 产出, null-filtered robust DEG (UP/DOWN per contrast)
                列: p_val, avg_log2FC, pct.1, pct.2, p_val_adj, gene, sample_consistency
                
                A1_UP_filtered.csv     2 gene  (Celf2, Gm42418) - 富集已跳过
                A1_DOWN_filtered.csv   977 gene
                A2_UP_filtered.csv     2486 gene  ⭐ invaded host signature
                A2_DOWN_filtered.csv   21 gene
                A3_UP_filtered.csv     2249 gene  ⭐ distant host signature
                A3_DOWN_filtered.csv   14 gene
                B3_UP_filtered.csv     2124 gene  (sample-driven, paper 不用)
                B3_DOWN_filtered.csv   168 gene
                C4_UP_filtered.csv     141 gene  ⭐⭐⭐ 脑端核心 finding
                C4_DOWN_filtered.csv   321 gene

signatures/     D3 step2 产出, top 50 UP gene 简版 (UCell scoring 用)
                每文件一行一个 gene symbol

EOF

# ---------- 4. 打包 ----------
echo ""
echo "打包中 ..."
tar -czf "$PKG_PATH" \
    -C "$PROJECT_ROOT" "D4_GO_output_v2" \
    -C "$TMP_REF" "D3_reference"

# 清理临时目录
rm -rf "$TMP_REF"

# ---------- 5. 完成 ----------
PKG_SIZE=$(du -h "$PKG_PATH" | cut -f1)
echo ""
echo "==============================================================================="
echo "打包完成"
echo "==============================================================================="
echo "包路径:    $PKG_PATH"
echo "大小:      $PKG_SIZE"
echo "MD5:       $(md5sum $PKG_PATH | cut -d' ' -f1)"
echo ""
echo "下载到本地 D:/SGS10 (Windows PowerShell 或 Git Bash):"
echo "  scp -P 8022 wangmingyuan@g2.cnsknowall.com:$PKG_PATH D:/SGS10/"
echo ""
echo "(或 cmd 下用反斜杠路径)"
echo "  scp -P 8022 wangmingyuan@g2.cnsknowall.com:$PKG_PATH D:\\SGS10\\"
echo ""
echo "解压 (Windows 用 7-Zip / WinRAR / tar 命令):"
echo "  tar -xzf $PKG_NAME -C D:/SGS10/"
