#!/bin/bash
# 30 秒微扫盘: 看 filtered.csv 的列结构 (第二轮没打印)
D3_DIR="/mnt/BioDisk/wangmingyuan/SGS10_analysis/_D3_FINAL/signature_output/filtered_DEG"

echo "===== 看 filtered.csv 列名格式 ====="
for f in A2_UP A2_DOWN C4_UP C4_DOWN A1_UP B3_UP; do
    fp="$D3_DIR/${f}_filtered.csv"
    if [ -f "$fp" ]; then
        echo ""
        echo "=== ${f}_filtered.csv (size: $(du -h $fp | cut -f1), rows: $(wc -l < $fp)) ==="
        head -3 "$fp"
        echo "..."
        echo "列数: $(head -1 $fp | awk -F',' '{print NF}')"
    fi
done

echo ""
echo "===== Null reference 表列名 ====="
NULL_LIV="/mnt/BioDisk/wangmingyuan/SGS10_analysis/_D3_FINAL/DEG_output/wilcoxon/_NULL_Bla_vs_Con_liver.csv"
if [ -f "$NULL_LIV" ]; then
    echo "rows: $(wc -l < $NULL_LIV)"
    head -3 "$NULL_LIV"
fi

echo ""
echo "===== 旧 D4_GO_output 的 enrichment csv 列名（看是否有 ID 转换问题）====="
OLD_GO="/mnt/BioDisk/wangmingyuan/SGS10_analysis/D4_GO_output/enrichment/A2_UP_GO.csv"
if [ -f "$OLD_GO" ]; then
    head -2 "$OLD_GO"
    echo "..."
fi
