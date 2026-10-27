#!/bin/bash
# D4 step1 准备扫盘: 看脑端 D1 输出 + D3 step1 脚本细节
# 用法: bash scan_d4_prep.sh > scan_d4_prep.txt 2>&1

set -u
PR="/mnt/BioDisk/wangmingyuan/SGS10_analysis"

echo "==============================================================================="
echo "D4 step1 准备扫盘"
echo "==============================================================================="

# 1. D3 step1 脚本看 grouping 怎么写的, 沿用其框架
echo ""
echo "### 1. D3_step1_per_cell_DEG.R 主体逻辑"
echo "-------------------------------------------------------------------------------"
SCRIPT="$PR/_D3_FINAL/scripts/D3_step1_per_cell_DEG.R"
if [ -f "$SCRIPT" ]; then
    echo "总行数: $(wc -l < $SCRIPT)"
    echo ""
    echo "--- grep contrast/zone/grouping 相关 ---"
    grep -n -E "(contrast|zone|grouping|FindMarkers|wilcox|pseudobulk|null|Bla.*Con)" "$SCRIPT" | head -50
    echo ""
    echo "--- 整体 head 100 行 (看输入加载逻辑) ---"
    head -100 "$SCRIPT"
fi

# 2. D1 脑端输出
echo ""
echo "### 2. D1 脑端 niche 数据"
echo "-------------------------------------------------------------------------------"
ls -la "$PR/D1_niche_output/_D1_FINAL/" 2>&1
echo ""
echo "--- niches_by_sample_brain.rds 结构 ---"
Rscript -e 'x <- readRDS("'$PR'/D1_niche_output/_D1_FINAL/niches_by_sample_brain.rds")
cat("class:", class(x), "\n")
cat("length:", length(x), "\n")
cat("names:", names(x), "\n")
if (is.list(x)) {
    cat("\n--- 第一个元素结构 (sample 1) ---\n")
    e <- x[[1]]
    cat("class:", class(e), "\n")
    if ("data.frame" %in% class(e) || is.data.frame(e)) {
        cat("dim:", dim(e), "\n")
        cat("colnames:", colnames(e), "\n")
        print(head(e, 5))
    } else if (inherits(e, "Seurat")) {
        cat("Seurat object\n")
        print(e)
        cat("\n--- meta.data colnames ---\n")
        print(colnames(e@meta.data))
        cat("\n--- meta.data head 3 ---\n")
        print(head(e@meta.data, 3))
    } else {
        cat("(其他类型)\n")
        str(e, max.level = 2)
    }
}' 2>&1 | head -80

# 3. 确认 spot_zone_stage_brain 在哪 (D2 的脑端 zone+stage 元数据)
echo ""
echo "### 3. 脑端 spot_zone_stage / D2 输出"
echo "-------------------------------------------------------------------------------"
find "$PR/D2_zone_output" -type f 2>/dev/null | head -10
echo ""
echo "--- spot_zone_stage_brain ---"
ls -la "$PR/D2_zone_output/_D2_FINAL/"*brain* 2>&1

# 4. 看 spatial_results/rctd 的 RDS 结构 (要从这里读 raw count)
echo ""
echo "### 4. RCTD 处理后的 spatial RDS"
echo "-------------------------------------------------------------------------------"
ls -la "$PR/spatial_results/rctd/enriched/" 2>&1 | head -10
echo ""
echo "--- 一个脑端 sample RDS 的结构 ---"
RDS_SAMPLE=$(ls "$PR/spatial_results/rctd/enriched/"*brain*.rds 2>/dev/null | head -1)
if [ -n "$RDS_SAMPLE" ]; then
    echo "示例: $RDS_SAMPLE"
    Rscript -e 'x <- readRDS("'$RDS_SAMPLE'")
cat("class:", class(x), "\n")
if (inherits(x, "Seurat")) {
    print(x)
    cat("\nAssays:", names(x@assays), "\n")
    cat("\nmeta.data colnames:", colnames(x@meta.data), "\n")
    cat("\nmeta.data head 3:\n"); print(head(x@meta.data, 3))
}' 2>&1 | head -40
fi

# 5. D3 step1 脚本剩余部分 (FindMarkers 调用 + null filter 逻辑)
echo ""
echo "### 5. D3 step1 脚本: FindMarkers + null filter"
echo "-------------------------------------------------------------------------------"
if [ -f "$SCRIPT" ]; then
    echo "--- grep null / Bla vs Con 段 ---"
    grep -n -A 15 "Bla.*Con\|_NULL_\|null_ref" "$SCRIPT" | head -60
    echo ""
    echo "--- grep FindMarkers / wilcox 调用段 ---"
    grep -n -A 5 "FindMarkers\|wilcox" "$SCRIPT" | head -40
fi

echo ""
echo "==============================================================================="
echo "扫盘完成"
echo "==============================================================================="
