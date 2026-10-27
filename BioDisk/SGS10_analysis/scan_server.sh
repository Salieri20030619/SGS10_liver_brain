#!/bin/bash
# ============================================================================
# SGS10 服务器扫盘脚本 (READ-ONLY)
# 目的: 定位 D3 输出结构, 为 D3 GO/KEGG/Reactome 脚本编写做准备
# 用法: bash scan_server.sh > scan_report.txt 2>&1
# 服务器: wangmingyuan@g2.cnsknowall.com:8022
# ============================================================================

set -u  # undefined var 报错（但不 set -e, 因为有些命令期望失败）

echo "==============================================================================="
echo "SGS10 服务器扫盘报告"
echo "时间: $(date)"
echo "用户: $(whoami)"
echo "主机: $(hostname)"
echo "工作目录: $(pwd)"
echo "==============================================================================="

# ---------- 1. 定位项目根目录 ----------
echo ""
echo "### 1. 项目根目录候选"
echo "-------------------------------------------------------------------------------"
CANDIDATES=(
    "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
    "$HOME/SGS10_analysis"
    "$HOME/projects/SGS10"
    "/data/SGS10_analysis"
    "/data/wangmingyuan/SGS10_analysis"
)

PROJECT_ROOT=""
for d in "${CANDIDATES[@]}"; do
    if [ -d "$d" ]; then
        echo "[FOUND] $d"
        if [ -z "$PROJECT_ROOT" ]; then
            PROJECT_ROOT="$d"
        fi
    else
        echo "[MISS]  $d"
    fi
done

if [ -z "$PROJECT_ROOT" ]; then
    echo ""
    echo "未找到候选路径, 尝试在 \$HOME 和 /mnt 下搜索 'SGS10' 关键字 (深度限 5)..."
    find "$HOME" -maxdepth 5 -type d -name '*SGS10*' 2>/dev/null | head -20
    find /mnt -maxdepth 5 -type d -name '*SGS10*' 2>/dev/null | head -20
    echo ""
    echo "请手动设置 PROJECT_ROOT 后重跑。脚本退出。"
    exit 1
fi

echo ""
echo "采用 PROJECT_ROOT = $PROJECT_ROOT"

# ---------- 2. 项目顶层结构 ----------
echo ""
echo "### 2. 项目顶层目录结构 (深度 2)"
echo "-------------------------------------------------------------------------------"
ls -la "$PROJECT_ROOT" 2>/dev/null
echo ""
echo "--- 一级子目录 ---"
find "$PROJECT_ROOT" -maxdepth 2 -type d 2>/dev/null | sort

# ---------- 3. 重点定位 D3 输出 ----------
echo ""
echo "### 3. D3 输出目录定位"
echo "-------------------------------------------------------------------------------"
D3_CANDIDATES=$(find "$PROJECT_ROOT" -maxdepth 3 -type d \( -name '*D3*' -o -name '*signature*' -o -name '*filtered_DEG*' \) 2>/dev/null)
if [ -z "$D3_CANDIDATES" ]; then
    echo "未发现 D3 / signature / filtered_DEG 目录, 在更深一层搜索..."
    D3_CANDIDATES=$(find "$PROJECT_ROOT" -maxdepth 5 -type d \( -name '*D3*' -o -name '*signature*' -o -name '*filtered_DEG*' \) 2>/dev/null | head -30)
fi
echo "$D3_CANDIDATES"

# ---------- 4. filtered_DEG 文件 ----------
echo ""
echo "### 4. filtered_DEG 文件枚举"
echo "-------------------------------------------------------------------------------"
FILTERED_FILES=$(find "$PROJECT_ROOT" -maxdepth 6 -type f \( -name '*robust_UP*' -o -name '*robust_DOWN*' -o -name '*filtered_UP*' -o -name '*filtered_DOWN*' \) 2>/dev/null)
echo "$FILTERED_FILES" | head -40
echo ""
N_FILTERED=$(echo "$FILTERED_FILES" | grep -c -v '^$')
echo "共找到 $N_FILTERED 个 filtered DEG 文件"

# ---------- 5. 抽 1 个 UP 文件查看列结构 ----------
echo ""
echo "### 5. 文件列结构抽样"
echo "-------------------------------------------------------------------------------"
SAMPLE_UP=$(echo "$FILTERED_FILES" | grep -E '(A2|C4).*UP' | head -1)
SAMPLE_UP=${SAMPLE_UP:-$(echo "$FILTERED_FILES" | grep -i UP | head -1)}

if [ -n "$SAMPLE_UP" ] && [ -f "$SAMPLE_UP" ]; then
    echo "示例文件: $SAMPLE_UP"
    echo "文件大小: $(du -h "$SAMPLE_UP" | cut -f1)"
    echo "总行数: $(wc -l < "$SAMPLE_UP")"
    echo ""
    echo "--- 前 3 行 ---"
    head -3 "$SAMPLE_UP"
    echo ""
    echo "--- 列分隔符检测 ---"
    head -1 "$SAMPLE_UP" | awk '{
        n_comma = gsub(",", ",", $0); 
        n_tab = gsub("\t", "\t", $0); 
        print "  逗号数: " n_comma "  制表符数: " n_tab
    }'
    echo ""
    echo "--- 列名 (按逗号拆分) ---"
    head -1 "$SAMPLE_UP" | tr ',' '\n' | nl
else
    echo "未找到合适的示例文件"
fi

# ---------- 6. signature gene list 文件 (top50) ----------
echo ""
echo "### 6. Signature top50 gene list 文件"
echo "-------------------------------------------------------------------------------"
SIG_FILES=$(find "$PROJECT_ROOT" -maxdepth 6 -type f \( -name '*top50*' -o -name '*signature*.txt' -o -name '*signature*.csv' \) 2>/dev/null | head -20)
echo "$SIG_FILES"

# ---------- 7. wilcoxon 原始 DEG 输出 (备用) ----------
echo ""
echo "### 7. Wilcoxon / pseudobulk 原始 DEG 文件 (备份用)"
echo "-------------------------------------------------------------------------------"
RAW_DEG=$(find "$PROJECT_ROOT" -maxdepth 6 -type d \( -name 'wilcoxon' -o -name 'pseudobulk*' -o -name 'combined' \) 2>/dev/null)
echo "$RAW_DEG"
echo ""
echo "--- 原始 DEG 目录内容样本 ---"
for d in $RAW_DEG; do
    echo "[$d]"
    ls -la "$d" 2>/dev/null | head -8
    echo ""
done

# ---------- 8. R 环境检查 ----------
echo ""
echo "### 8. R 环境与关键包检查"
echo "-------------------------------------------------------------------------------"
which R
R --version 2>/dev/null | head -1

if which R >/dev/null 2>&1; then
    R --vanilla --quiet -e '
    cat("\n--- R libPaths ---\n")
    print(.libPaths())
    cat("\n--- 关键包检查 ---\n")
    pkgs <- c("clusterProfiler", "ReactomePA", "org.Mm.eg.db", "org.Hs.eg.db",
              "AnnotationDbi", "enrichplot", "DOSE", "msigdbr", "ggplot2",
              "dplyr", "tidyr", "readr", "patchwork")
    for (p in pkgs) {
        installed <- requireNamespace(p, quietly = TRUE)
        if (installed) {
            ver <- as.character(packageVersion(p))
            cat(sprintf("  [OK]  %-22s  v%s\n", p, ver))
        } else {
            cat(sprintf("  [MISS] %-22s\n", p))
        }
    }
    cat("\n--- BiocManager 版本 ---\n")
    if (requireNamespace("BiocManager", quietly = TRUE)) {
        cat("  BiocManager v", as.character(packageVersion("BiocManager")), "\n", sep="")
        cat("  Bioconductor v", as.character(BiocManager::version()), "\n", sep="")
    } else {
        cat("  [MISS] BiocManager\n")
    }
    ' 2>&1
fi

# ---------- 9. CPU / 内存 ----------
echo ""
echo "### 9. 资源信息"
echo "-------------------------------------------------------------------------------"
echo "CPU 核数: $(nproc 2>/dev/null || echo unknown)"
free -h 2>/dev/null | head -3 || echo "free 命令不可用"
echo ""
echo "项目盘空间:"
df -h "$PROJECT_ROOT" 2>/dev/null | head -2

# ---------- 10. 写出推荐路径 ----------
echo ""
echo "==============================================================================="
echo "扫盘完成"
echo "==============================================================================="
echo ""
echo "推荐 D4 enrichment 输出路径:"
echo "  $PROJECT_ROOT/D4_enrichment_output/"
echo ""
echo "请把本扫盘报告 (scan_report.txt) 反馈给 Claude, 用于编写 GO/KEGG/Reactome 脚本。"
