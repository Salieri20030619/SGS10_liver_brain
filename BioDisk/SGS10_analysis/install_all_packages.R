#!/usr/bin/env Rscript

# ==================== 生物信息学R包安装脚本 ====================
# 用法: Rscript install_all_packages.R

cat("========================================\n")
cat("开始安装生物信息学R包\n")
cat("========================================\n\n")

# 设置镜像
cat("设置CRAN镜像...\n")
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")

# 安装BiocManager
cat("\n安装BiocManager...\n")
if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

# 定义要安装的Bioconductor包
bioc_packages <- c(
    "DESeq2",           # 差异表达分析
    "edgeR",            # 差异表达分析
    "limma",            # 线性模型分析
    "clusterProfiler",  # 富集分析
    "DOSE",             # 疾病本体富集
    "org.Hs.eg.db",     # 人类基因注释
    "org.Mm.eg.db",     # 小鼠基因注释
    "pathview",         # 通路可视化
    "EnhancedVolcano",  # 火山图
    "scater",           # 单细胞分析
    "scran",            # 单细胞分析
    "monocle3"          # 单细胞轨迹分析
)

# 安装Bioconductor包
cat("\n安装Bioconductor包...\n")
for(pkg in bioc_packages) {
    cat("  安装", pkg, "...\n")
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
}

# 定义要安装的CRAN包
cran_packages <- c(
    "ggplot2",      # 可视化
    "dplyr",        # 数据处理
    "tidyr",        # 数据整理
    "pheatmap",     # 热图
    "ggrepel",      # 标签防重叠
    "ggpubr",       # 统计可视化
    "plotly",       # 交互式绘图
    "cowplot",      # 组合图形
    "data.table",   # 高效数据处理
    "reshape2",     # 数据重塑
    "Seurat",       # 单细胞分析
    "patchwork",    # 图形组合
    "stringr",      # 字符串处理
    "purrr",        # 函数式编程
    "tibble",       # 现代数据框
    "forcats",      # 因子处理
    "RColorBrewer", # 颜色方案
    "viridis"       # 颜色方案
)

# 安装CRAN包
cat("\n安装CRAN包...\n")
for(pkg in cran_packages) {
    cat("  安装", pkg, "...\n")
    install.packages(pkg)
}

# ==================== 验证安装 ====================
cat("\n========================================\n")
cat("验证安装\n")
cat("========================================\n\n")

packages_to_check <- c(
    "DESeq2", "edgeR", "limma", "clusterProfiler", 
    "Seurat", "ggplot2", "dplyr", "pheatmap"
)

success <- TRUE
for(pkg in packages_to_check) {
    if(require(pkg, character.only = TRUE, quietly = TRUE)) {
        ver <- as.character(packageVersion(pkg))
        cat("✓", pkg, "版本:", ver, "\n")
    } else {
        cat("✗", pkg, "安装失败\n")
        success <- FALSE
    }
}

cat("\n========================================\n")
if(success) {
    cat("所有包安装成功！\n")
} else {
    cat("部分包安装失败，请检查错误信息\n")
}
cat("========================================\n")
