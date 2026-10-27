options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")

# 安装必需包
install.packages("BiocManager")
BiocManager::install(c("DESeq2", "edgeR", "limma", "EnhancedVolcano"), update = FALSE)
install.packages(c("ggplot2", "dplyr", "pheatmap"))

cat("完成！\n")
