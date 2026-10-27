# ==================== 肝脏手动注释脚本（调试版） ====================
# 用法: Rscript liver_manual_annotation_debug.R sample_name

args <- commandArgs(trailingOnly = TRUE)
sample_name <- args[1]

library(Seurat)
library(ggplot2)
library(dplyr)

cat("========================================\n")
cat("肝脏手动注释（调试版）:", sample_name, "\n")
cat("========================================\n\n")

# 工作目录
sample_dir <- file.path("/mnt/BioDisk/wangmingyuan/SGS10_analysis", sample_name)
setwd(sample_dir)

# 加载UMAP后的对象
if(file.exists(paste0(sample_name, "_umap.rds"))) {
  seurat_obj <- readRDS(paste0(sample_name, "_umap.rds"))
  cat("✓ 加载UMAP对象，细胞数:", ncol(seurat_obj), "\n")
} else {
  stop("未找到预处理对象")
}

# ==================== 肝脏细胞标记基因数据库 ====================
create_liver_marker_db <- function() {
  liver_markers <- list()
  
  liver_markers$Hepatocytes <- c("Alb", "Afp", "Ttr", "Hnf4a", "Hnf1a", "Foxa1", "Foxa2")
  liver_markers$Hepatocytes_Periportal <- c("Cyp2e1", "Cyp1a2", "G6pc", "Pck1", "Ass1")
  liver_markers$Hepatocytes_Perivenous <- c("Cyp3a11", "Cyp2f2", "Glul", "Oat", "Gstm1")
  liver_markers$Stellate <- c("Pdgfra", "Pdgfrb", "Des", "Rgs5", "Acta2", "Lrat", "Reln")
  liver_markers$HSC_Quiescent <- c("Rbp1", "Lrat", "Pparg")
  liver_markers$HSC_Activated <- c("Acta2", "Col1a1", "Col1a2", "Timp1", "Tgfb1")
  liver_markers$LSEC <- c("Clec4g", "Lyve1", "Stab1", "Stab2", "Fcgr2b", "Oit3")
  liver_markers$Endothelial <- c("Pecam1", "Cdh5", "Vwf", "Kdr", "Flt1")
  liver_markers$Kupffer <- c("Clec4f", "Vsig4", "Timd4", "C1qa", "C1qb", "C1qc", "Cd68", "Cd163")
  liver_markers$Macrophage_M2 <- c("Arg1", "Cd163", "Mrc1", "Il10")
  liver_markers$Macrophage_NAM <- c("Trem2", "Gpnmb", "Cd9", "Spp1")
  liver_markers$Cholangiocytes <- c("Krt7", "Krt8", "Krt18", "Krt19", "Epcam", "Sox9")
  liver_markers$T_cells <- c("Cd3d", "Cd3e", "Cd3g", "Cd4", "Cd8a")
  liver_markers$CD4_T <- c("Cd4", "Cd3e", "Il2ra")
  liver_markers$CD8_T <- c("Cd8a", "Cd8b", "Cd3e")
  liver_markers$Treg <- c("Foxp3", "Il2ra", "Ctla4")
  liver_markers$B_cells <- c("Cd19", "Ms4a1", "Cd79a", "Cd79b")
  liver_markers$NK_cells <- c("Nkg7", "Klrb1", "Klrd1", "Klrk1")
  liver_markers$Monocytes <- c("Cd14", "Ly6c2", "Ccr2", "S100a8", "S100a9")
  liver_markers$Neutrophils <- c("S100a8", "S100a9", "Cxcl2", "Cxcr2", "Mmp8")
  liver_markers$Pericytes <- c("Rgs5", "Pdgfrb", "Cspg4", "Kcnj8")
  
  return(liver_markers)
}

# 创建标记数据库
liver_markers <- create_liver_marker_db()
cat("定义了", length(liver_markers), "种细胞类型\n")

# 定义亚型和大类
subtype_names <- c(
  "Hepatocytes_Periportal", "Hepatocytes_Perivenous",
  "HSC_Quiescent", "HSC_Activated",
  "Macrophage_M2", "Macrophage_NAM",
  "CD4_T", "CD8_T", "Treg",
  "LSEC"
)

major_names <- c(
  "Hepatocytes", "Stellate", "Endothelial", "Kupffer", 
  "T_cells", "B_cells", "Monocytes", "Neutrophils", 
  "Cholangiocytes", "Pericytes", "NK_cells"
)

# 获取表达数据
cat("获取表达数据...\n")
expr_data <- GetAssayData(seurat_obj, layer = "data")

# 计算得分
cat("计算得分...\n")
all_scores <- matrix(0, nrow = ncol(seurat_obj), ncol = length(liver_markers))
colnames(all_scores) <- names(liver_markers)
rownames(all_scores) <- colnames(seurat_obj)

for(ct in names(liver_markers)) {
  markers <- liver_markers[[ct]]
  markers_present <- intersect(markers, rownames(expr_data))
  if(length(markers_present) > 0) {
    scores <- colMeans(expr_data[markers_present, , drop = FALSE])
    all_scores[, ct] <- scores
  }
}

# 提取得分矩阵
subtype_scores <- all_scores[, intersect(subtype_names, colnames(all_scores)), drop = FALSE]
major_scores <- all_scores[, intersect(major_names, colnames(all_scores)), drop = FALSE]

cat("亚型得分矩阵维度:", dim(subtype_scores), "\n")
cat("大类得分矩阵维度:", dim(major_scores), "\n")

# 计算每个细胞的最高得分
subtype_max <- apply(subtype_scores, 1, max, na.rm = TRUE)
major_max <- apply(major_scores, 1, max, na.rm = TRUE)
all_max <- apply(all_scores, 1, max, na.rm = TRUE)

# 查看前10个细胞的得分情况
cat("\n前10个细胞的得分情况:\n")
for(i in 1:10) {
  cell <- colnames(seurat_obj)[i]
  cat("\n细胞", i, ":", cell, "\n")
  cat("  亚型最高分:", round(subtype_max[i], 4), "\n")
  cat("  大类最高分:", round(major_max[i], 4), "\n")
  cat("  所有最高分:", round(all_max[i], 4), "\n")
  
  # 显示亚型得分
  sub_scores <- subtype_scores[i, ]
  valid_sub <- sub_scores[!is.na(sub_scores) & sub_scores > 0]
  if(length(valid_sub) > 0) {
    cat("  亚型得分:", paste(names(valid_sub), round(valid_sub, 3), collapse=", "), "\n")
  }
  
  # 显示大类得分
  maj_scores <- major_scores[i, ]
  valid_maj <- maj_scores[!is.na(maj_scores) & maj_scores > 0]
  if(length(valid_maj) > 0) {
    cat("  大类得分:", paste(names(valid_maj), round(valid_maj, 3), collapse=", "), "\n")
  }
}

# 计算阈值
positive_scores <- all_max[all_max > 0 & is.finite(all_max)]
subtype_threshold <- quantile(positive_scores, 0.05, na.rm = TRUE)
major_threshold <- quantile(positive_scores, 0.02, na.rm = TRUE)

cat("\n阈值:\n")
cat("  亚型阈值:", round(subtype_threshold, 4), "\n")
cat("  大类阈值:", round(major_threshold, 4), "\n")

# 检查有多少细胞能达到阈值
subtype_above <- sum(subtype_max >= subtype_threshold, na.rm = TRUE)
major_above <- sum(major_max >= major_threshold, na.rm = TRUE)
all_above <- sum(all_max >= major_threshold, na.rm = TRUE)

cat("\n达到阈值的细胞数:\n")
cat("  亚型得分 >= 阈值:", subtype_above, "/", ncol(seurat_obj), "\n")
cat("  大类得分 >= 阈值:", major_above, "/", ncol(seurat_obj), "\n")
cat("  所有得分 >= 大类阈值:", all_above, "/", ncol(seurat_obj), "\n")

# 直接分配：选择每个细胞的最高得分类型（不设阈值）
cat("\n直接分配（无阈值）:\n")
direct_assignment <- apply(all_scores, 1, function(x) {
  if(all(is.na(x)) || max(x, na.rm = TRUE) <= 0) return(NA)
  return(names(which.max(x)))
})

direct_counts <- table(direct_assignment)
cat("直接分配结果:\n")
print(sort(direct_counts, decreasing = TRUE))

