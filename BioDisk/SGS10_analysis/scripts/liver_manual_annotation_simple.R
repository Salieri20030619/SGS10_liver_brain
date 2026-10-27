# ==================== 肝脏手动注释脚本（简单有效版） ====================
# 用法: Rscript liver_manual_annotation_simple.R sample_name

args <- commandArgs(trailingOnly = TRUE)
sample_name <- args[1]

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

cat("========================================\n")
cat("肝脏手动注释（简单有效版）:", sample_name, "\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================\n\n")

# 工作目录
sample_dir <- file.path("/mnt/BioDisk/wangmingyuan/SGS10_analysis", sample_name)
setwd(sample_dir)

# 加载UMAP后的对象
if(file.exists(paste0(sample_name, "_umap.rds"))) {
  seurat_obj <- readRDS(paste0(sample_name, "_umap.rds"))
  cat("✓ 加载UMAP对象，细胞数:", ncol(seurat_obj), "\n")
} else if(file.exists(paste0(sample_name, "_singler.rds"))) {
  seurat_obj <- readRDS(paste0(sample_name, "_singler.rds"))
  cat("✓ 加载SingleR对象，细胞数:", ncol(seurat_obj), "\n")
} else {
  stop("未找到预处理对象")
}

# ==================== 肝脏细胞标记基因数据库 ====================
create_liver_marker_db <- function() {
  liver_markers <- list()
  
  # 肝细胞
  liver_markers$Hepatocytes <- c("Alb", "Afp", "Ttr", "Hnf4a", "Hnf1a", "Foxa1", "Foxa2")
  
  # 星状细胞
  liver_markers$Stellate <- c("Pdgfra", "Pdgfrb", "Des", "Rgs5", "Acta2", "Lrat", "Reln")
  liver_markers$HSC_Activated <- c("Acta2", "Col1a1", "Col1a2", "Timp1", "Tgfb1")
  
  # 内皮细胞
  liver_markers$Endothelial <- c("Pecam1", "Cdh5", "Vwf", "Kdr", "Flt1")
  liver_markers$LSEC <- c("Clec4g", "Lyve1", "Stab1", "Stab2", "Fcgr2b", "Oit3")
  
  # 枯否细胞
  liver_markers$Kupffer <- c("Clec4f", "Vsig4", "Timd4", "C1qa", "C1qb", "C1qc", "Cd68", "Cd163")
  
  # 巨噬细胞相关
  liver_markers$Macrophage_M2 <- c("Arg1", "Cd163", "Mrc1", "Il10")
  liver_markers$Macrophage_NAM <- c("Trem2", "Gpnmb", "Cd9", "Spp1")
  
  # 胆管细胞
  liver_markers$Cholangiocytes <- c("Krt7", "Krt8", "Krt18", "Krt19", "Epcam", "Sox9")
  
  # T细胞
  liver_markers$T_cells <- c("Cd3d", "Cd3e", "Cd3g", "Cd4", "Cd8a")
  liver_markers$CD4_T <- c("Cd4", "Cd3e", "Il2ra")
  liver_markers$CD8_T <- c("Cd8a", "Cd8b", "Cd3e")
  liver_markers$Treg <- c("Foxp3", "Il2ra", "Ctla4")
  
  # B细胞
  liver_markers$B_cells <- c("Cd19", "Ms4a1", "Cd79a", "Cd79b")
  
  # NK细胞
  liver_markers$NK_cells <- c("Nkg7", "Klrb1", "Klrd1", "Klrk1")
  
  # 单核细胞
  liver_markers$Monocytes <- c("Cd14", "Ly6c2", "Ccr2", "S100a8", "S100a9")
  
  # 中性粒细胞
  liver_markers$Neutrophils <- c("S100a8", "S100a9", "Cxcl2", "Cxcr2", "Mmp8")
  
  # 周细胞
  liver_markers$Pericytes <- c("Rgs5", "Pdgfrb", "Cspg4", "Kcnj8")
  
  return(liver_markers)
}

# 创建标记数据库
liver_markers <- create_liver_marker_db()
cat("定义了", length(liver_markers), "种细胞类型\n")

# 定义大类（主要使用大类）
major_names <- c(
  "Hepatocytes", "Stellate", "Endothelial", "Kupffer", 
  "T_cells", "B_cells", "Monocytes", "Neutrophils", 
  "Cholangiocytes", "Pericytes", "NK_cells",
  "Macrophage_M2", "Macrophage_NAM", "LSEC"
)

# 获取表达数据
cat("获取表达数据...\n")
expr_data <- GetAssayData(seurat_obj, layer = "data")

# 计算得分
cat("计算每个细胞类型的marker得分...\n")
all_scores <- matrix(0, nrow = ncol(seurat_obj), ncol = length(liver_markers))
colnames(all_scores) <- names(liver_markers)
rownames(all_scores) <- colnames(seurat_obj)

for(ct in names(liver_markers)) {
  markers <- liver_markers[[ct]]
  markers_present <- intersect(markers, rownames(expr_data))
  if(length(markers_present) > 0) {
    scores <- colMeans(expr_data[markers_present, , drop = FALSE])
    all_scores[, ct] <- scores
    cat("  ✓", ct, ": 找到", length(markers_present), "个marker\n")
  } else {
    cat("  ⚠️", ct, ": 未找到marker\n")
  }
}

# 只保留有marker的大类
major_scores <- all_scores[, intersect(major_names, colnames(all_scores)), drop = FALSE]
cat("\n大类得分矩阵维度:", dim(major_scores), "\n")

# 计算阈值
all_max_scores <- apply(all_scores, 1, max)
positive_scores <- all_max_scores[all_max_scores > 0 & is.finite(all_max_scores)]

# 使用2%分位数作为阈值
if(length(positive_scores) > 0) {
  threshold <- quantile(positive_scores, 0.02, na.rm = TRUE)
  if(is.na(threshold) || threshold < 0.05) threshold <- 0.05
} else {
  threshold <- 0.05
}

cat("\n阈值:", round(threshold, 4), "\n")

# 分配细胞类型
cat("\n分配细胞类型...\n")
final_assignment <- rep(NA, ncol(seurat_obj))
names(final_assignment) <- colnames(seurat_obj)

for(cell in colnames(seurat_obj)) {
  cell_scores <- major_scores[cell, , drop = FALSE]
  if(length(cell_scores) > 0) {
    valid_idx <- !is.na(cell_scores) & is.finite(cell_scores) & cell_scores > 0
    if(any(valid_idx)) {
      valid_scores <- cell_scores[valid_idx]
      max_score <- max(valid_scores)
      if(is.finite(max_score) && max_score >= threshold) {
        best_type <- names(valid_scores)[which.max(valid_scores)]
        final_assignment[cell] <- best_type
      }
    }
  }
}

# 分配结果统计
assigned_count <- sum(!is.na(final_assignment))
cat("\n分配结果统计:\n")
cat("  已分配细胞数:", assigned_count, "/", ncol(seurat_obj), "\n")
cat("  未分配细胞数:", sum(is.na(final_assignment)), "\n")

# 用SingleR结果填补NA
na_idx <- is.na(final_assignment)
if(sum(na_idx) > 0 && "celltype" %in% colnames(seurat_obj@meta.data)) {
  cat("用SingleR结果填补", sum(na_idx), "个未分配细胞...\n")
  final_assignment[na_idx] <- as.character(seurat_obj$celltype[na_idx])
}

# 确保所有细胞都有注释
na_idx <- is.na(final_assignment)
if(sum(na_idx) > 0) {
  final_assignment[na_idx] <- "Other"
}

# 应用注释
seurat_obj$final_celltype <- final_assignment
Idents(seurat_obj) <- seurat_obj$final_celltype

cat("\n最终细胞类型分布:\n")
celltype_counts <- table(Idents(seurat_obj))
print(celltype_counts)

# 保存比例
prop_df <- data.frame(
  CellType = names(celltype_counts),
  Count = as.numeric(celltype_counts),
  Percentage = as.numeric(prop.table(celltype_counts) * 100)
)
write.csv(prop_df, paste0(sample_name, "_manual_proportions.csv"), row.names = FALSE)

# 绘制对比图
pdf(paste0(sample_name, "_annotation_comparison.pdf"), width = 16, height = 8)
p1 <- DimPlot(seurat_obj, reduction = "umap", group.by = "final_celltype", 
              label = TRUE, label.size = 3, repel = TRUE) +
  ggtitle("Liver Manual Annotation")
if("celltype" %in% colnames(seurat_obj@meta.data)) {
  p2 <- DimPlot(seurat_obj, reduction = "umap", group.by = "celltype", 
                label = TRUE, label.size = 3, repel = TRUE) +
    ggtitle("SingleR Annotation")
  print(p1 | p2)
} else {
  print(p1)
}
dev.off()
cat("✓ 注释对比图已保存\n")

# 保存对象
saveRDS(seurat_obj, paste0(sample_name, "_manual_annotated.rds"))
cat("✓ 手动注释结果已保存\n")

cat("\n========== 肝脏注释完成 ==========\n")
