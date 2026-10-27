# ==================== 肝脏手动注释脚本（调试版2） ====================
args <- commandArgs(trailingOnly = TRUE)
sample_name <- args[1]

library(Seurat)

sample_dir <- file.path("/mnt/BioDisk/wangmingyuan/SGS10_analysis", sample_name)
setwd(sample_dir)

if(file.exists(paste0(sample_name, "_umap.rds"))) {
  seurat_obj <- readRDS(paste0(sample_name, "_umap.rds"))
  cat("✓ 加载UMAP对象，细胞数:", ncol(seurat_obj), "\n")
} else {
  stop("未找到预处理对象")
}

# 创建标记数据库
create_liver_marker_db <- function() {
  liver_markers <- list()
  liver_markers$Hepatocytes <- c("Alb", "Afp", "Ttr", "Hnf4a", "Hnf1a", "Foxa1", "Foxa2")
  liver_markers$Stellate <- c("Pdgfra", "Pdgfrb", "Des", "Rgs5", "Acta2", "Lrat", "Reln")
  liver_markers$HSC_Activated <- c("Acta2", "Col1a1", "Col1a2", "Timp1", "Tgfb1")
  liver_markers$Endothelial <- c("Pecam1", "Cdh5", "Vwf", "Kdr", "Flt1")
  liver_markers$LSEC <- c("Clec4g", "Lyve1", "Stab1", "Stab2", "Fcgr2b", "Oit3")
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

liver_markers <- create_liver_marker_db()
cat("定义了", length(liver_markers), "种细胞类型\n")

major_names <- c(
  "Hepatocytes", "Stellate", "Endothelial", "Kupffer", 
  "T_cells", "B_cells", "Monocytes", "Neutrophils", 
  "Cholangiocytes", "Pericytes", "NK_cells",
  "Macrophage_M2", "Macrophage_NAM", "LSEC", "HSC_Activated"
)

expr_data <- GetAssayData(seurat_obj, layer = "data")

# 计算得分
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

# 只保留有marker的大类
major_scores <- all_scores[, intersect(major_names, colnames(all_scores)), drop = FALSE]
cat("\n大类得分矩阵维度:", dim(major_scores), "\n")

# 计算每个细胞的最高得分
all_max_scores <- apply(all_scores, 1, max, na.rm = TRUE)
positive_scores <- all_max_scores[all_max_scores > 0 & is.finite(all_max_scores)]

if(length(positive_scores) > 0) {
  threshold <- quantile(positive_scores, 0.05, na.rm = TRUE)
  if(is.na(threshold) || threshold < 0.01) threshold <- 0.01
} else {
  threshold <- 0.01
}

cat("\n阈值:", round(threshold, 4), "\n")

# 测试前10个细胞
cat("\n=== 测试前10个细胞 ===\n")
for(i in 1:10) {
  cell <- colnames(seurat_obj)[i]
  cell_scores <- all_scores[cell, , drop = FALSE]
  
  # 找出所有正得分
  valid_idx <- !is.na(cell_scores) & is.finite(cell_scores) & cell_scores > 0
  cat("\n细胞", i, ":", cell, "\n")
  cat("  正得分类型数:", sum(valid_idx), "\n")
  
  if(any(valid_idx)) {
    valid_scores <- cell_scores[valid_idx]
    max_score <- max(valid_scores)
    best_type <- names(valid_scores)[which.max(valid_scores)]
    cat("  最高分:", round(max_score, 4), "类型:", best_type, "\n")
    cat("  是否达到阈值:", max_score >= threshold, "\n")
  } else {
    cat("  没有正得分\n")
  }
}

# 统计能分配到细胞的数量
cat("\n=== 统计分配情况 ===\n")
assigned_count <- 0
for(cell in colnames(seurat_obj)[1:1000]) {
  cell_scores <- all_scores[cell, , drop = FALSE]
  valid_idx <- !is.na(cell_scores) & is.finite(cell_scores) & cell_scores > 0
  if(any(valid_idx)) {
    valid_scores <- cell_scores[valid_idx]
    max_score <- max(valid_scores)
    if(is.finite(max_score) && max_score >= threshold) {
      assigned_count <- assigned_count + 1
    }
  }
}
cat("前1000个细胞中达到阈值的细胞数:", assigned_count, "/1000\n")

