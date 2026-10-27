# ==================== 肝脏手动注释脚本（脑样本阈值策略 + 只保留有marker的亚型） ====================
# 用法: Rscript liver_manual_annotation_fixed.R sample_name

args <- commandArgs(trailingOnly = TRUE)
sample_name <- args[1]

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

if(!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap", repos = "https://cloud.r-project.org/")
}
library(pheatmap)

cat("========================================\n")
cat("肝脏手动注释（脑样本阈值策略）:", sample_name, "\n")
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
  
  # 肝细胞（主要实质细胞）
  liver_markers$Hepatocytes <- c("Alb", "Afp", "Ttr", "Hnf4a", "Hnf1a", "Foxa1", "Foxa2")
  
  # 肝细胞亚型（基于区域）
  liver_markers$Hepatocytes_Periportal <- c("Cyp2e1", "Cyp1a2", "G6pc", "Pck1", "Ass1")
  liver_markers$Hepatocytes_Perivenous <- c("Cyp3a11", "Cyp2f2", "Glul", "Oat", "Gstm1")
  liver_markers$Hepatocytes_Proliferating <- c("Mki67", "Top2a", "Pcna", "Ccnb1")
  
  # 星状细胞
  liver_markers$Stellate <- c("Pdgfra", "Pdgfrb", "Des", "Rgs5", "Acta2", "Lrat", "Reln")
  liver_markers$HSC_Quiescent <- c("Rbp1", "Lrat", "Pparg")
  liver_markers$HSC_Activated <- c("Acta2", "Col1a1", "Col1a2", "Timp1", "Tgfb1")
  
  # 肝窦内皮细胞
  liver_markers$LSEC <- c("Clec4g", "Lyve1", "Stab1", "Stab2", "Fcgr2b", "Oit3")
  liver_markers$Endothelial <- c("Pecam1", "Cdh5", "Vwf", "Kdr", "Flt1")
  
  # 枯否细胞（肝脏巨噬细胞）
  liver_markers$Kupffer <- c("Clec4f", "Vsig4", "Timd4", "C1qa", "C1qb", "C1qc", "Cd68", "Cd163")
  
  # 巨噬细胞亚型
  liver_markers$Macrophage_M1 <- c("Il1b", "Tnfa", "Cd86", "Nos2")
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
  liver_markers$B_Plasma <- c("Igba", "Jchain", "Xbp1", "Irf4")
  
  # NK细胞
  liver_markers$NK_cells <- c("Nkg7", "Klrb1", "Klrd1", "Klrk1")
  
  # 单核细胞
  liver_markers$Monocytes <- c("Cd14", "Ly6c2", "Ccr2", "S100a8", "S100a9")
  
  # 树突状细胞
  liver_markers$Dendritic <- c("Flt3", "Zbtb46", "Clec9a", "Cd11c")
  
  # 中性粒细胞
  liver_markers$Neutrophils <- c("S100a8", "S100a9", "Cxcl2", "Cxcr2", "Mmp8")
  
  # 周细胞
  liver_markers$Pericytes <- c("Rgs5", "Pdgfrb", "Cspg4", "Kcnj8")
  
  # 红细胞前体
  liver_markers$Erythrocytes <- c("Hba-a1", "Hba-a2", "Hbb-bt", "Klf1")
  
  # 增殖细胞
  liver_markers$Proliferating <- c("Mki67", "Top2a", "Pcna", "Ccnb1")
  
  return(liver_markers)
}

# 创建标记数据库
liver_markers <- create_liver_marker_db()
cat("定义了", length(liver_markers), "种细胞类型\n")

# 定义所有可能的亚型和大类
all_subtype_names <- c(
  # 肝细胞亚型
  "Hepatocytes_Periportal", "Hepatocytes_Perivenous", "Hepatocytes_Proliferating",
  # 星状细胞亚型
  "HSC_Quiescent", "HSC_Activated",
  # 巨噬细胞亚型
  "Macrophage_M1", "Macrophage_M2", "Macrophage_NAM",
  # 免疫细胞亚型
  "CD4_T", "CD8_T", "Treg", "B_Plasma",
  # 内皮细胞亚型
  "LSEC"
)

all_major_names <- c(
  "Hepatocytes", "Stellate", "Endothelial", "Kupffer", 
  "T_cells", "B_cells", "Monocytes", "Neutrophils", 
  "Cholangiocytes", "Pericytes", "NK_cells"
)

# 获取表达数据
cat("获取表达数据...\n")
expr_data <- GetAssayData(seurat_obj, layer = "data")

# 计算每个细胞类型的marker得分
cat("计算每个细胞类型的marker得分...\n")

all_scores <- matrix(0, nrow = ncol(seurat_obj), ncol = length(liver_markers))
colnames(all_scores) <- names(liver_markers)
rownames(all_scores) <- colnames(seurat_obj)

# 记录哪些类型有marker
has_marker <- rep(FALSE, length(liver_markers))
names(has_marker) <- names(liver_markers)

for(ct in names(liver_markers)) {
  markers <- liver_markers[[ct]]
  markers_present <- intersect(markers, rownames(expr_data))
  if(length(markers_present) > 0) {
    scores <- colMeans(expr_data[markers_present, , drop = FALSE])
    all_scores[, ct] <- scores
    has_marker[ct] <- TRUE
    cat("  ✓", ct, ": 找到", length(markers_present), "个marker\n")
  } else {
    cat("  ⚠️", ct, ": 未找到marker\n")
  }
}

# ========== 只保留有marker的亚型和大类 ==========
subtype_names <- all_subtype_names[sapply(all_subtype_names, function(x) has_marker[x])]
major_names <- all_major_names[sapply(all_major_names, function(x) has_marker[x])]

cat("\n找到", length(subtype_names), "种亚型（有marker）:\n")
if(length(subtype_names) > 0) {
  for(ct in subtype_names) cat("  -", ct, "\n")
}

cat("\n找到", length(major_names), "种大类（有marker）:\n")
if(length(major_names) > 0) {
  for(ct in major_names) cat("  -", ct, "\n")
}

# 提取得分矩阵
if(length(subtype_names) > 0) {
  subtype_scores <- all_scores[, subtype_names, drop = FALSE]
  cat("\n亚型得分矩阵维度:", dim(subtype_scores), "\n")
} else {
  cat("\n警告: 没有找到任何亚型标记基因，将仅使用大类注释\n")
  subtype_scores <- matrix(0, nrow = ncol(seurat_obj), ncol = 1)
  subtype_names <- character(0)
}

if(length(major_names) > 0) {
  major_scores <- all_scores[, major_names, drop = FALSE]
  cat("大类得分矩阵维度:", dim(major_scores), "\n")
} else {
  stop("错误: 没有找到任何大类标记基因")
}

# ========== 使用脑样本的阈值策略 ==========
cat("\n计算阈值（脑样本策略）...\n")
all_max_scores <- apply(all_scores, 1, max)

# 亚型阈值：使用所有细胞最大得分的5%分位数
if(length(subtype_names) > 0) {
  valid_scores <- all_max_scores[all_max_scores > 0 & is.finite(all_max_scores)]
  if(length(valid_scores) > 0) {
    subtype_threshold <- quantile(valid_scores, 0.05, na.rm = TRUE)
    if(is.na(subtype_threshold) || subtype_threshold < 0.03) subtype_threshold <- 0.03
  } else {
    subtype_threshold <- 0.05
  }
} else {
  subtype_threshold <- Inf
}

# 大类阈值：使用所有细胞最大得分的2%分位数
valid_scores <- all_max_scores[all_max_scores > 0 & is.finite(all_max_scores)]
if(length(valid_scores) > 0) {
  major_threshold <- quantile(valid_scores, 0.02, na.rm = TRUE)
  if(is.na(major_threshold) || major_threshold < 0.02) major_threshold <- 0.02
} else {
  major_threshold <- 0.05
}

cat("  亚型阈值 (0.05分位数):", round(subtype_threshold, 4), "\n")
cat("  大类阈值 (0.02分位数):", round(major_threshold, 4), "\n")

# ========== 亚型优先分配策略（修复版） ==========
cat("\n亚型优先分配...\n")
final_assignment <- rep(NA, ncol(seurat_obj))
names(final_assignment) <- colnames(seurat_obj)

for(cell in colnames(seurat_obj)) {
  assigned <- FALSE
  
  # 第一步：检查亚型
  if(length(subtype_names) > 0) {
    cell_subtype_scores <- subtype_scores[cell, , drop = FALSE]
    if(length(cell_subtype_scores) > 0) {
      # 移除NA和无效值
      valid_idx <- !is.na(cell_subtype_scores) & is.finite(cell_subtype_scores)
      if(any(valid_idx)) {
        valid_scores <- cell_subtype_scores[valid_idx]
        max_subtype_score <- max(valid_scores)
        if(is.finite(max_subtype_score) && max_subtype_score >= subtype_threshold) {
          best_subtype <- names(valid_scores)[which.max(valid_scores)]
          final_assignment[cell] <- best_subtype
          assigned <- TRUE
        }
      }
    }
  }
  
  # 第二步：检查大类
  if(!assigned && length(major_names) > 0) {
    cell_major_scores <- major_scores[cell, , drop = FALSE]
    if(length(cell_major_scores) > 0) {
      # 移除NA和无效值
      valid_idx <- !is.na(cell_major_scores) & is.finite(cell_major_scores)
      if(any(valid_idx)) {
        valid_scores <- cell_major_scores[valid_idx]
        max_major_score <- max(valid_scores)
        if(is.finite(max_major_score) && max_major_score >= major_threshold) {
          best_major <- names(valid_scores)[which.max(valid_scores)]
          final_assignment[cell] <- best_major
          assigned <- TRUE
        }
      }
    }
  }
}

# 分配结果统计
cat("\n分配结果统计:\n")
assigned_count <- sum(!is.na(final_assignment))
cat("  已分配细胞数:", assigned_count, "/", ncol(seurat_obj), "\n")
cat("  未分配细胞数:", sum(is.na(final_assignment)), "\n")

# 用SingleR结果填补NA
na_idx <- is.na(final_assignment)
if(sum(na_idx) > 0 && "celltype" %in% colnames(seurat_obj@meta.data)) {
  cat("用SingleR结果填补", sum(na_idx), "个未分配细胞...\n")
  final_assignment[na_idx] <- as.character(seurat_obj$celltype[na_idx])
}

# 应用注释
seurat_obj$final_celltype <- final_assignment
Idents(seurat_obj) <- seurat_obj$final_celltype

cat("\n最终细胞类型分布:\n")
celltype_counts <- table(Idents(seurat_obj))
print(celltype_counts)

# 统计亚型和大类
final_types <- names(celltype_counts)
subtypes_found <- intersect(final_types, subtype_names)
majors_found <- intersect(final_types, major_names)

cat("\n发现的精细亚型 (", length(subtypes_found), "种):\n")
if(length(subtypes_found) > 0) print(subtypes_found)

cat("\n发现的大类 (", length(majors_found), "种):\n")
if(length(majors_found) > 0) print(majors_found)

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
seurat_obj$subtype_max_score <- all_max_scores
seurat_obj$subtype_threshold <- subtype_threshold
seurat_obj$major_threshold <- major_threshold

saveRDS(seurat_obj, paste0(sample_name, "_manual_annotated.rds"))
cat("✓ 手动注释结果已保存\n")

cat("\n========== 肝脏注释完成 ==========\n")
