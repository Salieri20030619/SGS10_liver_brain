# ==================== 肝脏手动注释脚本（保存修复版） ====================
# 用法: Rscript liver_manual_annotation_save_fix.R sample_name

args <- commandArgs(trailingOnly = TRUE)
sample_name <- args[1]

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

cat("========================================\n")
cat("肝脏手动注释:", sample_name, "\n")
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
  
  # 肝细胞亚型
  liver_markers$Hepatocytes_Periportal <- c("Cyp2e1", "Cyp1a2", "G6pc", "Pck1", "Ass1")
  liver_markers$Hepatocytes_Perivenous <- c("Cyp3a11", "Cyp2f2", "Glul", "Oat", "Gstm1")
  liver_markers$Hepatocytes_Proliferating <- c("Mki67", "Top2a", "Pcna", "Ccnb1")
  
  # 星状细胞
  liver_markers$Stellate <- c("Pdgfra", "Pdgfrb", "Des", "Rgs5", "Acta2", "Lrat", "Reln")
  liver_markers$HSC_Quiescent <- c("Rbp1", "Lrat", "Pparg")
  liver_markers$HSC_Activated <- c("Acta2", "Col1a1", "Col1a2", "Timp1", "Tgfb1")
  
  # 内皮细胞
  liver_markers$Endothelial <- c("Pecam1", "Cdh5", "Vwf", "Kdr", "Flt1")
  liver_markers$LSEC <- c("Clec4g", "Lyve1", "Stab1", "Stab2", "Fcgr2b", "Oit3")
  
  # 枯否细胞
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
  
  # 技术噪声
  liver_markers$TECHO <- c("Malat1", "Neat1", "Mt-rnr1")
  
  return(liver_markers)
}

# 创建标记数据库
liver_markers <- create_liver_marker_db()
cat("定义了", length(liver_markers), "种细胞类型\n")

# 定义大类
major_names <- c(
  "Hepatocytes", "Stellate", "Endothelial", "Kupffer", 
  "T_cells", "B_cells", "Monocytes", "Neutrophils", 
  "Cholangiocytes", "Pericytes", "NK_cells",
  "Macrophage_M2", "Macrophage_NAM", "LSEC", "HSC_Activated",
  "Dendritic", "Erythrocytes", "Proliferating", "TECHO"
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

# 计算每个细胞的最高得分
all_max_scores <- apply(all_scores, 1, max, na.rm = TRUE)
positive_scores <- all_max_scores[all_max_scores > 0 & is.finite(all_max_scores)]

if(length(positive_scores) > 0) {
  threshold <- quantile(positive_scores, 0.05, na.rm = TRUE)
  if(is.na(threshold) || threshold < 0.01) threshold <- 0.01
} else {
  threshold <- 0.01
}

cat("\n阈值 (正得分5%分位数):", round(threshold, 4), "\n")

# 分配细胞类型
cat("\n分配细胞类型...\n")
final_assignment <- rep("Other", ncol(seurat_obj))
names(final_assignment) <- colnames(seurat_obj)

for(cell in colnames(seurat_obj)) {
  cell_scores <- all_scores[cell, , drop = FALSE]
  pos_scores <- cell_scores[cell_scores > 0 & is.finite(cell_scores)]
  
  if(length(pos_scores) > 0) {
    max_score <- max(pos_scores)
    
    if(max_score >= threshold) {
      best_idx <- which.max(cell_scores)
      best_type <- colnames(all_scores)[best_idx]
      
      if(!is.null(best_type) && length(best_type) > 0 && best_type != "") {
        final_assignment[cell] <- best_type
      }
    }
  }
}

# 分配结果统计
assigned_count <- sum(final_assignment != "Other")
cat("\n分配结果统计:\n")
cat("  已分配细胞数:", assigned_count, "/", ncol(seurat_obj), "\n")
cat("  Other细胞数:", sum(final_assignment == "Other"), "\n")

# 用SingleR结果填补Other
other_idx <- final_assignment == "Other"
if(sum(other_idx) > 0 && "celltype" %in% colnames(seurat_obj@meta.data)) {
  cat("用SingleR结果填补", sum(other_idx), "个Other细胞...\n")
  final_assignment[other_idx] <- as.character(seurat_obj$celltype[other_idx])
}

# 确保所有细胞都有注释
na_idx <- is.na(final_assignment)
if(sum(na_idx) > 0) {
  final_assignment[na_idx] <- "Unknown"
}

# 应用注释 - 关键修复
seurat_obj$final_celltype <- final_assignment

# 确保所有细胞的final_celltype都有值
if(any(is.na(seurat_obj$final_celltype))) {
  cat("警告: 仍有NA值，用Unknown填充\n")
  seurat_obj$final_celltype[is.na(seurat_obj$final_celltype)] <- "Unknown"
}

# 设置Idents
Idents(seurat_obj) <- seurat_obj$final_celltype

cat("\n最终细胞类型分布:\n")
celltype_counts <- table(Idents(seurat_obj))
print(celltype_counts)

# 统计发现的细胞类型
cat("\n发现的细胞类型 (", length(celltype_counts), "种):\n")
print(sort(celltype_counts, decreasing = TRUE))

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
  ggtitle(paste0(sample_name, " - Manual Annotation"))
if("celltype" %in% colnames(seurat_obj@meta.data)) {
  p2 <- DimPlot(seurat_obj, reduction = "umap", group.by = "celltype", 
                label = TRUE, label.size = 3, repel = TRUE) +
    ggtitle(paste0(sample_name, " - SingleR Annotation"))
  print(p1 | p2)
} else {
  print(p1)
}
dev.off()
cat("✓ 注释对比图已保存\n")

# 保存对象 - 关键修复
cat("\n保存对象...\n")
# 确保在保存前设置正确的Idents
Idents(seurat_obj) <- seurat_obj$final_celltype

# 保存阈值信息
seurat_obj$annotation_threshold <- threshold

# 保存
saveRDS(seurat_obj, paste0(sample_name, "_manual_annotated.rds"))

# 验证保存
if(file.exists(paste0(sample_name, "_manual_annotated.rds"))) {
  file_size <- file.size(paste0(sample_name, "_manual_annotated.rds")) / 1024^3
  cat("✓ 手动注释结果已保存，大小:", round(file_size, 2), "GB\n")
  
  # 验证可以重新加载
  test_obj <- tryCatch({
    readRDS(paste0(sample_name, "_manual_annotated.rds"))
  }, error = function(e) {
    cat("  警告: 保存的文件无法读取 -", e$message, "\n")
    return(NULL)
  })
  
  if(!is.null(test_obj)) {
    cat("  验证: 文件可正常读取\n")
    if("final_celltype" %in% colnames(test_obj@meta.data)) {
      cat("  验证: final_celltype存在\n")
    }
  }
} else {
  cat("✗ 保存失败！文件未生成\n")
}

cat("\n========== 肝脏注释完成 ==========\n")
