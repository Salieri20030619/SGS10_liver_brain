# ==================== 肝脏手动注释脚本（完整调试版） ====================
# 用法: Rscript liver_manual_annotation_debug_full.R sample_name

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
  liver_markers$Hepatocytes_Proliferating <- c("Mki67", "Top2a", "Pcna", "Ccnb1")
  liver_markers$Stellate <- c("Pdgfra", "Pdgfrb", "Des", "Rgs5", "Acta2", "Lrat", "Reln")
  liver_markers$HSC_Quiescent <- c("Rbp1", "Lrat", "Pparg")
  liver_markers$HSC_Activated <- c("Acta2", "Col1a1", "Col1a2", "Timp1", "Tgfb1")
  liver_markers$Endothelial <- c("Pecam1", "Cdh5", "Vwf", "Kdr", "Flt1")
  liver_markers$LSEC <- c("Clec4g", "Lyve1", "Stab1", "Stab2", "Fcgr2b", "Oit3")
  liver_markers$Kupffer <- c("Clec4f", "Vsig4", "Timd4", "C1qa", "C1qb", "C1qc", "Cd68", "Cd163")
  liver_markers$Macrophage_M1 <- c("Il1b", "Tnfa", "Cd86", "Nos2")
  liver_markers$Macrophage_M2 <- c("Arg1", "Cd163", "Mrc1", "Il10")
  liver_markers$Macrophage_NAM <- c("Trem2", "Gpnmb", "Cd9", "Spp1")
  liver_markers$Cholangiocytes <- c("Krt7", "Krt8", "Krt18", "Krt19", "Epcam", "Sox9")
  liver_markers$T_cells <- c("Cd3d", "Cd3e", "Cd3g", "Cd4", "Cd8a")
  liver_markers$CD4_T <- c("Cd4", "Cd3e", "Il2ra")
  liver_markers$CD8_T <- c("Cd8a", "Cd8b", "Cd3e")
  liver_markers$Treg <- c("Foxp3", "Il2ra", "Ctla4")
  liver_markers$B_cells <- c("Cd19", "Ms4a1", "Cd79a", "Cd79b")
  liver_markers$B_Plasma <- c("Igba", "Jchain", "Xbp1", "Irf4")
  liver_markers$NK_cells <- c("Nkg7", "Klrb1", "Klrd1", "Klrk1")
  liver_markers$Monocytes <- c("Cd14", "Ly6c2", "Ccr2", "S100a8", "S100a9")
  liver_markers$Dendritic <- c("Flt3", "Zbtb46", "Clec9a", "Cd11c")
  liver_markers$Neutrophils <- c("S100a8", "S100a9", "Cxcl2", "Cxcr2", "Mmp8")
  liver_markers$Pericytes <- c("Rgs5", "Pdgfrb", "Cspg4", "Kcnj8")
  liver_markers$Erythrocytes <- c("Hba-a1", "Hba-a2", "Hbb-bt", "Klf1")
  liver_markers$Proliferating <- c("Mki67", "Top2a", "Pcna", "Ccnb1")
  liver_markers$TECHO <- c("Malat1", "Neat1", "Mt-rnr1")
  
  return(liver_markers)
}

# 创建标记数据库
liver_markers <- create_liver_marker_db()
cat("定义了", length(liver_markers), "种细胞类型\n")

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
final_assignment[is.na(final_assignment)] <- "Unknown"

# ========== 详细调试信息 ==========
cat("\n========== 调试信息 ==========\n")
cat("1. final_assignment 类型:", class(final_assignment), "\n")
cat("2. final_assignment 长度:", length(final_assignment), "\n")
cat("3. final_assignment 唯一值数量:", length(unique(final_assignment)), "\n")
cat("4. final_assignment 前10个值:\n")
print(head(final_assignment, 10))

cat("\n5. 元数据当前列名:\n")
print(colnames(seurat_obj@meta.data))

cat("\n6. 当前Idents:\n")
print(head(Idents(seurat_obj), 10))

# 尝试添加元数据
cat("\n7. 尝试添加 final_celltype 到元数据...\n")
tryCatch({
  seurat_obj$final_celltype <- final_assignment
  cat("  ✓ 成功添加 final_celltype\n")
}, error = function(e) {
  cat("  ✗ 添加失败:", e$message, "\n")
})

# 验证添加是否成功
cat("\n8. 验证添加结果:\n")
cat("  final_celltype 存在:", "final_celltype" %in% colnames(seurat_obj@meta.data), "\n")
if("final_celltype" %in% colnames(seurat_obj@meta.data)) {
  cat("  final_celltype 长度:", length(seurat_obj$final_celltype), "\n")
  cat("  final_celltype 前5个值:\n")
  print(head(seurat_obj$final_celltype, 5))
}

# 尝试设置Idents
cat("\n9. 尝试设置 Idents...\n")
tryCatch({
  Idents(seurat_obj) <- seurat_obj$final_celltype
  cat("  ✓ 成功设置 Idents\n")
}, error = function(e) {
  cat("  ✗ 设置 Idents 失败:", e$message, "\n")
})

# 验证Idents
cat("\n10. 验证 Idents:\n")
cat("  Idents 类型:", class(Idents(seurat_obj)), "\n")
cat("  Idents 长度:", length(Idents(seurat_obj)), "\n")
cat("  Idents 唯一值数量:", length(unique(Idents(seurat_obj))), "\n")
cat("  Idents 前10个值:\n")
print(head(Idents(seurat_obj), 10))

# 检查一致性
cat("\n11. 检查一致性:\n")
if("final_celltype" %in% colnames(seurat_obj@meta.data)) {
  ident_check <- all(as.character(Idents(seurat_obj)) == seurat_obj$final_celltype)
  cat("  Idents 与 final_celltype 一致:", ident_check, "\n")
  if(!ident_check) {
    cat("  差异数量:", sum(as.character(Idents(seurat_obj)) != seurat_obj$final_celltype), "\n")
  }
}

# 保存比例
prop_df <- data.frame(
  CellType = names(table(Idents(seurat_obj))),
  Count = as.numeric(table(Idents(seurat_obj))),
  Percentage = as.numeric(prop.table(table(Idents(seurat_obj))) * 100)
)
write.csv(prop_df, paste0(sample_name, "_manual_proportions.csv"), row.names = FALSE)
cat("\n✓ 比例已保存\n")

# 绘制对比图
cat("\n12. 绘制对比图...\n")
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

# 保存对象
cat("\n13. 保存对象...\n")
tryCatch({
  # 确保Idents有效
  Idents(seurat_obj) <- seurat_obj$final_celltype
  
  # 添加阈值信息
  seurat_obj$annotation_threshold <- threshold
  
  # 保存
  saveRDS(seurat_obj, paste0(sample_name, "_manual_annotated.rds"))
  cat("  ✓ 保存成功\n")
}, error = function(e) {
  cat("  ✗ 保存失败:", e$message, "\n")
  traceback()
})

# 验证保存
cat("\n14. 验证保存...\n")
if(file.exists(paste0(sample_name, "_manual_annotated.rds"))) {
  file_size <- file.size(paste0(sample_name, "_manual_annotated.rds")) / 1024^3
  cat("  ✓ 文件存在，大小:", round(file_size, 2), "GB\n")
  
  test_obj <- tryCatch({
    readRDS(paste0(sample_name, "_manual_annotated.rds"))
  }, error = function(e) {
    cat("  ✗ 无法读取:", e$message, "\n")
    return(NULL)
  })
  
  if(!is.null(test_obj)) {
    cat("  ✓ 文件可正常读取\n")
    if("final_celltype" %in% colnames(test_obj@meta.data)) {
      cat("  ✓ final_celltype存在\n")
    } else {
      cat("  ✗ final_celltype不存在\n")
    }
  }
} else {
  cat("  ✗ 文件不存在\n")
}

cat("\n========== 调试完成 ==========\n")
