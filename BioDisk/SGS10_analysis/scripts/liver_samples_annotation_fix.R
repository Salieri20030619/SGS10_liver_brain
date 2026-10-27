# ==================== 肝脏样本增强注释脚本 ====================
# 文件名: ~/SGS10_analysis/scripts/liver_samples_annotation_fix.R
# 用法: Rscript liver_samples_annotation_fix.R sample_name n_cores

args <- commandArgs(trailingOnly = TRUE)
sample_name <- args[1]
n_cores <- as.numeric(args[2])

cat("========================================\n")
cat("增强注释肝脏样本:", sample_name, "\n")
cat("========================================\n\n")

library(Seurat)
library(ggplot2)
library(dplyr)
library(SingleR)
library(celldex)

# 加载已完成的Seurat对象
sample_dir <- file.path("/mnt/BioDisk/wangmingyuan/SGS10_analysis", sample_name)
setwd(sample_dir)

# 加载UMAP后的对象
seurat_obj <- readRDS(paste0(sample_name, "_umap.rds"))
cat("加载对象完成，细胞数:", ncol(seurat_obj), "\n")
cat("当前cluster数:", length(unique(Idents(seurat_obj))), "\n\n")

# ==================== 基因名检查 ====================
cat("========== 基因名检查 ==========\n")
gene_names <- rownames(seurat_obj)
cat("基因名示例:\n")
print(head(gene_names, 30))

# 检查是否有人鼠混合标记
has_grch38 <- any(grepl("^GRCh38-", gene_names))
has_mm10 <- any(grepl("^mm10-", gene_names))
cat("\n检测到GRCh38前缀:", has_grch38, "\n")
cat("检测到mm10前缀:", has_mm10, "\n")

# ==================== 方案1: 直接使用人类marker手动注释 ====================
cat("\n========== 方案1: 肝脏专用marker手动注释 ==========\n")

# 定义肝脏各类细胞的marker（同时包含人类和小鼠同源基因）
liver_markers <- list(
  Hepatocytes = c("ALB", "AFP", "CYP3A4", "APOA1", "TTR", "HNF4A", 
                  "Alb", "Afp", "Cyp3a11", "Apoa1", "Ttr", "Hnf4a"),
  Endothelial = c("PECAM1", "CDH5", "VWF", "KDR", "CLEC4G",
                  "Pecam1", "Cdh5", "Vwf", "Kdr", "Clec4g"),
  HSCs = c("PDGFRA", "PDGFRB", "DES", "RGS5", "ACTA2", "LRAT",
           "Pdgfra", "Pdgfrb", "Des", "Rgs5", "Acta2", "Lrat"),
  Kupffer = c("CD68", "CD163", "VSIG4", "MARCO", "C1QA", "C1QB",
              "Cd68", "Cd163", "Vsig4", "Marco", "C1qa", "C1qb"),
  B_cells = c("CD79A", "MS4A1", "CD19", "BANK1",
              "Cd79a", "Ms4a1", "Cd19", "Bank1"),
  T_cells = c("CD3D", "CD3E", "CD4", "CD8A", "CD8B",
              "Cd3d", "Cd3e", "Cd4", "Cd8a", "Cd8b"),
  NK_cells = c("NKG7", "KLRB1", "KLRD1", "NCR1",
               "Nkg7", "Klrb1", "Klrd1", "Ncr1"),
  Monocytes = c("CD14", "FCGR3A", "LYZ", "S100A8", "S100A9",
                "Cd14", "Fcgr3a", "Lyz", "S100a8", "S100a9"),
  Cholangiocytes = c("KRT19", "EPCAM", "SOX9",
                     "Krt19", "Epcam", "Sox9")
)

# 获取表达数据
expr_data <- GetAssayData(seurat_obj, layer = "data")

# 计算每个cluster的marker得分
cluster_ids <- Idents(seurat_obj)
cluster_levels <- levels(cluster_ids)

# 创建得分矩阵
score_matrix <- matrix(0, nrow = length(cluster_levels), 
                       ncol = length(liver_markers))
rownames(score_matrix) <- cluster_levels
colnames(score_matrix) <- names(liver_markers)

for(cl in cluster_levels) {
  cells <- WhichCells(seurat_obj, idents = cl)
  if(length(cells) > 100) cells <- sample(cells, 100)
  
  for(celltype in names(liver_markers)) {
    markers <- liver_markers[[celltype]]
    markers_present <- intersect(markers, rownames(expr_data))
    
    if(length(markers_present) > 0) {
      # 计算平均表达
      avg_exp <- mean(expr_data[markers_present, cells])
      score_matrix[cl, celltype] <- avg_exp
    }
  }
}

# 为每个cluster分配细胞类型
cluster_assignments <- apply(score_matrix, 1, function(x) {
  max_val <- max(x)
  if(max_val > 0.1) {
    return(names(which.max(x)))
  } else {
    return("Unknown")
  }
})

cat("\n手动注释结果:\n")
print(table(cluster_assignments))

# 将注释应用到Seurat对象
seurat_obj$manual_celltype <- cluster_assignments[as.character(cluster_ids)]
seurat_obj$celltype <- seurat_obj$manual_celltype
Idents(seurat_obj) <- seurat_obj$celltype

# ==================== 方案2: 如果手动注释效果不好，使用参考数据 ====================

cat("\n========== 方案2: 使用参考数据验证 ==========\n")

# 清理基因名
expr_data_clean <- expr_data
rownames(expr_data_clean) <- gsub("^GRCh38-", "", rownames(expr_data_clean))
rownames(expr_data_clean) <- gsub("^mm10-", "", rownames(expr_data_clean))

# 去重
keep_idx <- !duplicated(rownames(expr_data_clean))
expr_data_clean <- expr_data_clean[keep_idx, ]

cat("清理后基因数:", nrow(expr_data_clean), "\n")

# 尝试使用多种参考数据
ref_datasets <- list(
  hpca = HumanPrimaryCellAtlasData(),
  monaco = MonacoImmuneData(),
  blueprint = BlueprintEncodeData()
)

singler_results <- list()

for(ref_name in names(ref_datasets)) {
  cat("\n尝试参考:", ref_name, "\n")
  ref <- ref_datasets[[ref_name]]
  
  common_genes <- intersect(rownames(expr_data_clean), rownames(ref))
  cat("共同基因数:", length(common_genes), "\n")
  
  if(length(common_genes) >= 100) {
    tryCatch({
      res <- SingleR(
        test = expr_data_clean[common_genes, ],
        ref = ref[common_genes, ],
        labels = ref$label.main,
        clusters = Idents(seurat_obj),
        assay.type.test = "logcounts"
      )
      singler_results[[ref_name]] <- res$labels
      cat("注释完成，标签分布:\n")
      print(table(res$labels))
    }, error = function(e) {
      cat("错误:", e$message, "\n")
    })
  }
}

# 投票决定最终标签
if(length(singler_results) > 0) {
  cat("\n========== 综合投票 ==========\n")
  
  cluster_levels <- levels(Idents(seurat_obj))
  vote_matrix <- matrix(NA, nrow = length(cluster_levels), ncol = length(singler_results))
  rownames(vote_matrix) <- cluster_levels
  colnames(vote_matrix) <- names(singler_results)
  
  for(i in 1:length(cluster_levels)) {
    cl <- cluster_levels[i]
    for(j in 1:length(singler_results)) {
      vote_matrix[i, j] <- singler_results[[j]][cl]
    }
  }
  
  final_labels <- apply(vote_matrix, 1, function(x) {
    tbl <- table(x)
    names(which.max(tbl))
  })
  
  cat("投票最终结果:\n")
  print(table(final_labels))
  
  # 映射到肝脏细胞类型
  label_mapping <- c(
    "Hepatocyte" = "Hepatocytes",
    "Hepatocytes" = "Hepatocytes",
    "Kupffer cell" = "Kupffer",
    "Macrophage" = "Kupffer",
    "Macrophages" = "Kupffer",
    "Endothelial cell" = "Endothelial",
    "Endothelial cells" = "Endothelial",
    "Stellate cell" = "Stellate",
    "Cholangiocyte" = "Cholangiocytes",
    "B cell" = "B_cells",
    "T cell" = "T_cells",
    "NK cell" = "NK_cells",
    "Monocyte" = "Monocytes"
  )
  
  mapped_labels <- final_labels
  for(old_name in names(label_mapping)) {
    mapped_labels[final_labels == old_name] <- label_mapping[old_name]
  }
  
  seurat_obj$singler_celltype <- mapped_labels[as.character(Idents(seurat_obj))]
  
  # 比较手动注释和参考注释
  cat("\n手动注释 vs 参考注释:\n")
  print(table(seurat_obj$manual_celltype, seurat_obj$singler_celltype))
  
  # 选择可信度高的注释
  # 如果参考注释不是"Unknown"且得分合理，使用参考注释
  seurat_obj$final_celltype <- seurat_obj$singler_celltype
  seurat_obj$final_celltype[is.na(seurat_obj$final_celltype)] <- seurat_obj$manual_celltype[is.na(seurat_obj$final_celltype)]
} else {
  seurat_obj$final_celltype <- seurat_obj$manual_celltype
}

# 设置最终细胞类型
Idents(seurat_obj) <- seurat_obj$final_celltype

# ==================== 可视化 ====================

cat("\n========== 最终细胞类型分布 ==========\n")
celltype_counts <- table(Idents(seurat_obj))
print(celltype_counts)

# 保存比例
prop_df <- data.frame(
  CellType = names(celltype_counts),
  Count = as.numeric(celltype_counts),
  Percentage = as.numeric(prop.table(celltype_counts) * 100)
)
write.csv(prop_df, paste0(sample_name, "_celltype_proportions_fixed.csv"), row.names = FALSE)

# 绘制UMAP
pdf(paste0(sample_name, "_celltypes_fixed.pdf"), width = 14, height = 10)
p <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, 
             label.size = 4, repel = TRUE) + 
  ggtitle(paste0(sample_name, " Cell Types (Fixed)\n", 
                 "Total cells: ", ncol(seurat_obj)))
print(p)
dev.off()

# ==================== 免疫细胞细分（如果存在多种免疫细胞） ====================

immune_types <- c("B_cells", "T_cells", "NK_cells", "Monocytes", "Kupffer")
immune_present <- immune_types[immune_types %in% unique(Idents(seurat_obj))]

if(length(immune_present) > 0) {
  cat("\n========== 免疫细胞精细分析 ==========\n")
  
  immune_cells <- WhichCells(seurat_obj, expression = final_celltype %in% immune_present)
  cat("免疫细胞数量:", length(immune_cells), "\n")
  
  if(length(immune_cells) >= 100) {
    immune_obj <- subset(seurat_obj, cells = immune_cells)
    
    # 重新聚类
    immune_obj <- NormalizeData(immune_obj, verbose = FALSE)
    immune_obj <- FindVariableFeatures(immune_obj, nfeatures = 1500, verbose = FALSE)
    immune_obj <- ScaleData(immune_obj, verbose = FALSE)
    immune_obj <- RunPCA(immune_obj, npcs = 30, verbose = FALSE)
    immune_obj <- FindNeighbors(immune_obj, dims = 1:20, verbose = FALSE)
    immune_obj <- FindClusters(immune_obj, resolution = 0.5, verbose = FALSE)
    immune_obj <- RunUMAP(immune_obj, dims = 1:20, verbose = FALSE)
    
    # 使用免疫细胞marker精细注释
    immune_markers <- list(
      CD4_T = c("CD4", "CD3D", "CD3E"),
      CD8_T = c("CD8A", "CD8B", "CD3D"),
      B = c("CD19", "MS4A1", "CD79A"),
      NK = c("NKG7", "KLRB1", "KLRD1"),
      Monocyte = c("CD14", "FCGR3A", "LYZ"),
      Macrophage = c("CD68", "CD163", "MRC1")
    )
    
    # 计算每个cluster的得分
    immune_clusters <- Idents(immune_obj)
    immune_cluster_levels <- levels(immune_clusters)
    immune_scores <- matrix(0, nrow = length(immune_cluster_levels), 
                            ncol = length(immune_markers))
    rownames(immune_scores) <- immune_cluster_levels
    colnames(immune_scores) <- names(immune_markers)
    
    for(cl in immune_cluster_levels) {
      cells <- WhichCells(immune_obj, idents = cl)
      if(length(cells) > 100) cells <- sample(cells, 100)
      
      for(celltype in names(immune_markers)) {
        markers <- immune_markers[[celltype]]
        markers_present <- intersect(markers, rownames(immune_obj))
        if(length(markers_present) > 0) {
          avg_exp <- mean(GetAssayData(immune_obj, layer = "data")[markers_present, cells])
          immune_scores[cl, celltype] <- avg_exp
        }
      }
    }
    
    immune_assignments <- apply(immune_scores, 1, function(x) {
      names(which.max(x))
    })
    
    immune_obj$immune_subtype <- immune_assignments[as.character(immune_clusters)]
    
    # 保存免疫细胞子集
    saveRDS(immune_obj, paste0(sample_name, "_immune_subset_fixed.rds"))
    
    # 绘图
    pdf(paste0(sample_name, "_immune_subtypes_fixed.pdf"), width = 14, height = 10)
    p1 <- DimPlot(immune_obj, reduction = "umap", group.by = "immune_subtype", 
                  label = TRUE, label.size = 3) +
      ggtitle(paste0(sample_name, " Immune Cell Subtypes"))
    p2 <- DimPlot(immune_obj, reduction = "umap", label = TRUE) +
      ggtitle("Immune Cell Clusters")
    print(p1 | p2)
    dev.off()
    
    # 免疫细胞比例
    immune_prop <- prop.table(table(immune_obj$immune_subtype)) * 100
    immune_df <- data.frame(
      Subtype = names(immune_prop),
      Percentage = as.numeric(immune_prop)
    )
    write.csv(immune_df, paste0(sample_name, "_immune_subtype_proportions_fixed.csv"), 
              row.names = FALSE)
  }
}

# ==================== 保存最终对象 ====================

saveRDS(seurat_obj, paste0(sample_name, "_complete_fixed.rds"))
cat("\n✓ 最终对象已保存\n")

# ==================== 生成报告 ====================

report_file <- paste0(sample_name, "_annotation_report.txt")
sink(report_file)

cat("========================================\n")
cat("肝脏样本增强注释报告\n")
cat("========================================\n\n")
cat("样本:", sample_name, "\n")
cat("总细胞数:", ncol(seurat_obj), "\n")
cat("总基因数:", nrow(seurat_obj), "\n\n")

cat("========================================\n")
cat("最终细胞类型分布\n")
cat("========================================\n")
for(i in 1:length(celltype_counts)) {
  cat(sprintf("  %-15s: %6d (%5.1f%%)\n", 
              names(celltype_counts)[i], 
              celltype_counts[i], 
              prop.table(celltype_counts)[i] * 100))
}

# 如果手动和参考注释都有，显示一致性
if(exists("singler_results") && length(singler_results) > 0) {
  cat("\n========================================\n")
  cat("注释方法对比\n")
  cat("========================================\n")
  cat("手动注释 vs 参考注释一致性:\n")
  print(table(seurat_obj$manual_celltype, seurat_obj$singler_celltype))
}

sink()

cat("\n✓ 报告已保存:", report_file, "\n")
cat("\n========== 注释完成 ==========\n")
