# ==================== 肝脏样本修复脚本（完全本地化，不依赖网络） ====================
# 文件名: ~/SGS10_analysis/scripts/liver_samples_fix_local.R
# 用法: Rscript liver_samples_fix_local.R sample_name n_cores

args <- commandArgs(trailingOnly = TRUE)
sample_name <- args[1]
n_cores <- as.numeric(args[2])

cat("========================================\n")
cat("肝脏样本修复脚本（完全本地化）:", sample_name, "\n")
cat("========================================\n\n")

library(Seurat)
library(ggplot2)
library(dplyr)
library(future)
library(Matrix)

# 设置并行
plan("multicore", workers = n_cores)
options(future.globals.maxSize = 50 * 1024^3)
set.seed(42)

sample_dir <- file.path("/mnt/BioDisk/wangmingyuan/SGS10_analysis", sample_name)
if(!dir.exists(sample_dir)) dir.create(sample_dir, recursive = TRUE)
setwd(sample_dir)
cat("工作目录:", getwd(), "\n\n")

# ==================== 数据路径 ====================
data_paths <- list(
  gan56 = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/gan56/gan_56_filtered_feature_bc_matrix",
  liver3 = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/liver_samples/25090559_Liver_3_expression/Liver_3_filtered_feature_bc_matrix",
  liver8 = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/liver_samples/25090559_Liver_8_expression/Liver_8_filtered_feature_bc_matrix",
  bla_liver = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/liver_samples/25090559_Bla_liver_1_expression/Bla_liver_1_filtered_feature_bc_matrix",
  con_liver = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/liver_samples/25090559_Con_liver_1_expression/Con_liver_1_filtered_feature_bc_matrix"
)

# ==================== 创建本地肝脏marker数据库 ====================
create_liver_marker_db <- function() {
  # 肝脏各类细胞的关键marker（人类）
  marker_db <- list(
    # 肝细胞 - 核心marker
    Hepatocytes = c("ALB", "AFP", "CYP3A4", "APOA1", "TTR", "HNF4A", "PCK1", "G6PC"),
    
    # 内皮细胞
    Endothelial = c("PECAM1", "CDH5", "VWF", "KDR", "CLEC4G", "LYVE1"),
    
    # 星状细胞
    HSCs = c("PDGFRA", "PDGFRB", "DES", "RGS5", "ACTA2", "LRAT"),
    
    # 枯否细胞（肝脏巨噬细胞）
    Kupffer = c("CD68", "CD163", "VSIG4", "MARCO", "C1QA", "C1QB"),
    
    # T细胞
    T_cells = c("CD3D", "CD3E", "CD4", "CD8A", "CD8B"),
    
    # B细胞
    B_cells = c("CD19", "MS4A1", "CD79A"),
    
    # NK细胞
    NK_cells = c("NKG7", "KLRB1", "KLRD1", "KLRK1"),
    
    # 单核细胞/巨噬细胞
    Monocytes = c("CD14", "FCGR3A", "LYZ", "S100A8", "S100A9"),
    
    # 胆管细胞
    Cholangiocytes = c("KRT19", "EPCAM", "SOX9", "KRT7"),
    
    # 中性粒细胞
    Neutrophils = c("FCGR3B", "CSF3R", "CXCL8", "S100A12"),
    
    # 肥大细胞
    Mast_cells = c("TPSAB1", "TPSB2", "CPA3", "MS4A2")
  )
  
  # 合并所有marker
  all_markers <- unique(unlist(marker_db))
  
  # 创建从基因名到细胞类型的映射
  gene_to_celltype <- list()
  for(celltype in names(marker_db)) {
    for(gene in marker_db[[celltype]]) {
      gene_to_celltype[[gene]] <- c(gene_to_celltype[[gene]], celltype)
    }
  }
  
  return(list(marker_db = marker_db, all_markers = all_markers, gene_to_celltype = gene_to_celltype))
}

liver_db <- create_liver_marker_db()
cat("肝脏marker数据库已创建，包含", length(liver_db$all_markers), "个关键基因\n")
cat("细胞类型:", paste(names(liver_db$marker_db), collapse = ", "), "\n\n")

# ==================== 步骤1: 读取features文件 ====================
cat("========== 步骤1: 读取features文件 ==========\n")

data_dir <- data_paths[[sample_name]]
features_file <- file.path(data_dir, "features.tsv.gz")

if(!file.exists(features_file)) {
  stop("features.tsv.gz not found at: ", features_file)
}

# 读取features文件
features <- read.table(features_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(features) <- c("gene_id", "gene_name", "feature_type")

cat("Features文件维度:", dim(features), "\n")
cat("Feature类型分布:\n")
print(table(features$feature_type))

# 只保留Gene Expression类型
gene_features <- features[features$feature_type == "Gene Expression", ]
cat("\nGene Expression基因数:", nrow(gene_features), "\n")

# ==================== 步骤2: 提取标准基因名 ====================
cat("\n========== 步骤2: 提取标准基因名 ==========\n")

# 方法1: 从gene_name中提取（去除GRCh38-前缀）
extract_standard_name <- function(gene_name) {
  # 去除GRCh38-和mm10-前缀
  cleaned <- gsub("^GRCh38-", "", gene_name)
  cleaned <- gsub("^mm10-", "", cleaned)
  return(cleaned)
}

# 方法2: 检查是否是已知的肝脏marker
is_liver_marker <- function(gene_name) {
  # 去除前缀后检查
  cleaned <- extract_standard_name(gene_name)
  return(cleaned %in% liver_db$all_markers)
}

# 处理每个基因
standard_names <- character(nrow(gene_features))
keep_idx <- logical(nrow(gene_features))

for(i in 1:nrow(gene_features)) {
  gene_id <- gene_features$gene_id[i]
  gene_name <- gene_features$gene_name[i]
  
  # 尝试从gene_name提取
  std_name <- extract_standard_name(gene_name)
  
  # 如果提取后是空字符串或太长（>50），尝试从gene_id提取
  if(std_name == "" || nchar(std_name) > 50) {
    # 从gene_id提取（如果是Ensembl格式）
    if(grepl("ENSG", gene_id)) {
      std_name <- gsub("\\..*$", "", gene_id)
    } else {
      std_name <- gene_id
    }
  }
  
  standard_names[i] <- std_name
  
  # 保留标准基因名不为空且不全是数字/符号的基因
  if(std_name != "" && !grepl("^[0-9_]+$", std_name)) {
    keep_idx[i] <- TRUE
  }
}

# 应用过滤
gene_features_filtered <- gene_features[keep_idx, ]
standard_names_filtered <- standard_names[keep_idx]

cat("过滤后基因数:", length(standard_names_filtered), "\n")
cat("标准基因名示例:\n")
print(head(standard_names_filtered, 30))

# 检查肝脏marker
found_markers <- intersect(standard_names_filtered, liver_db$all_markers)
cat("\n找到的肝脏marker数:", length(found_markers), "\n")
if(length(found_markers) > 0) {
  cat("找到的肝脏marker:", paste(head(found_markers, 20), collapse = ", "), "\n")
}

# ==================== 步骤3: 读取数据并应用标准基因名 ====================
cat("\n========== 步骤3: 读取10X数据 ==========\n")

# 读取原始数据
cat("读取10X数据...\n")
data <- Read10X(data_dir)
cat("原始数据维度:", dim(data), "\n")

# 确保数据与features文件顺序一致
if(nrow(data) == nrow(features)) {
  # 应用过滤
  data_filtered <- data[keep_idx, ]
  rownames(data_filtered) <- standard_names_filtered
  
  # 去重（保留表达量最高的）
  dup_genes <- unique(rownames(data_filtered)[duplicated(rownames(data_filtered))])
  if(length(dup_genes) > 0) {
    cat("处理重复基因:", length(dup_genes), "\n")
    # 对于重复基因，保留表达量总和最高的
    gene_max_exp <- data.frame(gene = rownames(data_filtered), 
                               max_exp = apply(data_filtered, 1, max))
    
    keep_dup <- c()
    for(gene in dup_genes) {
      idx <- which(rownames(data_filtered) == gene)
      # 保留表达量最高的
      max_idx <- idx[which.max(gene_max_exp$max_exp[idx])]
      keep_dup <- c(keep_dup, max_idx)
    }
    
    # 保留非重复基因和选中的重复基因
    non_dup_idx <- which(!duplicated(rownames(data_filtered)))
    final_idx <- sort(unique(c(non_dup_idx, keep_dup)))
    data_filtered <- data_filtered[final_idx, ]
    
    cat("去重后基因数:", nrow(data_filtered), "\n")
  }
  
  cat("过滤后数据维度:", dim(data_filtered), "\n")
  cat("示例基因名:\n")
  print(head(rownames(data_filtered), 30))
  
  # 再次检查肝脏marker
  found_markers_after <- intersect(rownames(data_filtered), liver_db$all_markers)
  cat("\n过滤后找到的肝脏marker:", paste(head(found_markers_after, 20), collapse = ", "), "\n")
  
  if(length(found_markers_after) < 5) {
    cat("\n⚠️ 警告：找到的肝脏marker较少，可能影响注释准确性\n")
    cat("建议：使用更宽松的过滤条件或检查数据质量\n")
  }
  
  # ==================== 步骤4: 创建Seurat对象 ====================
  cat("\n========== 步骤4: 创建Seurat对象 ==========\n")
  
  seurat_obj <- CreateSeuratObject(counts = data_filtered, project = sample_name)
  
  # 计算线粒体基因比例
  mt_genes <- grep("^MT-", rownames(seurat_obj), value = TRUE, ignore.case = TRUE)
  if(length(mt_genes) > 0) {
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  } else {
    seurat_obj[["percent.mt"]] <- 0
  }
  
  cat("✓ 创建成功！细胞数:", ncol(seurat_obj), "\n")
  cat("基因数:", nrow(seurat_obj), "\n")
  
  # 保存原始对象
  saveRDS(seurat_obj, paste0(sample_name, "_raw_fixed.rds"))
  cat("✓ 原始对象已保存\n")
  
  # ==================== 步骤5: 质控 ====================
  cat("\n========== 步骤5: 质控 ==========\n")
  
  cat("质控前统计:\n")
  cat("  细胞数:", ncol(seurat_obj), "\n")
  cat("  基因数范围:", range(seurat_obj$nFeature_RNA), "\n")
  cat("  UMI数范围:", range(seurat_obj$nCount_RNA), "\n")
  
  # 质控
  seurat_obj <- subset(seurat_obj, 
                       subset = nFeature_RNA > 200 & 
                                nFeature_RNA < 6000 & 
                                percent.mt < 20)
  
  cat("\n质控后统计:\n")
  cat("  细胞数:", ncol(seurat_obj), "\n")
  cat("  基因数范围:", range(seurat_obj$nFeature_RNA), "\n")
  
  if(ncol(seurat_obj) < 100) {
    stop("质控后细胞数不足100，无法继续分析")
  }
  
  # 绘制质控图
  pdf(paste0(sample_name, "_QC_fixed.pdf"), width = 12, height = 6)
  print(VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                ncol = 3, pt.size = 0.1))
  dev.off()
  cat("✓ 质控图已保存\n")
  
  # ==================== 步骤6: 标准化和降维 ====================
  cat("\n========== 步骤6: 标准化和降维 ==========\n")
  
  # 标准化
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 2000, verbose = FALSE)
  seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
  
  # PCA
  seurat_obj <- RunPCA(seurat_obj, npcs = 50, verbose = FALSE)
  
  # 肘部图
  pdf(paste0(sample_name, "_elbow_fixed.pdf"), width = 8, height = 6)
  print(ElbowPlot(seurat_obj, ndims = 50))
  dev.off()
  
  # 聚类
  pc_use <- 1:min(30, ncol(seurat_obj[["pca"]]@cell.embeddings))
  seurat_obj <- FindNeighbors(seurat_obj, dims = pc_use, verbose = FALSE)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.8, verbose = FALSE)
  cat("找到", length(unique(Idents(seurat_obj))), "个cluster\n")
  
  # UMAP
  seurat_obj <- RunUMAP(seurat_obj, dims = pc_use, verbose = FALSE)
  
  pdf(paste0(sample_name, "_clusters_fixed.pdf"), width = 10, height = 8)
  print(DimPlot(seurat_obj, reduction = "umap", label = TRUE) + 
          ggtitle(paste0(sample_name, " cell clusters (", ncol(seurat_obj), " cells)")))
  dev.off()
  cat("✓ UMAP图已保存\n")
  
  # ==================== 步骤7: 基于marker的手动注释 ====================
  cat("\n========== 步骤7: 基于marker的手动注释 ==========\n")
  
  expr_data <- GetAssayData(seurat_obj, layer = "data")
  cluster_ids <- Idents(seurat_obj)
  cluster_levels <- levels(cluster_ids)
  
  # 为每个cluster计算每种细胞类型的得分
  score_matrix <- matrix(0, nrow = length(cluster_levels), 
                         ncol = length(liver_db$marker_db))
  rownames(score_matrix) <- cluster_levels
  colnames(score_matrix) <- names(liver_db$marker_db)
  
  for(cl in cluster_levels) {
    cells <- WhichCells(seurat_obj, idents = cl)
    if(length(cells) > 200) cells <- sample(cells, 200)
    
    for(celltype in names(liver_db$marker_db)) {
      markers <- liver_db$marker_db[[celltype]]
      markers_present <- intersect(markers, rownames(expr_data))
      
      if(length(markers_present) > 0) {
        avg_exp <- mean(expr_data[markers_present, cells])
        score_matrix[cl, celltype] <- avg_exp
      }
    }
  }
  
  # 为每个cluster分配细胞类型
  assignments <- apply(score_matrix, 1, function(x) {
    max_val <- max(x)
    if(max_val > 0.1) {
      return(names(which.max(x)))
    } else {
      return("Unknown")
    }
  })
  
  cat("\n手动注释结果:\n")
  print(table(assignments))
  
  # 应用注释
  seurat_obj$celltype <- assignments[as.character(cluster_ids)]
  seurat_obj$celltype[is.na(seurat_obj$celltype)] <- "Unknown"
  Idents(seurat_obj) <- seurat_obj$celltype
  
  # 保存细胞类型比例
  celltype_counts <- table(Idents(seurat_obj))
  prop_df <- data.frame(
    CellType = names(celltype_counts),
    Count = as.numeric(celltype_counts),
    Percentage = as.numeric(prop.table(celltype_counts) * 100)
  )
  write.csv(prop_df, paste0(sample_name, "_celltype_proportions_fixed.csv"), row.names = FALSE)
  cat("✓ 细胞类型比例已保存\n")
  
  # 绘制细胞类型UMAP
  pdf(paste0(sample_name, "_celltypes_fixed.pdf"), width = 12, height = 8)
  print(DimPlot(seurat_obj, reduction = "umap", label = TRUE, 
                label.size = 4, repel = TRUE) + 
          ggtitle(paste0(sample_name, " Cell Types (", ncol(seurat_obj), " cells)")))
  dev.off()
  cat("✓ 细胞类型图已保存\n")
  
  # ==================== 步骤8: 保存最终对象 ====================
  cat("\n========== 步骤8: 保存最终对象 ==========\n")
  
  saveRDS(seurat_obj, paste0(sample_name, "_complete_fixed.rds"))
  cat("✓ 最终对象已保存\n")
  
  # ==================== 步骤9: 生成报告 ====================
  cat("\n========== 步骤9: 生成报告 ==========\n")
  
  report_file <- paste0(sample_name, "_fix_report.txt")
  sink(report_file)
  
  cat("========================================\n")
  cat("肝脏样本修复报告\n")
  cat("========================================\n\n")
  cat("样本:", sample_name, "\n")
  cat("修复时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("修复后细胞数:", ncol(seurat_obj), "\n")
  cat("修复后基因数:", nrow(seurat_obj), "\n")
  cat("找到的cluster数:", length(unique(Idents(seurat_obj))), "\n\n")
  
  cat("========================================\n")
  cat("细胞类型分布\n")
  cat("========================================\n")
  for(i in 1:nrow(prop_df)) {
    cat(sprintf("  %-15s: %6d (%5.1f%%)\n", 
                prop_df$CellType[i], 
                prop_df$Count[i], 
                prop_df$Percentage[i]))
  }
  
  cat("\n========================================\n")
  cat("关键肝脏marker基因存在性\n")
  cat("========================================\n")
  for(marker in c("ALB", "CD68", "PECAM1", "PDGFRA", "CD3D", "CD19")) {
    exists <- marker %in% rownames(seurat_obj)
    cat(sprintf("  %-15s: %s\n", marker, ifelse(exists, "✓ 存在", "✗ 不存在")))
  }
  
  cat("\n========================================\n")
  cat("生成的输出文件\n")
  cat("========================================\n")
  output_files <- list.files(pattern = paste0(sample_name, ".*fixed.*"))
  for(f in output_files) {
    file_size <- file.size(f)
    if(file_size > 1024^2) {
      size_str <- sprintf("%.2f MB", file_size / 1024^2)
    } else {
      size_str <- sprintf("%.2f KB", file_size / 1024)
    }
    cat("  -", f, "(", size_str, ")\n")
  }
  
  sink()
  cat("✓ 报告已保存:", report_file, "\n")
  
} else {
  stop("数据行数与features文件行数不匹配")
}

cat("\n========== 修复完成 ==========\n")
