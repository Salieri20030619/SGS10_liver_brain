# ==================== 肝脏样本修复脚本（参考脑样本成功经验） ====================
# 文件名: ~/SGS10_analysis/scripts/liver_samples_fix_like_brain.R
# 用法: Rscript liver_samples_fix_like_brain.R sample_name n_cores

args <- commandArgs(trailingOnly = TRUE)
sample_name <- args[1]
n_cores <- as.numeric(args[2])

cat("========================================\n")
cat("参考脑样本成功经验修复肝脏样本:", sample_name, "\n")
cat("========================================\n\n")

library(Seurat)
library(ggplot2)
library(dplyr)
library(future)
library(SingleR)
library(celldex)

# 设置并行
plan("multicore", workers = n_cores)
options(future.globals.maxSize = 50 * 1024^3)
set.seed(42)

sample_dir <- file.path("/mnt/BioDisk/wangmingyuan/SGS10_analysis", sample_name)
if(!dir.exists(sample_dir)) dir.create(sample_dir, recursive = TRUE)
setwd(sample_dir)

# ==================== 数据路径 ====================
data_paths <- list(
  gan56 = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/gan56/gan_56_filtered_feature_bc_matrix",
  liver3 = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/liver_samples/25090559_Liver_3_expression/Liver_3_filtered_feature_bc_matrix",
  liver8 = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/liver_samples/25090559_Liver_8_expression/Liver_8_filtered_feature_bc_matrix",
  bla_liver = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/liver_samples/25090559_Bla_liver_1_expression/Bla_liver_1_filtered_feature_bc_matrix",
  con_liver = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/liver_samples/25090559_Con_liver_1_expression/Con_liver_1_filtered_feature_bc_matrix"
)

# ==================== 关键修复：从features.tsv获取标准基因名 ====================
cat("========== 步骤1: 从原始数据提取标准基因名 ==========\n")

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

# 使用gene_name作为标准基因名
standard_gene_names <- gene_features$gene_name
# 去除空值
standard_gene_names <- standard_gene_names[standard_gene_names != "" & !is.na(standard_gene_names)]
cat("有效标准基因名数:", length(standard_gene_names), "\n")
cat("标准基因名示例:\n")
print(head(standard_gene_names, 30))

# 检查是否有蛋白编码基因（肝脏marker）
liver_markers_check <- c("ALB", "AFP", "CYP3A4", "APOA1", "PECAM1", "CDH5", 
                         "CD68", "CD163", "PDGFRA", "PDGFRB", "TTR", "HNF4A")
found_markers <- intersect(liver_markers_check, standard_gene_names)
cat("\n找到的肝脏marker:", paste(found_markers, collapse = ", "), "\n")

if(length(found_markers) == 0) {
  cat("\n⚠️ 警告：未找到标准肝脏marker，尝试使用Ensembl ID映射...\n")
  
  # 提取Ensembl ID
  has_ensembl <- grepl("ENSG", features$gene_id)
  if(sum(has_ensembl) > 0) {
    cat("找到", sum(has_ensembl), "个Ensembl ID\n")
    
    # 使用biomaRt映射
    library(biomaRt)
    ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    
    ensembl_ids <- gsub("\\..*$", "", features$gene_id[has_ensembl])
    gene_mapping <- getBM(
      attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype"),
      filters = "ensembl_gene_id",
      values = ensembl_ids,
      mart = ensembl
    )
    
    # 过滤蛋白编码基因
    protein_coding <- gene_mapping[gene_mapping$gene_biotype == "protein_coding", ]
    cat("蛋白编码基因数:", nrow(protein_coding), "\n")
    
    # 创建映射
    mapping_dict <- setNames(protein_coding$external_gene_name, protein_coding$ensembl_gene_id)
    standard_gene_names <- rep(NA, nrow(features))
    standard_gene_names[has_ensembl] <- mapping_dict[ensembl_ids]
    
    # 过滤NA
    valid_idx <- !is.na(standard_gene_names) & standard_gene_names != ""
    standard_gene_names <- standard_gene_names[valid_idx]
    cat("映射后的标准基因名数:", length(standard_gene_names), "\n")
  }
}

# ==================== 步骤2: 读取10X数据并应用标准基因名 ====================
cat("\n========== 步骤2: 读取数据并应用标准基因名 ==========\n")

# 读取原始数据
cat("读取10X数据...\n")
data <- Read10X(data_dir)
cat("原始数据维度:", dim(data), "\n")

# 确保基因名与features文件顺序一致
if(nrow(data) == nrow(features)) {
  # 只保留有标准基因名的基因
  keep_idx <- !is.na(standard_gene_names) & standard_gene_names != ""
  data_filtered <- data[keep_idx, ]
  rownames(data_filtered) <- standard_gene_names[keep_idx]
  
  # 去重（取平均值）
  dup_genes <- unique(rownames(data_filtered)[duplicated(rownames(data_filtered))])
  if(length(dup_genes) > 0) {
    cat("处理重复基因:", length(dup_genes), "\n")
    counts_unique <- matrix(0, nrow = length(unique(rownames(data_filtered))), 
                           ncol = ncol(data_filtered))
    rownames(counts_unique) <- unique(rownames(data_filtered))
    
    for(gene in unique(rownames(data_filtered))) {
      gene_idx <- which(rownames(data_filtered) == gene)
      if(length(gene_idx) > 1) {
        counts_unique[gene, ] <- colMeans(data_filtered[gene_idx, ])
      } else {
        counts_unique[gene, ] <- data_filtered[gene_idx, ]
      }
    }
    data_filtered <- counts_unique
  }
  
  cat("过滤后数据维度:", dim(data_filtered), "\n")
  cat("示例基因名:\n")
  print(head(rownames(data_filtered), 30))
  
  # 检查肝脏marker
  found_markers_after <- intersect(liver_markers_check, rownames(data_filtered))
  cat("\n过滤后找到的肝脏marker:", paste(found_markers_after, collapse = ", "), "\n")
  
  if(length(found_markers_after) == 0) {
    stop("错误：过滤后仍然没有找到肝脏marker，数据可能不是标准的基因表达数据")
  }
  
  # 创建Seurat对象
  cat("\n创建Seurat对象...\n")
  seurat_obj <- CreateSeuratObject(counts = data_filtered, project = sample_name)
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  cat("✓ 创建成功！细胞数:", ncol(seurat_obj), "\n")
  cat("基因数:", nrow(seurat_obj), "\n")
  
  # 保存原始对象
  saveRDS(seurat_obj, paste0(sample_name, "_raw_fixed.rds"))
  cat("✓ 原始对象已保存\n")
} else {
  stop("数据行数与features文件行数不匹配")
}

# ==================== 步骤3: 质控（参考脑样本） ====================
cat("\n========== 步骤3: 质控 ==========\n")

cat("质控前统计:\n")
cat("  细胞数:", ncol(seurat_obj), "\n")
cat("  基因数范围:", range(seurat_obj$nFeature_RNA), "\n")
cat("  线粒体比例范围:", range(seurat_obj$percent.mt), "\n")

seurat_obj <- subset(seurat_obj, 
                     subset = nFeature_RNA > 200 & 
                              nFeature_RNA < 6000 & 
                              percent.mt < 20)

cat("\n质控后统计:\n")
cat("  细胞数:", ncol(seurat_obj), "\n")
cat("  基因数范围:", range(seurat_obj$nFeature_RNA), "\n")

# ==================== 步骤4: 标准化和降维（完全参考脑样本） ====================
cat("\n========== 步骤4: 标准化和降维 ==========\n")

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
pc_use <- 1:30
seurat_obj <- FindNeighbors(seurat_obj, dims = pc_use, verbose = FALSE)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.8, verbose = FALSE)
cat("找到", length(unique(Idents(seurat_obj))), "个cluster\n")

# UMAP
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, verbose = FALSE)

pdf(paste0(sample_name, "_clusters_fixed.pdf"), width = 10, height = 8)
print(DimPlot(seurat_obj, reduction = "umap", label = TRUE) + 
        ggtitle(paste0(sample_name, " cell clusters (", ncol(seurat_obj), " cells)")))
dev.off()

# ==================== 步骤5: SingleR注释（参考脑样本） ====================
cat("\n========== 步骤5: SingleR自动注释 ==========\n")

expr_data <- GetAssayData(seurat_obj, layer = "data")

# 根据样本类型选择参考数据
if(sample_name %in% c("gan56", "liver3", "liver8")) {
  cat("人类样本，使用人类参考数据...\n")
  ref_hpca <- HumanPrimaryCellAtlasData()
  ref_monaco <- MonacoImmuneData()
  ref_blueprint <- BlueprintEncodeData()
  
  # 使用多个参考（类似脑样本的多参考策略）
  singler_results <- SingleR(
    test = expr_data,
    ref = list(HPCA = ref_hpca, Monaco = ref_monaco, Blueprint = ref_blueprint),
    labels = list(ref_hpca$label.main, ref_monaco$label.main, ref_blueprint$label.main),
    clusters = Idents(seurat_obj),
    assay.type.test = "logcounts",
    fine.tune = TRUE
  )
  
  # 投票机制
  vote_matrix <- do.call(cbind, lapply(singler_results$labels, function(x) x))
  final_labels <- apply(vote_matrix, 1, function(x) {
    tbl <- table(x)
    names(which.max(tbl))
  })
  
  seurat_obj$celltype <- final_labels[Idents(seurat_obj)]
  
} else {
  cat("小鼠样本，使用小鼠参考数据...\n")
  ref_immgen <- ImmGenData()
  ref_mouse <- MouseRNAseqData()
  
  singler_results <- SingleR(
    test = expr_data,
    ref = list(ImmGen = ref_immgen, Mouse = ref_mouse),
    labels = list(ref_immgen$label.main, ref_mouse$label.main),
    clusters = Idents(seurat_obj),
    assay.type.test = "logcounts"
  )
  
  seurat_obj$celltype <- singler_results$labels[Idents(seurat_obj)]
}

# 映射到肝脏细胞类型
liver_mapping <- c(
  "Hepatocyte" = "Hepatocytes",
  "Hepatocytes" = "Hepatocytes",
  "Kupffer_cell" = "Kupffer",
  "Macrophage" = "Kupffer",
  "Macrophages" = "Kupffer",
  "Endothelial_cell" = "Endothelial",
  "Endothelial cells" = "Endothelial",
  "Stellate_cell" = "Stellate"
)

for(old_name in names(liver_mapping)) {
  seurat_obj$celltype[seurat_obj$celltype == old_name] <- liver_mapping[old_name]
}

seurat_obj$celltype[is.na(seurat_obj$celltype)] <- "Other"
Idents(seurat_obj) <- seurat_obj$celltype

cat("\n最终细胞类型分布:\n")
print(table(Idents(seurat_obj)))

# 保存细胞类型比例
prop_df <- data.frame(
  CellType = names(table(Idents(seurat_obj))),
  Count = as.numeric(table(Idents(seurat_obj))),
  Percentage = as.numeric(prop.table(table(Idents(seurat_obj))) * 100)
)
write.csv(prop_df, paste0(sample_name, "_celltype_proportions_fixed.csv"), row.names = FALSE)

# 绘制细胞类型UMAP
pdf(paste0(sample_name, "_celltypes_fixed.pdf"), width = 12, height = 8)
print(DimPlot(seurat_obj, reduction = "umap", label = TRUE, 
              label.size = 4, repel = TRUE) + 
        ggtitle(paste0(sample_name, " Cell Types (Fixed)")))
dev.off()

# ==================== 步骤6: 保存最终对象 ====================
saveRDS(seurat_obj, paste0(sample_name, "_complete_fixed.rds"))
cat("\n✓ 最终对象已保存\n")

# ==================== 步骤7: 生成报告 ====================
report_file <- paste0(sample_name, "_fix_report.txt")
sink(report_file)

cat("========================================\n")
cat("肝脏样本修复报告\n")
cat("========================================\n\n")
cat("样本:", sample_name, "\n")
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

# 列出关键肝脏marker
cat("\n========================================\n")
cat("关键肝脏marker基因存在性\n")
cat("========================================\n")
for(marker in liver_markers_check) {
  exists <- marker %in% rownames(seurat_obj)
  cat(sprintf("  %-15s: %s\n", marker, ifelse(exists, "✓", "✗")))
}

sink()

cat("✓ 报告已保存:", report_file, "\n")
cat("\n========== 修复完成 ==========\n")
