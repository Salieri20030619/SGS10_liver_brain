# ==================== 肝脏样本修复脚本（转换为小鼠同源基因） ====================
# 文件名: ~/SGS10_analysis/scripts/liver_samples_fix_human2mouse.R
# 用法: Rscript liver_samples_fix_human2mouse.R sample_name n_cores

args <- commandArgs(trailingOnly = TRUE)
sample_name <- args[1]
n_cores <- as.numeric(args[2])

cat("========================================\n")
cat("肝脏样本修复脚本（人类→小鼠同源基因转换）:", sample_name, "\n")
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

# ==================== 创建人类-小鼠同源基因映射表 ====================
create_human_mouse_homologs <- function() {
  cat("创建人类-小鼠同源基因映射表...\n")
  
  # 肝脏关键基因的人类-小鼠同源映射（手动整理）
  homolog_map <- list(
    # 肝细胞标记
    "ALB" = "Alb", "AFP" = "Afp", "CYP3A4" = "Cyp3a11", "CYP2E1" = "Cyp2e1",
    "APOA1" = "Apoa1", "TTR" = "Ttr", "HNF4A" = "Hnf4a", "PCK1" = "Pck1",
    "G6PC" = "G6pc", "FOXA1" = "Foxa1", "FOXA2" = "Foxa2",
    
    # 内皮细胞标记
    "PECAM1" = "Pecam1", "CDH5" = "Cdh5", "VWF" = "Vwf", "KDR" = "Kdr",
    "FLT1" = "Flt1", "TEK" = "Tek", "CLEC4G" = "Clec4g", "LYVE1" = "Lyve1",
    
    # 星状细胞标记
    "PDGFRA" = "Pdgfra", "PDGFRB" = "Pdgfrb", "DES" = "Des", "RGS5" = "Rgs5",
    "ACTA2" = "Acta2", "COL1A1" = "Col1a1", "LRAT" = "Lrat", "RELN" = "Reln",
    
    # 枯否细胞/巨噬细胞标记
    "CD68" = "Cd68", "CD163" = "Cd163", "VSIG4" = "Vsig4", "MARCO" = "Marco",
    "C1QA" = "C1qa", "C1QB" = "C1qb", "MRC1" = "Mrc1", "ITGAM" = "Itgam",
    
    # T细胞标记
    "CD3D" = "Cd3d", "CD3E" = "Cd3e", "CD4" = "Cd4", "CD8A" = "Cd8a",
    "CD8B" = "Cd8b", "FOXP3" = "Foxp3",
    
    # B细胞标记
    "CD19" = "Cd19", "MS4A1" = "Ms4a1", "CD79A" = "Cd79a", "PAX5" = "Pax5",
    
    # NK细胞标记
    "NKG7" = "Nkg7", "KLRB1" = "Klrb1", "KLRD1" = "Klrd1", "KLRK1" = "Klrk1",
    "NCR1" = "Ncr1",
    
    # 单核细胞标记
    "CD14" = "Cd14", "FCGR3A" = "Fcgr3", "LYZ" = "Lyz", "S100A8" = "S100a8",
    "S100A9" = "S100a9", "CSF1R" = "Csf1r",
    
    # 胆管细胞标记
    "KRT19" = "Krt19", "EPCAM" = "Epcam", "SOX9" = "Sox9", "KRT7" = "Krt7",
    
    # 看家基因
    "GAPDH" = "Gapdh", "ACTB" = "Actb", "B2M" = "B2m", "HPRT1" = "Hprt"
  )
  
  # 返回映射表
  return(homolog_map)
}

homolog_map <- create_human_mouse_homologs()
cat("创建了", length(homolog_map), "个人类-小鼠同源基因映射\n")
cat("示例映射:\n")
print(head(names(homolog_map), 20))
print(head(unlist(homolog_map), 20))

# ==================== 读取features文件 ====================
cat("\n========== 步骤1: 读取features文件 ==========\n")

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

# ==================== 提取标准基因名 ====================
cat("\n========== 步骤2: 提取标准基因名 ==========\n")

extract_standard_name <- function(gene_name) {
  # 去除GRCh38-前缀
  cleaned <- gsub("^GRCh38-", "", gene_name)
  cleaned <- gsub("^mm10-", "", cleaned)
  # 去除版本号
  cleaned <- gsub("\\..*$", "", cleaned)
  return(cleaned)
}

# 提取所有基因名
gene_features$clean_name <- extract_standard_name(gene_features$gene_name)

cat("标准基因名示例:\n")
print(head(gene_features$clean_name, 30))

# ==================== 转换为小鼠同源基因 ====================
cat("\n========== 步骤3: 转换为小鼠同源基因 ==========\n")

# 创建映射向量
human_to_mouse <- unlist(homolog_map)
names(human_to_mouse) <- names(homolog_map)

# 转换基因名
gene_features$mouse_name <- gene_features$clean_name
for(human_gene in names(human_to_mouse)) {
  idx <- which(gene_features$clean_name == human_gene)
  if(length(idx) > 0) {
    gene_features$mouse_name[idx] <- human_to_mouse[human_gene]
  }
}

# 统计转换情况
converted_count <- sum(gene_features$clean_name != gene_features$mouse_name)
cat("转换的基因数:", converted_count, "\n")
cat("转换示例:\n")
converted_examples <- gene_features[gene_features$clean_name != gene_features$mouse_name, ]
if(nrow(converted_examples) > 0) {
  print(head(converted_examples[, c("clean_name", "mouse_name")], 20))
}

# ==================== 读取10X数据 ====================
cat("\n========== 步骤4: 读取10X数据 ==========\n")

# 读取原始数据
cat("读取10X数据...\n")
data <- Read10X(data_dir)
cat("原始数据维度:", dim(data), "\n")

# 确保数据与features文件顺序一致
if(nrow(data) == nrow(features)) {
  # 应用转换
  data_filtered <- data
  rownames(data_filtered) <- gene_features$mouse_name
  
  # 去重（保留表达量最高的）
  dup_genes <- unique(rownames(data_filtered)[duplicated(rownames(data_filtered))])
  if(length(dup_genes) > 0) {
    cat("处理重复基因:", length(dup_genes), "\n")
    
    # 对于重复基因，保留表达量最高的
    gene_max_exp <- data.frame(gene = rownames(data_filtered), 
                               max_exp = apply(data_filtered, 1, max))
    
    keep_idx <- rep(TRUE, nrow(data_filtered))
    for(gene in dup_genes) {
      idx <- which(rownames(data_filtered) == gene)
      # 只保留表达量最高的
      max_idx <- idx[which.max(gene_max_exp$max_exp[idx])]
      keep_idx[idx[-which(idx == max_idx)]] <- FALSE
    }
    data_filtered <- data_filtered[keep_idx, ]
    cat("去重后基因数:", nrow(data_filtered), "\n")
  }
  
  cat("转换后数据维度:", dim(data_filtered), "\n")
  cat("示例基因名:\n")
  print(head(rownames(data_filtered), 30))
  
  # 检查小鼠肝脏marker
  mouse_liver_markers <- c("Alb", "Cd68", "Pecam1", "Pdgfra", "Cd3d", "Cd19")
  found_markers <- intersect(mouse_liver_markers, rownames(data_filtered))
  cat("\n找到的小鼠肝脏marker:", paste(found_markers, collapse = ", "), "\n")
  
  if(length(found_markers) == 0) {
    cat("\n⚠️ 警告：未找到小鼠肝脏marker，尝试使用人类marker直接搜索...\n")
    human_liver_markers <- c("ALB", "CD68", "PECAM1", "PDGFRA", "CD3D", "CD19")
    found_human <- intersect(human_liver_markers, rownames(data_filtered))
    cat("找到的人类marker:", paste(found_human, collapse = ", "), "\n")
  }
  
  # ==================== 创建Seurat对象 ====================
  cat("\n========== 步骤5: 创建Seurat对象 ==========\n")
  
  # 检查是否有重复行名
  if(any(duplicated(rownames(data_filtered)))) {
    cat("移除重复行名...\n")
    data_filtered <- data_filtered[!duplicated(rownames(data_filtered)), ]
  }
  
  seurat_obj <- CreateSeuratObject(counts = data_filtered, project = sample_name)
  
  # 计算线粒体基因比例（小鼠格式：mt-）
  mt_genes <- grep("^mt-", rownames(seurat_obj), value = TRUE, ignore.case = TRUE)
  if(length(mt_genes) > 0) {
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
  } else {
    seurat_obj[["percent.mt"]] <- 0
  }
  
  cat("✓ 创建成功！细胞数:", ncol(seurat_obj), "\n")
  cat("基因数:", nrow(seurat_obj), "\n")
  
  # 保存原始对象
  saveRDS(seurat_obj, paste0(sample_name, "_raw_mouse_converted.rds"))
  cat("✓ 原始对象已保存\n")
  
  # ==================== 质控 ====================
  cat("\n========== 步骤6: 质控 ==========\n")
  
  cat("质控前统计:\n")
  cat("  细胞数:", ncol(seurat_obj), "\n")
  cat("  基因数范围:", range(seurat_obj$nFeature_RNA), "\n")
  
  # 质控（参考脑样本的阈值）
  seurat_obj <- subset(seurat_obj, 
                       subset = nFeature_RNA > 200 & 
                                nFeature_RNA < 6000 & 
                                percent.mt < 20)
  
  cat("\n质控后统计:\n")
  cat("  细胞数:", ncol(seurat_obj), "\n")
  
  if(ncol(seurat_obj) < 100) {
    stop("质控后细胞数不足100，无法继续分析")
  }
  
  # 绘制质控图
  pdf(paste0(sample_name, "_QC_mouse_converted.pdf"), width = 12, height = 6)
  print(VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                ncol = 3, pt.size = 0.1))
  dev.off()
  cat("✓ 质控图已保存\n")
  
  # ==================== 标准化和降维 ====================
  cat("\n========== 步骤7: 标准化和降维 ==========\n")
  
  # 标准化
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 2000, verbose = FALSE)
  seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
  
  # PCA
  seurat_obj <- RunPCA(seurat_obj, npcs = 50, verbose = FALSE)
  
  # 肘部图
  pdf(paste0(sample_name, "_elbow_mouse_converted.pdf"), width = 8, height = 6)
  print(ElbowPlot(seurat_obj, ndims = 50))
  dev.off()
  
  # 聚类
  pc_use <- 1:min(30, ncol(seurat_obj[["pca"]]@cell.embeddings))
  seurat_obj <- FindNeighbors(seurat_obj, dims = pc_use, verbose = FALSE)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.8, verbose = FALSE)
  cat("找到", length(unique(Idents(seurat_obj))), "个cluster\n")
  
  # UMAP
  seurat_obj <- RunUMAP(seurat_obj, dims = pc_use, verbose = FALSE)
  
  pdf(paste0(sample_name, "_clusters_mouse_converted.pdf"), width = 10, height = 8)
  print(DimPlot(seurat_obj, reduction = "umap", label = TRUE) + 
          ggtitle(paste0(sample_name, " cell clusters (", ncol(seurat_obj), " cells)")))
  dev.off()
  cat("✓ UMAP图已保存\n")
  
  # ==================== SingleR注释（使用小鼠参考数据，参考脑样本） ====================
  cat("\n========== 步骤8: SingleR自动注释（小鼠参考） ==========\n")
  
  library(SingleR)
  library(celldex)
  
  expr_data <- GetAssayData(seurat_obj, layer = "data")
  
  # 使用小鼠参考数据（与脑样本相同）
  cat("加载小鼠参考数据...\n")
  ref_immgen <- ImmGenData()
  ref_mouse <- MouseRNAseqData()
  
  cat("运行SingleR注释...\n")
  singler_results <- SingleR(
    test = expr_data,
    ref = list(ImmGen = ref_immgen, Mouse = ref_mouse),
    labels = list(ref_immgen$label.main, ref_mouse$label.main),
    clusters = Idents(seurat_obj),
    assay.type.test = "logcounts"
  )
  
  seurat_obj$singler_labels <- singler_results$labels[Idents(seurat_obj)]
  
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
    "B cell" = "B_cells",
    "T cell" = "T_cells",
    "NK cell" = "NK_cells",
    "Monocyte" = "Monocytes"
  )
  
  seurat_obj$celltype <- seurat_obj$singler_labels
  for(old_name in names(label_mapping)) {
    seurat_obj$celltype[seurat_obj$singler_labels == old_name] <- label_mapping[old_name]
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
  write.csv(prop_df, paste0(sample_name, "_celltype_proportions_mouse_converted.csv"), row.names = FALSE)
  cat("✓ 细胞类型比例已保存\n")
  
  # 绘制细胞类型UMAP
  pdf(paste0(sample_name, "_celltypes_mouse_converted.pdf"), width = 12, height = 8)
  print(DimPlot(seurat_obj, reduction = "umap", label = TRUE, 
                label.size = 4, repel = TRUE) + 
          ggtitle(paste0(sample_name, " Cell Types (Human→Mouse Converted)")))
  dev.off()
  cat("✓ 细胞类型图已保存\n")
  
  # ==================== 保存最终对象 ====================
  cat("\n========== 步骤9: 保存最终对象 ==========\n")
  
  saveRDS(seurat_obj, paste0(sample_name, "_complete_mouse_converted.rds"))
  cat("✓ 最终对象已保存\n")
  
  # ==================== 生成报告 ====================
  cat("\n========== 步骤10: 生成报告 ==========\n")
  
  report_file <- paste0(sample_name, "_mouse_conversion_report.txt")
  sink(report_file)
  
  cat("========================================\n")
  cat("肝脏样本修复报告（人类→小鼠同源转换）\n")
  cat("========================================\n\n")
  cat("样本:", sample_name, "\n")
  cat("修复时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("转换后细胞数:", ncol(seurat_obj), "\n")
  cat("转换后基因数:", nrow(seurat_obj), "\n")
  cat("找到的cluster数:", length(unique(Idents(seurat_obj))), "\n")
  cat("人类基因转小鼠基因数:", converted_count, "\n\n")
  
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
  cat("关键肝脏marker基因存在性（小鼠版本）\n")
  cat("========================================\n")
  for(marker in c("Alb", "Cd68", "Pecam1", "Pdgfra", "Cd3d", "Cd19")) {
    exists <- marker %in% rownames(seurat_obj)
    cat(sprintf("  %-15s: %s\n", marker, ifelse(exists, "✓ 存在", "✗ 不存在")))
  }
  
  cat("\n========================================\n")
  cat("生成的输出文件\n")
  cat("========================================\n")
  output_files <- list.files(pattern = paste0(sample_name, ".*mouse_converted.*"))
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
