# ==================== 肝脏样本完整分析脚本（最终版） ====================
# 文件名: ~/SGS10_analysis/scripts/liver_samples_analysis_final.R
# 用法: Rscript liver_samples_analysis_final.R sample_name n_cores

args <- commandArgs(trailingOnly = TRUE)
sample_name <- args[1]
n_cores <- as.numeric(args[2])

cat("========================================\n")
cat("肝脏样本完整分析（最终版）:", sample_name, "\n")
cat("========================================\n\n")

library(Seurat)
library(ggplot2)
library(dplyr)
library(future)
library(Matrix)
library(SingleR)
library(celldex)

# 设置并行
plan("multicore", workers = n_cores)
options(future.globals.maxSize = 50 * 1024^3)
options(future.rng.onMisuse = "ignore")
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
  
  homolog_map <- list(
    # 肝细胞标记
    "ALB" = "Alb", "AFP" = "Afp", "CYP3A4" = "Cyp3a11", "CYP2E1" = "Cyp2e1",
    "CYP1A2" = "Cyp1a2", "APOA1" = "Apoa1", "APOA2" = "Apoa2", "APOB" = "Apob",
    "TTR" = "Ttr", "HNF4A" = "Hnf4a", "HNF1A" = "Hnf1a", "PCK1" = "Pck1",
    "G6PC" = "G6pc", "FOXA1" = "Foxa1", "FOXA2" = "Foxa2", "FGG" = "Fgg",
    "FGA" = "Fga", "FGB" = "Fgb",
    
    # 内皮细胞标记
    "PECAM1" = "Pecam1", "CDH5" = "Cdh5", "VWF" = "Vwf", "KDR" = "Kdr",
    "FLT1" = "Flt1", "TEK" = "Tek", "ENG" = "Eng", "CLEC4G" = "Clec4g",
    "LYVE1" = "Lyve1", "ICAM1" = "Icam1", "VCAM1" = "Vcam1", "SELE" = "Sele",
    
    # 星状细胞标记
    "PDGFRA" = "Pdgfra", "PDGFRB" = "Pdgfrb", "DES" = "Des", "RGS5" = "Rgs5",
    "ACTA2" = "Acta2", "COL1A1" = "Col1a1", "COL1A2" = "Col1a2", "LRAT" = "Lrat",
    "RELN" = "Reln", "GFAP" = "Gfap",
    
    # 枯否细胞/巨噬细胞标记
    "CD68" = "Cd68", "CD163" = "Cd163", "VSIG4" = "Vsig4", "MARCO" = "Marco",
    "C1QA" = "C1qa", "C1QB" = "C1qb", "C1QC" = "C1qc", "MRC1" = "Mrc1",
    "ITGAM" = "Itgam", "CD86" = "Cd86",
    
    # T细胞标记
    "CD3D" = "Cd3d", "CD3E" = "Cd3e", "CD3G" = "Cd3g", "CD4" = "Cd4",
    "CD8A" = "Cd8a", "CD8B" = "Cd8b", "FOXP3" = "Foxp3", "CD28" = "Cd28",
    "CTLA4" = "Ctla4",
    
    # B细胞标记
    "CD19" = "Cd19", "MS4A1" = "Ms4a1", "CD79A" = "Cd79a", "CD79B" = "Cd79b",
    "PAX5" = "Pax5", "BLK" = "Blk", "BANK1" = "Bank1",
    
    # NK细胞标记
    "NKG7" = "Nkg7", "KLRB1" = "Klrb1", "KLRD1" = "Klrd1", "KLRK1" = "Klrk1",
    "NCR1" = "Ncr1", "NCR3" = "Ncr3", "FCGR3A" = "Fcgr3",
    
    # 单核细胞标记
    "CD14" = "Cd14", "FCGR3A" = "Fcgr3", "LYZ" = "Lyz", "S100A8" = "S100a8",
    "S100A9" = "S100a9", "CSF1R" = "Csf1r",
    
    # 胆管细胞标记
    "KRT7" = "Krt7", "KRT8" = "Krt8", "KRT18" = "Krt18", "KRT19" = "Krt19",
    "EPCAM" = "Epcam", "SOX9" = "Sox9", "HNF1B" = "Hnf1b",
    
    # 看家基因
    "GAPDH" = "Gapdh", "ACTB" = "Actb", "B2M" = "B2m", "HPRT1" = "Hprt",
    "RPLP0" = "Rplp0", "TUBB" = "Tubb", "TBP" = "Tbp",
    
    # 其他重要基因
    "PTPRC" = "Ptprc", "CD44" = "Cd44", "CD24" = "Cd24", "PROM1" = "Prom1",
    "EGFR" = "Egfr", "VEGFA" = "Vegfa", "TGFB1" = "Tgfb1", "TNF" = "Tnf",
    "IL6" = "Il6", "IL1B" = "Il1b", "CXCL8" = "Cxcl8", "MMP2" = "Mmp2",
    "MMP9" = "Mmp9", "TIMP1" = "Timp1"
  )
  
  return(homolog_map)
}

homolog_map <- create_human_mouse_homologs()
cat("创建了", length(homolog_map), "个人类-小鼠同源基因映射\n")

# ==================== 读取features文件 ====================
cat("\n========== 步骤1: 读取features文件 ==========\n")

data_dir <- data_paths[[sample_name]]
features_file <- file.path(data_dir, "features.tsv.gz")

if(!file.exists(features_file)) {
  stop("features.tsv.gz not found at: ", features_file)
}

features <- read.table(features_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(features) <- c("gene_id", "gene_name", "feature_type")
cat("Features文件维度:", dim(features), "\n")

# ==================== 智能基因名处理 ====================
cat("\n========== 步骤2: 智能基因名处理 ==========\n")

# 清理函数
clean_gene_name <- function(x) {
  x <- gsub("^GRCh38-", "", x)
  x <- gsub("^GRCh38_", "", x)
  x <- gsub("^mm10-", "", x)
  x <- gsub("^mm10_", "", x)
  x <- gsub("\\.[0-9]+$", "", x)
  return(x)
}

features$clean_name <- clean_gene_name(features$gene_name)

# 判断基因类型
features$is_mouse <- grepl("^__", features$clean_name)
features$is_human <- !features$is_mouse & toupper(features$clean_name) %in% names(homolog_map)
features$is_noncoding <- !features$is_mouse & !features$is_human

cat("小鼠基因数（带下划线）:", sum(features$is_mouse), "\n")
cat("人类基因数（可转换）:", sum(features$is_human), "\n")
cat("非编码RNA数:", sum(features$is_noncoding), "\n")

# 转换人类基因名
features$final_name <- features$clean_name
human_to_mouse <- unlist(homolog_map)
names(human_to_mouse) <- names(homolog_map)

for(human_gene in names(human_to_mouse)) {
  idx <- which(features$clean_name == human_gene)
  if(length(idx) > 0) {
    features$final_name[idx] <- human_to_mouse[human_gene]
  }
}

# 小鼠基因：移除开头的下划线
features$final_name[features$is_mouse] <- gsub("^__", "", features$final_name[features$is_mouse])

cat("\n转换后基因名示例（人类基因）:\n")
print(head(features[features$is_human, c("clean_name", "final_name")], 20))

cat("\n转换后基因名示例（小鼠基因）:\n")
print(head(features[features$is_mouse, c("clean_name", "final_name")], 20))

# ==================== 读取10X数据 ====================
cat("\n========== 步骤3: 读取10X数据 ==========\n")

data <- Read10X(data_dir)
cat("原始数据维度:", dim(data), "\n")

if(nrow(data) == nrow(features)) {
  # 应用最终基因名
  rownames(data) <- features$final_name
  
  # 过滤空基因名
  valid_idx <- rownames(data) != "" & !is.na(rownames(data))
  data <- data[valid_idx, ]
  cat("过滤空基因名后维度:", dim(data), "\n")
  
  # 去重
  dup_genes <- unique(rownames(data)[duplicated(rownames(data))])
  if(length(dup_genes) > 0) {
    cat("处理重复基因:", length(dup_genes), "\n")
    gene_max_exp <- data.frame(gene = rownames(data), 
                               max_exp = apply(data, 1, max))
    keep_idx <- rep(TRUE, nrow(data))
    for(gene in dup_genes) {
      idx <- which(rownames(data) == gene)
      if(length(idx) > 1) {
        max_idx <- idx[which.max(gene_max_exp$max_exp[idx])]
        keep_idx[idx[-which(idx == max_idx)]] <- FALSE
      }
    }
    data <- data[keep_idx, ]
    cat("去重后基因数:", nrow(data), "\n")
  }
  
  cat("最终数据维度:", dim(data), "\n")
  cat("示例基因名:\n")
  print(head(rownames(data), 50))
  
  # 检查关键marker
  key_markers <- c("Alb", "Cd68", "Pecam1", "Pdgfra", "Cd3d", "Cd19", "Gapdh", "Actb")
  found_markers <- intersect(key_markers, rownames(data))
  cat("\n找到的关键marker:", paste(found_markers, collapse = ", "), "\n")
  
  # ==================== 创建Seurat对象 ====================
  cat("\n========== 步骤4: 创建Seurat对象 ==========\n")
  
  seurat_obj <- CreateSeuratObject(counts = data, project = sample_name)
  
  # 计算线粒体基因比例
  mt_genes <- grep("^mt-", rownames(seurat_obj), value = TRUE, ignore.case = TRUE)
  if(length(mt_genes) > 0) {
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
  } else {
    seurat_obj[["percent.mt"]] <- 0
  }
  
  cat("✓ 创建成功！细胞数:", ncol(seurat_obj), "\n")
  cat("基因数:", nrow(seurat_obj), "\n")
  
  # ==================== 质控 ====================
  cat("\n========== 步骤5: 质控 ==========\n")
  
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
  
  if(ncol(seurat_obj) < 100) {
    stop("质控后细胞数不足100")
  }
  
  # 质控图
  pdf(paste0(sample_name, "_QC.pdf"), width = 12, height = 6)
  print(VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                ncol = 3, pt.size = 0.1))
  dev.off()
  cat("✓ 质控图已保存\n")
  
  # ==================== 标准化和降维 ====================
  cat("\n========== 步骤6: 标准化和降维 ==========\n")
  
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 2000, verbose = FALSE)
  seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
  seurat_obj <- RunPCA(seurat_obj, npcs = 50, verbose = FALSE)
  
  # 肘部图
  pdf(paste0(sample_name, "_elbow.pdf"), width = 8, height = 6)
  print(ElbowPlot(seurat_obj, ndims = 50))
  dev.off()
  
  # 聚类
  pc_use <- 1:min(30, ncol(seurat_obj[["pca"]]@cell.embeddings))
  seurat_obj <- FindNeighbors(seurat_obj, dims = pc_use, verbose = FALSE)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.8, verbose = FALSE)
  cat("找到", length(unique(Idents(seurat_obj))), "个cluster\n")
  
  # UMAP
  seurat_obj <- RunUMAP(seurat_obj, dims = pc_use, verbose = FALSE)
  
  pdf(paste0(sample_name, "_clusters.pdf"), width = 10, height = 8)
  print(DimPlot(seurat_obj, reduction = "umap", label = TRUE) + 
          ggtitle(paste0(sample_name, " cell clusters (", ncol(seurat_obj), " cells)")))
  dev.off()
  cat("✓ UMAP图已保存\n")
  
  # ==================== SingleR注释（小鼠参考） ====================
  cat("\n========== 步骤7: SingleR自动注释 ==========\n")
  
  expr_data <- GetAssayData(seurat_obj, layer = "data")
  
  # 加载小鼠参考数据
  cat("加载小鼠参考数据...\n")
  ref_immgen <- ImmGenData()
  ref_mouse <- MouseRNAseqData()
  
  cat("运行SingleR注释...\n")
  singler_results <- SingleR(
    test = expr_data,
    ref = list(ImmGen = ref_immgen, Mouse = ref_mouse),
    labels = list(ref_immgen$label.main, ref_mouse$label.main),
    clusters = Idents(seurat_obj),
    assay.type.test = "logcounts",
    fine.tune = TRUE
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
  write.csv(prop_df, paste0(sample_name, "_celltype_proportions.csv"), row.names = FALSE)
  cat("✓ 细胞类型比例已保存\n")
  
  # 细胞类型UMAP
  pdf(paste0(sample_name, "_celltypes.pdf"), width = 12, height = 8)
  print(DimPlot(seurat_obj, reduction = "umap", label = TRUE, 
                label.size = 4, repel = TRUE) + 
          ggtitle(paste0(sample_name, " Cell Types (", ncol(seurat_obj), " cells)")))
  dev.off()
  cat("✓ 细胞类型图已保存\n")
  
  # ==================== 保存最终对象 ====================
  cat("\n========== 步骤8: 保存最终对象 ==========\n")
  
  saveRDS(seurat_obj, paste0(sample_name, "_complete.rds"))
  cat("✓ 最终对象已保存\n")
  
  # ==================== 生成报告 ====================
  cat("\n========== 步骤9: 生成报告 ==========\n")
  
  report_file <- paste0(sample_name, "_analysis_report.txt")
  sink(report_file)
  
  cat("========================================\n")
  cat("肝脏样本分析报告\n")
  cat("========================================\n\n")
  cat("样本:", sample_name, "\n")
  cat("分析时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("总细胞数:", ncol(seurat_obj), "\n")
  cat("总基因数:", nrow(seurat_obj), "\n")
  cat("小鼠基因数:", sum(features$is_mouse), "\n")
  cat("转换的人类基因数:", sum(features$is_human), "\n")
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
  cat("关键marker基因存在性\n")
  cat("========================================\n")
  for(marker in c("Alb", "Cd68", "Pecam1", "Pdgfra", "Cd3d", "Cd19", "Gapdh")) {
    exists <- marker %in% rownames(seurat_obj)
    cat(sprintf("  %-15s: %s\n", marker, ifelse(exists, "✓ 存在", "✗ 不存在")))
  }
  
  sink()
  cat("✓ 报告已保存:", report_file, "\n")
  
} else {
  stop("数据行数与features文件行数不匹配")
}

# 关闭并行
plan("sequential")

cat("\n========== 分析完成 ==========\n")
cat("完成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
