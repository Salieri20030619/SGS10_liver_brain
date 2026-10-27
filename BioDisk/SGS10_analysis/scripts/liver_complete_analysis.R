# ==================== 肝脏样本完整分析脚本（整合版） ====================
# 文件名: ~/SGS10_analysis/scripts/liver_complete_analysis.R
# 用法: Rscript liver_complete_analysis.R sample_name n_cores

args <- commandArgs(trailingOnly = TRUE)
sample_name <- args[1]
n_cores <- as.numeric(args[2])

# 设置输出立即显示
options(echo = TRUE)

cat("========================================\n")
cat("肝脏样本完整分析:", sample_name, "\n")
cat("使用核心数:", n_cores, "\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================\n\n")

# 1. 加载包
cat("加载必要的R包...\n")
library(Seurat)
library(ggplot2)
library(dplyr)
library(future)
library(SingleR)
library(celldex)
library(patchwork)

# 检查并加载pheatmap
if(!requireNamespace("pheatmap", quietly = TRUE)) {
  cat("安装pheatmap包...\n")
  install.packages("pheatmap", repos = "https://cloud.r-project.org/")
}
library(pheatmap)

# 检查包版本
cat("✓ Seurat版本:", as.character(packageVersion("Seurat")), "\n")
cat("✓ SingleR版本:", as.character(packageVersion("SingleR")), "\n")
cat("✓ celldex版本:", as.character(packageVersion("celldex")), "\n\n")

# 2. 设置并行
cat("设置并行计算...\n")
plan("multicore", workers = n_cores)
options(future.globals.maxSize = 50 * 1024^3)
options(future.rng.onMisuse = "ignore")
set.seed(42)
cat("✓ 并行设置完成，使用", n_cores, "核心\n\n")

# 3. 创建工作目录
sample_dir <- file.path("/mnt/BioDisk/wangmingyuan/SGS10_analysis", sample_name)
if(!dir.exists(sample_dir)) {
  dir.create(sample_dir, recursive = TRUE)
}
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

# ==================== 检查点设置 ====================
checkpoint_file <- paste0(sample_name, "_checkpoint.txt")
completed_steps <- c()
if(file.exists(checkpoint_file)) {
  completed_steps <- readLines(checkpoint_file, warn = FALSE)
  cat("检测到已完成步骤:", paste(completed_steps, collapse=" → "), "\n\n")
}

save_checkpoint <- function(step) {
  write(c(completed_steps, step), checkpoint_file)
  completed_steps <<- c(completed_steps, step)
  cat("✓ 步骤已记录:", step, "\n")
}

check_step <- function(step, force = FALSE) {
  if(force) return(FALSE)
  if(step %in% completed_steps) {
    file_name <- switch(step,
      "data_loaded" = paste0(sample_name, "_raw.rds"),
      "qc_done" = paste0(sample_name, "_qc.rds"),
      "norm_done" = paste0(sample_name, "_norm.rds"),
      "pca_done" = paste0(sample_name, "_pca.rds"),
      "cluster_done" = paste0(sample_name, "_cluster.rds"),
      "umap_done" = paste0(sample_name, "_umap.rds"),
      "singler_done" = paste0(sample_name, "_singler.rds"),
      "manual_annotation_done" = paste0(sample_name, "_manual_annotated.rds")
    )
    if(!is.null(file_name) && file.exists(file_name)) {
      return(TRUE)
    }
  }
  return(FALSE)
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

# ==================== 创建人类-小鼠同源基因映射表 ====================
create_human_mouse_homologs <- function() {
  homolog_map <- list(
    # 肝细胞标记
    "ALB" = "Alb", "AFP" = "Afp", "TTR" = "Ttr", "HNF4A" = "Hnf4a", "HNF1A" = "Hnf1a",
    "CYP2E1" = "Cyp2e1", "CYP1A2" = "Cyp1a2", "CYP3A4" = "Cyp3a11",
    
    # 内皮细胞标记
    "PECAM1" = "Pecam1", "CDH5" = "Cdh5", "VWF" = "Vwf", "KDR" = "Kdr",
    "CLEC4G" = "Clec4g", "LYVE1" = "Lyve1",
    
    # 星状细胞标记
    "PDGFRA" = "Pdgfra", "PDGFRB" = "Pdgfrb", "DES" = "Des", "ACTA2" = "Acta2",
    "COL1A1" = "Col1a1", "LRAT" = "Lrat",
    
    # 枯否细胞/巨噬细胞标记
    "CD68" = "Cd68", "CD163" = "Cd163", "VSIG4" = "Vsig4", "C1QA" = "C1qa",
    "MRC1" = "Mrc1", "TREM2" = "Trem2", "GPNMB" = "Gpnmb",
    
    # T细胞标记
    "CD3D" = "Cd3d", "CD3E" = "Cd3e", "CD4" = "Cd4", "CD8A" = "Cd8a",
    "FOXP3" = "Foxp3",
    
    # B细胞标记
    "CD19" = "Cd19", "MS4A1" = "Ms4a1", "CD79A" = "Cd79a",
    
    # NK细胞标记
    "NKG7" = "Nkg7", "KLRB1" = "Klrb1",
    
    # 单核细胞标记
    "CD14" = "Cd14", "LYZ" = "Lyz", "S100A8" = "S100a8", "S100A9" = "S100a9",
    
    # 胆管细胞标记
    "KRT7" = "Krt7", "KRT8" = "Krt8", "KRT19" = "Krt19", "EPCAM" = "Epcam",
    
    # 看家基因
    "GAPDH" = "Gapdh", "ACTB" = "Actb"
  )
  return(homolog_map)
}

# ==================== 步骤1: 加载数据 ====================
cat("\n========== 步骤1: 加载数据 ==========\n")

if(check_step("data_loaded")) {
  cat("检测到原始数据文件存在，直接加载...\n")
  seurat_obj <- readRDS(paste0(sample_name, "_raw.rds"))
  cat("细胞数:", ncol(seurat_obj), "\n")
  cat("基因数:", nrow(seurat_obj), "\n")
} else {
  data_dir <- data_paths[[sample_name]]
  if(is.null(data_dir)) stop("未知的样本名: ", sample_name)
  
  cat("数据目录:", data_dir, "\n")
  if(!dir.exists(data_dir)) stop("数据目录不存在: ", data_dir)
  
  # 读取features文件进行基因名转换
  features_file <- file.path(data_dir, "features.tsv.gz")
  if(file.exists(features_file)) {
    features <- read.table(features_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    colnames(features) <- c("gene_id", "gene_name", "feature_type")
    cat("Features文件维度:", dim(features), "\n")
    
    clean_gene_name <- function(x) {
      x <- gsub("^GRCh38-", "", x)
      x <- gsub("^GRCh38_", "", x)
      x <- gsub("^mm10-", "", x)
      x <- gsub("^mm10_", "", x)
      x <- gsub("\\.[0-9]+$", "", x)
      return(x)
    }
    
    features$clean_name <- clean_gene_name(features$gene_name)
    homolog_map <- create_human_mouse_homologs()
    
    features$final_name <- features$clean_name
    for(human_gene in names(homolog_map)) {
      idx <- which(features$clean_name == human_gene)
      if(length(idx) > 0) {
        features$final_name[idx] <- homolog_map[[human_gene]]
      }
    }
    
    data <- Read10X(data_dir)
    cat("原始数据维度:", dim(data), "\n")
    
    if(nrow(data) == nrow(features)) {
      rownames(data) <- features$final_name
      
      dup_genes <- unique(rownames(data)[duplicated(rownames(data))])
      if(length(dup_genes) > 0) {
        gene_max_exp <- data.frame(gene = rownames(data), max_exp = apply(data, 1, max))
        keep_idx <- rep(TRUE, nrow(data))
        for(gene in dup_genes) {
          idx <- which(rownames(data) == gene)
          if(length(idx) > 1) {
            max_idx <- idx[which.max(gene_max_exp$max_exp[idx])]
            keep_idx[idx[-which(idx == max_idx)]] <- FALSE
          }
        }
        data <- data[keep_idx, ]
      }
    }
  } else {
    data <- Read10X(data_dir)
  }
  
  cat("最终数据维度:", dim(data), "\n")
  
  seurat_obj <- CreateSeuratObject(counts = data, project = sample_name)
  
  mt_genes <- grep("^mt-", rownames(seurat_obj), value = TRUE, ignore.case = TRUE)
  if(length(mt_genes) > 0) {
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
  } else {
    seurat_obj[["percent.mt"]] <- 0
  }
  
  saveRDS(seurat_obj, paste0(sample_name, "_raw.rds"))
  cat("✓ 原始数据已保存\n")
  save_checkpoint("data_loaded")
}

# ==================== 步骤2: 质控 ====================
cat("\n========== 步骤2: 质控 ==========\n")

if(check_step("qc_done")) {
  cat("检测到质控后数据存在，直接加载...\n")
  seurat_obj <- readRDS(paste0(sample_name, "_qc.rds"))
  cat("质控后细胞数:", ncol(seurat_obj), "\n")
} else {
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
  cat("  线粒体比例范围:", range(seurat_obj$percent.mt), "\n")
  
  pdf(paste0(sample_name, "_QC.pdf"), width = 12, height = 6)
  print(VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                ncol = 3, pt.size = 0.1))
  dev.off()
  cat("✓ 质控图已保存\n")
  
  saveRDS(seurat_obj, paste0(sample_name, "_qc.rds"))
  cat("✓ 质控后数据已保存\n")
  save_checkpoint("qc_done")
}

# ==================== 步骤3: 标准化 ====================
cat("\n========== 步骤3: 标准化 ==========\n")

if(check_step("norm_done")) {
  cat("检测到标准化后数据存在，直接加载...\n")
  seurat_obj <- readRDS(paste0(sample_name, "_norm.rds"))
  cat("高变基因数:", length(VariableFeatures(seurat_obj)), "\n")
} else {
  cat("标准化数据...\n")
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  
  cat("寻找高变基因...\n")
  seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 2000, verbose = FALSE)
  cat("高变基因数:", length(VariableFeatures(seurat_obj)), "\n")
  
  pdf(paste0(sample_name, "_variable_features.pdf"), width = 10, height = 6)
  top10 <- head(VariableFeatures(seurat_obj), 10)
  plot1 <- VariableFeaturePlot(seurat_obj)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  print(plot2)
  dev.off()
  cat("✓ 高变基因图已保存\n")
  
  cat("数据缩放...\n")
  seurat_obj <- ScaleData(seurat_obj, vars.to.regress = "percent.mt", verbose = FALSE)
  
  saveRDS(seurat_obj, paste0(sample_name, "_norm.rds"))
  cat("✓ 标准化后数据已保存\n")
  save_checkpoint("norm_done")
}

# ==================== 步骤4: PCA降维 ====================
cat("\n========== 步骤4: PCA降维 ==========\n")

if(check_step("pca_done")) {
  cat("检测到PCA结果存在，加载文件...\n")
  seurat_obj <- readRDS(paste0(sample_name, "_pca.rds"))
} else {
  seurat_obj <- RunPCA(seurat_obj, npcs = 50, verbose = FALSE)
  
  pdf(paste0(sample_name, "_elbow.pdf"), width = 8, height = 6)
  print(ElbowPlot(seurat_obj, ndims = 50))
  dev.off()
  cat("✓ 肘部图已保存\n")
  
  saveRDS(seurat_obj, paste0(sample_name, "_pca.rds"))
  cat("✓ PCA后数据已保存\n")
  save_checkpoint("pca_done")
}

# ==================== 步骤5: 聚类 ====================
cat("\n========== 步骤5: 聚类 ==========\n")

if(check_step("cluster_done")) {
  cat("检测到聚类结果存在，加载文件...\n")
  seurat_obj <- readRDS(paste0(sample_name, "_cluster.rds"))
  cat("聚类结果:\n")
  print(table(Idents(seurat_obj)))
} else {
  pc_use <- 1:30
  cat("使用前", length(pc_use), "个PC\n")
  
  cat("计算邻居关系...\n")
  seurat_obj <- FindNeighbors(seurat_obj, dims = pc_use, verbose = FALSE)
  
  cat("聚类...\n")
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.8, verbose = FALSE)
  
  cat("找到", length(unique(Idents(seurat_obj))), "个cluster\n")
  print(table(Idents(seurat_obj)))
  
  saveRDS(seurat_obj, paste0(sample_name, "_cluster.rds"))
  cat("✓ 聚类后数据已保存\n")
  save_checkpoint("cluster_done")
}

# ==================== 步骤6: UMAP降维 ====================
cat("\n========== 步骤6: UMAP降维 ==========\n")

if(check_step("umap_done")) {
  cat("检测到UMAP结果存在，加载文件...\n")
  seurat_obj <- readRDS(paste0(sample_name, "_umap.rds"))
} else {
  cat("细胞数:", ncol(seurat_obj), "\n")
  cat("运行UMAP...\n")
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, verbose = FALSE)
  
  pdf(paste0(sample_name, "_clusters.pdf"), width = 10, height = 8)
  print(DimPlot(seurat_obj, reduction = "umap", label = TRUE) + 
          ggtitle(paste0(sample_name, " cell clusters (", ncol(seurat_obj), " cells)")))
  dev.off()
  cat("✓ 聚类图已保存\n")
  
  saveRDS(seurat_obj, paste0(sample_name, "_umap.rds"))
  cat("✓ UMAP后数据已保存\n")
  save_checkpoint("umap_done")
}

# ==================== 步骤7: SingleR自动注释 ====================
cat("\n========== 步骤7: SingleR自动注释 ==========\n")

if(check_step("singler_done")) {
  cat("检测到SingleR注释结果存在，加载文件...\n")
  seurat_obj <- readRDS(paste0(sample_name, "_singler.rds"))
  if("celltype" %in% colnames(seurat_obj@meta.data)) {
    Idents(seurat_obj) <- seurat_obj$celltype
    cat("已将celltype设置为active identity\n")
  }
} else {
  cat("获取表达数据...\n")
  expr_data <- GetAssayData(seurat_obj, layer = "data")
  
  cat("加载参考数据集...\n")
  ref_immgen <- ImmGenData()
  ref_mouse <- MouseRNAseqData()
  
  cat("运行SingleR注释...\n")
  singler_results <- SingleR(test = expr_data,
                             ref = list(ImmGen = ref_immgen, Mouse = ref_mouse),
                             labels = list(ref_immgen$label.main, ref_mouse$label.main),
                             clusters = Idents(seurat_obj),
                             assay.type.test = "logcounts",
                             fine.tune = TRUE)
  
  seurat_obj$singler_labels <- singler_results$labels[Idents(seurat_obj)]
  
  # 映射肝脏细胞类型
  label_mapping <- c(
    "Hepatocyte" = "Hepatocytes", "Hepatocytes" = "Hepatocytes",
    "Kupffer cell" = "Kupffer", "Macrophage" = "Kupffer", "Macrophages" = "Kupffer",
    "Endothelial cell" = "Endothelial", "Endothelial cells" = "Endothelial",
    "Stellate cell" = "Stellate", "B cell" = "B_cells", "T cell" = "T_cells",
    "NK cell" = "NK_cells", "Monocyte" = "Monocytes", "Neutrophil" = "Neutrophils"
  )
  
  seurat_obj$celltype <- seurat_obj$singler_labels
  for(old_name in names(label_mapping)) {
    seurat_obj$celltype[seurat_obj$singler_labels == old_name] <- label_mapping[old_name]
  }
  
  other_cells <- is.na(seurat_obj$celltype)
  if(sum(other_cells) > 0) {
    seurat_obj$celltype[other_cells] <- "Other"
  }
  
  Idents(seurat_obj) <- seurat_obj$celltype
  
  cat("\nSingleR注释细胞类型分布:\n")
  print(table(Idents(seurat_obj)))
  
  prop_df <- data.frame(
    CellType = names(table(Idents(seurat_obj))),
    Count = as.numeric(table(Idents(seurat_obj))),
    Percentage = as.numeric(prop.table(table(Idents(seurat_obj))) * 100)
  )
  write.csv(prop_df, paste0(sample_name, "_singler_proportions.csv"), row.names = FALSE)
  
  pdf(paste0(sample_name, "_singler_celltypes.pdf"), width = 12, height = 8)
  print(DimPlot(seurat_obj, reduction = "umap", label = TRUE, 
                label.size = 4, repel = TRUE) + 
          ggtitle(paste0(sample_name, " SingleR Cell Types")))
  dev.off()
  
  saveRDS(seurat_obj, paste0(sample_name, "_singler.rds"))
  cat("✓ SingleR注释结果已保存\n")
  save_checkpoint("singler_done")
}

# ==================== 步骤8: 手动精细注释 ====================
cat("\n========== 步骤8: 手动精细注释 ==========\n")

if(check_step("manual_annotation_done")) {
  cat("检测到手动注释结果存在，加载文件...\n")
  seurat_obj <- readRDS(paste0(sample_name, "_manual_annotated.rds"))
} else {
  cat("创建肝脏细胞marker数据库...\n")
  liver_markers <- create_liver_marker_db()
  cat("定义了", length(liver_markers), "种细胞类型\n")
  
  # 定义大类
  major_names <- c(
    "Hepatocytes", "Stellate", "Endothelial", "Kupffer", 
    "T_cells", "B_cells", "Monocytes", "Neutrophils", 
    "Cholangiocytes", "Pericytes", "NK_cells",
    "Macrophage_M2", "Macrophage_NAM", "LSEC", "HSC_Activated"
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
  
  # 应用注释
  seurat_obj$final_celltype <- final_assignment
  Idents(seurat_obj) <- seurat_obj$final_celltype
  
  cat("\n最终细胞类型分布:\n")
  celltype_counts <- table(Idents(seurat_obj))
  print(celltype_counts)
  
  # 统计发现的细胞类型
  final_types <- names(celltype_counts)
  cat("\n发现的细胞类型 (", length(final_types), "种):\n")
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
  
  # 保存对象
  seurat_obj$subtype_max_score <- all_max_scores
  seurat_obj$annotation_threshold <- threshold
  
  saveRDS(seurat_obj, paste0(sample_name, "_manual_annotated.rds"))
  cat("✓ 手动注释结果已保存\n")
  save_checkpoint("manual_annotation_done")
}

# ==================== 步骤9: 生成报告 ====================
cat("\n========== 步骤9: 生成报告 ==========\n")

report_file <- paste0(sample_name, "_complete_report.txt")
sink(report_file)

cat("========================================\n")
cat("肝脏样本完整分析报告\n")
cat("========================================\n\n")
cat("样本:", sample_name, "\n")
cat("分析时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("总细胞数:", ncol(seurat_obj), "\n")
cat("总基因数:", nrow(seurat_obj), "\n\n")

cat("========================================\n")
cat("最终细胞类型分布\n")
cat("========================================\n")
celltype_counts <- table(seurat_obj$final_celltype)
for(i in 1:length(celltype_counts)) {
  cat(sprintf("  %-20s: %6d (%5.1f%%)\n", 
              names(celltype_counts)[i], 
              celltype_counts[i], 
              celltype_counts[i]/ncol(seurat_obj)*100))
}

cat("\n========================================\n")
cat("分析步骤完成情况\n")
cat("========================================\n")
cat("✓ 数据加载\n")
cat("✓ 质控\n")
cat("✓ 标准化\n")
cat("✓ PCA降维\n")
cat("✓ 聚类\n")
cat("✓ UMAP降维\n")
cat("✓ SingleR自动注释\n")
cat("✓ 手动精细注释\n")

sink()
cat("✓ 报告已保存:", report_file, "\n")

# ==================== 关闭并行 ====================
plan("sequential")
cat("\n并行计算已关闭\n")

cat("\n========== 肝脏样本分析完成 ==========\n")
cat("完成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

# 列出生成的文件
files <- list.files(pattern = paste0(sample_name, "_.*"))
cat("\n生成的文件:\n")
if(length(files) > 0) {
  pdf_files <- files[grep("\\.pdf$", files)]
  csv_files <- files[grep("\\.csv$", files)]
  rds_files <- files[grep("\\.rds$", files)]
  
  if(length(pdf_files) > 0) {
    cat("\n  PDF文件:\n")
    for(f in pdf_files) {
      file_size <- round(file.size(f) / 1024^2, 2)
      cat("    -", f, sprintf("(%.2f MB)", file_size), "\n")
    }
  }
  
  if(length(csv_files) > 0) {
    cat("\n  CSV文件:\n")
    for(f in csv_files) {
      file_size <- round(file.size(f) / 1024, 2)
      cat("    -", f, sprintf("(%.2f KB)", file_size), "\n")
    }
  }
  
  if(length(rds_files) > 0) {
    cat("\n  RDS文件:\n")
    for(f in rds_files) {
      file_size <- round(file.size(f) / 1024^3, 2)
      cat("    -", f, sprintf("(%.2f GB)", file_size), "\n")
    }
  }
}
