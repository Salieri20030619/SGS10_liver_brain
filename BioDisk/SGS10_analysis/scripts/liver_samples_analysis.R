# ==================== 肝脏样本统一分析脚本（人鼠混合数据修复版） ====================
# 文件名: ~/SGS10_analysis/scripts/liver_samples_analysis.R
# 用法: Rscript liver_samples_analysis.R sample_name n_cores [run_fine_annotation]
# 修复：将人鼠混合数据统一转换为人类标准格式，使用人类参考数据进行注释

args <- commandArgs(trailingOnly = TRUE)
sample_name <- args[1]
n_cores <- as.numeric(args[2])
run_fine_annotation <- ifelse(length(args) >= 3, as.logical(args[3]), TRUE)

# 设置输出立即显示
options(echo = TRUE)

cat("========================================\n")
cat("开始处理肝脏样本:", sample_name, "\n")
cat("使用核心数:", n_cores, "\n")
cat("精细注释:", ifelse(run_fine_annotation, "是", "否"), "\n")
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

# 检查包版本
cat("✓ Seurat版本:", as.character(packageVersion("Seurat")), "\n")
cat("✓ SingleR版本:", as.character(packageVersion("SingleR")), "\n")
cat("✓ celldex版本:", as.character(packageVersion("celldex")), "\n\n")

# 2. 设置并行
cat("设置并行计算...\n")
plan("multicore", workers = n_cores)
options(future.globals.maxSize = 30 * 1024^3)  # 30GB内存限制
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
  if(force) {
    cat("强制运行步骤:", step, "(忽略检查点)\n")
    return(FALSE)
  }
  
  if(step %in% completed_steps) {
    file_name <- switch(step,
      "data_loaded" = paste0(sample_name, "_raw.rds"),
      "qc_done" = paste0(sample_name, "_qc.rds"),
      "norm_done" = paste0(sample_name, "_norm.rds"),
      "pca_done" = paste0(sample_name, "_pca.rds"),
      "cluster_done" = paste0(sample_name, "_cluster.rds"),
      "umap_done" = paste0(sample_name, "_umap.rds"),
      "singler_done" = paste0(sample_name, "_singler.rds"),
      "marker_done" = paste0(sample_name, "_all_markers.csv"),
      "kupffer_done" = paste0(sample_name, "_kupffer_subset.rds"),
      "hepatocyte_done" = paste0(sample_name, "_hepatocyte_subset.rds"),
      "stellate_done" = paste0(sample_name, "_stellate_subset.rds")
    )
    if(!is.null(file_name) && file.exists(file_name)) {
      return(TRUE)
    }
  }
  return(FALSE)
}

# ==================== 判断样本类型 ====================
# 人类样本：gan56, liver3, liver8
is_human_sample <- function(sample) {
  human_samples <- c("gan56", "liver3", "liver8")
  return(sample %in% human_samples)
}

# ==================== 基因名清理函数 ====================
# 将人鼠混合数据的基因名统一转换为人类标准格式
convert_to_human_genes <- function(expr_data) {
  gene_names <- rownames(expr_data)
  original_count <- length(gene_names)
  
  cat("  原始基因名示例:\n")
  print(head(gene_names, 20))
  
  # 1. 移除 GRCh38- 前缀（人类基因）
  cleaned <- gsub("^GRCh38-", "", gene_names)
  
  # 2. 移除 mm10- 前缀（小鼠基因）并添加 "_mouse" 后缀避免冲突
  cleaned <- gsub("^mm10-", "", cleaned)
  
  # 3. 检查并处理重复
  dup_count <- sum(duplicated(cleaned))
  if(dup_count > 0) {
    cat("  发现", dup_count, "个重复基因名，正在处理...\n")
    # 为小鼠基因添加 "_mouse" 后缀
    mouse_genes <- grep("^mm10-", gene_names, value = TRUE)
    mouse_cleaned <- gsub("^mm10-", "", mouse_genes)
    
    # 重新处理：人类基因保持原名，小鼠基因加后缀
    cleaned <- gene_names
    cleaned <- gsub("^GRCh38-", "", cleaned)
    for(mg in mouse_genes) {
      idx <- which(gene_names == mg)
      cleaned[idx] <- paste0(gsub("^mm10-", "", mg), "_mouse")
    }
  }
  
  # 4. 去除空字符串
  valid_idx <- cleaned != ""
  expr_data <- expr_data[valid_idx, ]
  cleaned <- cleaned[valid_idx]
  
  # 5. 去重
  keep_idx <- !duplicated(cleaned)
  expr_data <- expr_data[keep_idx, ]
  cleaned <- cleaned[keep_idx]
  
  rownames(expr_data) <- cleaned
  
  cat("  清理后基因名示例:\n")
  print(head(rownames(expr_data), 20))
  cat("  基因数变化:", original_count, "→", nrow(expr_data), "\n")
  
  return(expr_data)
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
  
  cat("读取10X数据...\n")
  data <- Read10X(data_dir)
  cat("数据维度:", nrow(data), "基因 ×", ncol(data), "细胞\n")
  
  cat("创建Seurat对象...\n")
  seurat_obj <- CreateSeuratObject(counts = data, project = sample_name)
  cat("✓ 创建成功！细胞数:", ncol(seurat_obj), "\n")
  
  # 计算线粒体基因比例（使用统一的 mt- 模式）
  cat("计算线粒体基因比例...\n")
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "mt-")
  
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
  cat("  UMI数范围:", range(seurat_obj$nCount_RNA), "\n")
  cat("  线粒体比例范围:", range(seurat_obj$percent.mt), "\n")
  
  # 统一阈值
  seurat_obj <- subset(seurat_obj, 
                       subset = nFeature_RNA > 200 & 
                                nFeature_RNA < 6000 & 
                                percent.mt < 20)
  
  cat("\n质控后统计:\n")
  cat("  细胞数:", ncol(seurat_obj), "\n")
  cat("  基因数范围:", range(seurat_obj$nFeature_RNA), "\n")
  cat("  UMI数范围:", range(seurat_obj$nCount_RNA), "\n")
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
  
  pdf(paste0(sample_name, "_pca_dimplot.pdf"), width = 8, height = 6)
  print(DimPlot(seurat_obj, reduction = "pca") + ggtitle(paste0(sample_name, " PCA")))
  dev.off()
  cat("✓ PCA图已保存\n")
  
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
  cat("表达矩阵维度:", dim(expr_data), "\n")
  
  # 将基因名转换为人类标准格式
  cat("\n转换基因名为人类标准格式...\n")
  expr_data <- convert_to_human_genes(expr_data)
  
  # 根据样本类型选择参考数据
  if(is_human_sample(sample_name)) {
    cat("\n人类样本，使用人类参考数据 HumanPrimaryCellAtlasData...\n")
    ref <- HumanPrimaryCellAtlasData()
  } else {
    cat("\n小鼠样本，使用小鼠参考数据 MouseRNAseqData...\n")
    ref <- MouseRNAseqData()
  }
  
  cat("参考数据基因名示例:\n")
  print(head(rownames(ref), 20))
  
  # 检查共同基因
  common_genes <- intersect(rownames(expr_data), rownames(ref))
  cat("\n共同基因数:", length(common_genes), "\n")
  
  if(length(common_genes) < 100) {
    cat("⚠️ 警告：共同基因数较少，结果可能不准确\n")
  }
  
  cat("\n运行SingleR注释（这可能需要几分钟）...\n")
  
  singler_results <- SingleR(
    test = expr_data,
    ref = ref,
    labels = ref$label.main,
    clusters = Idents(seurat_obj),
    assay.type.test = "logcounts"
  )
  
  seurat_obj$singler_labels <- singler_results$labels[Idents(seurat_obj)]
  
  # 肝脏细胞类型映射
  label_mapping <- c(
    "Hepatocyte" = "Hepatocytes",
    "Hepatocytes" = "Hepatocytes",
    "Kupffer_cell" = "Kupffer",
    "Macrophage" = "Kupffer",
    "Macrophages" = "Kupffer",
    "Endothelial_cell" = "Endothelial",
    "Endothelial cells" = "Endothelial",
    "Stellate_cell" = "Stellate",
    "Cholangiocyte" = "Cholangiocytes",
    "B_cell" = "B_cells",
    "T_cell" = "T_cells",
    "NK_cell" = "NK_cells"
  )
  
  seurat_obj$celltype <- seurat_obj$singler_labels
  for (old_name in names(label_mapping)) {
    seurat_obj$celltype[seurat_obj$singler_labels == old_name] <- label_mapping[old_name]
  }
  
  other_cells <- is.na(seurat_obj$celltype)
  if(sum(other_cells) > 0) {
    seurat_obj$celltype[other_cells] <- "Other"
  }
  
  Idents(seurat_obj) <- seurat_obj$celltype
  
  cat("\n最终细胞类型分布:\n")
  print(table(Idents(seurat_obj)))
  
  prop_df <- data.frame(
    CellType = names(table(Idents(seurat_obj))),
    Count = as.numeric(table(Idents(seurat_obj))),
    Percentage = as.numeric(prop.table(table(Idents(seurat_obj))) * 100)
  )
  write.csv(prop_df, paste0(sample_name, "_celltype_proportions.csv"), row.names = FALSE)
  cat("✓ 细胞类型比例已保存\n")
  
  pdf(paste0(sample_name, "_celltypes.pdf"), width = 12, height = 8)
  print(DimPlot(seurat_obj, reduction = "umap", label = TRUE, 
                label.size = 4, repel = TRUE) + 
          ggtitle(paste0(sample_name, " Cell Types (", ncol(seurat_obj), " cells)")))
  dev.off()
  cat("✓ 细胞类型图已保存\n")
  
  saveRDS(seurat_obj, paste0(sample_name, "_singler.rds"))
  cat("✓ SingleR注释结果已保存\n")
  save_checkpoint("singler_done")
}

# ==================== 步骤8: marker分析 ====================

cat("\n========== 步骤8: marker分析 ==========\n")

if(check_step("marker_done")) {
  cat("检测到marker结果存在，跳过计算\n")
  all_markers <- read.csv(paste0(sample_name, "_all_markers.csv"))
  cat("已存在的marker基因数:", nrow(all_markers), "\n")
} else {
  cat("使用Seurat计算marker基因...\n")
  
  Idents(seurat_obj) <- seurat_obj$celltype
  
  start_time <- Sys.time()
  all_markers <- FindAllMarkers(seurat_obj, 
                                only.pos = TRUE,
                                min.pct = 0.25,
                                logfc.threshold = 0.25,
                                max.cells.per.ident = 200,
                                verbose = FALSE,
                                future.seed = TRUE)
  end_time <- Sys.time()
  
  cat("Seurat marker分析完成，用时:", round(difftime(end_time, start_time, units = "secs"), 2), "秒\n")
  cat("找到", nrow(all_markers), "个marker基因\n")
  
  if(nrow(all_markers) > 0) {
    write.csv(all_markers, paste0(sample_name, "_all_markers.csv"), row.names = FALSE)
    cat("✓ marker结果已保存\n")
    
    top_markers <- all_markers %>%
      group_by(cluster) %>%
      top_n(n = 10, wt = avg_log2FC)
    
    write.csv(top_markers, paste0(sample_name, "_top10_markers.csv"), row.names = FALSE)
    cat("✓ top10 marker已保存\n")
    
    top_features <- top_markers %>%
      group_by(cluster) %>%
      slice_head(n = 5) %>%
      pull(gene) %>%
      unique()
    
    if(length(top_features) > 0 && length(top_features) <= 100) {
      pdf(paste0(sample_name, "_marker_heatmap.pdf"), width = 14, height = 10)
      print(DoHeatmap(seurat_obj, features = top_features) + 
              scale_fill_gradientn(colors = c("blue", "white", "red")) +
              ggtitle(paste0(sample_name, " Top Marker Genes")))
      dev.off()
      cat("✓ marker热图已保存\n")
    }
  }
  save_checkpoint("marker_done")
}

# ==================== 步骤9: 枯否细胞亚群分析 ====================

if(run_fine_annotation) {
  cat("\n========== 步骤9: 枯否细胞(Kupffer)亚群分析 ==========\n")
  
  if(!check_step("kupffer_done", force = TRUE)) {
    kupffer_names <- c("Kupffer", "Macrophage", "Macrophages")
    kupffer_present <- kupffer_names[kupffer_names %in% unique(seurat_obj$celltype)]
    
    if(length(kupffer_present) > 0) {
      kupffer_name <- kupffer_present[1]
      kupffer_cells <- WhichCells(seurat_obj, expression = celltype == kupffer_name)
      
      cat("找到枯否细胞类型:", kupffer_name, "\n")
      cat("枯否细胞数量:", length(kupffer_cells), "\n")
      
      if(length(kupffer_cells) >= 100) {
        kupffer_obj <- subset(seurat_obj, cells = kupffer_cells)
        
        cat("\n重新处理枯否细胞...\n")
        kupffer_obj <- NormalizeData(kupffer_obj, verbose = FALSE)
        kupffer_obj <- FindVariableFeatures(kupffer_obj, nfeatures = 800, verbose = FALSE)
        kupffer_obj <- ScaleData(kupffer_obj, verbose = FALSE)
        kupffer_obj <- RunPCA(kupffer_obj, npcs = 15, verbose = FALSE)
        
        pdf(paste0(sample_name, "_kupffer_elbow.pdf"), width = 8, height = 6)
        print(ElbowPlot(kupffer_obj, ndims = 15))
        dev.off()
        cat("  肘部图已保存\n")
        
        kupffer_obj <- FindNeighbors(kupffer_obj, dims = 1:10, verbose = FALSE)
        kupffer_obj <- FindClusters(kupffer_obj, resolution = 0.2, verbose = FALSE)
        kupffer_obj <- RunUMAP(kupffer_obj, dims = 1:10, verbose = FALSE)
        
        pdf(paste0(sample_name, "_kupffer_clusters.pdf"), width = 10, height = 8)
        print(DimPlot(kupffer_obj, reduction = "umap", label = TRUE) + 
                ggtitle(paste0(sample_name, " Kupffer Subclusters")))
        dev.off()
        cat("  亚群UMAP图已保存\n")
        
        kupffer_counts <- table(Idents(kupffer_obj))
        cat("\n枯否细胞亚群分布:\n")
        print(kupffer_counts)
        
        # 枯否细胞功能状态marker
        kupffer_markers <- list(
          "Homeostatic" = c("CLEC4F", "VSIG4", "TIMD4", "CD163"),
          "M1_activated" = c("IL1B", "NOS2", "TNF", "CD86"),
          "M2_activated" = c("CD163", "MRC1", "IL10", "ARG1"),
          "Lipid_associated" = c("LPL", "CD36", "PPARG"),
          "Proliferating" = c("MKI67", "TOP2A", "PCNA")
        )
        
        # 计算每个cluster的功能状态分数
        Idents(kupffer_obj) <- kupffer_obj$seurat_clusters
        cluster_func <- c()
        
        for(i in levels(Idents(kupffer_obj))) {
          cells <- WhichCells(kupffer_obj, idents = i)
          if(length(cells) > 100) cells <- sample(cells, 100)
          
          avg_exp <- GetAssayData(kupffer_obj, layer = "data")[, cells] %>% rowMeans()
          
          scores <- sapply(names(kupffer_markers), function(st) {
            markers <- kupffer_markers[[st]]
            markers_present <- intersect(markers, names(avg_exp))
            if(length(markers_present) > 0) {
              return(mean(avg_exp[markers_present]))
            } else {
              return(0)
            }
          })
          cluster_func[i] <- names(scores)[which.max(scores)]
        }
        
        kupffer_obj$kupffer_state <- Idents(kupffer_obj)
        levels(kupffer_obj$kupffer_state) <- cluster_func[levels(kupffer_obj$kupffer_state)]
        
        pdf(paste0(sample_name, "_kupffer_states.pdf"), width = 12, height = 8)
        print(DimPlot(kupffer_obj, reduction = "umap", group.by = "kupffer_state", 
                      label = TRUE, label.size = 3, repel = TRUE) +
                ggtitle(paste0(sample_name, " Kupffer Cell States")))
        dev.off()
        cat("  枯否细胞状态图已保存\n")
        
        saveRDS(kupffer_obj, paste0(sample_name, "_kupffer_subset.rds"))
        cat("✓ 枯否细胞亚群对象已保存\n")
        save_checkpoint("kupffer_done")
      } else {
        cat("枯否细胞数量不足，跳过亚群分析\n")
      }
    } else {
      cat("未找到枯否细胞，跳过亚群分析\n")
    }
  }
}

# ==================== 步骤10: 肝细胞亚群精细注释 ====================

if(run_fine_annotation) {
  cat("\n========== 步骤10: 肝细胞(Hepatocytes)亚群精细注释 ==========\n")
  
  if(!check_step("hepatocyte_done", force = TRUE)) {
    if("Hepatocytes" %in% unique(seurat_obj$celltype)) {
      hep_cells <- WhichCells(seurat_obj, expression = celltype == "Hepatocytes")
      cat("肝细胞数量:", length(hep_cells), "\n")
      
      if(length(hep_cells) >= 500) {
        hep_obj <- subset(seurat_obj, cells = hep_cells)
        
        cat("\n重新处理肝细胞亚群...\n")
        hep_obj <- NormalizeData(hep_obj, verbose = FALSE)
        hep_obj <- FindVariableFeatures(hep_obj, nfeatures = 1000, verbose = FALSE)
        hep_obj <- ScaleData(hep_obj, verbose = FALSE)
        hep_obj <- RunPCA(hep_obj, npcs = 20, verbose = FALSE)
        
        pdf(paste0(sample_name, "_hepatocyte_elbow.pdf"), width = 8, height = 6)
        print(ElbowPlot(hep_obj, ndims = 20))
        dev.off()
        cat("  肘部图已保存\n")
        
        hep_obj <- FindNeighbors(hep_obj, dims = 1:15, verbose = FALSE)
        hep_obj <- FindClusters(hep_obj, resolution = 0.3, verbose = FALSE)
        hep_obj <- RunUMAP(hep_obj, dims = 1:15, verbose = FALSE)
        
        # 肝细胞区域/功能marker（人类）
        hep_markers <- list(
          "Periportal" = c("GLUL", "OAT", "ASS1"),
          "Pericentral" = c("CYP2E1", "GLS2", "CYP1A2"),
          "Metabolic_high" = c("ALB", "TTR", "APOA1"),
          "Proliferating" = c("MKI67", "TOP2A", "PCNA")
        )
        
        Idents(hep_obj) <- hep_obj$seurat_clusters
        cluster_annotations <- c()
        
        for(i in levels(Idents(hep_obj))) {
          cells <- WhichCells(hep_obj, idents = i)
          if(length(cells) > 100) cells <- sample(cells, 100)
          
          avg_exp <- GetAssayData(hep_obj, layer = "data")[, cells] %>% rowMeans()
          
          scores <- sapply(names(hep_markers), function(st) {
            markers <- hep_markers[[st]]
            markers_present <- intersect(markers, names(avg_exp))
            if(length(markers_present) > 0) {
              return(mean(avg_exp[markers_present]))
            } else {
              return(0)
            }
          })
          cluster_annotations[i] <- names(scores)[which.max(scores)]
        }
        
        hep_obj$hepatocyte_subtype <- Idents(hep_obj)
        levels(hep_obj$hepatocyte_subtype) <- cluster_annotations[levels(hep_obj$hepatocyte_subtype)]
        
        pdf(paste0(sample_name, "_hepatocyte_subtypes.pdf"), width = 12, height = 8)
        print(DimPlot(hep_obj, reduction = "umap", group.by = "hepatocyte_subtype", 
                      label = TRUE, label.size = 3, repel = TRUE) +
                ggtitle(paste0(sample_name, " Hepatocyte Subtypes")))
        dev.off()
        cat("  肝细胞亚型UMAP图已保存\n")
        
        saveRDS(hep_obj, paste0(sample_name, "_hepatocyte_subset.rds"))
        cat("✓ 肝细胞亚群对象已保存\n")
        save_checkpoint("hepatocyte_done")
      } else {
        cat("肝细胞数量不足，跳过精细注释\n")
      }
    } else {
      cat("未找到肝细胞，跳过精细注释\n")
    }
  }
}

# ==================== 步骤11: 星状细胞亚群精细注释 ====================

if(run_fine_annotation) {
  cat("\n========== 步骤11: 星状细胞(Stellate)亚群精细注释 ==========\n")
  
  if(!check_step("stellate_done", force = TRUE)) {
    stellate_names <- c("Stellate", "Hepatic stellate cells")
    stellate_present <- stellate_names[stellate_names %in% unique(seurat_obj$celltype)]
    
    if(length(stellate_present) > 0) {
      stellate_name <- stellate_present[1]
      stellate_cells <- WhichCells(seurat_obj, expression = celltype == stellate_name)
      
      cat("星状细胞数量:", length(stellate_cells), "\n")
      
      if(length(stellate_cells) >= 100) {
        stellate_obj <- subset(seurat_obj, cells = stellate_cells)
        
        cat("\n重新处理星状细胞亚群...\n")
        stellate_obj <- NormalizeData(stellate_obj, verbose = FALSE)
        stellate_obj <- FindVariableFeatures(stellate_obj, nfeatures = 800, verbose = FALSE)
        stellate_obj <- ScaleData(stellate_obj, verbose = FALSE)
        stellate_obj <- RunPCA(stellate_obj, npcs = 15, verbose = FALSE)
        
        pdf(paste0(sample_name, "_stellate_elbow.pdf"), width = 8, height = 6)
        print(ElbowPlot(stellate_obj, ndims = 15))
        dev.off()
        cat("  肘部图已保存\n")
        
        stellate_obj <- FindNeighbors(stellate_obj, dims = 1:10, verbose = FALSE)
        stellate_obj <- FindClusters(stellate_obj, resolution = 0.2, verbose = FALSE)
        stellate_obj <- RunUMAP(stellate_obj, dims = 1:10, verbose = FALSE)
        
        # 星状细胞功能状态marker（人类）
        stellate_markers <- list(
          "Quiescent" = c("LRAT", "RELN", "PPARG"),
          "Activated" = c("ACTA2", "DES", "PDGFRB", "COL1A1"),
          "Proliferating" = c("MKI67", "TOP2A", "PCNA")
        )
        
        Idents(stellate_obj) <- stellate_obj$seurat_clusters
        cluster_annotations <- c()
        
        for(i in levels(Idents(stellate_obj))) {
          cells <- WhichCells(stellate_obj, idents = i)
          if(length(cells) > 100) cells <- sample(cells, 100)
          
          avg_exp <- GetAssayData(stellate_obj, layer = "data")[, cells] %>% rowMeans()
          
          scores <- sapply(names(stellate_markers), function(st) {
            markers <- stellate_markers[[st]]
            markers_present <- intersect(markers, names(avg_exp))
            if(length(markers_present) > 0) {
              return(mean(avg_exp[markers_present]))
            } else {
              return(0)
            }
          })
          cluster_annotations[i] <- names(scores)[which.max(scores)]
        }
        
        stellate_obj$stellate_state <- Idents(stellate_obj)
        levels(stellate_obj$stellate_state) <- cluster_annotations[levels(stellate_obj$stellate_state)]
        
        pdf(paste0(sample_name, "_stellate_states.pdf"), width = 12, height = 8)
        print(DimPlot(stellate_obj, reduction = "umap", group.by = "stellate_state", 
                      label = TRUE, label.size = 3, repel = TRUE) +
                ggtitle(paste0(sample_name, " Stellate Cell States")))
        dev.off()
        cat("  星状细胞状态图已保存\n")
        
        saveRDS(stellate_obj, paste0(sample_name, "_stellate_subset.rds"))
        cat("✓ 星状细胞亚群对象已保存\n")
        save_checkpoint("stellate_done")
      } else {
        cat("星状细胞数量不足，跳过精细注释\n")
      }
    } else {
      cat("未找到星状细胞，跳过精细注释\n")
    }
  }
}

# ==================== 步骤12: 保存完整对象 ====================

cat("\n========== 步骤12: 保存完整结果 ==========\n")

if(!file.exists(paste0(sample_name, "_complete.rds"))) {
  saveRDS(seurat_obj, paste0(sample_name, "_complete.rds"))
  file_size <- round(file.size(paste0(sample_name, "_complete.rds")) / 1024^3, 2)
  cat("✓ 完整Seurat对象已保存，大小:", file_size, "GB\n")
}

# ==================== 关闭并行 ====================

plan("sequential")
cat("\n并行计算已关闭\n")

# ==================== 完成 ====================

cat("\n========== ", sample_name, "分析完成 ==========\n")
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
} else {
  cat("  暂无生成文件\n")
}
