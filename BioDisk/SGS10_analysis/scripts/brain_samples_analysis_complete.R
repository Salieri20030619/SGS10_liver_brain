# ==================== 脑样本统一分析脚本 ====================
# 文件名: ~/SGS10_analysis/scripts/brain_samples_analysis_complete.R
# 用法: Rscript brain_samples_analysis_complete.R sample_name n_cores [run_fine_annotation]
# 包含完整的亚群分析 + 神经元和星形胶质细胞精细注释

args <- commandArgs(trailingOnly = TRUE)
sample_name <- args[1]
n_cores <- as.numeric(args[2])
run_fine_annotation <- ifelse(length(args) >= 3, as.logical(args[3]), TRUE)

# 设置输出立即显示
options(echo = TRUE)

cat("========================================\n")
cat("开始处理脑样本:", sample_name, "\n")
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
  nao56 = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/nao_56_filtered_feature_bc_matrix",
  brain3 = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/brain_samples/Brain_3_filtered_feature_bc_matrix",
  brain8 = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/brain_samples/Brain_8_filtered_feature_bc_matrix",
  bla_brain = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/brain_samples/Bla_brain_1_filtered_feature_bc_matrix",
  con_brain = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/brain_samples/Con_brain_1_filtered_feature_bc_matrix"
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

check_step <- function(step) {
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
      "microglia_done" = paste0(sample_name, "_microglia_subset.rds"),
      "neuron_done" = paste0(sample_name, "_neuron_subset.rds"),
      "astro_done" = paste0(sample_name, "_astro_subset.rds")
    )
    if(!is.null(file_name) && file.exists(file_name)) {
      return(TRUE)
    }
  }
  return(FALSE)
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
  
  cat("计算线粒体基因比例...\n")
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
  
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
  
  cat("加载参考数据集...\n")
  ref_immgen <- celldex::ImmGenData()
  ref_mouse <- celldex::MouseRNAseqData()
  
  cat("运行SingleR注释...\n")
  singler_results <- SingleR(test = expr_data,
                             ref = list(ImmGen = ref_immgen, Mouse = ref_mouse),
                             labels = list(ref_immgen$label.main, ref_mouse$label.main),
                             clusters = Idents(seurat_obj),
                             assay.type.test = "logcounts")
  
  seurat_obj$singler_labels <- singler_results$labels[Idents(seurat_obj)]
  
  # 映射细胞类型
  label_mapping <- c(
    "Microglia" = "Microglia",
    "Macrophages" = "Microglia",
    "Monocytes" = "Microglia",
    "Astrocytes" = "Astrocytes",
    "Neurons" = "Neurons",
    "Oligodendrocytes" = "Oligodendrocytes",
    "Oligodendrocytes precursor cells" = "OPCs",
    "OPC" = "OPCs",
    "Endothelial cells" = "Endothelial",
    "Pericytes" = "Pericytes",
    "Ependymal cells" = "Ependymal",
    "Choroid plexus epithelial cells" = "CPECs",
    "Epithelial cells" = "CPECs"
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

# ==================== 步骤9: 小胶质细胞亚群分析 ====================

cat("\n========== 步骤9: 小胶质细胞亚群分析 ==========\n")

if(check_step("microglia_done")) {
  cat("检测到亚群分析结果已存在，跳过\n")
} else {
  microglia_names <- c("Microglia", "Mgl", "Macrophages")
  microglia_present <- microglia_names[microglia_names %in% unique(seurat_obj$celltype)]
  
  if(length(microglia_present) > 0) {
    microglia_name <- microglia_present[1]
    microglia_cells <- WhichCells(seurat_obj, expression = celltype == microglia_name)
    
    cat("找到小胶质细胞类型:", microglia_name, "\n")
    cat("小胶质细胞数量:", length(microglia_cells), "\n")
    cat("占总细胞比例:", round(length(microglia_cells)/ncol(seurat_obj)*100, 2), "%\n")
    
    if(length(microglia_cells) >= 100) {
      microglia_obj <- subset(seurat_obj, cells = microglia_cells)
      
      cat("\n重新处理小胶质细胞...\n")
      microglia_obj <- NormalizeData(microglia_obj, verbose = FALSE)
      microglia_obj <- FindVariableFeatures(microglia_obj, nfeatures = 800, verbose = FALSE)
      microglia_obj <- ScaleData(microglia_obj, verbose = FALSE)
      microglia_obj <- RunPCA(microglia_obj, npcs = 15, verbose = FALSE)
      
      pdf(paste0(sample_name, "_microglia_elbow.pdf"), width = 8, height = 6)
      print(ElbowPlot(microglia_obj, ndims = 15))
      dev.off()
      cat("  肘部图已保存\n")
      
      microglia_obj <- FindNeighbors(microglia_obj, dims = 1:10, verbose = FALSE)
      microglia_obj <- FindClusters(microglia_obj, resolution = 0.2, verbose = FALSE)
      microglia_obj <- RunUMAP(microglia_obj, dims = 1:10, verbose = FALSE)
      
      pdf(paste0(sample_name, "_microglia_clusters.pdf"), width = 10, height = 8)
      print(DimPlot(microglia_obj, reduction = "umap", label = TRUE) + 
              ggtitle(paste0(sample_name, " Microglia Subclusters")))
      dev.off()
      cat("  亚群UMAP图已保存\n")
      
      mg_counts <- table(Idents(microglia_obj))
      cat("\n小胶质细胞亚群分布:\n")
      print(mg_counts)
      
      prop_df <- data.frame(
        Subcluster = paste0("Cluster ", names(mg_counts)),
        Count = as.numeric(mg_counts),
        Percentage = as.numeric(prop.table(mg_counts) * 100)
      )
      write.csv(prop_df, paste0(sample_name, "_microglia_proportions.csv"), row.names = FALSE)
      cat("  亚群比例已保存\n")
      
      mg_markers <- FindAllMarkers(microglia_obj, 
                                   only.pos = TRUE,
                                   min.pct = 0.25,
                                   logfc.threshold = 0.25,
                                   max.cells.per.ident = 200,
                                   verbose = FALSE,
                                   future.seed = TRUE)
      
      if(nrow(mg_markers) > 0) {
        write.csv(mg_markers, paste0(sample_name, "_microglia_markers.csv"), row.names = FALSE)
        cat("  亚群marker已保存\n")
        
        top_mg <- mg_markers %>%
          group_by(cluster) %>%
          top_n(n = 10, wt = avg_log2FC)
        
        top_features_mg <- top_mg %>%
          group_by(cluster) %>%
          slice_head(n = 5) %>%
          pull(gene) %>%
          unique()
        
        if(length(top_features_mg) > 0 && length(top_features_mg) <= 100) {
          pdf(paste0(sample_name, "_microglia_heatmap.pdf"), width = 12, height = 8)
          print(DoHeatmap(microglia_obj, features = top_features_mg) + 
                  scale_fill_gradientn(colors = c("blue", "white", "red")) +
                  ggtitle("Microglia Subcluster Markers"))
          dev.off()
          cat("  亚群热图已保存\n")
        }
      }
      
      saveRDS(microglia_obj, paste0(sample_name, "_microglia_subset.rds"))
      cat("✓ 小胶质细胞亚群对象已保存\n")
      save_checkpoint("microglia_done")
    } else {
      cat("小胶质细胞数量不足，跳过亚群分析\n")
    }
  } else {
    cat("未找到小胶质细胞，跳过亚群分析\n")
  }
}

# ==================== 步骤10: 神经元亚群精细注释 ====================

if(run_fine_annotation) {
  cat("\n========== 步骤10: 神经元亚群精细注释 ==========\n")
  
  if(check_step("neuron_done")) {
    cat("检测到神经元亚群结果存在，跳过\n")
  } else {
    if("Neurons" %in% unique(seurat_obj$celltype)) {
      neuron_cells <- WhichCells(seurat_obj, expression = celltype == "Neurons")
      cat("神经元数量:", length(neuron_cells), "\n")
      
      if(length(neuron_cells) >= 500) {
        neuron_obj <- subset(seurat_obj, cells = neuron_cells)
        
        cat("\n重新处理神经元亚群...\n")
        neuron_obj <- NormalizeData(neuron_obj, verbose = FALSE)
        neuron_obj <- FindVariableFeatures(neuron_obj, nfeatures = 1000, verbose = FALSE)
        neuron_obj <- ScaleData(neuron_obj, verbose = FALSE)
        neuron_obj <- RunPCA(neuron_obj, npcs = 20, verbose = FALSE)
        
        # 肘部图
        pdf(paste0(sample_name, "_neuron_elbow.pdf"), width = 8, height = 6)
        print(ElbowPlot(neuron_obj, ndims = 20))
        dev.off()
        cat("  肘部图已保存\n")
        
        # 聚类（稍高分辨率以发现更多亚型）
        neuron_obj <- FindNeighbors(neuron_obj, dims = 1:15, verbose = FALSE)
        neuron_obj <- FindClusters(neuron_obj, resolution = 0.4, verbose = FALSE)
        neuron_obj <- RunUMAP(neuron_obj, dims = 1:15, verbose = FALSE)
        
        # 神经元亚型marker
        neuron_markers <- list(
          # 兴奋性神经元 - 层特异性
          "Excitatory_L2_3" = c("Cux1", "Cux2", "Rorb", "Satb2"),
          "Excitatory_L4" = c("Rorb", "Kcnk2", "Ntng1"),
          "Excitatory_L5" = c("Fezf2", "Bcl11b", "Etv1", "Tle4"),
          "Excitatory_L6" = c("Tle4", "Foxp2", "Sulf1", "Ctip2"),
          
          # 抑制性神经元
          "Inhibitory_PV" = c("Pvalb", "Gad1", "Gad2"),
          "Inhibitory_SST" = c("Sst", "Gad1", "Gad2"),
          "Inhibitory_VIP" = c("Vip", "Gad1", "Gad2"),
          "Inhibitory_ID2" = c("Id2", "Gad1"),
          
          # 区域特异性
          "Hippocampal" = c("Prox1", "Zbtb20", "Neurod2"),
          "Cortical" = c("Satb2", "Cux1", "Cux2"),
          "Striatal" = c("Isl1", "Foxp1", "Foxp2")
        )
        
        # 计算每个cluster的神经元亚型分数
        Idents(neuron_obj) <- neuron_obj$seurat_clusters
        cluster_annotations <- c()
        
        for(i in levels(Idents(neuron_obj))) {
          cells <- WhichCells(neuron_obj, idents = i)
          if(length(cells) > 100) cells <- sample(cells, 100)
          
          avg_exp <- GetAssayData(neuron_obj, layer = "data")[, cells] %>% rowMeans()
          
          scores <- sapply(names(neuron_markers), function(st) {
            markers <- neuron_markers[[st]]
            markers_present <- intersect(markers, names(avg_exp))
            if(length(markers_present) > 0) {
              return(mean(avg_exp[markers_present]))
            } else {
              return(0)
            }
          })
          cluster_annotations[i] <- names(scores)[which.max(scores)]
        }
        
        neuron_obj$neuron_subtype <- Idents(neuron_obj)
        levels(neuron_obj$neuron_subtype) <- cluster_annotations[levels(neuron_obj$neuron_subtype)]
        
        # 可视化神经元亚型
        pdf(paste0(sample_name, "_neuron_subtypes_umap.pdf"), width = 12, height = 8)
        print(DimPlot(neuron_obj, reduction = "umap", group.by = "neuron_subtype", 
                      label = TRUE, label.size = 3, repel = TRUE) +
                ggtitle(paste0(sample_name, " Neuron Subtypes (", ncol(neuron_obj), " cells)")))
        dev.off()
        cat("  神经元亚型UMAP图已保存\n")
        
        # 神经元亚型比例
        neuron_prop <- prop.table(table(neuron_obj$neuron_subtype)) * 100
        neuron_prop_df <- data.frame(
          Subtype = names(neuron_prop),
          Percentage = as.numeric(neuron_prop)
        )
        write.csv(neuron_prop_df, paste0(sample_name, "_neuron_proportions.csv"), row.names = FALSE)
        
        # 如果有空间信息，绘制空间分布图
        if("spatial" %in% names(seurat_obj@images)) {
          pdf(paste0(sample_name, "_neuron_subtypes_spatial.pdf"), width = 12, height = 10)
          print(SpatialDimPlot(neuron_obj, group.by = "neuron_subtype", label = FALSE) +
                  ggtitle("Neuron Subtypes Spatial Distribution"))
          dev.off()
          cat("  神经元亚型空间分布图已保存\n")
        }
        
        # 保存神经元亚群对象
        saveRDS(neuron_obj, paste0(sample_name, "_neuron_subset.rds"))
        cat("✓ 神经元亚群对象已保存\n")
        save_checkpoint("neuron_done")
      } else {
        cat("神经元数量不足，跳过精细注释\n")
      }
    } else {
      cat("未找到神经元，跳过精细注释\n")
    }
  }
}

# ==================== 步骤11: 星形胶质细胞亚群精细注释 ====================

if(run_fine_annotation) {
  cat("\n========== 步骤11: 星形胶质细胞亚群精细注释 ==========\n")
  
  if(check_step("astro_done")) {
    cat("检测到星形胶质细胞亚群结果存在，跳过\n")
  } else {
    if("Astrocytes" %in% unique(seurat_obj$celltype)) {
      astro_cells <- WhichCells(seurat_obj, expression = celltype == "Astrocytes")
      cat("星形胶质细胞数量:", length(astro_cells), "\n")
      
      if(length(astro_cells) >= 300) {
        astro_obj <- subset(seurat_obj, cells = astro_cells)
        
        cat("\n重新处理星形胶质细胞亚群...\n")
        astro_obj <- NormalizeData(astro_obj, verbose = FALSE)
        astro_obj <- FindVariableFeatures(astro_obj, nfeatures = 800, verbose = FALSE)
        astro_obj <- ScaleData(astro_obj, verbose = FALSE)
        astro_obj <- RunPCA(astro_obj, npcs = 15, verbose = FALSE)
        
        pdf(paste0(sample_name, "_astro_elbow.pdf"), width = 8, height = 6)
        print(ElbowPlot(astro_obj, ndims = 15))
        dev.off()
        cat("  肘部图已保存\n")
        
        astro_obj <- FindNeighbors(astro_obj, dims = 1:10, verbose = FALSE)
        astro_obj <- FindClusters(astro_obj, resolution = 0.3, verbose = FALSE)
        astro_obj <- RunUMAP(astro_obj, dims = 1:10, verbose = FALSE)
        
        # 星形胶质细胞亚型marker
        astro_markers <- list(
          # 区域特异性
          "Astro_Cortical" = c("S100b", "Gfap", "Aqp4", "Slc1a3"),
          "Astro_Hippocampal" = c("S100b", "Gfap", "Aqp4", "Aldh1l1"),
          "Astro_Striatal" = c("S100b", "Gfap", "Aqp4", "Crym"),
          "Astro_Cerebellar" = c("S100b", "Gfap", "Aqp4", "Gdf10"),
          
          # 功能状态
          "Astro_Homeostatic" = c("Slc1a3", "Slc1a2", "Glul", "Aldh1l1"),
          "Astro_Reactive_A1" = c("Gfap", "Serping1", "C3", "H2-D1"),
          "Astro_Reactive_A2" = c("Gfap", "Clcf1", "Tgm1", "S100a10"),
          
          # 发育阶段
          "Astro_Developing" = c("Fabp7", "Glast", "Blbp"),
          "Astro_Mature" = c("Gfap", "Aqp4", "Aldoc")
        )
        
        # 计算每个cluster的星形胶质细胞亚型分数
        Idents(astro_obj) <- astro_obj$seurat_clusters
        cluster_annotations <- c()
        
        for(i in levels(Idents(astro_obj))) {
          cells <- WhichCells(astro_obj, idents = i)
          if(length(cells) > 100) cells <- sample(cells, 100)
          
          avg_exp <- GetAssayData(astro_obj, layer = "data")[, cells] %>% rowMeans()
          
          scores <- sapply(names(astro_markers), function(st) {
            markers <- astro_markers[[st]]
            markers_present <- intersect(markers, names(avg_exp))
            if(length(markers_present) > 0) {
              return(mean(avg_exp[markers_present]))
            } else {
              return(0)
            }
          })
          cluster_annotations[i] <- names(scores)[which.max(scores)]
        }
        
        astro_obj$astro_subtype <- Idents(astro_obj)
        levels(astro_obj$astro_subtype) <- cluster_annotations[levels(astro_obj$astro_subtype)]
        
        # 可视化星形胶质细胞亚型
        pdf(paste0(sample_name, "_astro_subtypes_umap.pdf"), width = 12, height = 8)
        print(DimPlot(astro_obj, reduction = "umap", group.by = "astro_subtype", 
                      label = TRUE, label.size = 3, repel = TRUE) +
                ggtitle(paste0(sample_name, " Astrocyte Subtypes (", ncol(astro_obj), " cells)")))
        dev.off()
        cat("  星形胶质细胞亚型UMAP图已保存\n")
        
        # 星形胶质细胞亚型比例
        astro_prop <- prop.table(table(astro_obj$astro_subtype)) * 100
        astro_prop_df <- data.frame(
          Subtype = names(astro_prop),
          Percentage = as.numeric(astro_prop)
        )
        write.csv(astro_prop_df, paste0(sample_name, "_astro_proportions.csv"), row.names = FALSE)
        
        # 如果有空间信息，绘制空间分布图
        if("spatial" %in% names(seurat_obj@images)) {
          pdf(paste0(sample_name, "_astro_subtypes_spatial.pdf"), width = 12, height = 10)
          print(SpatialDimPlot(astro_obj, group.by = "astro_subtype", label = FALSE) +
                  ggtitle("Astrocyte Subtypes Spatial Distribution"))
          dev.off()
          cat("  星形胶质细胞亚型空间分布图已保存\n")
        }
        
        # 保存星形胶质细胞亚群对象
        saveRDS(astro_obj, paste0(sample_name, "_astro_subset.rds"))
        cat("✓ 星形胶质细胞亚群对象已保存\n")
        save_checkpoint("astro_done")
      } else {
        cat("星形胶质细胞数量不足，跳过精细注释\n")
      }
    } else {
      cat("未找到星形胶质细胞，跳过精细注释\n")
    }
  }
}

# ==================== 步骤12: 空间特征可视化 ====================

if(run_fine_annotation && "spatial" %in% names(seurat_obj@images)) {
  cat("\n========== 步骤12: 空间特征可视化 ==========\n")
  
  # 神经元亚型空间分布
  if(file.exists(paste0(sample_name, "_neuron_subset.rds"))) {
    neuron_obj <- readRDS(paste0(sample_name, "_neuron_subset.rds"))
    
    pdf(paste0(sample_name, "_neuron_subtypes_spatial.pdf"), width = 12, height = 10)
    print(SpatialDimPlot(neuron_obj, group.by = "neuron_subtype", label = FALSE) +
            ggtitle("Neuron Subtypes Spatial Distribution"))
    dev.off()
    cat("✓ 神经元亚型空间分布图已保存\n")
  }
  
  # 星形胶质细胞亚型空间分布
  if(file.exists(paste0(sample_name, "_astro_subset.rds"))) {
    astro_obj <- readRDS(paste0(sample_name, "_astro_subset.rds"))
    
    pdf(paste0(sample_name, "_astro_subtypes_spatial.pdf"), width = 12, height = 10)
    print(SpatialDimPlot(astro_obj, group.by = "astro_subtype", label = FALSE) +
            ggtitle("Astrocyte Subtypes Spatial Distribution"))
    dev.off()
    cat("✓ 星形胶质细胞亚型空间分布图已保存\n")
  }
  
  # 特征基因空间表达
  feature_genes <- c(
    "Cux1", "Cux2", "Rorb", "Fezf2", "Tle4", "Foxp2",
    "Gfap", "Aqp4", "Slc1a3", "C3", "S100a10"
  )
  
  feature_genes <- feature_genes[feature_genes %in% rownames(seurat_obj)]
  if(length(feature_genes) > 0) {
    pdf(paste0(sample_name, "_spatial_features.pdf"), width = 15, height = 12)
    for(i in seq(1, length(feature_genes), by = 4)) {
      genes <- feature_genes[i:min(i+3, length(feature_genes))]
      print(SpatialFeaturePlot(seurat_obj, features = genes, ncol = 2))
    }
    dev.off()
    cat("✓ 空间特征基因图已保存\n")
  }
}

# ==================== 步骤13: 保存完整对象 ====================

cat("\n========== 步骤13: 保存完整结果 ==========\n")

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
