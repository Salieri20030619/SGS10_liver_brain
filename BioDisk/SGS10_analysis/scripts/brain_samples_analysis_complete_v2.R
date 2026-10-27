# ==================== 脑样本完整分析脚本（SingleR + 手动精细注释） ====================
# 文件名: ~/SGS10_analysis/scripts/brain_samples_analysis_complete_v2.R
# 用法: Rscript brain_samples_analysis_complete_v2.R sample_name n_cores [run_fine_annotation]

args <- commandArgs(trailingOnly = TRUE)
sample_name <- args[1]
n_cores <- as.numeric(args[2])
run_fine_annotation <- ifelse(length(args) >= 3, as.logical(args[3]), TRUE)

# 设置输出立即显示
options(echo = TRUE)

cat("========================================\n")
cat("开始处理脑样本（增强版）:", sample_name, "\n")
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
library(patchwork)

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
      "manual_annotation_done" = paste0(sample_name, "_manual_annotated.rds"),
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

# ==================== 脑细胞精细注释marker库 ====================

create_brain_marker_db <- function() {
  # 基于参考图的脑细胞类型marker（小鼠）
  brain_markers <- list(
    # 星形胶质细胞 - 区域特异性
    "Astro_CTX" = c("Slc1a3", "Gfap", "Aqp4", "Aldh1l1", "S100b", "Ctx", "Cortex"),
    "Astro_HPC" = c("Slc1a3", "Gfap", "Aqp4", "Aldh1l1", "S100b", "Hpca", "Hippocampus"),
    "Astro_HYPO" = c("Slc1a3", "Gfap", "Aqp4", "Aldh1l1", "S100b", "Hypothalamus"),
    "Astro_THAL" = c("Slc1a3", "Gfap", "Aqp4", "Aldh1l1", "S100b", "Thalamus"),
    "Astro_WM" = c("Slc1a3", "Gfap", "Aqp4", "Aldh1l1", "S100b", "WhiteMatter", "Wm"),
    "Astro_AMY" = c("Slc1a3", "Gfap", "Aqp4", "Aldh1l1", "S100b", "Amygdala"),
    
    # 兴奋性神经元 - 层/区域特异性
    "Ext_L23" = c("Slc17a7", "Cux1", "Cux2", "Rorb", "Satb2"),
    "Ext_L4" = c("Slc17a7", "Rorb", "Kcnk2", "Ntng1"),
    "Ext_L5" = c("Slc17a7", "Fezf2", "Bcl11b", "Etv1", "Tle4"),
    "Ext_L6" = c("Slc17a7", "Tle4", "Foxp2", "Sulf1"),
    "Ext_Hpc_CA1" = c("Slc17a7", "Neurod6", "Prox1", "Camk2a", "Calb1"),
    "Ext_Hpc_CA3" = c("Slc17a7", "Neurod6", "Prox1", "Camk2a", "Calb2"),
    "Ext_Hpc_DG" = c("Slc17a7", "Prox1", "Dcx", "Calb1", "Pcp4"),
    "Ext_ClauPyr" = c("Slc17a7", "Claustrum", "Piriform", "Cortex"),
    
    # 抑制性神经元 - 分子类型
    "Inh_Sst" = c("Gad1", "Gad2", "Sst", "Cort"),
    "Inh_Pvalb" = c("Gad1", "Gad2", "Pvalb"),
    "Inh_Vip" = c("Gad1", "Gad2", "Vip", "Cck"),
    "Inh_Meis2" = c("Gad1", "Gad2", "Meis2", "Lhx6", "Lhx8"),
    "Inh_Lamp5" = c("Gad1", "Gad2", "Lamp5", "Ndnf"),
    "Inh_Npy" = c("Gad1", "Gad2", "Npy"),
    
    # 小胶质细胞
    "Mgl" = c("Tmem119", "P2ry12", "Cx3cr1", "Cd68", "Hexb", "Csf1r"),
    
    # 少突胶质细胞谱系
    "Oligo" = c("Mbp", "Mog", "Plp1", "Cnp", "Mobp"),
    "OPC" = c("Pdgfra", "Cspg4", "Olig1", "Olig2"),
    "OPCs" = c("Pdgfra", "Cspg4", "Olig1", "Olig2", "Nkx2-1"),
    
    # 内皮细胞
    "ECs" = c("Cldn5", "Pecam1", "Cdh5", "Vwf", "Flt1"),
    
    # 周细胞
    "Pericytes" = c("Rgs5", "Pdgfrb", "Cspg4", "Kcnj8", "Acta2"),
    
    # 室管膜细胞
    "Ependymal" = c("Foxj1", "Vim", "Ccdc153", "Cfap126"),
    
    # 脉络丛上皮细胞
    "CPECs" = c("Ttr", "Lcn2", "Clcn2", "Cldn3", "Slc2a1"),
    
    # 血管软脑膜细胞
    "VLMC" = c("Lum", "Dcn", "Col1a1", "Col3a1", "Pdgfra"),
    
    # 中棘神经元
    "Str_MSN" = c("Drd1", "Drd2", "Adora2a", "Ppp1r1b", "Penk"),
    
    # 技术噪声
    "TECHO" = c("Malat1", "Neat1", "Mtrnr2l2", "Mt-rnr1")
  )
  
  return(brain_markers)
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
                             assay.type.test = "logcounts",
                             fine.tune = TRUE)
  
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

# ==================== 步骤8: 手动精细注释（基于marker） ====================

cat("\n========== 步骤8: 手动精细注释 ==========\n")

if(check_step("manual_annotation_done")) {
  cat("检测到手动手册结果存在，加载文件...\n")
  seurat_obj <- readRDS(paste0(sample_name, "_manual_annotated.rds"))
  cat("手动注释细胞类型分布:\n")
  print(table(seurat_obj$manual_celltype))
} else {
  cat("创建脑细胞marker数据库...\n")
  brain_markers <- create_brain_marker_db()
  cat("定义了", length(brain_markers), "种细胞类型\n")
  
  cat("获取表达数据...\n")
  expr_data <- GetAssayData(seurat_obj, layer = "data")
  
  cat("计算每个细胞类型的marker得分...\n")
  celltype_scores <- matrix(0, nrow = ncol(seurat_obj), ncol = length(brain_markers))
  colnames(celltype_scores) <- names(brain_markers)
  rownames(celltype_scores) <- colnames(seurat_obj)
  
  for(ct in names(brain_markers)) {
    markers <- brain_markers[[ct]]
    markers_present <- intersect(markers, rownames(expr_data))
    if(length(markers_present) > 0) {
      scores <- colMeans(expr_data[markers_present, , drop = FALSE])
      celltype_scores[, ct] <- scores
    } else {
      cat("  ⚠️ 未找到", ct, "的marker，使用部分匹配\n")
      # 尝试部分匹配
      for(m in markers) {
        m_idx <- grep(m, rownames(expr_data), ignore.case = TRUE, value = TRUE)
        if(length(m_idx) > 0) {
          scores <- colMeans(expr_data[m_idx, , drop = FALSE])
          celltype_scores[, ct] <- celltype_scores[, ct] + scores
        }
      }
    }
  }
  
  # 分配细胞类型（基于最高得分）
  cat("分配细胞类型...\n")
  max_scores <- apply(celltype_scores, 1, max)
  primary_assignment <- apply(celltype_scores, 1, function(x) {
    if(max(x) > 0.1) {
      return(names(which.max(x)))
    } else {
      return("Unknown")
    }
  })
  
  seurat_obj$manual_celltype <- primary_assignment
  
  # 基于cluster的共识注释（更稳健）
  cat("基于cluster的共识注释...\n")
  cluster_ids <- Idents(seurat_obj)
  cluster_assignments <- table(cluster_ids, primary_assignment)
  cluster_majority <- apply(cluster_assignments, 1, function(x) {
    names(which.max(x))
  })
  seurat_obj$consensus_celltype <- cluster_majority[as.character(cluster_ids)]
  
  # 合并SingleR和手动注释
  cat("合并SingleR和手动注释...\n")
  seurat_obj$final_celltype <- seurat_obj$consensus_celltype
  # 对于Unknown的细胞，使用SingleR结果
  unknown_idx <- seurat_obj$final_celltype == "Unknown"
  if(sum(unknown_idx) > 0) {
    seurat_obj$final_celltype[unknown_idx] <- seurat_obj$celltype[unknown_idx]
  }
  
  Idents(seurat_obj) <- seurat_obj$final_celltype
  
  cat("\n最终手动注释细胞类型分布:\n")
  print(table(Idents(seurat_obj)))
  
  # 保存细胞类型比例
  prop_df_manual <- data.frame(
    CellType = names(table(Idents(seurat_obj))),
    Count = as.numeric(table(Idents(seurat_obj))),
    Percentage = as.numeric(prop.table(table(Idents(seurat_obj))) * 100)
  )
  write.csv(prop_df_manual, paste0(sample_name, "_manual_proportions.csv"), row.names = FALSE)
  
  # 绘制UMAP
  pdf(paste0(sample_name, "_manual_celltypes.pdf"), width = 14, height = 10)
  p1 <- DimPlot(seurat_obj, reduction = "umap", group.by = "final_celltype", 
                label = TRUE, label.size = 3, repel = TRUE) +
    ggtitle(paste0(sample_name, " Manual Annotation (", ncol(seurat_obj), " cells)"))
  
  p2 <- DimPlot(seurat_obj, reduction = "umap", group.by = "celltype", 
                label = TRUE, label.size = 3, repel = TRUE) +
    ggtitle("SingleR Annotation")
  
  print(p1 | p2)
  dev.off()
  
  # 保存手动注释对象
  saveRDS(seurat_obj, paste0(sample_name, "_manual_annotated.rds"))
  cat("✓ 手动注释结果已保存\n")
  save_checkpoint("manual_annotation_done")
}

# ==================== 步骤9: marker分析 ====================

cat("\n========== 步骤9: marker分析 ==========\n")

if(check_step("marker_done")) {
  cat("检测到marker结果存在，跳过计算\n")
  all_markers <- read.csv(paste0(sample_name, "_all_markers.csv"))
  cat("已存在的marker基因数:", nrow(all_markers), "\n")
} else {
  cat("使用Seurat计算marker基因...\n")
  
  Idents(seurat_obj) <- seurat_obj$final_celltype
  
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

# ==================== 步骤10: 小胶质细胞亚群分析 ====================

if(run_fine_annotation) {
  cat("\n========== 步骤10: 小胶质细胞亚群分析 ==========\n")
  
  if(check_step("microglia_done", force = TRUE)) {
    cat("检测到小胶质细胞亚群结果存在，跳过\n")
  } else {
    microglia_names <- c("Microglia", "Mgl", "Macrophages")
    microglia_present <- microglia_names[microglia_names %in% unique(seurat_obj$final_celltype)]
    
    if(length(microglia_present) > 0) {
      microglia_name <- microglia_present[1]
      microglia_cells <- WhichCells(seurat_obj, expression = final_celltype == microglia_name)
      
      cat("找到小胶质细胞类型:", microglia_name, "\n")
      cat("小胶质细胞数量:", length(microglia_cells), "\n")
      
      if(length(microglia_cells) >= 100) {
        microglia_obj <- subset(seurat_obj, cells = microglia_cells)
        
        cat("\n重新处理小胶质细胞...\n")
        microglia_obj <- NormalizeData(microglia_obj, verbose = FALSE)
        microglia_obj <- FindVariableFeatures(microglia_obj, nfeatures = 800, verbose = FALSE)
        microglia_obj <- ScaleData(microglia_obj, verbose = FALSE)
        microglia_obj <- RunPCA(microglia_obj, npcs = 15, verbose = FALSE)
        microglia_obj <- FindNeighbors(microglia_obj, dims = 1:10, verbose = FALSE)
        microglia_obj <- FindClusters(microglia_obj, resolution = 0.2, verbose = FALSE)
        microglia_obj <- RunUMAP(microglia_obj, dims = 1:10, verbose = FALSE)
        
        # 小胶质细胞功能状态marker
        microglia_state_markers <- list(
          "Homeostatic" = c("P2ry12", "Tmem119", "Cx3cr1", "Siglech"),
          "M1_Activated" = c("Il1b", "Tnfa", "Cd86", "Nos2"),
          "M2_Activated" = c("Arg1", "Cd163", "Mrc1", "Il10"),
          "Proliferating" = c("Mki67", "Top2a", "Pcna"),
          "Disease_Associated" = c("Apoe", "Clec7a", "Lpl", "Cst7")
        )
        
        # 计算每个cluster的状态
        Idents(microglia_obj) <- microglia_obj$seurat_clusters
        cluster_state <- c()
        
        for(i in levels(Idents(microglia_obj))) {
          cells <- WhichCells(microglia_obj, idents = i)
          if(length(cells) > 100) cells <- sample(cells, 100)
          avg_exp <- GetAssayData(microglia_obj, layer = "data")[, cells] %>% rowMeans()
          
          scores <- sapply(names(microglia_state_markers), function(st) {
            markers <- microglia_state_markers[[st]]
            markers_present <- intersect(markers, names(avg_exp))
            if(length(markers_present) > 0) {
              return(mean(avg_exp[markers_present]))
            } else {
              return(0)
            }
          })
          cluster_state[i] <- names(scores)[which.max(scores)]
        }
        
        microglia_obj$microglia_state <- cluster_state[as.character(Idents(microglia_obj))]
        
        pdf(paste0(sample_name, "_microglia_states.pdf"), width = 12, height = 8)
        p1 <- DimPlot(microglia_obj, reduction = "umap", group.by = "microglia_state", 
                      label = TRUE, label.size = 3) +
          ggtitle(paste0(sample_name, " Microglia States"))
        p2 <- DimPlot(microglia_obj, reduction = "umap", label = TRUE) +
          ggtitle("Microglia Clusters")
        print(p1 | p2)
        dev.off()
        
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
}

# ==================== 步骤11: 神经元亚群精细注释 ====================

if(run_fine_annotation) {
  cat("\n========== 步骤11: 神经元亚群精细注释 ==========\n")
  
  if(check_step("neuron_done", force = TRUE)) {
    cat("检测到神经元亚群结果存在，跳过\n")
  } else {
    neuron_names <- c("Neurons", "Ext_L23", "Ext_Hpc_CA1", "Ext_Hpc_CA3", "Ext_Hpc_DG")
    neuron_present <- neuron_names[neuron_names %in% unique(seurat_obj$final_celltype)]
    
    if(length(neuron_present) > 0) {
      neuron_cells <- WhichCells(seurat_obj, expression = final_celltype %in% neuron_present)
      cat("神经元数量:", length(neuron_cells), "\n")
      
      if(length(neuron_cells) >= 500) {
        neuron_obj <- subset(seurat_obj, cells = neuron_cells)
        
        cat("\n重新处理神经元亚群...\n")
        neuron_obj <- NormalizeData(neuron_obj, verbose = FALSE)
        neuron_obj <- FindVariableFeatures(neuron_obj, nfeatures = 1000, verbose = FALSE)
        neuron_obj <- ScaleData(neuron_obj, verbose = FALSE)
        neuron_obj <- RunPCA(neuron_obj, npcs = 20, verbose = FALSE)
        neuron_obj <- FindNeighbors(neuron_obj, dims = 1:15, verbose = FALSE)
        neuron_obj <- FindClusters(neuron_obj, resolution = 0.4, verbose = FALSE)
        neuron_obj <- RunUMAP(neuron_obj, dims = 1:15, verbose = FALSE)
        
        # 神经元亚型marker
        neuron_subtype_markers <- list(
          "Ext_L23" = c("Cux1", "Cux2", "Rorb", "Satb2"),
          "Ext_L4" = c("Rorb", "Kcnk2", "Ntng1"),
          "Ext_L5" = c("Fezf2", "Bcl11b", "Etv1", "Tle4"),
          "Ext_L6" = c("Tle4", "Foxp2", "Sulf1"),
          "Ext_Hpc_CA1" = c("Neurod6", "Prox1", "Camk2a", "Calb1"),
          "Ext_Hpc_CA3" = c("Neurod6", "Prox1", "Camk2a", "Calb2"),
          "Ext_Hpc_DG" = c("Prox1", "Dcx", "Calb1", "Pcp4"),
          "Inh_Sst" = c("Sst", "Gad1", "Gad2"),
          "Inh_Pvalb" = c("Pvalb", "Gad1", "Gad2"),
          "Inh_Vip" = c("Vip", "Gad1", "Gad2"),
          "Inh_Meis2" = c("Meis2", "Lhx6", "Gad1")
        )
        
        # 计算每个cluster的神经元亚型
        Idents(neuron_obj) <- neuron_obj$seurat_clusters
        cluster_annotations <- c()
        
        for(i in levels(Idents(neuron_obj))) {
          cells <- WhichCells(neuron_obj, idents = i)
          if(length(cells) > 100) cells <- sample(cells, 100)
          avg_exp <- GetAssayData(neuron_obj, layer = "data")[, cells] %>% rowMeans()
          
          scores <- sapply(names(neuron_subtype_markers), function(st) {
            markers <- neuron_subtype_markers[[st]]
            markers_present <- intersect(markers, names(avg_exp))
            if(length(markers_present) > 0) {
              return(mean(avg_exp[markers_present]))
            } else {
              return(0)
            }
          })
          cluster_annotations[i] <- names(scores)[which.max(scores)]
        }
        
        neuron_obj$neuron_subtype <- cluster_annotations[as.character(Idents(neuron_obj))]
        
        pdf(paste0(sample_name, "_neuron_subtypes.pdf"), width = 12, height = 8)
        print(DimPlot(neuron_obj, reduction = "umap", group.by = "neuron_subtype", 
                      label = TRUE, label.size = 3, repel = TRUE) +
                ggtitle(paste0(sample_name, " Neuron Subtypes")))
        dev.off()
        
        # 神经元亚型比例
        neuron_prop <- prop.table(table(neuron_obj$neuron_subtype)) * 100
        neuron_df <- data.frame(
          Subtype = names(neuron_prop),
          Count = as.numeric(table(neuron_obj$neuron_subtype)),
          Percentage = as.numeric(neuron_prop)
        )
        write.csv(neuron_df, paste0(sample_name, "_neuron_proportions.csv"), row.names = FALSE)
        
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

# ==================== 步骤12: 星形胶质细胞亚群精细注释 ====================

if(run_fine_annotation) {
  cat("\n========== 步骤12: 星形胶质细胞亚群精细注释 ==========\n")
  
  if(check_step("astro_done", force = TRUE)) {
    cat("检测到星形胶质细胞亚群结果存在，跳过\n")
  } else {
    astro_names <- c("Astrocytes", "Astro_CTX", "Astro_HPC", "Astro_THAL", "Astro_HYPO", "Astro_WM")
    astro_present <- astro_names[astro_names %in% unique(seurat_obj$final_celltype)]
    
    if(length(astro_present) > 0) {
      astro_cells <- WhichCells(seurat_obj, expression = final_celltype %in% astro_present)
      cat("星形胶质细胞数量:", length(astro_cells), "\n")
      
      if(length(astro_cells) >= 300) {
        astro_obj <- subset(seurat_obj, cells = astro_cells)
        
        cat("\n重新处理星形胶质细胞亚群...\n")
        astro_obj <- NormalizeData(astro_obj, verbose = FALSE)
        astro_obj <- FindVariableFeatures(astro_obj, nfeatures = 800, verbose = FALSE)
        astro_obj <- ScaleData(astro_obj, verbose = FALSE)
        astro_obj <- RunPCA(astro_obj, npcs = 15, verbose = FALSE)
        astro_obj <- FindNeighbors(astro_obj, dims = 1:10, verbose = FALSE)
        astro_obj <- FindClusters(astro_obj, resolution = 0.3, verbose = FALSE)
        astro_obj <- RunUMAP(astro_obj, dims = 1:10, verbose = FALSE)
        
        # 星形胶质细胞区域marker
        astro_region_markers <- list(
          "Astro_CTX" = c("Ctx", "Cortex", "Satb2"),
          "Astro_HPC" = c("Hpca", "Hippocampus", "Prox1"),
          "Astro_THAL" = c("Thalamus", "Thal", "Gbx2"),
          "Astro_HYPO" = c("Hypothalamus", "Npy", "Agrp"),
          "Astro_WM" = c("WhiteMatter", "CorpusCallosum", "Mbp")
        )
        
        # 计算每个cluster的区域
        Idents(astro_obj) <- astro_obj$seurat_clusters
        cluster_region <- c()
        
        for(i in levels(Idents(astro_obj))) {
          cells <- WhichCells(astro_obj, idents = i)
          if(length(cells) > 100) cells <- sample(cells, 100)
          avg_exp <- GetAssayData(astro_obj, layer = "data")[, cells] %>% rowMeans()
          
          scores <- sapply(names(astro_region_markers), function(st) {
            markers <- astro_region_markers[[st]]
            markers_present <- intersect(markers, names(avg_exp))
            if(length(markers_present) > 0) {
              return(mean(avg_exp[markers_present]))
            } else {
              return(0)
            }
          })
          cluster_region[i] <- names(scores)[which.max(scores)]
        }
        
        astro_obj$astro_region <- cluster_region[as.character(Idents(astro_obj))]
        
        pdf(paste0(sample_name, "_astro_regions.pdf"), width = 12, height = 8)
        print(DimPlot(astro_obj, reduction = "umap", group.by = "astro_region", 
                      label = TRUE, label.size = 3, repel = TRUE) +
                ggtitle(paste0(sample_name, " Astrocyte Regions")))
        dev.off()
        
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

# ==================== 步骤13: 空间特征可视化 ====================

if(run_fine_annotation && "spatial" %in% names(seurat_obj@images)) {
  cat("\n========== 步骤13: 空间特征可视化 ==========\n")
  
  # 细胞类型空间分布
  pdf(paste0(sample_name, "_spatial_celltypes.pdf"), width = 14, height = 10)
  p1 <- SpatialDimPlot(seurat_obj, group.by = "final_celltype", 
                       label = FALSE, pt.size.factor = 1.6) +
    ggtitle("Cell Type Spatial Distribution")
  print(p1)
  dev.off()
  
  # 关键marker空间表达
  key_markers <- c("Slc1a3", "Slc17a7", "Gad2", "Tmem119", "Mbp", "Pdgfra", "Cldn5")
  key_markers <- intersect(key_markers, rownames(seurat_obj))
  
  if(length(key_markers) > 0) {
    pdf(paste0(sample_name, "_spatial_key_markers.pdf"), width = 15, height = 12)
    for(i in seq(1, length(key_markers), by = 4)) {
      genes <- key_markers[i:min(i+3, length(key_markers))]
      print(SpatialFeaturePlot(seurat_obj, features = genes, ncol = 2, pt.size.factor = 1.6))
    }
    dev.off()
    cat("✓ 空间特征基因图已保存\n")
  }
}

# ==================== 步骤14: 保存完整对象 ====================

cat("\n========== 步骤14: 保存完整结果 ==========\n")

if(!file.exists(paste0(sample_name, "_complete.rds"))) {
  saveRDS(seurat_obj, paste0(sample_name, "_complete.rds"))
  file_size <- round(file.size(paste0(sample_name, "_complete.rds")) / 1024^3, 2)
  cat("✓ 完整Seurat对象已保存，大小:", file_size, "GB\n")
}

# ==================== 步骤15: 生成综合报告 ====================

cat("\n========== 步骤15: 生成综合报告 ==========\n")

report_file <- paste0(sample_name, "_complete_report.txt")
sink(report_file)

cat("========================================\n")
cat("脑样本完整分析报告\n")
cat("========================================\n\n")
cat("样本:", sample_name, "\n")
cat("分析时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("总细胞数:", ncol(seurat_obj), "\n")
cat("总基因数:", nrow(seurat_obj), "\n\n")

cat("========================================\n")
cat("SingleR注释结果\n")
cat("========================================\n")
singler_counts <- table(seurat_obj$celltype)
for(i in 1:length(singler_counts)) {
  cat(sprintf("  %-20s: %6d (%5.1f%%)\n", 
              names(singler_counts)[i], 
              singler_counts[i], 
              singler_counts[i]/ncol(seurat_obj)*100))
}

cat("\n========================================\n")
cat("手动精细注释结果\n")
cat("========================================\n")
manual_counts <- table(seurat_obj$final_celltype)
for(i in 1:length(manual_counts)) {
  cat(sprintf("  %-20s: %6d (%5.1f%%)\n", 
              names(manual_counts)[i], 
              manual_counts[i], 
              manual_counts[i]/ncol(seurat_obj)*100))
}

# 神经元亚群统计
if(file.exists(paste0(sample_name, "_neuron_proportions.csv"))) {
  cat("\n========================================\n")
  cat("神经元亚群分布\n")
  cat("========================================\n")
  neuron_df <- read.csv(paste0(sample_name, "_neuron_proportions.csv"))
  for(i in 1:nrow(neuron_df)) {
    cat(sprintf("  %-20s: %6d (%5.1f%%)\n", 
                neuron_df$Subtype[i], 
                neuron_df$Count[i], 
                neuron_df$Percentage[i]))
  }
}

sink()
cat("✓ 报告已保存:", report_file, "\n")

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
