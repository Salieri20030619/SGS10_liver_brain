# ==================== 脑样本完整分析脚本（修复版 - 亚型敏感分层共识） ====================
# 文件名: ~/SGS10_analysis/scripts/brain_samples_analysis_optimized_final.R
# 用法: Rscript brain_samples_analysis_optimized_final.R sample_name n_cores [run_fine_annotation]

args <- commandArgs(trailingOnly = TRUE)
sample_name <- args[1]
n_cores <- as.numeric(args[2])
run_fine_annotation <- ifelse(length(args) >= 3, as.logical(args[3]), TRUE)

# 设置输出立即显示
options(echo = TRUE)

cat("========================================\n")
cat("开始处理脑样本（优化版）:", sample_name, "\n")
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

# ==================== 脑细胞注释marker库 ====================

create_brain_marker_db <- function() {
  brain_markers <- list()
  
  # 神经元
  brain_markers$Neurons <- c("Snap25", "Syp", "Syn1", "Rbfox3", "Tubb3", "Map2")
  
  # 兴奋性神经元
  brain_markers$Excitatory_L23 <- c("Slc17a7", "Cux1", "Cux2", "Rorb", "Satb2")
  brain_markers$Excitatory_L45 <- c("Slc17a7", "Fezf2", "Bcl11b", "Etv1", "Tle4")
  brain_markers$Excitatory_L6 <- c("Slc17a7", "Tle4", "Foxp2", "Sulf1")
  
  # 海马体神经元
  brain_markers$Hpc_CA1 <- c("Slc17a7", "Neurod6", "Prox1", "Camk2a", "Calb1")
  brain_markers$Hpc_CA3 <- c("Slc17a7", "Neurod6", "Prox1", "Camk2a", "Calb2")
  brain_markers$Hpc_DG <- c("Slc17a7", "Prox1", "Dcx", "Calb1", "Pcp4")
  
  # 抑制性神经元
  brain_markers$Inhibitory_Sst <- c("Gad1", "Gad2", "Sst", "Cort")
  brain_markers$Inhibitory_Pvalb <- c("Gad1", "Gad2", "Pvalb")
  brain_markers$Inhibitory_Vip <- c("Gad1", "Gad2", "Vip", "Cck")
  brain_markers$Inhibitory_Meis2 <- c("Gad1", "Gad2", "Meis2", "Lhx6")
  
  # 星形胶质细胞
  brain_markers$Astrocytes <- c("Slc1a3", "Gfap", "Aqp4", "Aldh1l1", "S100b")
  brain_markers$Astro_CTX <- c("Slc1a3", "Gfap", "Aqp4", "Aldh1l1", "S100b", "Foxg1")
  brain_markers$Astro_HPC <- c("Slc1a3", "Gfap", "Aqp4", "Aldh1l1", "S100b", "Hpca")
  brain_markers$Astro_HYPO <- c("Slc1a3", "Gfap", "Aqp4", "Aldh1l1", "S100b", "Tll3")
  brain_markers$Astro_THAL <- c("Slc1a3", "Gfap", "Aqp4", "Aldh1l1", "S100b", "Lgr6")
  brain_markers$Astro_WM <- c("Slc1a3", "Gfap", "Aqp4", "Aldh1l1", "S100b", "WhiteMatter")
  
  # 小胶质细胞
  brain_markers$Microglia <- c("Tmem119", "P2ry12", "Cx3cr1", "Hexb", "Csf1r", "Cd68")
  brain_markers$Microglia_Homeostatic <- c("P2ry12", "Tmem119", "Cx3cr1", "Siglech")
  brain_markers$Microglia_M1 <- c("Il1b", "Tnfa", "Cd86", "Nos2")
  brain_markers$Microglia_M2 <- c("Arg1", "Cd163", "Mrc1", "Il10")
  brain_markers$Microglia_DAM <- c("Apoe", "Clec7a", "Lpl", "Cst7")
  
  # 少突胶质细胞
  brain_markers$Oligodendrocytes <- c("Mbp", "Mog", "Plp1", "Cnp", "Mobp")
  brain_markers$OPCs <- c("Pdgfra", "Cspg4", "Olig1", "Olig2")
  
  # 血管相关
  brain_markers$Endothelial <- c("Cldn5", "Pecam1", "Cdh5", "Vwf", "Flt1")
  brain_markers$Pericytes <- c("Rgs5", "Pdgfrb", "Cspg4", "Kcnj8", "Acta2")
  brain_markers$VLMC <- c("Lum", "Dcn", "Col1a1", "Col3a1")
  
  # 室管膜和脉络丛
  brain_markers$Ependymal <- c("Foxj1", "Vim", "Ccdc153", "Cfap126")
  brain_markers$CPECs <- c("Ttr", "Lcn2", "Clcn2", "Cldn3", "Slc2a1")
  
  # 免疫细胞
  brain_markers$T_Cells <- c("Cd3d", "Cd3e", "Cd3g", "Cd4", "Cd8a")
  brain_markers$B_Cells <- c("Cd19", "Ms4a1", "Cd79a", "Cd79b")
  brain_markers$NK_Cells <- c("Nkg7", "Klrb1", "Klrd1", "Klrk1")
  
  # 技术噪声
  brain_markers$TECHO <- c("Malat1", "Neat1", "Mtrnr2l2", "Mt-rnr1")
  
  return(brain_markers)
}

# ==================== 质量检查函数 ====================

quality_check_new_cells <- function(seurat_obj, all_scores, final_assignment, 
                                     old_assignment, subtype_threshold_old) {
  cat("\n========== 阈值调整质量检查 ==========\n")
  
  # 找出新增获得注释的细胞
  old_annotated <- names(old_assignment[!is.na(old_assignment)])
  new_annotated <- names(final_assignment[!is.na(final_assignment)])
  newly_annotated <- setdiff(new_annotated, old_annotated)
  
  cat("原有注释细胞数:", length(old_annotated), "\n")
  cat("新注释细胞数:", length(new_annotated), "\n")
  cat("新增注释细胞数:", length(newly_annotated), "\n")
  
  if(length(newly_annotated) == 0) {
    cat("没有新增细胞，跳过质量检查\n")
    return(NULL)
  }
  
  # 1. 新增细胞得分分布
  cat("\n--- 1. 得分分布检查 ---\n")
  new_scores <- all_scores[newly_annotated, ]
  new_max_scores <- apply(new_scores, 1, max)
  
  cat("新增细胞最高得分统计:\n")
  cat("  最小值:", round(min(new_max_scores), 4), "\n")
  cat("  中位数:", round(median(new_max_scores), 4), "\n")
  cat("  最大值:", round(max(new_max_scores), 4), "\n")
  cat("  25%分位数:", round(quantile(new_max_scores, 0.25), 4), "\n")
  cat("  75%分位数:", round(quantile(new_max_scores, 0.75), 4), "\n")
  
  # 2. 按亚型统计新增细胞
  cat("\n--- 2. 新增细胞亚型分布 ---\n")
  new_assignment <- final_assignment[newly_annotated]
  new_type_counts <- table(new_assignment)
  
  for(ct in names(sort(new_type_counts, decreasing = TRUE))) {
    cat("  ", ct, ": ", new_type_counts[ct], "个细胞\n")
  }
  
  # 3. 检查新增细胞的质量指标
  cat("\n--- 3. 新增细胞质量指标 ---\n")
  
  if("percent.mt" %in% colnames(seurat_obj@meta.data)) {
    original_mt <- seurat_obj@meta.data[old_annotated, "percent.mt"]
    new_mt <- seurat_obj@meta.data[newly_annotated, "percent.mt"]
    
    cat("线粒体比例:\n")
    cat("  原有细胞: 中位数=", round(median(original_mt), 2), 
        "%, 范围=", round(range(original_mt), 2), "\n")
    cat("  新增细胞: 中位数=", round(median(new_mt), 2), 
        "%, 范围=", round(range(new_mt), 2), "\n")
    
    mt_test <- wilcox.test(original_mt, new_mt)
    if(mt_test$p.value < 0.05) {
      cat("  ⚠️ 新增细胞线粒体比例显著不同 (p =", format(mt_test$p.value, digits = 3), ")\n")
    } else {
      cat("  ✓ 线粒体比例无显著差异\n")
    }
  }
  
  if("nFeature_RNA" %in% colnames(seurat_obj@meta.data)) {
    original_genes <- seurat_obj@meta.data[old_annotated, "nFeature_RNA"]
    new_genes <- seurat_obj@meta.data[newly_annotated, "nFeature_RNA"]
    
    cat("\n检测基因数:\n")
    cat("  原有细胞: 中位数=", round(median(original_genes), 0), 
        ", 范围=", round(range(original_genes), 0), "\n")
    cat("  新增细胞: 中位数=", round(median(new_genes), 0), 
        ", 范围=", round(range(new_genes), 0), "\n")
    
    genes_test <- wilcox.test(original_genes, new_genes)
    if(genes_test$p.value < 0.05) {
      cat("  ⚠️ 新增细胞基因数显著不同 (p =", format(genes_test$p.value, digits = 3), ")\n")
    } else {
      cat("  ✓ 基因数无显著差异\n")
    }
  }
  
  # 4. 生成可视化报告
  cat("\n--- 4. 生成可视化报告 ---\n")
  
  pdf(paste0(sample_name, "_threshold_quality_check.pdf"), width = 14, height = 10)
  
  p1 <- DimPlot(seurat_obj, reduction = "umap", 
                cells.highlight = newly_annotated, 
                sizes.highlight = 1, 
                cols.highlight = "red",
                cols = "grey") +
    ggtitle(paste0("Newly annotated cells (n=", length(newly_annotated), ")")) +
    theme(plot.title = element_text(size = 14))
  
  old_scores <- all_scores[old_annotated, ]
  old_max <- apply(old_scores, 1, max)
  
  score_df <- data.frame(
    Score = c(old_max, new_max_scores),
    Group = c(rep("Original", length(old_max)), rep("New", length(new_max_scores)))
  )
  
  p2 <- ggplot(score_df, aes(x = Score, fill = Group)) +
    geom_density(alpha = 0.5) +
    geom_vline(xintercept = subtype_threshold_old, linetype = "dashed", color = "blue") +
    annotate("text", x = subtype_threshold_old, y = 1, label = "Old Threshold", 
             hjust = -0.1, size = 3, color = "blue") +
    ggtitle("Score Distribution Comparison") +
    theme_minimal()
  
  print(p1)
  print(p2)
  dev.off()
  
  cat("✓ 质量检查报告已保存: ", sample_name, "_threshold_quality_check.pdf\n", sep="")
  
  # 5. 保存新增细胞列表
  new_cells_df <- data.frame(
    Cell = newly_annotated,
    Assigned_Type = new_assignment,
    Max_Score = new_max_scores,
    Percent_MT = seurat_obj@meta.data[newly_annotated, "percent.mt"],
    nFeature_RNA = seurat_obj@meta.data[newly_annotated, "nFeature_RNA"]
  )
  write.csv(new_cells_df, paste0(sample_name, "_newly_annotated_cells.csv"), row.names = FALSE)
  cat("✓ 新增细胞列表已保存\n")
  
  return(newly_annotated)
}

# ==================== 安全采样函数 ====================

safe_sample_cells <- function(obj, cluster_id, max_cells = 50) {
  cells <- WhichCells(obj, idents = cluster_id)
  n_cells <- length(cells)
  
  if(n_cells == 0) {
    return(character(0))
  }
  
  if(n_cells <= max_cells) {
    return(cells)
  } else {
    return(sample(cells, max_cells))
  }
}

# ==================== 绘制亚型热图的辅助函数（Seurat 5兼容版） ====================

plot_subtype_heatmap <- function(seurat_obj, sample_name, subtypes_found) {
  cat("\n绘制亚型热图...\n")
  
  tryCatch({
    # 设置身份
    Idents(seurat_obj) <- seurat_obj$final_celltype
    
    # 获取亚型细胞
    subtype_cells <- WhichCells(seurat_obj, idents = subtypes_found)
    
    if(length(subtype_cells) < 50) {
      cat("  警告: 亚型细胞数量不足50，跳过热图\n")
      return(FALSE)
    }
    
    cat("  亚型细胞数:", length(subtype_cells), "\n")
    
    # 创建临时子集对象
    temp_obj <- subset(seurat_obj, cells = subtype_cells)
    
    # 检查是否有足够的细胞和基因
    if(ncol(temp_obj) < 10 || nrow(temp_obj) < 10) {
      cat("  警告: 细胞或基因数量不足，跳过热图\n")
      return(FALSE)
    }
    
    # 计算marker基因
    cat("  计算差异表达基因...\n")
    markers <- tryCatch({
      FindAllMarkers(temp_obj, only.pos = TRUE, 
                     min.pct = 0.25, logfc.threshold = 0.25,
                     verbose = FALSE)
    }, error = function(e) {
      cat("  警告: 无法计算marker -", e$message, "\n")
      return(data.frame())
    })
    
    if(nrow(markers) == 0) {
      cat("  警告: 未找到差异表达基因，使用高变基因\n")
      if(length(VariableFeatures(temp_obj)) == 0) {
        VariableFeatures(temp_obj) <- rownames(temp_obj)[1:min(100, nrow(temp_obj))]
      }
      top_genes <- VariableFeatures(temp_obj)[1:min(20, length(VariableFeatures(temp_obj)))]
    } else {
      # 每个亚型取top 3个基因
      top_genes <- markers %>% 
        group_by(cluster) %>% 
        top_n(3, avg_log2FC) %>% 
        pull(gene) %>% 
        unique()
    }
    
    cat("  找到", length(top_genes), "个top marker基因\n")
    
    if(length(top_genes) < 2) {
      cat("  警告: 基因数量不足，跳过热图\n")
      return(FALSE)
    }
    
    # 限制基因数量避免热图过大
    if(length(top_genes) > 50) {
      top_genes <- top_genes[1:50]
      cat("  限制为前50个基因\n")
    }
    
    # 获取表达矩阵（使用layer参数，兼容Seurat 5）
    cat("  获取表达矩阵...\n")
    expr_matrix <- GetAssayData(temp_obj, assay = "RNA", layer = "data")
    expr_matrix_subset <- expr_matrix[top_genes, , drop = FALSE]
    
    # 确保矩阵不为空
    if(nrow(expr_matrix_subset) == 0 || ncol(expr_matrix_subset) == 0) {
      cat("  错误: 表达矩阵为空\n")
      return(FALSE)
    }
    
    # 计算z-score
    cat("  计算z-score...\n")
    expr_matrix_z <- t(scale(t(as.matrix(expr_matrix_subset))))
    expr_matrix_z[is.na(expr_matrix_z)] <- 0
    
    # 创建注释数据框
    annotation_df <- data.frame(
      CellType = Idents(temp_obj)[colnames(expr_matrix_z)]
    )
    rownames(annotation_df) <- colnames(expr_matrix_z)
    
    # 获取颜色
    cell_types <- unique(annotation_df$CellType)
    colors <- setNames(rainbow(length(cell_types)), cell_types)
    
    # 如果细胞数太多，随机采样
    if(ncol(expr_matrix_z) > 500) {
      cat("  细胞数过多，随机采样500个\n")
      set.seed(42)
      sample_cols <- sample(colnames(expr_matrix_z), 500)
      expr_matrix_z <- expr_matrix_z[, sample_cols]
      annotation_df <- annotation_df[sample_cols, , drop = FALSE]
    }
    
    # 打开PDF并绘制热图
    pdf(paste0(sample_name, "_subtype_heatmap.pdf"), width = 14, height = 10)
    
    cat("  绘制热图...\n")
    pheatmap::pheatmap(expr_matrix_z,
                       annotation_col = annotation_df,
                       annotation_colors = list(CellType = colors),
                       show_rownames = TRUE,
                       show_colnames = FALSE,
                       cluster_rows = TRUE,
                       cluster_cols = TRUE,
                       fontsize_row = 8,
                       main = paste0("Subtype Marker Genes (", length(subtypes_found), " subtypes)"),
                       color = colorRampPalette(c("blue", "white", "red"))(100))
    
    dev.off()
    cat("  ✓ 热图绘制成功\n")
    return(TRUE)
    
  }, error = function(e) {
    cat("  错误: 热图绘制失败 -", e$message, "\n")
    return(FALSE)
  })
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

# ==================== 步骤8: 手动精细注释（修复版） ====================

cat("\n========== 步骤8: 手动精细注释（修复版） ==========\n")

if(check_step("manual_annotation_done")) {
  cat("检测到手动注释结果存在，加载文件...\n")
  seurat_obj <- readRDS(paste0(sample_name, "_manual_annotated.rds"))
} else {
  cat("创建脑细胞marker数据库...\n")
  brain_markers <- create_brain_marker_db()
  cat("定义了", length(brain_markers), "种细胞类型\n")
  
  # 分离精细亚型和大类
  subtype_names <- c(
    # 兴奋性神经元亚型
    "Excitatory_L23", "Excitatory_L45", "Excitatory_L6",
    # 海马体神经元亚型
    "Hpc_CA1", "Hpc_CA3", "Hpc_DG",
    # 抑制性神经元亚型
    "Inhibitory_Sst", "Inhibitory_Pvalb", "Inhibitory_Vip", "Inhibitory_Meis2",
    # 星形胶质细胞亚型
    "Astro_CTX", "Astro_HPC", "Astro_HYPO", "Astro_THAL", "Astro_WM",
    # 小胶质细胞亚型
    "Microglia_Homeostatic", "Microglia_M1", "Microglia_M2", "Microglia_DAM",
    # 其他亚型
    "OPCs"
  )
  
  major_names <- c(
    "Neurons", "Astrocytes", "Microglia", "Oligodendrocytes",
    "Endothelial", "Pericytes", "VLMC", "Ependymal", "CPECs",
    "T_Cells", "B_Cells", "NK_Cells"
  )
  
  cat("精细亚型数:", length(subtype_names), "\n")
  cat("大类数:", length(major_names), "\n")
  
  cat("获取表达数据...\n")
  expr_data <- GetAssayData(seurat_obj, layer = "data")
  
  cat("计算每个细胞类型的marker得分...\n")
  
  # 计算所有细胞类型的得分
  all_scores <- matrix(0, nrow = ncol(seurat_obj), ncol = length(brain_markers))
  colnames(all_scores) <- names(brain_markers)
  rownames(all_scores) <- colnames(seurat_obj)
  
  for(ct in names(brain_markers)) {
    markers <- brain_markers[[ct]]
    markers_present <- intersect(markers, rownames(expr_data))
    if(length(markers_present) > 0) {
      scores <- colMeans(expr_data[markers_present, , drop = FALSE])
      all_scores[, ct] <- scores
      cat("  ✓", ct, ": 找到", length(markers_present), "个marker\n")
    } else {
      for(m in markers) {
        m_idx <- grep(m, rownames(expr_data), ignore.case = TRUE, value = TRUE)
        if(length(m_idx) > 0) {
          scores <- colMeans(expr_data[m_idx, , drop = FALSE])
          all_scores[, ct] <- all_scores[, ct] + scores
        }
      }
      cat("  ⚠️", ct, ": 使用模糊匹配\n")
    }
  }
  
  # 分离亚型得分和大类得分
  subtype_scores <- all_scores[, subtype_names, drop = FALSE]
  major_scores <- all_scores[, major_names, drop = FALSE]
  
  # 计算阈值
  cat("计算阈值...\n")
  all_max_scores <- apply(all_scores, 1, max)
  subtype_max_scores <- apply(subtype_scores, 1, max)
  major_max_scores <- apply(major_scores, 1, max)
  
  # 保存旧阈值用于对比
  subtype_threshold_old <- quantile(all_max_scores[all_max_scores > 0], 0.1, na.rm = TRUE)
  
  # 亚型阈值从0.1改为0.05
  subtype_threshold <- quantile(all_max_scores[all_max_scores > 0], 0.05, na.rm = TRUE)
  if(is.na(subtype_threshold) || subtype_threshold < 0.03) subtype_threshold <- 0.03
  
  # 大类阈值
  major_threshold <- quantile(all_max_scores[all_max_scores > 0], 0.02, na.rm = TRUE)
  if(is.na(major_threshold) || major_threshold < 0.02) major_threshold <- 0.02
  
  cat("  旧亚型阈值 (0.1分位数):", round(subtype_threshold_old, 4), "\n")
  cat("  新亚型阈值 (0.05分位数):", round(subtype_threshold, 4), "\n")
  cat("  大类阈值:", round(major_threshold, 4), "\n")
  
  # 亚型优先分配策略
  cat("亚型优先分配...\n")
  final_assignment <- rep(NA, ncol(seurat_obj))
  names(final_assignment) <- colnames(seurat_obj)
  
  # 同时保存旧阈值下的分配结果用于对比
  old_assignment <- rep(NA, ncol(seurat_obj))
  names(old_assignment) <- colnames(seurat_obj)
  
  for(cell in colnames(seurat_obj)) {
    # 第一步：检查亚型
    cell_subtype_scores <- subtype_scores[cell, ]
    max_subtype_score <- max(cell_subtype_scores)
    best_subtype <- names(which.max(cell_subtype_scores))
    
    # 第二步：检查大类
    cell_major_scores <- major_scores[cell, ]
    max_major_score <- max(cell_major_scores)
    best_major <- names(which.max(cell_major_scores))
    
    # 新阈值分配
    if(max_subtype_score >= subtype_threshold) {
      final_assignment[cell] <- best_subtype
    } else if(max_major_score >= major_threshold) {
      final_assignment[cell] <- best_major
    } else {
      final_assignment[cell] <- NA
    }
    
    # 旧阈值分配（用于对比）
    if(max_subtype_score >= subtype_threshold_old) {
      old_assignment[cell] <- best_subtype
    } else if(max_major_score >= major_threshold) {
      old_assignment[cell] <- best_major
    } else {
      old_assignment[cell] <- NA
    }
  }
  
  # 质量检查
  cat("\n执行阈值调整质量检查...\n")
  newly_annotated <- quality_check_new_cells(
    seurat_obj, all_scores, final_assignment, 
    old_assignment, subtype_threshold_old
  )
  
  # ==================== 亚型敏感分层共识注释 ====================
  cat("\n基于cluster的分层共识注释（亚型敏感）...\n")
  cluster_ids <- Idents(seurat_obj)
  cluster_levels <- levels(cluster_ids)
  
  # 初始化共识注释为final_assignment
  consensus_celltype <- final_assignment
  
  # 统计信息
  consensus_stats <- data.frame(
    Cluster = cluster_levels,
    Total_Cells = NA,
    Majority_Type = NA,
    Majority_Fraction = NA,
    Action = NA
  )
  
  for(cl in cluster_levels) {
    cells_in_cluster <- WhichCells(seurat_obj, idents = cl)
    total_cells <- length(cells_in_cluster)
    consensus_stats[consensus_stats$Cluster == cl, "Total_Cells"] <- total_cells
    
    if(length(cells_in_cluster) >= 5) {
      assignments <- final_assignment[cells_in_cluster]
      assignments_clean <- assignments[!is.na(assignments)]
      
      if(length(assignments_clean) >= 3) {
        # 分离亚型和大类注释
        subtype_assign <- assignments_clean[assignments_clean %in% subtype_names]
        major_assign <- assignments_clean[assignments_clean %in% major_names]
        
        # 统计各类型频率
        all_types <- c(subtype_assign, major_assign)
        freq_table <- table(all_types)
        majority_type <- names(freq_table)[which.max(freq_table)]
        majority_frac <- max(freq_table) / sum(freq_table)
        
        consensus_stats[consensus_stats$Cluster == cl, "Majority_Type"] <- majority_type
        consensus_stats[consensus_stats$Cluster == cl, "Majority_Fraction"] <- round(majority_frac, 3)
        
        # 决策逻辑
        if(length(subtype_assign) >= length(assignments_clean) * 0.5) {
          if(length(unique(subtype_assign)) == 1) {
            consensus_celltype[cells_in_cluster] <- unique(subtype_assign)
            consensus_stats[consensus_stats$Cluster == cl, "Action"] <- "Unified_Subtype"
            cat("  Cluster", cl, ": 统一为亚型", unique(subtype_assign), "\n")
          } else {
            consensus_stats[consensus_stats$Cluster == cl, "Action"] <- "Preserved_Diversity"
            cat("  Cluster", cl, ": 保留亚型多样性\n")
          }
        } else if(length(major_assign) >= 3 && majority_frac > 0.7) {
          consensus_celltype[cells_in_cluster] <- majority_type
          consensus_stats[consensus_stats$Cluster == cl, "Action"] <- "Unified_Major"
          cat("  Cluster", cl, ": 统一为大类", majority_type, "(", round(majority_frac*100, 1), "%)\n")
        } else if(length(subtype_assign) > 0) {
          consensus_stats[consensus_stats$Cluster == cl, "Action"] <- "Kept_Subtypes"
          cat("  Cluster", cl, ": 保留", length(unique(subtype_assign)), "种亚型\n")
        } else {
          consensus_stats[consensus_stats$Cluster == cl, "Action"] <- "Kept_Major"
          cat("  Cluster", cl, ": 保留大类多样性\n")
        }
      } else {
        consensus_stats[consensus_stats$Cluster == cl, "Action"] <- "Insufficient_Data"
        cat("  Cluster", cl, ": 注释不足\n")
      }
    } else {
      consensus_stats[consensus_stats$Cluster == cl, "Action"] <- "Small_Cluster"
      cat("  Cluster", cl, ": 小cluster\n")
    }
  }
  
  # 确保亚型不被大类覆盖
  subtype_to_major <- c(
    "Excitatory_L23" = "Neurons", "Excitatory_L45" = "Neurons", "Excitatory_L6" = "Neurons",
    "Hpc_CA1" = "Neurons", "Hpc_CA3" = "Neurons", "Hpc_DG" = "Neurons",
    "Inhibitory_Sst" = "Neurons", "Inhibitory_Pvalb" = "Neurons", 
    "Inhibitory_Vip" = "Neurons", "Inhibitory_Meis2" = "Neurons",
    "Astro_CTX" = "Astrocytes", "Astro_HPC" = "Astrocytes", 
    "Astro_HYPO" = "Astrocytes", "Astro_THAL" = "Astrocytes", "Astro_WM" = "Astrocytes",
    "Microglia_Homeostatic" = "Microglia", "Microglia_M1" = "Microglia",
    "Microglia_M2" = "Microglia", "Microglia_DAM" = "Microglia",
    "OPCs" = "Oligodendrocytes"
  )
  
  protection_count <- 0
  for(cell in names(final_assignment)) {
    if(!is.na(final_assignment[cell]) && final_assignment[cell] %in% names(subtype_to_major)) {
      if(!is.na(consensus_celltype[cell]) && 
         consensus_celltype[cell] == subtype_to_major[final_assignment[cell]]) {
        consensus_celltype[cell] <- final_assignment[cell]
        protection_count <- protection_count + 1
      }
    }
  }
  
  cat("\n亚型保护: 恢复了", protection_count, "个细胞的亚型注释\n")
  
  # 保存共识统计
  write.csv(consensus_stats, paste0(sample_name, "_consensus_stats.csv"), row.names = FALSE)
  cat("✓ 共识统计已保存\n")
  
  # 应用共识注释
  seurat_obj$manual_celltype <- consensus_celltype
  seurat_obj$final_celltype <- consensus_celltype
  
  # 填补完全未分配的细胞（NA）
  na_idx <- is.na(seurat_obj$final_celltype)
  if(sum(na_idx) > 0) {
    cat("\n填补", sum(na_idx), "个未分配细胞...\n")
    seurat_obj$final_celltype[na_idx] <- seurat_obj$celltype[na_idx]
  }
  
  Idents(seurat_obj) <- seurat_obj$final_celltype
  
  cat("\n最终细胞类型分布:\n")
  print(table(Idents(seurat_obj)))
  
  # 统计亚型和大类数量
  final_types <- names(table(Idents(seurat_obj)))
  subtypes_found <- intersect(final_types, subtype_names)
  majors_found <- intersect(final_types, major_names)
  
  cat("\n发现的精细亚型 (", length(subtypes_found), "种):\n")
  print(subtypes_found)
  
  cat("\n发现的大类 (", length(majors_found), "种):\n")
  print(majors_found)
  
  # 保存比例
  prop_df <- data.frame(
    CellType = names(table(Idents(seurat_obj))),
    Count = as.numeric(table(Idents(seurat_obj))),
    Percentage = as.numeric(prop.table(table(Idents(seurat_obj))) * 100)
  )
  write.csv(prop_df, paste0(sample_name, "_manual_proportions.csv"), row.names = FALSE)
  
  # 绘制对比图
  pdf(paste0(sample_name, "_annotation_comparison.pdf"), width = 16, height = 8)
  p1 <- DimPlot(seurat_obj, reduction = "umap", group.by = "final_celltype", 
                label = TRUE, label.size = 3, repel = TRUE) +
    ggtitle("Subtype-Sensitive Consensus Annotation")
  p2 <- DimPlot(seurat_obj, reduction = "umap", group.by = "celltype", 
                label = TRUE, label.size = 3, repel = TRUE) +
    ggtitle("SingleR Annotation")
  print(p1 | p2)
  dev.off()
  
  # 绘制亚型热图（使用修复后的函数）
  if(length(subtypes_found) > 1) {
    cat("\n绘制亚型热图...\n")
    tryCatch({
      plot_subtype_heatmap(seurat_obj, sample_name, subtypes_found)
    }, error = function(e) {
      cat("  警告: 热图绘制失败 -", e$message, "\n")
      cat("  继续保存对象...\n")
    })
  } else {
    cat("\n亚型数量不足，跳过热图绘制\n")
  }
  
  # ========== 关键修复：保存前的完整检查 ==========
  cat("\n=== 保存前验证 ===\n")
  
  # 1. 确保 final_celltype 存在
  if(!"final_celltype" %in% colnames(seurat_obj@meta.data)) {
    cat("警告: final_celltype不存在，从Idents创建\n")
    seurat_obj$final_celltype <- as.character(Idents(seurat_obj))
  }
  
  # 2. 确保没有NA
  na_count <- sum(is.na(seurat_obj$final_celltype))
  if(na_count > 0) {
    cat("警告: 发现", na_count, "个NA，使用celltype填充\n")
    na_idx <- is.na(seurat_obj$final_celltype)
    seurat_obj$final_celltype[na_idx] <- as.character(seurat_obj$celltype[na_idx])
  }
  
  # 3. 保存得分信息（向量）
  if(exists("all_max_scores") && length(all_max_scores) == ncol(seurat_obj)) {
    seurat_obj$subtype_max_score <- all_max_scores
    cat("✓ subtype_max_score已添加\n")
  } else {
    cat("警告: all_max_scores不存在或长度不匹配，跳过\n")
  }
  
  # 4. 保存阈值信息（存储到misc槽位，避免元数据错误）
  if(exists("subtype_threshold")) {
    seurat_obj@misc$subtype_threshold <- subtype_threshold
    cat("✓ subtype_threshold已保存到misc\n")
  }
  if(exists("major_threshold")) {
    seurat_obj@misc$major_threshold <- major_threshold
    cat("✓ major_threshold已保存到misc\n")
  }
  if(exists("subtype_threshold_old")) {
    seurat_obj@misc$subtype_threshold_old <- subtype_threshold_old
    cat("✓ subtype_threshold_old已保存到misc\n")
  }
  
  seurat_obj@misc$annotation_method <- "subtype_sensitive_consensus"
  
  # 5. 最终验证
  cat("\n最终验证:\n")
  cat("  final_celltype存在:", "final_celltype" %in% colnames(seurat_obj@meta.data), "\n")
  cat("  final_celltype非NA数量:", sum(!is.na(seurat_obj$final_celltype)), "/", ncol(seurat_obj), "\n")
  cat("  final_celltype唯一值:", length(unique(seurat_obj$final_celltype)), "\n")
  cat("  Idents唯一值:", length(unique(Idents(seurat_obj))), "\n")
  cat("  元数据列数:", ncol(seurat_obj@meta.data), "\n")
  cat("  misc槽位内容:", paste(names(seurat_obj@misc), collapse=", "), "\n")
  cat("==================\n")
  
  # 6. 保存对象
  cat("\n保存对象...\n")
  saveRDS(seurat_obj, paste0(sample_name, "_manual_annotated.rds"))
  
  # 7. 验证保存
  if(file.exists(paste0(sample_name, "_manual_annotated.rds"))) {
    file_size <- file.size(paste0(sample_name, "_manual_annotated.rds")) / 1024^3
    cat("✓ 手动注释结果已保存，大小:", round(file_size, 2), "GB\n")
    
    # 验证可以重新加载
    test_obj <- tryCatch({
      readRDS(paste0(sample_name, "_manual_annotated.rds"))
    }, error = function(e) {
      cat("  警告: 保存的文件无法读取 -", e$message, "\n")
      return(NULL)
    })
    
    if(!is.null(test_obj)) {
      cat("  验证: 文件可正常读取\n")
      if("final_celltype" %in% colnames(test_obj@meta.data)) {
        cat("  验证: final_celltype存在\n")
      }
      if("subtype_max_score" %in% colnames(test_obj@meta.data)) {
        cat("  验证: subtype_max_score存在\n")
      }
    }
    
    save_checkpoint("manual_annotation_done")
  } else {
    cat("✗ 保存失败！文件未生成\n")
  }
}

# ==================== 步骤9-13: 精细注释和报告 ====================

# 步骤9: 小胶质细胞亚群分析
if(run_fine_annotation) {
  cat("\n========== 步骤9: 小胶质细胞亚群分析 ==========\n")
  
  if(check_step("microglia_done", force = TRUE)) {
    cat("检测到小胶质细胞亚群结果存在，跳过\n")
  } else {
    microglia_names <- c("Microglia", "Microglia_Homeostatic", "Microglia_M1", "Microglia_M2", "Microglia_DAM")
    microglia_present <- microglia_names[microglia_names %in% unique(seurat_obj$final_celltype)]
    
    if(length(microglia_present) > 0) {
      microglia_cells <- WhichCells(seurat_obj, expression = final_celltype %in% microglia_present)
      cat("找到小胶质细胞数量:", length(microglia_cells), "\n")
      
      if(length(microglia_cells) >= 50) {
        microglia_obj <- subset(seurat_obj, cells = microglia_cells)
        
        cat("\n重新处理小胶质细胞...\n")
        microglia_obj <- NormalizeData(microglia_obj, verbose = FALSE)
        microglia_obj <- FindVariableFeatures(microglia_obj, nfeatures = 800, verbose = FALSE)
        microglia_obj <- ScaleData(microglia_obj, verbose = FALSE)
        microglia_obj <- RunPCA(microglia_obj, npcs = 15, verbose = FALSE)
        
        pdf(paste0(sample_name, "_microglia_elbow.pdf"), width = 8, height = 6)
        print(ElbowPlot(microglia_obj, ndims = 15))
        dev.off()
        
        n_pcs <- min(10, ncol(microglia_obj[["pca"]]@cell.embeddings))
        microglia_obj <- FindNeighbors(microglia_obj, dims = 1:n_pcs, verbose = FALSE)
        microglia_obj <- FindClusters(microglia_obj, resolution = 0.3, verbose = FALSE)
        microglia_obj <- RunUMAP(microglia_obj, dims = 1:n_pcs, verbose = FALSE)
        
        pdf(paste0(sample_name, "_microglia_states.pdf"), width = 12, height = 8)
        print(DimPlot(microglia_obj, reduction = "umap", label = TRUE) +
                ggtitle(paste0(sample_name, " Microglia Clusters")))
        dev.off()
        
        saveRDS(microglia_obj, paste0(sample_name, "_microglia_subset.rds"))
        cat("✓ 小胶质细胞亚群对象已保存\n")
        save_checkpoint("microglia_done")
      } else {
        cat("小胶质细胞数量不足50，跳过\n")
      }
    } else {
      cat("未找到小胶质细胞\n")
    }
  }
}

# 步骤10: 神经元亚群精细注释
if(run_fine_annotation) {
  cat("\n========== 步骤10: 神经元亚群精细注释 ==========\n")
  
  if(check_step("neuron_done", force = TRUE)) {
    cat("检测到神经元亚群结果存在，跳过\n")
  } else {
    neuron_names <- c("Neurons", "Excitatory_L23", "Excitatory_L45", "Excitatory_L6",
                      "Hpc_CA1", "Hpc_CA3", "Hpc_DG", "Inhibitory_Sst", "Inhibitory_Pvalb", "Inhibitory_Vip")
    neuron_present <- neuron_names[neuron_names %in% unique(seurat_obj$final_celltype)]
    
    if(length(neuron_present) > 0) {
      neuron_cells <- WhichCells(seurat_obj, expression = final_celltype %in% neuron_present)
      cat("神经元数量:", length(neuron_cells), "\n")
      
      if(length(neuron_cells) >= 200) {
        neuron_obj <- subset(seurat_obj, cells = neuron_cells)
        
        cat("\n重新处理神经元亚群...\n")
        neuron_obj <- NormalizeData(neuron_obj, verbose = FALSE)
        neuron_obj <- FindVariableFeatures(neuron_obj, nfeatures = 1000, verbose = FALSE)
        neuron_obj <- ScaleData(neuron_obj, verbose = FALSE)
        neuron_obj <- RunPCA(neuron_obj, npcs = 20, verbose = FALSE)
        
        n_pcs <- min(15, ncol(neuron_obj[["pca"]]@cell.embeddings))
        neuron_obj <- FindNeighbors(neuron_obj, dims = 1:n_pcs, verbose = FALSE)
        neuron_obj <- FindClusters(neuron_obj, resolution = 0.5, verbose = FALSE)
        neuron_obj <- RunUMAP(neuron_obj, dims = 1:n_pcs, verbose = FALSE)
        
        pdf(paste0(sample_name, "_neuron_subtypes.pdf"), width = 12, height = 8)
        print(DimPlot(neuron_obj, reduction = "umap", label = TRUE) +
                ggtitle(paste0(sample_name, " Neuron Clusters")))
        dev.off()
        
        saveRDS(neuron_obj, paste0(sample_name, "_neuron_subset.rds"))
        cat("✓ 神经元亚群对象已保存\n")
        save_checkpoint("neuron_done")
      } else {
        cat("神经元数量不足200，跳过\n")
      }
    } else {
      cat("未找到神经元\n")
    }
  }
}

# 步骤11: 星形胶质细胞亚群精细注释
if(run_fine_annotation) {
  cat("\n========== 步骤11: 星形胶质细胞亚群精细注释 ==========\n")
  
  if(check_step("astro_done", force = TRUE)) {
    cat("检测到星形胶质细胞亚群结果存在，跳过\n")
  } else {
    astro_names <- c("Astrocytes", "Astro_CTX", "Astro_HPC", "Astro_HYPO", "Astro_THAL", "Astro_WM")
    astro_present <- astro_names[astro_names %in% unique(seurat_obj$final_celltype)]
    
    if(length(astro_present) > 0) {
      astro_cells <- WhichCells(seurat_obj, expression = final_celltype %in% astro_present)
      cat("星形胶质细胞数量:", length(astro_cells), "\n")
      
      if(length(astro_cells) >= 150) {
        astro_obj <- subset(seurat_obj, cells = astro_cells)
        
        cat("\n重新处理星形胶质细胞亚群...\n")
        astro_obj <- NormalizeData(astro_obj, verbose = FALSE)
        astro_obj <- FindVariableFeatures(astro_obj, nfeatures = 800, verbose = FALSE)
        astro_obj <- ScaleData(astro_obj, verbose = FALSE)
        astro_obj <- RunPCA(astro_obj, npcs = 15, verbose = FALSE)
        
        n_pcs <- min(10, ncol(astro_obj[["pca"]]@cell.embeddings))
        astro_obj <- FindNeighbors(astro_obj, dims = 1:n_pcs, verbose = FALSE)
        astro_obj <- FindClusters(astro_obj, resolution = 0.3, verbose = FALSE)
        astro_obj <- RunUMAP(astro_obj, dims = 1:n_pcs, verbose = FALSE)
        
        pdf(paste0(sample_name, "_astro_regions.pdf"), width = 12, height = 8)
        print(DimPlot(astro_obj, reduction = "umap", label = TRUE) +
                ggtitle(paste0(sample_name, " Astrocyte Clusters")))
        dev.off()
        
        saveRDS(astro_obj, paste0(sample_name, "_astro_subset.rds"))
        cat("✓ 星形胶质细胞亚群对象已保存\n")
        save_checkpoint("astro_done")
      } else {
        cat("星形胶质细胞数量不足150，跳过\n")
      }
    } else {
      cat("未找到星形胶质细胞\n")
    }
  }
}

# ==================== 步骤12: 保存完整对象 ====================

cat("\n========== 步骤12: 保存完整结果 ==========\n")

if(!file.exists(paste0(sample_name, "_complete.rds"))) {
  saveRDS(seurat_obj, paste0(sample_name, "_complete.rds"))
  file_size <- round(file.size(paste0(sample_name, "_complete.rds")) / 1024^3, 2)
  cat("✓ 完整Seurat对象已保存，大小:", round(file_size, 2), "GB\n")
}

# ==================== 步骤13: 生成综合报告 ====================

cat("\n========== 步骤13: 生成综合报告 ==========\n")

report_file <- paste0(sample_name, "_complete_report.txt")
sink(report_file)

cat("========================================\n")
cat("脑样本完整分析报告（亚型敏感分层共识版）\n")
cat("========================================\n\n")
cat("样本:", sample_name, "\n")
cat("分析时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("总细胞数:", ncol(seurat_obj), "\n")
cat("总基因数:", nrow(seurat_obj), "\n")
if(!is.null(seurat_obj@misc$subtype_threshold)) {
  cat("亚型阈值 (0.05分位数):", round(seurat_obj@misc$subtype_threshold, 4), "\n")
}
if(!is.null(seurat_obj@misc$major_threshold)) {
  cat("大类阈值:", round(seurat_obj@misc$major_threshold, 4), "\n\n")
}

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
