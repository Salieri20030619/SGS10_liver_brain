# ==================== 肝脏样本完整分析脚本（从头开始版） ====================
# 文件名: ~/SGS10_analysis/scripts/liver_complete_analysis_clean.R
# 用法: Rscript liver_complete_analysis_clean.R sample_name n_cores

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

# ==================== 删除旧的中间文件和检查点 ====================
cat("清理旧的中间文件和检查点...\n")
old_files <- list.files(pattern = paste0(sample_name, "_.*\\.(rds|pdf|csv|txt)$"))
if(length(old_files) > 0) {
  for(f in old_files) {
    file.remove(f)
    cat("  删除:", f, "\n")
  }
}
if(file.exists(paste0(sample_name, "_checkpoint.txt"))) {
  file.remove(paste0(sample_name, "_checkpoint.txt"))
  cat("  删除检查点文件\n")
}
cat("✓ 清理完成\n\n")

# ==================== 数据路径 ====================
data_paths <- list(
  gan56 = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/gan56/gan_56_filtered_feature_bc_matrix",
  liver3 = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/liver_samples/25090559_Liver_3_expression/Liver_3_filtered_feature_bc_matrix",
  liver8 = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/liver_samples/25090559_Liver_8_expression/Liver_8_filtered_feature_bc_matrix",
  bla_liver = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/liver_samples/25090559_Bla_liver_1_expression/Bla_liver_1_filtered_feature_bc_matrix",
  con_liver = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/liver_samples/25090559_Con_liver_1_expression/Con_liver_1_filtered_feature_bc_matrix"
)

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
    "ALB" = "Alb", "AFP" = "Afp", "TTR" = "Ttr", "HNF4A" = "Hnf4a", "HNF1A" = "Hnf1a",
    "CYP2E1" = "Cyp2e1", "CYP1A2" = "Cyp1a2", "CYP3A4" = "Cyp3a11",
    "PECAM1" = "Pecam1", "CDH5" = "Cdh5", "VWF" = "Vwf", "KDR" = "Kdr",
    "CLEC4G" = "Clec4g", "LYVE1" = "Lyve1",
    "PDGFRA" = "Pdgfra", "PDGFRB" = "Pdgfrb", "DES" = "Des", "ACTA2" = "Acta2",
    "COL1A1" = "Col1a1", "LRAT" = "Lrat",
    "CD68" = "Cd68", "CD163" = "Cd163", "VSIG4" = "Vsig4", "C1QA" = "C1qa",
    "MRC1" = "Mrc1", "TREM2" = "Trem2", "GPNMB" = "Gpnmb",
    "CD3D" = "Cd3d", "CD3E" = "Cd3e", "CD4" = "Cd4", "CD8A" = "Cd8a",
    "FOXP3" = "Foxp3",
    "CD19" = "Cd19", "MS4A1" = "Ms4a1", "CD79A" = "Cd79a",
    "NKG7" = "Nkg7", "KLRB1" = "Klrb1",
    "CD14" = "Cd14", "LYZ" = "Lyz", "S100A8" = "S100a8", "S100A9" = "S100a9",
    "KRT7" = "Krt7", "KRT8" = "Krt8", "KRT19" = "Krt19", "EPCAM" = "Epcam",
    "GAPDH" = "Gapdh", "ACTB" = "Actb"
  )
  return(homolog_map)
}

# ==================== 步骤1: 加载数据 ====================
cat("\n========== 步骤1: 加载数据 ==========\n")

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
    
    # 去重
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
cat("✓ 原始数据加载完成\n\n")

# ==================== 步骤2: 质控 ====================
cat("\n========== 步骤2: 质控 ==========\n")
cat("质控前细胞数:", ncol(seurat_obj), "\n")

seurat_obj <- subset(seurat_obj, 
                     subset = nFeature_RNA > 200 & 
                              nFeature_RNA < 6000 & 
                              percent.mt < 20)

cat("质控后细胞数:", ncol(seurat_obj), "\n")

pdf(paste0(sample_name, "_QC.pdf"), width = 12, height = 6)
print(VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
              ncol = 3, pt.size = 0.1))
dev.off()
cat("✓ 质控图已保存\n\n")

# ==================== 步骤3: 标准化 ====================
cat("\n========== 步骤3: 标准化 ==========\n")
seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 2000, verbose = FALSE)
cat("高变基因数:", length(VariableFeatures(seurat_obj)), "\n")

pdf(paste0(sample_name, "_variable_features.pdf"), width = 10, height = 6)
top10 <- head(VariableFeatures(seurat_obj), 10)
plot1 <- VariableFeaturePlot(seurat_obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
print(plot2)
dev.off()
cat("✓ 高变基因图已保存\n")

seurat_obj <- ScaleData(seurat_obj, vars.to.regress = "percent.mt", verbose = FALSE)
cat("✓ 标准化完成\n\n")

# ==================== 步骤4: PCA降维 ====================
cat("\n========== 步骤4: PCA降维 ==========\n")
seurat_obj <- RunPCA(seurat_obj, npcs = 50, verbose = FALSE)

pdf(paste0(sample_name, "_elbow.pdf"), width = 8, height = 6)
print(ElbowPlot(seurat_obj, ndims = 50))
dev.off()
cat("✓ 肘部图已保存\n\n")

# ==================== 步骤5: 聚类 ====================
cat("\n========== 步骤5: 聚类 ==========\n")
pc_use <- 1:30
cat("使用前", length(pc_use), "个PC\n")

seurat_obj <- FindNeighbors(seurat_obj, dims = pc_use, verbose = FALSE)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.8, verbose = FALSE)
cat("找到", length(unique(Idents(seurat_obj))), "个cluster\n")

pdf(paste0(sample_name, "_clusters.pdf"), width = 10, height = 8)
print(DimPlot(seurat_obj, reduction = "pca", label = TRUE) + 
        ggtitle(paste0(sample_name, " PCA clusters")))
dev.off()
cat("✓ 聚类图已保存\n\n")

# ==================== 步骤6: UMAP降维 ====================
cat("\n========== 步骤6: UMAP降维 ==========\n")
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, verbose = FALSE)

pdf(paste0(sample_name, "_umap_clusters.pdf"), width = 10, height = 8)
print(DimPlot(seurat_obj, reduction = "umap", label = TRUE) + 
        ggtitle(paste0(sample_name, " UMAP clusters (", ncol(seurat_obj), " cells)")))
dev.off()
cat("✓ UMAP图已保存\n\n")

# ==================== 步骤7: SingleR自动注释 ====================
cat("\n========== 步骤7: SingleR自动注释 ==========\n")
expr_data <- GetAssayData(seurat_obj, layer = "data")

ref_immgen <- ImmGenData()
ref_mouse <- MouseRNAseqData()

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
cat("✓ SingleR注释完成\n\n")

# ==================== 步骤8: 手动精细注释 ====================
cat("\n========== 步骤8: 手动精细注释 ==========\n")

liver_markers <- create_liver_marker_db()
cat("定义了", length(liver_markers), "种细胞类型\n")

# 定义大类
major_names <- c(
  "Hepatocytes", "Stellate", "Endothelial", "Kupffer", 
  "T_cells", "B_cells", "Monocytes", "Neutrophils", 
  "Cholangiocytes", "Pericytes", "NK_cells",
  "Macrophage_M2", "Macrophage_NAM", "LSEC", "HSC_Activated",
  "Dendritic", "Erythrocytes", "Proliferating", "TECHO"
)

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

# 计算每个细胞的最高得分
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
if(sum(other_idx) > 0) {
  cat("用SingleR结果填补", sum(other_idx), "个Other细胞...\n")
  final_assignment[other_idx] <- as.character(seurat_obj$celltype[other_idx])
}

# 确保所有细胞都有注释
final_assignment[is.na(final_assignment)] <- "Unknown"

# 添加注释到对象
seurat_obj$final_celltype <- final_assignment
Idents(seurat_obj) <- seurat_obj$final_celltype

cat("\n最终细胞类型分布:\n")
celltype_counts <- table(Idents(seurat_obj))
print(celltype_counts)

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
p2 <- DimPlot(seurat_obj, reduction = "umap", group.by = "celltype", 
              label = TRUE, label.size = 3, repel = TRUE) +
  ggtitle(paste0(sample_name, " - SingleR Annotation"))
print(p1 | p2)
dev.off()
cat("✓ 注释对比图已保存\n")

# 保存对象
seurat_obj@misc$annotation_threshold <- threshold
saveRDS(seurat_obj, paste0(sample_name, "_manual_annotated.rds"))

# 验证保存
if(file.exists(paste0(sample_name, "_manual_annotated.rds"))) {
  file_size <- file.size(paste0(sample_name, "_manual_annotated.rds")) / 1024^3
  cat("\n✓ 手动注释结果已保存，大小:", round(file_size, 2), "GB\n")
} else {
  cat("\n✗ 保存失败！\n")
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
