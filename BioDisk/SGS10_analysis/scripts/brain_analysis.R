# ==================== 脑样本分析脚本 (纯鼠) ====================
# 文件名: ~/SGS10_analysis/scripts/brain_analysis.R
# 用法:   Rscript brain_analysis.R <sample_name> <n_cores>
# 说明:   全部脑样本均为纯鼠,无需 PDX 拆分。
#         保留亚型敏感分层共识注释 + 检查点机制。

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("用法: Rscript brain_analysis.R <sample_name> <n_cores>")
sample_name <- args[1]
n_cores     <- as.numeric(args[2])

options(echo = TRUE)
SCRIPT_DIR <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis/scripts"
source(file.path(SCRIPT_DIR, "pdx_utils.R"))

cat("========================================\n")
cat("脑样本分析:", sample_name, "\n")
cat("核心数:", n_cores, "\n")
cat("时间:", format(Sys.time(), "%F %T"), "\n")
cat("========================================\n\n")

suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(future)
  library(SingleR); library(celldex); library(patchwork); library(Matrix)
})

plan("multicore", workers = n_cores)
options(future.globals.maxSize = 50 * 1024^3)
options(future.rng.onMisuse = "ignore")
set.seed(42)

# ==================== 数据路径 ====================
data_paths <- list(
  nao56     = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/nao_56_filtered_feature_bc_matrix",
  brain3    = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/brain_samples/Brain_3_filtered_feature_bc_matrix",
  brain8    = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/brain_samples/Brain_8_filtered_feature_bc_matrix",
  bla_brain = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/brain_samples/Bla_brain_1_filtered_feature_bc_matrix",
  con_brain = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/brain_samples/Con_brain_1_filtered_feature_bc_matrix"
)

sample_dir <- file.path("/mnt/BioDisk/wangmingyuan/SGS10_analysis", sample_name, "mouse")
dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
setwd(sample_dir)
cat("工作目录:", getwd(), "\n\n")

# ==================== 检查点 ====================
checkpoint_file <- paste0(sample_name, "_checkpoint.txt")
completed_steps <- if (file.exists(checkpoint_file)) readLines(checkpoint_file, warn = FALSE) else character(0)
save_checkpoint <- function(step) {
  completed_steps <<- c(completed_steps, step)
  writeLines(completed_steps, checkpoint_file)
  cat("✓ 检查点:", step, "\n")
}
step_done <- function(step, file) {
  step %in% completed_steps && file.exists(file)
}

# ==================== 鼠脑 marker 库 ====================
create_brain_marker_db <- function() {
  list(
    Neurons              = c("Snap25","Syp","Syn1","Rbfox3","Tubb3","Map2"),
    Excitatory_L23       = c("Slc17a7","Cux1","Cux2","Rorb","Satb2"),
    Excitatory_L45       = c("Slc17a7","Fezf2","Bcl11b","Etv1","Tle4"),
    Excitatory_L6        = c("Slc17a7","Tle4","Foxp2","Sulf1"),
    Hpc_CA1              = c("Slc17a7","Neurod6","Camk2a","Calb1"),
    Hpc_CA3              = c("Slc17a7","Neurod6","Camk2a","Calb2"),
    Hpc_DG               = c("Prox1","Dcx","Calb1","Pcp4"),
    Inhibitory_Sst       = c("Gad1","Gad2","Sst","Cort"),
    Inhibitory_Pvalb     = c("Gad1","Gad2","Pvalb"),
    Inhibitory_Vip       = c("Gad1","Gad2","Vip","Cck"),
    Inhibitory_Meis2     = c("Gad1","Gad2","Meis2","Lhx6"),
    Astrocytes           = c("Slc1a3","Gfap","Aqp4","Aldh1l1","S100b"),
    Microglia            = c("Tmem119","P2ry12","Cx3cr1","Hexb","Csf1r","Cd68"),
    Microglia_Homeostatic= c("P2ry12","Tmem119","Cx3cr1","Siglech"),
    Microglia_DAM        = c("Apoe","Clec7a","Lpl","Cst7"),
    Microglia_M1         = c("Il1b","Tnf","Cd86","Nos2"),
    Microglia_M2         = c("Arg1","Cd163","Mrc1","Il10"),
    Oligodendrocytes     = c("Mbp","Mog","Plp1","Cnp","Mobp"),
    OPCs                 = c("Pdgfra","Cspg4","Olig1","Olig2"),
    Endothelial          = c("Cldn5","Pecam1","Cdh5","Vwf","Flt1"),
    Pericytes            = c("Rgs5","Pdgfrb","Cspg4","Kcnj8","Acta2"),
    VLMC                 = c("Lum","Dcn","Col1a1","Col3a1"),
    Ependymal            = c("Foxj1","Vim","Ccdc153","Cfap126"),
    CPECs                = c("Ttr","Lcn2","Clcn2","Cldn3","Slc2a1"),
    T_cells              = c("Cd3d","Cd3e","Cd3g","Cd4","Cd8a"),
    B_cells              = c("Cd19","Ms4a1","Cd79a","Cd79b"),
    NK_cells             = c("Nkg7","Klrb1c","Klrd1")
  )
}

SUBTYPES <- c("Excitatory_L23","Excitatory_L45","Excitatory_L6",
              "Hpc_CA1","Hpc_CA3","Hpc_DG",
              "Inhibitory_Sst","Inhibitory_Pvalb","Inhibitory_Vip","Inhibitory_Meis2",
              "Microglia_Homeostatic","Microglia_DAM","Microglia_M1","Microglia_M2",
              "OPCs")
MAJORS   <- c("Neurons","Astrocytes","Microglia","Oligodendrocytes",
              "Endothelial","Pericytes","VLMC","Ependymal","CPECs",
              "T_cells","B_cells","NK_cells")
SUBTYPE2MAJOR <- c(
  Excitatory_L23="Neurons", Excitatory_L45="Neurons", Excitatory_L6="Neurons",
  Hpc_CA1="Neurons", Hpc_CA3="Neurons", Hpc_DG="Neurons",
  Inhibitory_Sst="Neurons", Inhibitory_Pvalb="Neurons",
  Inhibitory_Vip="Neurons", Inhibitory_Meis2="Neurons",
  Microglia_Homeostatic="Microglia", Microglia_DAM="Microglia",
  Microglia_M1="Microglia", Microglia_M2="Microglia",
  OPCs="Oligodendrocytes"
)

# ==================== 步骤1: 加载并清理数据 ====================
raw_rds <- paste0(sample_name, "_raw.rds")
if (step_done("data_loaded", raw_rds)) {
  cat("加载已有原始数据...\n")
  seurat_obj <- readRDS(raw_rds)
} else {
  data_dir <- data_paths[[sample_name]]
  if (is.null(data_dir)) stop("未知样本: ", sample_name)

  # 用统一加载函数,自动剥离 mm10- 等前缀
  split_data <- load_and_split_pdx(data_dir, sample_name)
  if (is.null(split_data$mouse)) {
    stop("脑样本 ", sample_name, " 未能获得鼠源矩阵,请检查参考基因组")
  }
  if (split_data$is_mixed_ref && !is.null(split_data$human)) {
    n_h <- ncol(split_data$human); n_m <- ncol(split_data$mouse)
    cat("⚠️ 注意:脑样本检测到混合参考前缀,但按需求只取鼠源矩阵\n")
  }
  counts <- split_data$mouse

  seurat_obj <- CreateSeuratObject(counts = counts, project = sample_name)
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
  saveRDS(seurat_obj, raw_rds)
  save_checkpoint("data_loaded")
}

# ==================== 步骤2: QC ====================
qc_rds <- paste0(sample_name, "_qc.rds")
if (step_done("qc_done", qc_rds)) {
  seurat_obj <- readRDS(qc_rds)
} else {
  cat("QC 前细胞数:", ncol(seurat_obj), "\n")
  seurat_obj <- subset(seurat_obj,
                       subset = nFeature_RNA > 200 &
                                nFeature_RNA < 6000 &
                                percent.mt < 20)
  cat("QC 后细胞数:", ncol(seurat_obj), "\n")
  pdf(paste0(sample_name, "_QC.pdf"), width = 12, height = 5)
  print(VlnPlot(seurat_obj, features = c("nFeature_RNA","nCount_RNA","percent.mt"),
                ncol = 3, pt.size = 0.1))
  dev.off()
  saveRDS(seurat_obj, qc_rds)
  save_checkpoint("qc_done")
}

# ==================== 步骤3: 标准化 + Scale ====================
norm_rds <- paste0(sample_name, "_norm.rds")
if (step_done("norm_done", norm_rds)) {
  seurat_obj <- readRDS(norm_rds)
} else {
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 2000, verbose = FALSE)
  seurat_obj <- ScaleData(seurat_obj, vars.to.regress = "percent.mt", verbose = FALSE)
  saveRDS(seurat_obj, norm_rds)
  save_checkpoint("norm_done")
}

# ==================== 步骤4-6: PCA / Cluster / UMAP ====================
umap_rds <- paste0(sample_name, "_umap.rds")
if (step_done("umap_done", umap_rds)) {
  seurat_obj <- readRDS(umap_rds)
} else {
  seurat_obj <- RunPCA(seurat_obj, npcs = 50, verbose = FALSE)
  pdf(paste0(sample_name, "_elbow.pdf"), width = 7, height = 5)
  print(ElbowPlot(seurat_obj, ndims = 50)); dev.off()

  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30, verbose = FALSE)
  seurat_obj <- FindClusters(seurat_obj,  resolution = 0.8, verbose = FALSE)
  seurat_obj <- RunUMAP(seurat_obj,       dims = 1:30, verbose = FALSE)

  pdf(paste0(sample_name, "_clusters.pdf"), width = 9, height = 7)
  print(DimPlot(seurat_obj, reduction = "umap", label = TRUE) +
          ggtitle(paste0(sample_name, " clusters (n=", ncol(seurat_obj), ")")))
  dev.off()
  saveRDS(seurat_obj, umap_rds)
  save_checkpoint("umap_done")
}

# ==================== 步骤7: SingleR ====================
singler_rds <- paste0(sample_name, "_singler.rds")
if (step_done("singler_done", singler_rds)) {
  seurat_obj <- readRDS(singler_rds)
} else {
  cat("SingleR 注释...\n")
  expr <- GetAssayData(seurat_obj, layer = "data")
  ref1 <- celldex::ImmGenData()
  ref2 <- celldex::MouseRNAseqData()
  sr <- SingleR(test = expr,
                ref = list(ImmGen = ref1, Mouse = ref2),
                labels = list(ref1$label.main, ref2$label.main),
                clusters = Idents(seurat_obj),
                assay.type.test = "logcounts",
                fine.tune = TRUE)
  seurat_obj$singler_label <- sr$labels[as.character(Idents(seurat_obj))]

  label_map <- c("Microglia"="Microglia","Macrophages"="Microglia","Monocytes"="Microglia",
                 "Astrocytes"="Astrocytes","Neurons"="Neurons",
                 "Oligodendrocytes"="Oligodendrocytes","OPC"="OPCs",
                 "Endothelial cells"="Endothelial","Pericytes"="Pericytes",
                 "Ependymal cells"="Ependymal","Epithelial cells"="CPECs")
  ct <- seurat_obj$singler_label
  for (k in names(label_map)) ct[ct == k] <- label_map[k]
  ct[is.na(ct)] <- "Other"
  seurat_obj$celltype <- ct
  saveRDS(seurat_obj, singler_rds)
  save_checkpoint("singler_done")
}

# ==================== 步骤8: 亚型敏感共识注释 ====================
final_rds <- paste0(sample_name, "_manual_annotated.rds")
if (step_done("manual_annotation_done", final_rds)) {
  seurat_obj <- readRDS(final_rds)
} else {
  cat("\n亚型敏感共识注释...\n")
  brain_markers <- create_brain_marker_db()
  expr <- GetAssayData(seurat_obj, layer = "data")

  all_scores <- matrix(0, nrow = ncol(seurat_obj), ncol = length(brain_markers),
                       dimnames = list(colnames(seurat_obj), names(brain_markers)))
  for (ct in names(brain_markers)) {
    g <- intersect(brain_markers[[ct]], rownames(expr))
    if (length(g)) all_scores[, ct] <- Matrix::colMeans(expr[g, , drop = FALSE])
  }

  sub_scores <- all_scores[, intersect(SUBTYPES, colnames(all_scores)), drop = FALSE]
  maj_scores <- all_scores[, intersect(MAJORS,   colnames(all_scores)), drop = FALSE]

  pos <- all_scores[all_scores > 0]
  sub_thresh <- max(quantile(pos, 0.05, na.rm = TRUE), 0.03)
  maj_thresh <- max(quantile(pos, 0.02, na.rm = TRUE), 0.02)
  cat("  亚型阈值:", round(sub_thresh, 4), "  大类阈值:", round(maj_thresh, 4), "\n")

  # 每个细胞先找最强亚型,不够再找大类
  best_sub_idx <- max.col(sub_scores, ties.method = "first")
  best_sub     <- colnames(sub_scores)[best_sub_idx]
  best_sub_val <- sub_scores[cbind(seq_len(nrow(sub_scores)), best_sub_idx)]

  best_maj_idx <- max.col(maj_scores, ties.method = "first")
  best_maj     <- colnames(maj_scores)[best_maj_idx]
  best_maj_val <- maj_scores[cbind(seq_len(nrow(maj_scores)), best_maj_idx)]

  assignment <- rep(NA_character_, ncol(seurat_obj))
  assignment[best_sub_val >= sub_thresh] <- best_sub[best_sub_val >= sub_thresh]
  fill_maj <- is.na(assignment) & (best_maj_val >= maj_thresh)
  assignment[fill_maj] <- best_maj[fill_maj]

  # 基于 cluster 的共识 + 亚型保护
  cluster_ids <- as.character(Idents(seurat_obj))
  consensus <- assignment
  for (cl in unique(cluster_ids)) {
    in_cl <- which(cluster_ids == cl)
    if (length(in_cl) < 5) next
    a <- assignment[in_cl]
    a_clean <- a[!is.na(a)]
    if (length(a_clean) < 3) next
    sub_a <- a_clean[a_clean %in% SUBTYPES]
    maj_a <- a_clean[a_clean %in% MAJORS]
    if (length(sub_a) >= 0.5 * length(a_clean) && length(unique(sub_a)) == 1) {
      consensus[in_cl] <- unique(sub_a)
    } else if (length(maj_a) >= 3) {
      tt <- table(maj_a)
      if (max(tt) / sum(tt) > 0.7) consensus[in_cl] <- names(tt)[which.max(tt)]
    }
  }
  # 亚型保护:如果原本是亚型而 consensus 把它压成对应大类,保留亚型
  for (i in seq_along(assignment)) {
    if (!is.na(assignment[i]) && assignment[i] %in% names(SUBTYPE2MAJOR)) {
      if (!is.na(consensus[i]) && consensus[i] == SUBTYPE2MAJOR[assignment[i]]) {
        consensus[i] <- assignment[i]
      }
    }
  }
  na_idx <- is.na(consensus)
  if (any(na_idx)) consensus[na_idx] <- seurat_obj$celltype[na_idx]
  consensus[is.na(consensus)] <- "Unknown"

  seurat_obj$final_celltype <- consensus
  Idents(seurat_obj) <- seurat_obj$final_celltype

  cat("\n最终细胞类型分布:\n"); print(table(Idents(seurat_obj)))

  pdf(paste0(sample_name, "_celltypes_umap.pdf"), width = 12, height = 8)
  print(DimPlot(seurat_obj, reduction = "umap", group.by = "final_celltype",
                label = TRUE, repel = TRUE, label.size = 3) +
          ggtitle(paste0(sample_name, " - Manual annotation")))
  dev.off()

  pdf(paste0(sample_name, "_annotation_compare.pdf"), width = 16, height = 7)
  print(
    (DimPlot(seurat_obj, group.by = "final_celltype", label = TRUE, repel = TRUE) + ggtitle("Manual")) |
    (DimPlot(seurat_obj, group.by = "celltype",       label = TRUE, repel = TRUE) + ggtitle("SingleR"))
  )
  dev.off()

  ct_tab <- table(seurat_obj$final_celltype)
  prop_df <- data.frame(CellType = names(ct_tab),
                        Count = as.numeric(ct_tab),
                        Percentage = round(as.numeric(prop.table(ct_tab)) * 100, 2))
  write.csv(prop_df, paste0(sample_name, "_proportions.csv"), row.names = FALSE)

  seurat_obj@misc$sub_thresh <- sub_thresh
  seurat_obj@misc$maj_thresh <- maj_thresh
  saveRDS(seurat_obj, final_rds)
  save_checkpoint("manual_annotation_done")
}

# ==================== 报告 ====================
sink(paste0(sample_name, "_report.txt"))
cat("Sample:", sample_name, "\n")
cat("Time:", format(Sys.time(), "%F %T"), "\n")
cat("Cells:", ncol(seurat_obj), "  Genes:", nrow(seurat_obj), "\n\n")
cat("Cell type distribution:\n")
ct_tab <- table(seurat_obj$final_celltype)
for (i in seq_along(ct_tab)) {
  cat(sprintf("  %-22s %6d (%5.1f%%)\n",
              names(ct_tab)[i], ct_tab[i], 100*ct_tab[i]/sum(ct_tab)))
}
sink()

plan("sequential")
cat("\n========== ", sample_name, "脑分析完成 ==========\n")
cat("结束时间:", format(Sys.time(), "%F %T"), "\n")
