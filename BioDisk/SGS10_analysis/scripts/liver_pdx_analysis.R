# ==================== 肝脏样本分析脚本 v3 ====================
# 文件名: ~/SGS10_analysis/scripts/liver_pdx_analysis.R
# 用法:   Rscript liver_pdx_analysis.R <sample_name> <n_cores>
#
# v3 关键改动 (vs v2):
#   用 Seurat 的 AddModuleScore 替代 colMeans 做 marker 评分
#   AddModuleScore 会减去背景基因集的平均表达,能正确识别真正富集 marker 的细胞
#   这样不会被 ALB/TTR 这样的超高表达泛肝细胞 marker 主导

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("用法: Rscript liver_pdx_analysis.R <sample_name> <n_cores>")
sample_name <- args[1]
n_cores     <- as.numeric(args[2])

options(echo = TRUE)
SCRIPT_DIR <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis/scripts"
source(file.path(SCRIPT_DIR, "pdx_utils.R"))

cat("========================================\n")
cat("肝脏 PDX 分析 v3:", sample_name, "\n")
cat("核心数:", n_cores, "\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================\n\n")

suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(future)
  library(SingleR); library(celldex); library(patchwork); library(Matrix)
})

plan("multicore", workers = n_cores)
options(future.globals.maxSize = 50 * 1024^3)
options(future.rng.onMisuse = "ignore")
set.seed(42)

data_paths <- list(
  gan56     = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/gan56/gan_56_filtered_feature_bc_matrix",
  liver3    = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/liver_samples/25090559_Liver_3_expression/Liver_3_filtered_feature_bc_matrix",
  liver8    = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/liver_samples/25090559_Liver_8_expression/Liver_8_filtered_feature_bc_matrix",
  bla_liver = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/liver_samples/25090559_Bla_liver_1_expression/Bla_liver_1_filtered_feature_bc_matrix",
  con_liver = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/liver_samples/25090559_Con_liver_1_expression/Con_liver_1_filtered_feature_bc_matrix"
)
PDX_SAMPLES <- c("gan56", "liver3", "liver8")
is_pdx <- sample_name %in% PDX_SAMPLES

sample_dir <- file.path("/mnt/BioDisk/wangmingyuan/SGS10_analysis", sample_name)
dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
setwd(sample_dir)
cat("工作目录:", getwd(), "\n")
cat("是否 PDX 样本:", is_pdx, "\n\n")

# ==================== 鼠肝 marker 库 (去掉 Hepatocytes 总类) ====================
create_mouse_liver_markers <- function() {
  list(
    # 注意: 不再保留 "Hepatocytes" 泛肝细胞总类
    # 因为 Alb/Ttr 表达量太高会压制所有其他类别
    # 让分区 marker 自己区分真实肝细胞

    # Halpern 2017 分区 marker (加入少量泛肝标志以确保识别为肝细胞)
    Hep_Periportal     = c("Hal","Sds","Cyp2f2","Cdh1","Pck1","G6pc","Cps1","Sdhd","Arg1","Mup3"),
    Hep_Perivenous     = c("Glul","Cyp2e1","Cyp1a2","Oat","Axin2","Lgr5","Notum","Lhpp","Slc1a2","Cyp7a1"),
    Hep_Proliferating  = c("Mki67","Top2a","Pcna","Ccnb1","Mcm2","Ccnd1"),

    # HSC 亚型
    Stellate           = c("Pdgfrb","Des","Rgs5","Lrat","Reln","Colec11"),
    HSC_Quiescent      = c("Lrat","Angptl6","Vipr1","Rgs5","Reln","Pparg","Rbp1"),
    HSC_Activated_ECM  = c("Acta2","Col1a1","Col1a2","Col3a1","Timp1","Timp2","Aebp1","Bgn"),
    HSC_Activated_Chemo= c("Ccl2","Ccl6","Cxcl10","Ccl24","Ccr2"),
    HSC_Proliferating  = c("Mki67","Top2a","Mcm6","Ccnb1"),

    # 内皮分区 (Halpern 2018)
    Endothelial        = c("Pecam1","Cdh5","Vwf","Kdr","Flt1","Tek","Cldn5"),
    LSEC               = c("Clec4g","Lyve1","Stab1","Stab2","Fcgr2b","Oit3","Gpr182"),
    LSEC_Periportal    = c("Clec4g","Ednrb","Jag1","Lrg1","Efnb1","Ltbp4","Adgrg6"),
    LSEC_Pericentral   = c("Clec4g","Wnt9b","Rspo3","Cdh13","Wnt2","Kit","Thbd"),
    LSEC_Mid           = c("Clec4g","Lyve1","Ctsl","Hamp","Igfbp2"),

    # 巨噬细胞 (NAM 来自 Xiong 2019)
    Kupffer            = c("Clec4f","Vsig4","Timd4","Marco","C1qa","C1qb","C1qc","Cd68","Cd163","Adgre1"),
    KC_Homeostatic     = c("Clec4f","Vsig4","Timd4","Marco","Cdh5"),
    Macrophage_NAM     = c("Trem2","Gpnmb","Cd9","Spp1","Cd63","Lpl","Cstb","Fabp5"),
    Macrophage_M1      = c("Il1b","Tnf","Cd86","Nos2","Cxcl9","Cxcl10"),
    Macrophage_M2      = c("Arg1","Mrc1","Cd163","Il10","Retnla","Chil3"),
    Monocytes          = c("Cd14","Ly6c2","Ccr2","S100a8","S100a9","Fcn1"),
    LCM                = c("Cx3cr1","Itgax","S100a4"),

    Cholangiocytes     = c("Krt7","Krt8","Krt18","Krt19","Epcam","Sox9","Cftr","Spp1"),

    T_cells            = c("Cd3d","Cd3e","Cd3g"),
    CD4_T              = c("Cd4","Cd3e","Lef1"),
    CD8_T              = c("Cd8a","Cd8b1","Cd3e","Gzmk","Prf1"),
    Treg               = c("Foxp3","Il2ra","Ctla4"),
    NK_cells           = c("Nkg7","Klrb1c","Klrd1","Klrk1","Ncr1","Gnly"),
    NKT                = c("Cd3e","Nkg7","Klrb1c"),
    B_cells            = c("Cd19","Ms4a1","Cd79a","Cd79b","Ighd","Fcmr"),
    B_Plasma           = c("Jchain","Xbp1","Irf4","Mzb1","Igha","Sdc1"),

    Dendritic          = c("Flt3","Zbtb46","Itgax"),
    cDC1               = c("Clec9a","Xcr1","Batf3","Irf8"),
    cDC2               = c("Cd1c","Sirpa","Clec10a","Irf4"),
    pDC                = c("Siglech","Bst2","Tcf4","Irf7"),

    Neutrophils        = c("S100a8","S100a9","Cxcl2","Cxcr2","Mmp8","Csf3r","Ly6g"),
    Pericytes          = c("Rgs5","Pdgfrb","Cspg4","Kcnj8","Vtn","Higd1b"),
    Erythrocytes       = c("Hba-a1","Hba-a2","Hbb-bt","Hbb-bs","Klf1","Gypa")
  )
}

create_human_liver_markers <- function() {
  list(
    Tumor_Hepatocyte       = c("ALB","AFP","GPC3","HNF4A","HNF1A","TTR"),
    Tumor_Proliferating    = c("MKI67","TOP2A","PCNA","CCNB1","CDK1","MCM2"),
    Tumor_EMT              = c("VIM","SNAI1","ZEB1","CDH2","S100A4"),
    Tumor_Stem             = c("EPCAM","PROM1","CD44","KRT19","SOX9"),
    HCC_GPC3hi             = c("GPC3","AFP","ALB","SERPINA1"),

    Hep_Periportal         = c("HAL","SDS","CYP2F2","ARG1","PCK1","G6PC","CPS1"),
    Hep_Perivenous         = c("GLUL","CYP2E1","CYP1A2","OAT","AXIN2","LGR5"),

    Cholangiocytes         = c("KRT7","KRT19","EPCAM","SOX9","CFTR","SPP1"),

    Stellate               = c("ACTA2","COL1A1","COL1A2","PDGFRB","DES","LRAT","RGS5"),
    HSC_Activated          = c("ACTA2","COL1A1","TIMP1","TGFB1"),
    Endothelial            = c("PECAM1","CDH5","VWF","KDR","CLDN5"),
    LSEC                   = c("CLEC4G","LYVE1","STAB2","FCN3"),

    Kupffer_Macrophage     = c("CD68","CD163","VSIG4","C1QA","C1QB","MARCO","CD5L"),
    Macrophage_NAM         = c("TREM2","GPNMB","CD9","SPP1","CD63"),
    Monocytes              = c("CD14","LYZ","S100A8","S100A9","FCN1","VCAN"),
    Dendritic              = c("FLT3","CLEC9A","ITGAX","HLA-DRA","CD1C"),

    T_cells                = c("CD3D","CD3E","CD3G"),
    CD4_T                  = c("CD4","CD3E","IL7R"),
    CD8_T                  = c("CD8A","CD8B","CD3E","GZMK"),
    Treg                   = c("FOXP3","IL2RA","CTLA4"),
    NK_cells               = c("NKG7","KLRD1","KLRB1","GNLY","NCAM1","NCR1"),
    B_cells                = c("CD19","MS4A1","CD79A","CD79B"),
    Plasma                 = c("JCHAIN","XBP1","IRF4","MZB1","SDC1"),

    Neutrophils            = c("S100A8","S100A9","CXCL8","CXCR2","FCGR3B","CSF3R"),
    Erythrocytes           = c("HBA1","HBA2","HBB","KLF1")
  )
}

# ==================== 关键改进: 用 AddModuleScore 替代 colMeans ====================
# AddModuleScore 原理:
#   对每个 marker 基因集,计算每个细胞的平均表达
#   然后减去"背景"基因集(随机选取的表达水平匹配的基因)的平均表达
#   这样能正确识别"真正富集该 marker 组"的细胞,不被 Alb/Ttr 这类超高表达基因主导

score_markers_module <- function(obj, markers_db) {
  cat("  使用 AddModuleScore 评分 (", length(markers_db), "种类型)...\n")
  score_mat <- matrix(NA_real_, nrow = ncol(obj), ncol = length(markers_db),
                      dimnames = list(colnames(obj), names(markers_db)))

  for (ct in names(markers_db)) {
    g <- intersect(markers_db[[ct]], rownames(obj))
    if (length(g) >= 2) {
      obj <- tryCatch(
        AddModuleScore(obj, features = list(g), name = "tmpscore",
                       nbin = 24, ctrl = 100, seed = 42),
        error = function(e) {
          cat("    ", ct, ": AddModuleScore 失败 (", e$message, "), fallback 到均值\n")
          obj
        }
      )
      if ("tmpscore1" %in% colnames(obj@meta.data)) {
        score_mat[, ct] <- obj$tmpscore1
        obj$tmpscore1 <- NULL
      } else {
        # fallback: 用均值
        expr <- GetAssayData(obj, layer = "data")
        score_mat[, ct] <- Matrix::colMeans(expr[g, , drop = FALSE])
      }
    } else if (length(g) == 1) {
      expr <- GetAssayData(obj, layer = "data")
      score_mat[, ct] <- as.numeric(expr[g, ])
    } else {
      # 0 个 marker 匹配,填 -Inf 避免被选中
      score_mat[, ct] <- -Inf
    }
  }
  score_mat
}

# ==================== 通用单物种分析流程 ====================
run_single_species_pipeline <- function(counts_mat, sample_name, species, out_dir,
                                        n_cores = 4) {
  cat("\n========================================\n")
  cat("单物种流程:", species, "\n")
  cat("矩阵维度:", nrow(counts_mat), "基因 ×", ncol(counts_mat), "细胞\n")
  cat("输出目录:", out_dir, "\n")
  cat("========================================\n")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  setwd(out_dir)

  prefix <- paste0(sample_name, "_", species)

  obj <- CreateSeuratObject(counts = counts_mat, project = prefix)
  mt_pat <- if (species == "mouse") "^mt-" else "^MT-"
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mt_pat)
  cat("初始细胞数:", ncol(obj), " 平均 mt%:", round(mean(obj$percent.mt), 2), "\n")

  obj <- subset(obj,
                subset = nFeature_RNA > 200 &
                         nFeature_RNA < 6000 &
                         percent.mt   < 20)
  cat("QC 后细胞数:", ncol(obj), "\n")
  if (ncol(obj) < 50) {
    cat("⚠️ 细胞数过少 (<50),跳过此物种\n")
    return(NULL)
  }

  pdf(paste0(prefix, "_QC.pdf"), width = 12, height = 5)
  print(VlnPlot(obj, features = c("nFeature_RNA","nCount_RNA","percent.mt"),
                ncol = 3, pt.size = 0.1))
  dev.off()

  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
  obj <- ScaleData(obj, vars.to.regress = "percent.mt", verbose = FALSE)

  npcs <- min(50, ncol(obj) - 1)
  obj <- RunPCA(obj, npcs = npcs, verbose = FALSE)
  pdf(paste0(prefix, "_elbow.pdf"), width = 7, height = 5)
  print(ElbowPlot(obj, ndims = npcs))
  dev.off()

  use_dims <- 1:min(30, npcs)
  obj <- FindNeighbors(obj, dims = use_dims, verbose = FALSE)
  obj <- FindClusters(obj,  resolution = 0.8, verbose = FALSE)
  obj <- RunUMAP(obj,       dims = use_dims, verbose = FALSE)

  pdf(paste0(prefix, "_umap_clusters.pdf"), width = 9, height = 7)
  print(DimPlot(obj, reduction = "umap", label = TRUE) +
          ggtitle(paste0(prefix, " clusters (n=", ncol(obj), ")")))
  dev.off()

  # ----- SingleR (容错) -----
  cat("\nSingleR 注释 (", species, ")...\n")
  expr <- GetAssayData(obj, layer = "data")
  singler_labels <- tryCatch({
    if (species == "mouse") {
      ref1 <- celldex::ImmGenData()
      ref2 <- celldex::MouseRNAseqData()
      sr <- SingleR(test = expr,
                    ref = list(ImmGen = ref1, Mouse = ref2),
                    labels = list(ref1$label.main, ref2$label.main),
                    clusters = Idents(obj),
                    assay.type.test = "logcounts",
                    fine.tune = TRUE)
    } else {
      ref1 <- celldex::HumanPrimaryCellAtlasData()
      sr <- SingleR(test = expr,
                    ref = ref1,
                    labels = ref1$label.main,
                    clusters = Idents(obj),
                    assay.type.test = "logcounts",
                    fine.tune = TRUE)
    }
    sr$labels[as.character(Idents(obj))]
  }, error = function(e) {
    cat("  ⚠️ SingleR 失败:", e$message, "\n")
    cat("  ℹ️ 仅依赖 marker 评分注释\n")
    rep("Unknown", ncol(obj))
  })
  obj$singler_label <- singler_labels

  # ----- Marker 评分 (v3: AddModuleScore) -----
  cat("\nMarker 评分注释 (v3: AddModuleScore)...\n")
  markers_db <- if (species == "mouse") create_mouse_liver_markers() else create_human_liver_markers()
  score_mat <- score_markers_module(obj, markers_db)

  # AddModuleScore 评分可正可负,阈值设 0 (表示"高于背景")
  best_idx  <- max.col(score_mat, ties.method = "first")
  best_type <- colnames(score_mat)[best_idx]
  best_val  <- score_mat[cbind(seq_len(nrow(score_mat)), best_idx)]
  # 阈值: 必须高于背景基因集均值(0 是中位水平)
  best_type[!is.finite(best_val) | best_val < 0] <- NA_character_

  # SingleR 兜底
  fill_idx <- is.na(best_type)
  best_type[fill_idx] <- obj$singler_label[fill_idx]
  best_type[is.na(best_type)] <- "Unknown"

  obj$final_celltype <- best_type
  Idents(obj) <- obj$final_celltype

  cat("\n最终细胞类型分布:\n"); print(sort(table(Idents(obj)), decreasing = TRUE))

  pdf(paste0(prefix, "_celltypes_umap.pdf"), width = 12, height = 8)
  print(DimPlot(obj, reduction = "umap", group.by = "final_celltype",
                label = TRUE, repel = TRUE, label.size = 3) +
          ggtitle(paste0(prefix, " - Manual annotation v3")))
  dev.off()

  pdf(paste0(prefix, "_annotation_compare.pdf"), width = 16, height = 7)
  p1 <- DimPlot(obj, group.by = "final_celltype", label = TRUE, repel = TRUE) +
          ggtitle("Manual (v3 AddModuleScore)")
  p2 <- DimPlot(obj, group.by = "singler_label", label = TRUE, repel = TRUE) +
          ggtitle("SingleR")
  print(p1 | p2)
  dev.off()

  ct_tab <- table(obj$final_celltype)
  prop_df <- data.frame(CellType = names(ct_tab),
                        Count = as.numeric(ct_tab),
                        Percentage = round(as.numeric(prop.table(ct_tab)) * 100, 2))
  prop_df <- prop_df[order(-prop_df$Count), ]
  write.csv(prop_df, paste0(prefix, "_proportions.csv"), row.names = FALSE)

  saveRDS(obj, paste0(prefix, "_annotated.rds"))
  cat("✓ 已保存:", paste0(prefix, "_annotated.rds"), "\n")

  sink(paste0(prefix, "_report.txt"))
  cat("Sample:", prefix, "\n")
  cat("Species:", species, "\n")
  cat("Time:", format(Sys.time(), "%F %T"), "\n")
  cat("Cells:", ncol(obj), "  Genes:", nrow(obj), "\n")
  cat("Marker library: v3.2 (AddModuleScore, Halpern2017+Xiong2019+Bai2023)\n\n")
  cat("Cell type distribution (sorted by count):\n")
  for (i in seq_len(nrow(prop_df))) {
    cat(sprintf("  %-25s %6d (%5.1f%%)\n",
                prop_df$CellType[i], prop_df$Count[i], prop_df$Percentage[i]))
  }
  sink()

  obj
}

# ==================== 主流程 ====================
split_data <- load_and_split_pdx(data_paths[[sample_name]], sample_name)
# === 去除人鼠前缀 (CellRanger 混合参考: GRCh38_ 和 mm10___) ===
strip_prefix <- function(mat, prefix_re) {
  if (is.null(mat)) return(NULL)
  rn <- rownames(mat)
  new_rn <- sub(prefix_re, "", rn)
  rownames(mat) <- make.unique(new_rn)
  mat
}
if (!is.null(split_data$human)) {
  cat("\n清理人源基因前缀...\n")
  cat("  示例(原):", head(rownames(split_data$human), 3), "\n")
  split_data$human <- strip_prefix(split_data$human, "^GRCh38_")
  cat("  示例(新):", head(rownames(split_data$human), 3), "\n")
}
if (!is.null(split_data$mouse)) {
  cat("\n清理鼠源基因前缀...\n")
  cat("  示例(原):", head(rownames(split_data$mouse), 3), "\n")
  split_data$mouse <- strip_prefix(split_data$mouse, "^_+")
  cat("  示例(新):", head(rownames(split_data$mouse), 3), "\n")
}
if (is_pdx) {
  if (!split_data$is_mixed_ref) {
    warning("样本 ", sample_name, " 标记为 PDX 但未检测到混合参考基因前缀")
  }

  if (!is.null(split_data$human) && !is.null(split_data$mouse)) {
    class_df <- classify_barcodes_by_species(split_data$human, split_data$mouse,
                                             pure_threshold = 0.90,
                                             min_total_umi  = 200)
    write.csv(class_df, file.path(sample_dir, paste0(sample_name, "_species_classification.csv")),
              row.names = FALSE)
    plot_barnyard(class_df, sample_name,
                  file.path(sample_dir, paste0(sample_name, "_barnyard.pdf")))

    human_bc <- class_df$barcode[class_df$species == "Human"]
    mouse_bc <- class_df$barcode[class_df$species == "Mouse"]
    cat("\n纯人细胞:", length(human_bc), " 纯鼠细胞:", length(mouse_bc), "\n")

    if (length(human_bc) >= 50) {
      h_mat <- split_data$human[, human_bc, drop = FALSE]
      run_single_species_pipeline(h_mat, sample_name, "human",
                                  out_dir = file.path(sample_dir, "human"),
                                  n_cores = n_cores)
    } else {
      cat("⚠️ 人细胞数 <50,跳过人源分析\n")
    }

    if (length(mouse_bc) >= 50) {
      m_mat <- split_data$mouse[, mouse_bc, drop = FALSE]
      run_single_species_pipeline(m_mat, sample_name, "mouse",
                                  out_dir = file.path(sample_dir, "mouse"),
                                  n_cores = n_cores)
    } else {
      cat("⚠️ 鼠细胞数 <50,跳过鼠源分析\n")
    }
  } else {
    cat("⚠️ PDX 样本但无法获得双物种矩阵\n")
  }

} else {
  cat("\n纯鼠对照样本流程\n")
  m_mat <- if (!is.null(split_data$mouse)) split_data$mouse else split_data$human
  if (is.null(m_mat)) stop("无法读取计数矩阵")
  run_single_species_pipeline(m_mat, sample_name, "mouse",
                              out_dir = file.path(sample_dir, "mouse"),
                              n_cores = n_cores)
}

plan("sequential")
cat("\n========== ", sample_name, "肝脏分析 v3 完成 ==========\n")
cat("结束时间:", format(Sys.time(), "%F %T"), "\n")
