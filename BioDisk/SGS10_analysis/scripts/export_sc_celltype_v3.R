#!/usr/bin/env Rscript
# v3: 修复 mito 计算 + 新增 top30 marker 表导出

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

PROJECT_ROOT <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
setwd(PROJECT_ROOT)
source("scripts/diff_common.R")

OUT_DIR <- file.path(PROJECT_ROOT, "sc_celltype_export_v3")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

MARKER_DIR <- file.path(OUT_DIR, "top30_markers")
dir.create(MARKER_DIR, showWarnings = FALSE, recursive = TRUE)

RAW_SAMPLE_MAP <- list(
  liver3    = "25090559_Liver_3",
  liver8    = "25090559_Liver_8",
  bla_liver = "25090559_Bla_liver_1",
  con_liver = "25090559_Con_liver_1",
  brain3    = "25090559_Brain_3",
  brain8    = "25090559_Brain_8",
  bla_brain = "25090559_Bla_brain_1",
  con_brain = "25090559_Con_brain_1",
  gan56     = "24070264_gan_56",
  nao56     = "24070263_nao_56"
)

# ============================================================
# mito calc (修复: 返回 numeric vector 而不是索引)
# ============================================================
calc_mito <- function(obj) {
  genes <- rownames(obj)
  # mouse: ^mt-, human: ^MT-
  mito_genes <- grep("^(mt-|MT-)", genes, value = TRUE, ignore.case = FALSE)
  if (length(mito_genes) == 0) {
    return(list(pct = rep(0, ncol(obj)), n_genes = 0))
  }
  pct_result <- PercentageFeatureSet(obj, features = mito_genes)
  # 兼容 v4/v5: 结果可能是 data.frame 或 matrix 或 numeric
  if (is.data.frame(pct_result) || is.matrix(pct_result)) {
    pct_vec <- pct_result[, 1]
  } else {
    pct_vec <- as.numeric(pct_result)
  }
  # 确保长度对
  if (length(pct_vec) != ncol(obj)) {
    pct_vec <- rep(0, ncol(obj))
  }
  list(pct = pct_vec, n_genes = length(mito_genes))
}

# ============================================================
# 导出单个对象的 metadata (按客户格式)
# ============================================================
export_obj <- function(obj, sample_name, species_label, raw_sample) {
  obj <- tryCatch(JoinLayers(obj), error = function(e) obj)
  md <- obj@meta.data

  mito_info <- calc_mito(obj)
  mito_pct <- round(mito_info$pct, 3)

  if (species_label == "mouse") {
    pct_human <- rep(0, ncol(obj))
    pct_mouse <- rep(100, ncol(obj))
    species_call <- rep("mouse", ncol(obj))
  } else if (species_label == "human") {
    pct_human <- rep(100, ncol(obj))
    pct_mouse <- rep(0, ncol(obj))
    species_call <- rep("human", ncol(obj))
  } else {
    pct_human <- rep(NA, ncol(obj))
    pct_mouse <- rep(NA, ncol(obj))
    species_call <- rep("undefine", ncol(obj))
  }

  if ("celltype_refined" %in% colnames(md)) {
    cell_anno <- as.character(md$celltype_refined)
  } else if ("celltype_v32" %in% colnames(md)) {
    cell_anno <- as.character(md$celltype_v32)
  } else if ("final_celltype" %in% colnames(md)) {
    cell_anno <- as.character(md$final_celltype)
  } else {
    cell_anno <- as.character(Idents(obj))
  }

  orig_ident <- if ("orig.ident" %in% colnames(md)) {
    as.character(md$orig.ident)
  } else {
    rep(sample_name, ncol(obj))
  }

  data.frame(
    barcode = colnames(obj),
    orig.ident = orig_ident,
    nCount_RNA = md$nCount_RNA,
    nFeature_RNA = md$nFeature_RNA,
    mito = mito_pct,
    Sample = sample_name,
    percent.human = pct_human,
    percent.mouse = pct_mouse,
    species_call = species_call,
    raw_Sample = raw_sample,
    percent.mito = mito_pct,
    mitorelatedgenes = mito_info$n_genes,
    CellAnnotation = cell_anno,
    stringsAsFactors = FALSE
  )
}

# ============================================================
# Step 1: 导出 metadata CSV
# ============================================================
cat("========================================\n")
cat("STEP 1: Export cell metadata (client format)\n")
cat("========================================\n\n")

samples_mouse <- c("liver3", "liver8", "gan56", "bla_liver", "con_liver",
                   "brain3", "brain8", "nao56", "bla_brain", "con_brain")
samples_human <- c("liver3", "liver8", "gan56")

all_dfs <- list()
objects_cache <- list()  # 缓存对象供 Step 2 使用

cat("--- Mouse cells ---\n")
for (s in samples_mouse) {
  fp <- tryCatch(resolve_rds_path(s, "v4"), error = function(e) NULL)
  if (is.null(fp) || !file.exists(fp)) {
    cat(sprintf("  ⚠ %s: RDS not found\n", s))
    next
  }
  obj <- readRDS(fp)
  obj <- tryCatch(JoinLayers(obj), error = function(e) obj)
  raw_s <- RAW_SAMPLE_MAP[[s]] %||% s

  df <- export_obj(obj, s, "mouse", raw_s)
  out_csv <- file.path(OUT_DIR, sprintf("%s_mouse_celltype.csv", s))
  write.csv(df, out_csv, row.names = FALSE)
  cat(sprintf("  ✓ %s mouse: %d cells, %d cell types\n",
              s, nrow(df), length(unique(df$CellAnnotation))))
  all_dfs[[paste0(s, "_mouse")]] <- df
  objects_cache[[paste0(s, "_mouse")]] <- obj
  # 不 rm, 下面要用
}

cat("\n--- Human cells (PDX) ---\n")
for (s in samples_human) {
  candidates <- c(
    file.path(PROJECT_ROOT, s, "human", sprintf("%s_human_annotated_v4.rds", s)),
    file.path(PROJECT_ROOT, s, "human", sprintf("%s_human_manual_annotated_v4.rds", s))
  )
  human_rds <- candidates[file.exists(candidates)]
  if (length(human_rds) == 0) {
    cat(sprintf("  ⚠ %s human: RDS not found\n", s))
    next
  }
  obj <- readRDS(human_rds[1])
  obj <- tryCatch(JoinLayers(obj), error = function(e) obj)
  raw_s <- RAW_SAMPLE_MAP[[s]] %||% s

  df <- export_obj(obj, s, "human", raw_s)
  out_csv <- file.path(OUT_DIR, sprintf("%s_human_celltype.csv", s))
  write.csv(df, out_csv, row.names = FALSE)
  cat(sprintf("  ✓ %s human: %d cells, %d cell types\n",
              s, nrow(df), length(unique(df$CellAnnotation))))
  all_dfs[[paste0(s, "_human")]] <- df
  objects_cache[[paste0(s, "_human")]] <- obj
}

cat("\n--- Merging ---\n")
merged <- do.call(rbind, all_dfs)
rownames(merged) <- NULL
write.csv(merged, file.path(OUT_DIR, "all_samples_celltype.csv"), row.names = FALSE)
cat(sprintf("✓ merged: %d cells total\n", nrow(merged)))
cat("\ncells by Sample × species:\n")
print(table(merged$Sample, merged$species_call))

# ============================================================
# Step 2: Top 30 markers per cell type
# ============================================================
cat("\n\n========================================\n")
cat("STEP 2: Compute top 30 markers per cell type\n")
cat("========================================\n\n")

# 策略: 把同物种 + 同器官的样本合并, 然后 FindAllMarkers
# 4 个合并对象:
#   (a) liver_mouse (5 个样本 host liver)
#   (b) brain_mouse (5 个样本 host brain)
#   (c) liver_human (3 个 PDX HCC)

compute_markers <- function(obj_list, group_name, ident_col = "celltype_refined") {
  cat(sprintf("\n--- %s ---\n", group_name))
  if (length(obj_list) == 0) {
    cat("  no objects, skip\n")
    return(NULL)
  }

  cat(sprintf("  merging %d objects...\n", length(obj_list)))

  # merge
  if (length(obj_list) == 1) {
    merged_obj <- obj_list[[1]]
  } else {
    merged_obj <- merge(obj_list[[1]],
                        y = obj_list[-1],
                        add.cell.ids = names(obj_list))
    merged_obj <- tryCatch(JoinLayers(merged_obj), error = function(e) merged_obj)
  }

  cat(sprintf("  merged: %d cells × %d genes\n",
              ncol(merged_obj), nrow(merged_obj)))

  # 设置 Idents 为 celltype_refined
  md <- merged_obj@meta.data
  if (ident_col %in% colnames(md)) {
    Idents(merged_obj) <- md[[ident_col]]
  } else if ("celltype_v32" %in% colnames(md)) {
    Idents(merged_obj) <- md$celltype_v32
  } else {
    cat("  ⚠ no celltype column, using seurat_clusters\n")
  }

  # 过滤掉太小的 cell type (< 10 cells 算不出稳定 marker)
  ident_counts <- table(Idents(merged_obj))
  valid_idents <- names(ident_counts[ident_counts >= 10])
  cat(sprintf("  cell types: %d total, %d with >=10 cells\n",
              length(ident_counts), length(valid_idents)))

  if (length(valid_idents) == 0) {
    cat("  ⚠ no valid cell types\n")
    return(NULL)
  }

  merged_obj <- subset(merged_obj, idents = valid_idents)

  # Normalize if needed
  if (!"data" %in% Layers(merged_obj[["RNA"]])) {
    merged_obj <- NormalizeData(merged_obj, verbose = FALSE)
  }

  # FindAllMarkers
  cat("  Running FindAllMarkers (this may take 2-5 min)...\n")
  t0 <- proc.time()
  markers <- tryCatch({
    FindAllMarkers(merged_obj,
                   only.pos = TRUE,
                   min.pct = 0.25,
                   logfc.threshold = 0.25,
                   verbose = FALSE)
  }, error = function(e) {
    cat(sprintf("  ✗ FindAllMarkers failed: %s\n", e$message))
    NULL
  })

  if (is.null(markers) || nrow(markers) == 0) {
    cat("  ⚠ no markers found\n")
    return(NULL)
  }

  elapsed <- (proc.time() - t0)["elapsed"] / 60
  cat(sprintf("  ✓ %d marker rows in %.1f min\n", nrow(markers), elapsed))

  # Top 30 per cluster
  top30 <- markers %>%
    group_by(cluster) %>%
    arrange(desc(avg_log2FC)) %>%
    slice_head(n = 30) %>%
    ungroup() %>%
    select(cluster, gene, avg_log2FC, pct.1, pct.2, p_val, p_val_adj)

  # 排序: cluster 名 + log2FC
  top30 <- top30 %>%
    arrange(cluster, desc(avg_log2FC))

  top30_csv <- file.path(MARKER_DIR, sprintf("%s_top30_markers.csv", group_name))
  write.csv(top30, top30_csv, row.names = FALSE)
  cat(sprintf("  ✓ wrote %d rows to %s\n", nrow(top30), basename(top30_csv)))

  # 也写一份 wide format (每个 cell type 一列, 30 行)
  top30_wide <- top30 %>%
    group_by(cluster) %>%
    mutate(rank = row_number()) %>%
    ungroup() %>%
    select(cluster, gene, rank) %>%
    tidyr::pivot_wider(names_from = cluster, values_from = gene)
  wide_csv <- file.path(MARKER_DIR, sprintf("%s_top30_markers_wide.csv", group_name))
  write.csv(top30_wide, wide_csv, row.names = FALSE)

  # 也写完整 markers 表 (不只 top30) 作为附录
  full_csv <- file.path(MARKER_DIR, sprintf("%s_all_markers.csv", group_name))
  write.csv(markers %>% select(cluster, gene, avg_log2FC, pct.1, pct.2,
                               p_val, p_val_adj) %>%
              arrange(cluster, desc(avg_log2FC)),
            full_csv, row.names = FALSE)
  cat(sprintf("  ✓ wrote full markers (%d rows) to %s\n",
              nrow(markers), basename(full_csv)))

  rm(merged_obj); gc(verbose = FALSE)
  return(top30)
}

# (a) liver mouse (host 微环境)
liver_mouse_objs <- objects_cache[c("liver3_mouse", "liver8_mouse",
                                     "gan56_mouse", "bla_liver_mouse",
                                     "con_liver_mouse")]
liver_mouse_objs <- liver_mouse_objs[!sapply(liver_mouse_objs, is.null)]
compute_markers(liver_mouse_objs, "liver_mouse_host")

# (b) brain mouse
brain_mouse_objs <- objects_cache[c("brain3_mouse", "brain8_mouse",
                                     "nao56_mouse", "bla_brain_mouse",
                                     "con_brain_mouse")]
brain_mouse_objs <- brain_mouse_objs[!sapply(brain_mouse_objs, is.null)]
compute_markers(brain_mouse_objs, "brain_mouse_host")

# (c) liver human (PDX HCC)
liver_human_objs <- objects_cache[c("liver3_human", "liver8_human",
                                     "gan56_human")]
liver_human_objs <- liver_human_objs[!sapply(liver_human_objs, is.null)]
compute_markers(liver_human_objs, "liver_human_PDX")

cat(sprintf("\n\n========================================\n"))
cat(sprintf("ALL DONE\n"))
cat(sprintf("Output: %s\n", OUT_DIR))
cat(sprintf("  - 10 mouse + 3 human CSV (cell metadata)\n"))
cat(sprintf("  - all_samples_celltype.csv (merged)\n"))
cat(sprintf("  - top30_markers/ (top30 + full markers × 3 groups)\n"))
cat(sprintf("========================================\n"))
