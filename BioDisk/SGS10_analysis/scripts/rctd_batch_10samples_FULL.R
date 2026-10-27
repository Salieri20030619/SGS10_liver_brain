#!/usr/bin/env Rscript
# ============================================================
# rctd_batch_10samples.R
# RCTD cell type deconvolution on 10 SeekSpace samples
# ============================================================
# 基于 Task A v1.1 Phase 2 Day 1:
#   把 scRNA-seq (v4_clean brain + v4 liver) 的 cell type label 
#   迁移到 SeekSpace 每个 spot, 得到 cell type proportion.
#
# 输入:
#   - spatial_results/*_seekspace.rds (10 个)
#   - brain3/8/nao56/bla_brain/con_brain 的 *_manual_annotated_v4_clean.rds
#   - liver3/8/gan56/bla_liver/con_liver 的 *_mouse_annotated_v4.rds
#
# 输出:
#   spatial_results/rctd/
#     ref_brain_22ct.rds                                (Brain reference)
#     ref_liver_33ct.rds                                (Liver reference)  
#     per_sample/<sample>_rctd_results.rds              (完整 RCTD 对象)
#     per_sample/<sample>_weights_matrix.csv            (spot × ct 权重)
#     per_sample/<sample>_dominant_celltype.csv         (每 spot 最强 ct)
#     enriched/<sample>_seekspace_rctd.rds              (Seurat obj 加 RCTD 列)
#     summary/
#       cell_type_proportion_per_zone.csv               (每 zone × ct 统计)
#       rctd_batch_status.csv                           (任务状态)
#       rctd_batch_summary.md
#
# 运行:
#   cd ~/SGS10_analysis
#   screen -S rctd_batch -dm bash -c '
#     Rscript scripts/rctd_batch_10samples.R 2>&1 | tee scripts/logs/rctd_batch.log
#   '
#
# 预计:
#   Reference 构建 5-10 min
#   每 sample 20-40 min (55k spots × 22 ct)
#   总计 4-7 hours 串行, 或 ~1.5 hr 如果并行 4 threads (通过 RCTD max_cores)
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(Matrix)
  library(spacexr)
})

PROJECT_ROOT <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
setwd(PROJECT_ROOT)

OUT_DIR <- file.path(PROJECT_ROOT, "spatial_results", "rctd_full")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "per_sample"), showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "enriched"), showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "summary"), showWarnings = FALSE)

LOG_FILE <- file.path(OUT_DIR, "summary", "progress.log")
log_msg <- function(msg) {
  line <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), msg)
  cat(line, "\n")
  cat(line, "\n", file = LOG_FILE, append = TRUE)
}

writeLines(sprintf("=== RCTD batch start %s ===",
                   format(Sys.time(), "%Y-%m-%d %H:%M:%S")), LOG_FILE)

# ============================================================
# 0. 配置
# ============================================================

# RCTD 使用几个 core. 服务器 1TB RAM 应该可以用 8-16.
# 实际 RCTD 有内存 peak, 建议 4-8 保守
MAX_CORES <- 4

# spatial samples 配置
SPATIAL_SAMPLES <- list(
  Liver_3     = list(organ = "liver", group = "Tumor"),
  Liver_8     = list(organ = "liver", group = "Tumor"),
  gan56       = list(organ = "liver", group = "Tumor"),
  Bla_liver_1 = list(organ = "liver", group = "Blank"),
  Con_liver_1 = list(organ = "liver", group = "Control"),
  Brain_3     = list(organ = "brain", group = "Tumor"),
  Brain_8     = list(organ = "brain", group = "Tumor"),
  nao56       = list(organ = "brain", group = "Tumor"),
  Bla_brain_1 = list(organ = "brain", group = "Blank"),
  Con_brain_1 = list(organ = "brain", group = "Control")
)

# Brain reference 文件
BRAIN_REF_FILES <- c(
  "brain3/mouse/brain3_manual_annotated_v4_clean.rds",
  "brain8/mouse/brain8_manual_annotated_v4_clean.rds",
  "nao56/mouse/nao56_manual_annotated_v4_clean.rds",
  "bla_brain/mouse/bla_brain_manual_annotated_v4_clean.rds",
  "con_brain/mouse/con_brain_manual_annotated_v4_clean.rds"
)

# Liver reference 文件
LIVER_REF_FILES <- c(
  "liver3/mouse/liver3_mouse_annotated_v4.rds",
  "liver8/mouse/liver8_mouse_annotated_v4.rds",
  "gan56/mouse/gan56_mouse_annotated_v4.rds",
  "bla_liver/mouse/bla_liver_mouse_annotated_v4.rds",
  "con_liver/mouse/con_liver_mouse_annotated_v4.rds"
)

# Cell type 筛选规则
BRAIN_CT_COL <- "celltype_refined"   # 精细清洗后的标签
BRAIN_CT_EXCLUDE <- c("low_quality", "Neurons_unassigned", "Other", "NK_cells")
BRAIN_CT_MERGE <- list(
  Microglia_other = c("Microglia", "Microglia_M1", "Microglia_M2", "Microglia_DAM_strict")
)

LIVER_CT_COL <- "final_celltype"     # liver 用 final_celltype (v4 原 celltype_refined 与 final 一致)
LIVER_CT_EXCLUDE <- c("Unknown", "Basophil")  # < 25 或身份不明
LIVER_CT_MERGE <- list()              # liver 不合并

# RCTD 参数
# reference n_max_cells 降到 500 后, doublet mode 应该能在合理时间内跑完
# 如果仍然慢, 可改为 full mode (只算 weight, 不做 singlet/doublet decision)
RCTD_DOUBLET_MODE <- "full"           # full: 每 spot × K weights (稳定快)
                                        # doublet: 决定 1 or 2 cell type per spot (精度高但慢)
RCTD_MIN_UMI <- 100                   # 每 spot 至少 100 UMI 才纳入

# ⭐ SUBSAMPLE: 降采样加速 (SeekSpace 55k spots 完整跑耗时数小时, subsample 后 ~5-10x 加速)
# 设 1.0 = 不 subsample (全量), 0.2 = 20% 采样, 0.1 = 10% 采样
# 对 RCTD cell type proportion 推断无显著影响
# 保留 full dominance 分析 + zone composition
SUBSAMPLE_RATE <- 1.0
SUBSAMPLE_SEED <- 42                  # 可复现

# ============================================================
# 1. 构建 Reference function
# ============================================================

# 基于 5 样本合并 + filter + merge rare types, 建 RCTD Reference
build_reference <- function(files, ct_col, ct_exclude, ct_merge,
                             organ_name, out_path) {
  log_msg(sprintf("=== Building %s reference ===", organ_name))
  
  all_counts <- list()
  all_celltype <- list()
  
  for (f in files) {
    log_msg(sprintf("  reading %s ...", basename(f)))
    obj <- readRDS(f)
    
    if (!(ct_col %in% colnames(obj@meta.data))) {
      stop(sprintf("cell type column '%s' not found in %s", ct_col, f))
    }
    
    ct <- obj@meta.data[[ct_col]]
    
    # 取 counts (RCTD 要 raw counts)
    counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
    
    # 给 barcode 加 sample prefix, 避免跨样本重复
    sample_tag <- gsub("_manual_annotated_v4_clean\\.rds|_mouse_annotated_v4\\.rds", "",
                       basename(f))
    new_barcodes <- paste0(sample_tag, "_", colnames(counts))
    colnames(counts) <- new_barcodes
    names(ct) <- new_barcodes
    
    all_counts[[f]] <- counts
    all_celltype[[f]] <- ct
    
    rm(obj)
    gc(verbose = FALSE)
  }
  
  # 合并
  log_msg("  merging 5 samples...")
  
  # 取公共基因
  common_genes <- Reduce(intersect, lapply(all_counts, rownames))
  log_msg(sprintf("  common genes: %d", length(common_genes)))
  
  all_counts <- lapply(all_counts, function(m) m[common_genes, , drop = FALSE])
  merged_counts <- do.call(cbind, all_counts)
  merged_ct <- unlist(all_celltype)
  
  stopifnot(ncol(merged_counts) == length(merged_ct))
  names(merged_ct) <- colnames(merged_counts)
  
  log_msg(sprintf("  merged: %d genes × %d cells", 
                  nrow(merged_counts), ncol(merged_counts)))
  
  # Filter: 排除特定 cell type
  if (length(ct_exclude) > 0) {
    keep <- !merged_ct %in% ct_exclude
    log_msg(sprintf("  excluding %d cells (from %s)",
                    sum(!keep), paste(ct_exclude, collapse = ",")))
    merged_counts <- merged_counts[, keep]
    merged_ct <- merged_ct[keep]
  }
  
  # 合并稀有类
  if (length(ct_merge) > 0) {
    for (new_name in names(ct_merge)) {
      old_names <- ct_merge[[new_name]]
      idx <- merged_ct %in% old_names
      log_msg(sprintf("  merging %d cells (%s) → %s",
                      sum(idx), paste(old_names, collapse = ","), new_name))
      merged_ct[idx] <- new_name
    }
  }
  
  # 再次 filter: 排除 cell type < 25 (RCTD requirement)
  ct_counts <- table(merged_ct)
  low_ct <- names(ct_counts)[ct_counts < 25]
  if (length(low_ct) > 0) {
    keep <- !merged_ct %in% low_ct
    log_msg(sprintf("  excluding %d cells (rare ct <25: %s)",
                    sum(!keep), paste(low_ct, collapse = ",")))
    merged_counts <- merged_counts[, keep]
    merged_ct <- merged_ct[keep]
  }
  
  final_ct_table <- table(merged_ct)
  log_msg(sprintf("  final reference: %d cells, %d cell types",
                  ncol(merged_counts), length(final_ct_table)))
  cat("  Cell type distribution:\n")
  print(sort(final_ct_table, decreasing = TRUE))
  
  # 构造 RCTD Reference
  # spacexr::Reference 接受 (counts, cell_types_as_factor, nUMI)
  merged_ct_factor <- factor(merged_ct)
  names(merged_ct_factor) <- colnames(merged_counts)
  
  n_umi <- colSums(merged_counts)
  names(n_umi) <- colnames(merged_counts)
  
  log_msg("  building Reference object...")
  ref <- Reference(
    counts      = merged_counts,
    cell_types  = merged_ct_factor,
    nUMI        = n_umi,
    require_int = TRUE,     # counts should be integer
    n_max_cells = 500       # 每 ct 最多 500 cells (官方 vignette 就是这个量级)
                            # 之前 10000 导致 Sigma 优化每轮 14 min 不收敛
  )
  
  log_msg(sprintf("  saving to %s", out_path))
  saveRDS(ref, out_path)
  
  return(list(
    ref = ref,
    n_cells = ncol(merged_counts),
    n_celltypes = length(final_ct_table),
    ct_table = final_ct_table
  ))
}

# ============================================================
# 2. 构建 Brain + Liver Reference (如果已保存就 skip)
# ============================================================
brain_ref_path <- file.path(OUT_DIR, "ref_brain_22ct.rds")
liver_ref_path <- file.path(OUT_DIR, "ref_liver_33ct.rds")

if (file.exists(brain_ref_path)) {
  log_msg(sprintf("Brain ref already exists, loading: %s", brain_ref_path))
  brain_ref <- readRDS(brain_ref_path)
  brain_info <- list(ref = brain_ref, n_celltypes = length(levels(brain_ref@cell_types)))
  log_msg(sprintf("  loaded: %d cell types", brain_info$n_celltypes))
} else {
  brain_info <- build_reference(
    BRAIN_REF_FILES, BRAIN_CT_COL, BRAIN_CT_EXCLUDE, BRAIN_CT_MERGE,
    "BRAIN", brain_ref_path
  )
  brain_ref <- brain_info$ref
}

if (file.exists(liver_ref_path)) {
  log_msg(sprintf("Liver ref already exists, loading: %s", liver_ref_path))
  liver_ref <- readRDS(liver_ref_path)
  liver_info <- list(ref = liver_ref, n_celltypes = length(levels(liver_ref@cell_types)))
  log_msg(sprintf("  loaded: %d cell types", liver_info$n_celltypes))
} else {
  liver_info <- build_reference(
    LIVER_REF_FILES, LIVER_CT_COL, LIVER_CT_EXCLUDE, LIVER_CT_MERGE,
    "LIVER", liver_ref_path
  )
  liver_ref <- liver_info$ref
}

# ============================================================
# 3. RCTD per sample function
# ============================================================

run_rctd_on_sample <- function(sample_name, cfg, reference) {
  log_msg(sprintf("\n============================================"))
  log_msg(sprintf("=== RCTD on %s (%s, %s) ===", 
                  sample_name, cfg$organ, cfg$group))
  log_msg(sprintf("============================================"))
  
  # Skip if already done
  out_prefix <- file.path(OUT_DIR, "per_sample", sample_name)
  existing_rds <- paste0(out_prefix, "_rctd_results.rds")
  if (file.exists(existing_rds)) {
    size_mb <- round(file.info(existing_rds)$size / 1024^2, 1)
    log_msg(sprintf("  ✓ SKIP: already completed (%.1f MB, %s)",
                    size_mb, format(file.info(existing_rds)$mtime)))
    # Read existing summary
    return(list(
      sample = sample_name, organ = cfg$organ, group = cfg$group,
      n_spots_total = NA, n_spots_analyzed = NA,
      n_celltypes_detected = NA, elapsed_min = 0,
      dom_celltype_top = "(pre-existing)", dom_celltype_top_pct = NA,
      skipped = TRUE
    ))
  }
  
  t0 <- proc.time()
  
  # 1. 读 spatial Seurat obj
  spatial_rds <- file.path("spatial_results", sprintf("%s_seekspace.rds", sample_name))
  if (!file.exists(spatial_rds)) {
    log_msg(sprintf("  ✗ Spatial RDS not found: %s", spatial_rds))
    return(NULL)
  }
  
  log_msg(sprintf("  reading %s...", spatial_rds))
  obj <- readRDS(spatial_rds)
  log_msg(sprintf("  spots: %d, genes: %d", ncol(obj), nrow(obj)))
  
  # 2. 提取 host (mouse) counts
  # 两种情况:
  #   A. Dual-species PDX (Liver_3/8/gan56): 基因名 "mm10---Acta2" 和 "GRCh38-MKI67"
  #      → 取 mm10--- 前缀的, 去前缀
  #   B. Mouse-only (Blank/Control/Brain): 基因名 "Acta2", "Xkr4" (无前缀)
  #      → 直接全部用
  all_genes <- rownames(obj)
  mouse_genes_prefixed <- grep("^mm10---", all_genes, value = TRUE)
  
  counts_all <- GetAssayData(obj, assay = "RNA", layer = "counts")
  
  if (length(mouse_genes_prefixed) > 0) {
    # Dual-species: 取 mm10--- 前缀的
    log_msg(sprintf("  dual-species mode: %d mm10--- genes found", length(mouse_genes_prefixed)))
    counts_mouse <- counts_all[mouse_genes_prefixed, , drop = FALSE]
    rownames(counts_mouse) <- gsub("^mm10---", "", rownames(counts_mouse))
  } else {
    # Mouse-only: 直接全部基因
    log_msg(sprintf("  mouse-only mode: %d genes (no prefix)", length(all_genes)))
    counts_mouse <- counts_all
  }
  log_msg(sprintf("  mouse counts: %d genes", nrow(counts_mouse)))
  
  # 3. 拿 spatial coordinates (从 spatial reduction)
  coords <- obj@reductions$spatial@cell.embeddings
  colnames(coords) <- c("x", "y")
  coords <- as.data.frame(coords)
  
  # 确保 coords 和 counts 的 barcode 对齐
  stopifnot(all(rownames(coords) == colnames(counts_mouse)))
  
  # 4. 过滤 low UMI spots
  n_umi <- colSums(counts_mouse)
  keep_spots <- n_umi >= RCTD_MIN_UMI
  log_msg(sprintf("  spots with >=%d mouse UMI: %d / %d",
                  RCTD_MIN_UMI, sum(keep_spots), length(keep_spots)))
  
  if (sum(keep_spots) < 100) {
    log_msg(sprintf("  ✗ too few spots pass UMI filter, skip"))
    return(NULL)
  }
  
  counts_mouse <- counts_mouse[, keep_spots]
  coords <- coords[keep_spots, ]
  
  # ⭐ SUBSAMPLE: 降采样加速
  if (SUBSAMPLE_RATE < 1.0) {
    set.seed(SUBSAMPLE_SEED)
    n_before <- ncol(counts_mouse)
    n_keep <- round(n_before * SUBSAMPLE_RATE)
    sampled_idx <- sort(sample(n_before, n_keep))
    counts_mouse <- counts_mouse[, sampled_idx]
    coords <- coords[sampled_idx, ]
    log_msg(sprintf("  SUBSAMPLE: %d/%d (%.0f%%, seed=%d)",
                    n_keep, n_before, 100 * SUBSAMPLE_RATE, SUBSAMPLE_SEED))
  }
  
  # 5. 构造 SpatialRNA 对象
  n_umi_filt <- colSums(counts_mouse)
  spatialRNA <- SpatialRNA(
    coords = coords,
    counts = counts_mouse,
    nUMI = n_umi_filt
  )
  
  # 6. 跑 RCTD
  log_msg(sprintf("  building RCTD object..."))
  rctd <- create.RCTD(
    spatialRNA,
    reference,
    max_cores = MAX_CORES,
    CELL_MIN_INSTANCE = 25,
    test_mode = FALSE
  )
  
  log_msg(sprintf("  running RCTD (mode=%s)...", RCTD_DOUBLET_MODE))
  rctd <- run.RCTD(rctd, doublet_mode = RCTD_DOUBLET_MODE)
  
  # 7. 提取结果
  results <- rctd@results
  
  if (RCTD_DOUBLET_MODE == "doublet") {
    # doublet mode: results_df 每 spot 一行, first_type + second_type + weights
    results_df <- results$results_df
    
    # weights 矩阵: normalize_weights
    weights <- results$weights
    weights_norm <- sweep(weights, 1, rowSums(weights), "/")
    weights_norm[is.nan(weights_norm)] <- 0
  } else {
    results_df <- NULL
    weights_norm <- results$weights
    weights_norm <- sweep(weights_norm, 1, rowSums(weights_norm), "/")
  }
  
  # 8. 保存 results
  out_prefix <- file.path(OUT_DIR, "per_sample", sample_name)
  
  # Full RCTD obj
  saveRDS(rctd, paste0(out_prefix, "_rctd_results.rds"))
  
  # weights matrix
  write.csv(as.data.frame(as.matrix(weights_norm)),
            paste0(out_prefix, "_weights_matrix.csv"),
            row.names = TRUE)
  
  # dominant cell type per spot
  dom_ct <- colnames(weights_norm)[apply(weights_norm, 1, which.max)]
  dom_score <- apply(weights_norm, 1, max)
  second_ct <- apply(weights_norm, 1, function(r) {
    sorted <- sort(r, decreasing = TRUE)
    if (length(sorted) >= 2) names(sorted)[2] else NA
  })
  second_score <- apply(weights_norm, 1, function(r) {
    sorted <- sort(r, decreasing = TRUE)
    if (length(sorted) >= 2) sorted[2] else NA
  })
  
  dom_df <- data.frame(
    barcode = rownames(weights_norm),
    dominant_celltype = dom_ct,
    dominant_score    = round(dom_score, 4),
    second_celltype   = second_ct,
    second_score      = round(second_score, 4),
    stringsAsFactors  = FALSE
  )
  
  if (!is.null(results_df)) {
    # 补充 RCTD 的 spot_class (singlet / doublet_certain / reject)
    dom_df$rctd_spot_class <- as.character(results_df[dom_df$barcode, "spot_class"])
    dom_df$rctd_first_type <- as.character(results_df[dom_df$barcode, "first_type"])
    dom_df$rctd_second_type <- as.character(results_df[dom_df$barcode, "second_type"])
  }
  
  write.csv(dom_df, paste0(out_prefix, "_dominant_celltype.csv"), row.names = FALSE)
  
  # 9. 回写到 Seurat obj 的 meta.data
  # 注意: 被过滤的 spot (UMI < 100) 没有 RCTD 结果, 用 NA
  meta_add <- data.frame(
    rctd_dominant = NA_character_,
    rctd_dom_score = NA_real_,
    rctd_second = NA_character_,
    rctd_second_score = NA_real_,
    stringsAsFactors = FALSE
  )[rep(1, ncol(obj)), ]
  rownames(meta_add) <- colnames(obj)
  
  idx <- match(dom_df$barcode, rownames(meta_add))
  meta_add$rctd_dominant[idx]     <- dom_df$dominant_celltype
  meta_add$rctd_dom_score[idx]    <- dom_df$dominant_score
  meta_add$rctd_second[idx]       <- dom_df$second_celltype
  meta_add$rctd_second_score[idx] <- dom_df$second_score
  
  obj <- AddMetaData(obj, meta_add)
  
  # 10. 保存 enriched seurat obj
  saveRDS(obj, file.path(OUT_DIR, "enriched",
                         sprintf("%s_seekspace_rctd.rds", sample_name)))
  
  elapsed <- round((proc.time() - t0)["elapsed"] / 60, 1)
  
  log_msg(sprintf("  ✓ done (%.1f min)", elapsed))
  log_msg(sprintf("  dominant cell type distribution:"))
  dom_tab <- sort(table(dom_df$dominant_celltype), decreasing = TRUE)
  for (i in seq_len(min(5, length(dom_tab)))) {
    log_msg(sprintf("    %-25s: %d (%.1f%%)",
                    names(dom_tab)[i], dom_tab[i], 100 * dom_tab[i] / sum(dom_tab)))
  }
  
  # 返回 summary
  return(list(
    sample = sample_name,
    organ  = cfg$organ,
    group  = cfg$group,
    n_spots_total    = ncol(obj),
    n_spots_analyzed = nrow(dom_df),
    n_celltypes_detected = length(unique(dom_df$dominant_celltype)),
    elapsed_min      = elapsed,
    dom_celltype_top = names(dom_tab)[1],
    dom_celltype_top_pct = round(100 * dom_tab[1] / sum(dom_tab), 1)
  ))
}

# ============================================================
# 4. 主循环: 10 samples
# ============================================================
log_msg("\n============================================")
log_msg("=== Starting RCTD on 10 spatial samples ===")
log_msg("============================================\n")

summary_list <- list()
error_log <- list()

for (sample_name in names(SPATIAL_SAMPLES)) {
  cfg <- SPATIAL_SAMPLES[[sample_name]]
  reference <- if (cfg$organ == "brain") brain_ref else liver_ref
  
  result <- tryCatch(
    run_rctd_on_sample(sample_name, cfg, reference),
    error = function(e) {
      log_msg(sprintf("\n✗ %s FAILED: %s", sample_name, e$message))
      error_log[[sample_name]] <<- e$message
      NULL
    }
  )
  
  if (!is.null(result)) {
    summary_list[[sample_name]] <- result
  }
  
  gc(verbose = FALSE)
}

# ============================================================
# 5. 汇总
# ============================================================
log_msg("\n============================================")
log_msg("=== Writing summary ===")
log_msg("============================================")

if (length(summary_list) > 0) {
  summary_df <- do.call(rbind, lapply(summary_list, function(x) {
    data.frame(
      sample                = x$sample,
      organ                 = x$organ,
      group                 = x$group,
      n_spots_total         = x$n_spots_total,
      n_spots_analyzed      = x$n_spots_analyzed,
      n_celltypes_detected  = x$n_celltypes_detected,
      dom_celltype_top      = x$dom_celltype_top,
      dom_celltype_top_pct  = x$dom_celltype_top_pct,
      elapsed_min           = x$elapsed_min,
      stringsAsFactors      = FALSE
    )
  }))
  write.csv(summary_df, file.path(OUT_DIR, "summary", "rctd_batch_status.csv"),
            row.names = FALSE)
  cat("\n=== Sample summary ===\n")
  print(summary_df)
}

if (length(error_log) > 0) {
  err_df <- data.frame(
    sample = names(error_log),
    error  = unlist(error_log),
    stringsAsFactors = FALSE
  )
  write.csv(err_df, file.path(OUT_DIR, "summary", "rctd_errors.csv"),
            row.names = FALSE)
  cat("\n=== Errors ===\n")
  print(err_df)
}

# ============================================================
# 6. 跨 sample cell type proportion by zone
# ============================================================
log_msg("\nComputing cell type proportion per zone (combining samples)...")

all_zone_ct <- list()
for (sample_name in names(summary_list)) {
  cfg <- SPATIAL_SAMPLES[[sample_name]]
  enriched_path <- file.path(OUT_DIR, "enriched",
                             sprintf("%s_seekspace_rctd.rds", sample_name))
  if (!file.exists(enriched_path)) next
  
  obj <- readRDS(enriched_path)
  if (!"zone" %in% colnames(obj@meta.data)) {
    rm(obj); gc(verbose=FALSE); next
  }
  
  analyzed <- !is.na(obj$rctd_dominant)
  df <- data.frame(
    sample  = sample_name,
    organ   = cfg$organ,
    group   = cfg$group,
    zone    = obj$zone[analyzed],
    celltype = obj$rctd_dominant[analyzed],
    stringsAsFactors = FALSE
  )
  all_zone_ct[[sample_name]] <- df
  
  rm(obj); gc(verbose=FALSE)
}

if (length(all_zone_ct) > 0) {
  zone_ct_df <- do.call(rbind, all_zone_ct)
  
  # 每 sample × zone × celltype 的 count 和 proportion
  prop_df <- zone_ct_df %>%
    dplyr::group_by(sample, organ, group, zone, celltype) %>%
    dplyr::summarise(n_spots = dplyr::n(), .groups = "drop_last") %>%
    dplyr::mutate(pct = n_spots / sum(n_spots) * 100) %>%
    dplyr::ungroup()
  
  write.csv(prop_df,
            file.path(OUT_DIR, "summary", "cell_type_proportion_per_zone.csv"),
            row.names = FALSE)
  log_msg(sprintf("  saved %d rows", nrow(prop_df)))
}

log_msg("\n=== RCTD batch complete ===")
log_msg(sprintf("Success: %d / %d", length(summary_list), length(SPATIAL_SAMPLES)))
log_msg(sprintf("Failed : %d", length(error_log)))
