#!/usr/bin/env Rscript
# ============================================================
# cellchat_baseline_control_v1.R   --  T6 supplement
# ============================================================
#
# 目的: 补跑 baseline + control 的 spatial CellChat
#       用与 4/20 cellchat_niche_communication.R 完全一致的方法
#       让 3 组 (Tumor / Blank / Control) 形成对称 + 可比较的数据集
#
# 输入 (与 4/20 完全一致):
#   - spatial_results/<sample>_seekspace.rds (gene expression)
#   - spatial_results/rctd/per_sample/<sample>_dominant_celltype.csv
#
# 设计:
#   liver_blank:   Bla_liver_1 (1 sample, 14074 spots, 25 cell types)
#   liver_control: Con_liver_1 (1 sample, 10201 spots, 25 cell types)
#   brain_blank:   Bla_brain_1 (1 sample, 12248 spots, 22 cell types)
#   brain_control: Con_brain_1 (1 sample, 12731 spots, 22 cell types)
#
# 物种处理 (dual-species liver):
#   mm10--- prefix gene rename to bare name
#   drop GRCh38- (PDX human reads)
#
# 输出:
#   taskB/cellchat/
#     liver_blank_cellchat.rds
#     liver_blank_lr_pairs.csv
#     liver_blank_pathway_pairs.csv
#     liver_blank_aggregate_network.pdf
#     liver_blank_count_matrix.csv          ⭐ 数量表 (客户要求)
#     liver_blank_weight_matrix.csv         ⭐ 权重表 (客户要求)
#     liver_blank_top_pairs.csv             top cell-cell pair 排序
#     ... (4 group × 7 文件 = 28 文件)
#
# 与 4/20 一致的关键参数:
#   - CellChat 创建: createCellChat (no spatial coords)
#   - DB: CellChatDB.mouse, all categories
#   - computeCommunProb: type = "triMean", raw.use = TRUE
#   - filterCommunication: min.cells = 10
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(dplyr)
  library(patchwork)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)
future::plan("multisession", workers = 4)
options(future.globals.maxSize = 64 * 1024^3)

PROJECT_ROOT  <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
SPATIAL_DIR   <- file.path(PROJECT_ROOT, "spatial_results")
RCTD_DIR      <- file.path(SPATIAL_DIR, "rctd", "per_sample")
OUT_DIR       <- file.path(PROJECT_ROOT, "taskB", "cellchat")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

setwd(PROJECT_ROOT)

cat("============================================================\n")
cat("T6 supplement: CellChat for Baseline + Control (Liver + Brain)\n")
cat(sprintf("Start: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("============================================================\n")

# ============================================================
# Groups (4 个新组, 与 4/20 tumor 互补)
# ============================================================
GROUPS <- list(
  liver_blank   = c("Bla_liver_1"),
  liver_control = c("Con_liver_1"),
  brain_blank   = c("Bla_brain_1"),
  brain_control = c("Con_brain_1")
)

MIN_CELLS_PER_CT <- 10  # 跟 4/20 完全一致

# ============================================================
# Function: load + prep one sample (复用 4/20 逻辑)
# ============================================================
load_for_cellchat <- function(sample) {
  rds_path  <- file.path(SPATIAL_DIR, paste0(sample, "_seekspace.rds"))
  rctd_path <- file.path(RCTD_DIR, paste0(sample, "_dominant_celltype.csv"))
  
  if (!file.exists(rds_path)) {
    cat(sprintf("  ✗ MISSING rds: %s\n", rds_path))
    return(NULL)
  }
  if (!file.exists(rctd_path)) {
    cat(sprintf("  ✗ MISSING rctd: %s\n", rctd_path))
    return(NULL)
  }
  
  cat(sprintf("  Loading %s ... ", sample))
  t0 <- Sys.time()
  obj <- readRDS(rds_path)
  
  # Rename rownames: mm10---Gene -> Gene, drop GRCh38- (human)
  old_names <- rownames(obj)
  has_mm10  <- grepl("^mm10---", old_names)
  has_grch  <- grepl("^GRCh38-", old_names)
  
  if (any(has_mm10) || any(has_grch)) {
    mouse_features <- old_names[has_mm10 | (!has_mm10 & !has_grch)]
    new_names <- sub("^mm10---", "", mouse_features)
    
    counts <- GetAssayData(obj, layer = "counts")[mouse_features, ]
    rownames(counts) <- new_names
    
    if (any(duplicated(new_names))) {
      cat(sprintf("\n    ⚠ Duplicated bare names after rename, keeping first\n"))
      keep_idx <- !duplicated(new_names)
      counts <- counts[keep_idx, ]
    }
    
    obj2 <- CreateSeuratObject(counts = counts, meta.data = obj@meta.data,
                               project = sample)
    obj2 <- NormalizeData(obj2, verbose = FALSE)
    rm(obj)
    obj <- obj2
    cat(sprintf("rebuilt mouse-only (%d features) ... ", nrow(obj)))
  } else {
    # brain (mouse only) doesn't need rename, but still NormalizeData
    obj <- NormalizeData(obj, verbose = FALSE)
  }
  
  # Merge RCTD dominant_celltype
  rctd_df <- read.csv(rctd_path, stringsAsFactors = FALSE)
  
  # 看 csv 列名 (一般是 spot_id + dominant_celltype)
  if (!"dominant_celltype" %in% colnames(rctd_df)) {
    cat(sprintf("\n    ⚠ no dominant_celltype col, found: %s\n",
                paste(colnames(rctd_df), collapse=",")))
    return(NULL)
  }
  
  # 第一列是 barcode/spot_id
  rownames(rctd_df) <- rctd_df[[1]]
  
  md <- obj@meta.data
  md$dominant_celltype <- rctd_df[rownames(md), "dominant_celltype"]
  obj@meta.data <- md
  
  # Filter: keep only spots with RCTD label
  keep <- !is.na(obj$dominant_celltype)
  obj <- obj[, keep]
  
  obj$sample_id <- sample
  
  cat(sprintf("%d spots kept (%.1fs)\n",
              ncol(obj), as.numeric(Sys.time() - t0)))
  obj
}

# ============================================================
# Run CellChat for one group
# ============================================================
run_cellchat_group <- function(group_name, samples) {
  cat(sprintf("\n\n#### === %s (samples: %s) ===\n",
              group_name, paste(samples, collapse = ", ")))
  
  t_total <- Sys.time()
  
  # ---- Load + merge samples ----
  obj_list <- lapply(samples, load_for_cellchat)
  obj_list <- obj_list[!sapply(obj_list, is.null)]
  
  if (length(obj_list) == 0) {
    cat(sprintf("  ✗ No valid samples for %s, skip\n", group_name))
    return(NULL)
  }
  
  if (length(obj_list) == 1) {
    obj <- obj_list[[1]]
  } else {
    # multi-sample merge (虽然 baseline 都是单样本, 留 future-proof)
    obj <- merge(obj_list[[1]], obj_list[-1], add.cell.ids = samples)
    obj <- JoinLayers(obj)
    obj <- NormalizeData(obj, verbose = FALSE)
  }
  
  cat(sprintf("  Total spots: %d\n", ncol(obj)))
  cat(sprintf("  Cell types: %d\n", length(unique(obj$dominant_celltype))))
  
  ct_counts <- table(obj$dominant_celltype)
  cat("  Cell type spots:\n")
  for (ct in names(sort(ct_counts, decreasing = TRUE))) {
    cat(sprintf("    %-30s %d\n", ct, ct_counts[ct]))
  }
  
  # ---- CellChat workflow (与 4/20 完全一致) ----
  cat("  createCellChat -> preprocessing -> computeCommunProb...\n")
  
  data_input <- GetAssayData(obj, layer = "data")
  meta <- data.frame(
    labels  = obj$dominant_celltype,
    samples = obj$sample_id,
    row.names = colnames(obj)
  )
  
  cellchat <- createCellChat(
    object   = data_input,
    meta     = meta,
    group.by = "labels"
  )
  
  CellChatDB <- CellChatDB.mouse
  CellChatDB.use <- CellChatDB
  cellchat@DB <- CellChatDB.use
  
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  n_over_inter <- nrow(cellchat@LR$LRsig)
  cat(sprintf("  over-expressed LR pairs: %d\n", n_over_inter))
  
  # 关键: 与 4/20 完全一致的 CellChat 参数
  t_comm <- proc.time()
  cellchat <- computeCommunProb(
    cellchat,
    type            = "triMean",
    raw.use         = TRUE,
    population.size = FALSE
  )
  cat(sprintf("  computeCommunProb done in %.1f min\n",
              (proc.time() - t_comm)["elapsed"] / 60))
  
  cellchat <- filterCommunication(cellchat, min.cells = MIN_CELLS_PER_CT)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  # ---- 提取数据 ----
  df_lr      <- subsetCommunication(cellchat)
  df_pathway <- subsetCommunication(cellchat, slot.name = "netP")
  
  cat(sprintf("  ✓ significant LR pairs: %d, pathways: %d\n",
              nrow(df_lr), length(unique(df_pathway$pathway_name))))
  
  # ---- 保存 ----
  prefix <- file.path(OUT_DIR, group_name)
  
  # 1. cellchat object
  saveRDS(cellchat, paste0(prefix, "_cellchat.rds"))
  
  # 2. LR pairs
  write.csv(df_lr,      paste0(prefix, "_lr_pairs.csv"),      row.names = FALSE)
  write.csv(df_pathway, paste0(prefix, "_pathway_pairs.csv"), row.names = FALSE)
  
  # 3. ⭐ 数量表 + 权重表 (客户重点要求)
  count_mat  <- cellchat@net$count
  weight_mat <- cellchat@net$weight
  write.csv(count_mat,  paste0(prefix, "_count_matrix.csv"),  row.names = TRUE)
  write.csv(weight_mat, paste0(prefix, "_weight_matrix.csv"), row.names = TRUE)
  
  # 4. top cell-cell pairs (long format, sorted by weight)
  pair_df <- expand.grid(
    sender   = rownames(count_mat),
    receiver = colnames(count_mat),
    stringsAsFactors = FALSE
  )
  pair_df$count  <- count_mat[as.matrix(pair_df[, c("sender", "receiver")])]
  pair_df$weight <- weight_mat[as.matrix(pair_df[, c("sender", "receiver")])]
  pair_df <- pair_df[pair_df$count > 0, ]
  pair_df <- pair_df[order(-pair_df$weight), ]
  write.csv(pair_df, paste0(prefix, "_top_pairs.csv"), row.names = FALSE)
  
  # 5. aggregate network plot (CellChat 默认圆图, 客户也要)
  tryCatch({
    pdf(paste0(prefix, "_aggregate_network.pdf"), width = 14, height = 7)
    par(mfrow = c(1, 2), xpd = TRUE)
    netVisual_circle(
      count_mat, vertex.weight = as.numeric(table(meta$labels)),
      weight.scale = TRUE, label.edge = FALSE,
      title.name = "Number of interactions"
    )
    netVisual_circle(
      weight_mat, vertex.weight = as.numeric(table(meta$labels)),
      weight.scale = TRUE, label.edge = FALSE,
      title.name = "Interaction weight/strength"
    )
    dev.off()
    cat(sprintf("  ✓ aggregate_network.pdf saved\n"))
  }, error = function(e) {
    cat(sprintf("  ⚠ plot fail: %s\n", e$message))
  })
  
  cat(sprintf("\n  ✓✓ DONE %s in %.1f min\n",
              group_name,
              as.numeric(Sys.time() - t_total, units = "mins")))
  
  # 返回 summary
  list(
    group         = group_name,
    n_spots       = ncol(obj),
    n_cell_types  = length(unique(obj$dominant_celltype)),
    n_lr_pairs    = nrow(df_lr),
    n_pathways    = length(unique(df_pathway$pathway_name)),
    total_count   = sum(count_mat),
    total_weight  = sum(weight_mat)
  )
}

# ============================================================
# 主循环
# ============================================================
all_summary <- list()
for (gname in names(GROUPS)) {
  s <- tryCatch(
    run_cellchat_group(gname, GROUPS[[gname]]),
    error = function(e) {
      cat(sprintf("\n  ✗✗ FAIL %s: %s\n", gname, e$message))
      NULL
    }
  )
  if (!is.null(s)) all_summary[[gname]] <- s
  
  # 释放内存
  gc(verbose = FALSE)
}

# ============================================================
# 写 summary
# ============================================================
cat("\n\n============================================================\n")
cat("FINAL SUMMARY\n")
cat("============================================================\n")

if (length(all_summary) > 0) {
  summary_df <- do.call(rbind, lapply(all_summary, function(s) {
    data.frame(
      group        = s$group,
      n_spots      = s$n_spots,
      n_cell_types = s$n_cell_types,
      n_lr_pairs   = s$n_lr_pairs,
      n_pathways   = s$n_pathways,
      total_count  = s$total_count,
      total_weight = round(s$total_weight, 4),
      stringsAsFactors = FALSE
    )
  }))
  
  write.csv(summary_df, file.path(OUT_DIR, "summary_baseline_control.csv"),
            row.names = FALSE)
  
  print(summary_df)
}

cat(sprintf("\nDone: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(sprintf("Output: %s\n", OUT_DIR))
cat("\nNew files in taskB/cellchat/:\n")
system(sprintf("ls -la %s/{liver_blank,liver_control,brain_blank,brain_control}_*", OUT_DIR))
