#!/usr/bin/env Rscript
# ============================================================
# cellchat_liver_expanded.R --- liver CellChat 扩展版本
# ============================================================
# 用法:
#   cd ~/SGS10_analysis
#   Rscript scripts/cellchat_liver_expanded.R 2>&1 | tee scripts/logs/cellchat_liver_expanded.log
#
# 改进 vs cellchat_v4_analysis.R:
#   1. 扩展 liver 虚拟细类映射 (5 → 22), 加入 T/NK/DC/B/各 macrophage 亚群/LSEC 等
#   2. 修复 HSC bug: HSC_Activated_ECM 从 HSC 范畴移除 (避免 for 循环覆盖)
#   3. 仍然用 NAM clean CSV (NAM 仅 72 cells, 主要在 Tumor)
#   4. 只跑 liver 3 个 task (brain 不变, 因为 brain 没有这个映射问题)
#
# 输出:
#   cellchat_results_expanded/
#     liver_{Tumor,Blank,Control}/
#       (同 cellchat_v4_analysis.R 的输出结构)
#     comparison/
#       liver_*_diff.csv
#     summary/
#       task_status.csv
#       cellchat_summary.md
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(patchwork)
  library(dplyr)
  library(ggplot2)
})

select <- dplyr::select; rename <- dplyr::rename; filter <- dplyr::filter
mutate <- dplyr::mutate; arrange <- dplyr::arrange; bind_rows <- dplyr::bind_rows

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 100 * 1024^3)

# ============================================================
# 0. 路径
# ============================================================
PROJECT_ROOT <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
SCRIPT_DIR   <- file.path(PROJECT_ROOT, "scripts")
setwd(PROJECT_ROOT)
source(file.path(SCRIPT_DIR, "diff_common.R"))

OUT_DIR <- file.path(PROJECT_ROOT, "cellchat_results_expanded")
for (d in c("comparison", "summary")) {
  dir.create(file.path(OUT_DIR, d), recursive = TRUE, showWarnings = FALSE)
}

ERROR_LOG <- file.path(OUT_DIR, "summary", "error_log.txt")
writeLines(sprintf("== CellChat liver expanded - %s ==",
                   format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ERROR_LOG)
write_error <- function(context, msg) {
  line <- sprintf("[%s] %s => %s",
                  format(Sys.time(), "%H:%M:%S"), context, msg)
  cat("  ⚠", line, "\n")
  cat(line, "\n", file = ERROR_LOG, append = TRUE)
}

# ============================================================
# 1. 扩展的 liver 细类映射
# ============================================================
# 原 LIVER_CELL_TYPES 只有 8 个虚拟类 (其中 All_Cells/HSC 因 bug 不参与 CellChat),
# 实际只 5 个细类进 CellChat. 现在扩展到 22 个.
#
# CRITICAL: HSC_Activated_ECM 从 HSC 范畴移除, 避免 for 循环覆盖 bug

LIVER_CELL_TYPES_EXPANDED <- list(
  # 原有 (修复 HSC)
  Hepatocyte_PP        = "Hep_Periportal",
  Hepatocyte_PV        = "Hep_Perivenous",
  HSC                  = c("HSC_Quiescent", "HSC_Activated_Chemo",
                            "HSC_Activated_Remodel", "HSC_Proliferating"),
  HSC_Activated_ECM    = "HSC_Activated_ECM",
  Kupffer              = "Kupffer_KC",
  Macrophage_NAM       = "Macrophage_NAM",  # 用 NAM clean CSV 覆盖

  # 新增 myeloid (之前漏掉)
  Macrophage_M1        = "Macrophage_M1",
  Macrophage_M2        = "Macrophage_M2",
  Macrophage_MDM       = "Macrophage_MDM",
  LCM                  = "LCM",  # liver capsular macrophage

  # 新增 DC
  DC_cDC1              = "DC_cDC1",
  DC_cDC2              = "DC_cDC2",
  DC_pDC               = "DC_pDC",

  # 新增 T cells
  T_CD4                = "T_CD4",
  T_CD8                = "T_CD8",
  T_NKT                = "T_NKT",
  T_Treg               = "T_Treg",

  # 新增 NK
  NK                   = "NK",

  # 新增 B cells (合并三个亚型)
  B_cell               = c("B_Memory", "B_Naive", "B_Plasma"),

  # 新增 Neutrophil
  Neutrophil           = "Neutrophil",

  # 新增 endothelial 亚群
  LSEC                 = c("LSEC", "LSEC_Pericentral", "LSEC_Periportal"),
  LVEC                 = "LVEC",

  # 新增 Cholangiocyte
  Cholangiocyte        = "Cholangiocyte"
)

cat("========================================\n")
cat("CellChat liver expanded analysis\n")
cat("Start:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================\n\n")

cat(sprintf("Expanded LIVER_CELL_TYPES: %d virtual cell types\n",
            length(LIVER_CELL_TYPES_EXPANDED)))
for (vct in names(LIVER_CELL_TYPES_EXPANDED)) {
  raw <- LIVER_CELL_TYPES_EXPANDED[[vct]]
  cat(sprintf("  %-22s ← %s\n", vct, paste(raw, collapse = " + ")))
}
cat("\n")

# ============================================================
# 2. 参数
# ============================================================
COMM_PROB_TYPE     <- "triMean"
RAW_USE            <- TRUE
POPULATION_SIZE    <- TRUE
MIN_CELLS_PER_TYPE <- 10

MIN_FREE_GB_TO_START <- 200
MIN_FREE_GB_ABORT    <- 100
MAX_R_PROCS_ALLOWED  <- 2

LIVER_SAMPLES_BY_GROUP <- list(
  Tumor   = c("liver3", "liver8", "gan56"),
  Blank   = c("bla_liver"),
  Control = c("con_liver")
)
GROUPS <- c("Tumor", "Blank", "Control")

# ============================================================
# 3. 安全检查
# ============================================================
get_free_mem_gb <- function() {
  tryCatch({
    out <- system("free -g | awk 'NR==2 {print $7}'", intern = TRUE)
    as.numeric(out)
  }, error = function(e) NA)
}
get_r_proc_count <- function() {
  tryCatch({
    out <- system("pgrep -u $USER -c '^(R|Rscript)$' 2>/dev/null", intern = TRUE)
    n <- suppressWarnings(as.integer(out)); if (is.na(n)) 0 else n
  }, error = function(e) 0)
}
preflight_check <- function() {
  cat("========== Pre-run safety check ==========\n")
  n_r <- get_r_proc_count()
  cat(sprintf("  Residual R processes: %d (limit: %d)\n",
              n_r, MAX_R_PROCS_ALLOWED))
  if (n_r > MAX_R_PROCS_ALLOWED) {
    stop(sprintf("ABORT: %d R 进程残留", n_r))
  }
  free_gb <- get_free_mem_gb()
  cat(sprintf("  Available memory: %.0f GB\n",
              if (is.na(free_gb)) 0 else free_gb))
  if (is.na(free_gb) || free_gb < MIN_FREE_GB_TO_START) {
    stop("ABORT: 内存不足")
  }
  cat("  ✓ Safety check passed\n\n")
}
preflight_check()

# ============================================================
# 4. 加载 NAM clean CSV
# ============================================================
NAM_CLEAN_CSV <- file.path(PROJECT_ROOT, "nam_clean_results",
                            "v4_refined_clean_full.csv")
NAM_CLEAN_ANNO <- NULL
if (file.exists(NAM_CLEAN_CSV)) {
  cat("========== Loading NAM clean annotation ==========\n")
  NAM_CLEAN_ANNO <- read.csv(NAM_CLEAN_CSV, stringsAsFactors = FALSE)
  cat(sprintf("  Loaded %d cells\n\n", nrow(NAM_CLEAN_ANNO)))
} else {
  cat("⚠ NAM clean CSV not found\n\n")
}

# ============================================================
# 5. CellChat database
# ============================================================
cat("========== Loading CellChatDB.mouse ==========\n")
CellChatDB <- CellChatDB.mouse
cat(sprintf("  total LR pairs: %d\n\n", nrow(CellChatDB$interaction)))
CellChatDB.use <- CellChatDB

# ============================================================
# 6. 加载 liver samples
# ============================================================
cat("========== Loading liver samples ==========\n")
samples <- unlist(LIVER_SAMPLES_BY_GROUP)
objs <- list()
for (s in samples) {
  fp <- resolve_rds_path(s, "v4")
  if (!file.exists(fp)) { write_error(s, paste0("missing: ", fp)); next }
  o <- tryCatch({
    x <- readRDS(fp)
    if (inherits(x, "Seurat")) x <- tryCatch(JoinLayers(x), error = function(e) x)
    x <- RenameCells(x, add.cell.id = s)
    x[["sample"]] <- s
    grp <- NA
    for (g in names(LIVER_SAMPLES_BY_GROUP)) {
      if (s %in% LIVER_SAMPLES_BY_GROUP[[g]]) { grp <- g; break }
    }
    x[["group"]] <- grp
    x
  }, error = function(e) { write_error(s, e$message); NULL })
  if (!is.null(o)) {
    objs[[s]] <- o
    cat(sprintf("  ✓ %-12s n=%d  group=%s\n", s, ncol(o), o$group[1]))
  }
}

cat("  merging...\n")
seurat_liver <- merge(objs[[1]], objs[-1])
seurat_liver <- tryCatch(JoinLayers(seurat_liver), error = function(e) seurat_liver)
rm(objs); gc(verbose = FALSE)
cat(sprintf("  merged: %d cells\n\n", ncol(seurat_liver)))

# ============================================================
# 7. build_celltype_labels (扩展版 + NAM clean join + HSC fix)
# ============================================================
build_celltype_labels_expanded <- function(seurat_obj) {
  ct_col <- NULL
  for (cand in c("celltype_refined", "v4_refined", "final_celltype")) {
    if (cand %in% colnames(seurat_obj@meta.data)) { ct_col <- cand; break }
  }
  if (is.null(ct_col)) stop("No cell type column")

  raw_labels <- as.character(seurat_obj@meta.data[[ct_col]])
  cell_ids   <- colnames(seurat_obj)

  # NAM clean join
  if (!is.null(NAM_CLEAN_ANNO)) {
    n_v4_nam <- sum(raw_labels == "Macrophage_NAM", na.rm = TRUE)
    match_idx <- match(cell_ids, NAM_CLEAN_ANNO$cell_id)
    has_clean <- !is.na(match_idx)
    is_v4_nam <- !is.na(raw_labels) & raw_labels == "Macrophage_NAM"
    apply_mask <- is_v4_nam & has_clean
    if (sum(apply_mask) > 0) {
      new_label <- NAM_CLEAN_ANNO$v4_refined_clean[match_idx[apply_mask]]
      new_label[new_label == "Macrophage_NAM_clean"] <- "Macrophage_NAM"
      raw_labels[apply_mask] <- new_label
    }
    n_after <- sum(raw_labels == "Macrophage_NAM", na.rm = TRUE)
    n_misannot <- sum(raw_labels == "Macrophage_NAM_misannotated", na.rm = TRUE)
    cat(sprintf("    [NAM clean] before %d v4 NAM, after %d clean kept, %d misannotated\n",
                n_v4_nam, n_after, n_misannot))
  }

  # 反向映射: raw label → virtual label (用扩展版 ct_map)
  # CRITICAL: 把 HSC_Activated_ECM 放在 HSC 之后, 即使覆盖也能保证 ECM 正确
  # (实际上扩展版的 HSC 已经不包含 ECM, 不会冲突)
  virtual_labels <- rep("Other", length(raw_labels))
  for (vct in names(LIVER_CELL_TYPES_EXPANDED)) {
    matching_raw <- LIVER_CELL_TYPES_EXPANDED[[vct]]
    virtual_labels[raw_labels %in% matching_raw] <- vct
  }
  return(virtual_labels)
}

# ============================================================
# 8. 单组 CellChat 分析
# ============================================================
run_one_cellchat <- function(seurat_obj, group_name, task_idx) {
  task_key <- paste0("liver_", group_name)
  task_dir <- file.path(OUT_DIR, task_key)
  plot_dir <- file.path(task_dir, "plots")
  dir.create(task_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  result <- list(
    key = task_key, group = group_name, status = "error",
    n_cells = NA, n_celltypes = NA, n_lr_pairs_sig = NA,
    n_pathways = NA, elapsed_min = NA
  )
  t0 <- proc.time()

  cat(sprintf("\n### Task %d: %s ###\n", task_idx, task_key))
  free_gb <- get_free_mem_gb()
  cat(sprintf("  free mem: %.0f GB\n",
              if (is.na(free_gb)) 0 else free_gb))
  if (!is.na(free_gb) && free_gb < MIN_FREE_GB_ABORT) {
    write_error(task_key, sprintf("free %.0f < %d", free_gb, MIN_FREE_GB_ABORT))
    return(result)
  }

  tryCatch({
    cell_idx <- which(seurat_obj$group == group_name)
    if (length(cell_idx) < 100) stop("group too small")
    sub <- seurat_obj[, cell_idx]
    sub <- tryCatch(JoinLayers(sub), error = function(e) sub)
    cat(sprintf("  n_cells=%d\n", ncol(sub)))

    virtual_labels <- build_celltype_labels_expanded(sub)
    sub$cellchat_label <- virtual_labels

    label_counts <- table(virtual_labels)
    cat("  cell type counts (sorted):\n")
    sorted_counts <- sort(label_counts, decreasing = TRUE)
    for (lab in names(sorted_counts)) {
      mark <- if (sorted_counts[[lab]] >= MIN_CELLS_PER_TYPE) "✓" else "✗"
      cat(sprintf("    %s %-25s : %d\n", mark, lab, sorted_counts[[lab]]))
    }

    valid_types <- names(label_counts)[label_counts >= MIN_CELLS_PER_TYPE]
    valid_types <- setdiff(valid_types, "Other")
    if (length(valid_types) < 3) stop(sprintf("valid types < 3 (%d)", length(valid_types)))

    keep_idx <- which(sub$cellchat_label %in% valid_types)
    sub <- sub[, keep_idx]
    cat(sprintf("  filtered: %d cells, %d types\n", ncol(sub), length(valid_types)))
    result$n_cells <- ncol(sub)
    result$n_celltypes <- length(valid_types)

    # CellChat
    if (length(sub@assays$RNA@layers) > 0 &&
        "data" %in% names(sub@assays$RNA@layers) &&
        sum(sub@assays$RNA@layers$data) > 0) {
      data_input <- GetAssayData(sub, assay = "RNA", layer = "data")
    } else {
      sub <- NormalizeData(sub, verbose = FALSE)
      data_input <- GetAssayData(sub, assay = "RNA", layer = "data")
    }

    meta <- data.frame(
      labels = sub$cellchat_label,
      row.names = colnames(sub),
      stringsAsFactors = FALSE
    )

    cat("  createCellChat → preprocessing...\n")
    cellchat <- createCellChat(object = data_input, meta = meta, group.by = "labels")
    cellchat@DB <- CellChatDB.use
    cellchat <- subsetData(cellchat)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)

    n_over_inter <- nrow(cellchat@LR$LRsig)
    cat(sprintf("  over-expressed LR pairs: %d\n", n_over_inter))

    cat("  computeCommunProb...\n")
    t_comm <- proc.time()
    cellchat <- computeCommunProb(cellchat,
                                   type = COMM_PROB_TYPE,
                                   raw.use = RAW_USE,
                                   population.size = POPULATION_SIZE)
    cat(sprintf("  done in %.1f min\n",
                (proc.time() - t_comm)["elapsed"] / 60))

    cellchat <- filterCommunication(cellchat, min.cells = MIN_CELLS_PER_TYPE)
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)

    df_lr <- subsetCommunication(cellchat)
    df_pathway <- subsetCommunication(cellchat, slot.name = "netP")

    result$n_lr_pairs_sig <- nrow(df_lr)
    result$n_pathways <- length(unique(df_pathway$pathway_name))
    cat(sprintf("  ✓ LR pairs: %d, pathways: %d\n",
                result$n_lr_pairs_sig, result$n_pathways))

    # 保存
    saveRDS(cellchat, file.path(task_dir, "cellchat_obj.rds"))
    write.csv(df_lr, file.path(task_dir, "lr_pairs_filtered.csv"), row.names = FALSE)
    write.csv(df_pathway, file.path(task_dir, "pathway_summary.csv"), row.names = FALSE)
    write.csv(cellchat@net$count, file.path(task_dir, "interactions_count.csv"))
    write.csv(cellchat@net$weight, file.path(task_dir, "interactions_weight.csv"))

    # 关键的可视化
    tryCatch({
      pdf(file.path(plot_dir, "circle_count.pdf"), width = 11, height = 11)
      par(xpd = TRUE)
      netVisual_circle(cellchat@net$count,
                       vertex.weight = as.numeric(table(cellchat@idents)),
                       weight.scale = TRUE, label.edge = FALSE,
                       title.name = sprintf("%s · interaction count", task_key))
      dev.off()
    }, error = function(e) write_error(paste0(task_key, "_circle_count"), e$message))

    tryCatch({
      pdf(file.path(plot_dir, "circle_weight.pdf"), width = 11, height = 11)
      par(xpd = TRUE)
      netVisual_circle(cellchat@net$weight,
                       vertex.weight = as.numeric(table(cellchat@idents)),
                       weight.scale = TRUE, label.edge = FALSE,
                       title.name = sprintf("%s · interaction weight", task_key))
      dev.off()
    }, error = function(e) write_error(paste0(task_key, "_circle_weight"), e$message))

    if (length(cellchat@netP$pathways) > 0) {
      tryCatch({
        ht <- netAnalysis_signalingRole_heatmap(cellchat,
                                                 pattern = "outgoing",
                                                 height = max(6, length(cellchat@netP$pathways) * 0.2),
                                                 width = 10, font.size = 7)
        pdf(file.path(plot_dir, "heatmap_pathway_outgoing.pdf"),
            width = 11, height = max(6, length(cellchat@netP$pathways) * 0.2))
        print(ht)
        dev.off()
      }, error = function(e) write_error(paste0(task_key, "_ht_out"), e$message))

      tryCatch({
        ht <- netAnalysis_signalingRole_heatmap(cellchat,
                                                 pattern = "incoming",
                                                 height = max(6, length(cellchat@netP$pathways) * 0.2),
                                                 width = 10, font.size = 7)
        pdf(file.path(plot_dir, "heatmap_pathway_incoming.pdf"),
            width = 11, height = max(6, length(cellchat@netP$pathways) * 0.2))
        print(ht)
        dev.off()
      }, error = function(e) write_error(paste0(task_key, "_ht_in"), e$message))
    }

    result$status <- "ok"
    result$elapsed_min <- round((proc.time() - t0)["elapsed"] / 60, 1)
    cat(sprintf("  ✓ DONE %s: cells=%d, types=%d, LR=%d, pathways=%d (%.1f min)\n",
                task_key, result$n_cells, result$n_celltypes,
                result$n_lr_pairs_sig, result$n_pathways, result$elapsed_min))

    rm(cellchat, sub, data_input, meta)
    for (j in 1:3) gc(verbose = FALSE)

  }, error = function(e) {
    result$elapsed_min <<- round((proc.time() - t0)["elapsed"] / 60, 1)
    write_error(task_key, e$message)
    cat("  ✗ ERROR:", e$message, "\n")
  })

  return(result)
}

# ============================================================
# 9. 跑 3 个 liver tasks
# ============================================================
t_all_start <- proc.time()
all_results <- list()

for (i in seq_along(GROUPS)) {
  grp <- GROUPS[i]
  res <- run_one_cellchat(seurat_liver, grp, i)
  all_results[[paste0("liver_", grp)]] <- res
}

t_all_elapsed <- (proc.time() - t_all_start)["elapsed"]
cat(sprintf("\n\n3 liver tasks elapsed: %.1f min\n", t_all_elapsed / 60))

# 释放
rm(seurat_liver)
for (j in 1:5) gc(verbose = FALSE)

# ============================================================
# 10. task_status.csv
# ============================================================
ts_df <- bind_rows(lapply(all_results, function(r) {
  data.frame(
    key = r$key, group = r$group, status = r$status,
    n_cells = if (is.null(r$n_cells)) NA else r$n_cells,
    n_celltypes = if (is.null(r$n_celltypes)) NA else r$n_celltypes,
    n_lr_pairs_sig = if (is.null(r$n_lr_pairs_sig)) NA else r$n_lr_pairs_sig,
    n_pathways = if (is.null(r$n_pathways)) NA else r$n_pathways,
    elapsed_min = if (is.null(r$elapsed_min)) NA else r$elapsed_min,
    stringsAsFactors = FALSE
  )
}))
write.csv(ts_df, file.path(OUT_DIR, "summary", "task_status.csv"),
          row.names = FALSE)
cat("\n========== task_status.csv ==========\n")
print(ts_df)

# ============================================================
# 11. 跨组对比 (简化版, 只做 LR diff 不做 mergeCellChat)
# ============================================================
cat("\n========== Liver cross-group LR diff ==========\n")
cc_list <- list()
for (grp in GROUPS) {
  fp <- file.path(OUT_DIR, paste0("liver_", grp), "cellchat_obj.rds")
  if (file.exists(fp)) {
    cc_list[[grp]] <- tryCatch(readRDS(fp), error = function(e) NULL)
  }
}
cc_list <- cc_list[!sapply(cc_list, is.null)]

if (length(cc_list) >= 2) {
  comparisons <- list(
    c("Tumor", "Control"),
    c("Tumor", "Blank"),
    c("Blank", "Control")
  )

  for (cmp in comparisons) {
    if (!all(cmp %in% names(cc_list))) next
    cmp_name <- paste0(cmp[1], "_vs_", cmp[2])
    cat(sprintf("  %s\n", cmp_name))

    df1 <- subsetCommunication(cc_list[[cmp[1]]])
    df2 <- subsetCommunication(cc_list[[cmp[2]]])
    df1$lr_key <- paste(df1$source, df1$target, df1$interaction_name, sep = "__")
    df2$lr_key <- paste(df2$source, df2$target, df2$interaction_name, sep = "__")

    only_in_1 <- setdiff(df1$lr_key, df2$lr_key)
    only_in_2 <- setdiff(df2$lr_key, df1$lr_key)
    shared    <- intersect(df1$lr_key, df2$lr_key)

    diff_df <- data.frame(
      lr_key = c(only_in_1, only_in_2, shared),
      group1 = cmp[1], group2 = cmp[2],
      category = c(
        rep(paste0("only_in_", cmp[1]), length(only_in_1)),
        rep(paste0("only_in_", cmp[2]), length(only_in_2)),
        rep("shared", length(shared))
      ),
      stringsAsFactors = FALSE
    )
    df1_lookup <- df1[match(diff_df$lr_key, df1$lr_key),
                      c("source", "target", "interaction_name", "pathway_name", "prob")]
    df2_lookup <- df2[match(diff_df$lr_key, df2$lr_key), c("prob")]
    diff_df$source <- df1_lookup$source
    diff_df$target <- df1_lookup$target
    diff_df$interaction_name <- df1_lookup$interaction_name
    diff_df$pathway_name <- df1_lookup$pathway_name
    diff_df$prob_group1 <- df1_lookup$prob
    diff_df$prob_group2 <- df2_lookup
    diff_df$prob_diff <- diff_df$prob_group1 - diff_df$prob_group2

    only2_idx <- which(diff_df$category == paste0("only_in_", cmp[2]))
    if (length(only2_idx) > 0) {
      df2_for_only2 <- df2[match(diff_df$lr_key[only2_idx], df2$lr_key),
                           c("source", "target", "interaction_name", "pathway_name", "prob")]
      diff_df$source[only2_idx] <- df2_for_only2$source
      diff_df$target[only2_idx] <- df2_for_only2$target
      diff_df$interaction_name[only2_idx] <- df2_for_only2$interaction_name
      diff_df$pathway_name[only2_idx] <- df2_for_only2$pathway_name
      diff_df$prob_group2[only2_idx] <- df2_for_only2$prob
    }

    diff_df <- diff_df[order(abs(diff_df$prob_diff), decreasing = TRUE,
                              na.last = TRUE), ]
    write.csv(diff_df,
              file.path(OUT_DIR, "comparison",
                        paste0("liver_", cmp_name, "_diff.csv")),
              row.names = FALSE)
    cat(sprintf("    %s: %d total (only_%s=%d, only_%s=%d, shared=%d)\n",
                cmp_name, nrow(diff_df),
                cmp[1], length(only_in_1),
                cmp[2], length(only_in_2),
                length(shared)))
  }
}

# ============================================================
# 12. summary md
# ============================================================
summary_lines <- c(
  "# CellChat liver expanded 分析总结",
  "",
  sprintf("**生成时间:** %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  sprintf("**总耗时:** %.1f min", t_all_elapsed / 60),
  "",
  "## 改进 vs cellchat_v4_analysis.R",
  "",
  "1. **扩展虚拟细类**: 5 → 22 (加入 T/NK/DC/B/Macrophage 各亚群/LSEC 等)",
  "2. **修复 HSC bug**: HSC_Activated_ECM 不再被 for 循环覆盖, HSC 包含 4 个剩余亚群",
  "3. **NAM clean**: 使用 nam_clean_v2.R 输出的清洗后 NAM (72 cells)",
  "",
  "## 任务状态",
  "",
  "| 任务 | 细胞数 | 细类数 | LR pairs | Pathways | 时间 (min) |",
  "|---|---:|---:|---:|---:|---:|"
)
for (i in seq_len(nrow(ts_df))) {
  r <- ts_df[i, ]
  summary_lines <- c(summary_lines,
    sprintf("| %s | %s | %s | %s | %s | %s |",
            r$key,
            if (is.na(r$n_cells)) "-" else as.character(r$n_cells),
            if (is.na(r$n_celltypes)) "-" else as.character(r$n_celltypes),
            if (is.na(r$n_lr_pairs_sig)) "-" else as.character(r$n_lr_pairs_sig),
            if (is.na(r$n_pathways)) "-" else as.character(r$n_pathways),
            if (is.na(r$elapsed_min)) "-" else sprintf("%.1f", r$elapsed_min)))
}
summary_lines <- c(summary_lines, "",
  "## 关键问题",
  "",
  "1. **TGFB1-TGFBR pairs**: 扩展细类后是否能检测到? 哪个细类是 source?",
  "2. **IL6-IL6R pairs**: 同上",
  "3. **HSC ↔ T cell / NK / DC 通讯**: 之前看不到",
  "4. **NAM (cleaned 72) 真实通讯网络**: 之前 NAM 在 5 cell types 里通讯模式")

writeLines(summary_lines, file.path(OUT_DIR, "summary", "cellchat_summary.md"))

cat("\n========================================\n")
cat("CellChat liver expanded done\n")
cat(sprintf("  Output: %s\n", OUT_DIR))
cat(sprintf("  Total: %.1f min\n", t_all_elapsed / 60))
cat("========================================\n")
