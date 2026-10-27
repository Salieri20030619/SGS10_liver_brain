#!/usr/bin/env Rscript
# ============================================================
# cellchat_v4_analysis.R --- CellChat 细胞通讯分析
# ============================================================
# 用法:
#   cd ~/SGS10_analysis
#   Rscript scripts/cellchat_v4_analysis.R 2>&1 | tee scripts/logs/cellchat.log
#
# 设计:
#   - 6 个子任务: liver × {Tumor, Blank, Control} + brain × {Tumor, Blank, Control}
#   - 每个子任务独立: subset → CellChat → computeCommunProb → 保存
#   - 6 个子任务完成后: 肝 3 组 mergeCellChat + 对比, 脑 3 组同样
#   - 复用 diff_common.R 的 sample/cell type 定义
#   - 内存防护: 三层 (preflight + runtime + 任务间释放)
#
# 输出:
#   cellchat_results/
#     {liver,brain}_{Tumor,Blank,Control}/
#       cellchat_obj.rds        ← CellChat 对象
#       lr_pairs_filtered.csv   ← 显著的 LR 配对
#       pathway_summary.csv     ← 每条通路的汇总
#       interactions_count.csv  ← cell-cell 交互数矩阵
#       interactions_weight.csv ← cell-cell 交互强度矩阵
#       plots/
#         circle_count.pdf
#         circle_weight.pdf
#         heatmap_pathway.pdf
#         bubble_lr.pdf
#     comparison/
#       liver_tumor_vs_control_diff.csv
#       liver_tumor_vs_blank_diff.csv
#       brain_tumor_vs_control_diff.csv
#       brain_tumor_vs_blank_diff.csv
#       rankNet_*.pdf
#       merged_circle_*.pdf
#     summary/
#       task_status.csv
#       cellchat_summary.md
#       error_log.txt
#
# 资源:
#   单器官单组 50-90 GB 内存, 30-60 min CPU
#   6 任务总计 3-6 h
#   服务器 1 TB 内存足够 (单任务峰值远低于阈值)
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(patchwork)
  library(NMF)
  library(ggalluvial)
  library(dplyr)
  library(ggplot2)
})

# dplyr namespace lock
select    <- dplyr::select
rename    <- dplyr::rename
filter    <- dplyr::filter
mutate    <- dplyr::mutate
arrange   <- dplyr::arrange
bind_rows <- dplyr::bind_rows

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 100 * 1024^3)  # 100 GB future limit

# ============================================================
# 0. 路径 & 共用模块
# ============================================================
PROJECT_ROOT <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
SCRIPT_DIR   <- file.path(PROJECT_ROOT, "scripts")
setwd(PROJECT_ROOT)

source(file.path(SCRIPT_DIR, "diff_common.R"))

OUT_DIR <- file.path(PROJECT_ROOT, "cellchat_results")
for (d in c("comparison", "summary")) {
  dir.create(file.path(OUT_DIR, d), recursive = TRUE, showWarnings = FALSE)
}

ERROR_LOG <- file.path(OUT_DIR, "summary", "error_log.txt")
writeLines(sprintf("== CellChat error log - %s ==",
                   format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ERROR_LOG)

write_error <- function(context, msg) {
  line <- sprintf("[%s] %s => %s",
                  format(Sys.time(), "%H:%M:%S"), context, msg)
  cat("  ⚠", line, "\n")
  cat(line, "\n", file = ERROR_LOG, append = TRUE)
}

# ============================================================
# 1. 参数
# ============================================================

# 内存防护
MIN_FREE_GB_TO_START <- 200
MIN_FREE_GB_ABORT    <- 100
MAX_R_PROCS_ALLOWED  <- 2

# CellChat 参数
COMM_PROB_TYPE  <- "triMean"   # 默认, 比 "truncatedMean" 更严格
RAW_USE         <- TRUE        # 用原始 counts 计算 (推荐)
POPULATION_SIZE <- TRUE        # 考虑细胞数量加权
MIN_CELLS_PER_TYPE <- 10       # 每个细类最少细胞 (低于此被剔除)

# 组定义 (复用 diff_common.R 的样本)
LIVER_SAMPLES_BY_GROUP <- list(
  Tumor   = c("liver3", "liver8", "gan56"),
  Blank   = c("bla_liver"),
  Control = c("con_liver")
)
BRAIN_SAMPLES_BY_GROUP <- list(
  Tumor   = c("brain3", "brain8", "nao56"),
  Blank   = c("bla_brain"),
  Control = c("con_brain")
)

GROUPS <- c("Tumor", "Blank", "Control")
ORGANS <- c("liver", "brain")

cat("========================================\n")
cat("CellChat v4 cell-cell communication analysis\n")
cat("Start:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================\n\n")

# ============================================================
# 1.5 安全检查函数 (复用拟时序脚本逻辑)
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
    n <- suppressWarnings(as.integer(out))
    if (is.na(n)) 0 else n
  }, error = function(e) 0)
}

preflight_check <- function() {
  cat("========== Pre-run safety check ==========\n")
  n_r <- get_r_proc_count()
  cat(sprintf("  Residual R processes: %d (limit: %d)\n",
              n_r, MAX_R_PROCS_ALLOWED))
  if (n_r > MAX_R_PROCS_ALLOWED) {
    stop(sprintf("ABORT: %d R 进程残留. pkill -9 -u $USER R", n_r))
  }
  free_gb <- get_free_mem_gb()
  cat(sprintf("  Available memory: %.0f GB (min %d GB)\n",
              if (is.na(free_gb)) 0 else free_gb, MIN_FREE_GB_TO_START))
  if (is.na(free_gb) || free_gb < MIN_FREE_GB_TO_START) {
    stop(sprintf("ABORT: 可用内存不足"))
  }
  cat("  ✓ Safety check passed\n\n")
}

runtime_check <- function(task_idx, task_key) {
  free_gb <- get_free_mem_gb()
  if (is.na(free_gb)) {
    cat("  ⚠ runtime_check: 无法获取内存状态\n")
    return(invisible(NULL))
  }
  cat(sprintf("  [runtime_check] task %d (%s) starting, free: %.0f GB\n",
              task_idx, task_key, free_gb))
  if (free_gb < MIN_FREE_GB_ABORT) {
    stop(sprintf("ABORT at task %d: free %.0f GB < %d GB",
                 task_idx, free_gb, MIN_FREE_GB_ABORT))
  }
  invisible(NULL)
}

preflight_check()

# ============================================================
# 1.6 加载 NAM 清洗后注释 (nam_clean_v2.R 产出)
# ============================================================
NAM_CLEAN_CSV <- file.path(PROJECT_ROOT, "nam_clean_results",
                            "v4_refined_clean_full.csv")
NAM_CLEAN_ANNO <- NULL
if (file.exists(NAM_CLEAN_CSV)) {
  cat(sprintf("\n========== Loading NAM clean annotation ==========\n"))
  NAM_CLEAN_ANNO <- read.csv(NAM_CLEAN_CSV, stringsAsFactors = FALSE)
  cat(sprintf("  Loaded %d cells from %s\n",
              nrow(NAM_CLEAN_ANNO), basename(NAM_CLEAN_CSV)))
  cat(sprintf("  v4 NAM cleaning impact:\n"))
  print(table(NAM_CLEAN_ANNO$v4_refined_clean[
    grepl("Macrophage_NAM", NAM_CLEAN_ANNO$v4_refined_clean)]))
} else {
  cat(sprintf("\n⚠ NAM clean CSV not found at: %s\n", NAM_CLEAN_CSV))
  cat("  Will use original v4_refined NAM labels (with hepatocyte contamination)\n")
}

# ============================================================
# 2. 加载 CellChat mouse 数据库
# ============================================================
cat("========== Loading CellChat mouse database ==========\n")

CellChatDB <- CellChatDB.mouse
cat(sprintf("  database total LR pairs: %d\n", nrow(CellChatDB$interaction)))

# 通常分为 Secreted Signaling, ECM-Receptor, Cell-Cell Contact
# 默认全部用 (用 'Secreted Signaling' 也可以但会丢失 ECM 相关)
CellChatDB.use <- CellChatDB
cat(sprintf("  using all categories: %s\n",
            paste(unique(CellChatDB$interaction$annotation), collapse = ", ")))
cat("\n")

# ============================================================
# 3. 加载样本 (按 organ 分段, 避免同时持有两个 organ)
# ============================================================

load_organ_samples <- function(organ) {
  grp_map <- if (organ == "liver") LIVER_SAMPLES_BY_GROUP else BRAIN_SAMPLES_BY_GROUP
  samples <- unlist(grp_map)

  cat(sprintf("\n[%s] loading %d samples\n", organ, length(samples)))

  objs <- list()
  for (s in samples) {
    fp <- resolve_rds_path(s, "v4")
    if (!file.exists(fp)) {
      write_error(s, paste0("missing v4 RDS: ", fp))
      next
    }
    o <- tryCatch({
      x <- readRDS(fp)
      if (inherits(x, "Seurat")) x <- tryCatch(JoinLayers(x), error = function(e) x)
      x <- RenameCells(x, add.cell.id = s)
      x$sample <- s
      # 找 group
      grp <- NA
      for (g in names(grp_map)) {
        if (s %in% grp_map[[g]]) { grp <- g; break }
      }
      x$group <- grp
      x$organ <- organ
      x
    }, error = function(e) { write_error(s, e$message); NULL })
    if (!is.null(o)) {
      objs[[s]] <- o
      cat(sprintf("  ✓ %-12s  cells=%6d  group=%s\n", s, ncol(o), o$group[1]))
    }
  }

  if (length(objs) < 2) {
    write_error(organ, "fewer than 2 samples loaded")
    return(NULL)
  }

  cat(sprintf("  merging %d %s samples...\n", length(objs), organ))
  merged <- merge(objs[[1]], objs[-1])
  merged <- tryCatch(JoinLayers(merged), error = function(e) merged)
  cat(sprintf("  [%s] merged cells: %d\n", organ, ncol(merged)))

  rm(objs); gc(verbose = FALSE)
  return(merged)
}

# ============================================================
# 4. 将 v4 注释展开成实际 CellChat 用的标签
# ============================================================
# diff_common.R 提供 LIVER_CELL_TYPES / BRAIN_CELL_TYPES 这样的 list
# 每个 virtual cell type 对应一组 v4_refined 实际 label
# CellChat 需要一个扁平的 label 向量, 所以我们用 virtual_ct 作为 CellChat label
# (对应原始细胞 → 虚拟细类的映射)

build_celltype_labels <- function(seurat_obj, organ) {
  ct_map <- if (organ == "liver") LIVER_CELL_TYPES else BRAIN_CELL_TYPES

  # 探测注释列
  ct_col <- NULL
  for (candidate in c("celltype_refined", "v4_refined", "final_celltype")) {
    if (candidate %in% colnames(seurat_obj@meta.data)) {
      ct_col <- candidate; break
    }
  }
  if (is.null(ct_col)) {
    stop("没找到细胞类型列")
  }

  raw_labels <- as.character(seurat_obj@meta.data[[ct_col]])
  cell_ids   <- colnames(seurat_obj)

  # === NAM 清洗 join (仅 liver) ===
  # 用 nam_clean_v2.R 输出的 v4_refined_clean 覆盖原 v4 的 Macrophage_NAM 标签
  # - 真 NAM (72) → 保留为 "Macrophage_NAM" (映射到虚拟细类 Macrophage_NAM)
  # - 假 NAM (226) → 改成 "Macrophage_NAM_misannotated" → 自动落入 "Other" 而被排除
  if (organ == "liver" && !is.null(NAM_CLEAN_ANNO)) {
    n_before <- sum(raw_labels == "Macrophage_NAM", na.rm = TRUE)
    match_idx <- match(cell_ids, NAM_CLEAN_ANNO$cell_id)
    has_clean <- !is.na(match_idx)
    n_matched <- sum(has_clean)
    cat(sprintf("    [NAM clean join] %d / %d cells matched in clean CSV\n",
                n_matched, length(cell_ids)))

    # 只对 v4 NAM 细胞应用 clean (其他细胞标签不动)
    is_v4_nam <- raw_labels == "Macrophage_NAM"
    is_v4_nam[is.na(is_v4_nam)] <- FALSE
    apply_mask <- is_v4_nam & has_clean
    n_apply <- sum(apply_mask)

    if (n_apply > 0) {
      new_label_for_nam <- NAM_CLEAN_ANNO$v4_refined_clean[match_idx[apply_mask]]
      # Map clean label back to v4-compatible:
      #   Macrophage_NAM_clean       → Macrophage_NAM (真 NAM, 保留虚拟映射)
      #   Macrophage_NAM_misannotated→ Macrophage_NAM_misannotated (新 label, 落入 Other)
      new_label_for_nam[new_label_for_nam == "Macrophage_NAM_clean"] <-
        "Macrophage_NAM"
      raw_labels[apply_mask] <- new_label_for_nam

      n_clean_kept <- sum(raw_labels == "Macrophage_NAM", na.rm = TRUE)
      n_misannot   <- sum(raw_labels == "Macrophage_NAM_misannotated", na.rm = TRUE)
      cat(sprintf("    [NAM clean join] before: %d v4 NAM\n", n_before))
      cat(sprintf("    [NAM clean join] after:  %d clean NAM kept, %d misannotated → Other\n",
                  n_clean_kept, n_misannot))
    }
  }

  # 反向映射: raw label → virtual label
  virtual_labels <- rep("Other", length(raw_labels))
  for (vct in names(ct_map)) {
    matching_raw <- ct_map[[vct]]
    virtual_labels[raw_labels %in% matching_raw] <- vct
  }

  # 排除 All_Cells (这个虚拟类是"全部细胞", 不参与 CellChat)
  if (any(virtual_labels == "All_Cells")) {
    virtual_labels[virtual_labels == "All_Cells"] <- "Other"
  }

  return(virtual_labels)
}

# ============================================================
# 5. 单组 CellChat 分析的核心函数
# ============================================================
run_one_cellchat <- function(seurat_obj, organ, group_name, task_idx) {
  task_key <- paste0(organ, "_", group_name)
  task_dir <- file.path(OUT_DIR, task_key)
  plot_dir <- file.path(task_dir, "plots")
  dir.create(task_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  result <- list(
    organ = organ, group = group_name, key = task_key,
    status = "error", n_cells = NA, n_celltypes = NA,
    n_lr_pairs_sig = NA, n_pathways = NA, elapsed_min = NA
  )
  t0 <- proc.time()

  cat(sprintf("\n### Task %d: %s ###\n", task_idx, task_key))
  runtime_check(task_idx, task_key)

  tryCatch({
    # 5.1 抽出该组细胞
    cell_idx <- which(seurat_obj$group == group_name)
    if (length(cell_idx) < 100) {
      stop(sprintf("group %s 细胞数 %d < 100", group_name, length(cell_idx)))
    }
    sub <- seurat_obj[, cell_idx]
    sub <- tryCatch(JoinLayers(sub), error = function(e) sub)
    cat(sprintf("  [%s] n_cells=%d\n", task_key, ncol(sub)))

    # 5.2 构造细类标签
    virtual_labels <- build_celltype_labels(sub, organ)
    sub$cellchat_label <- virtual_labels

    label_counts <- table(virtual_labels)
    cat(sprintf("  cell type counts:\n"))
    for (lab in names(label_counts)) {
      cat(sprintf("    %-25s : %d\n", lab, label_counts[[lab]]))
    }

    # 排除 cell_count < MIN_CELLS_PER_TYPE 的细类
    valid_types <- names(label_counts)[label_counts >= MIN_CELLS_PER_TYPE]
    valid_types <- setdiff(valid_types, "Other")  # 也排除 Other
    if (length(valid_types) < 3) {
      stop(sprintf("有效细类 < 3 (%d), 无法做 CellChat",
                   length(valid_types)))
    }

    keep_idx <- which(sub$cellchat_label %in% valid_types)
    sub <- sub[, keep_idx]
    cat(sprintf("  after filtering: %d cells, %d cell types\n",
                ncol(sub), length(valid_types)))
    result$n_cells <- ncol(sub)
    result$n_celltypes <- length(valid_types)

    # 5.3 准备 CellChat 输入
    # data.input: log-normalized matrix (CellChat 需要 normalized data)
    if (length(sub@assays$RNA@layers) > 0 &&
        "data" %in% names(sub@assays$RNA@layers) &&
        sum(sub@assays$RNA@layers$data) > 0) {
      data_input <- GetAssayData(sub, assay = "RNA", layer = "data")
    } else {
      cat("  Normalizing data (RNA layer 'data' not found)...\n")
      sub <- NormalizeData(sub, verbose = FALSE)
      data_input <- GetAssayData(sub, assay = "RNA", layer = "data")
    }

    meta <- data.frame(
      labels = sub$cellchat_label,
      sample = sub$sample,
      row.names = colnames(sub),
      stringsAsFactors = FALSE
    )

    # 5.4 创建 CellChat 对象
    cat(sprintf("  Creating CellChat object...\n"))
    cellchat <- createCellChat(
      object   = data_input,
      meta     = meta,
      group.by = "labels"
    )
    cellchat@DB <- CellChatDB.use

    # 5.5 预处理: 子集 + 识别过表达基因 / interactions
    cat("  subsetData → identifyOverExpressedGenes → identifyOverExpressedInteractions\n")
    cellchat <- subsetData(cellchat)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)

    n_over_genes <- length(cellchat@var.features$features)
    n_over_inter <- nrow(cellchat@LR$LRsig)
    cat(sprintf("  over-expressed genes: %d, over-expressed LR: %d\n",
                n_over_genes, n_over_inter))

    # 5.6 计算通讯概率 (这是最贵的一步)
    cat(sprintf("  computeCommunProb (type=%s)...\n", COMM_PROB_TYPE))
    t_comm <- proc.time()
    cellchat <- computeCommunProb(
      cellchat,
      type            = COMM_PROB_TYPE,
      raw.use         = RAW_USE,
      population.size = POPULATION_SIZE
    )
    cat(sprintf("  computeCommunProb done in %.1f min\n",
                (proc.time() - t_comm)["elapsed"] / 60))

    # 过滤掉 cell counts 过少的细类的通讯
    cellchat <- filterCommunication(cellchat, min.cells = MIN_CELLS_PER_TYPE)

    # 5.7 通路水平计算
    cat("  computeCommunProbPathway → aggregateNet\n")
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)

    # 5.8 提取结果
    df_lr <- subsetCommunication(cellchat)  # LR-level
    df_pathway <- subsetCommunication(cellchat, slot.name = "netP")  # pathway-level

    result$n_lr_pairs_sig <- nrow(df_lr)
    result$n_pathways <- length(unique(df_pathway$pathway_name))

    cat(sprintf("  ✓ significant LR pairs: %d, pathways: %d\n",
                result$n_lr_pairs_sig, result$n_pathways))

    # 5.9 保存数据
    saveRDS(cellchat, file.path(task_dir, "cellchat_obj.rds"))
    write.csv(df_lr, file.path(task_dir, "lr_pairs_filtered.csv"),
              row.names = FALSE)
    write.csv(df_pathway, file.path(task_dir, "pathway_summary.csv"),
              row.names = FALSE)

    # interaction count / weight 矩阵
    write.csv(cellchat@net$count,
              file.path(task_dir, "interactions_count.csv"))
    write.csv(cellchat@net$weight,
              file.path(task_dir, "interactions_weight.csv"))

    # 5.10 可视化
    tryCatch({
      pdf(file.path(plot_dir, "circle_count.pdf"), width = 9, height = 9)
      par(xpd=TRUE)
      netVisual_circle(
        cellchat@net$count, vertex.weight = as.numeric(table(cellchat@idents)),
        weight.scale = TRUE, label.edge = FALSE,
        title.name = sprintf("%s · interaction count", task_key)
      )
      dev.off()
    }, error = function(e) write_error(paste0(task_key, "_circle_count"), e$message))

    tryCatch({
      pdf(file.path(plot_dir, "circle_weight.pdf"), width = 9, height = 9)
      par(xpd=TRUE)
      netVisual_circle(
        cellchat@net$weight, vertex.weight = as.numeric(table(cellchat@idents)),
        weight.scale = TRUE, label.edge = FALSE,
        title.name = sprintf("%s · interaction weight", task_key)
      )
      dev.off()
    }, error = function(e) write_error(paste0(task_key, "_circle_weight"), e$message))

    # heatmap of pathways
    if (length(cellchat@netP$pathways) > 0) {
      tryCatch({
        ht <- netAnalysis_signalingRole_heatmap(
          cellchat,
          pattern = "outgoing",
          height = max(6, length(cellchat@netP$pathways) * 0.2),
          width  = 8,
          font.size = 7
        )
        pdf(file.path(plot_dir, "heatmap_pathway_outgoing.pdf"),
            width = 9, height = max(6, length(cellchat@netP$pathways) * 0.2))
        print(ht)
        dev.off()
      }, error = function(e) write_error(paste0(task_key, "_ht_out"), e$message))

      tryCatch({
        ht <- netAnalysis_signalingRole_heatmap(
          cellchat,
          pattern = "incoming",
          height = max(6, length(cellchat@netP$pathways) * 0.2),
          width  = 8,
          font.size = 7
        )
        pdf(file.path(plot_dir, "heatmap_pathway_incoming.pdf"),
            width = 9, height = max(6, length(cellchat@netP$pathways) * 0.2))
        print(ht)
        dev.off()
      }, error = function(e) write_error(paste0(task_key, "_ht_in"), e$message))
    }

    # bubble plot of top LR pairs
    if (nrow(df_lr) > 0) {
      tryCatch({
        # top 50 by prob
        top_lr_pairs <- df_lr %>%
          dplyr::arrange(dplyr::desc(prob)) %>%
          dplyr::slice_head(n = 50)

        if (nrow(top_lr_pairs) > 0) {
          # netVisual_bubble 需要 sources.use / targets.use
          all_celltypes <- levels(cellchat@idents)
          p_bubble <- netVisual_bubble(
            cellchat,
            sources.use = all_celltypes,
            targets.use = all_celltypes,
            remove.isolate = TRUE,
            angle.x = 45,
            return.data = FALSE
          )
          ggsave(file.path(plot_dir, "bubble_lr.pdf"), p_bubble,
                 width = 14, height = 10, dpi = 200, limitsize = FALSE)
        }
      }, error = function(e) write_error(paste0(task_key, "_bubble"), e$message))
    }

    result$status <- "ok"
    result$elapsed_min <- round((proc.time() - t0)["elapsed"] / 60, 1)
    cat(sprintf("  ✓ DONE  %s  n_cells=%d  LR=%d  pathways=%d  (%.1f min)\n",
                task_key, result$n_cells,
                result$n_lr_pairs_sig, result$n_pathways,
                result$elapsed_min))

    # 立刻清理
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
# 6. 跑 6 个子任务: liver × 3 组 + brain × 3 组
# ============================================================

t_all_start <- proc.time()
all_results <- list()
task_idx <- 0

for (organ in ORGANS) {
  cat(sprintf("\n\n========== Loading %s samples ==========\n", organ))
  organ_obj <- load_organ_samples(organ)

  if (is.null(organ_obj)) {
    write_error(organ, "load_organ_samples returned NULL, skip")
    next
  }

  free_after_load <- get_free_mem_gb()
  cat(sprintf("  [%s loaded] free mem: %.0f GB\n",
              organ,
              if (is.na(free_after_load)) -1 else free_after_load))

  for (grp in GROUPS) {
    task_idx <- task_idx + 1
    res <- run_one_cellchat(organ_obj, organ, grp, task_idx)
    all_results[[paste0(organ, "_", grp)]] <- res

    free_after_task <- get_free_mem_gb()
    cat(sprintf("  [after task %d] free mem: %.0f GB\n",
                task_idx,
                if (is.na(free_after_task)) -1 else free_after_task))
  }

  # 释放 organ 对象, 准备下一个 organ
  cat(sprintf("\n-- Releasing %s merged object --\n", organ))
  rm(organ_obj)
  for (j in 1:5) gc(verbose = FALSE)

  free_after_release <- get_free_mem_gb()
  cat(sprintf("  [%s released] free mem: %.0f GB\n",
              organ,
              if (is.na(free_after_release)) -1 else free_after_release))
}

t_all_elapsed <- (proc.time() - t_all_start)["elapsed"]
cat(sprintf("\n\nAll 6 tasks elapsed: %.1f min\n", t_all_elapsed / 60))

# ============================================================
# 7. task_status.csv
# ============================================================
cat("\n========== Writing task_status.csv ==========\n")

ts_df <- dplyr::bind_rows(lapply(all_results, function(r) {
  data.frame(
    key            = r$key,
    organ          = r$organ,
    group          = r$group,
    status         = r$status,
    n_cells        = if (is.null(r$n_cells)) NA else r$n_cells,
    n_celltypes    = if (is.null(r$n_celltypes)) NA else r$n_celltypes,
    n_lr_pairs_sig = if (is.null(r$n_lr_pairs_sig)) NA else r$n_lr_pairs_sig,
    n_pathways     = if (is.null(r$n_pathways)) NA else r$n_pathways,
    elapsed_min    = if (is.null(r$elapsed_min)) NA else r$elapsed_min,
    stringsAsFactors = FALSE
  )
}))
write.csv(ts_df, file.path(OUT_DIR, "summary", "task_status.csv"),
          row.names = FALSE)
print(ts_df)

# ============================================================
# 8. mergeCellChat 跨组对比 (按 organ)
# ============================================================
cat("\n\n========== Running cross-group comparisons ==========\n")

run_comparison <- function(organ) {
  cat(sprintf("\n[%s comparison]\n", organ))

  # 加载 3 个 cellchat 对象
  cc_list <- list()
  for (grp in GROUPS) {
    fp <- file.path(OUT_DIR, paste0(organ, "_", grp), "cellchat_obj.rds")
    if (!file.exists(fp)) {
      write_error(paste0(organ, "_compare"), paste0("missing: ", fp))
      next
    }
    cc_list[[grp]] <- tryCatch(readRDS(fp),
                               error = function(e) { write_error(grp, e$message); NULL })
  }
  cc_list <- cc_list[!sapply(cc_list, is.null)]

  if (length(cc_list) < 2) {
    write_error(paste0(organ, "_compare"), "fewer than 2 valid CellChat objects")
    return(invisible(NULL))
  }

  # mergeCellChat
  merged <- tryCatch(
    mergeCellChat(cc_list, add.names = names(cc_list), cell.prefix = TRUE),
    error = function(e) {
      write_error(paste0(organ, "_merge"), e$message)
      NULL
    }
  )
  if (is.null(merged)) return(invisible(NULL))

  saveRDS(merged, file.path(OUT_DIR, "comparison",
                             paste0(organ, "_merged.rds")))

  # 8.1 总 interaction 数对比 bar plot
  tryCatch({
    p1 <- compareInteractions(merged, show.legend = FALSE,
                              group = seq_along(cc_list))
    p2 <- compareInteractions(merged, show.legend = FALSE,
                              group = seq_along(cc_list), measure = "weight")
    pp <- p1 + p2
    ggsave(file.path(OUT_DIR, "comparison",
                     paste0(organ, "_compare_interactions.pdf")),
           pp, width = 10, height = 5, dpi = 200)
  }, error = function(e) write_error(paste0(organ, "_compareInter"), e$message))

  # 8.2 rankNet (通路水平差异)
  tryCatch({
    p <- rankNet(merged, mode = "comparison",
                 stacked = TRUE, do.stat = TRUE)
    ggsave(file.path(OUT_DIR, "comparison",
                     paste0(organ, "_rankNet_stacked.pdf")),
           p, width = 8, height = 14, dpi = 200, limitsize = FALSE)
  }, error = function(e) write_error(paste0(organ, "_rankNet"), e$message))

  # 8.3 对每对组合做差异分析
  comparisons <- list(
    c("Tumor", "Control"),
    c("Tumor", "Blank"),
    c("Blank", "Control")
  )

  for (cmp in comparisons) {
    if (!all(cmp %in% names(cc_list))) next

    cmp_name <- paste0(cmp[1], "_vs_", cmp[2])
    cat(sprintf("  comparison: %s %s\n", organ, cmp_name))

    # 提取两组的 LR 通讯, 计算差异
    df1 <- tryCatch(subsetCommunication(cc_list[[cmp[1]]]),
                    error = function(e) NULL)
    df2 <- tryCatch(subsetCommunication(cc_list[[cmp[2]]]),
                    error = function(e) NULL)
    if (is.null(df1) || is.null(df2)) next

    # LR pair 唯一 key: source__target__interaction_name
    df1$lr_key <- paste(df1$source, df1$target, df1$interaction_name, sep = "__")
    df2$lr_key <- paste(df2$source, df2$target, df2$interaction_name, sep = "__")

    # diff: in tumor only / in control only / shared with prob diff
    only_in_1 <- setdiff(df1$lr_key, df2$lr_key)
    only_in_2 <- setdiff(df2$lr_key, df1$lr_key)
    shared    <- intersect(df1$lr_key, df2$lr_key)

    diff_df <- data.frame(
      lr_key = c(only_in_1, only_in_2, shared),
      group1 = cmp[1],
      group2 = cmp[2],
      category = c(
        rep(paste0("only_in_", cmp[1]), length(only_in_1)),
        rep(paste0("only_in_", cmp[2]), length(only_in_2)),
        rep("shared", length(shared))
      ),
      stringsAsFactors = FALSE
    )

    # 加上 prob 信息
    df1_lookup <- df1[match(diff_df$lr_key, df1$lr_key), c("source", "target",
                                                           "interaction_name",
                                                           "pathway_name", "prob")]
    df2_lookup <- df2[match(diff_df$lr_key, df2$lr_key), c("prob")]
    diff_df$source <- df1_lookup$source
    diff_df$target <- df1_lookup$target
    diff_df$interaction_name <- df1_lookup$interaction_name
    diff_df$pathway_name <- df1_lookup$pathway_name
    diff_df$prob_group1 <- df1_lookup$prob
    diff_df$prob_group2 <- df2_lookup
    diff_df$prob_diff <- diff_df$prob_group1 - diff_df$prob_group2

    # 对于 only_in_1 / only_in_2, source/target/interaction 用对应组的
    only2_idx <- which(diff_df$category == paste0("only_in_", cmp[2]))
    if (length(only2_idx) > 0) {
      df2_for_only2 <- df2[match(diff_df$lr_key[only2_idx], df2$lr_key),
                           c("source", "target", "interaction_name",
                             "pathway_name", "prob")]
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
                        paste0(organ, "_", cmp_name, "_diff.csv")),
              row.names = FALSE)
    cat(sprintf("    %s: total %d (only_%s=%d, only_%s=%d, shared=%d)\n",
                cmp_name, nrow(diff_df),
                cmp[1], length(only_in_1),
                cmp[2], length(only_in_2),
                length(shared)))

    # rankNet for this pair specifically
    pair_idx <- match(cmp, names(cc_list))
    if (length(pair_idx) == 2 && all(!is.na(pair_idx))) {
      tryCatch({
        merged_pair <- mergeCellChat(cc_list[pair_idx],
                                     add.names = cmp,
                                     cell.prefix = TRUE)
        p <- rankNet(merged_pair, mode = "comparison",
                     stacked = TRUE, do.stat = TRUE,
                     comparison = c(1, 2))
        ggsave(file.path(OUT_DIR, "comparison",
                         paste0(organ, "_rankNet_", cmp_name, ".pdf")),
               p, width = 8, height = 12, dpi = 200, limitsize = FALSE)
      }, error = function(e) write_error(paste0(organ, "_rank_", cmp_name), e$message))
    }
  }

  rm(cc_list, merged); for (j in 1:3) gc(verbose = FALSE)
  invisible(NULL)
}

run_comparison("liver")
run_comparison("brain")

# ============================================================
# 9. 写 summary md
# ============================================================
cat("\n========== Writing cellchat_summary.md ==========\n")

summary_lines <- c(
  "# CellChat 细胞通讯分析总结",
  "",
  sprintf("**生成时间:** %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  sprintf("**总耗时:** %.1f min", t_all_elapsed / 60),
  "",
  "## 方法",
  "",
  "- **工具:** CellChat (mouse database, all categories)",
  "- **设计:** 按 organ × group 拆成 6 个子任务独立运行",
  "  - liver × {Tumor, Blank, Control}",
  "  - brain × {Tumor, Blank, Control}",
  sprintf("- **通讯概率类型:** %s", COMM_PROB_TYPE),
  sprintf("- **最少细胞数 / 细类:** %d", MIN_CELLS_PER_TYPE),
  "- **细胞类型来源:** v4_refined → 虚拟细类映射 (复用 diff_common.R 的 LIVER_CELL_TYPES / BRAIN_CELL_TYPES)",
  sprintf("- **NAM 清洗:** %s",
          if (!is.null(NAM_CLEAN_ANNO))
            sprintf("YES (用 %s, 真 NAM 72 个保留, 226 个 hepatocyte 污染排除)",
                    basename(NAM_CLEAN_CSV))
          else "NO (使用原 v4 NAM, 包含 76%% 肝细胞污染 - 不推荐)"),
  "- **后续对比:** mergeCellChat → compareInteractions / rankNet / 自定义 LR diff",
  "",
  "## 6 个子任务状态",
  "",
  "| 任务 | 细胞数 | 细类数 | LR pairs (sig) | Pathways | 时间 (min) |",
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
  "## 关键产出",
  "",
  "### 单组结果 (`<organ>_<group>/`)",
  "- `cellchat_obj.rds` — CellChat 对象",
  "- `lr_pairs_filtered.csv` — 显著 LR 配对 (source, target, interaction_name, pathway_name, prob, pval)",
  "- `pathway_summary.csv` — 通路水平汇总",
  "- `interactions_count.csv` / `interactions_weight.csv` — 细胞间交互数 / 强度矩阵",
  "- `plots/circle_count.pdf` / `plots/circle_weight.pdf` — 细胞间通讯网络环图",
  "- `plots/heatmap_pathway_outgoing.pdf` / `_incoming.pdf` — 通路 outgoing/incoming heatmap",
  "- `plots/bubble_lr.pdf` — top LR 配对的 bubble plot",
  "",
  "### 跨组对比 (`comparison/`)",
  "- `<organ>_merged.rds` — mergeCellChat 合并对象",
  "- `<organ>_compare_interactions.pdf` — 总交互数对比 bar plot",
  "- `<organ>_rankNet_stacked.pdf` — 通路差异 rankNet (stacked)",
  "- `<organ>_<grp1>_vs_<grp2>_diff.csv` — LR 配对差异表 (only_in_X / shared / prob_diff)",
  "- `<organ>_rankNet_<cmp>.pdf` — 各对组合的 rankNet",
  "",
  "## 关键问题待回答",
  "",
  "1. **轴 3 验证**: HSC 是否通过 TGFB1-TGFBR1/2 → Microglia/Astrocyte 的配对在 tumor 组特异性激活",
  "2. **NAM 通讯角色**: Macrophage_NAM 在 liver tumor 组的 outgoing 通讯是否以 IL6/CXCL/TNF 为主",
  "3. **轴 1 IFN 共激活**: 脑内 IFN-γ 信号的 source 是哪些细胞 (Microglia / Endothelial / NK)",
  "4. **HSC 自分泌环**: HSC 自身是否通过 TGFb 自分泌强化激活",
  "",
  "## 局限",
  "",
  "- CellChat 假设细胞在同一组织内空间接触, **不能直接验证肝→脑跨器官通讯**",
  "- Blank 组样本数 1, Control 组样本数 1, 单样本带来的 bias 不可避免",
  "- 细类标签依赖 v4_refined 注释的准确性 (NAM 内部异质性等已知问题)"
)

writeLines(summary_lines, file.path(OUT_DIR, "summary", "cellchat_summary.md"))

cat("\n========================================\n")
cat("CellChat analysis complete\n")
cat(sprintf("  Output dir : %s\n", OUT_DIR))
cat(sprintf("  Summary    : %s\n", file.path(OUT_DIR, "summary", "cellchat_summary.md")))
cat(sprintf("  Total time : %.1f min\n", t_all_elapsed / 60))
cat("========================================\n")
