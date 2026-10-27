#!/usr/bin/env Rscript
# ============================================================
# diff_v4_core.R  --  阶段 1: 鼠源 DEG 核心分析 (无 viz/无 enrichment)
# ============================================================
# 用法:
#   cd ~/SGS10_analysis
#   Rscript scripts/diff_v4_core.R 2>&1 | tee scripts/logs/diff_core.log
#
# 设计:
#   - 只做 DEG + 写 CSV, 不做火山图/热图/富集
#     (失败率最低, 后续 viz/enrichment 脚本读 CSV)
#   - 三个 annotation 版本并行: v4_refined, v4_broad, v32
#     v4_* 读 <sample>/mouse/<sample>_mouse|manual_annotated_v4.rds
#     v32  读 <sample>/<sample>_manual_annotated.rds (顶层旧文件)
#   - 鼠源 only. PDX 肝样本只用 mouse RDS, 人源留给后续
#   - 6 个对比 (肝 3 + 脑 3)
#   - 并行: 32 并发 × 4 核/任务 = 128 核
#
# 输出:
#   diff_stage1_results/
#     deg_v4_refined/{liver,brain}/<comp>_<virtual_ct>.csv
#     deg_v4_broad/{liver,brain}/<comp>_<virtual_ct>.csv
#     deg_v32/{liver,brain}/<comp>_<virtual_ct>.csv
#     summary/
#       task_status.csv       每个任务 status/n_cells_g1/n_cells_g2/up/down
#       deg_summary.csv       Up/Down 计数 + Il6st 专项追踪
#       error_log.txt
# ============================================================

suppressPackageStartupMessages({
  library(igraph)
  library(Seurat)
  library(presto)
  library(parallel)
  library(future)
  library(dplyr)
})

# 锁定 dplyr 动词防 namespace 冲突
select    <- dplyr::select
rename    <- dplyr::rename
filter    <- dplyr::filter
mutate    <- dplyr::mutate
arrange   <- dplyr::arrange
summarise <- dplyr::summarise
group_by  <- dplyr::group_by
bind_rows <- dplyr::bind_rows

# ============================================================
# 0. 路径 & 共用模块
# ============================================================
PROJECT_ROOT <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
SCRIPT_DIR   <- file.path(PROJECT_ROOT, "scripts")
setwd(PROJECT_ROOT)

source(file.path(SCRIPT_DIR, "diff_common.R"))

OUT_DIR <- file.path(PROJECT_ROOT, "diff_stage1_results")
for (d in c(
  "deg_v4_refined/liver", "deg_v4_refined/brain",
  "deg_v4_broad/liver",   "deg_v4_broad/brain",
  "deg_v32/liver",         "deg_v32/brain",
  "summary"
)) dir.create(file.path(OUT_DIR, d), recursive = TRUE, showWarnings = FALSE)

ERROR_LOG   <- file.path(OUT_DIR, "summary", "error_log.txt")
TASK_STATUS <- file.path(OUT_DIR, "summary", "task_status.csv")
DEG_SUMMARY <- file.path(OUT_DIR, "summary", "deg_summary.csv")

writeLines(sprintf(
  "========================================\nDEG core error log - %s\n========================================",
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ERROR_LOG)

write_error <- function(context, msg) {
  line <- sprintf("[%s] %s => %s",
                  format(Sys.time(), "%H:%M:%S"), context, msg)
  cat("  ⚠", line, "\n")
  cat(line, "\n", file = ERROR_LOG, append = TRUE)
}

# ============================================================
# 1. 参数 & 并行
# ============================================================
P <- DEFAULT_PARAMS
FC_THRESHOLD <- P$FC_THRESHOLD
PADJ_CUT     <- P$PADJ_CUT
MIN_PCT      <- P$MIN_PCT
MAX_CELLS    <- P$MAX_CELLS
MIN_CELLS    <- P$MIN_CELLS
SEED         <- P$RANDOM_SEED

CORES_PER_TASK   <- 4
CORES_TOTAL      <- 128
CONCURRENT_TASKS <- floor(CORES_TOTAL / CORES_PER_TASK)

cat("========================================\n")
cat("Stage 1: DEG Core (mouse cells only)\n")
cat("Start:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================\n")
cat(sprintf("Parallel: %d concurrent × %d cores = %d\n",
            CONCURRENT_TASKS, CORES_PER_TASK,
            CONCURRENT_TASKS * CORES_PER_TASK))
cat(sprintf("Params: FC>%.2f  padj<%.2f  min_pct>%.2f  max_cells=%d  min_cells=%d\n\n",
            FC_THRESHOLD, PADJ_CUT, MIN_PCT, MAX_CELLS, MIN_CELLS))

options(future.globals.maxSize = 128 * 1024^3)

# ============================================================
# 2. 加载所有 RDS (v4 + v32)
# ============================================================
cat("========== Loading RDS ==========\n")

load_all_samples <- function(version) {
  out <- list()
  for (s in ALL_SAMPLES) {
    fp <- resolve_rds_path(s, version)
    if (!file.exists(fp)) {
      write_error(paste0(s, "/", version), paste0("missing: ", fp))
      next
    }
    obj <- tryCatch({
      o <- readRDS(fp)
      if (inherits(o, "Seurat")) {
        o <- tryCatch(JoinLayers(o), error = function(e) o)
      }
      RenameCells(o, add.cell.id = s)
    }, error = function(e) {
      write_error(paste0(s, "/", version), e$message)
      NULL
    })
    if (!is.null(obj)) {
      out[[s]] <- obj
      cat(sprintf("  [%-3s] %-12s cells=%6d\n",
                  version, s, ncol(obj)))
    }
  }
  out
}

cat("-- v4 (refined + broad share the same RDS) --\n")
seurat_v4 <- load_all_samples("v4")
cat(sprintf("  loaded %d/%d\n\n", length(seurat_v4), length(ALL_SAMPLES)))

cat("-- v3.2 (legacy top-level RDS) --\n")
seurat_v32 <- load_all_samples("v32")
cat(sprintf("  loaded %d/%d\n\n", length(seurat_v32), length(ALL_SAMPLES)))

# ============================================================
# 3. Helpers: 提取细胞, 合并多样本, 二组合并
# ============================================================

# 根据 virtual_ct + version 从一个对象取细胞
get_cells_for_vct <- function(obj, virtual_ct, version, cell_type_map) {
  tryCatch(
    extract_cells_v2(obj, virtual_ct, version, cell_type_map),
    error = function(e) character(0)
  )
}

# 把多个样本合并成一个 "group 对象" (只保留指定 virtual_ct 的细胞)
merge_samples_to_group <- function(obj_list, sample_names, virtual_ct,
                                   version, cell_type_map, group_label) {
  subs <- list()
  for (s in sample_names) {
    obj <- obj_list[[s]]
    if (is.null(obj)) next
    cells <- get_cells_for_vct(obj, virtual_ct, version, cell_type_map)
    if (length(cells) < MIN_CELLS) next
    sub <- subset(obj, cells = cells)
    sub$group <- group_label
    sub$sample_id <- s
    subs[[s]] <- sub
  }
  if (length(subs) == 0) return(NULL)
  if (length(subs) == 1) return(subs[[1]])
  merged <- Reduce(function(a, b) merge(a, b), subs)
  tryCatch(JoinLayers(merged), error = function(e) merged)
}

# 从一个对比定义构建合并对象
build_merged_for_comparison <- function(comp, virtual_ct, obj_list,
                                        version, cell_type_map) {
  g1 <- merge_samples_to_group(obj_list, comp$g1_samples, virtual_ct,
                               version, cell_type_map, comp$g1_label)
  g2 <- merge_samples_to_group(obj_list, comp$g2_samples, virtual_ct,
                               version, cell_type_map, comp$g2_label)
  if (is.null(g1) || is.null(g2)) return(NULL)
  if (ncol(g1) < MIN_CELLS || ncol(g2) < MIN_CELLS) return(NULL)
  merged <- merge(g1, g2)
  tryCatch(JoinLayers(merged), error = function(e) merged)
}

# ============================================================
# 4. DEG 核心函数
# ============================================================
run_one_deg <- function(merged_obj, comp_name, group1, group2, virtual_ct,
                        fm_cores = CORES_PER_TASK) {
  n1 <- sum(merged_obj$group == group1)
  n2 <- sum(merged_obj$group == group2)
  if (n1 < MIN_CELLS || n2 < MIN_CELLS) {
    return(list(status = "insufficient_cells", n1 = n1, n2 = n2, res = NULL))
  }

  # downsample
  set.seed(SEED)
  cg1 <- colnames(merged_obj)[merged_obj$group == group1]
  cg2 <- colnames(merged_obj)[merged_obj$group == group2]
  n1_used <- length(cg1); n2_used <- length(cg2)
  if (n1_used > MAX_CELLS) { cg1 <- sample(cg1, MAX_CELLS); n1_used <- MAX_CELLS }
  if (n2_used > MAX_CELLS) { cg2 <- sample(cg2, MAX_CELLS); n2_used <- MAX_CELLS }

  obj_sub <- subset(merged_obj, cells = c(cg1, cg2))
  Idents(obj_sub) <- obj_sub$group

  plan("multicore", workers = fm_cores)
  res <- tryCatch(
    FindMarkers(obj_sub,
                ident.1 = group1, ident.2 = group2,
                test.use = "wilcox",
                min.pct = MIN_PCT,
                logfc.threshold = 0,
                only.pos = FALSE,
                verbose = FALSE),
    error = function(e) {
      stop(sprintf("FindMarkers failed [%s|%s]: %s",
                   comp_name, virtual_ct, e$message))
    }
  )
  plan("sequential")

  if (is.null(res) || nrow(res) == 0) {
    return(list(status = "no_deg", n1 = n1, n2 = n2,
                n1_used = n1_used, n2_used = n2_used, res = NULL))
  }

  res$gene <- rownames(res)
  res <- res %>%
    dplyr::select(gene, avg_log2FC, pct.1, pct.2, p_val, p_val_adj) %>%
    dplyr::rename(pct1 = pct.1, pct2 = pct.2) %>%
    dplyr::arrange(p_val_adj, dplyr::desc(abs(avg_log2FC)))

  res$sig <- ifelse(
    abs(res$avg_log2FC) > FC_THRESHOLD & res$p_val_adj < PADJ_CUT,
    ifelse(res$avg_log2FC > 0, "Up", "Down"),
    "NS"
  )

  list(
    status = "ok", n1 = n1, n2 = n2,
    n1_used = n1_used, n2_used = n2_used,
    res = res
  )
}

# ============================================================
# 5. 单任务完整执行 (mclapply 的 worker)
# ============================================================
# 每个 worker:
#   1. 通过 fork 共享的 seurat_v4/seurat_v32 查到 cell type map
#   2. 在 worker 内部做 subset -> merge -> JoinLayers (只在本 worker 消耗内存)
#   3. 跑 DEG
#   4. 写 CSV
#   5. 返回结果元数据 (不返回 merged 对象, 不返回 res 全部内容, 减少 fork 回传开销)
#   6. 函数返回时 worker 进程退出, 所有内存释放
run_one_task <- function(task) {
  # 子进程重新加载包 (mclapply fork 后包环境仍然共享, 但 future plan 可能被重置)
  suppressPackageStartupMessages({
    library(igraph); library(Seurat); library(presto)
    library(future); library(dplyr)
  })
  # 子进程内同样锁 dplyr
  select    <- dplyr::select
  rename    <- dplyr::rename
  filter    <- dplyr::filter
  mutate    <- dplyr::mutate
  arrange   <- dplyr::arrange
  summarise <- dplyr::summarise
  group_by  <- dplyr::group_by
  bind_rows <- dplyr::bind_rows

  result_meta <- list(
    key       = task$key,
    comp_name = task$comp_name,
    virtual_ct= task$virtual_ct,
    version   = task$version,
    organ     = task$organ,
    status    = "error",
    n1 = NA, n2 = NA, n1_used = NA, n2_used = NA,
    n_up = NA, n_down = NA,
    il6st_fc = NA, il6st_padj = NA, il6st_sig = "Not Found",
    il6st_pct1 = NA, il6st_pct2 = NA,
    il6_core_fc  = setNames(as.list(rep(NA, length(IL6_PATHWAY_CORE))),
                            IL6_PATHWAY_CORE),
    il6_core_sig = setNames(as.list(rep("Not Found", length(IL6_PATHWAY_CORE))),
                            IL6_PATHWAY_CORE)
  )

  tryCatch({
    # ---- Step 1: 在 worker 内构建 merged 对象 ----
    # obj_list 通过 fork 共享, 不需要传参
    obj_list <- if (task$version == "v32") seurat_v32 else seurat_v4

    merged <- build_merged_for_comparison(
      task$comp, task$virtual_ct, obj_list,
      task$version,
      if (task$organ == "liver") LIVER_CELL_TYPES else BRAIN_CELL_TYPES
    )
    if (is.null(merged)) {
      result_meta$status <- "insufficient_cells"
      return(result_meta)
    }

    # ---- Step 2: 跑 DEG ----
    out <- run_one_deg(merged, task$comp_name,
                       task$group1, task$group2, task$virtual_ct,
                       fm_cores = task$fm_cores)

    # 立即释放 merged
    rm(merged); gc(verbose = FALSE)

    result_meta$status  <- out$status
    result_meta$n1      <- out$n1
    result_meta$n2      <- out$n2
    if (!is.null(out$n1_used)) result_meta$n1_used <- out$n1_used
    if (!is.null(out$n2_used)) result_meta$n2_used <- out$n2_used

    # ---- Step 3: 写 CSV + 提取 Il6st 等元数据 ----
    if (out$status == "ok" && !is.null(out$res)) {
      write.csv(out$res,
                file.path(task$deg_dir, paste0(task$key, ".csv")),
                row.names = FALSE)

      result_meta$n_up   <- sum(out$res$sig == "Up",   na.rm = TRUE)
      result_meta$n_down <- sum(out$res$sig == "Down", na.rm = TRUE)

      # Il6st 专项
      il6_row <- out$res[out$res$gene == "Il6st", , drop = FALSE]
      if (nrow(il6_row) > 0) {
        result_meta$il6st_fc   <- round(il6_row$avg_log2FC[1], 4)
        result_meta$il6st_padj <- il6_row$p_val_adj[1]
        result_meta$il6st_sig  <- il6_row$sig[1]
        result_meta$il6st_pct1 <- il6_row$pct1[1]
        result_meta$il6st_pct2 <- il6_row$pct2[1]
      }

      # IL6 核心基因面板
      for (g in IL6_PATHWAY_CORE) {
        row <- out$res[out$res$gene == g, , drop = FALSE]
        if (nrow(row) > 0) {
          result_meta$il6_core_fc[[g]]  <- round(row$avg_log2FC[1], 3)
          result_meta$il6_core_sig[[g]] <- row$sig[1]
        }
      }

      # 释放 res (只保留已经提取的元数据)
      rm(out); gc(verbose = FALSE)
    }

    # 打印到 stderr 让主进程看到进度
    message(sprintf("  ✓ [%s|%s|%s] n1=%d n2=%d Up=%s Down=%s",
                    task$version, task$comp_name, task$virtual_ct,
                    result_meta$n1, result_meta$n2,
                    result_meta$n_up, result_meta$n_down))

    return(result_meta)

  }, error = function(e) {
    msg <- sprintf("[%s] task[%s|%s|%s] FAILED: %s",
                   format(Sys.time(), "%H:%M:%S"),
                   task$comp_name, task$virtual_ct, task$version, e$message)
    # worker 里不能写同一个文件 (race condition), 只 print 到 stderr
    message(msg)
    result_meta$status <- "error"
    return(result_meta)
  })
}

# ============================================================
# 6. 任务构建 (只记录元数据, 不实体化 merged 对象)
# ============================================================
# 每个 task 只保留:
#   - comp, virtual_ct, version, organ, key, deg_dir
# merge 延后到 run_one_task 里由每个 worker 进程各自做
# 这样主进程只保留原始 seurat_v4/seurat_v32, 不会爆内存
build_task_list <- function(comparisons, cell_type_map, organ,
                            version, deg_dir_root) {
  deg_dir <- file.path(deg_dir_root, organ)
  tasks <- list()
  for (comp in comparisons) {
    for (virtual_ct in names(cell_type_map)) {
      # 检查该版本是否支持这个 virtual_ct
      if (is.null(cell_type_map[[virtual_ct]][[version]])) next

      tasks[[length(tasks) + 1]] <- list(
        comp       = comp,           # 整个对比定义 (含 g1/g2 samples 和 labels)
        virtual_ct = virtual_ct,
        version    = version,
        organ      = organ,
        comp_name  = comp$comp,
        group1     = comp$g1_label,
        group2     = comp$g2_label,
        key        = paste0(comp$comp, "_", virtual_ct),
        fm_cores   = CORES_PER_TASK,
        deg_dir    = deg_dir
      )
    }
  }
  tasks
}

# ============================================================
# 7. 构建全部任务 (3 版本 × 肝+脑)
# ============================================================
cat("========== Building task list (metadata only) ==========\n")

VERSION_TO_DIR <- list(
  v4_refined = file.path(OUT_DIR, "deg_v4_refined"),
  v4_broad   = file.path(OUT_DIR, "deg_v4_broad"),
  v32        = file.path(OUT_DIR, "deg_v32")
)

all_tasks <- list()
for (version in names(VERSION_TO_DIR)) {
  deg_root <- VERSION_TO_DIR[[version]]

  liver_tasks <- build_task_list(LIVER_COMPARISONS, LIVER_CELL_TYPES,
                                 "liver", version, deg_root)
  brain_tasks <- build_task_list(BRAIN_COMPARISONS, BRAIN_CELL_TYPES,
                                 "brain", version, deg_root)

  cat(sprintf("  [%-10s] liver=%d brain=%d total=%d\n",
              version, length(liver_tasks), length(brain_tasks),
              length(liver_tasks) + length(brain_tasks)))
  all_tasks <- c(all_tasks, liver_tasks, brain_tasks)
}
cat(sprintf("  TOTAL tasks: %d (merge will happen inside each worker)\n\n",
            length(all_tasks)))

if (length(all_tasks) == 0) {
  stop("No tasks to run, check cell type mapping and data loading")
}

# ============================================================
# 8. 执行所有任务 (并行)
# ============================================================
cat("========== Running DEG (parallel) ==========\n")
plan("sequential")
t_start <- proc.time()

n_parallel <- min(CONCURRENT_TASKS, length(all_tasks))
cat(sprintf("  %d tasks, %d concurrent workers\n\n",
            length(all_tasks), n_parallel))

results_raw <- mclapply(
  all_tasks, run_one_task,
  mc.cores = n_parallel,
  mc.preschedule = FALSE,
  mc.cleanup = TRUE
)

t_elapsed <- proc.time() - t_start
cat(sprintf("\n  Elapsed: %.1f seconds (%.1f min)\n\n",
            t_elapsed["elapsed"], t_elapsed["elapsed"] / 60))

# ============================================================
# 9. 收集结果 + 写 task_status.csv
# ============================================================
cat("========== Collecting results ==========\n")

status_rows <- list()
summary_rows <- list()
n_ok <- 0; n_insufficient <- 0; n_no_deg <- 0; n_error <- 0

for (r in results_raw) {
  if (is.null(r) || !is.list(r) || is.null(r$key)) {
    n_error <- n_error + 1
    next
  }

  status_rows[[length(status_rows) + 1]] <- data.frame(
    key        = r$key,
    comp_name  = r$comp_name,
    virtual_ct = r$virtual_ct,
    version    = r$version,
    organ      = r$organ,
    status     = r$status,
    n1         = if (is.null(r$n1)) NA else r$n1,
    n2         = if (is.null(r$n2)) NA else r$n2,
    n1_used    = if (is.null(r$n1_used)) NA else r$n1_used,
    n2_used    = if (is.null(r$n2_used)) NA else r$n2_used,
    n_up       = if (is.null(r$n_up)) NA else r$n_up,
    n_down     = if (is.null(r$n_down)) NA else r$n_down,
    stringsAsFactors = FALSE
  )

  if (r$status == "ok") {
    n_ok <- n_ok + 1

    summary_row <- data.frame(
      version    = r$version,
      organ      = r$organ,
      key        = r$key,
      n_up       = r$n_up,
      n_down     = r$n_down,
      n_total    = r$n_up + r$n_down,
      Il6st_FC   = r$il6st_fc,
      Il6st_padj = if (is.na(r$il6st_padj)) NA else sprintf("%.2e", r$il6st_padj),
      Il6st_sig  = r$il6st_sig,
      Il6st_pct1 = if (is.na(r$il6st_pct1)) NA else sprintf("%.1f%%", r$il6st_pct1*100),
      Il6st_pct2 = if (is.na(r$il6st_pct2)) NA else sprintf("%.1f%%", r$il6st_pct2*100),
      stringsAsFactors = FALSE
    )
    # 拼接 7 个核心基因列
    for (g in IL6_PATHWAY_CORE) {
      summary_row[[paste0(g, "_FC")]]  <- r$il6_core_fc[[g]]
      summary_row[[paste0(g, "_sig")]] <- r$il6_core_sig[[g]]
    }
    summary_rows[[length(summary_rows) + 1]] <- summary_row

  } else if (r$status == "insufficient_cells") {
    n_insufficient <- n_insufficient + 1
  } else if (r$status == "no_deg") {
    n_no_deg <- n_no_deg + 1
  } else {
    n_error <- n_error + 1
  }
}

task_status_df <- dplyr::bind_rows(status_rows)
write.csv(task_status_df, TASK_STATUS, row.names = FALSE)

cat(sprintf("  ok                 : %d\n", n_ok))
cat(sprintf("  insufficient_cells : %d\n", n_insufficient))
cat(sprintf("  no_deg             : %d\n", n_no_deg))
cat(sprintf("  error              : %d\n", n_error))
cat(sprintf("  task_status.csv -> %s\n\n", TASK_STATUS))

# ============================================================
# 10. 生成 DEG summary (含 Il6st 专项追踪)
# ============================================================
cat("========== Generating DEG summary ==========\n")

if (length(summary_rows) > 0) {
  deg_summary_df <- dplyr::bind_rows(summary_rows) %>%
    dplyr::arrange(version, organ, key)
  write.csv(deg_summary_df, DEG_SUMMARY, row.names = FALSE)
  cat(sprintf("  DEG summary: %d rows -> %s\n", nrow(deg_summary_df), DEG_SUMMARY))

  # 打印 Il6st 专项 (仅打印 Il6st 被检测到的行)
  cat("\n【Il6st 专项追踪】\n")
  cat("──────────────────────────────────────────────────────────\n")
  il6 <- deg_summary_df %>%
    dplyr::filter(!is.na(Il6st_FC), Il6st_sig != "Not Found") %>%
    dplyr::select(version, organ, key, Il6st_FC, Il6st_padj, Il6st_sig,
                  Il6st_pct1, Il6st_pct2)
  if (nrow(il6) > 0) {
    print(il6, row.names = FALSE)
  } else {
    cat("  (Il6st 未在任何成功 DEG 中检测到显著性)\n")
  }
  cat("──────────────────────────────────────────────────────────\n")
} else {
  cat("  ⚠ No ok results to summarize\n")
}

# ============================================================
# 11. 完成
# ============================================================
cat("\n========================================\n")
cat("Stage 1 DEG core complete\n")
cat(sprintf("  Total tasks  : %d\n", length(all_tasks)))
cat(sprintf("  Success (ok) : %d\n", n_ok))
cat(sprintf("  Elapsed      : %.1f min\n", t_elapsed["elapsed"] / 60))
cat(sprintf("  Output dir   : %s\n", OUT_DIR))
cat(sprintf("  Task status  : %s\n", TASK_STATUS))
cat(sprintf("  DEG summary  : %s\n", DEG_SUMMARY))
cat(sprintf("  Error log    : %s\n", ERROR_LOG))
cat("========================================\n")
