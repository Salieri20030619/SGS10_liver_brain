#!/usr/bin/env Rscript
# ============================================================
# diff_v4_paired.R  --  阶段 3: 跨器官配对分析
# ============================================================
# 用法:
#   cd ~/SGS10_analysis
#   Rscript scripts/diff_v4_paired.R 2>&1 | tee scripts/logs/diff_paired.log
#
# 前置: Stage 1 (diff_v4_core.R) 已完成
#
# 配对假设: 肝癌通过循环影响脑 -> 同一只 PDX 鼠内部, 肝响应越强脑响应也越强
# 配对样本:
#   mouse3  = liver3 <-> brain3
#   mouse8  = liver8 <-> brain8
#   mouse56 = gan56  <-> nao56
#
# 4 个子分析:
#   1. Per-sample DEG: 每只 PDX vs 对应器官的 con (n=1 对照)
#      13 virtual_ct x 3 mice x 2 organs = 39 基础任务, 实际约 39-42 个
#   2. Overlap 分析: 每只鼠的肝 DEG <-> 脑 DEG gene symbol 重叠 + hypergeometric test
#   3. FC 相关性: 机制基因在肝脑 FC 散点图 (3 张)
#   4. 响应强度 ranking: 3 只鼠的肝和脑响应强度排序一致性
#
# 输出:
#   diff_stage3_results/
#     per_sample_deg/{liver,brain}/<mouse>_<virtual_ct>.csv
#     overlap/mouse<id>_overlap_summary.csv
#     fc_correlation/mouse<id>_mechanism_scatter.{pdf,png}
#     response_ranking/
#       per_sample_response_scores.csv
#       ranking_plot.{pdf,png}
#     summary/
#       task_status.csv
#       paired_summary.md
#       error_log.txt
# ============================================================

suppressPackageStartupMessages({
  library(igraph)
  library(Seurat)
  library(presto)
  library(parallel)
  library(future)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(tidyr)
})

# dplyr namespace lock
select    <- dplyr::select
rename    <- dplyr::rename
filter    <- dplyr::filter
mutate    <- dplyr::mutate
arrange   <- dplyr::arrange
summarise <- dplyr::summarise
group_by  <- dplyr::group_by
bind_rows <- dplyr::bind_rows
pull      <- dplyr::pull

# ============================================================
# 0. 路径 & 共用模块
# ============================================================
PROJECT_ROOT <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
SCRIPT_DIR   <- file.path(PROJECT_ROOT, "scripts")
setwd(PROJECT_ROOT)

source(file.path(SCRIPT_DIR, "diff_common.R"))

OUT_DIR <- file.path(PROJECT_ROOT, "diff_stage3_results")
for (d in c(
  "per_sample_deg/liver", "per_sample_deg/brain",
  "overlap",
  "fc_correlation",
  "response_ranking",
  "summary"
)) dir.create(file.path(OUT_DIR, d), recursive = TRUE, showWarnings = FALSE)

ERROR_LOG   <- file.path(OUT_DIR, "summary", "error_log.txt")
TASK_STATUS <- file.path(OUT_DIR, "summary", "task_status.csv")
PAIRED_MD   <- file.path(OUT_DIR, "summary", "paired_summary.md")

writeLines(sprintf("== paired analysis error log - %s ==",
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
cat("Stage 3: Cross-organ paired analysis\n")
cat("Start:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================\n\n")

options(future.globals.maxSize = 128 * 1024^3)

# ============================================================
# 2. 配对定义
# ============================================================
PAIRED_MICE <- list(
  list(id = "mouse3",  liver = "liver3", brain = "brain3"),
  list(id = "mouse8",  liver = "liver8", brain = "brain8"),
  list(id = "mouse56", liver = "gan56",  brain = "nao56")
)

# 对照样本 (n=1 每器官, 和 Stage 1 保持一致)
CONTROL_LIVER <- "con_liver"
CONTROL_BRAIN <- "con_brain"

# ============================================================
# 3. 加载 v4 RDS (只读, worker 通过 fork 共享)
# ============================================================
cat("========== Loading v4 RDS ==========\n")

load_v4 <- function() {
  out <- list()
  for (s in ALL_SAMPLES) {
    fp <- resolve_rds_path(s, "v4")
    if (!file.exists(fp)) { write_error(s, paste0("missing: ", fp)); next }
    obj <- tryCatch({
      o <- readRDS(fp)
      if (inherits(o, "Seurat")) o <- tryCatch(JoinLayers(o), error = function(e) o)
      RenameCells(o, add.cell.id = s)
    }, error = function(e) { write_error(s, e$message); NULL })
    if (!is.null(obj)) {
      out[[s]] <- obj
      cat(sprintf("  [v4] %-12s cells=%6d\n", s, ncol(obj)))
    }
  }
  out
}

seurat_v4 <- load_v4()
cat(sprintf("  loaded %d/%d\n\n", length(seurat_v4), length(ALL_SAMPLES)))

# ============================================================
# 4. Helper: 单样本 vs 单样本 的 merge (只取指定 virtual_ct 的细胞)
# ============================================================
build_single_pair_merged <- function(pdx_sample, ctrl_sample, virtual_ct,
                                     cell_type_map) {
  if (is.null(seurat_v4[[pdx_sample]]) || is.null(seurat_v4[[ctrl_sample]])) {
    return(NULL)
  }
  cells_pdx <- tryCatch(
    extract_cells_v2(seurat_v4[[pdx_sample]], virtual_ct, "v4_refined", cell_type_map),
    error = function(e) character(0)
  )
  cells_ctrl <- tryCatch(
    extract_cells_v2(seurat_v4[[ctrl_sample]], virtual_ct, "v4_refined", cell_type_map),
    error = function(e) character(0)
  )
  if (length(cells_pdx) < MIN_CELLS || length(cells_ctrl) < MIN_CELLS) return(NULL)

  sub_pdx <- subset(seurat_v4[[pdx_sample]], cells = cells_pdx)
  sub_pdx$group <- "PDX"
  sub_ctrl <- subset(seurat_v4[[ctrl_sample]], cells = cells_ctrl)
  sub_ctrl$group <- "Control"

  merged <- merge(sub_pdx, sub_ctrl)
  tryCatch(JoinLayers(merged), error = function(e) merged)
}

# ============================================================
# 5. DEG 核心 (复用 Stage 1 逻辑)
# ============================================================
run_one_deg <- function(merged_obj, key, fm_cores = CORES_PER_TASK) {
  n1 <- sum(merged_obj$group == "PDX")
  n2 <- sum(merged_obj$group == "Control")
  if (n1 < MIN_CELLS || n2 < MIN_CELLS) {
    return(list(status = "insufficient_cells", n1 = n1, n2 = n2, res = NULL))
  }

  set.seed(SEED)
  cg1 <- colnames(merged_obj)[merged_obj$group == "PDX"]
  cg2 <- colnames(merged_obj)[merged_obj$group == "Control"]
  n1_used <- length(cg1); n2_used <- length(cg2)
  if (n1_used > MAX_CELLS) { cg1 <- sample(cg1, MAX_CELLS); n1_used <- MAX_CELLS }
  if (n2_used > MAX_CELLS) { cg2 <- sample(cg2, MAX_CELLS); n2_used <- MAX_CELLS }

  obj_sub <- subset(merged_obj, cells = c(cg1, cg2))
  Idents(obj_sub) <- obj_sub$group

  plan("multicore", workers = fm_cores)
  res <- tryCatch(
    FindMarkers(obj_sub,
                ident.1 = "PDX", ident.2 = "Control",
                test.use = "wilcox",
                min.pct = MIN_PCT,
                logfc.threshold = 0,
                only.pos = FALSE,
                verbose = FALSE),
    error = function(e) stop(sprintf("FindMarkers [%s]: %s", key, e$message))
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
    ifelse(res$avg_log2FC > 0, "Up", "Down"), "NS"
  )

  list(status = "ok", n1 = n1, n2 = n2,
       n1_used = n1_used, n2_used = n2_used, res = res)
}

# ============================================================
# 6. Worker: 单任务执行
# ============================================================
run_one_task <- function(task) {
  suppressPackageStartupMessages({
    library(igraph); library(Seurat); library(presto)
    library(future); library(dplyr)
  })
  select    <- dplyr::select
  rename    <- dplyr::rename
  filter    <- dplyr::filter
  arrange   <- dplyr::arrange
  bind_rows <- dplyr::bind_rows

  result <- list(
    key       = task$key,
    mouse_id  = task$mouse_id,
    pdx_sample= task$pdx_sample,
    virtual_ct= task$virtual_ct,
    organ     = task$organ,
    status    = "error",
    n1 = NA, n2 = NA, n1_used = NA, n2_used = NA,
    n_up = NA, n_down = NA,
    il6st_fc = NA, smad3_fc = NA, elmo1_fc = NA, h2d1_fc = NA
  )

  tryCatch({
    # 在 worker 内做 subset + merge (主进程只保留原始 seurat_v4)
    merged <- build_single_pair_merged(
      task$pdx_sample, task$ctrl_sample, task$virtual_ct,
      if (task$organ == "liver") LIVER_CELL_TYPES else BRAIN_CELL_TYPES
    )
    if (is.null(merged)) {
      result$status <- "insufficient_cells"
      return(result)
    }

    out <- run_one_deg(merged, task$key, fm_cores = CORES_PER_TASK)
    rm(merged); gc(verbose = FALSE)

    result$status  <- out$status
    result$n1      <- out$n1
    result$n2      <- out$n2
    if (!is.null(out$n1_used)) result$n1_used <- out$n1_used
    if (!is.null(out$n2_used)) result$n2_used <- out$n2_used

    if (out$status == "ok" && !is.null(out$res)) {
      write.csv(out$res,
                file.path(task$deg_dir, paste0(task$key, ".csv")),
                row.names = FALSE)

      result$n_up   <- sum(out$res$sig == "Up",   na.rm = TRUE)
      result$n_down <- sum(out$res$sig == "Down", na.rm = TRUE)

      # 关键基因追踪 (配对分析的 4 个 marker gene)
      get_fc <- function(res, gene) {
        row <- res[res$gene == gene, , drop = FALSE]
        if (nrow(row) > 0) row$avg_log2FC[1] else NA
      }
      result$il6st_fc <- get_fc(out$res, "Il6st")
      result$smad3_fc <- get_fc(out$res, "Smad3")
      result$elmo1_fc <- get_fc(out$res, "Elmo1")
      result$h2d1_fc  <- get_fc(out$res, "H2-D1")

      rm(out); gc(verbose = FALSE)
    }

    message(sprintf("  ✓ [%s|%s|%s] n1=%d n2=%d Up=%s Down=%s",
                    task$mouse_id, task$organ, task$virtual_ct,
                    result$n1, result$n2, result$n_up, result$n_down))
    return(result)

  }, error = function(e) {
    msg <- sprintf("[%s] task[%s] FAILED: %s",
                   format(Sys.time(), "%H:%M:%S"), task$key, e$message)
    message(msg)
    result$status <- "error"
    return(result)
  })
}

# ============================================================
# 7. 子分析 1: 构建 per-sample DEG 任务列表
# ============================================================
cat("========== Sub-analysis 1: Building per-sample DEG tasks ==========\n")

all_tasks <- list()
for (mouse in PAIRED_MICE) {
  # 肝: 8 个 virtual_ct
  for (vct in names(LIVER_CELL_TYPES)) {
    if (is.null(LIVER_CELL_TYPES[[vct]][["v4_refined"]])) next
    all_tasks[[length(all_tasks) + 1]] <- list(
      mouse_id   = mouse$id,
      pdx_sample = mouse$liver,
      ctrl_sample= CONTROL_LIVER,
      virtual_ct = vct,
      organ      = "liver",
      key        = paste0(mouse$id, "_liver_", vct),
      deg_dir    = file.path(OUT_DIR, "per_sample_deg", "liver")
    )
  }
  # 脑: 5 个 virtual_ct
  for (vct in names(BRAIN_CELL_TYPES)) {
    if (is.null(BRAIN_CELL_TYPES[[vct]][["v4_refined"]])) next
    all_tasks[[length(all_tasks) + 1]] <- list(
      mouse_id   = mouse$id,
      pdx_sample = mouse$brain,
      ctrl_sample= CONTROL_BRAIN,
      virtual_ct = vct,
      organ      = "brain",
      key        = paste0(mouse$id, "_brain_", vct),
      deg_dir    = file.path(OUT_DIR, "per_sample_deg", "brain")
    )
  }
}
cat(sprintf("  total tasks: %d (3 mice x (%d liver + %d brain) = %d)\n",
            length(all_tasks),
            length(LIVER_CELL_TYPES), length(BRAIN_CELL_TYPES),
            3 * (length(LIVER_CELL_TYPES) + length(BRAIN_CELL_TYPES))))

# ============================================================
# 8. 执行 per-sample DEG
# ============================================================
cat("\n========== Running per-sample DEG (parallel) ==========\n")
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
# 9. 收集任务状态
# ============================================================
cat("========== Collecting task status ==========\n")

status_rows <- list()
n_ok <- 0; n_insufficient <- 0; n_no_deg <- 0; n_error <- 0

for (r in results_raw) {
  if (is.null(r) || !is.list(r)) { n_error <- n_error + 1; next }

  status_rows[[length(status_rows) + 1]] <- data.frame(
    key        = r$key,
    mouse_id   = r$mouse_id,
    pdx_sample = r$pdx_sample,
    virtual_ct = r$virtual_ct,
    organ      = r$organ,
    status     = r$status,
    n1         = if (is.null(r$n1)) NA else r$n1,
    n2         = if (is.null(r$n2)) NA else r$n2,
    n1_used    = if (is.null(r$n1_used)) NA else r$n1_used,
    n2_used    = if (is.null(r$n2_used)) NA else r$n2_used,
    n_up       = if (is.null(r$n_up)) NA else r$n_up,
    n_down     = if (is.null(r$n_down)) NA else r$n_down,
    il6st_fc   = if (is.null(r$il6st_fc)) NA else r$il6st_fc,
    smad3_fc   = if (is.null(r$smad3_fc)) NA else r$smad3_fc,
    elmo1_fc   = if (is.null(r$elmo1_fc)) NA else r$elmo1_fc,
    h2d1_fc    = if (is.null(r$h2d1_fc)) NA else r$h2d1_fc,
    stringsAsFactors = FALSE
  )

  if (r$status == "ok") n_ok <- n_ok + 1
  else if (r$status == "insufficient_cells") n_insufficient <- n_insufficient + 1
  else if (r$status == "no_deg") n_no_deg <- n_no_deg + 1
  else n_error <- n_error + 1
}

task_status_df <- dplyr::bind_rows(status_rows)
write.csv(task_status_df, TASK_STATUS, row.names = FALSE)

cat(sprintf("  ok                 : %d\n", n_ok))
cat(sprintf("  insufficient_cells : %d\n", n_insufficient))
cat(sprintf("  no_deg             : %d\n", n_no_deg))
cat(sprintf("  error              : %d\n\n", n_error))

# ============================================================
# 10. Helper: 读回 DEG CSV
# ============================================================
read_deg_csv <- function(mouse_id, organ, virtual_ct) {
  fp <- file.path(OUT_DIR, "per_sample_deg", organ,
                  paste0(mouse_id, "_", organ, "_", virtual_ct, ".csv"))
  if (!file.exists(fp)) return(NULL)
  tryCatch(read.csv(fp, stringsAsFactors = FALSE), error = function(e) NULL)
}

# ============================================================
# 11. 子分析 2: Liver ↔ Brain DEG 重叠 (每只鼠)
# ============================================================
cat("========== Sub-analysis 2: Liver <-> Brain DEG overlap ==========\n")

# 对每只鼠, 合并所有细类的 Up 和 Down 基因作为"整体响应 gene set"
collect_all_sig_genes <- function(mouse_id, organ) {
  cell_type_map <- if (organ == "liver") LIVER_CELL_TYPES else BRAIN_CELL_TYPES
  up_set <- character(0); down_set <- character(0)
  universe <- character(0)
  for (vct in names(cell_type_map)) {
    if (is.null(cell_type_map[[vct]][["v4_refined"]])) next
    df <- read_deg_csv(mouse_id, organ, vct)
    if (is.null(df) || nrow(df) == 0) next
    universe <- union(universe, df$gene)
    up_set   <- union(up_set,   df$gene[df$sig == "Up"])
    down_set <- union(down_set, df$gene[df$sig == "Down"])
  }
  list(up = up_set, down = down_set, universe = universe)
}

overlap_rows <- list()
for (mouse in PAIRED_MICE) {
  cat(sprintf("  [%s] %s(liver) <-> %s(brain)\n",
              mouse$id, mouse$liver, mouse$brain))

  liver_sets <- collect_all_sig_genes(mouse$id, "liver")
  brain_sets <- collect_all_sig_genes(mouse$id, "brain")

  # universe = 肝脑两个 DEG 表里出现过的基因并集
  universe <- union(liver_sets$universe, brain_sets$universe)
  N <- length(universe)

  # 计算 4 种重叠
  overlaps <- list(
    list(name = "Up_Up",     liver = liver_sets$up,   brain = brain_sets$up),
    list(name = "Down_Down", liver = liver_sets$down, brain = brain_sets$down),
    list(name = "Up_Down",   liver = liver_sets$up,   brain = brain_sets$down),
    list(name = "Down_Up",   liver = liver_sets$down, brain = brain_sets$up)
  )

  for (ov in overlaps) {
    k <- length(intersect(ov$liver, ov$brain))  # overlap size
    K <- length(ov$liver)                        # liver set size
    n <- length(ov$brain)                        # brain set size

    # Hypergeometric test: 在 N 个基因里, 抽 n 个, 抽到至少 k 个属于 liver set 的概率
    if (K > 0 && n > 0 && N > 0) {
      pval <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)
      expected <- K * n / N
    } else {
      pval <- NA; expected <- NA
    }

    overlap_rows[[length(overlap_rows) + 1]] <- data.frame(
      mouse_id      = mouse$id,
      overlap_type  = ov$name,
      liver_n       = K,
      brain_n       = n,
      overlap_n     = k,
      expected      = round(expected, 2),
      fold_enrich   = if (is.na(expected) || expected == 0) NA else round(k / expected, 2),
      p_hyper       = pval,
      universe_N    = N,
      shared_genes  = paste(head(intersect(ov$liver, ov$brain), 20), collapse = ";"),
      stringsAsFactors = FALSE
    )
    cat(sprintf("    %-10s : %d / %d (liver) vs %d (brain), expected=%.1f, p=%.2e\n",
                ov$name, k, K, n,
                if (is.na(expected)) 0 else expected,
                if (is.na(pval)) 1 else pval))
  }

  # 每只鼠单独保存一份文件
  write.csv(
    dplyr::bind_rows(overlap_rows[(length(overlap_rows) - 3):length(overlap_rows)]),
    file.path(OUT_DIR, "overlap", paste0(mouse$id, "_overlap_summary.csv")),
    row.names = FALSE
  )
}

overlap_all <- dplyr::bind_rows(overlap_rows)
write.csv(overlap_all,
          file.path(OUT_DIR, "overlap", "all_mice_overlap.csv"),
          row.names = FALSE)
cat("  done\n\n")

# ============================================================
# 12. 子分析 3: 机制基因的跨器官 FC 散点图
# ============================================================
cat("========== Sub-analysis 3: Mechanism gene FC correlation ==========\n")

# 对每只鼠, 对每个机制基因取"该鼠所有肝细类中的最大 |FC|" 和 "该鼠所有脑细类中的最大 |FC|"
# (因为同一个基因在不同细类下 FC 不同, 取代表值)
get_max_fc_per_gene <- function(mouse_id, organ, gene_list) {
  cell_type_map <- if (organ == "liver") LIVER_CELL_TYPES else BRAIN_CELL_TYPES
  out <- data.frame(gene = gene_list,
                    max_abs_fc = 0, best_fc = 0, best_ct = NA_character_,
                    stringsAsFactors = FALSE)

  for (vct in names(cell_type_map)) {
    if (is.null(cell_type_map[[vct]][["v4_refined"]])) next
    df <- read_deg_csv(mouse_id, organ, vct)
    if (is.null(df) || nrow(df) == 0) next
    sub <- df[df$gene %in% gene_list, , drop = FALSE]
    if (nrow(sub) == 0) next
    for (i in seq_len(nrow(sub))) {
      g <- sub$gene[i]; fc <- sub$avg_log2FC[i]
      idx <- which(out$gene == g)
      if (abs(fc) > out$max_abs_fc[idx]) {
        out$max_abs_fc[idx] <- abs(fc)
        out$best_fc[idx]    <- fc
        out$best_ct[idx]    <- vct
      }
    }
  }
  out
}

for (mouse in PAIRED_MICE) {
  tryCatch({
    liver_fc <- get_max_fc_per_gene(mouse$id, "liver", ALL_MECH_GENES)
    brain_fc <- get_max_fc_per_gene(mouse$id, "brain", ALL_MECH_GENES)

    merged <- data.frame(
      gene       = ALL_MECH_GENES,
      liver_fc   = liver_fc$best_fc,
      liver_ct   = liver_fc$best_ct,
      brain_fc   = brain_fc$best_fc,
      brain_ct   = brain_fc$best_ct,
      stringsAsFactors = FALSE
    )

    # 通路标注
    pw_map <- dplyr::bind_rows(lapply(names(MECHANISM_GENES), function(nm) {
      data.frame(gene = MECHANISM_GENES[[nm]], pathway = nm,
                 stringsAsFactors = FALSE)
    }))
    pw_map <- pw_map[!duplicated(pw_map$gene), ]
    merged$pathway <- pw_map$pathway[match(merged$gene, pw_map$gene)]

    # 只画有信号的基因 (任意一侧 |FC| > 0.25)
    plot_df <- merged[abs(merged$liver_fc) > 0.25 | abs(merged$brain_fc) > 0.25, ]

    write.csv(merged,
              file.path(OUT_DIR, "fc_correlation",
                        paste0(mouse$id, "_mechanism_fc.csv")),
              row.names = FALSE)

    if (nrow(plot_df) < 3) {
      cat(sprintf("  [%s] too few genes with |FC| > 0.25, skip plot\n", mouse$id))
      next
    }

    # Pearson 相关
    cor_val <- tryCatch(
      cor(plot_df$liver_fc, plot_df$brain_fc, method = "pearson", use = "complete.obs"),
      error = function(e) NA
    )

    # 散点图
    label_df <- plot_df[abs(plot_df$liver_fc) > 0.5 | abs(plot_df$brain_fc) > 0.5, ]

    p <- ggplot(plot_df, aes(x = liver_fc, y = brain_fc, color = pathway)) +
      geom_hline(yintercept = 0, linewidth = 0.3, color = "grey60") +
      geom_vline(xintercept = 0, linewidth = 0.3, color = "grey60") +
      geom_abline(slope = 1, intercept = 0, linewidth = 0.3,
                  linetype = "dashed", color = "grey40") +
      geom_point(size = 2.5, alpha = 0.8) +
      geom_text_repel(data = label_df, aes(label = gene),
                      size = 3, max.overlaps = 30, segment.size = 0.2,
                      show.legend = FALSE) +
      scale_color_brewer(palette = "Set1", na.value = "grey50") +
      labs(title = sprintf("Mechanism gene FC: %s (%s liver <-> %s brain)",
                           mouse$id, mouse$liver, mouse$brain),
           subtitle = sprintf("Pearson r = %.3f, n = %d genes (|FC| > 0.25 in at least one organ)",
                              if (is.na(cor_val)) 0 else cor_val, nrow(plot_df)),
           x = "Liver max |FC| (with sign)",
           y = "Brain max |FC| (with sign)",
           color = "Pathway") +
      theme_bw(base_size = 11) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
            legend.position = "right",
            legend.text = element_text(size = 7))

    base <- file.path(OUT_DIR, "fc_correlation",
                      paste0(mouse$id, "_mechanism_scatter"))
    ggsave(paste0(base, ".pdf"), p, width = 10, height = 7, dpi = 300)
    ggsave(paste0(base, ".png"), p, width = 10, height = 7, dpi = 200)

    cat(sprintf("  [%s] Pearson r = %.3f (n=%d genes plotted)\n",
                mouse$id, if (is.na(cor_val)) 0 else cor_val, nrow(plot_df)))
  }, error = function(e) {
    write_error(paste0("fc_cor_", mouse$id), e$message)
  })
}
cat("\n")

# ============================================================
# 13. 子分析 4: 响应强度 ranking
# ============================================================
cat("========== Sub-analysis 4: Response strength ranking ==========\n")

# 每只鼠 × 每个细类的响应强度分数
# score_a: sum(|FC|) over significant genes
# score_b: n_up + n_down
# score_c: mean(|FC|) for key 10 mechanism genes

KEY_10_GENES <- c("Il6st", "Smad3", "Stat3", "Jak1", "Jak2",
                  "Socs3", "Tgfb1", "Il6", "Nfkb1", "Rela")

response_rows <- list()
for (mouse in PAIRED_MICE) {
  for (organ in c("liver", "brain")) {
    cell_type_map <- if (organ == "liver") LIVER_CELL_TYPES else BRAIN_CELL_TYPES
    organ_sum_fc <- 0
    organ_n_sig  <- 0
    key_fc_vals  <- numeric(0)

    for (vct in names(cell_type_map)) {
      if (is.null(cell_type_map[[vct]][["v4_refined"]])) next
      df <- read_deg_csv(mouse$id, organ, vct)
      if (is.null(df) || nrow(df) == 0) next

      sig_df <- df[df$sig != "NS", , drop = FALSE]
      score_a <- sum(abs(sig_df$avg_log2FC), na.rm = TRUE)
      score_b <- nrow(sig_df)
      organ_sum_fc <- organ_sum_fc + score_a
      organ_n_sig  <- organ_n_sig + score_b

      key_fc_vals <- c(key_fc_vals,
                       abs(df$avg_log2FC[df$gene %in% KEY_10_GENES]))
    }

    response_rows[[length(response_rows) + 1]] <- data.frame(
      mouse_id     = mouse$id,
      organ        = organ,
      sample       = if (organ == "liver") mouse$liver else mouse$brain,
      score_a_sumfc= round(organ_sum_fc, 2),
      score_b_nsig = organ_n_sig,
      score_c_keyfc= round(mean(key_fc_vals, na.rm = TRUE), 3),
      stringsAsFactors = FALSE
    )
  }
}
response_df <- dplyr::bind_rows(response_rows)
write.csv(response_df,
          file.path(OUT_DIR, "response_ranking", "per_sample_response_scores.csv"),
          row.names = FALSE)

cat("  Response scores:\n")
print(response_df, row.names = FALSE)

# Ranking 图
tryCatch({
  # 整理为 wide
  wide_a <- response_df %>%
    dplyr::select(mouse_id, organ, score_a_sumfc) %>%
    tidyr::pivot_wider(names_from = organ, values_from = score_a_sumfc)

  cor_val_a <- tryCatch(
    cor(wide_a$liver, wide_a$brain, method = "spearman"),
    error = function(e) NA
  )

  p <- ggplot(wide_a, aes(x = liver, y = brain, label = mouse_id)) +
    geom_point(size = 4, color = "#E41A1C") +
    geom_text_repel(size = 4, fontface = "bold") +
    geom_smooth(method = "lm", se = FALSE, color = "grey50",
                linetype = "dashed", linewidth = 0.5) +
    labs(title = "Paired response strength: Liver vs Brain",
         subtitle = sprintf("Spearman rho = %.3f (n=3 mice, descriptive only)",
                            if (is.na(cor_val_a)) 0 else cor_val_a),
         x = "Liver response score (sum |FC| across cell types)",
         y = "Brain response score (sum |FC| across cell types)") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, color = "grey40"))

  base <- file.path(OUT_DIR, "response_ranking", "ranking_plot")
  ggsave(paste0(base, ".pdf"), p, width = 8, height = 7, dpi = 300)
  ggsave(paste0(base, ".png"), p, width = 8, height = 7, dpi = 200)
  cat(sprintf("\n  Spearman rho (liver response vs brain response) = %.3f\n",
              if (is.na(cor_val_a)) 0 else cor_val_a))
}, error = function(e) write_error("ranking_plot", e$message))

cat("\n")

# ============================================================
# 14. 生成 summary markdown
# ============================================================
cat("========== Writing paired_summary.md ==========\n")

summary_lines <- c(
  "# Stage 3: Cross-organ paired analysis summary",
  "",
  sprintf("**生成时间：** %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## 任务统计",
  sprintf("- 总任务: %d", length(all_tasks)),
  sprintf("- 成功 (ok): %d", n_ok),
  sprintf("- insufficient_cells: %d", n_insufficient),
  sprintf("- no_deg: %d", n_no_deg),
  sprintf("- error: %d", n_error),
  sprintf("- 运行时间: %.1f min", t_elapsed["elapsed"] / 60),
  "",
  "## 配对鼠",
  "",
  "| mouse_id | liver sample | brain sample |",
  "|---|---|---|"
)
for (mouse in PAIRED_MICE) {
  summary_lines <- c(summary_lines,
    sprintf("| %s | %s | %s |", mouse$id, mouse$liver, mouse$brain))
}

summary_lines <- c(summary_lines, "",
  "## 子分析 2: DEG overlap (同一只鼠的肝 ↔ 脑基因重叠)",
  "",
  "以下表按 mouse_id × overlap_type 排列。每行是一只鼠的肝 DEG 基因集和脑 DEG 基因集在某个方向组合下的重叠统计。",
  "",
  "| Mouse | Type | Liver N | Brain N | Overlap | Expected | Fold | p_hyper |",
  "|---|---|---:|---:|---:|---:|---:|---:|"
)
for (i in seq_len(nrow(overlap_all))) {
  r <- overlap_all[i, ]
  summary_lines <- c(summary_lines,
    sprintf("| %s | %s | %d | %d | %d | %.1f | %.2f | %.2e |",
            r$mouse_id, r$overlap_type, r$liver_n, r$brain_n, r$overlap_n,
            if (is.na(r$expected)) 0 else r$expected,
            if (is.na(r$fold_enrich)) 0 else r$fold_enrich,
            if (is.na(r$p_hyper)) 1 else r$p_hyper))
}

summary_lines <- c(summary_lines, "",
  "## 子分析 4: 响应强度 ranking",
  "",
  "| Mouse | Organ | Sample | Sum\\|FC\\| | n_sig_genes | Key10 mean\\|FC\\| |",
  "|---|---|---|---:|---:|---:|"
)
for (i in seq_len(nrow(response_df))) {
  r <- response_df[i, ]
  summary_lines <- c(summary_lines,
    sprintf("| %s | %s | %s | %.2f | %d | %.3f |",
            r$mouse_id, r$organ, r$sample,
            r$score_a_sumfc, r$score_b_nsig, r$score_c_keyfc))
}

summary_lines <- c(summary_lines, "",
  "## 核心基因在 3 只鼠的跨器官响应",
  "",
  "每只鼠每个细类的 Il6st / Smad3 / Elmo1 / H2-D1 FC (NA = 该基因没出现在 DEG 表/min.pct 过滤)",
  "",
  "见 `summary/task_status.csv` 的 `il6st_fc / smad3_fc / elmo1_fc / h2d1_fc` 列",
  "",
  "## 文件位置",
  "",
  "- Per-sample DEG: `per_sample_deg/{liver,brain}/<mouse>_<organ>_<virtual_ct>.csv`",
  "- Overlap: `overlap/<mouse>_overlap_summary.csv`, `overlap/all_mice_overlap.csv`",
  "- FC correlation: `fc_correlation/<mouse>_mechanism_scatter.{pdf,png}`, `<mouse>_mechanism_fc.csv`",
  "- Response ranking: `response_ranking/per_sample_response_scores.csv`, `ranking_plot.{pdf,png}`",
  "- Task status: `summary/task_status.csv`",
  "- Error log: `summary/error_log.txt`"
)

writeLines(summary_lines, PAIRED_MD)

cat("\n========================================\n")
cat("Stage 3 complete\n")
cat(sprintf("  Total tasks  : %d\n", length(all_tasks)))
cat(sprintf("  Success (ok) : %d\n", n_ok))
cat(sprintf("  Elapsed      : %.1f min\n", t_elapsed["elapsed"] / 60))
cat(sprintf("  Output dir   : %s\n", OUT_DIR))
cat(sprintf("  Summary      : %s\n", PAIRED_MD))
cat("========================================\n")
