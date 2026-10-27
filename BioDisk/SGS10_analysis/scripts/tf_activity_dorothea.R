#!/usr/bin/env Rscript
# ============================================================
# tf_activity_dorothea.R --- DoRothEA + decoupleR TF 活性分析 (模式 A)
# ============================================================
# 用法:
#   cd ~/SGS10_analysis
#   Rscript scripts/tf_activity_dorothea.R 2>&1 | tee scripts/logs/tf_activity.log
#
# 前置:
#   - Stage 1 已完成, diff_stage1_results/deg_v4_refined/ 下有 DEG CSV
#   - 安装 decoupleR (见末尾的依赖说明)
#
# 做什么:
#   1. 加载 DoRothEA mouse 数据库 (A+B+C 置信度)
#   2. 对每个 Stage 1 v4_refined DEG CSV 读取 avg_log2FC 作为 gene signature
#   3. 用 decoupleR::run_ulm 计算 TF 活性 (score + p value)
#   4. 汇总成跨对比跨细类的 TF 活性矩阵
#   5. 生成核心图: TF × 对比热图 + 三轴核心 TF 单独柱图
#   6. 输出 summary md
#
# 输出:
#   dorothea_results/
#     per_contrast/*_tf_activity.csv          (39 个文件, 每个对比一个)
#     summary/
#       tf_activity_matrix.csv               (TF × 对比 的大宽表)
#       top_active_tfs_per_contrast.csv      (每对比 top 20 激活 TF)
#       core_axis_tf_mapping.csv             (三轴核心 TF 汇总)
#       task_status.csv
#       dorothea_summary.md
#     plots/
#       tf_activity_heatmap_brain.pdf
#       tf_activity_heatmap_liver.pdf
#       core_tf_by_contrast.pdf
#
# 资源:
#   内存 ~2 GB, CPU 单核, 时间 15-30 min
#   可以和拟时序并行跑
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
})

# 尝试加载 decoupleR (核心包)
if (!requireNamespace("decoupleR", quietly = TRUE)) {
  stop("ERROR: decoupleR 未安装。请先运行:\n",
       "  BiocManager::install('decoupleR')\n",
       "或:\n",
       "  remotes::install_github('saezlab/decoupleR')")
}
library(decoupleR)

# dplyr namespace lock
select    <- dplyr::select
rename    <- dplyr::rename
filter    <- dplyr::filter
mutate    <- dplyr::mutate
arrange   <- dplyr::arrange
bind_rows <- dplyr::bind_rows

# ============================================================
# 0. 路径 & 参数
# ============================================================
PROJECT_ROOT <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
setwd(PROJECT_ROOT)

DEG_INPUT_DIR <- file.path(PROJECT_ROOT, "diff_stage1_results", "deg_v4_refined")
OUT_DIR       <- file.path(PROJECT_ROOT, "dorothea_results")

for (d in c("per_contrast", "summary", "plots")) {
  dir.create(file.path(OUT_DIR, d), recursive = TRUE, showWarnings = FALSE)
}

ERROR_LOG <- file.path(OUT_DIR, "summary", "error_log.txt")
writeLines(sprintf("== DoRothEA error log - %s ==",
                   format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ERROR_LOG)

write_error <- function(context, msg) {
  line <- sprintf("[%s] %s => %s",
                  format(Sys.time(), "%H:%M:%S"), context, msg)
  cat("  ⚠", line, "\n")
  cat(line, "\n", file = ERROR_LOG, append = TRUE)
}

# 置信度过滤
CONFIDENCE_LEVELS <- c("A", "B", "C")  # 不用 D, E (噪声多)

# 三轴核心 TF (预先感兴趣的 TF, 报告里专门列出)
# 注意: DoRothEA mouse 数据库用鼠源 Symbol 格式 (首字母大写), 不是人源全大写
CORE_AXIS_TFS <- list(
  # 轴 1: BBB + 脑内皮 MHC-I
  Axis1_BBB = c("Klf4", "Klf2", "Sox17", "Erg", "Foxf2", "Mef2c"),
  Axis1_MHC = c("Stat1", "Irf1", "Irf2", "Irf3", "Nlrc5", "Ciita"),
  # 轴 2: Efferocytosis / Microglia
  Axis2_Microglia = c("Spi1", "Irf8", "Runx1", "Cebpb", "Batf3", "Mafb"),
  # 轴 3: TGFb/IL-6
  Axis3_TGFb = c("Smad2", "Smad3", "Smad4", "Smad7", "Foxp3"),
  Axis3_IL6  = c("Stat3", "Stat1", "Nfkb1", "Rela", "Rel", "Nfkb2"),
  # 常见炎症
  Inflammation = c("Nfkb1", "Rela", "Nfe2l2", "Atf3", "Fos", "Jun", "Junb",
                    "Jund", "Egr1", "Egr2", "Atf4")
)
# 合并成 flat list
CORE_TFS_FLAT <- unique(unlist(CORE_AXIS_TFS))

# ULM 的 minimum target 参数
MIN_N_TARGETS <- 5

cat("========================================\n")
cat("DoRothEA + decoupleR TF activity analysis\n")
cat("Start:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================\n\n")

# ============================================================
# 1. 加载 DoRothEA mouse 数据库
# ============================================================
cat("========== Loading DoRothEA mouse database ==========\n")

# 使用 dorothea 包的内嵌数据 dorothea_mm (不需要网络)
# decoupleR::get_dorothea 会尝试从 UniProt/OmnipathR 下载, 在无网服务器上失败
if (!requireNamespace("dorothea", quietly = TRUE)) {
  stop("dorothea 包未安装。请先运行:\n",
       "  BiocManager::install('dorothea')")
}

data("dorothea_mm", package = "dorothea")

# 过滤置信度 A+B+C
dorothea_db <- dorothea_mm %>%
  dplyr::filter(confidence %in% CONFIDENCE_LEVELS)

# decoupleR::run_ulm 要求列名: source, target, mor
# dorothea 包原始列名: tf, target, confidence, mor
dorothea_db <- dorothea_db %>%
  dplyr::select(source = tf, target, mor, confidence)

cat(sprintf("  database rows: %d\n", nrow(dorothea_db)))
cat(sprintf("  unique TFs: %d\n", length(unique(dorothea_db$source))))
cat(sprintf("  unique targets: %d\n", length(unique(dorothea_db$target))))
cat(sprintf("  confidence levels kept: %s\n",
            paste(CONFIDENCE_LEVELS, collapse = ",")))
cat("  head of database:\n")
print(head(dorothea_db, 5))

# 清理原始对象 (dorothea_mm 有 380K 行, 保留 filter 后的即可)
rm(dorothea_mm)
gc(verbose = FALSE)
cat("\n")

# ============================================================
# 2. 发现所有 DEG CSV 文件
# ============================================================
cat("========== Discovering DEG CSV files ==========\n")

deg_files <- list.files(
  DEG_INPUT_DIR,
  pattern  = "\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

# 只保留 tumor_vs_control 和 tumor_vs_blank 的文件 (跳过 blank_vs_control)
# 因为 blank_vs_control 是 "手术对照效应" 不是主要兴趣
deg_files_keep <- deg_files[
  grepl("tumor_vs_control|tumor_vs_blank", basename(deg_files))
]
cat(sprintf("  total DEG files: %d\n", length(deg_files)))
cat(sprintf("  tumor contrasts: %d\n", length(deg_files_keep)))

if (length(deg_files_keep) == 0) {
  stop("没找到任何 tumor_vs_* DEG 文件")
}

# ============================================================
# 3. 解析文件名, 提取元信息
# ============================================================
parse_deg_filename <- function(fp) {
  bn <- tools::file_path_sans_ext(basename(fp))
  # 文件名格式: brain_tumor_vs_control_Microglia_Homeostatic
  #           : liver_tumor_vs_blank_Macrophage_NAM
  parts <- strsplit(bn, "_")[[1]]
  organ <- parts[1]  # brain or liver
  # 找 tumor_vs_control 或 tumor_vs_blank
  if (parts[3] == "vs" && parts[2] == "tumor") {
    contrast <- paste0(parts[2], "_vs_", parts[4])
    ct_start_idx <- 5
  } else {
    contrast <- "unknown"
    ct_start_idx <- 2
  }
  virtual_ct <- paste(parts[ct_start_idx:length(parts)], collapse = "_")
  list(
    organ      = organ,
    contrast   = contrast,
    virtual_ct = virtual_ct,
    key        = paste(organ, contrast, virtual_ct, sep = "__")
  )
}

# ============================================================
# 4. 对每个 DEG 文件跑 decoupleR::run_ulm
# ============================================================
cat("\n========== Running decoupleR::run_ulm on each contrast ==========\n")

run_one_contrast <- function(fp) {
  info <- parse_deg_filename(fp)
  key  <- info$key
  cat(sprintf("  [%s]\n", key))

  result <- list(
    key        = key,
    organ      = info$organ,
    contrast   = info$contrast,
    virtual_ct = info$virtual_ct,
    status     = "error",
    n_input_genes  = NA,
    n_tfs_output   = NA,
    n_sig_strict   = NA,
    n_sig_relaxed  = NA,
    elapsed_sec    = NA
  )
  t0 <- proc.time()

  tryCatch({
    deg <- read.csv(fp, stringsAsFactors = FALSE)
    if (!all(c("gene", "avg_log2FC") %in% colnames(deg))) {
      stop("DEG CSV 缺少必要列 (gene, avg_log2FC)")
    }

    # 构造 gene-level signature: 行 = 基因, 1 列 = avg_log2FC
    # decoupleR::run_ulm 要求 matrix 输入
    sig_mat <- matrix(
      deg$avg_log2FC,
      ncol = 1,
      dimnames = list(deg$gene, "log2fc")
    )
    # 移除 NA 和 Inf
    sig_mat <- sig_mat[is.finite(sig_mat[, 1]), , drop = FALSE]
    result$n_input_genes <- nrow(sig_mat)

    if (nrow(sig_mat) < 100) {
      result$status <- "too_few_genes"
      result$elapsed_sec <- round((proc.time() - t0)["elapsed"], 1)
      return(result)
    }

    # 运行 ULM
    ulm_res <- decoupleR::run_ulm(
      mat     = sig_mat,
      network = dorothea_db,
      .source = "source",
      .target = "target",
      .mor    = "mor",
      minsize = MIN_N_TARGETS
    )

    # ulm_res 是长格式: statistic, source (TF), condition, score, p_value
    # 我们只关心 ulm 统计量
    ulm_main <- ulm_res %>%
      dplyr::filter(statistic == "ulm") %>%
      dplyr::select(TF = source, score, p_value) %>%
      dplyr::mutate(
        padj = p.adjust(p_value, method = "BH"),
        # 两档显著性 (strict 给 paper, relaxed 给探索性分析)
        sig_strict = ifelse(padj < 0.05,
                            ifelse(score > 0, "Up", "Down"), "NS"),
        sig_relaxed = ifelse(padj < 0.10,
                             ifelse(score > 0, "Up", "Down"), "NS"),
        # 保留 sig 列作为旧代码兼容 (用 strict)
        sig = sig_strict
      ) %>%
      dplyr::arrange(padj, dplyr::desc(abs(score)))

    result$n_tfs_output <- nrow(ulm_main)
    result$n_sig_strict <- sum(ulm_main$sig_strict != "NS", na.rm = TRUE)
    result$n_sig_relaxed <- sum(ulm_main$sig_relaxed != "NS", na.rm = TRUE)

    # 保存到 per_contrast/
    out_fp <- file.path(OUT_DIR, "per_contrast",
                        paste0(key, "_tf_activity.csv"))
    write.csv(ulm_main, out_fp, row.names = FALSE)

    result$status <- "ok"
    result$elapsed_sec <- round((proc.time() - t0)["elapsed"], 1)
    cat(sprintf("    ✓ n_genes=%d, n_tfs=%d, strict=%d relaxed=%d (%.1fs)\n",
                result$n_input_genes, result$n_tfs_output,
                result$n_sig_strict, result$n_sig_relaxed, result$elapsed_sec))

  }, error = function(e) {
    write_error(key, e$message)
    result$elapsed_sec <- round((proc.time() - t0)["elapsed"], 1)
  })

  return(result)
}

# 并行还是顺序? 每个对比只用几秒, 单核串行就行
t_start <- proc.time()
all_results <- lapply(deg_files_keep, run_one_contrast)
t_elapsed <- (proc.time() - t_start)["elapsed"]
cat(sprintf("\n  all contrasts done in %.1f min\n\n", t_elapsed / 60))

# ============================================================
# 5. 汇总 task_status
# ============================================================
cat("========== Building task_status and summary ==========\n")

ts_df <- dplyr::bind_rows(lapply(all_results, function(r) {
  data.frame(
    key            = r$key,
    organ          = r$organ,
    contrast       = r$contrast,
    virtual_ct     = r$virtual_ct,
    status         = r$status,
    n_input_genes  = if (is.null(r$n_input_genes)) NA else r$n_input_genes,
    n_tfs_output   = if (is.null(r$n_tfs_output))  NA else r$n_tfs_output,
    n_sig_strict   = if (is.null(r$n_sig_strict))  NA else r$n_sig_strict,
    n_sig_relaxed  = if (is.null(r$n_sig_relaxed)) NA else r$n_sig_relaxed,
    elapsed_sec    = if (is.null(r$elapsed_sec))   NA else r$elapsed_sec,
    stringsAsFactors = FALSE
  )
}))
write.csv(ts_df, file.path(OUT_DIR, "summary", "task_status.csv"),
          row.names = FALSE)
cat("  task status:\n")
print(ts_df)

# ============================================================
# 6. 读回所有 per_contrast 结果, 做大宽表
# ============================================================
cat("\n========== Aggregating TF activity matrix ==========\n")

all_tf_rows <- list()
for (r in all_results) {
  if (r$status != "ok") next
  fp <- file.path(OUT_DIR, "per_contrast",
                  paste0(r$key, "_tf_activity.csv"))
  if (!file.exists(fp)) next
  df <- tryCatch(read.csv(fp, stringsAsFactors = FALSE),
                 error = function(e) NULL)
  if (is.null(df)) next
  df$key        <- r$key
  df$organ      <- r$organ
  df$contrast   <- r$contrast
  df$virtual_ct <- r$virtual_ct
  all_tf_rows[[length(all_tf_rows) + 1]] <- df
}

if (length(all_tf_rows) == 0) {
  stop("没有成功的 TF activity 结果可聚合")
}

tf_all <- dplyr::bind_rows(all_tf_rows)

# 长格式 → 宽格式 (TF 作为行, 对比作为列, cell = score)
tf_mat <- tf_all %>%
  dplyr::select(TF, key, score) %>%
  tidyr::pivot_wider(names_from = key, values_from = score, values_fill = NA) %>%
  as.data.frame()
rownames(tf_mat) <- tf_mat$TF
tf_mat_data <- as.matrix(tf_mat[, -1])

write.csv(tf_mat, file.path(OUT_DIR, "summary", "tf_activity_matrix.csv"),
          row.names = FALSE)
cat(sprintf("  TF matrix: %d TFs × %d contrasts\n",
            nrow(tf_mat_data), ncol(tf_mat_data)))

# ============================================================
# 7. Top TFs per contrast (用 relaxed 标准, 保留两档列)
# ============================================================
top_per_contrast <- tf_all %>%
  dplyr::filter(sig_relaxed != "NS") %>%
  dplyr::group_by(key, organ, contrast, virtual_ct) %>%
  dplyr::arrange(dplyr::desc(abs(score))) %>%
  dplyr::slice_head(n = 20) %>%
  dplyr::ungroup() %>%
  dplyr::select(organ, contrast, virtual_ct, TF, score, p_value, padj,
                sig_strict, sig_relaxed)

write.csv(top_per_contrast,
          file.path(OUT_DIR, "summary", "top_active_tfs_per_contrast.csv"),
          row.names = FALSE)
cat(sprintf("  top TFs per contrast (relaxed): %d rows\n", nrow(top_per_contrast)))

# ============================================================
# 8. 核心 TF mapping (三轴相关)
# ============================================================
cat("\n========== Core axis TF extraction ==========\n")

core_tf_rows <- tf_all %>%
  dplyr::filter(TF %in% CORE_TFS_FLAT) %>%
  dplyr::mutate(
    axis_category = sapply(TF, function(t) {
      cats <- c()
      for (cat_name in names(CORE_AXIS_TFS)) {
        if (t %in% CORE_AXIS_TFS[[cat_name]]) cats <- c(cats, cat_name)
      }
      if (length(cats) == 0) "Other" else paste(cats, collapse = ";")
    })
  ) %>%
  dplyr::select(organ, contrast, virtual_ct, TF, axis_category,
                score, p_value, padj, sig_strict, sig_relaxed, sig)

write.csv(core_tf_rows,
          file.path(OUT_DIR, "summary", "core_axis_tf_mapping.csv"),
          row.names = FALSE)
cat(sprintf("  core axis TF rows: %d\n", nrow(core_tf_rows)))

# ============================================================
# 8.5 Focused test: 只对三轴核心 TF 做 BH 校正
# ============================================================
# 动机: 默认 BH 校正 (run_ulm 之后) 是对所有 ~200 TF 一起做的, 太严格
# 对预先定义的 hypothesis set (30 个三轴核心 TF) 做 focused correction 更合理
# 同一 raw p value 下, focused correction 的 padj 会比 global padj 小很多
cat("\n========== Focused test on core axis TFs ==========\n")

core_focused <- tf_all %>%
  dplyr::filter(TF %in% CORE_TFS_FLAT) %>%
  dplyr::select(organ, contrast, virtual_ct, TF, score, p_value) %>%
  dplyr::mutate(
    axis_category = sapply(TF, function(t) {
      cats <- c()
      for (cat_name in names(CORE_AXIS_TFS)) {
        if (t %in% CORE_AXIS_TFS[[cat_name]]) cats <- c(cats, cat_name)
      }
      if (length(cats) == 0) "Other" else paste(cats, collapse = ";")
    })
  )

# 按每个对比独立做 focused BH (只对 core TF 校正)
core_focused <- core_focused %>%
  dplyr::group_by(organ, contrast, virtual_ct) %>%
  dplyr::mutate(
    padj_focused = p.adjust(p_value, method = "BH"),
    sig_focused_strict  = ifelse(padj_focused < 0.05,
                                 ifelse(score > 0, "Up", "Down"), "NS"),
    sig_focused_relaxed = ifelse(padj_focused < 0.10,
                                 ifelse(score > 0, "Up", "Down"), "NS")
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(organ, contrast, virtual_ct, padj_focused)

write.csv(core_focused,
          file.path(OUT_DIR, "summary", "core_tf_focused_test.csv"),
          row.names = FALSE)

n_focused_strict  <- sum(core_focused$sig_focused_strict  != "NS", na.rm = TRUE)
n_focused_relaxed <- sum(core_focused$sig_focused_relaxed != "NS", na.rm = TRUE)
cat(sprintf("  focused test rows: %d\n", nrow(core_focused)))
cat(sprintf("  focused strict sig: %d (padj_focused < 0.05)\n", n_focused_strict))
cat(sprintf("  focused relaxed sig: %d (padj_focused < 0.10)\n", n_focused_relaxed))

# 打印 focused test 最显著的 TF (strict 下)
if (n_focused_strict > 0) {
  cat("\n  Top focused-strict-sig core TFs:\n")
  top_focused <- core_focused %>%
    dplyr::filter(sig_focused_strict != "NS") %>%
    dplyr::arrange(padj_focused) %>%
    dplyr::slice_head(n = 20)
  print(as.data.frame(top_focused[, c("organ", "virtual_ct", "TF",
                                       "axis_category", "score", "padj_focused",
                                       "sig_focused_strict")]))
}

# ============================================================
# 9. 热图: TF × 对比 (分 organ)
# ============================================================
cat("\n========== Plotting TF activity heatmaps ==========\n")

plot_tf_heatmap <- function(organ_name) {
  organ_cols <- grep(paste0("^", organ_name, "__"), colnames(tf_mat_data), value = TRUE)
  if (length(organ_cols) == 0) return(invisible(NULL))

  mat_sub <- tf_mat_data[, organ_cols, drop = FALSE]

  # 过滤: 至少一列 |score| > 1.5 的 TF
  # (阈值从 2 降到 1.5, 配合 sig 定义从 padj+score 双门槛降为仅 padj)
  keep_rows <- apply(mat_sub, 1, function(r) {
    any(abs(r) > 1.5, na.rm = TRUE)
  })
  mat_sub <- mat_sub[keep_rows, , drop = FALSE]

  if (nrow(mat_sub) < 5) {
    cat(sprintf("  [%s] too few significant TFs (%d), skip heatmap\n",
                organ_name, nrow(mat_sub)))
    return(invisible(NULL))
  }

  # 裁剪色标 + NA 处理 (NA → 0, 因为 pheatmap clustering 无法处理 NA)
  mat_clip <- mat_sub
  mat_clip[is.na(mat_clip)] <- 0   # 关键: NA 替换为 0 避免 pheatmap 聚类失败
  mat_clip[mat_clip >  5] <-  5
  mat_clip[mat_clip < -5] <- -5

  # 简化列名 (去掉 organ__ 前缀)
  colnames(mat_clip) <- gsub(paste0("^", organ_name, "__"), "", colnames(mat_clip))

  cat(sprintf("  [%s] heatmap: %d TFs × %d contrasts\n",
              organ_name, nrow(mat_clip), ncol(mat_clip)))

  tryCatch({
    pheatmap(
      mat_clip,
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
      breaks = seq(-5, 5, length.out = 101),
      main = sprintf("%s TF activity (decoupleR ULM score, NA→0)\n%d TFs × %d contrasts",
                     toupper(organ_name), nrow(mat_clip), ncol(mat_clip)),
      fontsize_row = 6,
      fontsize_col = 7,
      filename = file.path(OUT_DIR, "plots",
                            paste0("tf_activity_heatmap_", organ_name, ".pdf")),
      width = max(8, 1 + ncol(mat_clip) * 0.6),
      height = max(6, 1 + nrow(mat_clip) * 0.16)
    )
  }, error = function(e) write_error(paste0("heatmap_", organ_name), e$message))
}

plot_tf_heatmap("liver")
plot_tf_heatmap("brain")

# ============================================================
# 10. 核心 TF 条形图
# ============================================================
cat("\n========== Plotting core axis TF bar plot ==========\n")

if (nrow(core_tf_rows) > 0) {
  # 构造三级显著性 level (for fill):
  #   Strict Up/Down: padj<0.05 激活/抑制
  #   Relaxed Up/Down: padj<0.10 激活/抑制 (strict 之外的)
  #   NS: padj>=0.10
  core_tf_rows$sig_level <- dplyr::case_when(
    core_tf_rows$sig_strict  == "Up"   ~ "Strict Up",
    core_tf_rows$sig_strict  == "Down" ~ "Strict Down",
    core_tf_rows$sig_relaxed == "Up"   ~ "Relaxed Up",
    core_tf_rows$sig_relaxed == "Down" ~ "Relaxed Down",
    TRUE                                ~ "NS"
  )
  core_tf_rows$sig_level <- factor(
    core_tf_rows$sig_level,
    levels = c("Strict Up", "Relaxed Up", "NS", "Relaxed Down", "Strict Down")
  )

  # 每个 core TF 在所有对比里的 score 柱状图
  tryCatch({
    p_core <- ggplot(
      core_tf_rows,
      aes(x = virtual_ct, y = score, fill = sig_level)
    ) +
      geom_bar(stat = "identity", position = "dodge") +
      facet_grid(TF ~ organ + contrast, scales = "free", space = "free") +
      scale_fill_manual(
        values = c("Strict Up"    = "#B2182B",
                   "Relaxed Up"   = "#F4A582",
                   "NS"           = "grey85",
                   "Relaxed Down" = "#92C5DE",
                   "Strict Down"  = "#2166AC")
      ) +
      geom_hline(yintercept = 0, linewidth = 0.3, color = "grey60") +
      geom_hline(yintercept = c(-1.5, 1.5), linewidth = 0.3,
                 color = "grey40", linetype = "dashed") +
      labs(title = "Core axis TF activity",
           subtitle = "Strict: padj<0.05 | Relaxed: padj<0.10 | dashed: |score| 1.5",
           x = NULL, y = "ULM score (TF activity)") +
      theme_bw(base_size = 9) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, size = 8),
            axis.text.x = element_text(angle = 60, hjust = 1, size = 6),
            strip.text.y = element_text(angle = 0, size = 6),
            legend.position = "bottom")

    n_tfs_core <- length(unique(core_tf_rows$TF))
    n_contrasts_core <- length(unique(core_tf_rows$virtual_ct))
    ggsave(
      file.path(OUT_DIR, "plots", "core_tf_by_contrast.pdf"),
      p_core,
      width  = max(10, 2 + n_contrasts_core * 0.5),
      height = max(8, 1 + n_tfs_core * 0.6),
      dpi = 300, limitsize = FALSE
    )
    cat(sprintf("  ✓ core TF bar plot saved\n"))
  }, error = function(e) write_error("core_tf_plot", e$message))
}

# ============================================================
# 11. 生成 summary md
# ============================================================
cat("\n========== Writing dorothea_summary.md ==========\n")

# 找几个特别重要的 TF 在关键细类的 score
extract_tf_score <- function(tf_name, contrast_pattern, ct_pattern) {
  row_idx <- which(
    grepl(tf_name, tf_all$TF, fixed = TRUE) &
    grepl(contrast_pattern, tf_all$contrast) &
    grepl(ct_pattern, tf_all$virtual_ct)
  )
  if (length(row_idx) == 0) return("N/A")
  sprintf("%.2f (padj=%.2g)",
          tf_all$score[row_idx[1]],
          tf_all$padj[row_idx[1]])
}

# 关键查询
key_tf_summary <- list(
  list(tf = "Smad3",  ct = "Microglia_Homeostatic", contrast = "tumor_vs_control"),
  list(tf = "Smad2",  ct = "Microglia_Homeostatic", contrast = "tumor_vs_control"),
  list(tf = "Stat3",  ct = "Microglia_Homeostatic", contrast = "tumor_vs_control"),
  list(tf = "Stat1",  ct = "Endothelial",           contrast = "tumor_vs_control"),
  list(tf = "Irf1",   ct = "Endothelial",           contrast = "tumor_vs_control"),
  list(tf = "Klf4",   ct = "Endothelial",           contrast = "tumor_vs_control"),
  list(tf = "Spi1",   ct = "Microglia_Homeostatic", contrast = "tumor_vs_control"),
  list(tf = "Irf8",   ct = "Microglia_Homeostatic", contrast = "tumor_vs_control"),
  list(tf = "Nfkb1",  ct = "Macrophage_NAM",        contrast = "tumor_vs_control"),
  list(tf = "Rela",   ct = "HSC_Activated_ECM",     contrast = "tumor_vs_control")
)

summary_lines <- c(
  "# DoRothEA + decoupleR TF 活性分析总结",
  "",
  sprintf("**生成时间：** %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  sprintf("**运行耗时：** %.1f min", t_elapsed / 60),
  "",
  "## 方法",
  "",
  "- **数据库：** DoRothEA mouse (CollecTRI)，置信度 A+B+C",
  "- **算法：** `decoupleR::run_ulm` (Univariate Linear Model on gene signatures)",
  "- **输入：** Stage 1 `v4_refined` DEG CSV 的 `avg_log2FC` 作为 gene-level signature",
  sprintf("- **对比范围：** %d 个 tumor_vs_* 对比 (跳过 blank_vs_control)",
          length(deg_files_keep)),
  sprintf("- **最小 target 数：** %d", MIN_N_TARGETS),
  "",
  "## 任务统计",
  "",
  sprintf("- 成功 (ok): %d", sum(ts_df$status == "ok")),
  sprintf("- 基因数不足: %d", sum(ts_df$status == "too_few_genes")),
  sprintf("- 错误: %d", sum(ts_df$status == "error")),
  sprintf("- 每对比平均耗时: %.1f 秒",
          mean(ts_df$elapsed_sec, na.rm = TRUE)),
  sprintf("- 平均显著 TF 数 (strict padj<0.05): %.1f",
          mean(ts_df$n_sig_strict, na.rm = TRUE)),
  sprintf("- 平均显著 TF 数 (relaxed padj<0.10): %.1f",
          mean(ts_df$n_sig_relaxed, na.rm = TRUE)),
  "",
  "## Focused test 统计 (只对三轴核心 TF 校正)",
  "",
  sprintf("- 核心 TF 数: %d", length(CORE_TFS_FLAT)),
  sprintf("- Focused test 总行数: %d", nrow(core_focused)),
  sprintf("- Focused strict sig (padj_focused<0.05): %d", n_focused_strict),
  sprintf("- Focused relaxed sig (padj_focused<0.10): %d", n_focused_relaxed),
  "",
  "## 关键 TF 活性速查",
  "",
  "| 轴 | TF | 细类 | 对比 | ULM score (padj) |",
  "|---|---|---|---|---|"
)

axis_labels <- list(
  Smad3 = "轴 3 TGFb", Smad2 = "轴 3 TGFb",
  Stat3 = "轴 3 IL-6", Stat1 = "轴 1 MHC-I", Irf1 = "轴 1 MHC-I",
  Klf4  = "轴 1 BBB",  Spi1 = "轴 2 Microglia", Irf8 = "轴 2 Microglia",
  Nfkb1 = "炎症通用", Rela = "炎症通用"
)

for (q in key_tf_summary) {
  score_str <- extract_tf_score(q$tf, q$contrast, q$ct)
  axis_str  <- axis_labels[[q$tf]]
  if (is.null(axis_str)) axis_str <- "-"
  summary_lines <- c(summary_lines,
    sprintf("| %s | %s | %s | %s | %s |",
            axis_str, q$tf, q$ct, q$contrast, score_str))
}

summary_lines <- c(summary_lines, "",
  "## 解读提示",
  "",
  "ULM score 的解读：",
  "- **score > 0**: TF target 基因集整体上调 → TF 活性增加",
  "- **score < 0**: TF target 基因集整体下调 → TF 活性降低",
  "- **padj < 0.05 (strict)**: 严格显著, paper 主图可用",
  "- **padj < 0.10 (relaxed)**: 宽松显著, 趋势级别, 附录/讨论可用",
  "- **padj_focused**: 只对 30 个三轴核心 TF 做 BH 校正, 比 global padj 宽松 5-6 倍",
  "",
  "Focused test 的使用场景: 你的 hypothesis 只关心三轴核心 TF (Smad3/Stat3/Klf4 等),",
  "对这个预先定义的 set 做 focused correction 比对全部 200 TF 做更合理.",
  "如果 padj_focused<0.05 但 global padj>0.05, 说明 TF 在核心 set 里显著但在全 TF 层面被稀释.",
  "",
  "这些 TF score 是**从 DEG 的 avg_log2FC 反推**的, 对应 Stage 1 的结果方向.",
  "如果某个 TF score 和其基因表达方向不一致, 这个 TF 可能是受翻译后修饰调控.",
  "",
  "## 文件位置",
  "",
  "- `per_contrast/` — 每个对比的 TF 活性表 (包含 score, padj, sig_strict, sig_relaxed)",
  "- `summary/tf_activity_matrix.csv` — TF × 对比 大宽表 (score only)",
  "- `summary/top_active_tfs_per_contrast.csv` — 每对比 relaxed 显著的 top 20 TF",
  "- `summary/core_axis_tf_mapping.csv` — 三轴核心 TF 汇总 (global padj)",
  "- `summary/core_tf_focused_test.csv` — **关键**: 只对核心 TF 做 BH 的 focused test 结果",
  "- `plots/tf_activity_heatmap_{liver,brain}.pdf` — 主热图",
  "- `plots/core_tf_by_contrast.pdf` — 核心 TF 条形图 (三级显著性着色)"
)

writeLines(summary_lines, file.path(OUT_DIR, "summary", "dorothea_summary.md"))

cat("\n========================================\n")
cat("DoRothEA analysis complete\n")
cat(sprintf("  Output dir : %s\n", OUT_DIR))
cat(sprintf("  Summary    : %s\n", file.path(OUT_DIR, "summary", "dorothea_summary.md")))
cat(sprintf("  Elapsed    : %.1f min\n", t_elapsed / 60))
cat("========================================\n")
