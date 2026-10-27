#!/usr/bin/env Rscript
# ============================================================
# nam_reannotate.R --- NAM marker-based 重注释
# ============================================================
# 用法:
#   cd ~/SGS10_analysis
#   Rscript scripts/nam_reannotate.R 2>&1 | tee scripts/logs/nam_reannotate.log
#
# 背景:
#   nam_subtype_check 揭示 v4_refined 的 "Macrophage_NAM" 标签严重污染了肝细胞:
#     - 298 个被标为 NAM 的细胞中, 173 个 (主要 blank+control) 实际是 hepatocyte
#     - 真 NAM 只有 ~125 个 (主要 tumor 组)
#     - 污染机制: 可能是 doublet / ambient RNA 没清干净
#
# 方案:
#   1. 在 liver 全细胞 (~136K) 上 marker-based 重新识别真 NAM
#   2. 用 "正向 marker AND 反向 marker" 双重筛选
#   3. 不修改原 RDS, 只输出 cell-level CSV (cell_id, v4_refined_clean)
#   4. 下游脚本只需要 left-join 这个 CSV 到 metadata
#
# 严格 NAM 定义:
#   PositiveMarkers (myeloid/macrophage 谱系): 至少 3 个 > 0
#     - Itgb2 (CD18, 白细胞 integrin)
#     - Cd74 (MHC-II invariant chain)
#     - Csf1r (M-CSF receptor)
#     - Adgre1 (F4/80)
#     - Lgmn (asparaginyl endopeptidase)
#     - Ctss (cathepsin S)
#     - Cd68 (pan-macrophage)
#     - Lyz2 (lysozyme)
#     - Aif1 (Iba1)
#   NegativeMarkers (hepatocyte): 必须全部 == 0
#     - Alb (白蛋白)
#     - Itih4
#     - Mat1a
#     - Pzp
#     - Fga
#     - Apoa1
#   ExcludeMarkers (其他细胞): 必须全部 == 0
#     - Cdh5 / Pecam1 (内皮)
#     - Acta2 (HSC/肌纤维)
#     - Krt19 / Krt7 (胆管)
#     - Epcam (上皮)
#
# 输出 (nam_reannotate_results/):
#   v4_refined_clean.csv         ← cell_id × new_label, 用于下游 join
#   reannotation_summary.csv     ← 各样本/组的真 NAM 数量
#   misannotated_cells.csv       ← 被踢出的 NAM 细胞列表 + 它们高表达的污染 marker
#   verification_markers.csv     ← 重注释后真 NAM 的 marker 验证
#   summary.md                   ← 报告
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(Matrix)
})

select <- dplyr::select
filter <- dplyr::filter

# ============================================================
# 0. 路径
# ============================================================
PROJECT_ROOT <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
SCRIPT_DIR   <- file.path(PROJECT_ROOT, "scripts")
setwd(PROJECT_ROOT)

source(file.path(SCRIPT_DIR, "diff_common.R"))

OUT_DIR <- file.path(PROJECT_ROOT, "nam_reannotate_results")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("========================================\n")
cat("NAM Reannotation: marker-based strict filter\n")
cat("Start:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================\n\n")

t0 <- proc.time()

# ============================================================
# 1. Marker 定义
# ============================================================
POSITIVE_MARKERS <- c(
  "Itgb2",   # CD18, 白细胞 integrin
  "Cd74",    # MHC-II invariant chain
  "Csf1r",   # M-CSF receptor
  "Adgre1",  # F4/80
  "Lgmn",    # asparaginyl endopeptidase
  "Ctss",    # cathepsin S
  "Cd68",    # pan-macrophage
  "Lyz2",    # lysozyme
  "Aif1"     # Iba1
)

# 需要至少 3 个 positive markers 表达 (>0)
MIN_POSITIVE <- 3

NEGATIVE_HEPATOCYTE <- c(
  "Alb",     # 白蛋白
  "Itih4",   # IT-H4
  "Mat1a",   # 甲硫氨酸腺苷转移酶
  "Pzp",     # alpha-2-巨球蛋白
  "Fga",     # 纤维蛋白原 A
  "Apoa1"    # 载脂蛋白
)

# 反向 markers: 加和必须 == 0 (严格), 或 < THRESHOLD (宽松)
NEGATIVE_THRESHOLD <- 0.5  # log1p 后的总和阈值, 0 或 0.5 都可以

EXCLUDE_OTHER <- c(
  "Cdh5", "Pecam1",     # 内皮
  "Acta2",              # HSC/肌纤维
  "Krt19", "Krt7",      # 胆管
  "Epcam"               # 上皮
)

cat("Marker panels:\n")
cat(sprintf("  Positive (need >= %d expressed): %s\n",
            MIN_POSITIVE, paste(POSITIVE_MARKERS, collapse = ", ")))
cat(sprintf("  Negative hepatocyte (sum must be < %.1f): %s\n",
            NEGATIVE_THRESHOLD, paste(NEGATIVE_HEPATOCYTE, collapse = ", ")))
cat(sprintf("  Negative other (sum must be < %.1f): %s\n\n",
            NEGATIVE_THRESHOLD, paste(EXCLUDE_OTHER, collapse = ", ")))

# ============================================================
# 2. 加载 liver 样本
# ============================================================
LIVER_SAMPLES_BY_GROUP <- list(
  Tumor   = c("liver3", "liver8", "gan56"),
  Blank   = c("bla_liver"),
  Control = c("con_liver")
)

get_sample_group <- function(sample) {
  for (g in names(LIVER_SAMPLES_BY_GROUP)) {
    if (sample %in% LIVER_SAMPLES_BY_GROUP[[g]]) return(g)
  }
  return(NA)
}

cat("========== Loading liver samples ==========\n")
samples <- unlist(LIVER_SAMPLES_BY_GROUP)
objs <- list()
for (s in samples) {
  fp <- resolve_rds_path(s, "v4")
  if (!file.exists(fp)) {
    cat(sprintf("  ⚠ missing %s\n", s)); next
  }
  o <- tryCatch({
    x <- readRDS(fp)
    if (inherits(x, "Seurat")) x <- tryCatch(JoinLayers(x), error = function(e) x)
    x <- RenameCells(x, add.cell.id = s)
    x$sample <- s
    x$group  <- get_sample_group(s)
    x
  }, error = function(e) { cat(sprintf("  ⚠ %s: %s\n", s, e$message)); NULL })
  if (!is.null(o)) {
    objs[[s]] <- o
    cat(sprintf("  ✓ %-12s n=%d\n", s, ncol(o)))
  }
}

cat("  merging...\n")
seurat_liver <- merge(objs[[1]], objs[-1])
seurat_liver <- tryCatch(JoinLayers(seurat_liver), error = function(e) seurat_liver)
cat(sprintf("  merged n=%d cells\n\n", ncol(seurat_liver)))
rm(objs); gc(verbose = FALSE)

# ============================================================
# 3. 提取表达矩阵
# ============================================================
cat("========== Extracting expression matrix ==========\n")

# 用 normalized data (RNA layer 'data')
seurat_liver <- NormalizeData(seurat_liver, verbose = FALSE)
data_mat <- GetAssayData(seurat_liver, assay = "RNA", layer = "data")

cat(sprintf("  expression matrix: %d genes × %d cells\n",
            nrow(data_mat), ncol(data_mat)))

# 检查 marker 是否在数据里
pos_in_data <- intersect(POSITIVE_MARKERS, rownames(data_mat))
neg_in_data <- intersect(NEGATIVE_HEPATOCYTE, rownames(data_mat))
exc_in_data <- intersect(EXCLUDE_OTHER, rownames(data_mat))

cat(sprintf("  Positive markers in data: %d/%d (%s)\n",
            length(pos_in_data), length(POSITIVE_MARKERS),
            paste(setdiff(POSITIVE_MARKERS, pos_in_data), collapse = ",")))
cat(sprintf("  Negative hepatocyte markers in data: %d/%d (%s)\n",
            length(neg_in_data), length(NEGATIVE_HEPATOCYTE),
            paste(setdiff(NEGATIVE_HEPATOCYTE, neg_in_data), collapse = ",")))
cat(sprintf("  Exclude other markers in data: %d/%d (%s)\n\n",
            length(exc_in_data), length(EXCLUDE_OTHER),
            paste(setdiff(EXCLUDE_OTHER, exc_in_data), collapse = ",")))

if (length(pos_in_data) < MIN_POSITIVE) {
  stop("Positive markers 数量不足 MIN_POSITIVE, 调整 marker 列表")
}

# ============================================================
# 4. 严格筛选: 全 liver 细胞上找真 NAM
# ============================================================
cat("========== Strict NAM identification ==========\n")

# 4.1 计算每个细胞的 positive marker 表达计数 (有多少个 marker > 0)
pos_mat <- data_mat[pos_in_data, , drop = FALSE]
n_pos_expressed <- colSums(pos_mat > 0)  # 每个细胞表达的 positive marker 数

# 4.2 计算每个细胞的 negative hepatocyte 加和
neg_mat <- data_mat[neg_in_data, , drop = FALSE]
neg_sum <- colSums(neg_mat)  # log1p 后的加和

# 4.3 计算每个细胞的 exclude other 加和
exc_mat <- data_mat[exc_in_data, , drop = FALSE]
exc_sum <- colSums(exc_mat)

# 4.4 综合判定
is_real_nam <- (n_pos_expressed >= MIN_POSITIVE) &
               (neg_sum < NEGATIVE_THRESHOLD) &
               (exc_sum < NEGATIVE_THRESHOLD)

n_real_nam <- sum(is_real_nam)
cat(sprintf("  Total real NAM cells (strict): %d\n", n_real_nam))

# 4.5 分布到样本和组
real_nam_meta <- seurat_liver@meta.data[is_real_nam, c("sample", "group")]
cat("\n  Distribution by sample:\n")
print(table(real_nam_meta$sample))
cat("\n  Distribution by group:\n")
print(table(real_nam_meta$group))

# 4.6 与原 v4_refined NAM 标签对比
ct_col <- NULL
for (candidate in c("celltype_refined", "v4_refined", "final_celltype")) {
  if (candidate %in% colnames(seurat_liver@meta.data)) {
    ct_col <- candidate; break
  }
}
if (is.null(ct_col)) stop("找不到细胞类型列")

cat(sprintf("\n  using cell type column: %s\n", ct_col))
old_nam_labels <- LIVER_CELL_TYPES[["Macrophage_NAM"]]
cat(sprintf("  v4_refined NAM labels: %s\n",
            paste(old_nam_labels, collapse = ", ")))

was_old_nam <- seurat_liver@meta.data[[ct_col]] %in% old_nam_labels
n_was_old_nam <- sum(was_old_nam)
cat(sprintf("\n  Was v4 NAM: %d\n", n_was_old_nam))
cat(sprintf("  Is now real NAM: %d\n", sum(is_real_nam & was_old_nam)))
cat(sprintf("  Was v4 NAM but now hepatocyte/other: %d (kicked out)\n",
            sum(was_old_nam & !is_real_nam)))
cat(sprintf("  NEW NAM not in v4 NAM: %d (newly added)\n",
            sum(is_real_nam & !was_old_nam)))

# ============================================================
# 5. 构造新注释列 v4_refined_clean
# ============================================================
cat("\n========== Building v4_refined_clean ==========\n")

# 默认: 完全沿用原注释
new_labels <- as.character(seurat_liver@meta.data[[ct_col]])

# 把 was_old_nam 但 !is_real_nam 的细胞标记为 'Macrophage_NAM_misannotated'
# (这些细胞 173 个左右, 主要是 hepatocyte contamination)
to_kick <- was_old_nam & !is_real_nam
new_labels[to_kick] <- "Macrophage_NAM_misannotated"

# 把 newly identified real NAM 标记为 'Macrophage_NAM_clean'
# (这包括:
#   1) was_old_nam & is_real_nam (原 NAM 中真的部分)
#   2) !was_old_nam & is_real_nam (新发现的真 NAM)
# 都统一叫 Macrophage_NAM_clean
to_promote <- is_real_nam
new_labels[to_promote] <- "Macrophage_NAM_clean"

cat(sprintf("  Final v4_refined_clean distribution for NAM-related labels:\n"))
nam_related_idx <- new_labels %in% c("Macrophage_NAM_clean",
                                      "Macrophage_NAM_misannotated",
                                      old_nam_labels)
print(table(new_labels[nam_related_idx]))

# ============================================================
# 6. 输出 cell-level CSV
# ============================================================
cat("\n========== Writing v4_refined_clean.csv ==========\n")

cell_anno <- data.frame(
  cell_id = colnames(seurat_liver),
  sample  = seurat_liver$sample,
  group   = seurat_liver$group,
  organ   = "liver",
  v4_refined_old   = as.character(seurat_liver@meta.data[[ct_col]]),
  v4_refined_clean = new_labels,
  n_pos_markers_expr = n_pos_expressed,
  neg_hepatocyte_sum = round(neg_sum, 3),
  neg_other_sum      = round(exc_sum, 3),
  stringsAsFactors = FALSE
)

# 只保存有注释变更或属于 NAM-related 的细胞 (减小文件)
cell_anno_changed <- cell_anno[cell_anno$v4_refined_old != cell_anno$v4_refined_clean |
                                cell_anno$v4_refined_clean %in%
                                  c("Macrophage_NAM_clean", "Macrophage_NAM_misannotated"), ]

write.csv(cell_anno, file.path(OUT_DIR, "v4_refined_clean_full.csv"),
          row.names = FALSE)
write.csv(cell_anno_changed, file.path(OUT_DIR, "v4_refined_clean.csv"),
          row.names = FALSE)
cat(sprintf("  Full annotation: %d cells → v4_refined_clean_full.csv\n",
            nrow(cell_anno)))
cat(sprintf("  Changed/NAM-related: %d cells → v4_refined_clean.csv\n",
            nrow(cell_anno_changed)))

# ============================================================
# 7. 输出 misannotated cells 详情
# ============================================================
cat("\n========== Writing misannotated_cells.csv ==========\n")

mis_idx <- which(was_old_nam & !is_real_nam)
if (length(mis_idx) > 0) {
  mis_df <- data.frame(
    cell_id = colnames(seurat_liver)[mis_idx],
    sample  = seurat_liver$sample[mis_idx],
    group   = seurat_liver$group[mis_idx],
    n_pos_markers_expr = n_pos_expressed[mis_idx],
    neg_hepatocyte_sum = round(neg_sum[mis_idx], 3),
    neg_other_sum      = round(exc_sum[mis_idx], 3),
    Alb     = round(if ("Alb" %in% rownames(data_mat)) data_mat["Alb", mis_idx] else 0, 3),
    Itih4   = round(if ("Itih4" %in% rownames(data_mat)) data_mat["Itih4", mis_idx] else 0, 3),
    Mat1a   = round(if ("Mat1a" %in% rownames(data_mat)) data_mat["Mat1a", mis_idx] else 0, 3),
    Pzp     = round(if ("Pzp" %in% rownames(data_mat)) data_mat["Pzp", mis_idx] else 0, 3),
    Fga     = round(if ("Fga" %in% rownames(data_mat)) data_mat["Fga", mis_idx] else 0, 3),
    Csf1r   = round(if ("Csf1r" %in% rownames(data_mat)) data_mat["Csf1r", mis_idx] else 0, 3),
    Cd68    = round(if ("Cd68" %in% rownames(data_mat)) data_mat["Cd68", mis_idx] else 0, 3),
    Adgre1  = round(if ("Adgre1" %in% rownames(data_mat)) data_mat["Adgre1", mis_idx] else 0, 3),
    stringsAsFactors = FALSE
  )
  mis_df <- mis_df[order(-mis_df$neg_hepatocyte_sum), ]
  write.csv(mis_df, file.path(OUT_DIR, "misannotated_cells.csv"),
            row.names = FALSE)
  cat(sprintf("  %d misannotated cells written\n", nrow(mis_df)))

  # 摘要
  cat("\n  Misannotated 细胞的肝细胞 marker 表达汇总:\n")
  cat(sprintf("    Mean Alb expression: %.2f\n", mean(mis_df$Alb)))
  cat(sprintf("    Mean Itih4 expression: %.2f\n", mean(mis_df$Itih4)))
  cat(sprintf("    Mean Mat1a expression: %.2f\n", mean(mis_df$Mat1a)))
  cat(sprintf("    Mean Csf1r expression: %.2f (should be low)\n", mean(mis_df$Csf1r)))
  cat(sprintf("    Mean Cd68 expression: %.2f (should be low)\n", mean(mis_df$Cd68)))
}

# ============================================================
# 8. 验证: 重注释后的 NAM 是不是真巨噬?
# ============================================================
cat("\n========== Verification: real NAM marker expression ==========\n")

real_nam_idx <- which(is_real_nam)
if (length(real_nam_idx) >= 10) {
  # 计算真 NAM 上 marker 的平均表达
  marker_panel <- c(POSITIVE_MARKERS, NEGATIVE_HEPATOCYTE, EXCLUDE_OTHER)
  marker_in_data <- intersect(marker_panel, rownames(data_mat))

  marker_expr_in_real_nam <- data.frame(
    gene = marker_in_data,
    mean_in_real_nam = round(rowMeans(data_mat[marker_in_data, real_nam_idx,
                                                drop = FALSE]), 3),
    pct_in_real_nam = round(rowMeans(data_mat[marker_in_data, real_nam_idx,
                                                drop = FALSE] > 0) * 100, 1),
    category = ifelse(marker_in_data %in% POSITIVE_MARKERS, "Positive_myeloid",
               ifelse(marker_in_data %in% NEGATIVE_HEPATOCYTE, "Negative_hepatocyte",
                      "Negative_other")),
    stringsAsFactors = FALSE
  )
  marker_expr_in_real_nam <- marker_expr_in_real_nam[
    order(marker_expr_in_real_nam$category, -marker_expr_in_real_nam$mean_in_real_nam), ]
  write.csv(marker_expr_in_real_nam,
            file.path(OUT_DIR, "verification_markers.csv"),
            row.names = FALSE)

  cat("  Marker expression in cleaned real NAM:\n")
  print(marker_expr_in_real_nam)
}

# ============================================================
# 9. 各样本/组 NAM 数量汇总
# ============================================================
cat("\n========== Reannotation summary by sample ==========\n")

summary_df <- data.frame(
  sample = character(),
  group = character(),
  v4_NAM_count = integer(),
  real_NAM_count = integer(),
  misannotated_count = integer(),
  newly_added_NAM = integer(),
  stringsAsFactors = FALSE
)

for (s in unique(seurat_liver$sample)) {
  cell_in_sample <- which(seurat_liver$sample == s)
  v4_nam_in_s <- sum(was_old_nam[cell_in_sample])
  real_nam_in_s <- sum(is_real_nam[cell_in_sample])
  mis_in_s <- sum((was_old_nam & !is_real_nam)[cell_in_sample])
  new_in_s <- sum((is_real_nam & !was_old_nam)[cell_in_sample])
  grp <- get_sample_group(s)

  summary_df <- rbind(summary_df, data.frame(
    sample = s, group = grp,
    v4_NAM_count = v4_nam_in_s,
    real_NAM_count = real_nam_in_s,
    misannotated_count = mis_in_s,
    newly_added_NAM = new_in_s,
    stringsAsFactors = FALSE
  ))
}
write.csv(summary_df, file.path(OUT_DIR, "reannotation_summary.csv"),
          row.names = FALSE)
print(summary_df)

# ============================================================
# 10. 写 summary md
# ============================================================
cat("\n========== Writing summary.md ==========\n")

elapsed_min <- (proc.time() - t0)["elapsed"] / 60

n_v4_nam <- sum(was_old_nam)
n_kept   <- sum(was_old_nam & is_real_nam)
n_kicked <- sum(was_old_nam & !is_real_nam)
n_added  <- sum(!was_old_nam & is_real_nam)
n_final  <- sum(is_real_nam)

md_lines <- c(
  "# NAM 重注释报告",
  "",
  sprintf("**生成时间:** %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  sprintf("**耗时:** %.1f min", elapsed_min),
  "",
  "## 背景",
  "",
  "v4_refined 的 'Macrophage_NAM' 标签被发现严重污染了肝细胞 (见 nam_subtype_check_results/verdict.md)。",
  "本脚本用 marker-based 严格筛选重新识别真 NAM, 不修改原 RDS, 输出 cell-level 注释 CSV。",
  "",
  "## Marker 定义",
  "",
  sprintf("**Positive (myeloid 谱系)**: 至少 %d 个 > 0", MIN_POSITIVE),
  paste0("  ", paste(POSITIVE_MARKERS, collapse = ", ")),
  "",
  sprintf("**Negative (hepatocyte)**: 加和 < %.1f", NEGATIVE_THRESHOLD),
  paste0("  ", paste(NEGATIVE_HEPATOCYTE, collapse = ", ")),
  "",
  sprintf("**Negative (其他)**: 加和 < %.1f", NEGATIVE_THRESHOLD),
  paste0("  ", paste(EXCLUDE_OTHER, collapse = ", ")),
  "",
  "## 数量变化",
  "",
  sprintf("- **原 v4 NAM**: %d", n_v4_nam),
  sprintf("- **保留 (真 NAM)**: %d", n_kept),
  sprintf("- **踢出 (实际是 hepatocyte 等)**: %d", n_kicked),
  sprintf("- **新发现 (v4 没标但实际是 NAM)**: %d", n_added),
  sprintf("- **最终真 NAM 总数**: %d", n_final),
  "",
  sprintf("**踢出比例**: %.1f%% (%d/%d)",
          100 * n_kicked / max(1, n_v4_nam), n_kicked, n_v4_nam),
  "",
  "## 各样本/组分布",
  "",
  "| sample | group | v4 NAM | real NAM | misannot | new added |",
  "|---|---|---:|---:|---:|---:|"
)
for (i in seq_len(nrow(summary_df))) {
  r <- summary_df[i, ]
  md_lines <- c(md_lines,
    sprintf("| %s | %s | %d | %d | %d | %d |",
            r$sample, r$group,
            r$v4_NAM_count, r$real_NAM_count,
            r$misannotated_count, r$newly_added_NAM))
}

md_lines <- c(md_lines, "",
  "## 输出文件",
  "",
  "- `v4_refined_clean.csv` — 包含变更和 NAM-related 细胞的注释 (用于下游 join)",
  "- `v4_refined_clean_full.csv` — 完整的 cell-level 注释 (备份)",
  "- `misannotated_cells.csv` — 被踢出的细胞 + 它们的肝细胞 marker 表达",
  "- `verification_markers.csv` — 真 NAM 上的 marker 表达验证",
  "- `reannotation_summary.csv` — 各样本数量汇总",
  "",
  "## 下游脚本怎么用",
  "",
  "下游脚本 (DoRothEA, NAM 拟时序, CellChat) 需要在加载 RDS 后做这一步:",
  "",
  "```r",
  "# 1. 加载 cell-level 注释",
  "anno <- read.csv('nam_reannotate_results/v4_refined_clean_full.csv')",
  "",
  "# 2. 把新注释 join 到 metadata",
  "match_idx <- match(colnames(seurat_obj), anno$cell_id)",
  "seurat_obj$v4_refined_clean <- anno$v4_refined_clean[match_idx]",
  "",
  "# 3. 用 v4_refined_clean 替代 v4_refined / celltype_refined",
  "# 4. 把 'Macrophage_NAM' 替换为 'Macrophage_NAM_clean'",
  "# 5. 'Macrophage_NAM_misannotated' 视为待排除或独立分析",
  "```",
  "",
  "## 接下来的工作",
  "",
  "1. **重跑 DoRothEA**: NAM 相关对比, 用 Macrophage_NAM_clean",
  "2. **重跑 NAM 拟时序**: 验证 monocyte→成熟方向是否依然成立",
  "3. **跑 CellChat**: 用 v4_refined_clean 注释",
  "4. **更新 v2.1 报告 §2.5**: 修订 NAM 部分的所有结论"
)

writeLines(md_lines, file.path(OUT_DIR, "summary.md"))

# ============================================================
# 11. 完成
# ============================================================
cat(sprintf("\n========================================\n"))
cat(sprintf("NAM reannotation complete\n"))
cat(sprintf("  Output dir : %s\n", OUT_DIR))
cat(sprintf("  Real NAM   : %d (was %d, kicked %d, added %d)\n",
            n_final, n_v4_nam, n_kicked, n_added))
cat(sprintf("  Elapsed    : %.1f min\n", elapsed_min))
cat(sprintf("========================================\n"))
