#!/usr/bin/env Rscript
# ============================================================
# nam_clean_v2.R --- 在 v4_refined 的 298 个 NAM 内部清洗肝细胞污染
# ============================================================
# 用法:
#   cd ~/SGS10_analysis
#   Rscript scripts/nam_clean_v2.R 2>&1 | tee scripts/logs/nam_clean_v2.log
#
# 设计:
#   1. 加载 5 个 liver 样本, merge
#   2. 找原 v4_refined 里的 Macrophage_NAM 细胞 (应该 298 个)
#   3. 用 hepatocyte 反向 marker 加和 < threshold 清洗
#   4. 不在全 liver 找新 NAM (避免上次的错误)
#   5. 输出 cell-level CSV + marker 验证
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(Matrix)
})

PROJECT_ROOT <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
setwd(PROJECT_ROOT)
source("scripts/diff_common.R")

OUT_DIR <- file.path(PROJECT_ROOT, "nam_clean_results")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("========================================\n")
cat("NAM Clean v2: 298 internal hepatocyte filter\n")
cat("========================================\n\n")

# ============================================================
# 1. 加载
# ============================================================
LIVER_SAMPLES <- c("liver3", "liver8", "gan56", "bla_liver", "con_liver")
LIVER_GROUP   <- c(liver3 = "Tumor", liver8 = "Tumor", gan56 = "Tumor",
                   bla_liver = "Blank", con_liver = "Control")

cat("Loading liver samples...\n")
objs <- list()
for (s in LIVER_SAMPLES) {
  fp <- resolve_rds_path(s, "v4")
  o <- readRDS(fp)
  if (inherits(o, "Seurat")) {
    o <- tryCatch(JoinLayers(o), error = function(e) o)
  }
  o <- RenameCells(o, add.cell.id = s)
  # 直接用 [[<- 操作而不是 $<- (避免对齐问题)
  o[["sample"]] <- s
  o[["group"]]  <- LIVER_GROUP[[s]]
  objs[[s]] <- o
  cat(sprintf("  %-12s n=%d\n", s, ncol(o)))
}

cat("merging...\n")
sl <- merge(objs[[1]], objs[-1])
sl <- tryCatch(JoinLayers(sl), error = function(e) sl)
rm(objs); gc(verbose = FALSE)
cat(sprintf("merged: %d cells\n", ncol(sl)))

# ============================================================
# 2. 找 v4 NAM 细胞
# ============================================================
ct_col <- NULL
for (cand in c("celltype_refined", "v4_refined", "final_celltype")) {
  if (cand %in% colnames(sl@meta.data)) { ct_col <- cand; break }
}
if (is.null(ct_col)) stop("No celltype column found")
cat(sprintf("\nUsing celltype column: %s\n", ct_col))

nam_mask <- sl@meta.data[[ct_col]] == "Macrophage_NAM"
nam_mask[is.na(nam_mask)] <- FALSE
n_v4_nam <- sum(nam_mask)
cat(sprintf("v4 NAM cells: %d\n", n_v4_nam))

if (n_v4_nam == 0) stop("No v4 NAM cells found")

# ============================================================
# 3. 提取表达矩阵 + 反向筛选
# ============================================================
cat("\nNormalizing...\n")
sl <- NormalizeData(sl, verbose = FALSE)
data_mat <- GetAssayData(sl, assay = "RNA", layer = "data")
cat(sprintf("data matrix: %d genes × %d cells\n", nrow(data_mat), ncol(data_mat)))

NEG_HEP <- c("Alb", "Itih4", "Mat1a", "Pzp", "Fga", "Apoa1")
NEG_HEP_in <- intersect(NEG_HEP, rownames(data_mat))
cat(sprintf("hepatocyte markers used: %s\n",
            paste(NEG_HEP_in, collapse = ",")))

# 提取 NAM 细胞的列索引
nam_idx <- which(nam_mask)

# 计算 negative hepatocyte 加和 (只在 NAM 细胞上)
neg_sum <- colSums(data_mat[NEG_HEP_in, nam_idx, drop = FALSE])
cat(sprintf("\nneg_hep_sum stats on %d NAM cells:\n", length(neg_sum)))
print(summary(neg_sum))

# ============================================================
# 4. 阈值扫描
# ============================================================
cat("\n=== Threshold scan ===\n")
cat("threshold | n_keep | pct_keep\n")
cat("----------|--------|----------\n")
for (thresh in c(0, 0.1, 0.3, 0.5, 1.0, 2.0, 5.0)) {
  n_k <- sum(neg_sum < thresh)
  cat(sprintf("  %.2f    | %5d  | %5.1f%%\n",
              thresh, n_k, 100 * n_k / length(neg_sum)))
}

# ============================================================
# 5. 用 < 0.5 阈值生成最终列表
# ============================================================
THRESH <- 0.5
keep_mask <- neg_sum < THRESH
keep_cells_global_idx <- nam_idx[keep_mask]
kick_cells_global_idx <- nam_idx[!keep_mask]

n_keep <- length(keep_cells_global_idx)
n_kick <- length(kick_cells_global_idx)

cat(sprintf("\n=== Final at threshold %.1f ===\n", THRESH))
cat(sprintf("kept (clean NAM):     %d (%.1f%%)\n",
            n_keep, 100 * n_keep / n_v4_nam))
cat(sprintf("kicked (hepatocyte):  %d (%.1f%%)\n",
            n_kick, 100 * n_kick / n_v4_nam))

# ============================================================
# 6. 各样本/组分布
# ============================================================
cat("\n=== Distribution of CLEAN NAM ===\n")
keep_meta <- sl@meta.data[keep_cells_global_idx, c("sample", "group")]
cat("\nBy sample:\n")
print(table(keep_meta$sample))
cat("\nBy group:\n")
print(table(keep_meta$group))

cat("\n=== Distribution of KICKED (hepatocyte contamination) ===\n")
kick_meta <- sl@meta.data[kick_cells_global_idx, c("sample", "group")]
cat("\nBy sample:\n")
print(table(kick_meta$sample))
cat("\nBy group:\n")
print(table(kick_meta$group))

# ============================================================
# 7. Marker 验证
# ============================================================
cat("\n=== Marker verification on CLEAN NAM ===\n")
verify_genes <- c(
  # 核心 myeloid
  "Cd68", "Adgre1", "Csf1r", "Cd74", "Ctss", "Lyz2", "Itgb2",
  "Aif1", "Lgmn", "Spi1", "Cebpb",
  # 反向 hepatocyte
  "Alb", "Itih4", "Mat1a", "Pzp", "Fga", "Apoa1",
  # NAM-relevant TF (DoRothEA hits)
  "Runx1", "Mafb"
)
verify_in <- intersect(verify_genes, rownames(data_mat))
cat(sprintf("verify markers in data: %d\n", length(verify_in)))

if (n_keep >= 5) {
  vmat <- data_mat[verify_in, keep_cells_global_idx, drop = FALSE]
  verify_df <- data.frame(
    gene = verify_in,
    mean = round(rowMeans(vmat), 3),
    pct  = round(rowMeans(vmat > 0) * 100, 1),
    stringsAsFactors = FALSE
  )
  verify_df$category <- ifelse(verify_df$gene %in% NEG_HEP_in,
                                "Negative_hep", "Positive_myeloid")
  verify_df <- verify_df[order(verify_df$category, -verify_df$mean), ]
  cat("\n")
  print(verify_df)
} else {
  cat("Too few cells to verify\n")
}

# ============================================================
# 8. 也看 KICKED 细胞的 marker (确认是肝细胞)
# ============================================================
if (n_kick >= 5) {
  cat("\n=== Marker check on KICKED cells (should be hepatocyte) ===\n")
  vmat_k <- data_mat[verify_in, kick_cells_global_idx, drop = FALSE]
  verify_kdf <- data.frame(
    gene = verify_in,
    mean = round(rowMeans(vmat_k), 3),
    pct  = round(rowMeans(vmat_k > 0) * 100, 1),
    stringsAsFactors = FALSE
  )
  verify_kdf$category <- ifelse(verify_kdf$gene %in% NEG_HEP_in,
                                 "Negative_hep", "Positive_myeloid")
  verify_kdf <- verify_kdf[order(verify_kdf$category, -verify_kdf$mean), ]
  print(verify_kdf)
}

# ============================================================
# 9. 输出 CSV
# ============================================================
cat("\n=== Writing output ===\n")
cell_names <- colnames(sl)
out_df <- data.frame(
  cell_id = cell_names[nam_idx],
  sample  = sl@meta.data$sample[nam_idx],
  group   = sl@meta.data$group[nam_idx],
  v4_refined_old = "Macrophage_NAM",
  v4_refined_clean = ifelse(keep_mask,
                             "Macrophage_NAM_clean",
                             "Macrophage_NAM_misannotated"),
  neg_hep_sum = round(neg_sum, 3),
  stringsAsFactors = FALSE
)

write.csv(out_df,
          file.path(OUT_DIR, "nam_298_classification.csv"),
          row.names = FALSE)
cat(sprintf("  → %s\n", file.path(OUT_DIR, "nam_298_classification.csv")))

# 也写一个完整的 cell-level 注释 (所有细胞), 用于下游 join
full_anno <- data.frame(
  cell_id = cell_names,
  v4_refined_old = as.character(sl@meta.data[[ct_col]]),
  stringsAsFactors = FALSE
)
full_anno$v4_refined_clean <- full_anno$v4_refined_old
full_anno$v4_refined_clean[nam_idx[keep_mask]] <- "Macrophage_NAM_clean"
full_anno$v4_refined_clean[nam_idx[!keep_mask]] <- "Macrophage_NAM_misannotated"

write.csv(full_anno,
          file.path(OUT_DIR, "v4_refined_clean_full.csv"),
          row.names = FALSE)
cat(sprintf("  → %s (%d cells)\n",
            file.path(OUT_DIR, "v4_refined_clean_full.csv"),
            nrow(full_anno)))

cat("\n========================================\n")
cat("Done\n")
cat(sprintf("  v4 NAM:    %d\n", n_v4_nam))
cat(sprintf("  Clean:     %d\n", n_keep))
cat(sprintf("  Kicked:    %d\n", n_kick))
cat("========================================\n")
