#!/usr/bin/env Rscript
# ============================================================
# pseudotime_nam_only.R --- Macrophage_NAM 单独拟时序分析
# ============================================================
# 用法:
#   cd ~/SGS10_analysis
#   Rscript scripts/pseudotime_nam_only.R 2>&1 | tee scripts/logs/pseudotime_nam.log
#
# 动机:
#   原 pseudotime_v4_core.R 排除了 Macrophage_NAM (n 不足判断错误)
#   v2 整合报告显示 NAM 是核心细类:
#     - DoRothEA Spi1 score 8.63, padj 1.8e-16 (全数据集最强 TF 信号)
#     - Smad3 score -2.92, padj 0.017 (排除 NAM 作为 TGFb 源)
#     - Nfkb1 +3.39, Irf1 +2.94 (NFκB/IFN 双驱动)
#   补跑 NAM 拟时序验证: Spi1 程序激活是否沿伪时间渐变
#
# 参数与原脚本的关键区别:
#   1. 只处理 NAM 单一细类
#   2. min_cells 阈值从 500 降到 200 (NAM 在某些样本可能较少)
#   3. Root marker 用 Spi1 (NAM 主 regulator, 高活性细胞作为终点; 起点用低 Spi1)
#      实际上我们用 Ly6c2 / Ccr2 (monocyte 起点 marker) 作为起点
#   4. Core panel 加入 NAM 特异基因 (骨髓单核谱系 TF + chemokine + IFN 通路)
#   5. ENABLE_TRADESEQ = FALSE (复用 v4 脚本的降级方案, 避免 fitGAM 卡顿)
#
# 输出:
#   pseudotime_nam_results/
#     Macrophage_NAM/
#       trajectory.{pdf,png}
#       curve_assignment.{pdf,png}
#       pseudotime_by_group.{pdf,png}      关键
#       spearman_correlation.csv
#       heatmap_top_genes.pdf
#       dynamic_core_genes.{pdf,png}        关键
#       pseudotime_values.csv
#     summary/
#       task_status.csv
#       pseudotime_nam_summary.md
#       error_log.txt
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(slingshot)
  library(harmony)
  library(FNN)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(RColorBrewer)
  library(viridis)
  library(pheatmap)
})

# dplyr namespace lock
select    <- dplyr::select
rename    <- dplyr::rename
filter    <- dplyr::filter
mutate    <- dplyr::mutate
arrange   <- dplyr::arrange
bind_rows <- dplyr::bind_rows

# ============================================================
# 0. 路径 & 共用模块
# ============================================================
PROJECT_ROOT <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
SCRIPT_DIR   <- file.path(PROJECT_ROOT, "scripts")
setwd(PROJECT_ROOT)

source(file.path(SCRIPT_DIR, "diff_common.R"))

OUT_DIR <- file.path(PROJECT_ROOT, "pseudotime_nam_results")
for (d in c("Macrophage_NAM", "summary")) {
  dir.create(file.path(OUT_DIR, d), recursive = TRUE, showWarnings = FALSE)
}

ERROR_LOG <- file.path(OUT_DIR, "summary", "error_log.txt")
writeLines(sprintf("== pseudotime NAM error log - %s ==",
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

VIRTUAL_CT <- "Macrophage_NAM"
ORGAN      <- "liver"

HARMONY_DIMS <- 30
UMAP_DIMS    <- 30
KMEANS_K     <- 12       # 比 10 细类脚本的 15 略小, 适合中等细胞数
KNN_K        <- 10       # NA KNN 填充

MIN_CELLS_FOR_RUN <- 200   # 比原脚本的 500 更宽松

# 内存硬限制 (复用原脚本的防护机制)
MIN_FREE_GB_TO_START <- 200
MIN_FREE_GB_ABORT    <- 100
MAX_R_PROCS_ALLOWED  <- 2

# tradeSeq 开关 (强制 FALSE, 避免再次卡顿)
ENABLE_TRADESEQ <- FALSE

# 核心基因面板 (NAM 特异 + 三轴相关)
CORE_PANEL <- unique(c(
  # ===== NAM 谱系 TF (DoRothEA 强激活的) =====
  "Spi1",     # PU.1, NAM 主 regulator (DoRothEA score 8.63)
  "Runx1",    # padj 3e-7
  "Ikzf1",    # padj 1e-3
  "Pou2f2",   # padj 2e-5
  "Lyl1",     # padj 1e-15 (DoRothEA 最强之一)
  "Cebpb",    # 巨噬谱系
  "Mafb",     # 巨噬成熟 TF

  # ===== 单核细胞起点 marker (作为 trajectory 起点判定) =====
  "Ly6c2",    # 经典单核
  "Ccr2",     # 单核浸润 marker
  "Plac8",    # 单核 marker
  "Ms4a4c",
  "Treml4",

  # ===== 组织驻留巨噬终点 marker =====
  "Adgre1",   # F4/80
  "Csf1r",    # M-CSF receptor
  "Mertk",    # 吞噬受体
  "Cd68",
  "Marco",
  "Vsig4",    # Kupffer-like

  # ===== NFκB / IFN 通路 (DoRothEA 显示在 NAM 激活) =====
  "Nfkb1", "Rela", "Relb", "Nfkbia",
  "Irf1", "Irf7", "Stat1", "Stat2",

  # ===== 分泌型炎症因子 =====
  "Il6", "Il1b", "Tnf", "Cxcl9", "Cxcl10", "Cxcl16",
  "Ccl2", "Ccl5", "Ccl7", "Ccl8",

  # ===== TGFb 通路 (验证 NAM 不是 TGFb 源) =====
  "Tgfb1", "Smad3", "Smad7", "Tgfbr1", "Tgfbr2",

  # ===== Efferocytosis 通路 (轴 2 相关) =====
  "Elmo1", "Dock4", "Rac1", "Itgb5",
  "Lrp1", "Pld4", "Marcks", "Gas6", "Axl",

  # ===== MHC-I / 抗原呈递 =====
  "B2m", "H2-D1", "H2-K1", "H2-Aa", "H2-Ab1",

  # ===== Macrophage 极化 (M1 / M2) =====
  "Nos2", "Arg1", "Mrc1", "Cd163", "Il10"
))

# Root marker: 用 Ly6c2 (单核细胞起点)
# 高 Ly6c2 = 新近浸润的单核; 低 Ly6c2 = 已分化的组织巨噬 / NAM
ROOT_MARKER <- "Ly6c2"

# 组颜色
GROUP_COLORS <- c("Tumor"="#E41A1C", "Blank"="#FF7F00", "Control"="#377EB8")

cat("========================================\n")
cat("Pseudotime: Macrophage_NAM only\n")
cat("Start:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================\n\n")

options(future.globals.maxSize = 256 * 1024^3)

# ============================================================
# 1.5 安全检查 (复用 v4 脚本的逻辑)
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
    stop(sprintf("ABORT: %d R 进程残留 (限制 %d)。请先: pkill -9 -u $USER R",
                 n_r, MAX_R_PROCS_ALLOWED))
  }

  free_gb <- get_free_mem_gb()
  cat(sprintf("  Available memory: %.0f GB (min required: %d GB)\n",
              if (is.na(free_gb)) 0 else free_gb, MIN_FREE_GB_TO_START))
  if (is.na(free_gb) || free_gb < MIN_FREE_GB_TO_START) {
    stop(sprintf("ABORT: 可用内存不足"))
  }

  cat("  ✓ Safety check passed\n\n")
}

preflight_check()

# ============================================================
# 2. sample → group 映射
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

# ============================================================
# 3. 加载 liver v4 RDS 并 merge
# ============================================================
cat("========== Loading liver v4 RDS ==========\n")

samples <- unlist(LIVER_SAMPLES_BY_GROUP)
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
    x$group  <- get_sample_group(s)
    x$organ  <- "liver"
    x
  }, error = function(e) { write_error(s, e$message); NULL })
  if (!is.null(o)) {
    objs[[s]] <- o
    cat(sprintf("  ✓ %-12s  cells=%6d  group=%s\n", s, ncol(o), o$group[1]))
  }
}

if (length(objs) < 2) stop("加载样本数 < 2, 无法 merge")

cat(sprintf("  merging %d liver samples...\n", length(objs)))
seurat_liver <- merge(objs[[1]], objs[-1])
seurat_liver <- tryCatch(JoinLayers(seurat_liver), error = function(e) seurat_liver)
cat(sprintf("  [liver] merged cells: %d\n", ncol(seurat_liver)))
rm(objs); gc(verbose = FALSE)

# ============================================================
# 4. KNN 填充 NA 函数
# ============================================================
fill_pt_na_knn <- function(umap_mat, pseudotime_vec, k = KNN_K) {
  na_idx <- which(is.na(pseudotime_vec))
  nona_idx <- which(!is.na(pseudotime_vec))
  if (length(na_idx) == 0) return(pseudotime_vec)
  if (length(nona_idx) < k) k <- max(1, floor(length(nona_idx) / 2))

  knn_res <- FNN::get.knnx(
    data  = umap_mat[nona_idx, , drop = FALSE],
    query = umap_mat[na_idx,   , drop = FALSE],
    k     = min(k, length(nona_idx))
  )
  filled <- numeric(length(na_idx))
  for (i in seq_along(na_idx)) {
    neighbors <- knn_res$nn.index[i, ]
    dists     <- knn_res$nn.dist[i, ]
    if (all(dists == 0)) dists <- rep(1, length(dists))
    else dists[dists == 0] <- min(dists[dists > 0]) / 10
    w <- 1 / dists; w <- w / sum(w)
    filled[i] <- sum(pseudotime_vec[nona_idx[neighbors]] * w)
  }
  pseudotime_vec[na_idx] <- filled
  pseudotime_vec
}

# ============================================================
# 5. 主分析: 单细类拟时序
# ============================================================
cat("\n========== Running pseudotime on Macrophage_NAM ==========\n")

t_begin <- proc.time()
ct_dir <- file.path(OUT_DIR, "Macrophage_NAM")

result <- list(
  organ = ORGAN, virtual_ct = VIRTUAL_CT,
  status = "error", n_cells = NA, n_curves = NA,
  n_spearman_sig = NA, elapsed_min = NA
)

tryCatch({

  # 5.1 抽出 NAM 细胞
  cells <- extract_cells_v2(seurat_liver, VIRTUAL_CT, "v4_refined", LIVER_CELL_TYPES)
  cat(sprintf("  Cells extracted: %d\n", length(cells)))

  if (length(cells) < MIN_CELLS_FOR_RUN) {
    stop(sprintf("Cells (%d) < MIN_CELLS_FOR_RUN (%d)",
                 length(cells), MIN_CELLS_FOR_RUN))
  }

  sub <- subset(seurat_liver, cells = cells)
  sub <- tryCatch(JoinLayers(sub), error = function(e) sub)
  result$n_cells <- ncol(sub)

  cat(sprintf("  [NAM] n_cells=%d, groups: %s\n",
              ncol(sub),
              paste0(names(table(sub$group)), "=",
                     as.vector(table(sub$group)), collapse = ", ")))

  # 释放 parent (NAM 不需要 parent 了)
  rm(seurat_liver); gc(verbose = FALSE)

  # 5.2 子集内 PCA / Harmony / UMAP
  cat("  Normalize → VariableFeatures → Scale → PCA → Harmony → UMAP\n")
  sub <- NormalizeData(sub, verbose = FALSE)
  sub <- FindVariableFeatures(sub, nfeatures = 2000, verbose = FALSE)
  sub <- ScaleData(sub, verbose = FALSE)
  sub <- RunPCA(sub, npcs = 30, verbose = FALSE)
  sub <- tryCatch(
    harmony::RunHarmony(sub, group.by.vars = "sample", verbose = FALSE),
    error = function(e) { write_error("NAM_harmony", e$message); sub }
  )
  red_use <- if ("harmony" %in% names(sub@reductions)) "harmony" else "pca"
  sub <- RunUMAP(sub, reduction = red_use, dims = 1:min(20, UMAP_DIMS),
                 verbose = FALSE)

  umap_mat <- Embeddings(sub, "umap")
  colnames(umap_mat) <- c("UMAP1", "UMAP2")

  # 5.3 kmeans clustering
  set.seed(42)
  n_clust <- min(KMEANS_K, max(3, floor(ncol(sub) / 100)))
  kc <- kmeans(umap_mat, centers = n_clust, iter.max = 50)
  sub$sling_cluster <- as.factor(kc$cluster)
  cat(sprintf("  kmeans clusters: %d\n", n_clust))

  # 5.4 转 SCE
  counts <- GetAssayData(sub, assay = "RNA", layer = "counts")
  if (!inherits(counts, "dgCMatrix")) counts <- as(counts, "dgCMatrix")
  sce <- SingleCellExperiment(
    assays  = list(counts = counts),
    colData = sub@meta.data
  )
  reducedDim(sce, "UMAP") <- umap_mat

  # 5.5 找 root cluster (基于 Ly6c2)
  start_cluster <- NULL
  if (ROOT_MARKER %in% rownames(counts)) {
    marker_expr <- as.numeric(counts[ROOT_MARKER, ])
    cluster_mean <- tapply(marker_expr, sub$sling_cluster, mean)
    start_cluster <- as.character(which.max(cluster_mean))
    cat(sprintf("  root marker=%s, start_cluster=%s (mean_expr=%.3f)\n",
                ROOT_MARKER, start_cluster, max(cluster_mean)))
  } else {
    cat(sprintf("  root marker %s 不在数据里, slingshot 自动选\n", ROOT_MARKER))
  }

  # 5.6 运行 slingshot
  sce <- slingshot::slingshot(
    sce,
    clusterLabels = sub$sling_cluster,
    reducedDim    = "UMAP",
    start.clus    = start_cluster,
    allow.breaks  = TRUE
  )

  pt_mat <- slingPseudotime(sce)
  cw_mat <- slingCurveWeights(sce)
  if (is.null(pt_mat)) stop("slingshot 未生成 pseudotime")

  n_curves <- ncol(pt_mat)
  result$n_curves <- n_curves
  cat(sprintf("  detected %d curves\n", n_curves))

  # 5.7 NA 填充 (主曲线)
  pt_main <- pt_mat[, 1]
  na_pct <- 100 * sum(is.na(pt_main)) / length(pt_main)
  cat(sprintf("  main curve NA: %.1f%%\n", na_pct))
  if (na_pct > 0) {
    pt_main <- fill_pt_na_knn(umap_mat, pt_main, k = KNN_K)
  }
  sub$pseudotime <- pt_main
  sce$pseudotime <- pt_main

  # 5.8 可视化: trajectory
  ct_col_sub <- NULL
  for (candidate in c("celltype_refined", "v4_refined", "final_celltype")) {
    if (candidate %in% colnames(sub@meta.data)) { ct_col_sub <- candidate; break }
  }

  df_traj <- data.frame(
    UMAP1      = umap_mat[, 1],
    UMAP2      = umap_mat[, 2],
    pseudotime = pt_main,
    group      = sub$group,
    sample     = sub$sample,
    celltype   = if (!is.null(ct_col_sub)) as.character(sub@meta.data[[ct_col_sub]])
                 else rep(VIRTUAL_CT, ncol(sub)),
    sling_cluster = as.character(sub$sling_cluster)
  )

  p_traj <- ggplot(df_traj, aes(UMAP1, UMAP2, color = pseudotime)) +
    geom_point(size = 0.6, alpha = 0.85) +
    scale_color_viridis_c(option = "plasma") +
    labs(title = sprintf("liver · Macrophage_NAM · pseudotime"),
         subtitle = sprintf("n_cells=%d, n_curves=%d, root=%s",
                            ncol(sub), n_curves, ROOT_MARKER)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5))

  curves <- slingshot::slingCurves(sce)
  for (cv in curves) {
    cdf <- as.data.frame(cv$s[cv$tag, ])
    colnames(cdf) <- c("UMAP1", "UMAP2")
    p_traj <- p_traj + geom_path(data = cdf,
                                  aes(UMAP1, UMAP2),
                                  color = "black", linewidth = 0.8,
                                  inherit.aes = FALSE)
  }
  ggsave(file.path(ct_dir, "trajectory.pdf"), p_traj, width = 9, height = 7, dpi = 300)
  ggsave(file.path(ct_dir, "trajectory.png"), p_traj, width = 9, height = 7, dpi = 200)

  # 5.9 pseudotime by group (关键!)
  p_by_group <- ggplot(df_traj, aes(x = pseudotime, fill = group)) +
    geom_density(alpha = 0.55) +
    scale_fill_manual(values = GROUP_COLORS) +
    labs(title = "liver · Macrophage_NAM · pseudotime by group",
         subtitle = "tumor vs blank vs control 的伪时间分布") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5))
  ggsave(file.path(ct_dir, "pseudotime_by_group.pdf"),
         p_by_group, width = 9, height = 6, dpi = 300)
  ggsave(file.path(ct_dir, "pseudotime_by_group.png"),
         p_by_group, width = 9, height = 6, dpi = 200)

  # 5.10 curve assignment
  df_traj$curve_assign <- as.character(apply(cw_mat, 1, which.max))
  p_curve <- ggplot(df_traj, aes(UMAP1, UMAP2, color = curve_assign)) +
    geom_point(size = 0.6, alpha = 0.85) +
    scale_color_brewer(palette = "Set1") +
    labs(title = "liver · Macrophage_NAM · curve assignment") +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5))
  ggsave(file.path(ct_dir, "curve_assignment.pdf"),
         p_curve, width = 9, height = 7, dpi = 300)
  ggsave(file.path(ct_dir, "curve_assignment.png"),
         p_curve, width = 9, height = 7, dpi = 200)

  # 导出 pseudotime 值
  df_pt_out <- df_traj
  for (i in seq_len(n_curves)) {
    df_pt_out[[paste0("curve_", i, "_pt")]] <- pt_mat[, i]
  }
  write.csv(df_pt_out, file.path(ct_dir, "pseudotime_values.csv"),
            row.names = FALSE)

  # 5.11 Spearman 全基因组相关
  cat("  Spearman correlation on main curve...\n")
  log_counts <- log1p(as.matrix(counts))
  min_cells_expr <- max(20, ncol(sub) * 0.05)
  genes_keep <- rowSums(counts > 0) >= min_cells_expr
  log_counts_f <- log_counts[genes_keep, ]
  cat(sprintf("  Spearman on %d genes × %d cells\n",
              nrow(log_counts_f), ncol(log_counts_f)))

  cor_df <- data.frame(
    gene = rownames(log_counts_f),
    rho  = NA_real_,
    pval = NA_real_,
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(log_counts_f))) {
    tryCatch({
      ct <- cor.test(log_counts_f[i, ], pt_main, method = "spearman",
                     exact = FALSE)
      cor_df$rho[i]  <- ct$estimate
      cor_df$pval[i] <- ct$p.value
    }, error = function(e) NULL)
  }
  cor_df$padj <- p.adjust(cor_df$pval, method = "BH")
  cor_df <- cor_df[order(abs(cor_df$rho), decreasing = TRUE), ]
  write.csv(cor_df, file.path(ct_dir, "spearman_correlation.csv"),
            row.names = FALSE)

  n_sig <- sum(cor_df$padj < 0.05 & abs(cor_df$rho) > 0.3, na.rm = TRUE)
  result$n_spearman_sig <- n_sig
  cat(sprintf("  Spearman sig (padj<0.05, |rho|>0.3): %d\n", n_sig))

  # 5.12 Top gene heatmap
  top_genes_heatmap <- head(cor_df$gene[cor_df$padj < 0.05 &
                                        abs(cor_df$rho) > 0.3], 80)
  if (length(top_genes_heatmap) >= 10) {
    n_cells_hm <- min(2000, ncol(sub))
    set.seed(42)
    cell_idx_hm <- sample(ncol(sub), n_cells_hm)
    pt_hm <- pt_main[cell_idx_hm]
    ord_hm <- order(pt_hm)
    mat_hm <- log_counts_f[top_genes_heatmap, cell_idx_hm][, ord_hm]
    mat_hm <- t(scale(t(mat_hm)))
    mat_hm[mat_hm > 2] <- 2
    mat_hm[mat_hm < -2] <- -2
    mat_hm <- mat_hm[complete.cases(mat_hm), ]
    if (nrow(mat_hm) >= 5) {
      col_anno <- data.frame(
        pseudotime = pt_hm[ord_hm],
        group      = df_traj$group[cell_idx_hm][ord_hm]
      )
      rownames(col_anno) <- colnames(mat_hm)
      tryCatch({
        pheatmap(mat_hm,
                 cluster_rows = TRUE, cluster_cols = FALSE,
                 annotation_col = col_anno,
                 show_colnames = FALSE,
                 main = sprintf("liver · Macrophage_NAM · top %d Spearman genes",
                                nrow(mat_hm)),
                 filename = file.path(ct_dir, "heatmap_top_genes.pdf"),
                 width = 11, height = 12)
      }, error = function(e) write_error("NAM_heatmap", e$message))
    }
  }

  # 5.13 Core panel dynamic (lowess by group, 关键图)
  core_in_data <- intersect(CORE_PANEL, rownames(log_counts_f))
  cat(sprintf("  core panel genes in data: %d / %d\n",
              length(core_in_data), length(CORE_PANEL)))

  if (length(core_in_data) >= 3) {
    panel_dfs <- list()
    for (g in core_in_data) {
      expr <- log_counts_f[g, ]
      for (grp in unique(df_traj$group)) {
        idx <- which(df_traj$group == grp)
        if (length(idx) < 20) next
        ord <- idx[order(pt_main[idx])]
        sm <- tryCatch(
          stats::lowess(pt_main[ord], expr[ord], f = 0.2),
          error = function(e) NULL
        )
        if (!is.null(sm)) {
          panel_dfs[[length(panel_dfs) + 1]] <- data.frame(
            gene = g, group = grp,
            pseudotime = sm$x, expression = sm$y
          )
        }
      }
    }
    if (length(panel_dfs) > 0) {
      pdf_panel <- dplyr::bind_rows(panel_dfs)
      p_dyn <- ggplot(pdf_panel,
                      aes(pseudotime, expression, color = group)) +
        geom_line(linewidth = 0.8) +
        facet_wrap(~ gene, scales = "free_y", ncol = 5) +
        scale_color_manual(values = GROUP_COLORS) +
        labs(title = "liver · Macrophage_NAM · core panel dynamics",
             subtitle = "NAM 谱系 TF + 炎症通路 + 三轴相关基因",
             x = "Pseudotime", y = "Log1p(Expression) lowess") +
        theme_bw(base_size = 9) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              plot.subtitle = element_text(hjust = 0.5),
              strip.text = element_text(size = 7),
              legend.position = "bottom")
      n_genes_pdyn <- length(unique(pdf_panel$gene))
      ncol_dyn <- 5
      nrow_dyn <- ceiling(n_genes_pdyn / ncol_dyn)
      ggsave(file.path(ct_dir, "dynamic_core_genes.pdf"),
             p_dyn,
             width = ncol_dyn * 3.0,
             height = nrow_dyn * 1.8, dpi = 300, limitsize = FALSE)
      ggsave(file.path(ct_dir, "dynamic_core_genes.png"),
             p_dyn,
             width = ncol_dyn * 3.0,
             height = nrow_dyn * 1.8, dpi = 200, limitsize = FALSE)
    }
  }

  result$status <- "ok"
  result$elapsed_min <- round((proc.time() - t_begin)["elapsed"] / 60, 1)
  cat(sprintf("\n  ✓ DONE  n=%d  curves=%d  sp_sig=%s  (%.1f min)\n",
              result$n_cells, result$n_curves,
              result$n_spearman_sig, result$elapsed_min))

}, error = function(e) {
  result$elapsed_min <<- round((proc.time() - t_begin)["elapsed"] / 60, 1)
  write_error("NAM", e$message)
  cat("  ✗ ERROR:", e$message, "\n")
})

gc(verbose = FALSE)

# ============================================================
# 6. 写 task_status + summary
# ============================================================
cat("\n========== Writing summary ==========\n")

ts_df <- data.frame(
  key            = "liver_Macrophage_NAM",
  organ          = "liver",
  virtual_ct     = "Macrophage_NAM",
  status         = result$status,
  n_cells        = if (is.null(result$n_cells)) NA else result$n_cells,
  n_curves       = if (is.null(result$n_curves)) NA else result$n_curves,
  n_spearman_sig = if (is.null(result$n_spearman_sig)) NA else result$n_spearman_sig,
  elapsed_min    = if (is.null(result$elapsed_min)) NA else result$elapsed_min,
  stringsAsFactors = FALSE
)
write.csv(ts_df, file.path(OUT_DIR, "summary", "task_status.csv"),
          row.names = FALSE)
print(ts_df)

summary_lines <- c(
  "# Pseudotime: Macrophage_NAM 单独分析",
  "",
  sprintf("**Generated:** %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  sprintf("**Elapsed:** %s min",
          if (is.na(result$elapsed_min)) "-" else as.character(result$elapsed_min)),
  "",
  "## 动机",
  "",
  "原 `pseudotime_v4_core.R` 排除了 Macrophage_NAM (n 不足判断)。",
  "v2 整合报告显示 NAM 是核心细类:",
  "",
  "- DoRothEA Spi1 score 8.63, padj 1.8e-16 (全数据集最强 TF 信号)",
  "- Smad3 score -2.92, padj 0.017 (排除 NAM 作为 TGFb 源)",
  "- Nfkb1 +3.39, Irf1 +2.94 (NFκB/IFN 双驱动)",
  "",
  "本分析补跑 NAM 拟时序, 验证 Spi1 程序激活是否沿伪时间渐变。",
  "",
  "## 方法",
  "",
  "- 5 个肝样本 merge → Harmony (by sample) → subset 出 NAM",
  "- 子集内重新 PCA / Harmony / UMAP",
  "- kmeans clustering (k = 12)",
  "- slingshot trajectory (root = Ly6c2 高表达 cluster, 单核细胞起点)",
  "- KNN (k=10) 填充主曲线 NA",
  "- Spearman 全基因组相关",
  "- Lowess 分组核心基因动态曲线",
  "- 跳过 tradeSeq fitGAM (ENABLE_TRADESEQ=FALSE)",
  "",
  "## 任务状态",
  "",
  "| 项目 | 值 |",
  "|---|---|",
  sprintf("| n_cells | %s |", if (is.na(result$n_cells)) "-" else as.character(result$n_cells)),
  sprintf("| n_curves | %s |", if (is.na(result$n_curves)) "-" else as.character(result$n_curves)),
  sprintf("| Spearman sig (padj<0.05, |rho|>0.3) | %s |",
          if (is.na(result$n_spearman_sig)) "-" else as.character(result$n_spearman_sig)),
  sprintf("| Elapsed | %s min |",
          if (is.na(result$elapsed_min)) "-" else as.character(result$elapsed_min)),
  sprintf("| Status | %s |", result$status),
  "",
  "## Core panel 关键基因 (动态图)",
  "",
  "**NAM 谱系 TF**: Spi1 (PU.1), Runx1, Ikzf1, Pou2f2, Lyl1, Cebpb, Mafb",
  "",
  "**单核起点 marker** (用于 root): Ly6c2, Ccr2, Plac8",
  "",
  "**驻留巨噬终点**: Adgre1 (F4/80), Csf1r, Cd68, Mertk, Vsig4",
  "",
  "**NFκB / IFN**: Nfkb1, Rela, Relb, Irf1, Irf7, Stat1",
  "",
  "**炎症因子**: Il6, Il1b, Tnf, Cxcl9, Cxcl10, Ccl2, Ccl5",
  "",
  "**TGFb 通路** (验证 NAM 不是 TGFb 源): Tgfb1, Smad3, Smad7",
  "",
  "**Efferocytosis**: Elmo1, Dock4, Lrp1, Mertk, Gas6, Axl",
  "",
  "**MHC-I**: B2m, H2-D1, H2-K1, H2-Aa, H2-Ab1",
  "",
  "**M1/M2 极化**: Nos2, Arg1, Mrc1, Cd163, Il10",
  "",
  "## 关键产出",
  "",
  "- `Macrophage_NAM/trajectory.pdf` — UMAP + slingshot 曲线",
  "- `Macrophage_NAM/pseudotime_by_group.pdf` — **关键** tumor/blank/control 在同一轨迹上的分布",
  "- `Macrophage_NAM/dynamic_core_genes.pdf` — **关键** 核心基因沿伪时间分组曲线",
  "- `Macrophage_NAM/curve_assignment.pdf` — 多曲线分配",
  "- `Macrophage_NAM/heatmap_top_genes.pdf` — top 80 Spearman 基因热图",
  "- `Macrophage_NAM/spearman_correlation.csv` — 全基因 Spearman 结果",
  "- `Macrophage_NAM/pseudotime_values.csv` — 每个细胞的 pseudotime 值"
)

writeLines(summary_lines, file.path(OUT_DIR, "summary", "pseudotime_nam_summary.md"))

cat("\n========================================\n")
cat("Pseudotime NAM analysis complete\n")
cat(sprintf("  Output dir : %s\n", OUT_DIR))
cat(sprintf("  Status     : %s\n", result$status))
cat("========================================\n")
