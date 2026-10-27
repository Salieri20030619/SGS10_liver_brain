#!/usr/bin/env Rscript
# =============================================================================
# SGS10 Task B - BayesSpace Test (Single Sample, Brain_3)
# Date: 2026-04-20
# Author: Claude
#
# 目的: 验证 BayesSpace 在 SeekSpace 稀疏数据 + RCTD doublet weights 上
#       是否能复现客户参考图 1 (彩色大块分层) 的视觉效果
#
# 策略:
#   1. 用 RCTD weights (spot × celltype) 当作 "低维 expression"
#   2. 走 BayesSpace 的 spatial clustering (q=10, platform=ST)
#   3. 画 spatial scatter 比较 v2 (k-means + KNN) vs v4 (BayesSpace)
#
# Input:
#   ~/SGS10_analysis/spatial_results/rctd/per_sample/Brain_3_rctd_results.rds
#
# Output:
#   ~/SGS10_analysis/taskB_bayesspace_test/
#     Brain_3_bayesspace_domain.png        <- 主图
#     Brain_3_kmeans_vs_bayesspace.png     <- 对照图
#     Brain_3_bayesspace_labels.csv        <- spot → domain 标签
#     run.log
#
# Usage (server):
#   mkdir -p ~/SGS10_analysis/taskB_bayesspace_test
#   # scp this script there
#   cd ~/SGS10_analysis/
#   # 检查依赖 (第一次跑前做)
#   Rscript -e 'if(!requireNamespace("BayesSpace")) BiocManager::install("BayesSpace")'
#   # 跑脚本
#   screen -S bayesspace_test
#   Rscript taskB_bayesspace_test/test_bayesspace_brain3.R 2>&1 | tee run.log
#   # Ctrl+A D
# =============================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(SingleCellExperiment)
  library(BayesSpace)
  library(RColorBrewer)
  library(scales)
})

# ---------- 配置 ----------
SAMPLE      <- "Brain_3"
RCTD_FILE   <- path.expand(sprintf(
  "~/SGS10_analysis/spatial_results/rctd/per_sample/%s_rctd_results.rds", SAMPLE))
OUT_DIR     <- path.expand("~/SGS10_analysis/taskB_bayesspace_test")
DPI         <- 300
Q           <- 10       # 目标 domain 数
NREP        <- 10000    # BayesSpace MCMC 迭代数 (推荐 10000+)
SET_SEED    <- 42

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
LOG_FILE <- file.path(OUT_DIR, "run.log")
writeLines(c(paste("# BayesSpace test run log -", Sys.time()), ""),
           LOG_FILE)

log_msg <- function(msg) {
  line <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", msg)
  cat(line, "\n")
  cat(line, "\n", file = LOG_FILE, append = TRUE, sep = "")
}

# ---------- 读 RCTD ----------
log_msg(sprintf("Reading RCTD: %s", basename(RCTD_FILE)))
rctd <- readRDS(RCTD_FILE)

w <- NULL
if (!is.null(rctd@results$weights)) {
  w <- as.matrix(rctd@results$weights)
} else if (!is.null(rctd@results$weights_doublet)) {
  if (!is.null(rctd@results$norm_weights)) {
    w <- as.matrix(rctd@results$norm_weights)
  } else {
    w <- as.matrix(rctd@results$weights_doublet)
  }
}
if (is.null(w)) stop("No weights found in RCTD")

coords <- rctd@spatialRNA@coords
common <- intersect(rownames(w), rownames(coords))
w <- w[common, , drop = FALSE]
coords <- coords[common, , drop = FALSE]

# Normalize rows
row_sums <- rowSums(w)
row_sums[row_sums == 0] <- 1
w <- w / row_sums

n_spots <- nrow(w)
n_ct <- ncol(w)
log_msg(sprintf("Loaded: %d spots × %d celltypes", n_spots, n_ct))

# ---------- 构造 SingleCellExperiment ----------
# BayesSpace 需要: counts (genes × cells) + colData (x, y, col, row)
# Trick: 把 RCTD weights 转置作为 "counts" (CT × spots)
# 这样 PCA 会在 CT 空间里做降维
log_msg("Constructing SingleCellExperiment for BayesSpace...")

# 坐标 - BayesSpace 要求整数 col/row
# 把 SeekSpace 浮点坐标离散化为格点 (col = round(x / step), row = round(y / step))
x <- as.numeric(coords[, 1])
y <- as.numeric(coords[, 2])

# 估算 step (相邻 spot 最小距离的 median)
# 避免太细, 用分位数的粗估
rg_x <- diff(range(x))
rg_y <- diff(range(y))
# SeekSpace spot 大约 ~10000 个, 组织尺寸 ~100 units -> step ~1
# 用相邻 spot 距离中位数作为步长
set.seed(SET_SEED)
sample_idx <- sample(seq_len(n_spots), min(2000, n_spots))
d_sub <- as.matrix(dist(cbind(x[sample_idx], y[sample_idx])))
diag(d_sub) <- Inf
step <- median(apply(d_sub, 1, min))
log_msg(sprintf("Estimated grid step = %.3f", step))

col_int <- as.integer(round((x - min(x)) / step)) + 1L
row_int <- as.integer(round((y - min(y)) / step)) + 1L

# 如果有 col/row 冲突 (多个 spot 映到同一格), 加扰动
dup_check <- duplicated(paste(col_int, row_int))
n_dup <- sum(dup_check)
if (n_dup > 0) {
  log_msg(sprintf("  [warn] %d duplicated col/row pairs, resolving...", n_dup))
  while (any(duplicated(paste(col_int, row_int)))) {
    dup_idx <- which(duplicated(paste(col_int, row_int)))
    # 微调行
    row_int[dup_idx] <- row_int[dup_idx] + sample(c(-1L, 1L), length(dup_idx),
                                                   replace = TRUE)
  }
}

col_data <- DataFrame(
  row    = row_int,
  col    = col_int,
  imagerow = y,
  imagecol = x
)
rownames(col_data) <- rownames(w)

# counts: rows=features (= cell types), cols=spots
# 放大 weights 到整数 counts scale (× 1000) 以满足 BayesSpace 内部检查
counts_mat <- t(w) * 1000
storage.mode(counts_mat) <- "integer"
counts_mat <- as(counts_mat, "CsparseMatrix")

sce <- SingleCellExperiment(
  assays = list(counts = counts_mat, logcounts = log1p(counts_mat)),
  colData = col_data
)
log_msg(sprintf("SCE constructed: %d × %d", nrow(sce), ncol(sce)))

# ---------- 跑 BayesSpace ----------
# spatialPreprocess 是标准前置 (PCA + log)
log_msg("BayesSpace spatialPreprocess (PCA)...")
sce <- spatialPreprocess(sce,
                         platform = "ST",
                         n.PCs = min(15, n_ct - 1),
                         n.HVGs = n_ct,   # 全用 (CT 就 22 个, 没必要筛)
                         log.normalize = FALSE,
                         skip.PCA = FALSE)

log_msg(sprintf("BayesSpace spatialCluster q=%d nrep=%d (5-15 min)...", Q, NREP))
t_cluster <- Sys.time()
set.seed(SET_SEED)
sce <- spatialCluster(sce,
                      q = Q,
                      platform = "ST",
                      d = min(15, n_ct - 1),
                      init.method = "mclust",
                      model = "t",
                      gamma = 2,          # spatial smoothing strength (ST 默认 2)
                      nrep = NREP,
                      burn.in = 1000,
                      save.chain = FALSE)
dt_cluster <- as.numeric(difftime(Sys.time(), t_cluster, units = "mins"))
log_msg(sprintf("spatialCluster done in %.1f min", dt_cluster))

# 提取 cluster labels
bs_labels <- colData(sce)$spatial.cluster
log_msg(sprintf("BayesSpace cluster distribution:"))
print(table(bs_labels))

# ---------- 自动命名每 domain ----------
# 每 domain 的 top 2 CT (按 mean weight)
name_domains_bs <- function(w, labels, k) {
  out <- character(k)
  for (i in seq_len(k)) {
    idx <- which(labels == i)
    if (length(idx) == 0) {
      out[i] <- sprintf("D%02d [empty]", i)
      next
    }
    mean_w <- colMeans(w[idx, , drop = FALSE])
    top2 <- sort(mean_w, decreasing = TRUE)[1:2]
    out[i] <- sprintf("D%02d [%s+%s]", i,
                      names(top2)[1], names(top2)[2])
  }
  out
}

domain_names <- name_domains_bs(w, bs_labels, Q)
log_msg("Domain names:")
for (dn in domain_names) log_msg(sprintf("  %s", dn))

# ---------- 画主图: BayesSpace spatial domain ----------
log_msg("Plotting BayesSpace spatial domain...")

pal10 <- c(brewer.pal(9, "Set1"), "#666666")

df_plot <- data.frame(x = x, y = y,
                      Domain = factor(domain_names[bs_labels],
                                      levels = domain_names))

p_bs <- ggplot(df_plot, aes(x = x, y = y, color = Domain)) +
  geom_point(size = 1.2, alpha = 0.85) +
  scale_color_manual(values = pal10) +
  scale_y_reverse() + coord_fixed() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1),
                              ncol = 1)) +
  labs(title = sprintf("BayesSpace Spatial Domain - %s (tumor)", SAMPLE),
       subtitle = sprintf("%d spots, q=%d, nrep=%d, gamma=2 (ST platform)",
                          n_spots, Q, NREP),
       x = "Spatial x", y = "Spatial y", color = "Domain") +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank(),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.4, "cm"),
        plot.title = element_text(face = "bold"),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        aspect.ratio = 1)

ggsave(file.path(OUT_DIR, sprintf("%s_bayesspace_domain.png", SAMPLE)),
       p_bs, width = 9, height = 7, dpi = DPI, bg = "white")

# ---------- 保存 spot 级标签 CSV ----------
out_csv <- file.path(OUT_DIR,
                     sprintf("%s_bayesspace_labels.csv", SAMPLE))
write.csv(data.frame(barcode = rownames(w),
                     x = x, y = y,
                     col = col_int, row = row_int,
                     bayesspace_cluster = bs_labels,
                     domain = domain_names[bs_labels]),
          out_csv, row.names = FALSE)
log_msg(sprintf("Saved labels: %s", basename(out_csv)))

# ---------- 对比图: v2 k-means+KNN vs BayesSpace ----------
# 顺便跑一遍 v2 的 k-means + KNN, 放一起对比
log_msg("Running v2 k-means+KNN for comparison...")

set.seed(SET_SEED)
km <- kmeans(w, centers = Q, iter.max = 50, nstart = 10)

# KNN k=15 smoothing
if (requireNamespace("FNN", quietly = TRUE)) {
  library(FNN)
  xy <- cbind(x, y)
  nn <- FNN::get.knn(xy, k = 15)$nn.index
  smoothed <- km$cluster
  for (it in 1:2) {
    new_l <- integer(length(smoothed))
    for (i in seq_along(smoothed)) {
      tab <- table(c(smoothed[i], smoothed[nn[i, ]]))
      new_l[i] <- as.integer(names(tab)[which.max(tab)])
    }
    smoothed <- new_l
  }
  km_names <- name_domains_bs(w, smoothed, Q)

  df_km <- data.frame(x = x, y = y,
                      Domain = factor(km_names[smoothed], levels = km_names))

  p_km <- ggplot(df_km, aes(x = x, y = y, color = Domain)) +
    geom_point(size = 1.2, alpha = 0.85) +
    scale_color_manual(values = pal10) +
    scale_y_reverse() + coord_fixed() +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1),
                                ncol = 1)) +
    labs(title = sprintf("v2 k-means + KNN k=15 - %s", SAMPLE),
         subtitle = sprintf("%d spots, k=%d (no spatial prior, just smoothing)",
                           n_spots, Q),
         x = "Spatial x", y = "Spatial y", color = "Domain") +
    theme_bw(base_size = 11) +
    theme(panel.grid = element_blank(),
          legend.text = element_text(size = 7),
          legend.key.size = unit(0.4, "cm"),
          plot.title = element_text(face = "bold"),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          aspect.ratio = 1)

  # 拼图
  if (requireNamespace("patchwork", quietly = TRUE)) {
    library(patchwork)
    p_combined <- p_km | p_bs
    ggsave(file.path(OUT_DIR,
                     sprintf("%s_kmeans_vs_bayesspace.png", SAMPLE)),
           p_combined, width = 18, height = 7, dpi = DPI, bg = "white")
    log_msg("Saved comparison plot.")
  }
}

log_msg("")
log_msg("=== DONE ===")
log_msg(sprintf("Check outputs in: %s", OUT_DIR))
system(sprintf("ls -lh %s", OUT_DIR))
