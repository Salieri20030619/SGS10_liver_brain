#!/usr/bin/env Rscript
# =============================================================================
# SGS10 Task B - HE Overlay v2 (fix Y-flip + use cell_locations for extent)
# Date: 2026-04-20
# Author: Claude
#
# v2 vs v1 fixes:
#   1. 读 cell_locations.tsv.gz 获取真正的坐标 extent (相当于官方 size_x/size_y)
#   2. 修复 Y 轴翻转: RCTD 里 y 被 flip 过了 (y = y_max - y_original)
#      必须把 RCTD.y 转回 HE 坐标系: y_correct = y_max - RCTD.y
#   3. Brain_3/Brain_8 暂不做 (HE 路径还没确认, 先跑 Liver_3 + nao56 验证对齐)
# =============================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(FNN)
  library(png)
  library(grid)
  library(RColorBrewer)
})

RCTD_DIR  <- "~/SGS10_analysis/spatial_results/rctd/per_sample"
OUT_ROOT  <- "~/SGS10_analysis/taskB_he_overlay_v2"
DPI       <- 300
K_CLUSTER <- 10
KNN_SMOOTH <- 15
SET_SEED  <- 42
TOP_N_CT  <- 10
WEIGHT_THR <- 0.1

dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)
LOG_FILE <- file.path(OUT_ROOT, "run.log")
writeLines(c(paste("# HE overlay v2 run -", Sys.time()), ""), LOG_FILE)

log_msg <- function(msg) {
  line <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", msg)
  cat(line, "\n")
  cat(line, "\n", file = LOG_FILE, append = TRUE, sep = "")
}

# ---------- 样本 -> HE 路径 + cell_locations 路径 ----------
SAMPLE_PATHS <- list(
  "nao56" = list(
    he_timg = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/SGS24070263-nao56/data/matrix/24070263_nao_56_expression/nao_56_aligned_HE_TIMG.png",
    cell_locations = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/nao_56_filtered_feature_bc_matrix/cell_locations.tsv.gz"
  ),
  "Liver_3" = list(
    he_timg = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/liver_samples/25090559_Liver_3_expression/Liver_3_aligned_HE_TIMG.png",
    cell_locations = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/liver_samples/25090559_Liver_3_expression/Liver_3_filtered_feature_bc_matrix/cell_locations.tsv.gz"
  ),
  "Liver_8" = list(
    he_timg = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/liver_samples/25090559_Liver_8_expression/Liver_8_aligned_HE_TIMG.png",
    cell_locations = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/liver_samples/25090559_Liver_8_expression/Liver_8_filtered_feature_bc_matrix/cell_locations.tsv.gz"
  ),
  "Con_liver_1" = list(
    he_timg = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/liver_samples/25090559_Con_liver_1_expression/Con_liver_1_aligned_HE_TIMG.png",
    cell_locations = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/liver_samples/25090559_Con_liver_1_expression/Con_liver_1_filtered_feature_bc_matrix/cell_locations.tsv.gz"
  ),
  "Bla_liver_1" = list(
    he_timg = NA_character_,   # find 结果里 bla_liver 目录没 HE, 需确认
    cell_locations = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/bla_liver/Bla_liver_1_filtered_feature_bc_matrix/cell_locations.tsv.gz"
  ),
  "gan56" = list(
    he_timg = NA_character_,   # find 结果里 gan56 没 HE, 需确认
    cell_locations = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/gan56/gan_56_filtered_feature_bc_matrix/cell_locations.tsv.gz"
  )
)

GROUPS <- c(
  "Brain_3" = "tumor", "Liver_3" = "tumor", "nao56" = "tumor",
  "Liver_8" = "tumor", "gan56" = "tumor",
  "Con_liver_1" = "baseline", "Bla_liver_1" = "baseline",
  "Brain_8" = "tumor", "Con_brain_1" = "baseline", "Bla_brain_1" = "baseline"
)

# ---------- 读 cell_locations, 返回 extent + barcode->coord 映射 ----------
load_cell_locations <- function(path) {
  df <- read.table(gzfile(path), header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE)
  # 列: Cell_Barcode, X, Y  (按 head 结果)
  stopifnot(all(c("Cell_Barcode", "X", "Y") %in% colnames(df)) ||
            ncol(df) >= 3)
  if (!"Cell_Barcode" %in% colnames(df)) {
    colnames(df)[1:3] <- c("Cell_Barcode", "X", "Y")
  }
  list(
    df = df,
    x_min = min(df$X), x_max = max(df$X),
    y_min = min(df$Y), y_max = max(df$Y)
  )
}

extract_rctd <- function(rctd_file) {
  rctd <- readRDS(rctd_file)
  w <- NULL
  if (!is.null(rctd@results$weights)) {
    w <- as.matrix(rctd@results$weights)
  } else if (!is.null(rctd@results$weights_doublet)) {
    w <- if (!is.null(rctd@results$norm_weights))
           as.matrix(rctd@results$norm_weights)
         else as.matrix(rctd@results$weights_doublet)
  }
  if (is.null(w)) stop("No weights")
  coords <- rctd@spatialRNA@coords
  common <- intersect(rownames(w), rownames(coords))
  w <- w[common, , drop = FALSE]
  coords <- coords[common, , drop = FALSE]
  rs <- rowSums(w); rs[rs == 0] <- 1
  w <- w / rs
  list(weights = w, coords = coords,
       n_spots = nrow(w), celltypes = colnames(w))
}

smooth_labels <- function(coords, labels, k = 15, iters = 2) {
  xy <- as.matrix(coords[, 1:2])
  if (!is.numeric(xy[1, 1])) xy <- apply(xy, 2, as.numeric)
  nn <- FNN::get.knn(xy, k = k)$nn.index
  current <- labels
  for (it in seq_len(iters)) {
    new_l <- if (is.character(current)) character(length(current))
             else integer(length(current))
    for (i in seq_along(current)) {
      tab <- table(c(current[i], current[nn[i, ]]))
      best <- names(tab)[which.max(tab)]
      new_l[i] <- if (is.character(current)) best else as.integer(best)
    }
    current <- new_l
  }
  current
}

name_domains <- function(w, cluster, k) {
  out <- character(k)
  for (i in seq_len(k)) {
    idx <- which(cluster == i)
    if (length(idx) == 0) { out[i] <- sprintf("D%02d [empty]", i); next }
    mw <- colMeans(w[idx, , drop = FALSE])
    top2 <- sort(mw, decreasing = TRUE)[1:2]
    out[i] <- sprintf("D%02d [%s+%s]", i, names(top2)[1], names(top2)[2])
  }
  out
}

load_he <- function(he_path) {
  if (is.na(he_path) || !file.exists(he_path)) {
    return(NULL)
  }
  img <- png::readPNG(he_path)
  grob <- grid::rasterGrob(img, interpolate = TRUE,
                           width = grid::unit(1, "npc"),
                           height = grid::unit(1, "npc"))
  list(grob = grob, width = dim(img)[2], height = dim(img)[1])
}

# ---------- 主处理函数 ----------
process_sample <- function(sample_name) {
  t0 <- Sys.time()
  log_msg(sprintf("=== %s ===", sample_name))

  paths <- SAMPLE_PATHS[[sample_name]]
  if (is.null(paths)) {
    log_msg(sprintf("  [SKIP] no paths for %s", sample_name))
    return(FALSE)
  }

  # 1. 读 cell_locations 拿 extent
  log_msg("  reading cell_locations...")
  cl <- tryCatch(load_cell_locations(paths$cell_locations),
                 error = function(e) { log_msg(sprintf("  [ERR] %s", e$message)); NULL })
  if (is.null(cl)) return(FALSE)
  log_msg(sprintf("  cell_locations extent: x[%.0f,%.0f] y[%.0f,%.0f] n=%d",
                  cl$x_min, cl$x_max, cl$y_min, cl$y_max, nrow(cl$df)))

  # 2. 读 RCTD
  rctd_file <- path.expand(file.path(RCTD_DIR,
                                     sprintf("%s_rctd_results.rds", sample_name)))
  if (!file.exists(rctd_file)) {
    log_msg(sprintf("  [SKIP] RCTD not found"))
    return(FALSE)
  }
  log_msg("  reading RCTD...")
  d <- extract_rctd(rctd_file)
  log_msg(sprintf("  %d RCTD spots × %d CT", d$n_spots, length(d$celltypes)))

  # 3. ★ 关键: 重新对齐 RCTD 坐标到 HE 坐标系 ★
  # 通过匹配 barcode, 用 cell_locations 的真实 XY 替换 RCTD 的 coords
  log_msg("  re-aligning RCTD coords via barcode match to cell_locations...")
  bc_rctd <- rownames(d$coords)
  bc_cl <- cl$df$Cell_Barcode
  common_bc <- intersect(bc_rctd, bc_cl)
  log_msg(sprintf("    matched barcodes: %d / %d RCTD, %d cell_locations",
                  length(common_bc), length(bc_rctd), length(bc_cl)))

  if (length(common_bc) < 0.5 * length(bc_rctd)) {
    log_msg("  [WARN] <50% barcode overlap, abort - barcode mismatch?")
    log_msg(sprintf("    RCTD bc sample: %s", bc_rctd[1]))
    log_msg(sprintf("    cell_loc bc sample: %s", bc_cl[1]))
    return(FALSE)
  }

  # 对齐: 按 common_bc 顺序, 同时索引 d$weights 和 cl$df
  idx_rctd <- match(common_bc, bc_rctd)
  idx_cl   <- match(common_bc, bc_cl)

  weights_aligned <- d$weights[idx_rctd, , drop = FALSE]
  x_aligned <- cl$df$X[idx_cl]
  y_aligned <- cl$df$Y[idx_cl]

  log_msg(sprintf("  aligned: %d spots", length(common_bc)))
  log_msg(sprintf("  aligned coord: x[%.0f,%.0f] y[%.0f,%.0f]",
                  min(x_aligned), max(x_aligned),
                  min(y_aligned), max(y_aligned)))

  # 4. 读 HE
  he <- load_he(paths$he_timg)
  if (is.null(he)) {
    log_msg("  [SKIP] HE image not available")
    return(FALSE)
  }
  log_msg(sprintf("  HE image: %d x %d px", he$width, he$height))

  # 5. annotation extent: 用 cell_locations 的 xmax/ymax (SeekSpace 官方 size_x/y)
  #    HE 图会被 annotation_custom 自动拉伸到这个范围
  xmin_ann <- 0
  xmax_ann <- cl$x_max * 1.02
  ymin_ann <- 0
  ymax_ann <- cl$y_max * 1.02

  # 6. T4 domain (k-means + KNN on aligned coords)
  log_msg("  k-means + KNN smoothing...")
  set.seed(SET_SEED)
  km <- kmeans(weights_aligned, centers = K_CLUSTER, iter.max = 50, nstart = 10)
  coords_aligned <- data.frame(x = x_aligned, y = y_aligned)
  smoothed <- smooth_labels(coords_aligned, km$cluster,
                            k = KNN_SMOOTH, iters = 2)
  dom_names <- name_domains(weights_aligned, smoothed, K_CLUSTER)

  # 7. top CT
  log_msg("  top-CT argmax...")
  ct_sums <- colSums(weights_aligned)
  top_cts <- names(sort(ct_sums, decreasing = TRUE))[1:min(TOP_N_CT, length(ct_sums))]
  w_top <- weights_aligned[, top_cts, drop = FALSE]
  argmax_idx <- max.col(w_top, ties.method = "first")
  argmax_w <- w_top[cbind(seq_len(nrow(w_top)), argmax_idx)]
  top_ct_label <- ifelse(argmax_w >= WEIGHT_THR, top_cts[argmax_idx], "Other")

  # 8. 绘图 df
  df_dom <- data.frame(spatial_1 = x_aligned, spatial_2 = y_aligned,
                       Domain = factor(dom_names[smoothed], levels = dom_names))
  df_ct  <- data.frame(spatial_1 = x_aligned, spatial_2 = y_aligned,
                       TopCT = factor(top_ct_label, levels = c(top_cts, "Other")))

  out_dir <- file.path(OUT_ROOT, sample_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # ===== 图 1: HE + T4 domain =====
  log_msg("  plot 1: HE + domain")
  pal10 <- c(brewer.pal(9, "Set1"), "#999999")

  p_dom <- ggplot() +
    annotation_custom(grob = he$grob,
                      xmin = xmin_ann, xmax = xmax_ann,
                      ymin = ymin_ann, ymax = ymax_ann) +
    geom_point(data = df_dom,
               aes(x = spatial_1, y = spatial_2, color = Domain),
               size = 0.8, alpha = 0.75) +
    scale_color_manual(values = pal10) +
    scale_x_continuous(limits = c(xmin_ann, xmax_ann), expand = c(0, 0)) +
    scale_y_continuous(limits = c(ymin_ann, ymax_ann), expand = c(0, 0)) +
    coord_fixed() +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1),
                                ncol = 1)) +
    labs(title = sprintf("HE + T4 Spatial Domain - %s (%s)",
                         sample_name, GROUPS[sample_name]),
         subtitle = sprintf("%d spots aligned, k-means k=%d + KNN k=%d",
                            length(common_bc), K_CLUSTER, KNN_SMOOTH),
         x = NULL, y = NULL, color = "Niche") +
    theme_void(base_size = 11) +
    theme(legend.text = element_text(size = 7),
          legend.key.size = unit(0.4, "cm"),
          plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 9, color = "grey30"),
          legend.position = "right",
          plot.background = element_rect(fill = "white", color = NA))

  ggsave(file.path(out_dir, "01_HE_overlay_domain.png"),
         p_dom, width = 12, height = 6, dpi = DPI, bg = "white")

  # ===== 图 2: HE + top CT =====
  log_msg("  plot 2: HE + top CT")
  pal_ct <- c(brewer.pal(9, "Set1"), brewer.pal(7, "Dark2"))[1:TOP_N_CT]
  pal_ct <- setNames(c(pal_ct, "grey80"), c(top_cts, "Other"))

  p_ct <- ggplot() +
    annotation_custom(grob = he$grob,
                      xmin = xmin_ann, xmax = xmax_ann,
                      ymin = ymin_ann, ymax = ymax_ann) +
    geom_point(data = df_ct %>% dplyr::filter(TopCT == "Other"),
               aes(x = spatial_1, y = spatial_2),
               color = "grey75", size = 0.5, alpha = 0.4) +
    geom_point(data = df_ct %>% dplyr::filter(TopCT != "Other"),
               aes(x = spatial_1, y = spatial_2, color = TopCT),
               size = 0.9, alpha = 0.85) +
    scale_color_manual(values = pal_ct, drop = FALSE) +
    scale_x_continuous(limits = c(xmin_ann, xmax_ann), expand = c(0, 0)) +
    scale_y_continuous(limits = c(ymin_ann, ymax_ann), expand = c(0, 0)) +
    coord_fixed() +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1),
                                ncol = 1)) +
    labs(title = sprintf("HE + Top CT (argmax) - %s (%s)",
                         sample_name, GROUPS[sample_name]),
         subtitle = sprintf("top-%d CT, weight >= %.2f (else grey)",
                            TOP_N_CT, WEIGHT_THR),
         x = NULL, y = NULL, color = "Top cell type") +
    theme_void(base_size = 11) +
    theme(legend.text = element_text(size = 7),
          legend.key.size = unit(0.4, "cm"),
          plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 9, color = "grey30"),
          legend.position = "right",
          plot.background = element_rect(fill = "white", color = NA))

  ggsave(file.path(out_dir, "02_HE_overlay_topCT.png"),
         p_ct, width = 12, height = 6, dpi = DPI, bg = "white")

  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  log_msg(sprintf("  done in %.1f s", dt))
  TRUE
}

# ---------- 调度 ----------
log_msg("HE overlay v2: Liver_3 + nao56 (test)")
log_msg("")

for (s in c("Liver_3", "nao56")) {
  process_sample(s)
}

log_msg("")
log_msg("=== DONE ===")
system(sprintf("find %s -name '*.png' | xargs ls -lh", path.expand(OUT_ROOT)))
