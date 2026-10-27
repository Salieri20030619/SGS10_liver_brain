#!/usr/bin/env Rscript
# ============================================================
# fix_he_overlay_BC_v2.R
# 仅重跑 Panel B + Panel C，修复:
#   B: domain 配色在 H&E 上不可区分 → 加粗点+提高 alpha+确认 palette
#   C: 图例标题过长 → 缩短; 留白过大 → 紧凑裁切
#
# 数据读取方式与原 he_overlay_RCTD100.R 完全一致:
#   坐标: rctd@spatialRNA@coords
#   权重: rctd@results$weights (归一化)
#   H&E:  official RDS → seu@misc$info[[1]]$size_x/y + aligned PNG
#   域:   cosine similarity assign (原脚本逻辑复制)
#
# 用法: Rscript scripts/fix_he_overlay_BC_v2.R
#       Rscript scripts/fix_he_overlay_BC_v2.R --sample Liver_3
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(png)
  library(grid)
  library(dplyr)
  library(RANN)
  library(optparse)
})

option_list <- list(
  make_option("--sample", type = "character", default = "all",
              help = "sample name or 'all'"),
  make_option("--out_dir", type = "character",
              default = path.expand("~/SGS10_analysis/taskB_he_overlay_RCTD100"))
)
opt <- parse_args(OptionParser(option_list = option_list))

# ======================== Config ========================
PR       <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
OUT_DIR  <- path.expand(opt$out_dir)
RCTD_DIR <- file.path(PR, "spatial_results/rctd_full/per_sample")
SD_DIR   <- file.path(PR, "taskB/spatial_domain_RCTD100")

KNN_K       <- 15
POINT_SIZE  <- 0.8     # 比原脚本 0.6 稍大, 让颜色更醒目
POINT_ALPHA <- 0.85
PLOT_WIDTH  <- 16

log_msg <- function(msg) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))
}

# ======================== Paths ========================
BD <- "/mnt/BioData/bio/SGS10例"
MATRIX_LIVER <- file.path(BD, "SGS25090559-8例样本数据/data/matrix")
MATRIX_NAO56 <- file.path(BD, "SGS24070263-nao56/data/matrix/24070263_nao_56_expression")
MATRIX_GAN56 <- file.path(BD, "SGS24070264-gan56/data/matrix/24070264_gan_56_expression")
RDS_LIVER    <- file.path(BD, "SGS25090559-8例样本数据/data/rds")
RDS_NAO56    <- file.path(BD, "SGS24070263-nao56/data/rds/24070263_nao_56_expression")
RDS_GAN56    <- file.path(BD, "SGS24070264-gan56/data/rds/24070264_gan_56_expression")

SAMPLES <- list(
  nao56       = list(tissue = "brain", group = "tumor",
                     official = file.path(RDS_NAO56, "nao_56.rds"),
                     he = file.path(MATRIX_NAO56, "nao_56_aligned_HE_TIMG.png")),
  gan56       = list(tissue = "liver", group = "tumor",
                     official = file.path(RDS_GAN56, "gan_56.rds"),
                     he = file.path(MATRIX_GAN56, "gan_56_aligned_HE_TIMG.png")),
  Brain_3     = list(tissue = "brain", group = "tumor",
                     official = file.path(RDS_LIVER, "25090559_Brain_3_expression/Brain_3.rds"),
                     he = file.path(MATRIX_LIVER, "25090559_Brain_3_expression/Brain_3_aligned_HE_TIMG.png")),
  Brain_8     = list(tissue = "brain", group = "tumor",
                     official = file.path(RDS_LIVER, "25090559_Brain_8_expression/Brain_8.rds"),
                     he = file.path(MATRIX_LIVER, "25090559_Brain_8_expression/Brain_8_aligned_HE_TIMG.png")),
  Bla_brain_1 = list(tissue = "brain", group = "baseline",
                     official = file.path(RDS_LIVER, "25090559_Bla_brain_1_expression/Bla_brain_1.rds"),
                     he = file.path(MATRIX_LIVER, "25090559_Bla_brain_1_expression/Bla_brain_1_aligned_HE_TIMG.png")),
  Con_brain_1 = list(tissue = "brain", group = "control",
                     official = file.path(RDS_LIVER, "25090559_Con_brain_1_expression/Con_brain_1.rds"),
                     he = file.path(MATRIX_LIVER, "25090559_Con_brain_1_expression/Con_brain_1_aligned_HE_TIMG.png")),
  Liver_3     = list(tissue = "liver", group = "tumor",
                     official = file.path(RDS_LIVER, "25090559_Liver_3_expression/Liver_3.rds"),
                     he = file.path(MATRIX_LIVER, "25090559_Liver_3_expression/Liver_3_aligned_HE_TIMG.png")),
  Liver_8     = list(tissue = "liver", group = "tumor",
                     official = file.path(RDS_LIVER, "25090559_Liver_8_expression/Liver_8.rds"),
                     he = file.path(MATRIX_LIVER, "25090559_Liver_8_expression/Liver_8_aligned_HE_TIMG.png")),
  Bla_liver_1 = list(tissue = "liver", group = "baseline",
                     official = file.path(RDS_LIVER, "25090559_Bla_liver_1_expression/Bla_liver_1.rds"),
                     he = file.path(MATRIX_LIVER, "25090559_Bla_liver_1_expression/Bla_liver_1_aligned_HE_TIMG.png")),
  Con_liver_1 = list(tissue = "liver", group = "control",
                     official = file.path(RDS_LIVER, "25090559_Con_liver_1_expression/Con_liver_1.rds"),
                     he = file.path(MATRIX_LIVER, "25090559_Con_liver_1_expression/Con_liver_1_aligned_HE_TIMG.png"))
)

if (opt$sample != "all") {
  if (!opt$sample %in% names(SAMPLES)) stop("Unknown sample: ", opt$sample)
  SAMPLES <- SAMPLES[opt$sample]
}

# ======================== Load domain profiles + naming ========================
load_domain_profile <- function(f) {
  d <- read.csv(f, check.names = FALSE)
  rownames(d) <- d[[1]]; d[[1]] <- NULL
  as.matrix(d)
}
load_naming <- function(f) {
  d <- read.csv(f, stringsAsFactors = FALSE)
  setNames(d$auto_name, d$domain)
}

prof_brain <- load_domain_profile(file.path(SD_DIR, "brain_domain_ct_profile.csv"))
prof_liver <- load_domain_profile(file.path(SD_DIR, "liver_domain_ct_profile.csv"))
name_brain <- load_naming(file.path(SD_DIR, "brain_domain_naming.csv"))
name_liver <- load_naming(file.path(SD_DIR, "liver_domain_naming.csv"))

log_msg(sprintf("Profiles loaded: brain %dx%d, liver %dx%d",
                nrow(prof_brain), ncol(prof_brain),
                nrow(prof_liver), ncol(prof_liver)))

# ======================== Color palettes ========================
# FIX B: 高饱和度 10 色 palette, 在 H&E 粉紫底图上能清晰区分
# 原脚本的 DOMAIN_COLORS 配色理论上没问题, 但实际渲染时
# 可能因为 alpha + 点太小 + H&E 底色 导致看不清
# → 加大 POINT_SIZE, 提高 POINT_ALPHA, 并微调配色使其更鲜艳
DOMAIN_COLORS <- c(
  "D01" = "#E31A1C",   # 鲜红
  "D02" = "#1F78B4",   # 深蓝
  "D03" = "#33A02C",   # 深绿
  "D04" = "#6A3D9A",   # 深紫
  "D05" = "#FF7F00",   # 鲜橙
  "D06" = "#FFD700",   # 金黄
  "D07" = "#8B4513",   # 棕
  "D08" = "#E91E90",   # 亮粉 (原 #f781bf 太淡)
  "D09" = "#00CED1",   # 青绿 (原 #66c2a5 太淡)
  "D10" = "#708090"    # 石板灰 (原 #bebada 与粉紫色 H&E 太接近)
)

# ======================== Helpers (与原脚本一致) ========================
smooth_weights <- function(coords, weights, k = KNN_K) {
  nn <- RANN::nn2(coords[, c("x", "y")], coords[, c("x", "y")], k = k + 1)
  nn_idx <- nn$nn.idx
  out <- matrix(0, nrow = nrow(weights), ncol = ncol(weights),
                dimnames = dimnames(weights))
  for (j in seq_len(ncol(weights))) {
    out[, j] <- rowMeans(matrix(weights[nn_idx, j], nrow = nrow(weights)))
  }
  rs <- rowSums(out)
  out / pmax(rs, 1e-9)
}

assign_domain <- function(weights_smoothed, profile_matrix) {
  common_ct <- intersect(colnames(weights_smoothed), colnames(profile_matrix))
  W <- weights_smoothed[, common_ct, drop = FALSE]
  P <- profile_matrix[, common_ct, drop = FALSE]
  W_norm <- sqrt(rowSums(W^2))
  P_norm <- sqrt(rowSums(P^2))
  cos_sim <- (W %*% t(P)) / (pmax(W_norm, 1e-9) %o% pmax(P_norm, 1e-9))
  rownames(profile_matrix)[max.col(cos_sim, ties.method = "first")]
}

# ======================== Main loop ========================
log_msg(sprintf("=== fix_he_overlay_BC_v2 start (sample=%s) ===", opt$sample))

for (s in names(SAMPLES)) {
  t0 <- Sys.time()
  cfg <- SAMPLES[[s]]
  sample_out <- file.path(OUT_DIR, s)
  dir.create(sample_out, showWarnings = FALSE, recursive = TRUE)
  log_msg(sprintf("---- [%s] %s/%s ----", s, cfg$tissue, cfg$group))

  # --- Select organ-specific config ---
  if (cfg$tissue == "brain") {
    dom_profile <- prof_brain
    dom_naming  <- name_brain
  } else {
    dom_profile <- prof_liver
    dom_naming  <- name_liver
  }

  # --- Load official RDS for H&E size ---
  if (!file.exists(cfg$official)) { log_msg("  SKIP: official missing"); next }
  seu <- readRDS(cfg$official)
  info <- seu@misc$info[[1]]
  size_x <- info$size_x; size_y <- info$size_y
  rm(seu); gc(verbose = FALSE)
  log_msg(sprintf("  HE size: %d x %d", size_x, size_y))

  # --- Load H&E image ---
  if (!file.exists(cfg$he)) { log_msg("  SKIP: HE missing"); next }
  he_img <- png::readPNG(cfg$he)
  he_grob <- rasterGrob(he_img, interpolate = FALSE,
                        width = unit(1, "npc"), height = unit(1, "npc"))
  rm(he_img); gc(verbose = FALSE)

  # --- Load RCTD results ---
  rctd_path <- file.path(RCTD_DIR, paste0(s, "_rctd_results.rds"))
  if (!file.exists(rctd_path)) { log_msg("  SKIP: RCTD missing"); next }
  rctd <- readRDS(rctd_path)
  coords  <- rctd@spatialRNA@coords
  weights <- as.matrix(rctd@results$weights)
  rs <- rowSums(weights); weights <- weights / pmax(rs, 1e-9)
  rm(rctd); gc(verbose = FALSE)

  common <- intersect(rownames(coords), rownames(weights))
  coords  <- coords[common, , drop = FALSE]
  weights <- weights[common, , drop = FALSE]
  log_msg(sprintf("  RCTD: %d spots, %d CT", nrow(coords), ncol(weights)))

  # --- Smooth + assign domain ---
  weights_smoothed <- smooth_weights(coords, weights, k = KNN_K)
  domain <- assign_domain(weights_smoothed, dom_profile)

  # --- Plot settings ---
  aspect      <- size_y / size_x
  plot_height <- PLOT_WIDTH * aspect * 1.05

  plot_theme <- theme_void(base_size = 12) +
    theme(plot.title    = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(color = "grey30", size = 10),
          legend.position = "right",
          legend.title = element_text(face = "bold", size = 11),
          legend.text  = element_text(size = 9),
          plot.margin  = margin(5, 5, 5, 5))

  # ===== Figure B: domain (FIXED) =====
  plot_df_B <- data.frame(x = coords$x, y = coords$y, domain = domain)
  all_doms <- sprintf("D%02d", 1:10)
  plot_df_B$domain <- factor(plot_df_B$domain, levels = all_doms)

  # 构建 legend labels: "D01 name" 格式
  # dom_naming 的 key 是 D_new_XX, 需要映射到 D01-D10
  # 但原脚本直接用 dom_naming[all_doms] — 说明 dom_naming 的 key 就是 D01-D10
  domain_legend <- dom_naming[all_doms]
  # 如果某些 NA (naming 文件里 key 不匹配), fallback 到 D01 等
  domain_legend[is.na(domain_legend)] <- all_doms[is.na(domain_legend)]

  color_map_B <- setNames(DOMAIN_COLORS[all_doms], all_doms)

  pB <- ggplot(plot_df_B, aes(x = x, y = y, color = domain)) +
    annotation_custom(he_grob, xmin = 0, xmax = size_x,
                      ymin = size_y, ymax = 0) +
    geom_point(size = POINT_SIZE, alpha = POINT_ALPHA) +
    scale_color_manual(values = color_map_B,
                       labels = domain_legend,
                       name   = "Spatial domain (10)",
                       drop   = FALSE) +
    coord_fixed(xlim = c(0, size_x), ylim = c(size_y, 0), expand = FALSE) +
    labs(title    = sprintf("HE + Spatial Domain (10) — %s (%s)", s, cfg$group),
         subtitle = sprintf("%d spots assigned via cosine similarity",
                            nrow(plot_df_B))) +
    plot_theme +
    guides(color = guide_legend(override.aes = list(size = 4, alpha = 1),
                                ncol = 1))

  ggsave(file.path(sample_out, "B_spatial_domain_on_HE.png"), pB,
         width = PLOT_WIDTH, height = plot_height, dpi = 150,
         limitsize = FALSE, bg = "white")
  log_msg("  B saved (fixed palette)")
  rm(pB, plot_df_B); gc(verbose = FALSE)

  # ===== Figure C: tumor-induced score (FIXED) =====
  if (cfg$tissue == "liver") {
    score_cts <- intersect(c("Macrophage_M1", "HSC_Activated_ECM", "cDC2",
                             "cDC1", "Macrophage_NAM", "Neutrophils"),
                           colnames(weights))
  } else {
    score_cts <- intersect(c("CPECs", "Ependymal", "VLMC",
                             "Microglia_Homeostatic"),
                           colnames(weights))
  }

  if (length(score_cts) > 0) {
    score <- rowSums(weights[, score_cts, drop = FALSE])
    lim_q <- quantile(score, c(0.05, 0.95), na.rm = TRUE)
    score_clip <- pmin(pmax(score, lim_q[1]), lim_q[2])

    plot_df_C <- data.frame(x = coords$x, y = coords$y, score = score_clip)

    # FIX C: 图例标题缩短, 不再列出所有 cell type 名
    legend_title_C <- "Tumor-induced\nniche score"

    pC <- ggplot(plot_df_C, aes(x = x, y = y, color = score)) +
      annotation_custom(he_grob, xmin = 0, xmax = size_x,
                        ymin = size_y, ymax = 0) +
      geom_point(size = POINT_SIZE * 1.2, alpha = 0.95) +
      scale_color_gradient2(
        low      = "#1b7837",
        mid      = "#f7f7f7",
        high     = "#b35806",
        midpoint = median(score, na.rm = TRUE),
        name     = legend_title_C    # ← FIXED: 短标题
      ) +
      coord_fixed(xlim = c(0, size_x), ylim = c(size_y, 0), expand = FALSE) +
      labs(
        title    = sprintf("HE + Tumor-induced niche score — %s (%s)",
                           s, cfg$group),
        # cell type 细节移到 subtitle
        subtitle = sprintf("Sum of %d cell type weights (%s)",
                           length(score_cts),
                           paste(score_cts, collapse = " + "))
      ) +
      plot_theme

    ggsave(file.path(sample_out, "C_tumor_induced_score_on_HE.png"), pC,
           width = PLOT_WIDTH, height = plot_height, dpi = 150,
           limitsize = FALSE, bg = "white")
    log_msg(sprintf("  C saved (fixed legend, score range [%.3f, %.3f])",
                    min(score, na.rm = TRUE), max(score, na.rm = TRUE)))
    rm(pC, plot_df_C, score, score_clip); gc(verbose = FALSE)
  } else {
    log_msg("  C SKIP: no score cell types")
  }

  rm(coords, weights, weights_smoothed, domain); gc(verbose = FALSE)
  log_msg(sprintf("  done in %.1f min", difftime(Sys.time(), t0, units = "mins")))
}

log_msg("=== fix_he_overlay_BC_v2 ALL DONE ===")
