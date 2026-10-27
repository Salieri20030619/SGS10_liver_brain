#!/usr/bin/env Rscript
# ============================================================
# fix_he_overlay_BC_v3.R
# 仅重跑 Panel B + C + D + E，修复:
#   B: domain 名 D_new_XX 被强制转 D01 → 全 NA → 灰色
#      → 用 D_new_XX 原名做 factor, 手动映射到高对比色
#   C: 图例标题过长 → 缩短
#   D/E: 脑端新增 M04 (Vascular_BBB) 圈标注
#
# 用法: Rscript scripts/fix_he_overlay_BC_v3.R
#       Rscript scripts/fix_he_overlay_BC_v3.R --sample Liver_3
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(png)
  library(grid)
  library(dplyr)
  library(RANN)
  library(optparse)
  library(hexbin)
})

option_list <- list(
  make_option("--sample", type = "character", default = "all"),
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
POINT_SIZE  <- 0.6
POINT_ALPHA <- 0.9
PLOT_WIDTH  <- 16
HEX_GRID_N  <- 200

log_msg <- function(msg) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))
}

# ======================== Sample paths ========================
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
load_microenv <- function(organ) {
  f <- file.path(SD_DIR, paste0(organ, "_microenv_naming.csv"))
  d <- read.csv(f, stringsAsFactors = FALSE)
  list(
    domain_to_microenv = setNames(d$microenv_id, d$domain),
    microenv_to_label  = setNames(d$microenv_label, d$microenv_id) |>
                          (\(x) x[!duplicated(names(x))])()
  )
}

prof_brain <- load_domain_profile(file.path(SD_DIR, "brain_domain_ct_profile.csv"))
prof_liver <- load_domain_profile(file.path(SD_DIR, "liver_domain_ct_profile.csv"))
name_brain <- load_naming(file.path(SD_DIR, "brain_domain_naming.csv"))
name_liver <- load_naming(file.path(SD_DIR, "liver_domain_naming.csv"))
me_brain   <- load_microenv("brain")
me_liver   <- load_microenv("liver")

log_msg(sprintf("Profiles: brain %dx%d, liver %dx%d",
                nrow(prof_brain), ncol(prof_brain),
                nrow(prof_liver), ncol(prof_liver)))
log_msg(sprintf("Microenv: brain %d, liver %d",
                length(me_brain$microenv_to_label),
                length(me_liver$microenv_to_label)))

# ======================== Color palettes ========================
# FIX B 核心: palette 用 D_new_XX 作 key (与 profile CSV 行名一致)
DOMAIN_COLORS_10 <- c(
  "#E31A1C", "#1F78B4", "#33A02C", "#6A3D9A", "#FF7F00",
  "#FFD700", "#8B4513", "#E91E90", "#00CED1", "#708090"
)

# 动态构建: 从 profile 行名读取实际 domain 名
build_domain_palette <- function(profile_matrix) {
  doms <- rownames(profile_matrix)
  n <- length(doms)
  cols <- DOMAIN_COLORS_10[seq_len(n)]
  setNames(cols, doms)
}

LIVER_MICROENV_COLORS <- c(
  M01 = "#2C7FB8", M02 = "#7FCDBB", M03 = "#41B6C4",
  M04 = "#A1D99B", M05 = "#E31A1C", M06 = "#FDAE6B",
  M07 = "#9E9AC8", M08 = "#cab2d6", M09 = "#fb9a99"
)
BRAIN_MICROENV_COLORS <- c(
  M01 = "#E31A1C", M02 = "#33A02C", M03 = "#A6761D",
  M04 = "#1F78B4", M05 = "#FF7F00", M06 = "#984EA3"
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

hex_bin_majority <- function(coords, label_vec, n_bins = HEX_GRID_N,
                             size_x, size_y) {
  hb <- hexbin::hexbin(coords$x, coords$y, xbins = n_bins,
                       xbnds = c(0, size_x), ybnds = c(0, size_y), IDs = TRUE)
  spot_hex_id <- hb@cID
  bin_df <- data.frame(hex_id = spot_hex_id, label = label_vec, stringsAsFactors = FALSE)
  bin_majority <- bin_df %>%
    group_by(hex_id) %>%
    summarise(n_spots = n(),
              label = names(sort(table(label), decreasing = TRUE))[1],
              .groups = "drop") %>%
    filter(!is.na(label))
  hex_centers <- hexbin::hcell2xy(hb)
  hex_centers_df <- data.frame(hex_id = hb@cell, x_c = hex_centers$x, y_c = hex_centers$y)
  bin_majority <- bin_majority %>% left_join(hex_centers_df, by = "hex_id")
  dx <- (size_x) / n_bins / 2
  dy <- dx * 2 / sqrt(3)
  hex_offsets_x <- c(0, dx, dx, 0, -dx, -dx)
  hex_offsets_y <- c(dy, dy/2, -dy/2, -dy, -dy/2, dy/2)
  poly_list <- vector("list", nrow(bin_majority))
  for (i in seq_len(nrow(bin_majority))) {
    poly_list[[i]] <- data.frame(
      hex_id = bin_majority$hex_id[i], vertex_id = 1:6,
      x = bin_majority$x_c[i] + hex_offsets_x,
      y = bin_majority$y_c[i] + hex_offsets_y,
      label = bin_majority$label[i], stringsAsFactors = FALSE)
  }
  poly_df <- do.call(rbind, poly_list)
  poly_df$group_id <- paste0("hex_", poly_df$hex_id)
  list(bin_df = bin_majority, poly_df = poly_df,
       hex_dx = dx, hex_dy = dy, n_hex = nrow(bin_majority))
}

detect_enrichment_centers <- function(bin_df, target_id, hex_dx = 100,
                                      n_top = 3, min_hex = 5,
                                      neighbor_factor = 3) {
  target_hex <- bin_df[bin_df$label == target_id & !is.na(bin_df$label), , drop = FALSE]
  if (nrow(target_hex) < min_hex) return(NULL)
  cut_h <- hex_dx * neighbor_factor
  if (nrow(target_hex) <= 2) {
    return(data.frame(cluster = 1, x_c = mean(target_hex$x_c),
                      y_c = mean(target_hex$y_c), radius = hex_dx * 5,
                      n_hex = nrow(target_hex), stringsAsFactors = FALSE))
  }
  d <- dist(target_hex[, c("x_c", "y_c")])
  hc <- hclust(d, method = "single")
  cluster_ids <- cutree(hc, h = cut_h)
  target_hex$cluster <- cluster_ids
  centers <- target_hex %>%
    group_by(cluster) %>%
    summarise(x_c = mean(x_c), y_c = mean(y_c), n_hex = n(),
              x_min = min(x_c), x_max = max(x_c),
              y_min = min(y_c), y_max = max(y_c), .groups = "drop") %>%
    mutate(width = x_max - x_min, height = y_max - y_min,
           radius = pmax(width, height) / 2 + hex_dx * 1.5,
           radius = pmin(radius, hex_dx * 30)) %>%
    filter(n_hex >= min_hex) %>%
    arrange(desc(n_hex)) %>%
    slice_head(n = n_top)
  if (nrow(centers) == 0) return(NULL)
  centers
}

# FIX D/E: 肝端标 M05, 脑端标 M04 (BBB)
get_target_microenv <- function(tissue) {
  if (tissue == "liver") "M05" else "M04"
}
get_target_label <- function(tissue) {
  if (tissue == "liver") "Tumor_immune_niche" else "Vascular_BBB"
}

panel_d_subtitle <- function(s, group, tissue, n_microenv) {
  if (group == "tumor") {
    if (tissue == "liver") {
      sprintf("\u2B50 关键发现：红色 (M05 Tumor_immune_niche) 在 tumor 切片中富集，与肿瘤区域空间共定位 [%d 微环境]", n_microenv)
    } else {
      sprintf("\u2B50 关键发现：蓝色 (M04 Vascular_BBB) 在 tumor 切片中的分布变化 [%d 微环境]", n_microenv)
    }
  } else if (group == "baseline") {
    if (tissue == "liver") {
      sprintf("对照 (baseline)：肝端正常微环境分布，无 M05 富集 [%d 微环境]", n_microenv)
    } else {
      sprintf("对照 (baseline)：脑端正常微环境分布 [%d 微环境]", n_microenv)
    }
  } else {
    if (tissue == "liver") {
      sprintf("对照 (control)：肝端正常微环境分布，无 M05 富集 [%d 微环境]", n_microenv)
    } else {
      sprintf("对照 (control)：脑端正常微环境分布 [%d 微环境]", n_microenv)
    }
  }
}

# ======================== Main loop ========================
log_msg(sprintf("=== fix_he_overlay_BC_v3 start (sample=%s) ===", opt$sample))

for (s in names(SAMPLES)) {
  t0 <- Sys.time()
  cfg <- SAMPLES[[s]]
  sample_out <- file.path(OUT_DIR, s)
  dir.create(sample_out, showWarnings = FALSE, recursive = TRUE)
  log_msg(sprintf("---- [%s] %s/%s ----", s, cfg$tissue, cfg$group))

  if (cfg$tissue == "brain") {
    dom_profile <- prof_brain; dom_naming <- name_brain
    me_map <- me_brain; me_colors <- BRAIN_MICROENV_COLORS
  } else {
    dom_profile <- prof_liver; dom_naming <- name_liver
    me_map <- me_liver; me_colors <- LIVER_MICROENV_COLORS
  }

  # --- Official RDS for H&E size ---
  if (!file.exists(cfg$official)) { log_msg("  SKIP: official missing"); next }
  seu <- readRDS(cfg$official)
  info <- seu@misc$info[[1]]
  size_x <- info$size_x; size_y <- info$size_y
  rm(seu); gc(verbose = FALSE)
  log_msg(sprintf("  HE size: %d x %d", size_x, size_y))

  # --- H&E image ---
  if (!file.exists(cfg$he)) { log_msg("  SKIP: HE missing"); next }
  he_img <- png::readPNG(cfg$he)
  he_grob <- rasterGrob(he_img, interpolate = FALSE,
                        width = unit(1, "npc"), height = unit(1, "npc"))
  rm(he_img); gc(verbose = FALSE)

  # --- RCTD results ---
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

  top_ct <- colnames(weights)[max.col(weights, ties.method = "first")]
  weights_smoothed <- smooth_weights(coords, weights, k = KNN_K)
  domain <- assign_domain(weights_smoothed, dom_profile)
  microenv <- me_map$domain_to_microenv[domain]

  log_msg(sprintf("  Domain values: %s",
                  paste(head(sort(unique(domain)), 5), collapse = ", ")))
  log_msg(sprintf("  Microenv dist: %s",
                  paste(names(sort(table(microenv), decreasing = TRUE)), collapse = ",")))

  aspect      <- size_y / size_x
  plot_height <- PLOT_WIDTH * aspect * 1.05
  plot_theme  <- theme_void(base_size = 12) +
    theme(plot.title    = element_text(face = "bold"),
          plot.subtitle = element_text(color = "grey30"),
          legend.position = "right",
          legend.title = element_text(face = "bold"),
          plot.margin  = margin(5, 5, 5, 5))

  # ===== Figure B: domain (FIXED — 用 D_new_XX 原名) =====
  all_doms <- sort(rownames(dom_profile))  # D_new_01 .. D_new_10
  domain_palette <- build_domain_palette(dom_profile)
  domain_legend  <- dom_naming[all_doms]
  domain_legend[is.na(domain_legend)] <- all_doms[is.na(domain_legend)]

  plot_df_B <- data.frame(x = coords$x, y = coords$y, domain = domain,
                          stringsAsFactors = FALSE)
  plot_df_B$domain <- factor(plot_df_B$domain, levels = all_doms)

  log_msg(sprintf("  B factor levels: %s", paste(levels(plot_df_B$domain), collapse = ",")))
  log_msg(sprintf("  B NA count: %d / %d", sum(is.na(plot_df_B$domain)), nrow(plot_df_B)))

  pB <- ggplot(plot_df_B, aes(x = x, y = y, color = domain)) +
    annotation_custom(he_grob, xmin = 0, xmax = size_x, ymin = size_y, ymax = 0) +
    geom_point(size = POINT_SIZE, alpha = POINT_ALPHA) +
    scale_color_manual(values = domain_palette, labels = domain_legend,
                       name = "Spatial domain (10)", drop = FALSE) +
    coord_fixed(xlim = c(0, size_x), ylim = c(size_y, 0), expand = FALSE) +
    labs(title = sprintf("HE + Spatial Domain (10) \u2014 %s (%s)", s, cfg$group),
         subtitle = sprintf("%d spots assigned via cosine similarity", nrow(plot_df_B))) +
    plot_theme +
    guides(color = guide_legend(override.aes = list(size = 4, alpha = 1), ncol = 1))

  ggsave(file.path(sample_out, "B_spatial_domain_on_HE.png"), pB,
         width = PLOT_WIDTH, height = plot_height, dpi = 150,
         limitsize = FALSE, bg = "white")
  log_msg("  B saved (FIXED palette)")
  rm(pB, plot_df_B); gc(verbose = FALSE)

  # ===== Figure C: tumor-induced score (FIXED legend) =====
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

    pC <- ggplot(plot_df_C, aes(x = x, y = y, color = score)) +
      annotation_custom(he_grob, xmin = 0, xmax = size_x, ymin = size_y, ymax = 0) +
      geom_point(size = POINT_SIZE * 1.2, alpha = 0.95) +
      scale_color_gradient2(low = "#1b7837", mid = "#f7f7f7", high = "#b35806",
                            midpoint = median(score, na.rm = TRUE),
                            name = "Tumor-induced\nniche score") +
      coord_fixed(xlim = c(0, size_x), ylim = c(size_y, 0), expand = FALSE) +
      labs(title = sprintf("HE + Tumor-induced niche score \u2014 %s (%s)", s, cfg$group),
           subtitle = sprintf("Sum of %d cell type weights (%s)",
                              length(score_cts),
                              paste(score_cts, collapse = " + "))) +
      plot_theme

    ggsave(file.path(sample_out, "C_tumor_induced_score_on_HE.png"), pC,
           width = PLOT_WIDTH, height = plot_height, dpi = 150,
           limitsize = FALSE, bg = "white")
    log_msg(sprintf("  C saved (score [%.3f, %.3f])",
                    min(score, na.rm = TRUE), max(score, na.rm = TRUE)))
    rm(pC, plot_df_C, score, score_clip); gc(verbose = FALSE)
  } else {
    log_msg("  C SKIP: no score cell types")
  }

  # ===== Figure D: microenvironment (with BBB annotation for brain) =====
  all_me <- sort(unique(microenv))
  all_me <- all_me[!is.na(all_me)]
  me_legend_full <- paste0(all_me, " ", me_map$microenv_to_label[all_me])
  names(me_legend_full) <- all_me
  color_map_D <- me_colors[all_me]
  n_me <- length(all_me)

  plot_df_D <- data.frame(x = coords$x, y = coords$y,
                          microenv_id = unname(microenv), stringsAsFactors = FALSE)
  plot_df_D$microenv_id <- factor(plot_df_D$microenv_id, levels = all_me)

  target_me    <- get_target_microenv(cfg$tissue)
  target_label <- get_target_label(cfg$tissue)
  d_subtitle   <- panel_d_subtitle(s, cfg$group, cfg$tissue, n_me)

  # 检测富集中心 (tumor 样本: liver=M05, brain=M04)
  enrich_centers <- NULL
  if (cfg$group == "tumor" && target_me %in% all_me) {
    hex_tmp <- hex_bin_majority(coords, microenv, n_bins = HEX_GRID_N,
                                size_x = size_x, size_y = size_y)
    enrich_centers <- detect_enrichment_centers(hex_tmp$bin_df, target_me,
                                                hex_dx = hex_tmp$hex_dx, n_top = 3)
    if (!is.null(enrich_centers)) {
      log_msg(sprintf("  D: %d %s center(s), n_hex=[%s]",
                      nrow(enrich_centers), target_me,
                      paste(enrich_centers$n_hex, collapse = ",")))
    }
    rm(hex_tmp); gc(verbose = FALSE)
  }

  pD <- ggplot(plot_df_D, aes(x = x, y = y, color = microenv_id)) +
    annotation_custom(he_grob, xmin = 0, xmax = size_x, ymin = size_y, ymax = 0) +
    geom_point(size = POINT_SIZE, alpha = POINT_ALPHA) +
    scale_color_manual(values = color_map_D, labels = me_legend_full,
                       name = sprintf("Microenvironment\n(%d total)", n_me), drop = FALSE) +
    coord_fixed(xlim = c(0, size_x), ylim = c(size_y, 0), expand = FALSE) +
    labs(title = sprintf("HE + \u5FAE\u73AF\u5883 (%d) \u2014 %s (%s)", n_me, s, cfg$group),
         subtitle = d_subtitle) +
    plot_theme +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1))

  if (!is.null(enrich_centers)) {
    n_pts <- 64; angles <- seq(0, 2 * pi, length.out = n_pts)
    circle_df <- do.call(rbind, lapply(seq_len(nrow(enrich_centers)), function(i) {
      data.frame(circle_id = i,
                 x = enrich_centers$x_c[i] + cos(angles) * enrich_centers$radius[i],
                 y = enrich_centers$y_c[i] + sin(angles) * enrich_centers$radius[i])
    }))
    label_df <- data.frame(
      circle_id = seq_len(nrow(enrich_centers)),
      x_lab = pmin(enrich_centers$x_c + enrich_centers$radius * 0.7, size_x * 0.92),
      y_lab = pmax(enrich_centers$y_c - enrich_centers$radius * 0.85, size_y * 0.05),
      label = if (nrow(enrich_centers) == 1) {
        sprintf("\u2B50 %s\n(%s)", target_me, target_label)
      } else {
        sprintf("\u2B50 %s #%d", target_me, seq_len(nrow(enrich_centers)))
      }
    )
    pD <- pD +
      geom_path(data = circle_df, aes(x = x, y = y, group = circle_id),
                color = "black", linewidth = 1.4, linetype = "dashed",
                inherit.aes = FALSE) +
      geom_text(data = label_df, aes(x = x_lab, y = y_lab, label = label),
                color = "black", fontface = "bold", size = 5,
                hjust = 0, vjust = 0, inherit.aes = FALSE)
  } else if (cfg$group %in% c("baseline", "control")) {
    no_label <- if (cfg$tissue == "liver") {
      "\u2713 \u65E0 M05 \u5BCC\u96C6 (Tumor_immune_niche absent)"
    } else {
      "\u2713 M04 Vascular_BBB baseline distribution"
    }
    pD <- pD + annotate("text", x = size_x * 0.05, y = size_y * 0.05,
                        label = no_label, color = "darkgreen",
                        fontface = "bold", size = 4.5, hjust = 0)
  }

  ggsave(file.path(sample_out, "D_microenvironment_on_HE.png"), pD,
         width = PLOT_WIDTH, height = plot_height, dpi = 150,
         limitsize = FALSE, bg = "white")
  log_msg("  D saved")
  rm(pD); gc(verbose = FALSE)

  # ===== Figure E: hex microenv (with BBB annotation for brain) =====
  hex_res <- hex_bin_majority(coords, microenv, n_bins = HEX_GRID_N,
                              size_x = size_x, size_y = size_y)
  poly_df <- hex_res$poly_df
  poly_df$microenv_id <- factor(poly_df$label, levels = all_me)

  enrich_centers_E <- NULL
  if (cfg$group == "tumor" && target_me %in% all_me) {
    enrich_centers_E <- detect_enrichment_centers(hex_res$bin_df, target_me,
                                                  hex_dx = hex_res$hex_dx, n_top = 3)
  }

  pE <- ggplot(poly_df, aes(x = x, y = y, group = group_id, fill = microenv_id)) +
    annotation_custom(he_grob, xmin = 0, xmax = size_x, ymin = size_y, ymax = 0) +
    geom_polygon(alpha = 0.78, color = NA) +
    scale_fill_manual(values = color_map_D, labels = me_legend_full,
                      name = sprintf("Microenvironment\n(%d total)", n_me), drop = FALSE) +
    coord_fixed(xlim = c(0, size_x), ylim = c(size_y, 0), expand = FALSE) +
    labs(title = sprintf("HE + Hex \u5FAE\u73AF\u5883\u5730\u56FE \u2014 %s (%s)", s, cfg$group),
         subtitle = d_subtitle) +
    plot_theme +
    guides(fill = guide_legend(override.aes = list(alpha = 1), ncol = 1))

  if (!is.null(enrich_centers_E)) {
    n_pts <- 64; angles <- seq(0, 2 * pi, length.out = n_pts)
    circle_df <- do.call(rbind, lapply(seq_len(nrow(enrich_centers_E)), function(i) {
      data.frame(circle_id = i,
                 x = enrich_centers_E$x_c[i] + cos(angles) * enrich_centers_E$radius[i],
                 y = enrich_centers_E$y_c[i] + sin(angles) * enrich_centers_E$radius[i])
    }))
    label_df <- data.frame(
      circle_id = seq_len(nrow(enrich_centers_E)),
      x_lab = pmin(enrich_centers_E$x_c + enrich_centers_E$radius * 0.7, size_x * 0.92),
      y_lab = pmax(enrich_centers_E$y_c - enrich_centers_E$radius * 0.85, size_y * 0.05),
      label = if (nrow(enrich_centers_E) == 1) {
        sprintf("\u2B50 %s\n(%s)", target_me, target_label)
      } else {
        sprintf("\u2B50 %s #%d", target_me, seq_len(nrow(enrich_centers_E)))
      }
    )
    pE <- pE +
      geom_path(data = circle_df, aes(x = x, y = y, group = circle_id),
                color = "black", linewidth = 1.4, linetype = "dashed",
                inherit.aes = FALSE) +
      geom_text(data = label_df, aes(x = x_lab, y = y_lab, label = label),
                color = "black", fontface = "bold", size = 5,
                hjust = 0, vjust = 0, inherit.aes = FALSE)
  } else if (cfg$group %in% c("baseline", "control")) {
    no_label <- if (cfg$tissue == "liver") {
      "\u2713 \u65E0 M05 \u5BCC\u96C6"
    } else {
      "\u2713 M04 BBB baseline"
    }
    pE <- pE + annotate("text", x = size_x * 0.05, y = size_y * 0.05,
                        label = no_label, color = "darkgreen",
                        fontface = "bold", size = 4.5, hjust = 0)
  }

  ggsave(file.path(sample_out, "E_hex_microenv_on_HE.png"), pE,
         width = PLOT_WIDTH, height = plot_height, dpi = 150,
         limitsize = FALSE, bg = "white")
  log_msg(sprintf("  E saved (%d hex)", hex_res$n_hex))

  rm(pE, poly_df, hex_res, plot_df_D, coords, weights, weights_smoothed,
     he_grob, domain, microenv)
  gc(verbose = FALSE)
  log_msg(sprintf("  [%s] done in %.1f min",
                  s, difftime(Sys.time(), t0, units = "mins")))
}

log_msg("=== fix_he_overlay_BC_v3 ALL DONE ===")
