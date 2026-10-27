#!/usr/bin/env Rscript
source("/mnt/BioDisk/wangmingyuan/SGS10_analysis/scripts/fig_save.R")  # [restyle] 三格式输出器
# ============================================================
# he_overlay_RCTD100.R  (restyle: 100% RCTD, 删标题/放大字号/三格式)
# ============================================================
# Purpose:
#   Generate paper-grade H&E + spatial overlay images per sample.
#   5 panels per sample (A: top CT, B: domain, C: tumor score,
#                       D: microenv overlay, E: hex microenv map)
#
# Inputs:
#   - Official Seurat .rds with HE info
#   - HE aligned PNG
#   - 100% RCTD: spatial_results/rctd_full/per_sample/<s>_rctd_results.rds
#   - Domain profile + naming: taskB/spatial_domain_RCTD100/
#   - Microenv naming: taskB/spatial_domain_RCTD100/<organ>_microenv_naming.csv
#
# Outputs: figures/restyle_v1/<sample>/  (PDF + TIFF + PNG@300)
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

# ---------- CLI ----------
option_list <- list(
  make_option("--sample", type="character", default="all",
              help="sample name or 'all'"),
  make_option("--out_dir", type="character",
              default=path.expand("~/SGS10_analysis/taskB_he_overlay_RCTD100"))
)
opt <- parse_args(OptionParser(option_list=option_list))

# ---------- Config ----------
PR <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
OUT_DIR <- file.path(PR, "figures/restyle_v1")   # [restyle] 新图统一目录
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
LOG_FILE <- file.path(OUT_DIR, "he_overlay_RCTD100_run.log")
if (file.exists(LOG_FILE) && opt$sample == "all") file.remove(LOG_FILE)

log_msg <- function(msg) {
  ts <- format(Sys.time(), "%H:%M:%S")
  line <- sprintf("[%s] %s", ts, msg)
  cat(line, "\n")
  cat(line, "\n", file = LOG_FILE, append = TRUE)
}

log_msg(sprintf("=== HE overlay RCTD100 start (sample=%s) ===", opt$sample))

KNN_K       <- 15
POINT_SIZE  <- 0.6
POINT_ALPHA <- 0.9
PLOT_WIDTH  <- 16
HEX_GRID_N  <- 200

# ---------- Paths ----------
BD <- "/mnt/BioData/bio/SGS10例"
MATRIX_LIVER <- file.path(BD, "SGS25090559-8例样本数据/data/matrix")
MATRIX_NAO56 <- file.path(BD, "SGS24070263-nao56/data/matrix/24070263_nao_56_expression")
MATRIX_GAN56 <- file.path(BD, "SGS24070264-gan56/data/matrix/24070264_gan_56_expression")
RDS_LIVER    <- file.path(BD, "SGS25090559-8例样本数据/data/rds")
RDS_NAO56    <- file.path(BD, "SGS24070263-nao56/data/rds/24070263_nao_56_expression")
RDS_GAN56    <- file.path(BD, "SGS24070264-gan56/data/rds/24070264_gan_56_expression")

SAMPLES <- list(   # [restyle] 只保留 docx 需要的两个样本 (img5/6 Top-CT, img7 微环境, img17 Hex)
  Brain_3 = list(tissue = "brain", group = "tumor", official = file.path(RDS_LIVER, "25090559_Brain_3_expression/Brain_3.rds"), he = file.path(MATRIX_LIVER, "25090559_Brain_3_expression/Brain_3_aligned_HE_TIMG.png")),
  Liver_3 = list(tissue = "liver", group = "tumor", official = file.path(RDS_LIVER, "25090559_Liver_3_expression/Liver_3.rds"), he = file.path(MATRIX_LIVER, "25090559_Liver_3_expression/Liver_3_aligned_HE_TIMG.png"))
)

if (opt$sample != "all") {
  if (!opt$sample %in% names(SAMPLES)) stop("Unknown sample: ", opt$sample)
  SAMPLES <- SAMPLES[opt$sample]
}

RCTD_DIR     <- file.path(PR, "spatial_results/rctd_full/per_sample")
SD_DIR       <- file.path(PR, "taskB/spatial_domain_RCTD100")

# ---------- Load domain profile + naming ----------
load_domain_profile <- function(f) {
  d <- read.csv(f, check.names = FALSE)
  rownames(d) <- d[[1]]
  d[[1]] <- NULL
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

prof_brain  <- load_domain_profile(file.path(SD_DIR, "brain_domain_ct_profile.csv"))
prof_liver  <- load_domain_profile(file.path(SD_DIR, "liver_domain_ct_profile.csv"))
name_brain  <- load_naming(file.path(SD_DIR, "brain_domain_naming.csv"))
name_liver  <- load_naming(file.path(SD_DIR, "liver_domain_naming.csv"))
me_brain    <- load_microenv("brain")
me_liver    <- load_microenv("liver")

log_msg(sprintf("Profiles: brain %dx%d, liver %dx%d",
                nrow(prof_brain), ncol(prof_brain),
                nrow(prof_liver), ncol(prof_liver)))
log_msg(sprintf("Microenv: brain %d, liver %d",
                length(me_brain$microenv_to_label),
                length(me_liver$microenv_to_label)))

# ---------- Color palettes ----------
CT_COLORS_BRAIN <- c(
  "Hpc_CA1"="#e41a1c","Hpc_CA3"="#a65628","Hpc_DG"="#ff7f00",
  "Excitatory_L23"="#4daf4a","Excitatory_L45"="#66a61e","Excitatory_L6"="#006d2c",
  "Neurons"="#984ea3","Inhibitory_Meis2"="#f781bf","Inhibitory_Pvalb"="#b3de69",
  "Inhibitory_Sst"="#fdb462","Inhibitory_Vip"="#80b1d3",
  "Oligo_Mature"="#377eb8","OPCs"="#1f78b4",
  "Astrocyte_Resting"="#ffff33","Astrocyte_Gfap"="#bebada",
  "Microglia_Homeostatic"="#8c510a","Microglia_other"="#dfc27d",
  "Endothelial"="#fb8072","Pericytes"="#bf5b17","VLMC"="#737373",
  "CPECs"="#2b8cbe","Ependymal"="#01665e"
)
CT_COLORS_LIVER <- c(
  "Hep_Periportal"="#377eb8","Hep_Perivenous"="#984ea3","Hep_Proliferating"="#6a3d9a",
  "Pericytes"="#e41a1c","Stellate"="#b15928","HSC_Quiescent"="#33a02c",
  "HSC_Activated_ECM"="#a65628","HSC_Activated_Chemo"="#c2a5cf","HSC_Proliferating"="#762a83",
  "Kupffer"="#ff7f00","KC_Homeostatic"="#fdbf6f","Macrophage_M1"="#b35806",
  "Macrophage_M2"="#fee090","Macrophage_NAM"="#fdae61","Monocytes"="#f46d43",
  "Neutrophils"="#abd9e9","cDC1"="#4575b4","cDC2"="#5e4fa2","pDC"="#3288bd",
  "Dendritic"="#66c2a5","B_cells"="#d6604d","B_Plasma"="#4daf4a",
  "CD4_T"="#f781bf","CD8_T"="#e78ac3","T_cells"="#fccde5","Treg"="#bc80bd",
  "NK_cells"="#8dd3c7","NKT"="#80cdc1",
  "LSEC"="#ffff33","LSEC_Mid"="#ffeda0","LSEC_Pericentral"="#fed976","LSEC_Periportal"="#feb24c",
  "Endothelial"="#1b9e77","Cholangiocytes"="#e5f5e0","Erythrocytes"="#7f0000","LCM"="#c49c94"
)
OTHER_COLOR <- "#BDBDBD"

DOMAIN_COLORS <- c(
  "D01"="#e41a1c","D02"="#377eb8","D03"="#4daf4a","D04"="#984ea3",
  "D05"="#ff7f00","D06"="#ffd92f","D07"="#a65628","D08"="#f781bf",
  "D09"="#66c2a5","D10"="#bebada"
)
LIVER_MICROENV_COLORS <- c(
  M01 = "#2C7FB8",  M02 = "#7FCDBB",  M03 = "#41B6C4",
  M04 = "#B15928",  M05 = "#E31A1C",  # paper-key niche
  M06 = "#FDAE6B",  M07 = "#9E9AC8",  M08 = "#cab2d6",
  M09 = "#fb9a99"
)
BRAIN_MICROENV_COLORS <- c(
  M01 = "#E31A1C",  M02 = "#33A02C",  M03 = "#A6761D",
  M04 = "#1F78B4",  M05 = "#FF7F00",  M06 = "#984EA3"
)

# ---------- Helpers ----------
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
  hb <- hexbin::hexbin(coords$x, coords$y,
                       xbins = n_bins,
                       xbnds = c(0, size_x), ybnds = c(0, size_y),
                       IDs = TRUE)
  spot_hex_id <- hb@cID
  bin_df <- data.frame(hex_id = spot_hex_id, label = label_vec,
                       stringsAsFactors = FALSE)
  bin_majority <- bin_df %>%
    group_by(hex_id) %>%
    summarise(n_spots = n(),
              label = names(sort(table(label), decreasing = TRUE))[1],
              .groups = "drop") %>%
    mutate(label = ifelse(is.na(label), "Unassigned", label))

  hex_centers <- hexbin::hcell2xy(hb)
  hex_centers_df <- data.frame(hex_id = hb@cell, x_c = hex_centers$x, y_c = hex_centers$y)
  bin_majority <- bin_majority %>% left_join(hex_centers_df, by = "hex_id")

  dx <- (size_x - 0) / n_bins / 2
  dy <- dx * 2 / sqrt(3)
  hex_offsets_x <- c(0, dx, dx, 0, -dx, -dx)
  hex_offsets_y <- c(dy, dy/2, -dy/2, -dy, -dy/2, dy/2)

  poly_list <- vector("list", nrow(bin_majority))
  for (i in seq_len(nrow(bin_majority))) {
    poly_list[[i]] <- data.frame(
      hex_id = bin_majority$hex_id[i],
      vertex_id = 1:6,
      x = bin_majority$x_c[i] + hex_offsets_x,
      y = bin_majority$y_c[i] + hex_offsets_y,
      label = bin_majority$label[i],
      stringsAsFactors = FALSE)
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
    return(data.frame(cluster = 1,
                      x_c = mean(target_hex$x_c), y_c = mean(target_hex$y_c),
                      radius = hex_dx * 5, n_hex = nrow(target_hex),
                      stringsAsFactors = FALSE))
  }

  d <- dist(target_hex[, c("x_c", "y_c")])
  hc <- hclust(d, method = "single")
  cluster_ids <- cutree(hc, h = cut_h)
  target_hex$cluster <- cluster_ids

  centers <- target_hex %>%
    group_by(cluster) %>%
    summarise(x_c = mean(x_c), y_c = mean(y_c),
              n_hex = n(),
              x_min = min(x_c), x_max = max(x_c),
              y_min = min(y_c), y_max = max(y_c),
              .groups = "drop") %>%
    mutate(width = x_max - x_min, height = y_max - y_min,
           span = pmax(width, height),
           radius = pmax(span / 2 * 1.3, hex_dx * 8),
           radius = pmin(radius, hex_dx * 25)) %>%
    filter(n_hex >= min_hex) %>%
    arrange(desc(n_hex)) %>%
    slice_head(n = n_top)

  if (nrow(centers) == 0) return(NULL)
  centers
}

get_target_microenv <- function(tissue) {
  if (tissue == "liver") "M05" else NA
}

# ============================================================
# Main loop
# ============================================================
for (s in names(SAMPLES)) {
  sample_t0 <- Sys.time()
  cfg <- SAMPLES[[s]]
  sample_out <- file.path(OUT_DIR, s)
  dir.create(sample_out, showWarnings = FALSE, recursive = TRUE)
  log_msg(sprintf("---- [%s] %s/%s ----", s, cfg$tissue, cfg$group))

  if (cfg$tissue == "brain") {
    ct_colors    <- CT_COLORS_BRAIN
    dom_profile  <- prof_brain
    dom_naming   <- name_brain
    me_map       <- me_brain
    me_colors    <- BRAIN_MICROENV_COLORS
  } else {
    ct_colors    <- CT_COLORS_LIVER
    dom_profile  <- prof_liver
    dom_naming   <- name_liver
    me_map       <- me_liver
    me_colors    <- LIVER_MICROENV_COLORS
  }

  if (!file.exists(cfg$official)) { log_msg("  SKIP: official missing"); next }
  seu <- readRDS(cfg$official)
  info <- seu@misc$info[[1]]
  size_x <- info$size_x; size_y <- info$size_y
  rm(seu); gc(verbose = FALSE)
  log_msg(sprintf("  size: %d x %d", size_x, size_y))
  ROTATE180 <- (cfg$tissue == "brain")   # [client] 脑 HE 系列整体转 180度
  if (ROTATE180) log_msg("  ROTATE180: brain HE + points mirrored")

  if (!file.exists(cfg$he)) { log_msg("  SKIP: HE missing"); next }
  he_img <- png::readPNG(cfg$he)
  if (ROTATE180) {
    .d <- dim(he_img)
    he_img <- if (length(.d) == 3) he_img[.d[1]:1, .d[2]:1, , drop = FALSE] else he_img[.d[1]:1, .d[2]:1]
  }
  he_grob <- rasterGrob(he_img, interpolate = FALSE,
                        width = unit(1, "npc"), height = unit(1, "npc"))
  rm(he_img); gc(verbose = FALSE)

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
  if (ROTATE180) { coords$x <- size_x - coords$x; coords$y <- size_y - coords$y }  # [client] 点坐标镜像

  top_ct_idx <- max.col(weights, ties.method = "first")
  top_ct     <- colnames(weights)[top_ct_idx]

  weights_smoothed <- smooth_weights(coords, weights, k = KNN_K)
  domain <- assign_domain(weights_smoothed, dom_profile)
  log_msg(sprintf("  Domain dist: %s",
                  paste(names(sort(table(domain), decreasing=TRUE)), collapse=",")))

  microenv <- me_map$domain_to_microenv[domain]
  log_msg(sprintf("  Microenv dist: %s",
                  paste(names(sort(table(microenv), decreasing=TRUE)), collapse=",")))

  # [restyle] 标题/副标题(含中文)移除; 字号放大
  plot_theme <- theme_void(base_size = 14) +
    theme(plot.title    = element_blank(),
          plot.subtitle = element_blank(),
          legend.position = "right",
          legend.title  = element_text(size = 15, face = "bold"),
          legend.text   = element_text(size = 13, face = "bold"),
          plot.margin = margin(5, 5, 5, 5))
  aspect <- size_y / size_x
  plot_height <- PLOT_WIDTH * aspect * 1.05

  # ===== Figure A: top CT =====
  plot_df_A <- data.frame(x = coords$x, y = coords$y, ct = top_ct)
  present_ct <- sort(unique(top_ct))
  present_ct <- intersect(present_ct, names(ct_colors))
  plot_df_A$ct_label <- ifelse(plot_df_A$ct %in% present_ct, plot_df_A$ct, "Other")
  plot_df_A$ct_label <- factor(plot_df_A$ct_label, levels = c(present_ct, "Other"))
  color_map_A <- c(ct_colors[present_ct], Other = OTHER_COLOR)

  pA <- ggplot(plot_df_A, aes(x = x, y = y, color = ct_label)) +
    annotation_custom(he_grob, xmin = 0, xmax = size_x, ymin = size_y, ymax = 0) +
    geom_point(size = POINT_SIZE, alpha = POINT_ALPHA) +
    scale_color_manual(values = color_map_A, name = "Top cell type", drop = FALSE) +
    coord_fixed(xlim = c(0, size_x), ylim = c(size_y, 0), expand = FALSE) +
    plot_theme +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = if (cfg$tissue == "brain") 6 else 7, byrow = TRUE, title.position = "top")) +
    theme(legend.position = "bottom", legend.justification = if (cfg$tissue == "brain") "center" else c(0, 0.5), legend.key.spacing.x = unit(if (cfg$tissue == "brain") 8 else 20, "pt"))  # LEGEND_BOTTOM_A LEGEND_FIX2 LEGEND_SPREAD LEGEND_FONTUP

  save_fig_all(pA, file.path(sample_out, "A_top_CT_on_HE"),
               width = PLOT_WIDTH, height = plot_height + 3)   # [restyle] PDF+TIFF+PNG@300
  log_msg("  A saved")
  rm(pA, plot_df_A); gc(verbose = FALSE)

  # ===== Figure HE: 纯 H&E 底图(无叠加), 与 A/D/E 同几何, 供 docx img3/4 拼图对齐 =====
  pHE <- ggplot() +
    annotation_custom(he_grob, xmin = 0, xmax = size_x, ymin = size_y, ymax = 0) +
    geom_blank(data = data.frame(x = c(0, size_x), y = c(0, size_y)),
               aes(x = x, y = y)) +
    coord_fixed(xlim = c(0, size_x), ylim = c(size_y, 0), expand = FALSE) +
    theme_void()
  save_fig_all(pHE, file.path(sample_out, "HE_only"),
               width = PLOT_WIDTH, height = plot_height)   # [restyle] PDF+TIFF+PNG@300
  log_msg("  HE saved")
  rm(pHE); gc(verbose = FALSE)

  # [restyle] Figure B(spatial domain) 与 C(tumor score) 已移除: docx 未使用, 且 B 配色异常

  # ===== Figure D: microenvironment =====
  plot_df_D <- data.frame(x = coords$x, y = coords$y,
                          microenv_id = unname(microenv), stringsAsFactors = FALSE)
  all_me <- sort(unique(plot_df_D$microenv_id))
  plot_df_D$microenv_id <- factor(plot_df_D$microenv_id, levels = all_me)
  me_legend_full <- paste0(all_me, " ", me_map$microenv_to_label[all_me])
  names(me_legend_full) <- all_me
  color_map_D <- me_colors[all_me]
  n_me <- length(all_me)

  target_me <- get_target_microenv(cfg$tissue)

  enrich_centers <- NULL
  if (cfg$group == "tumor" && !is.na(target_me) && target_me %in% all_me) {
    hex_for_center <- hex_bin_majority(coords, microenv, n_bins = HEX_GRID_N,
                                        size_x = size_x, size_y = size_y)
    enrich_centers <- detect_enrichment_centers(hex_for_center$bin_df, target_me,
                                                 hex_dx = hex_for_center$hex_dx, n_top = 3)
    if (!is.null(enrich_centers)) {
      log_msg(sprintf("  D: %d %s enrichment center(s), n_hex=[%s]",
                      nrow(enrich_centers), target_me,
                      paste(enrich_centers$n_hex, collapse=",")))
    }
    rm(hex_for_center); gc(verbose = FALSE)
  }

  pD <- ggplot(plot_df_D, aes(x = x, y = y, color = microenv_id)) +
    annotation_custom(he_grob, xmin = 0, xmax = size_x, ymin = size_y, ymax = 0) +
    geom_point(size = POINT_SIZE, alpha = POINT_ALPHA) +
    scale_color_manual(values = color_map_D, labels = me_legend_full,
                       name = sprintf("Microenvironment (%d total)", n_me),
                       drop = FALSE) +
    coord_fixed(xlim = c(0, size_x), ylim = c(size_y, 0), expand = FALSE) +
    plot_theme +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 6, byrow = TRUE, title.position = "top")) +
    theme(legend.position = "bottom", legend.justification = "center", legend.key.spacing.x = unit(5, "pt"))  # LEGEND_BOTTOM_DE LEGEND_SPREAD LEGEND_FONTUP

  if (!is.null(enrich_centers)) {
    n_pts <- 64
    angles <- seq(0, 2*pi, length.out = n_pts)
    circle_df <- do.call(rbind, lapply(seq_len(nrow(enrich_centers)), function(i) {
      data.frame(circle_id = i,
                 x = enrich_centers$x_c[i] + cos(angles) * enrich_centers$radius[i],
                 y = enrich_centers$y_c[i] + sin(angles) * enrich_centers$radius[i])
    }))
    label_df <- data.frame(
      circle_id = seq_len(nrow(enrich_centers)),
      x_lab = pmin(enrich_centers$x_c + enrich_centers$radius * 0.7, size_x * 0.92),
      y_lab = pmax(enrich_centers$y_c - enrich_centers$radius * 0.85, size_y * 0.05),
      label = if (nrow(enrich_centers) == 1) sprintf("⭐ %s\n(Tumor-induced niche)", target_me)
              else sprintf("⭐ %s #%d", target_me, seq_len(nrow(enrich_centers)))
    )
    pD <- pD +
      geom_path(data = circle_df, aes(x = x, y = y, group = circle_id),
                color = "black", linewidth = 1.4, linetype = "dashed",
                inherit.aes = FALSE) +
      geom_text(data = label_df, aes(x = x_lab, y = y_lab, label = label),
                color = "black", fontface = "bold", size = 5,
                hjust = 0, vjust = 0, inherit.aes = FALSE)
  } else if (cfg$group %in% c("baseline", "control") && cfg$tissue == "liver") {
    pD <- pD + annotate("text",
                        x = size_x * 0.05, y = size_y * 0.05,
                        label = "✓ 无 M05 富集 (Tumor-immune niche absent)",
                        color = "darkgreen", fontface = "bold", size = 4.5, hjust = 0)
  }

  save_fig_all(pD, file.path(sample_out, "D_microenvironment_on_HE"),
               width = PLOT_WIDTH, height = plot_height + 2)   # [restyle] PDF+TIFF+PNG@300
  log_msg("  D saved")
  rm(pD); gc(verbose = FALSE)

  # ===== Figure E: hex microenv map =====
  hex_res <- hex_bin_majority(coords, microenv, n_bins = HEX_GRID_N,
                               size_x = size_x, size_y = size_y)
  poly_df <- hex_res$poly_df
  poly_df$microenv_id <- factor(poly_df$label, levels = all_me)
  # [client] M04 等"散点有、hex 无"的微环境：补不可见多边形(画外/零面积)，使图例每项都取到色
  .missing_me <- setdiff(as.character(all_me), as.character(unique(poly_df$label)))
  if (length(.missing_me) > 0) {
    .ph <- do.call(rbind, lapply(.missing_me, function(m) {
      data.frame(hex_id = paste0("ph_", m), vertex_id = 1:3,
                 x = -1e5, y = -1e5, label = m,
                 group_id = paste0("ph_", m),
                 microenv_id = factor(m, levels = all_me),
                 stringsAsFactors = FALSE)
    }))
    poly_df <- rbind(poly_df, .ph[, names(poly_df)])
  }

  enrich_centers_E <- NULL
  if (cfg$group == "tumor" && !is.na(target_me) && target_me %in% all_me) {
    enrich_centers_E <- detect_enrichment_centers(hex_res$bin_df, target_me,
                                                    hex_dx = hex_res$hex_dx, n_top = 3)
    if (!is.null(enrich_centers_E)) {
      log_msg(sprintf("  E: %d %s enrichment center(s)",
                      nrow(enrich_centers_E), target_me))
    }
  }

  pE <- ggplot(poly_df, aes(x = x, y = y, group = group_id, fill = microenv_id)) +
    annotation_custom(he_grob, xmin = 0, xmax = size_x, ymin = size_y, ymax = 0) +
    geom_polygon(alpha = 0.78, color = NA) +
    scale_fill_manual(values = color_map_D, labels = me_legend_full,
                      name = sprintf("Microenvironment (%d total)", n_me),
                      drop = FALSE, na.value = "#D9D9D9") +
    coord_fixed(xlim = c(0, size_x), ylim = c(size_y, 0), expand = FALSE) +
    plot_theme +
    guides(fill = guide_legend(override.aes = list(alpha = 1), ncol = 6, byrow = TRUE, title.position = "top")) +
    theme(legend.position = "bottom", legend.justification = "center", legend.key.spacing.x = unit(5, "pt"))

  if (!is.null(enrich_centers_E)) {
    n_pts <- 64
    angles <- seq(0, 2*pi, length.out = n_pts)
    circle_df <- do.call(rbind, lapply(seq_len(nrow(enrich_centers_E)), function(i) {
      data.frame(circle_id = i,
                 x = enrich_centers_E$x_c[i] + cos(angles) * enrich_centers_E$radius[i],
                 y = enrich_centers_E$y_c[i] + sin(angles) * enrich_centers_E$radius[i])
    }))
    label_df <- data.frame(
      circle_id = seq_len(nrow(enrich_centers_E)),
      x_lab = pmin(enrich_centers_E$x_c + enrich_centers_E$radius * 0.7, size_x * 0.92),
      y_lab = pmax(enrich_centers_E$y_c - enrich_centers_E$radius * 0.85, size_y * 0.05),
      label = if (nrow(enrich_centers_E) == 1) sprintf("⭐ %s\n(Tumor-induced niche)", target_me)
              else sprintf("⭐ %s #%d", target_me, seq_len(nrow(enrich_centers_E)))
    )
    pE <- pE +
      geom_path(data = circle_df, aes(x = x, y = y, group = circle_id),
                color = "black", linewidth = 1.4, linetype = "dashed",
                inherit.aes = FALSE) +
      geom_text(data = label_df, aes(x = x_lab, y = y_lab, label = label),
                color = "black", fontface = "bold", size = 5,
                hjust = 0, vjust = 0, inherit.aes = FALSE)
  } else if (cfg$group %in% c("baseline", "control") && cfg$tissue == "liver") {
    pE <- pE + annotate("text",
                        x = size_x * 0.05, y = size_y * 0.05,
                        label = "✓ 无 M05 富集 (Tumor-immune niche absent)",
                        color = "darkgreen", fontface = "bold", size = 4.5, hjust = 0)
  }

  save_fig_all(pE, file.path(sample_out, "E_hex_microenv_on_HE"),
               width = PLOT_WIDTH, height = plot_height + 2)   # [restyle] PDF+TIFF+PNG@300
  log_msg(sprintf("  E saved (%d hex bins)", hex_res$n_hex))
  rm(pE, poly_df, hex_res, plot_df_D, coords, weights, weights_smoothed, he_grob)
  gc(verbose = FALSE)

  elapsed <- as.numeric(difftime(Sys.time(), sample_t0, units = "secs"))
  log_msg(sprintf("  [%s] done in %.1fs", s, elapsed))
}

log_msg("=== All done ===")
log_msg(sprintf("Output: %s", OUT_DIR))
