#!/usr/bin/env Rscript
# ============================================================
# HE overlay v5 — 5 plots per sample with consistent colors:
#   A: per-spot top CT (复用 v4)
#   B: per-spot 10 spatial domain (复用 v4 cosine assign)
#   C: per-spot tumor-induced score (复用 v4)
#   D: per-spot microenvironment (新增, 7 肝 / 6 脑)
#   E: hex-binned microenvironment majority (新增, 修复 RCTD 20% 抽样导致的散点感)
# 
# Inputs:
#   原 v4 全部输入 +
#   D1_v2_RCTD_composition/{liver,brain}_microenv_naming.csv  (manually curated)
#
# Outputs:
#   taskB_he_overlay_v5/{sample}/{A,B,C,D,E}_*.png
#
# 用法:
#   单 sample 测试: Rscript he_overlay_v5_all10.R --sample Bla_brain_1
#   全量:        Rscript he_overlay_v5_all10.R
#
# 时间: 单 sample 1-2 分钟; 全 10 sample 15-25 分钟
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(png)
  library(grid)
  library(dplyr)
  library(RANN)
  library(optparse)
})

# ---------- CLI ----------
option_list <- list(
  make_option("--sample", type="character", default="all",
              help="sample name or 'all'"),
  make_option("--out_dir", type="character",
              default=path.expand("~/SGS10_analysis/taskB_he_overlay_v5"))
)
opt <- parse_args(OptionParser(option_list=option_list))

# ---------- Config ----------
OUT_DIR <- path.expand(opt$out_dir)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
LOG_FILE <- file.path(OUT_DIR, "run.log")
if (file.exists(LOG_FILE) && opt$sample == "all") file.remove(LOG_FILE)

log_msg <- function(msg) {
  ts <- format(Sys.time(), "%H:%M:%S")
  line <- sprintf("[%s] %s", ts, msg)
  cat(line, "\n")
  cat(line, "\n", file = LOG_FILE, append = TRUE)
}

log_msg(sprintf("=== HE overlay v5 start (sample=%s) ===", opt$sample))

KNN_K       <- 15
POINT_SIZE  <- 0.6
POINT_ALPHA <- 0.9
PLOT_WIDTH  <- 16
HEX_GRID_N  <- 200   # panel E hex grid resolution (200×200, fine hex)

# ---------- Paths (复用 v4) ----------
BD <- "/mnt/BioData/bio/SGS10例"
MATRIX_LIVER <- file.path(BD, "SGS25090559-8例样本数据/data/matrix")
MATRIX_NAO56 <- file.path(BD, "SGS24070263-nao56/data/matrix/24070263_nao_56_expression")
MATRIX_GAN56 <- file.path(BD, "SGS24070264-gan56/data/matrix/24070264_gan_56_expression")
RDS_LIVER    <- file.path(BD, "SGS25090559-8例样本数据/data/rds")
RDS_NAO56    <- file.path(BD, "SGS24070263-nao56/data/rds/24070263_nao_56_expression")
RDS_GAN56    <- file.path(BD, "SGS24070264-gan56/data/rds/24070264_gan_56_expression")

SAMPLES <- list(
  nao56       = list(tissue = "brain", group = "tumor",    official = file.path(RDS_NAO56, "nao_56.rds"),                                                           he = file.path(MATRIX_NAO56, "nao_56_aligned_HE_TIMG.png")),
  gan56       = list(tissue = "liver", group = "tumor",    official = file.path(RDS_GAN56, "gan_56.rds"),                                                           he = file.path(MATRIX_GAN56, "gan_56_aligned_HE_TIMG.png")),
  Brain_3     = list(tissue = "brain", group = "tumor",    official = file.path(RDS_LIVER, "25090559_Brain_3_expression/Brain_3.rds"),                              he = file.path(MATRIX_LIVER, "25090559_Brain_3_expression/Brain_3_aligned_HE_TIMG.png")),
  Brain_8     = list(tissue = "brain", group = "tumor",    official = file.path(RDS_LIVER, "25090559_Brain_8_expression/Brain_8.rds"),                              he = file.path(MATRIX_LIVER, "25090559_Brain_8_expression/Brain_8_aligned_HE_TIMG.png")),
  Bla_brain_1 = list(tissue = "brain", group = "baseline", official = file.path(RDS_LIVER, "25090559_Bla_brain_1_expression/Bla_brain_1.rds"),                      he = file.path(MATRIX_LIVER, "25090559_Bla_brain_1_expression/Bla_brain_1_aligned_HE_TIMG.png")),
  Con_brain_1 = list(tissue = "brain", group = "control",  official = file.path(RDS_LIVER, "25090559_Con_brain_1_expression/Con_brain_1.rds"),                      he = file.path(MATRIX_LIVER, "25090559_Con_brain_1_expression/Con_brain_1_aligned_HE_TIMG.png")),
  Liver_3     = list(tissue = "liver", group = "tumor",    official = file.path(RDS_LIVER, "25090559_Liver_3_expression/Liver_3.rds"),                              he = file.path(MATRIX_LIVER, "25090559_Liver_3_expression/Liver_3_aligned_HE_TIMG.png")),
  Liver_8     = list(tissue = "liver", group = "tumor",    official = file.path(RDS_LIVER, "25090559_Liver_8_expression/Liver_8.rds"),                              he = file.path(MATRIX_LIVER, "25090559_Liver_8_expression/Liver_8_aligned_HE_TIMG.png")),
  Bla_liver_1 = list(tissue = "liver", group = "baseline", official = file.path(RDS_LIVER, "25090559_Bla_liver_1_expression/Bla_liver_1.rds"),                      he = file.path(MATRIX_LIVER, "25090559_Bla_liver_1_expression/Bla_liver_1_aligned_HE_TIMG.png")),
  Con_liver_1 = list(tissue = "liver", group = "control",  official = file.path(RDS_LIVER, "25090559_Con_liver_1_expression/Con_liver_1.rds"),                      he = file.path(MATRIX_LIVER, "25090559_Con_liver_1_expression/Con_liver_1_aligned_HE_TIMG.png"))
)

if (opt$sample != "all") {
  if (!opt$sample %in% names(SAMPLES)) {
    stop("Unknown sample: ", opt$sample)
  }
  SAMPLES <- SAMPLES[opt$sample]
}

RCTD_DIR <- path.expand("~/SGS10_analysis/spatial_results/rctd/per_sample")
T4_DIR   <- path.expand("~/SGS10_analysis/SGS10_deliverables_2026-04-20/08_spatial_Task_B/T4_spatial_domain")
T1_DIR   <- path.expand("~/SGS10_analysis/SGS10_deliverables_2026-04-20/08_spatial_Task_B/T1_niche_10sample")
MICROENV_DIR <- path.expand("~/SGS10_analysis/D1_v2_RCTD_composition")

# ---------- Load T4 domain profiles (复用 v4) ----------
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
load_delta <- function(f) {
  d <- read.csv(f, check.names = FALSE); rownames(d) <- d[[1]]; d[[1]] <- NULL; as.matrix(d)
}

prof_brain <- load_domain_profile(file.path(T4_DIR, "brain_domain_ct_profile.csv"))
prof_liver <- load_domain_profile(file.path(T4_DIR, "liver_domain_ct_profile.csv"))
name_brain <- load_naming(file.path(T4_DIR, "brain_domain_naming.csv"))
name_liver <- load_naming(file.path(T4_DIR, "liver_domain_naming.csv"))
delta_brain <- load_delta(file.path(T1_DIR, "brain_niche_delta_matrix.csv"))
delta_liver <- load_delta(file.path(T1_DIR, "liver_niche_delta_matrix.csv"))

log_msg(sprintf("T4 profiles: brain %dx%d, liver %dx%d",
                nrow(prof_brain), ncol(prof_brain), nrow(prof_liver), ncol(prof_liver)))

# ---------- Load microenv naming (v5 新增) ----------
load_microenv <- function(organ) {
  f <- file.path(MICROENV_DIR, paste0(organ, "_microenv_naming.csv"))
  d <- read.csv(f, stringsAsFactors = FALSE)
  list(
    domain_to_microenv = setNames(d$microenv_id, d$domain),
    microenv_to_label  = setNames(d$microenv_label, d$microenv_id) |> (\(x) x[!duplicated(names(x))])()
  )
}
me_brain <- load_microenv("brain")
me_liver <- load_microenv("liver")

log_msg(sprintf("Microenv loaded: brain %d microenv, liver %d microenv",
                length(me_brain$microenv_to_label),
                length(me_liver$microenv_to_label)))

# ---------- Color palettes (v4 全部 + v5 新增 microenv) ----------
ALL_CT_BRAIN <- c("Astrocyte_Gfap", "Astrocyte_Resting", "CPECs", "Endothelial", "Ependymal",
                  "Excitatory_L23", "Excitatory_L45", "Excitatory_L6",
                  "Hpc_CA1", "Hpc_CA3", "Hpc_DG",
                  "Inhibitory_Meis2", "Inhibitory_Pvalb", "Inhibitory_Sst", "Inhibitory_Vip",
                  "Microglia_Homeostatic", "Microglia_other", "Neurons",
                  "Oligo_Mature", "OPCs", "Pericytes", "VLMC")
ALL_CT_LIVER <- c("B_cells", "B_Plasma", "CD4_T", "CD8_T", "cDC1", "cDC2",
                  "Cholangiocytes", "Dendritic", "Endothelial", "Erythrocytes",
                  "Hep_Periportal", "Hep_Perivenous", "Hep_Proliferating",
                  "HSC_Activated_Chemo", "HSC_Activated_ECM", "HSC_Proliferating", "HSC_Quiescent",
                  "KC_Homeostatic", "Kupffer",
                  "LCM", "LSEC", "LSEC_Mid", "LSEC_Pericentral", "LSEC_Periportal",
                  "Macrophage_M1", "Macrophage_M2", "Macrophage_NAM",
                  "Monocytes", "Neutrophils", "NK_cells", "NKT",
                  "pDC", "Pericytes", "Stellate", "T_cells", "Treg")

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

# v5 新增: microenv 配色 (与 make_microenvironment_v2.R 保持一致)
LIVER_MICROENV_COLORS <- c(
  M01 = "#2C7FB8",  M02 = "#7FCDBB",  M03 = "#41B6C4",
  M04 = "#A1D99B",  M05 = "#E31A1C",  # ⭐ paper key niche
  M06 = "#FDAE6B",  M07 = "#9E9AC8"
)
BRAIN_MICROENV_COLORS <- c(
  M01 = "#E31A1C",  M02 = "#33A02C",  M03 = "#A6761D",
  M04 = "#1F78B4",  M05 = "#FF7F00",  M06 = "#984EA3"
)

# ---------- Helper: KNN smooth (v4) ----------
smooth_weights <- function(coords, weights, k = KNN_K) {
  nn <- RANN::nn2(coords[, c("x", "y")], coords[, c("x", "y")], k = k + 1)
  nn_idx <- nn$nn.idx
  out <- matrix(0, nrow = nrow(weights), ncol = ncol(weights),
                dimnames = dimnames(weights))
  for (j in seq_len(ncol(weights))) {
    out[, j] <- rowMeans(matrix(weights[nn_idx, j], nrow = nrow(weights)))
  }
  rs <- rowSums(out)
  out <- out / pmax(rs, 1e-9)
  out
}

# ---------- Helper: cosine similarity assign domain (v4) ----------
assign_domain <- function(weights_smoothed, profile_matrix) {
  common_ct <- intersect(colnames(weights_smoothed), colnames(profile_matrix))
  W <- weights_smoothed[, common_ct, drop = FALSE]
  P <- profile_matrix[, common_ct, drop = FALSE]
  W_norm <- sqrt(rowSums(W^2))
  P_norm <- sqrt(rowSums(P^2))
  dot    <- W %*% t(P)
  cos_sim <- dot / (pmax(W_norm, 1e-9) %o% pmax(P_norm, 1e-9))
  dom_idx <- max.col(cos_sim, ties.method = "first")
  rownames(profile_matrix)[dom_idx]
}

# ---------- Helper: TRUE hexagon binning + majority vote (v5 升级) ----------
# 用 hexbin 包构造真六边形 grid, 每 hex 内 majority vote, 输出 hex 顶点供 geom_polygon
hex_bin_majority <- function(coords, label_vec, n_bins = HEX_GRID_N,
                             size_x, size_y) {
  if (!requireNamespace("hexbin", quietly = TRUE)) {
    stop("Please install hexbin: install.packages('hexbin')")
  }
  
  # 用 hexbin 把 spot 聚到 hex
  # IDs=TRUE 让我们能查每 spot 落在哪个 hex
  hb <- hexbin::hexbin(coords$x, coords$y,
                       xbins = n_bins,
                       xbnds = c(0, size_x), ybnds = c(0, size_y),
                       IDs = TRUE)
  
  # 每 spot 的 hex ID
  spot_hex_id <- hb@cID
  
  # 在每 hex 内 majority vote
  bin_df <- data.frame(
    hex_id = spot_hex_id,
    label = label_vec,
    stringsAsFactors = FALSE
  )
  
  bin_majority <- bin_df %>%
    group_by(hex_id) %>%
    summarise(
      n_spots = n(),
      label = names(sort(table(label), decreasing = TRUE))[1],
      .groups = "drop"
    ) %>%
    filter(!is.na(label))
  
  # 获取每个 hex 的中心坐标
  hex_centers <- hexbin::hcell2xy(hb)
  hex_centers_df <- data.frame(
    hex_id = hb@cell,        # cell IDs (only nonzero hex)
    x_c = hex_centers$x,
    y_c = hex_centers$y
  )
  
  bin_majority <- bin_majority %>%
    left_join(hex_centers_df, by = "hex_id")
  
  # 构造每个 hex 的 6 个顶点 (用于 geom_polygon)
  # hexbin 单元的 横向间距 = (size_x / n_bins), 纵向间距 = 横向 * sqrt(3)/2
  dx <- (size_x - 0) / n_bins / 2  # half-width
  dy <- dx * 2 / sqrt(3)            # height of hex (vertical radius)
  
  # 6 个顶点 (顺时针, pointy-top hexagon)
  hex_offsets_x <- c(0, dx, dx, 0, -dx, -dx)
  hex_offsets_y <- c(dy, dy/2, -dy/2, -dy, -dy/2, dy/2)
  
  # 展开成 polygon 顶点 data.frame
  poly_list <- vector("list", nrow(bin_majority))
  for (i in seq_len(nrow(bin_majority))) {
    poly_list[[i]] <- data.frame(
      hex_id   = bin_majority$hex_id[i],
      vertex_id = 1:6,
      x = bin_majority$x_c[i] + hex_offsets_x,
      y = bin_majority$y_c[i] + hex_offsets_y,
      label = bin_majority$label[i],
      stringsAsFactors = FALSE
    )
  }
  poly_df <- do.call(rbind, poly_list)
  poly_df$group_id <- paste0("hex_", poly_df$hex_id)  # geom_polygon 需要 group
  
  list(
    bin_df = bin_majority,    # 每 hex 中心 + label (legacy)
    poly_df = poly_df,         # 每 hex 6 个顶点 (用于 geom_polygon)
    hex_dx = dx,
    hex_dy = dy,
    n_hex = nrow(bin_majority)
  )
}

# ---------- Validate CT colors ----------
missing_brain <- setdiff(ALL_CT_BRAIN, names(CT_COLORS_BRAIN))
missing_liver <- setdiff(ALL_CT_LIVER, names(CT_COLORS_LIVER))
if (length(missing_brain) > 0) log_msg(sprintf("WARN: brain CT missing color: %s", paste(missing_brain, collapse=",")))
if (length(missing_liver) > 0) log_msg(sprintf("WARN: liver CT missing color: %s", paste(missing_liver, collapse=",")))

# ============================================================
# Main loop
# ============================================================
for (s in names(SAMPLES)) {
  sample_t0 <- Sys.time()
  cfg <- SAMPLES[[s]]
  sample_out <- file.path(OUT_DIR, s)
  dir.create(sample_out, showWarnings = FALSE, recursive = TRUE)
  log_msg(sprintf("---- [%s] %s/%s ----", s, cfg$tissue, cfg$group))

  # Tissue-specific refs
  if (cfg$tissue == "brain") {
    ct_colors    <- CT_COLORS_BRAIN
    dom_profile  <- prof_brain
    dom_naming   <- name_brain
    delta_mat    <- delta_brain
    me_map       <- me_brain
    me_colors    <- BRAIN_MICROENV_COLORS
  } else {
    ct_colors    <- CT_COLORS_LIVER
    dom_profile  <- prof_liver
    dom_naming   <- name_liver
    delta_mat    <- delta_liver
    me_map       <- me_liver
    me_colors    <- LIVER_MICROENV_COLORS
  }

  # Load size
  if (!file.exists(cfg$official)) { log_msg("  SKIP: official missing"); next }
  seu <- readRDS(cfg$official)
  info <- seu@misc$info[[1]]
  size_x <- info$size_x; size_y <- info$size_y
  rm(seu); gc(verbose = FALSE)
  log_msg(sprintf("  size: %d x %d", size_x, size_y))

  # HE
  if (!file.exists(cfg$he)) { log_msg("  SKIP: HE missing"); next }
  he_img <- png::readPNG(cfg$he)
  he_grob <- rasterGrob(he_img, interpolate = FALSE,
                        width = unit(1, "npc"), height = unit(1, "npc"))
  rm(he_img); gc(verbose = FALSE)

  # RCTD
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

  # Per-spot top CT
  top_ct_idx <- max.col(weights, ties.method = "first")
  top_ct     <- colnames(weights)[top_ct_idx]

  # KNN-smoothed weights
  weights_smoothed <- smooth_weights(coords, weights, k = KNN_K)

  # Assign domain (cosine similarity)
  domain <- assign_domain(weights_smoothed, dom_profile)
  log_msg(sprintf("  Domain dist: %s",
                  paste(names(sort(table(domain), decreasing=TRUE)),
                        collapse=",")))
  
  # v5 新增: 把 domain 映射到 microenv
  microenv <- me_map$domain_to_microenv[domain]
  log_msg(sprintf("  Microenv dist: %s",
                  paste(names(sort(table(microenv), decreasing=TRUE)),
                        collapse=",")))

  # Common ggplot helpers
  plot_theme <- theme_void(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(color = "grey30"),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      plot.margin = margin(5, 5, 5, 5)
    )
  aspect <- size_y / size_x
  plot_height <- PLOT_WIDTH * aspect * 1.05

  # ===== Figure A: top CT (复用 v4) =====
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
    labs(title = sprintf("HE + Top CT (argmax) — %s (%s)", s, cfg$group),
         subtitle = sprintf("%d spots, CT colors consistent across samples", nrow(plot_df_A))) +
    plot_theme +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 2))

  ggsave(file.path(sample_out, "A_top_CT_on_HE.png"), pA,
         width = PLOT_WIDTH, height = plot_height, dpi = 150,
         limitsize = FALSE, bg = "white")
  log_msg("  A saved")
  rm(pA, plot_df_A); gc(verbose = FALSE)

  # ===== Figure B: 10 spatial domain (复用 v4) =====
  plot_df_B <- data.frame(x = coords$x, y = coords$y, domain = domain)
  all_doms <- sprintf("D%02d", 1:10)
  plot_df_B$domain <- factor(plot_df_B$domain, levels = all_doms)
  domain_legend <- dom_naming[all_doms]
  color_map_B <- setNames(DOMAIN_COLORS[all_doms], all_doms)

  pB <- ggplot(plot_df_B, aes(x = x, y = y, color = domain)) +
    annotation_custom(he_grob, xmin = 0, xmax = size_x, ymin = size_y, ymax = 0) +
    geom_point(size = POINT_SIZE, alpha = POINT_ALPHA) +
    scale_color_manual(values = color_map_B, labels = domain_legend,
                       name = "Spatial domain (10)", drop = FALSE) +
    coord_fixed(xlim = c(0, size_x), ylim = c(size_y, 0), expand = FALSE) +
    labs(title = sprintf("HE + Spatial Domain (10 domains) — %s (%s)", s, cfg$group),
         subtitle = sprintf("%d spots assigned via cosine similarity to T4 CT profiles",
                            nrow(plot_df_B))) +
    plot_theme +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1))

  ggsave(file.path(sample_out, "B_spatial_domain_on_HE.png"), pB,
         width = PLOT_WIDTH, height = plot_height, dpi = 150,
         limitsize = FALSE, bg = "white")
  log_msg("  B saved")
  rm(pB, plot_df_B); gc(verbose = FALSE)

  # ===== Figure C: tumor-induced score (复用 v4) =====
  nn <- RANN::nn2(coords[, c("x", "y")], coords[, c("x", "y")], k = KNN_K + 1)
  nn_idx <- nn$nn.idx[, -1, drop = FALSE]

  anchor_row_names <- rownames(delta_mat)
  ct_col_names     <- colnames(delta_mat)
  n_spots <- nrow(coords)
  score <- rep(NA_real_, n_spots)
  has_anchor <- top_ct %in% anchor_row_names
  neigh_ct <- matrix(top_ct[nn_idx], nrow = n_spots)

  for (i in seq_len(n_spots)) {
    if (!has_anchor[i]) next
    neigh <- neigh_ct[i, ]
    valid <- neigh %in% ct_col_names
    if (!any(valid)) { score[i] <- 0; next }
    neigh_v <- neigh[valid]
    d_vals <- delta_mat[top_ct[i], neigh_v, drop = TRUE]
    score[i] <- sum(d_vals) / KNN_K
  }

  valid_score <- score[!is.na(score)]
  if (length(valid_score) > 0) {
    lo <- quantile(valid_score, 0.05, na.rm = TRUE)
    hi <- quantile(valid_score, 0.95, na.rm = TRUE)
    lim_abs <- max(abs(lo), abs(hi))
    if (lim_abs < 0.1) lim_abs <- 0.1
  } else lim_abs <- 1
  score_clip <- pmin(pmax(score, -lim_abs), lim_abs)

  plot_df_C <- data.frame(x = coords$x, y = coords$y,
                          score = score_clip, has_score = !is.na(score))

  pC <- ggplot() +
    annotation_custom(he_grob, xmin = 0, xmax = size_x, ymin = size_y, ymax = 0) +
    geom_point(data = plot_df_C[!plot_df_C$has_score, ],
               aes(x = x, y = y), color = OTHER_COLOR,
               size = POINT_SIZE * 0.5, alpha = 0.3) +
    geom_point(data = plot_df_C[plot_df_C$has_score, ],
               aes(x = x, y = y, color = score),
               size = POINT_SIZE * 1.5, alpha = 1.0) +
    scale_color_gradient2(low = "#1b7837", mid = "#f7f7f7", high = "#b35806",
                          midpoint = 0, limits = c(-lim_abs, lim_abs),
                          name = "Tumor-induced\nscore") +
    coord_fixed(xlim = c(0, size_x), ylim = c(size_y, 0), expand = FALSE) +
    labs(title = sprintf("HE + Tumor-induced niche score — %s (%s)", s, cfg$group),
         subtitle = sprintf("KNN k=%d (orange = tumor-like, green = baseline-like)", KNN_K)) +
    plot_theme

  ggsave(file.path(sample_out, "C_tumor_induced_score_on_HE.png"), pC,
         width = PLOT_WIDTH, height = plot_height, dpi = 150,
         limitsize = FALSE, bg = "white")
  log_msg(sprintf("  C saved (score range: [%.2f, %.2f], %d/%d have score)",
                  min(valid_score, na.rm=TRUE), max(valid_score, na.rm=TRUE),
                  sum(!is.na(score)), n_spots))
  rm(pC, plot_df_C, nn, nn_idx, neigh_ct, score, score_clip, valid_score)
  gc(verbose = FALSE)

  # ===== Figure D: microenvironment (v5 新增) =====
  plot_df_D <- data.frame(x = coords$x, y = coords$y,
                          microenv_id = unname(microenv),
                          stringsAsFactors = FALSE)
  
  all_me <- sort(unique(plot_df_D$microenv_id))
  plot_df_D$microenv_id <- factor(plot_df_D$microenv_id, levels = all_me)
  
  # Legend: "M01 Hippocampus"
  me_legend_full <- paste0(all_me, " ",
                           me_map$microenv_to_label[all_me])
  names(me_legend_full) <- all_me
  
  color_map_D <- me_colors[all_me]
  
  n_me <- length(all_me)
  pD <- ggplot(plot_df_D, aes(x = x, y = y, color = microenv_id)) +
    annotation_custom(he_grob, xmin = 0, xmax = size_x, ymin = size_y, ymax = 0) +
    geom_point(size = POINT_SIZE, alpha = POINT_ALPHA) +
    scale_color_manual(values = color_map_D, labels = me_legend_full,
                       name = sprintf("Microenvironment\n(%d total)", n_me),
                       drop = FALSE) +
    coord_fixed(xlim = c(0, size_x), ylim = c(size_y, 0), expand = FALSE) +
    labs(title = sprintf("HE + Microenvironment (%d) — %s (%s)",
                         n_me, s, cfg$group),
         subtitle = sprintf("%d spots aggregated from 10 spatial domains via biological merge",
                            nrow(plot_df_D))) +
    plot_theme +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1))

  ggsave(file.path(sample_out, "D_microenvironment_on_HE.png"), pD,
         width = PLOT_WIDTH, height = plot_height, dpi = 150,
         limitsize = FALSE, bg = "white")
  log_msg("  D saved")
  rm(pD); gc(verbose = FALSE)

  # ===== Figure E: TRUE hex-binned microenv majority (v5 升级) =====
  hex_res <- hex_bin_majority(coords, microenv, n_bins = HEX_GRID_N,
                              size_x = size_x, size_y = size_y)
  poly_df <- hex_res$poly_df
  poly_df$microenv_id <- factor(poly_df$label, levels = all_me)
  
  # 用 geom_polygon 画真六边形
  pE <- ggplot(poly_df,
               aes(x = x, y = y, group = group_id, fill = microenv_id)) +
    annotation_custom(he_grob, xmin = 0, xmax = size_x, ymin = size_y, ymax = 0) +
    geom_polygon(alpha = 0.78, color = NA) +
    scale_fill_manual(values = color_map_D, labels = me_legend_full,
                      name = sprintf("Microenvironment\n(%d total)", n_me),
                      drop = FALSE) +
    coord_fixed(xlim = c(0, size_x), ylim = c(size_y, 0), expand = FALSE) +
    labs(title = sprintf("HE + Hexagon-binned microenvironment — %s (%s)",
                         s, cfg$group),
         subtitle = sprintf("%d xbins (true hexagons), %d hex bins, majority vote",
                            HEX_GRID_N, hex_res$n_hex)) +
    plot_theme +
    guides(fill = guide_legend(override.aes = list(alpha = 1), ncol = 1))

  ggsave(file.path(sample_out, "E_hex_microenv_on_HE.png"), pE,
         width = PLOT_WIDTH, height = plot_height, dpi = 150,
         limitsize = FALSE, bg = "white")
  log_msg(sprintf("  E saved (%d hex bins, %d polygon vertices)",
                  hex_res$n_hex, nrow(poly_df)))
  rm(pE, poly_df, hex_res, plot_df_D, coords, weights, weights_smoothed, he_grob)
  gc(verbose = FALSE)

  elapsed <- as.numeric(difftime(Sys.time(), sample_t0, units = "secs"))
  log_msg(sprintf("  [%s] done in %.1fs", s, elapsed))
}

log_msg("=== All done ===")
log_msg(sprintf("Output: %s", OUT_DIR))
