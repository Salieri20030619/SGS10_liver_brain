#!/usr/bin/env Rscript
# ============================================================================
# fig16_19_spatial_overlay.R
# H&E + microenv overlay, paper-grade (Batch 6).
# Based on he_overlay_RCTD100.R panel D, with paper changes:
#   - English title removed (no inline title/subtitle, paper caption only)
#   - paper_theme fonts; legend right side
#   - Triple-format output to /figures/paper_v2/{PDF,TIFF,PNG}
#   - Only 4 representative samples:
#       fig16 Liver_3      (tumor liver,  M05 Tumor_immune_niche)
#       fig17 Brain_3      (tumor brain,  M04 Vascular_BBB)
#       fig18 Con_liver_1  (control liver, no M05)
#       fig19 Con_brain_1  (control brain, normal microenv)
# ============================================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(png)
  library(grid)
  library(dplyr)
  library(RANN)
})
source("/mnt/BioDisk/wangmingyuan/SGS10_analysis/scripts/paper_theme.R")

PR <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
BD <- "/mnt/BioData/bio/SGS10例"
SD_DIR <- file.path(PR, "taskB/spatial_domain_RCTD100")
RCTD_DIR <- file.path(PR, "spatial_results/rctd_full/per_sample")
OUT_DIR <- file.path(PR, "figures/paper_v2")
for (sub in c("PDF","TIFF","PNG")) dir.create(file.path(OUT_DIR, sub),
                                              showWarnings = FALSE, recursive = TRUE)

KNN_K       <- 15
POINT_SIZE  <- 0.65
POINT_ALPHA <- 0.92

# ---------- Sample registry ----------
MATRIX_LIVER <- file.path(BD, "SGS25090559-8例样本数据/data/matrix")
RDS_LIVER    <- file.path(BD, "SGS25090559-8例样本数据/data/rds")
SAMPLES <- list(
  Liver_3 = list(
    tissue = "liver", group = "tumor", fig = "fig16",
    official = file.path(RDS_LIVER, "25090559_Liver_3_expression/Liver_3.rds"),
    he       = file.path(MATRIX_LIVER, "25090559_Liver_3_expression/Liver_3_aligned_HE_TIMG.png")
  ),
  Brain_3 = list(
    tissue = "brain", group = "tumor", fig = "fig17",
    official = file.path(RDS_LIVER, "25090559_Brain_3_expression/Brain_3.rds"),
    he       = file.path(MATRIX_LIVER, "25090559_Brain_3_expression/Brain_3_aligned_HE_TIMG.png")
  ),
  Con_liver_1 = list(
    tissue = "liver", group = "control", fig = "fig18",
    official = file.path(RDS_LIVER, "25090559_Con_liver_1_expression/Con_liver_1.rds"),
    he       = file.path(MATRIX_LIVER, "25090559_Con_liver_1_expression/Con_liver_1_aligned_HE_TIMG.png")
  ),
  Con_brain_1 = list(
    tissue = "brain", group = "control", fig = "fig19",
    official = file.path(RDS_LIVER, "25090559_Con_brain_1_expression/Con_brain_1.rds"),
    he       = file.path(MATRIX_LIVER, "25090559_Con_brain_1_expression/Con_brain_1_aligned_HE_TIMG.png")
  )
)

# ---------- Load profiles / naming / microenv (copied from original) ----------
load_domain_profile <- function(f) {
  d <- read.csv(f, check.names = FALSE)
  rownames(d) <- d[[1]]; d[[1]] <- NULL
  as.matrix(d)
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
me_brain   <- load_microenv("brain")
me_liver   <- load_microenv("liver")

# ---------- Microenv colors (copied) ----------
LIVER_MICROENV_COLORS <- c(
  M01 = "#2C7FB8",  M02 = "#7FCDBB",  M03 = "#41B6C4",
  M04 = "#A1D99B",  M05 = "#E31A1C",   # paper-key
  M06 = "#FDAE6B",  M07 = "#9E9AC8",
  M08 = "#cab2d6",  M09 = "#fb9a99"
)
BRAIN_MICROENV_COLORS <- c(
  M01 = "#E31A1C",  M02 = "#33A02C",  M03 = "#A6761D",
  M04 = "#1F78B4",  M05 = "#FF7F00",  M06 = "#984EA3"
)

# ---------- Helpers (copied) ----------
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

# ---------- Paper-grade theme for spatial overlay ----------
plot_theme_paper <- function() {
  theme_void(base_size = 12) +
    theme(
      plot.title    = element_blank(),
      plot.subtitle = element_blank(),
      legend.position = "right",
      legend.title    = element_text(size = 11, face = "bold", color = "black"),
      legend.text     = element_text(size = 9.5, color = "black"),
      legend.key      = element_blank(),
      legend.background = element_blank(),
      plot.margin = margin(5, 5, 5, 5)
    )
}

# ---------- Main render function ----------
render_sample <- function(sname, cfg) {
  cat(sprintf("\n[%s] %s/%s  (%s)\n", cfg$fig, cfg$tissue, cfg$group, sname))
  t0 <- Sys.time()

  if (!file.exists(cfg$official)) stop("Missing official RDS: ", cfg$official)
  if (!file.exists(cfg$he))       stop("Missing HE PNG:      ", cfg$he)

  # 1) Seurat RDS -> size_x/y from misc info
  seu <- readRDS(cfg$official)
  info <- seu@misc$info[[1]]
  size_x <- info$size_x; size_y <- info$size_y
  rm(seu); gc(verbose = FALSE)
  cat(sprintf("  size: %d x %d\n", size_x, size_y))

  # 2) HE PNG -> rasterGrob
  he_img <- png::readPNG(cfg$he)
  he_grob <- rasterGrob(he_img, interpolate = FALSE,
                        width = unit(1, "npc"), height = unit(1, "npc"))
  rm(he_img); gc(verbose = FALSE)

  # 3) RCTD weights + coords
  rctd_path <- file.path(RCTD_DIR, paste0(sname, "_rctd_results.rds"))
  if (!file.exists(rctd_path)) stop("Missing RCTD: ", rctd_path)
  rctd <- readRDS(rctd_path)
  coords  <- rctd@spatialRNA@coords
  weights <- as.matrix(rctd@results$weights)
  rs <- rowSums(weights); weights <- weights / pmax(rs, 1e-9)
  rm(rctd); gc(verbose = FALSE)
  common <- intersect(rownames(coords), rownames(weights))
  coords  <- coords[common, , drop = FALSE]
  weights <- weights[common, , drop = FALSE]
  cat(sprintf("  RCTD: %d spots, %d CT\n", nrow(coords), ncol(weights)))

  # 4) Domain + microenv inference
  if (cfg$tissue == "brain") {
    dom_profile <- prof_brain
    me_map      <- me_brain
    me_colors   <- BRAIN_MICROENV_COLORS
  } else {
    dom_profile <- prof_liver
    me_map      <- me_liver
    me_colors   <- LIVER_MICROENV_COLORS
  }

  weights_smoothed <- smooth_weights(coords, weights, k = KNN_K)
  domain   <- assign_domain(weights_smoothed, dom_profile)
  microenv <- me_map$domain_to_microenv[domain]
  cat(sprintf("  Microenv dist: %s\n",
              paste(sprintf("%s=%d", names(table(microenv)),
                            as.integer(table(microenv))), collapse=", ")))

  # 5) Build plot dataframe
  plot_df <- data.frame(x = coords$x, y = coords$y,
                        microenv_id = unname(microenv),
                        stringsAsFactors = FALSE)
  all_me <- sort(unique(plot_df$microenv_id))
  plot_df$microenv_id <- factor(plot_df$microenv_id, levels = all_me)
  me_legend <- paste0(all_me, " - ", me_map$microenv_to_label[all_me])
  names(me_legend) <- all_me
  color_map <- me_colors[all_me]

  aspect <- size_y / size_x
  PLOT_W <- 10
  plot_h <- PLOT_W * aspect * 1.10  # 10% padding for legend

  # 6) Render plot
  p <- ggplot(plot_df, aes(x = x, y = y, color = microenv_id)) +
    annotation_custom(he_grob, xmin = 0, xmax = size_x,
                      ymin = size_y, ymax = 0) +
    geom_point(size = POINT_SIZE, alpha = POINT_ALPHA) +
    scale_color_manual(values = color_map, labels = me_legend,
                       name = "Microenvironment", drop = FALSE) +
    coord_fixed(xlim = c(0, size_x), ylim = c(size_y, 0), expand = FALSE) +
    plot_theme_paper() +
    guides(color = guide_legend(
      override.aes = list(size = 3.5, alpha = 1),
      ncol = 1, byrow = TRUE
    ))

  # 7) Triple-format save
  basename <- sprintf("%s_%s_microenv_HE", cfg$fig, sname)
  pdf_path <- file.path(OUT_DIR, "PDF",  paste0(basename, ".pdf"))
  tif_path <- file.path(OUT_DIR, "TIFF", paste0(basename, ".tif"))
  png_path <- file.path(OUT_DIR, "PNG",  paste0(basename, ".png"))

  ggsave(pdf_path, p, width = PLOT_W, height = plot_h,
         device = cairo_pdf, family = "Helvetica",
         limitsize = FALSE)
  ggsave(tif_path, p, width = PLOT_W, height = plot_h,
         dpi = 300, compression = "lzw", device = "tiff",
         limitsize = FALSE)
  ggsave(png_path, p, width = PLOT_W, height = plot_h,
         dpi = 200, device = "png", bg = "white",
         limitsize = FALSE)

  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("[OK] %s  %.0f x %.0f in  ->  PDF + TIFF + PNG  (%.0fs)\n",
              basename, PLOT_W, plot_h, dt))

  rm(p, plot_df, coords, weights, weights_smoothed, he_grob)
  gc(verbose = FALSE)
}

# ============================================================================
# Run all 4 samples
# ============================================================================
cat("##############################################\n")
cat("# fig16-19 Spatial overlay paper-grade\n")
cat(sprintf("# Started: %s\n", Sys.time()))
cat("##############################################\n")

for (sname in names(SAMPLES)) {
  tryCatch(render_sample(sname, SAMPLES[[sname]]),
           error = function(e) {
             cat(sprintf("[ERROR] %s: %s\n", sname, conditionMessage(e)))
           })
}

cat(sprintf("\n[FINISHED] %s\n", Sys.time()))
