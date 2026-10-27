#!/usr/bin/env Rscript
# ============================================================================
# compose_figures_1_5.R  (v5: full layout, fig03/04 restored from PDF)
# Panel labels via ggplot PNG (magick freetype broken on this server).
# ============================================================================
suppressPackageStartupMessages({
  library(magick)
  library(ggplot2)
})

PNG_DIR  <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis/figures/paper_v2/PNG"
OUT_DIR  <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis/figures/paper_v2/COMPOSED"
TMP_DIR  <- tempdir()
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

make_label_png <- function(letter, out_path, w_px = 100, h_px = 100) {
  p <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = letter,
             fontface = "bold", size = 35) +
    scale_x_continuous(limits = c(0, 1), expand = c(0,0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
    theme_void() +
    theme(plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent", color = NA))
  ggsave(out_path, p, width = w_px/100, height = h_px/100, dpi = 100,
         bg = "transparent")
  invisible(out_path)
}

LABEL_CACHE <- list()
for (lt in LETTERS[1:8]) {
  pp <- file.path(TMP_DIR, sprintf("label_%s.png", lt))
  make_label_png(lt, pp, w_px = 100, h_px = 100)
  LABEL_CACHE[[lt]] <- pp
}
cat(sprintf("[label] cached %d label PNGs\n", length(LABEL_CACHE)))

add_label <- function(img, letter, pad_x = 15, pad_y = 15) {
  label_img <- image_read(LABEL_CACHE[[letter]])
  image_composite(img, label_img, operator = "over",
                  offset = sprintf("+%d+%d", pad_x, pad_y))
}

load_panel <- function(name, target_w = 2000) {
  path <- file.path(PNG_DIR, paste0(name, ".png"))
  if (!file.exists(path)) stop("Missing: ", path)
  img <- image_read(path)
  info <- image_info(img)
  cat(sprintf("  loaded %s  (%d x %d)\n", name, info$width, info$height))
  if (info$width != target_w) img <- image_scale(img, sprintf("%dx", target_w))
  img
}

pad_to <- function(img, target_w, target_h) {
  info <- image_info(img)
  if (info$width == target_w && info$height == target_h) return(img)
  image_extent(img, sprintf("%dx%d", target_w, target_h),
               color = "white", gravity = "center")
}
row_combine <- function(imgs) image_append(image_join(imgs), stack = FALSE)
col_combine <- function(imgs) image_append(image_join(imgs), stack = TRUE)

save_composed <- function(img, basename) {
  out_png <- file.path(OUT_DIR, paste0(basename, ".png"))
  image_write(img, out_png, format = "png", quality = 95)
  info <- image_info(img)
  cat(sprintf("[OK] %s  %d x %d\n", basename, info$width, info$height))
}

# ============================================================================
# Figure 1: Cell atlas (A liver UMAP, B brain UMAP)
# ============================================================================
cat("\n=== Figure 1: Cell atlas ===\n")
W <- 2000
f1A <- add_label(load_panel("fig20_liver_umap", W), "A")
f1B <- add_label(load_panel("fig21_brain_umap", W), "B")
h_max <- max(image_info(f1A)$height, image_info(f1B)$height)
f1A <- pad_to(f1A, W, h_max); f1B <- pad_to(f1B, W, h_max)
save_composed(row_combine(list(f1A, f1B)), "Figure_1_cell_atlas")

# ============================================================================
# Figure 2: Liver host response (2x2)
#   A volcano  B HSC pseudobulk lollipop  C HSC->Kupffer hub  D silent inflam
# ============================================================================
cat("\n=== Figure 2: Liver host response ===\n")
f2A <- add_label(load_panel("fig01_HSC_volcano", W), "A")
f2B <- add_label(load_panel("fig04_HSC_DEG_lollipop_pb", W), "B")
f2C <- add_label(load_panel("fig05_HSC_Kupffer_hub", W), "C")
f2D <- add_label(load_panel("fig10_silent_inflammation", W), "D")

h_top <- max(image_info(f2A)$height, image_info(f2B)$height)
f2A <- pad_to(f2A, W, h_top); f2B <- pad_to(f2B, W, h_top)
row1 <- row_combine(list(f2A, f2B))

h_bot <- max(image_info(f2C)$height, image_info(f2D)$height)
f2C <- pad_to(f2C, W, h_bot); f2D <- pad_to(f2D, W, h_bot)
row2 <- row_combine(list(f2C, f2D))

w_total <- max(image_info(row1)$width, image_info(row2)$width)
row1 <- pad_to(row1, w_total, image_info(row1)$height)
row2 <- pad_to(row2, w_total, image_info(row2)$height)
save_composed(col_combine(list(row1, row2)), "Figure_2_liver_host_response")

# ============================================================================
# Figure 3: Brain dual-mode TF paradigm
#   A TF dual mode (wide top, full row)
#   B TF Tumor vs Control  C TF Tumor vs Blank  D TGFb pair
# ============================================================================
cat("\n=== Figure 3: Brain dual-mode TF paradigm ===\n")
f3A <- add_label(load_panel("fig09_TF_dual_mode", W * 2), "A")
f3B <- add_label(load_panel("fig03_Excit_L45_TF_lollipop", W), "B")
f3C <- add_label(load_panel("fig08_Excit_L45_TF_tvsb", W), "C")
f3D <- add_label(load_panel("fig11_X1_TGFb_pair", W), "D")

h_bot3 <- max(image_info(f3B)$height, image_info(f3C)$height,
              image_info(f3D)$height)
f3B <- pad_to(f3B, W, h_bot3)
f3C <- pad_to(f3C, W, h_bot3)
f3D <- pad_to(f3D, W, h_bot3)
row_bot <- row_combine(list(f3B, f3C, f3D))

w_total <- max(image_info(f3A)$width, image_info(row_bot)$width)
f3A     <- pad_to(f3A, w_total, image_info(f3A)$height)
row_bot <- pad_to(row_bot, w_total, image_info(row_bot)$height)
save_composed(col_combine(list(f3A, row_bot)), "Figure_3_brain_dual_mode_TF")

# ============================================================================
# Figure 4: Cross-organ axes (A thyroid cascade, B efferocytosis 3 layers)
# ============================================================================
cat("\n=== Figure 4: Cross-organ axes ===\n")
f4A <- add_label(load_panel("fig06_thyroid_cascade", W), "A")
f4B <- add_label(load_panel("fig07_efferocytosis_3layers", W), "B")
h_max <- max(image_info(f4A)$height, image_info(f4B)$height)
f4A <- pad_to(f4A, W, h_max); f4B <- pad_to(f4B, W, h_max)
save_composed(row_combine(list(f4A, f4B)), "Figure_4_cross_organ_axes")

# ============================================================================
# Figure 5: Niche-level integrative (2x3)
# ============================================================================
cat("\n=== Figure 5: Niche-level integrative ===\n")
W_5 <- 1800
f5A <- add_label(load_panel("fig12_liver_microenv_heatmap", W_5), "A")
f5B <- add_label(load_panel("fig13_brain_microenv_heatmap", W_5), "B")
f5C <- add_label(load_panel("fig16_Liver_3_microenv_HE",  W_5 * 2), "C")
f5D <- add_label(load_panel("fig17_Brain_3_microenv_HE",  W_5 * 2), "D")
f5E <- add_label(load_panel("fig14_liver_cellchat_chord", W_5 * 2), "E")
f5F <- add_label(load_panel("fig15_brain_cellchat_chord", W_5 * 2), "F")

h_r1 <- max(image_info(f5A)$height, image_info(f5B)$height)
f5A <- pad_to(f5A, W_5, h_r1); f5B <- pad_to(f5B, W_5, h_r1)
row1_5 <- row_combine(list(f5A, f5B))

f5C_r <- image_scale(f5C, sprintf("%dx", W_5))
f5D_r <- image_scale(f5D, sprintf("%dx", W_5))
h_r2 <- max(image_info(f5C_r)$height, image_info(f5D_r)$height)
f5C_r <- pad_to(f5C_r, W_5, h_r2); f5D_r <- pad_to(f5D_r, W_5, h_r2)
row2_5 <- row_combine(list(f5C_r, f5D_r))

f5E_r <- image_scale(f5E, sprintf("%dx", W_5))
f5F_r <- image_scale(f5F, sprintf("%dx", W_5))
h_r3 <- max(image_info(f5E_r)$height, image_info(f5F_r)$height)
f5E_r <- pad_to(f5E_r, W_5, h_r3); f5F_r <- pad_to(f5F_r, W_5, h_r3)
row3_5 <- row_combine(list(f5E_r, f5F_r))

w_total <- max(image_info(row1_5)$width,
               image_info(row2_5)$width,
               image_info(row3_5)$width)
row1_5 <- pad_to(row1_5, w_total, image_info(row1_5)$height)
row2_5 <- pad_to(row2_5, w_total, image_info(row2_5)$height)
row3_5 <- pad_to(row3_5, w_total, image_info(row3_5)$height)
save_composed(col_combine(list(row1_5, row2_5, row3_5)),
              "Figure_5_niche_integrative")

cat("\n[FINISHED] 5 figures composed in:\n  ", OUT_DIR, "\n")
