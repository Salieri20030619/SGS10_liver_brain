#!/usr/bin/env Rscript
# ============================================================
# compose_liver_microenv.R  (紧凑版; B/C 等高填满, 角标留窄白栏)
#   A=HE叠加  B=marker热图  C=占比堆叠柱  D=组成热图
# 输出: figures/restyle_v1/COMPOSED/Figure_liver_microenv_composed.{pdf,tiff,png}
# ============================================================
suppressPackageStartupMessages({ library(magick); library(cowplot); library(ggplot2) })

DIR <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis/figures/restyle_v1"

# ---- 可调参数 ----
W      <- 2600    # 整图宽(px)
pad    <- 22      # 行间留白(px)
GL     <- 110     # 左侧角标白栏宽(px) -- A/B/D 角标
GC     <- 110     # B 与 C 之间的缝(px) -- C 角标
lab_sz <- 12      # 角标字号(pt)
outer  <- 6       # 外边白(px)

cat("step1: 读图 + 裁白边\n")
FZ <- 2
rd <- function(f) image_trim(image_read(file.path(DIR, f)), fuzz = FZ)
A  <- rd("Liver_3/D_microenvironment_on_HE.png")
B  <- rd("liver_microenv_markers_heatmap.png")
Cc <- rd("liver_microenv_stacked.png")
D  <- rd("liver_microenv_celltype_heatmap.png")

cat("step2: 堆叠(整宽行, 左留窄白栏; B/C 等高填满)\n")
Wc <- W - GL
mkfull <- function(img) {                       # 内容靠右, 左补 GL 白 -> 宽 W
  r <- image_resize(img, paste0(Wc, "x"))
  image_extent(r, paste0(W, "x", image_info(r)$height), gravity = "east", color = "white")
}
rA <- mkfull(A)
rD <- mkfull(D)

# B/C: 按各自宽高比反算, 使两者等高且 wB+GC+wC = Wc (底部不留空)
aB <- image_info(B)$width  / image_info(B)$height
aC <- image_info(Cc)$width / image_info(Cc)$height
avail <- Wc - GC
h2 <- round(avail / (aB + aC))
Bs <- image_resize(B,  paste0("x", h2))
Cs <- image_resize(Cc, paste0("x", h2))
wB <- image_info(Bs)$width
rBC_c <- image_append(c(Bs, image_blank(GC, h2, "white"), Cs))
rBC   <- image_extent(rBC_c, paste0(W, "x", h2), gravity = "east", color = "white")

fig <- image_append(c(rA, image_blank(W, pad, "white"), rBC,
                      image_blank(W, pad, "white"), rD), stack = TRUE)
if (outer > 0) fig <- image_border(fig, "white", paste0(outer, "x", outer))

cat("step3: 写临时 PNG\n")
tmp <- tempfile(fileext = ".png")
image_write(fig, tmp, format = "png", density = 300)

cat("step4: cowplot 叠角标\n")
HA <- image_info(rA)$height; HD <- image_info(rD)$height
total_h <- HA + pad + h2 + pad + HD + 2 * outer
total_w <- W + 2 * outer
xL <- (outer + 14)           / total_w
xC <- (outer + GL + wB + 20) / total_w
yA <- 1 - (outer + 10)                       / total_h
yB <- 1 - (outer + HA + pad + 10)            / total_h
yD <- 1 - (outer + HA + pad + h2 + pad + 10) / total_h
g <- ggdraw() +
  draw_image(tmp, x = 0, y = 0, width = 1, height = 1) +
  draw_label("A", x = xL, y = yA, fontface = "bold", size = lab_sz, hjust = 0, vjust = 1) +
  draw_label("B", x = xL, y = yB, fontface = "bold", size = lab_sz, hjust = 0, vjust = 1) +
  draw_label("C", x = xC, y = yB, fontface = "bold", size = lab_sz, hjust = 0, vjust = 1) +
  draw_label("D", x = xL, y = yD, fontface = "bold", size = lab_sz, hjust = 0, vjust = 1) +
  theme(plot.background = element_rect(fill = "white", colour = NA))

cat("step5: 存盘\n")
OUT  <- file.path(DIR, "COMPOSED"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
base <- file.path(OUT, "Figure_liver_microenv_composed")
Win <- total_w / 300; Hin <- total_h / 300
sav <- function(ext, ...) {
  f <- paste0(base, ext)
  ok <- tryCatch({ ggsave(f, g, width = Win, height = Hin, dpi = 300, limitsize = FALSE, ...); TRUE },
                 error = function(e) { cat("  [FAIL]", ext, "->", conditionMessage(e), "\n"); FALSE })
  if (ok) cat("  [OK]  ", f, "\n")
}
sav(".pdf", device = cairo_pdf)
sav(".png")
sav(".tiff", compression = "lzw")
cat(sprintf("DONE  %.0f x %.0f px  (%.1f x %.1f in)  | B高=C高=%d px\n", total_w, total_h, Win, Hin, h2))
