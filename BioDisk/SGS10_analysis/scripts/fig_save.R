#!/usr/bin/env Rscript
# ============================================================
# fig_save.R  —  统一出图器 (PDF矢量 + TIFF 300dpi LZW + PNG 300dpi)
# 目的: 所有脚本走同一个保存入口, 保证拼图时分辨率/格式一致
# 用法:
#   source("/mnt/BioDisk/wangmingyuan/SGS10_analysis/scripts/fig_save.R")
#   # ggplot / patchwork:
#   save_fig_all(p, "/path/fig01", width = 10, height = 8)
#   # pheatmap 等 grob/gtable:
#   save_fig_all(ph$gtable, "/path/heatmap", width = 12, height = 6)
#   # base 图形 (CellChat netVisual_* 等), 用 drawfun:
#   save_fig_all(drawfun = function() {
#       par(mfrow = c(1,2)); netVisual_circle(...); netVisual_circle(...)
#   }, file = "/path/cellchat", width = 14, height = 7)
# 说明: width/height 单位均为英寸; TIFF/PNG 用 dpi 控分辨率, PDF 为矢量
# ============================================================

save_fig_all <- function(obj = NULL, file, width, height,
                         dpi = 300, drawfun = NULL,
                         formats = c("pdf", "tiff", "png"),
                         tiff_compression = "lzw", bg = "white") {

  # 统一的"绘制动作": ggplot/patchwork -> print; grob/gtable -> grid.draw; 或自定义 drawfun
  .draw <- function() {
    if (!is.null(drawfun)) {
      drawfun()
    } else if (inherits(obj, "ggplot") || inherits(obj, "patchwork")) {
      print(obj)
    } else if (inherits(obj, c("gtable", "grob", "gTree"))) {
      grid::grid.newpage(); grid::grid.draw(obj)
    } else {
      stop("save_fig_all: obj 类型未知, 请改用 drawfun=")
    }
  }

  if ("pdf" %in% formats) {
    if (capabilities("cairo")) {
      grDevices::cairo_pdf(paste0(file, ".pdf"), width = width, height = height,
                           family = "Helvetica", bg = bg)
    } else {
      grDevices::pdf(paste0(file, ".pdf"), width = width, height = height, bg = bg)
    }
    .draw(); grDevices::dev.off()
  }
  if ("tiff" %in% formats) {
    grDevices::tiff(paste0(file, ".tiff"), width = width, height = height,
                    units = "in", res = dpi, compression = tiff_compression, bg = bg)
    .draw(); grDevices::dev.off()
  }
  if ("png" %in% formats) {
    grDevices::png(paste0(file, ".png"), width = width, height = height,
                   units = "in", res = dpi, bg = bg)
    .draw(); grDevices::dev.off()
  }
  cat(sprintf("[save_fig_all] %s -> %s\n", basename(file),
              paste(formats, collapse = "/")))
  invisible(file)
}
