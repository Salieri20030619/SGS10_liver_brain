#!/usr/bin/env Rscript
source("/mnt/BioDisk/wangmingyuan/SGS10_analysis/scripts/fig_save.R")
# ============================================================
# microenv_marker_genes_RCTD100_replot.R   (纯重绘, 不重算)
# 只读已存盘的 combined_seurat.rds + markers_all.csv,
# 只跑 DotPlot / DoHeatmap; 跳过 merge / SCTransform /
# PrepSCTFindMarkers / FindAllMarkers (这些是跑一天的部分).
#   输入(原脚本已生成):
#     D1_v2_RCTD_composition/microenv_markers_RCTD100/
#       <organ>_combined_seurat.rds
#       <organ>_microenv_markers_all.csv
#   输出: figures/restyle_v1/  (PDF + TIFF + PNG@300)
# ============================================================
options(future.globals.maxSize = 50 * 1024^3)
suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(ggplot2); library(tibble)
  library(ComplexHeatmap); library(circlize)
})

PR          <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
IN_DIR      <- file.path(PR, "D1_v2_RCTD_composition/microenv_markers_RCTD100")
RESTYLE_DIR <- file.path(PR, "figures/restyle_v1")
dir.create(RESTYLE_DIR, showWarnings = FALSE, recursive = TRUE)
set.seed(42)

ORGANS <- c("liver", "brain")

replot_organ <- function(organ) {
  cat(sprintf("\n##### %s #####\n", toupper(organ)))
  rds_path <- file.path(IN_DIR, paste0(organ, "_combined_seurat.rds"))
  mk_path  <- file.path(IN_DIR, paste0(organ, "_microenv_markers_all.csv"))
  if (!file.exists(rds_path) || !file.exists(mk_path)) {
    cat(sprintf("WARN 缺少 rds/csv, 跳过 %s\n  %s\n  %s\n", organ, rds_path, mk_path))
    return(invisible())
  }
  cat("  载入 combined_seurat.rds (已SCT, 不重算) ...\n")
  combined <- readRDS(rds_path)
  markers  <- read.csv(mk_path, stringsAsFactors = FALSE)
  markers$cluster <- as.character(markers$cluster)

  # ---- DotPlot (每微环境 top5) ----
  top5_genes <- markers %>% group_by(cluster) %>%
    slice_max(avg_log2FC, n = 5, with_ties = FALSE) %>% pull(gene) %>% unique()
  p_dot <- DotPlot(combined, features = top5_genes,
                   group.by = "microenv_id", scale = TRUE) +
    RotatedAxis() +
    scale_color_gradient2(low = "navy", mid = "white", high = "firebrick3", midpoint = 0) +
    theme(axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 13))
  .dw <- max(10, length(top5_genes) * 0.18 + 2)
  .dh <- max(5, length(unique(markers$cluster)) * 0.5 + 2)
  save_fig_all(p_dot,
               file.path(RESTYLE_DIR, paste0(organ, "_microenv_markers_dotplot")),
               width = .dw, height = .dh)

  # ---- DoHeatmap (每微环境 top10, 每组抽样200 spot) ----
  top10_genes <- markers %>% group_by(cluster) %>%
    slice_max(avg_log2FC, n = 5, with_ties = FALSE) %>% pull(gene) %>% unique()
  set.seed(42)
  cells_subset <- combined@meta.data %>%
    rownames_to_column("cell") %>%
    group_by(microenv_id) %>% slice_sample(n = 200) %>% pull(cell)
  combined_sub <- subset(combined, cells = cells_subset)
  combined_sub <- ScaleData(combined_sub, features = top10_genes, verbose = FALSE)

  # ===== [client] ComplexHeatmap 复刻 + anno_mark 连线标注 =====
  .mat <- as.matrix(suppressWarnings(SeuratObject::LayerData(combined_sub, layer = "scale.data")))
  .feat <- top10_genes[top10_genes %in% rownames(.mat)]
  .mat <- .mat[.feat, , drop = FALSE]
  .me  <- factor(as.character(combined_sub$microenv_id[colnames(.mat)]))
  .me  <- factor(.me, levels = sort(levels(.me)))            # M01..M0N
  .ord <- order(.me); .mat <- .mat[, .ord, drop = FALSE]; .me <- .me[.ord]
  .mat[.mat >  2.5] <-  2.5; .mat[.mat < -2.5] <- -2.5
  .pal <- Seurat::PurpleAndYellow()
  .col_fun <- circlize::colorRamp2(seq(min(.mat), max(.mat), length.out = length(.pal)), .pal)

  # 重点基因: liver=fig12 key / brain=微环境M04 marker
  .top1  <- markers %>% dplyr::group_by(cluster) %>%
              dplyr::slice_max(avg_log2FC, n = 1, with_ties = FALSE) %>% dplyr::pull(gene)
  .top1c <- gsub("^mm10---|^GRCh38-", "", .top1)
  if (organ == "brain") {
    .keyc <- gsub("^mm10---|^GRCh38-", "",
                  markers$gene[as.character(markers$cluster) %in% c("M04","4","M4","04")])
  } else {
    .keyc <- c("Cd74","H2-Aa","Psap","Inpp5d","Ctsb","Zeb2","Runx1","Laptm5")
  }
  .rnc <- gsub("^mm10---|^GRCh38-", "", rownames(.mat))
  .is_key <- .rnc %in% .keyc; .is_top1 <- .rnc %in% .top1c
  .idx <- which(.is_key | .is_top1)
  .labs <- rownames(.mat)[.idx]
  .lcol <- ifelse(.is_key[.idx], "#C0143B", "black")
  cat(sprintf("[%s] mark: red=%d black=%d (of %d rows)\n",
              organ, sum(.is_key), sum(.is_top1 & !.is_key), nrow(.mat)))

  .melv   <- levels(.me)
  .barcol <- setNames(scales::hue_pal()(length(.melv)), .melv)
  .topann <- ComplexHeatmap::HeatmapAnnotation(
               Microenv = .me, col = list(Microenv = .barcol),
               show_legend = FALSE, show_annotation_name = FALSE,
               simple_anno_size = grid::unit(4, "mm"))
  .rowann <- ComplexHeatmap::rowAnnotation(
               mark = ComplexHeatmap::anno_mark(
                 at = .idx, labels = .labs, which = "row",
                 labels_gp = grid::gpar(col = .lcol, fontsize = 20),
                 link_width = grid::unit(10, "mm")))
  .ht <- ComplexHeatmap::Heatmap(
           .mat, name = "scaled", col = .col_fun,
           cluster_rows = FALSE, cluster_columns = FALSE,
           show_row_names = FALSE, show_column_names = FALSE,
           column_split = .me, column_gap = grid::unit(1.5, "mm"),
           column_title_gp = grid::gpar(fontsize = 18),
           top_annotation = .topann, right_annotation = .rowann,
           use_raster = TRUE, raster_quality = 4,
           show_heatmap_legend = FALSE, border = TRUE)

  .hmw <- 14; .hmh <- 10
  .base <- file.path(RESTYLE_DIR, paste0(organ, "_microenv_markers_heatmap"))
  grDevices::cairo_pdf(paste0(.base, ".pdf"), width = .hmw, height = .hmh); ComplexHeatmap::draw(.ht); grDevices::dev.off()
  grDevices::tiff(paste0(.base, ".tiff"), width = .hmw, height = .hmh, units = "in", res = 300, compression = "lzw"); ComplexHeatmap::draw(.ht); grDevices::dev.off()
  grDevices::png(paste0(.base, ".png"),  width = .hmw, height = .hmh, units = "in", res = 300); ComplexHeatmap::draw(.ht); grDevices::dev.off()
  rm(combined, combined_sub); gc(verbose = FALSE)
  cat(sprintf("  %s done -> dotplot + heatmap\n", organ))
}

for (o in ORGANS) replot_organ(o)
cat("\nAll done. 图已写入 figures/restyle_v1/\n")
