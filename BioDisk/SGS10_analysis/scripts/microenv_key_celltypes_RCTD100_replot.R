#!/usr/bin/env Rscript
source("/mnt/BioDisk/wangmingyuan/SGS10_analysis/scripts/fig_save.R")
# ============================================================
# microenv_key_celltypes_RCTD100_replot.R   (纯重绘, 不重算)
# 直接读已存盘的汇总 csv 画图, 完全不碰 RCTD 权重.
#   输入(原脚本已生成):
#     D1_v2_RCTD_composition/microenv_key_celltypes_RCTD100/
#       <organ>_microenv_celltypes_full.csv   -> 热图
#       <organ>_microenv_celltypes_top5.csv   -> 气泡
#   输出: figures/restyle_v1/  (PDF + TIFF + PNG@300)
# ============================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
  library(pheatmap); library(RColorBrewer); library(tibble)
})

PR          <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
IN_DIR      <- file.path(PR, "D1_v2_RCTD_composition/microenv_key_celltypes_RCTD100")
RESTYLE_DIR <- file.path(PR, "figures/restyle_v1")
dir.create(RESTYLE_DIR, showWarnings = FALSE, recursive = TRUE)

ORGANS <- c("liver", "brain")

replot_organ <- function(organ) {
  cat(sprintf("\n##### %s #####\n", toupper(organ)))
  full_csv <- file.path(IN_DIR, paste0(organ, "_microenv_celltypes_full.csv"))
  top5_csv <- file.path(IN_DIR, paste0(organ, "_microenv_celltypes_top5.csv"))
  if (!file.exists(full_csv) || !file.exists(top5_csv)) {
    cat(sprintf("WARN 缺少 csv, 跳过 %s\n  %s\n  %s\n", organ, full_csv, top5_csv))
    return(invisible())
  }
  me_ct <- read.csv(full_csv, stringsAsFactors = FALSE, check.names = FALSE)
  top5  <- read.csv(top5_csv, stringsAsFactors = FALSE, check.names = FALSE)

  lab_df      <- dplyr::distinct(me_ct, microenv_id, microenv_label)
  me_to_label <- setNames(lab_df$microenv_label, lab_df$microenv_id)

  # ---- 热图 (microenv × top25 cell type) ----
  heat_mat <- me_ct %>%
    pivot_wider(id_cols = microenv_id, names_from = celltype, values_from = mean_weight) %>%
    column_to_rownames("microenv_id") %>% as.matrix()
  heat_mat <- heat_mat[sort(rownames(heat_mat)), , drop = FALSE]
  col_max  <- apply(heat_mat, 2, max)
  top_cts  <- names(sort(col_max, decreasing = TRUE))[1:min(25, length(col_max))]
  heat_top <- heat_mat[, top_cts, drop = FALSE]
  rownames(heat_top) <- paste0(rownames(heat_top), " ", me_to_label[rownames(heat_top)])

  .hm_h <- max(4, nrow(heat_top) * 0.45 + 2)
  .ph <- pheatmap(heat_top,
           color = colorRampPalette(c("white","yellow","red","darkred"))(50),
           cluster_rows = FALSE, cluster_cols = FALSE,
           main = NA, silent = TRUE,
           fontsize_row = 12, fontsize_col = 11,
           display_numbers = TRUE, number_format = "%.2f", number_color = "black",
           fontsize_number = 9)
  save_fig_all(.ph$gtable,
               file.path(RESTYLE_DIR, sprintf("%s_microenv_celltype_heatmap", organ)),
               width = 12, height = .hm_h)

  # ---- 气泡 (top5 per microenv) ----
  top5_for_dot <- top5 %>% mutate(microenv_full = paste0(microenv_id, " ", microenv_label))
  p_dot <- ggplot(top5_for_dot,
                  aes(x = celltype, y = microenv_full,
                      size = n_spots, color = mean_weight)) +
    geom_point() +
    scale_color_gradient(low = "lightyellow", high = "red", name = "Mean weight") +
    scale_size_continuous(range = c(3, 12), name = "n spots") +
    labs(x = "Cell type", y = "Microenvironment") +
    theme_classic(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          axis.text.y = element_text(size = 12),
          legend.title = element_text(size = 13),
          legend.text  = element_text(size = 11))
  .dot_w <- max(8, length(unique(top5_for_dot$celltype)) * 0.3 + 4)
  .dot_h <- max(5, length(unique(top5_for_dot$microenv_full)) * 0.5 + 2)
  save_fig_all(p_dot,
               file.path(RESTYLE_DIR, sprintf("%s_microenv_celltype_dotplot", organ)),
               width = .dot_w, height = .dot_h)
  cat(sprintf("  %s done -> heatmap + dotplot\n", organ))
}

for (o in ORGANS) replot_organ(o)
cat("\nAll done. 图已写入 figures/restyle_v1/\n")
