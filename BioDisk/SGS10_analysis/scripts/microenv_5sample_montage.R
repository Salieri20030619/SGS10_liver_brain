#!/usr/bin/env Rscript
# ============================================================================
# microenv_5sample_montage.R   [restyle, replot-only]
# docx img8(liver) / img18(brain): 每器官 5 样本微环境空间图横排(1x5) + 右侧共用图例.
#   - 读已存的 spot->microenv 分配表 (taskB/spatial_domain_RCTD100/<organ>_spot_microenv_assignment.csv)
#   - 坐标按 barcode 去 rctd_full/per_sample/<sample>_rctd_results.rds 取 (不重算微环境)
#   - 无 HE 背景; y 轴翻转与 HE/解剖方向一致
#   - 先不画 ★ niche 标记 (需要可加回)
#   输出: figures/restyle_v1/<organ>_microenv_5sample.{pdf,tiff,png}
# ============================================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})
source("/mnt/BioDisk/wangmingyuan/SGS10_analysis/scripts/fig_save.R")

PR       <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
SD_DIR   <- file.path(PR, "taskB", "spatial_domain_RCTD100")
RCTD_DIR <- file.path(PR, "spatial_results", "rctd_full", "per_sample")
OUT_DIR  <- file.path(PR, "figures", "restyle_v1")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- 可调参数 ----
POINT_SIZE <- 0.40
FIG_W      <- 22     # 整张宽(英寸); 太挤/太空可调
FIG_H      <- 3.2    # 整张高; 上下留白多就调小, 标题被切就调大

# ---- 微环境配色 (沿用 fig16_19) ----
LIVER_MICROENV_COLORS <- c(M01="#2C7FB8", M02="#7FCDBB", M03="#41B6C4",
  M04="#A1D99B", M05="#E31A1C", M06="#FDAE6B", M07="#9E9AC8",
  M08="#cab2d6", M09="#fb9a99")
BRAIN_MICROENV_COLORS <- c(M01="#E31A1C", M02="#33A02C", M03="#A6761D",
  M04="#1F78B4", M05="#FF7F00", M06="#984EA3")

# ---- 样本顺序 (control, blank, 3x tumor) ----
SAMPLES <- list(
  liver = c("Con_liver_1", "Bla_liver_1", "Liver_3", "Liver_8", "gan56"),
  brain = c("Con_brain_1", "Bla_brain_1", "Brain_3", "Brain_8", "nao56")
)

# Bla_brain_1 原始坐标朝向已正确, 其余4个脑样本需转180度
ROTATE_BRAIN_SAMPLES <- c("Con_brain_1", "Brain_3", "Brain_8", "nao56")

# ---- 取某样本坐标 (barcode, x, y) ----
load_coords <- function(sample) {
  f <- file.path(RCTD_DIR, paste0(sample, "_rctd_results.rds"))
  if (!file.exists(f)) stop("缺 rctd_results: ", f)
  rctd <- readRDS(f)
  co <- as.data.frame(rctd@spatialRNA@coords)
  co$barcode <- rownames(co)
  rm(rctd); gc(verbose = FALSE)
  co[, c("barcode", "x", "y")]
}

# ---- 单样本面板 ----
build_panel <- function(sample, assign_df, me_colors, me_levels, me_labels, me_title) {
  sub <- assign_df[assign_df$sample == sample, ]
  sub$barcode <- sub("^.*__", "", sub$spot_id)
  co <- load_coords(sample)
  df <- merge(co, sub[, c("barcode", "microenv_id")], by = "barcode")
  if (sample %in% ROTATE_BRAIN_SAMPLES) {  # [client] 脑样本坐标绕各样本中心镜像 -> 整体转180度
    df$x <- max(df$x) + min(df$x) - df$x
    df$y <- max(df$y) + min(df$y) - df$y
  }
  cat(sprintf("  [%s] %d spots (assign %d, coords %d)\n",
              sample, nrow(df), nrow(sub), nrow(co)))
  df$microenv_id <- factor(df$microenv_id, levels = me_levels)

  ggplot(df, aes(x = x, y = y, color = microenv_id)) +
    geom_point(size = POINT_SIZE, alpha = 0.9) +
    scale_color_manual(values = me_colors, labels = me_labels,
                       name = me_title, drop = FALSE) +
    coord_fixed(xlim = range(df$x), ylim = rev(range(df$y)), expand = FALSE) +
    labs(title = sample) +
    theme_void(base_size = 13) +
    theme(
      plot.title    = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.title  = element_text(size = 12, face = "bold"),
      legend.text   = element_text(size = 11),
      plot.margin   = margin(4, 4, 4, 4)
    ) +
    guides(color = guide_legend(override.aes = list(size = 4, alpha = 1), ncol = 1))
}

# ---- 一个器官的 1x5 拼图 ----
make_montage <- function(organ) {
  cat(sprintf("\n==== %s ====\n", organ))
  assign_df <- read.csv(file.path(SD_DIR, paste0(organ, "_spot_microenv_assignment.csv")),
                        stringsAsFactors = FALSE)
  me_colors_all <- if (organ == "liver") LIVER_MICROENV_COLORS else BRAIN_MICROENV_COLORS

  present   <- sort(unique(assign_df$microenv_id))
  lab_map   <- assign_df %>% distinct(microenv_id, microenv_name)
  me_labels <- setNames(lab_map$microenv_name, lab_map$microenv_id)[present]
  me_colors <- me_colors_all[present]
  me_title  <- sprintf("Microenvironment\n(%d total)", length(present))

  panels <- lapply(SAMPLES[[organ]], build_panel,
                   assign_df = assign_df, me_colors = me_colors,
                   me_levels = present, me_labels = me_labels, me_title = me_title)

  combined <- wrap_plots(panels, nrow = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "right")

  save_fig_all(combined,
               file   = file.path(OUT_DIR, sprintf("%s_microenv_5sample", organ)),
               width  = FIG_W, height = FIG_H, dpi = 300)
  cat(sprintf("  [OK] %s_microenv_5sample -> restyle_v1\n", organ))
}

make_montage("liver")
make_montage("brain")

cat("\nDone.\n")
