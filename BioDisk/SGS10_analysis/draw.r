# ============================================================
# Microenv 堆叠柱状图 v14 (final)
# 客户最终需求:
#   1. 配色 = NEW10 (用户提供的范例 10 色) + paper-key M05/M04 = #E66D50 橙红
#   2. 纵坐标 = "Proportion"
#   3. 删 "Sample" x 轴标题
#   4. 样本名 30° 斜排, 右下角对齐柱子右边缘
#   5. baseline/tumor 分组线 + 标签放在样本名下方
# 用 patchwork 上下拼接, 保证 xtick 与柱子严格对齐
# ============================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(patchwork)
})

DATA_DIR <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis/taskB/spatial_domain_RCTD100"
OUT_DIR  <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis/figures/microenv_stacked"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- 范例 NEW10 (从用户提供的范例图采样) ----
NEW10 <- c(
  D01 = "#6BBAA2", D02 = "#D09C77", D03 = "#B99998",
  D04 = "#AA9BC2", D05 = "#D999B3", D06 = "#ACCF67",
  D07 = "#DED549", D08 = "#ECD361", D09 = "#D6BDA1",
  D10 = "#B2B3B7"
)
KEY_COLOR <- "#E66D50"   # paper-key niche 用的橙红 (黄昏绮景)

# ---- liver 7 微环境配色 ----
liver_cols <- c(
  "M01 Hepatocyte_zonation"     = NEW10[["D01"]],
  "M02 Vascular_perisinusoidal" = NEW10[["D03"]],
  "M03 Sinusoidal_endothelium"  = NEW10[["D04"]],
  "M04 Stellate_quiescent"      = NEW10[["D06"]],
  "M05 Tumor_immune_niche"      = KEY_COLOR,     # ★ paper-key
  "M06 B_plasma"                = NEW10[["D08"]],
  "M07 Kupffer"                 = NEW10[["D09"]]
)
# ---- brain 6 微环境配色 ----
brain_cols <- c(
  "M01 Hippocampus"        = NEW10[["D01"]],
  "M02 Cortex_excitatory"  = NEW10[["D03"]],
  "M03 Cortex_inhibitory"  = NEW10[["D04"]],
  "M04 Vascular_BBB"       = KEY_COLOR,          # ★ paper-key
  "M05 Glia"               = NEW10[["D08"]],
  "M06 General_neurons"    = NEW10[["D09"]]
)

# ---- 分组工具 ----
sample_to_group <- function(s) {
  ifelse(grepl("^Con_", s), "control",
  ifelse(grepl("^Bla_", s), "blank", "tumor"))
}
sample_to_macro <- function(s) {
  ifelse(grepl("^(Con_|Bla_)", s), "baseline", "tumor")
}

liver_sample_order <- c("Con_liver_1","Bla_liver_1","Liver_3","Liver_8","gan56")
brain_sample_order <- c("Con_brain_1","Bla_brain_1","Brain_3","Brain_8","nao56")
group_order <- c("control","blank","tumor")

BAR_W <- 0.62   # 柱宽 (用于计算样本名右移量)

# ---- 主图 (per_sample 用, 不显示 x text 和 title, 留给下方分组条) ----
plot_main_persample <- function(df, fill_levels, cols, x_levels) {
  df$sample <- factor(df$sample, levels = x_levels)
  df$microenv_name <- factor(df$microenv_name, levels = fill_levels)
  
  agg <- df |>
    group_by(sample, microenv_name) |>
    summarise(n = n(), .groups = "drop") |>
    group_by(sample) |>
    mutate(prop = n / sum(n)) |>
    ungroup()
  
  ggplot(agg, aes(x = sample, y = prop, fill = microenv_name)) +
    geom_col(width = BAR_W, color = NA) +
    scale_fill_manual(values = cols, name = NULL, drop = FALSE) +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       expand = c(0, 0),
                       breaks = c(0, 0.25, 0.50, 0.75, 1.0),
                       limits = c(0, 1)) +
    labs(y = "Proportion") +
    theme_bw(base_size = 12) +
    theme(
      panel.grid       = element_blank(),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.8),
      axis.ticks.y     = element_line(color = "black", linewidth = 0.4),
      axis.ticks.x     = element_blank(),
      axis.ticks.length = unit(3, "pt"),
      axis.text.x      = element_blank(),
      axis.text.y      = element_text(color = "black", size = 9),
      axis.title.x     = element_blank(),
      axis.title.y     = element_text(size = 11),
      legend.key.size  = unit(11, "pt"),
      legend.text      = element_text(size = 9),
      plot.margin      = margin(8, 6, 0, 6)   # 底部 margin=0, 紧贴下方分组条
    )
}

# ---- 下方分组条: 30° 斜排样本名 (右下角对齐柱子右边缘) + 分组线 + baseline/tumor 文字 ----
plot_grouptrack <- function(x_levels, macro_groups) {
  macros <- macro_groups[x_levels]
  grp_runs <- rle(macros)
  pos <- cumsum(c(0, grp_runs$lengths))
  
  seg_df <- data.frame(
    label  = grp_runs$values,
    xstart = pos[-length(pos)] + 1 - 0.45,
    xend   = pos[-1] + 0.45,
    xmid   = (pos[-length(pos)] + 1 + pos[-1]) / 2
  )
  
  # 样本名: 用 numeric x_right (柱子中心 + 半柱宽 = 右边缘), 不用 factor x
  # 这样名字的右下角对齐到柱子右边缘, 而不是中心
  label_df <- data.frame(
    x_right = seq_along(x_levels) + BAR_W / 2,
    label   = x_levels
  )
  
  # 用 factor dummy 锚定 x scale (与主图一致), 这样 patchwork 对齐时不偏
  dummy <- data.frame(x = factor(x_levels, levels = x_levels), y = 0)
  
  # 布局 (y 范围 0~1):
  #   y = 0.92: 样本名顶部 (30° 斜排, hjust=1 右对齐到柱子右边缘)
  #   y = 0.32: 分组线
  #   y = 0.10: baseline / tumor 文字
  ggplot(dummy, aes(x = x, y = y)) +
    geom_point(alpha = 0, show.legend = FALSE) +
    # 样本名 (numeric x, 偏移到柱子右边缘)
    geom_text(data = label_df, inherit.aes = FALSE,
              aes(x = x_right, y = 0.92, label = label),
              angle = 30, hjust = 1, vjust = 1,
              size = 3.6, color = "black") +
    # 分组线
    geom_segment(data = seg_df, inherit.aes = FALSE,
                 aes(x = xstart, xend = xend, y = 0.32, yend = 0.32),
                 linewidth = 0.7, color = "black") +
    # baseline / tumor 文字
    geom_text(data = seg_df, inherit.aes = FALSE,
              aes(x = xmid, y = 0.10, label = label),
              size = 4.2, fontface = "bold", color = "black",
              hjust = 0.5, vjust = 0.5) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme_void(base_size = 12) +
    theme(plot.margin = margin(0, 6, 4, 6))
}

# ---- per_group 单图 (control / blank / tumor, 不需要分组条) ----
plot_per_group <- function(df, fill_levels, cols, x_levels) {
  df$group <- factor(df$group, levels = x_levels)
  df$microenv_name <- factor(df$microenv_name, levels = fill_levels)
  
  agg <- df |>
    group_by(group, microenv_name) |>
    summarise(n = n(), .groups = "drop") |>
    group_by(group) |>
    mutate(prop = n / sum(n)) |>
    ungroup()
  
  ggplot(agg, aes(x = group, y = prop, fill = microenv_name)) +
    geom_col(width = BAR_W, color = NA) +
    scale_fill_manual(values = cols, name = NULL, drop = FALSE) +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       expand = c(0, 0),
                       breaks = c(0, 0.25, 0.50, 0.75, 1.0),
                       limits = c(0, 1)) +
    labs(y = "Proportion", x = NULL) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid       = element_blank(),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.8),
      axis.ticks       = element_line(color = "black", linewidth = 0.4),
      axis.ticks.length = unit(3, "pt"),
      axis.text.x      = element_text(color = "black", size = 10,
                                      angle = 30, hjust = 1),
      axis.text.y      = element_text(color = "black", size = 9),
      axis.title.y     = element_text(size = 11),
      legend.key.size  = unit(11, "pt"),
      legend.text      = element_text(size = 9),
      plot.margin      = margin(8, 6, 8, 6)
    )
}

# ---- 主流程 ----
process_organ <- function(organ, cols, sample_order) {
  csv <- file.path(DATA_DIR, sprintf("%s_spot_microenv_assignment.csv", organ))
  cat("Reading", csv, "...\n")
  df <- read.csv(csv, stringsAsFactors = FALSE)
  df$group <- sample_to_group(df$sample)
  
  fill_lvls <- names(cols)
  macro_map <- setNames(sample_to_macro(sample_order), sample_order)
  
  # --- per_sample: patchwork 上下拼接 ---
  p_main <- plot_main_persample(df, fill_lvls, cols, sample_order)
  p_grp  <- plot_grouptrack(sample_order, macro_map)
  
  p_combined <- p_main / p_grp +
    plot_layout(heights = c(10, 2.8), guides = "collect") &
    theme(legend.position = "right")
  
  ggsave(file.path(OUT_DIR, sprintf("%s_per_sample.pdf", organ)), p_combined,
         width = 7.5, height = 6.0, device = cairo_pdf)
  ggsave(file.path(OUT_DIR, sprintf("%s_per_sample.png", organ)), p_combined,
         width = 7.5, height = 6.0, dpi = 300)
  
  # --- per_group: 单图 ---
  p2 <- plot_per_group(df, fill_lvls, cols, group_order)
  ggsave(file.path(OUT_DIR, sprintf("%s_per_group.pdf", organ)), p2,
         width = 5.4, height = 4.8, device = cairo_pdf)
  ggsave(file.path(OUT_DIR, sprintf("%s_per_group.png", organ)), p2,
         width = 5.4, height = 4.8, dpi = 300)
  
  cat("Saved 4 files for", organ, "\n\n")
}

process_organ("liver", liver_cols, liver_sample_order)
process_organ("brain", brain_cols, brain_sample_order)

cat("=== ALL DONE ===\n")
cat("Output dir:", OUT_DIR, "\n")
