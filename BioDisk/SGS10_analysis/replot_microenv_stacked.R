# ============================================================
# Microenv 堆叠柱状图 v11 (final, 稳定版)
# - 配色: NEW10 + paper-key M05/M04 = #E66D50 橙红
# - 纵坐标 = "Proportion"
# - 关键修复: 用 patchwork 拼接, 分组条独立 ggplot, 不再用 y<0 负坐标 hack
#   保证 baseline/tumor 分组线/文字 / 样本名 / Sample 标题 三层完全分离
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
  D01 = "#6BBAA2",  D02 = "#D09C77",  D03 = "#B99998",
  D04 = "#AA9BC2",  D05 = "#D999B3",  D06 = "#ACCF67",
  D07 = "#DED549",  D08 = "#ECD361",  D09 = "#D6BDA1",
  D10 = "#B2B3B7"
)
KEY_COLOR <- "#E66D50"

liver_cols <- c(
  "M01 Hepatocyte_zonation"     = NEW10[["D01"]],
  "M02 Vascular_perisinusoidal" = NEW10[["D03"]],
  "M03 Sinusoidal_endothelium"  = NEW10[["D04"]],
  "M04 Stellate_quiescent"      = NEW10[["D06"]],
  "M05 Tumor_immune_niche"      = KEY_COLOR,
  "M06 B_plasma"                = NEW10[["D08"]],
  "M07 Kupffer"                 = NEW10[["D09"]]
)
brain_cols <- c(
  "M01 Hippocampus"        = NEW10[["D01"]],
  "M02 Cortex_excitatory"  = NEW10[["D03"]],
  "M03 Cortex_inhibitory"  = NEW10[["D04"]],
  "M04 Vascular_BBB"       = KEY_COLOR,
  "M05 Glia"               = NEW10[["D08"]],
  "M06 General_neurons"    = NEW10[["D09"]]
)

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

# ---- 主柱状图 (无 xtext 版, 留给下方分组条) ----
plot_main <- function(df, x_col, fill_levels, cols, x_levels, show_xtext = TRUE) {
  df[[x_col]] <- factor(df[[x_col]], levels = x_levels)
  df$microenv_name <- factor(df$microenv_name, levels = fill_levels)
  
  agg <- df |>
    group_by(.data[[x_col]], microenv_name) |>
    summarise(n = n(), .groups = "drop") |>
    group_by(.data[[x_col]]) |>
    mutate(prop = n / sum(n)) |>
    ungroup()
  
  p <- ggplot(agg, aes(x = .data[[x_col]], y = prop, fill = microenv_name)) +
    geom_col(width = 0.62, color = NA) +
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
      axis.ticks       = element_line(color = "black", linewidth = 0.4),
      axis.ticks.length = unit(3, "pt"),
      axis.text.y      = element_text(color = "black", size = 9),
      axis.title.y     = element_text(size = 11),
      legend.key.size  = unit(11, "pt"),
      legend.text      = element_text(size = 9),
      plot.margin      = margin(8, 6, 4, 6)
    )
  
  if (show_xtext) {
    # per_group 图: 直接显示 group 标签, 不需要分组条
    p <- p + labs(x = "Group") +
      theme(axis.text.x = element_text(color = "black", size = 10,
                                       angle = 30, hjust = 1),
            axis.title.x = element_text(size = 11))
  } else {
    # per_sample 图: 上半部分主图, xtext 隐藏 (移到下方分组条)
    p <- p + theme(axis.title.x = element_blank(),
                   axis.text.x  = element_blank(),
                   axis.ticks.x = element_blank())
  }
  
  p
}

# ---- 下方分组条: 单独 ggplot, 上半画 baseline/tumor 线+文字, 下半放样本名 ----
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
  
  # 关键: x 轴用与主图相同的 factor (确保 patchwork align 时 bar 中心对齐)
  dummy <- data.frame(x = factor(x_levels, levels = x_levels), y = 0)
  
  ggplot(dummy, aes(x = x, y = y)) +
    # 占位 geom (透明) 把 x scale 锚定到 factor
    geom_point(alpha = 0, show.legend = FALSE) +
    # 在 y=0.85 画分组线, y=0.50 放文字
    geom_segment(data = seg_df, inherit.aes = FALSE,
                 aes(x = xstart, xend = xend, y = 0.85, yend = 0.85),
                 linewidth = 0.7, color = "black") +
    geom_text(data = seg_df, inherit.aes = FALSE,
              aes(x = xmid, y = 0.45, label = label),
              size = 4.0, fontface = "bold", color = "black", vjust = 0.5) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(x = "Sample") +
    theme_void(base_size = 12) +
    theme(
      axis.text.x  = element_text(color = "black", size = 10,
                                  angle = 30, hjust = 1,
                                  margin = margin(t = 2)),
      axis.title.x = element_text(size = 11, margin = margin(t = 8)),
      plot.margin  = margin(0, 6, 8, 6)
    )
}

process_organ <- function(organ, cols, sample_order) {
  csv <- file.path(DATA_DIR, sprintf("%s_spot_microenv_assignment.csv", organ))
  cat("Reading", csv, "...\n")
  df <- read.csv(csv, stringsAsFactors = FALSE)
  df$group <- sample_to_group(df$sample)
  
  fill_lvls <- names(cols)
  macro_map <- setNames(sample_to_macro(sample_order), sample_order)
  
  # --- per_sample: patchwork 上下拼接 ---
  p_main <- plot_main(df, "sample", fill_lvls, cols, sample_order, show_xtext = FALSE)
  p_grp  <- plot_grouptrack(sample_order, macro_map)
  
  # 拼接: heights = c(14, 2.2) 比例; collect 图例放右侧
  p_combined <- p_main / p_grp +
    plot_layout(heights = c(14, 2.2), guides = "collect") &
    theme(legend.position = "right")
  
  ggsave(file.path(OUT_DIR, sprintf("%s_per_sample.pdf", organ)), p_combined,
         width = 7.5, height = 5.8, device = cairo_pdf)
  ggsave(file.path(OUT_DIR, sprintf("%s_per_sample.png", organ)), p_combined,
         width = 7.5, height = 5.8, dpi = 300)
  
  # --- per_group: 单一 ggplot, 不需要分组条 ---
  p2 <- plot_main(df, "group", fill_lvls, cols, group_order, show_xtext = TRUE)
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
