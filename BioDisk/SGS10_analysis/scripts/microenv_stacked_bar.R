#!/usr/bin/env Rscript
# ============================================================================
# microenv_stacked_bar.R   [restyle, fresh]  v2: 客户统一配色 + patchwork 分组条
# docx img9 (liver) 微环境构成堆叠柱 + 对称 brain 版.
#   每器官: 5 样本 x 微环境占比(0-100%); 主图全封闭外框;
#           下方独立 strip 画 baseline/tumor 横线+标签(patchwork 对齐, 不跑偏).
#   配色: 单一客户调色板(按 M 序号), liver/brain 共用.
#   输出: figures/restyle_v1/<organ>_microenv_stacked.{pdf,tiff,png}
# ============================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(scales); library(patchwork)
})
source("/mnt/BioDisk/wangmingyuan/SGS10_analysis/scripts/fig_save.R")

PR      <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
SD_DIR  <- file.path(PR, "taskB", "spatial_domain_RCTD100")
OUT_DIR <- file.path(PR, "figures", "restyle_v1")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- 客户统一配色 (从 docx img9 取色; liver/brain 按 M 序号共用) ----
MICROENV_COLORS <- c(
  M01 = "#68B8A0", M02 = "#B89898", M03 = "#A898C0", M04 = "#A8C860",
  M05 = "#E06850", M06 = "#E8D060", M07 = "#D0B8A0"
)

# 样本顺序: 前2个=baseline, 其余=tumor
SAMPLES <- list(
  liver = c("Con_liver_1", "Bla_liver_1", "Liver_3", "Liver_8", "gan56"),
  brain = c("Con_brain_1", "Bla_brain_1", "Brain_3", "Brain_8", "nao56")
)

make_stacked <- function(organ) {
  cat(sprintf("\n==== %s ====\n", organ))
  df <- read.csv(file.path(SD_DIR, paste0(organ, "_spot_microenv_assignment.csv")),
                 stringsAsFactors = FALSE)
  samples <- SAMPLES[[organ]]
  df <- df[df$sample %in% samples, ]

  present   <- sort(unique(df$microenv_id))
  lab_map   <- df %>% distinct(microenv_id, microenv_name)
  me_labels <- setNames(lab_map$microenv_name, lab_map$microenv_id)[present]
  me_colors <- MICROENV_COLORS[present]

  prop <- df %>%
    count(sample, microenv_id) %>%
    group_by(sample) %>% mutate(frac = n / sum(n)) %>% ungroup()
  prop$sample      <- factor(prop$sample, levels = samples)
  prop$microenv_id <- factor(prop$microenv_id, levels = present)

  n_base <- 2; n_all <- length(samples)

  # ---- 主图: 堆叠柱 + 全封闭外框, x 轴样本名(45°) ----
  p <- ggplot(prop, aes(x = sample, y = frac, fill = microenv_id)) +
    geom_col(width = 0.8) +
    scale_fill_manual(values = me_colors, labels = me_labels,
                      name = "Microenvironment") +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       expand = expansion(mult = c(0, 0.02))) +
    scale_x_discrete(expand = expansion(add = 0.5)) +
    labs(x = NULL, y = "Proportion") +
    theme_classic(base_size = 18) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      axis.line    = element_blank(),
      axis.text.x  = element_text(angle = 45, hjust = 1, size = 16),
      axis.ticks.x = element_line(linewidth = 0.4),
      axis.title.y = element_text(size = 19),
      axis.text.y  = element_text(size = 16),
      legend.title = element_text(size = 17, face = "bold"),
      legend.text  = element_text(size = 16),
      plot.margin  = margin(6, 6, 2, 6)
    )

  # ---- 下方分组条: baseline/tumor 横线 + 标签 (与主图 x 对齐) ----
  seg_df <- data.frame(
    x   = c(0.6, n_base + 0.6),
    xend= c(n_base + 0.4, n_all + 0.4),
    xc  = c((1 + n_base) / 2, (n_base + 1 + n_all) / 2),
    lab = c("baseline", "tumor")
  )
  pb <- ggplot(seg_df) +
    geom_segment(aes(x = x, xend = xend, y = 0.78, yend = 0.78), linewidth = 0.7) +
    geom_text(aes(x = xc, y = 0.30, label = lab), fontface = "bold", size = 6.5) +
    scale_x_continuous(limits = c(0.5, n_all + 0.5), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    theme_void() +
    theme(plot.margin = margin(0, 6, 4, 6))

  final <- p / pb + plot_layout(heights = c(11, 1), guides = "collect") &
    theme(legend.position = "right")

  save_fig_all(final, file = file.path(OUT_DIR, sprintf("%s_microenv_stacked", organ)),
               width = 7.5, height = 7.5, dpi = 300)
  cat(sprintf("  [OK] %s_microenv_stacked -> restyle_v1\n", organ))
}

make_stacked("liver")
make_stacked("brain")

cat("\nDone.\n")
