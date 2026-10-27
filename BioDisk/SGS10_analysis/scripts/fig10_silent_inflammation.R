#!/usr/bin/env Rscript
# ============================================================================
# fig10_silent_inflammation.R  (v2: tighter y-axis, broken effect)
# Silent inflammation cross-modality negative finding.
# 4 cytokines (Il1b/Tnf/Il6/Ifng) max pct vs sepsis reference line.
# Uses ylim trick: zoom to 0-2% panel, with sepsis ref. as inset annotation.
# ============================================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})
source("/mnt/BioDisk/wangmingyuan/SGS10_analysis/scripts/paper_theme.R")

df <- tibble::tribble(
  ~cytokine, ~metric,                ~pct,
  "Il1b",    "Max pct (10 spatial)", 0.72,
  "Il1b",    "Final sc report",      0.57,
  "Tnf",     "Max pct (10 spatial)", 0.68,
  "Tnf",     "Final sc report",      0.26,
  "Il6",     "Max pct (10 spatial)", 0.08,
  "Il6",     "Final sc report",      0.02,
  "Ifng",    "Max pct (10 spatial)", 0.06,
  "Ifng",    "Final sc report",      0.06
)
df$cytokine <- factor(df$cytokine, levels = c("Il1b", "Tnf", "Il6", "Ifng"))
df$metric   <- factor(df$metric,
                      levels = c("Max pct (10 spatial)", "Final sc report"))
cat("[fig10]\n"); print(df)

p <- ggplot(df, aes(x = cytokine, y = pct, fill = metric)) +
  geom_col(position = position_dodge(0.75), width = 0.65, color = NA) +
  geom_text(aes(label = sprintf("%.2f", pct)),
            position = position_dodge(0.75),
            vjust = -0.5, size = 2.8, color = "grey25") +
  # Sepsis reference as annotation (out of plot range, shown as text only)
  annotate("text", x = 2.5, y = 1.55,
           label = "Sepsis ref. TNF 5-20% (Hotchkiss 2016)",
           color = "#C0143B", size = 3.2, fontface = "bold") +
  annotate("segment", x = 0.5, xend = 4.5, y = 1.40, yend = 1.40,
           color = "#C0143B", linewidth = 0.5, linetype = "dashed") +
  annotate("text", x = 2.5, y = 1.30,
           label = "(off scale, 50x higher than max sample)",
           color = "#C0143B", size = 2.6, fontface = "italic") +
  scale_fill_manual(values = c("Max pct (10 spatial)" = "#7E3F98",
                               "Final sc report"      = "#D6BDA1"),
                    name = NULL) +
  scale_y_continuous(breaks = c(0, 0.5, 1.0),
                     limits = c(0, 1.7),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "Expression pct (%)") +
  theme_paper() +
  theme(
    axis.text.x  = element_text(size = FS$axis_text, face = "bold.italic"),
    legend.position = "right"
  )

save_paper_fig(p, "fig10_silent_inflammation", canvas = "WIDE")
