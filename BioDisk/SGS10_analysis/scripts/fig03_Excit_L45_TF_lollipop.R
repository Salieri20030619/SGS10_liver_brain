#!/usr/bin/env Rscript
# ============================================================================
# fig03_Excit_L45_TF_lollipop.R
# Excitatory_L45 transcription factor activity (tumor vs control).
# Shows Tead4/Pou2f2/Myc/Cebpa/Nr2c2 systematic DOWN; Tfap2a UP.
# Data: dorothea_results_brain_clean/per_contrast/
#       brain__tumor_vs_control__Excitatory_L45_tf_activity.csv
# Cols: TF, score, p_value, padj, sig_strict, sig_relaxed, sig
# ============================================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})
source("/mnt/BioDisk/wangmingyuan/SGS10_analysis/scripts/paper_theme.R")

# ---- Load ----
csv_path <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis/dorothea_results_brain_clean/per_contrast/brain__tumor_vs_control__Excitatory_L45_tf_activity.csv"
df <- read.csv(csv_path, stringsAsFactors = FALSE)

# ---- Pick TFs of interest (from previous Fig 3B) ----
focus <- c("Tfap2a", "Myc", "Pou2f2", "Tead4", "Cebpa", "Nr2c2")

sub <- df %>%
  filter(TF %in% focus) %>%
  mutate(
    direction = ifelse(score > 0, "Up", "Down"),
    direction = factor(direction, levels = c("Down", "Up"))
  ) %>%
  arrange(score) %>%                          # bottom = most negative
  mutate(TF = factor(TF, levels = TF))

missing <- setdiff(focus, df$TF)
if (length(missing) > 0) {
  cat(sprintf("[fig03] WARNING: missing TFs: %s\n",
              paste(missing, collapse = ", ")))
}

# ---- Plot ----
p <- ggplot(sub, aes(x = score, y = TF, color = direction)) +
  geom_segment(aes(xend = 0, yend = TF), linewidth = 0.7) +
  geom_point(size = 5) +
  geom_text(aes(label = sprintf("padj=%.0e", padj)),
            hjust = ifelse(sub$score > 0, -0.15, 1.15),
            size = 2.8, color = "grey25") +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.4) +
  scale_color_manual(values = c("Down" = PALETTE$down, "Up" = PALETTE$up),
                     labels = c("Down (padj<0.05)", "Up (padj<0.05)"),
                     name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.20, 0.20))) +
  labs(x = "TF activity score", y = NULL) +
  theme_paper() +
  theme(
    axis.text.y = element_text(size = FS$axis_text, face = "bold"),
    legend.position = "right"
  )

save_paper_fig(p, "fig03_Excit_L45_TF_lollipop", canvas = "WIDE")
