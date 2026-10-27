#!/usr/bin/env Rscript
# ============================================================================
# fig08_Excit_L45_TF_tvsb.R  (v2: ns shown in grey)
# Excitatory_L45 TF activity (Tumor vs Blank, complement to fig03 Tumor vs Control).
# Cross-contrast comparison: T_vs_C shows strong sig TFs, T_vs_B shows attenuated.
# Data: dorothea_results_brain_clean/per_contrast/
# ============================================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})
source("/mnt/BioDisk/wangmingyuan/SGS10_analysis/scripts/paper_theme.R")

csv_path <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis/dorothea_results_brain_clean/per_contrast/brain__tumor_vs_blank__Excitatory_L45_tf_activity.csv"
df <- read.csv(csv_path, stringsAsFactors = FALSE)
cat(sprintf("[fig08] %d TFs in table\n", nrow(df)))

focus <- c("Tfap2a", "Myc", "Pou2f2", "Tead4", "Cebpa", "Nr2c2")
sub <- df %>%
  filter(TF %in% focus) %>%
  arrange(score) %>%
  mutate(TF = factor(TF, levels = TF),
         is_sig = padj < 0.05,
         color_grp = ifelse(!is_sig, "ns",
                            ifelse(score >= 0, "up", "down")),
         color_grp = factor(color_grp, levels = c("down", "ns", "up")))

print(sub %>% select(TF, score, padj, color_grp))

hj <- ifelse(sub$score >= 0, -0.15, 1.15)

p <- ggplot(sub, aes(x = score, y = TF, color = color_grp)) +
  geom_segment(aes(xend = 0, yend = TF), linewidth = 0.7) +
  geom_point(size = 5) +
  geom_text(aes(label = sprintf("padj=%.0e", padj)),
            hjust = hj, size = 2.8, color = "grey25") +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.4) +
  scale_color_manual(
    values = c("down" = PALETTE$down, "ns" = PALETTE$ns, "up" = PALETTE$up),
    labels = c("down" = "Down (sig)", "ns" = "ns", "up" = "Up (sig)"),
    name = NULL, drop = FALSE,
    limits = c("up", "down", "ns")
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.20, 0.20))) +
  labs(x = "TF activity score (Tumor vs Blank)", y = NULL) +
  theme_paper() +
  theme(
    axis.text.y = element_text(size = FS$axis_text, face = "bold"),
    legend.position = "right"
  )

save_paper_fig(p, "fig08_Excit_L45_TF_tvsb", canvas = "WIDE")
