#!/usr/bin/env Rscript
# ============================================================================
# figS2_species_separation.R
# Supplementary Fig S2: HepG2 xenograft species separation.
# 3 tumor liver samples (gan56, Liver_3, Liver_8) x Human/Mouse spot proportion.
# Data: species_assignment/species_summary.csv
# Note: Only 2 classes (human/mouse) in the assignment table, paper draft
#       mentioned 4 (Human/Mouse/Doublet/LowUMI) — Doublet/LowUMI was a
#       pre-processing filter, removed before this table.
# ============================================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})
source("/mnt/BioDisk/wangmingyuan/SGS10_analysis/scripts/paper_theme.R")

df <- read.csv("/mnt/BioDisk/wangmingyuan/SGS10_analysis/species_assignment/species_summary.csv",
               stringsAsFactors = FALSE)
cat("[figS2] species_summary raw:\n")
print(df)

# Normalize sample names to paper convention: gan56 / Liver_3 / Liver_8
sample_rename <- c("gan56" = "gan56",
                   "liver3" = "Liver_3",
                   "liver8" = "Liver_8")
df$sample_disp <- sample_rename[df$sample]

# Long format
df_long <- df %>%
  pivot_longer(cols = c(human, mouse),
               names_to = "species", values_to = "n_spots") %>%
  mutate(
    species = factor(species, levels = c("mouse", "human"),
                     labels = c("Mouse (mm10)", "Human (GRCh38)")),
    sample_disp = factor(sample_disp,
                         levels = c("gan56", "Liver_3", "Liver_8"))
  )

# Pct labels per bar segment
df_long <- df_long %>%
  group_by(sample_disp) %>%
  mutate(total = sum(n_spots),
         pct = n_spots / total * 100) %>%
  ungroup()
cat("\n[figS2] Long form:\n")
print(df_long)

p <- ggplot(df_long, aes(x = sample_disp, y = pct, fill = species)) +
  geom_col(width = 0.65, color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", pct, n_spots)),
            position = position_stack(vjust = 0.5),
            color = "white", fontface = "bold", size = 3.4) +
  scale_fill_manual(values = c("Mouse (mm10)" = "#3B82C4",
                               "Human (GRCh38)" = "#C0143B"),
                    name = NULL) +
  scale_y_continuous(limits = c(0, 100),
                     breaks = c(0, 25, 50, 75, 100),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "Spot proportion (%)") +
  theme_paper() +
  theme(
    axis.text.x = element_text(face = "bold", size = 11),
    legend.position = "right"
  )

save_paper_fig(p, "figS2_species_separation",
               canvas = NULL, width = 7, height = 5)
