#!/usr/bin/env Rscript
# ============================================================================
# fig09_TF_dual_mode.R  (v2: ggrepel labels, top-only labeling)
# Brain dual-mode paradigm: non-neuron TF-decoupled vs Excit_L45 TF-driven.
# Each panel shows TF activity score vs -log10(TF padj); red dashed = 0.05 cutoff.
# Data: dorothea_results_brain_clean/per_contrast/*tumor_vs_control*_tf_activity.csv
# ============================================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
})
source("/mnt/BioDisk/wangmingyuan/SGS10_analysis/scripts/paper_theme.R")

TF_DIR <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis/dorothea_results_brain_clean/per_contrast"

celltypes <- list(
  "Microglia_Homeostatic" = "Microglia (non-neuron)",
  "Astrocyte_Gfap"        = "Astrocyte (non-neuron)",
  "Endothelial"           = "Endothelial (non-neuron)",
  "Pericytes"             = "Pericytes (non-neuron)",
  "Excitatory_L45"        = "Excitatory_L45 (NEURON)"
)

dfs <- list()
for (ct in names(celltypes)) {
  fn <- sprintf("brain__tumor_vs_control__%s_tf_activity.csv", ct)
  full <- file.path(TF_DIR, fn)
  if (!file.exists(full)) {
    cat(sprintf("[fig09] MISSING %s\n", full)); next
  }
  d <- read.csv(full, stringsAsFactors = FALSE)
  d$cell_type <- celltypes[[ct]]
  d$is_neuron <- grepl("NEURON", celltypes[[ct]])
  dfs[[ct]] <- d
  cat(sprintf("[fig09] %s: %d TFs, %d strict sig\n",
              ct, nrow(d), sum(d$padj < 0.05, na.rm = TRUE)))
}
all_df <- bind_rows(dfs)
all_df$neglog10_padj <- -log10(pmax(all_df$padj, 1e-300))
all_df$cell_type <- factor(all_df$cell_type, levels = unname(unlist(celltypes)))

# Only label sig TFs (padj < 0.05), top by padj per panel
sig_label <- all_df %>%
  filter(padj < 0.05) %>%
  group_by(cell_type) %>%
  slice_min(padj, n = 5) %>%
  ungroup()
cat(sprintf("[fig09] %d labels total\n", nrow(sig_label)))

# Compute uniform y limit
y_max <- max(all_df$neglog10_padj) * 1.10

p <- ggplot(all_df, aes(x = score, y = neglog10_padj, color = is_neuron)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = "red", linewidth = 0.4) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.4) +
  geom_point(size = 1.4, alpha = 0.55, shape = 16) +
  geom_text_repel(
    data = sig_label,
    aes(label = TF),
    size = 2.8, fontface = "bold",
    max.overlaps = Inf,
    box.padding = 0.5, point.padding = 0.3,
    force = 4, force_pull = 0.3,
    min.segment.length = 0,
    segment.size = 0.25, segment.color = "grey60", segment.alpha = 0.7,
    seed = 42, show.legend = FALSE
  ) +
  scale_color_manual(values = c("FALSE" = PALETTE$baseline,
                                "TRUE"  = PALETTE$tumor),
                     labels = c("FALSE" = "non-neuron (TF-decoupled)",
                                "TRUE"  = "neuron (TF-driven)"),
                     name = NULL) +
  scale_y_continuous(limits = c(0, y_max), expand = expansion(mult = c(0, 0))) +
  facet_wrap(~ cell_type, ncol = 5, scales = "free_x") +
  labs(x = "TF activity score (Tumor vs Control)",
       y = expression(-log[10]~"(TF padj)")) +
  theme_paper() +
  theme(
    strip.text = element_text(face = "bold", size = 9),
    legend.position = "bottom",
    panel.spacing = unit(0.6, "lines")
  )

save_paper_fig(p, "fig09_TF_dual_mode", canvas = "LARGE",
               width = 11, height = 4.5)
