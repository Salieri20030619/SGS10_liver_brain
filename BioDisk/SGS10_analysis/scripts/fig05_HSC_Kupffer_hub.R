#!/usr/bin/env Rscript
# ============================================================================
# fig05_HSC_Kupffer_hub.R  (v2: clean panels, no redundant legends)
# HSC -> Kupffer signaling hub (two-panel lollipop, Tumor vs Control).
# Top:    HSC_Activated_ECM ligand source (Gas6, C3)
# Bottom: Kupffer receptor / marker (Cd74, Csf1r, Axl)
# ============================================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})
source("/mnt/BioDisk/wangmingyuan/SGS10_analysis/scripts/paper_theme.R")

hsc <- load_deg("liver", "tumor_vs_control", "HSC_Activated_ECM")
cat("[fig05] HSC panel:\n")
p_hsc <- build_lollipop(
  hsc, c("Gas6", "C3"),
  x_lab = "avg_log2FC (HSC_Activated_ECM, Tumor vs Control)",
  color_up = "#7E3F98",
  show_legend = FALSE, show_size_legend = FALSE
)

kc <- load_deg("liver", "tumor_vs_control", "Kupffer")
cat("[fig05] Kupffer panel:\n")
p_kc <- build_lollipop(
  kc, c("Cd74", "Csf1r", "Axl"),
  x_lab = "avg_log2FC (Kupffer, Tumor vs Control)",
  color_up = "#3B82C4",
  show_legend = FALSE, show_size_legend = FALSE
)

p_combined <- p_hsc / p_kc + plot_layout(heights = c(2, 3))

save_paper_fig(p_combined, "fig05_HSC_Kupffer_hub", canvas = "LARGE")
