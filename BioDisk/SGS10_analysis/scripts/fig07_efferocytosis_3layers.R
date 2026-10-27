#!/usr/bin/env Rscript
# ============================================================================
# fig07_efferocytosis_3layers.R  (v2: clean panels)
# X3 efferocytosis parallel 3-layer response.
#   Layer 1: Liver HSC_Activated_ECM ligand (Gas6, C3)        - brown
#   Layer 2: Liver Kupffer professional phagocyte (Cd74/Csf1r/Axl) - orange
#   Layer 3: Brain Pericytes non-prof uptake (Pecam1/Bsg/Calr/Mfge8) - purple
# ============================================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})
source("/mnt/BioDisk/wangmingyuan/SGS10_analysis/scripts/paper_theme.R")

hsc <- load_deg("liver", "tumor_vs_control", "HSC_Activated_ECM")
cat("[fig07] Layer 1 - HSC ligand:\n")
p1 <- build_lollipop(
  hsc, c("Gas6", "C3"),
  x_lab = "Liver HSC ligand (Tumor vs Control)",
  color_up = "#8B4513",
  show_legend = FALSE, show_size_legend = FALSE
)

kc <- load_deg("liver", "tumor_vs_control", "Kupffer")
cat("[fig07] Layer 2 - Kupffer phagocyte:\n")
p2 <- build_lollipop(
  kc, c("Axl", "Cd74", "Csf1r"),
  x_lab = "Liver Kupffer phagocyte (Tumor vs Control)",
  color_up = "#E67E22",
  show_legend = FALSE, show_size_legend = FALSE
)

peri <- load_deg("brain", "Tumor_vs_Control", "Pericytes", source = "deliverables_brain")
cat("[fig07] Layer 3 - Brain Pericytes uptake:\n")
p3 <- build_lollipop(
  peri, c("Pecam1", "Bsg", "Calr", "Mfge8"),
  x_lab = "Brain Pericytes (Tumor vs Control)",
  color_up = "#6A5ACD",
  show_legend = FALSE, show_size_legend = FALSE
)

p_combined <- p1 / p2 / p3 + plot_layout(heights = c(2, 3, 4))

save_paper_fig(p_combined, "fig07_efferocytosis_3layers", canvas = "LARGE")
