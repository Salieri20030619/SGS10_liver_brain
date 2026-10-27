#!/usr/bin/env Rscript
# ============================================================================
# fig06_thyroid_cascade.R  (v2: clean panels)
# X2 thyroid cross-organ cascade (3 layers).
#   Layer 1: Liver Hep_PV Ttr        padj 1.8e-2  (orange)
#   Layer 2: Brain CPECs Ttr          padj 9.6e-10 (key orange)
#   Layer 3: Brain Endothelial Thra   padj 2.5e-6  (brain purple)
# ============================================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})
source("/mnt/BioDisk/wangmingyuan/SGS10_analysis/scripts/paper_theme.R")

hep_pv <- load_deg("liver", "tumor_vs_control", "Hepatocyte_PV")
cat("[fig06] Layer 1 - Liver Hep_PV:\n")
p1 <- build_lollipop(
  hep_pv, c("Ttr", "Thrb", "Thra"),
  x_lab = "Liver Hep_PV (Tumor vs Control)",
  color_up = "#D2691E",
  show_legend = FALSE, show_size_legend = FALSE
)

cpec <- load_deg("brain", "Tumor_vs_Control", "CPECs", source = "deliverables_brain")
cat("[fig06] Layer 2 - Brain CPECs:\n")
p2 <- build_lollipop(
  cpec, c("Ttr", "Atp1a1"),
  x_lab = "Brain CPECs (Tumor vs Control)",
  color_up = "#E66D50",
  show_legend = FALSE, show_size_legend = FALSE
)

endo <- load_deg("brain", "Tumor_vs_Control", "Endothelial", source = "deliverables_brain")
cat("[fig06] Layer 3 - Brain Endothelial:\n")
p3 <- build_lollipop(
  endo, c("Thra", "Atp1a2", "Ctnnb1", "Cdh5"),
  x_lab = "Brain Endothelial (Tumor vs Control)",
  color_up = "#6A5ACD",
  show_legend = FALSE, show_size_legend = FALSE
)

p_combined <- p1 / p2 / p3 + plot_layout(heights = c(3, 2, 4))

save_paper_fig(p_combined, "fig06_thyroid_cascade", canvas = "LARGE")
