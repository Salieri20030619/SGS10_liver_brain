#!/usr/bin/env Rscript
# ============================================================================
# fig11_X1_TGFb_pair.R  (v2: brain side uses clean DEG)
# X1 TGFb axis: parallel autocrine activation in Liver HSC + Brain Microglia.
# Top:    Liver HSC_Activated_ECM     - Tgfb1 / Smad3 / Tgfbr2 (host autocrine)
# Bottom: Brain Microglia_Homeostatic - Smad3 / Apoe / Lrp1 / Spred2 / Tgfbr2 / Cx3cr1
#
# Data:
#   Liver: diff_stage1_results/deg_v4_refined/liver        (source="diff_stage1")
#   Brain: SGS10_deliverables_2026-04-20/.../brain_clean   (source="deliverables_brain")
# README expected values:
#   Liver HSC_ECM Tgfb1, Smad3, Tgfbr2: most NS in this contrast (host autocrine
#     hard to detect at this cell type; expected signal at host-zone spatial layer).
#   Brain MG (clean): Smad3 3.7e-13, Apoe 2.3e-18, Lrp1 7.9e-19, Tgfbr2 2.5e-5,
#                     Spred2 6.3e-10, Cx3cr1 2.4e-22.
# ============================================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})
source("/mnt/BioDisk/wangmingyuan/SGS10_analysis/scripts/paper_theme.R")

# Liver HSC (diff_stage1)
hsc <- load_deg("liver", "tumor_vs_control", "HSC_Activated_ECM")
cat("[fig11] Liver HSC panel:\n")
p1 <- build_lollipop(
  hsc, c("Tgfb1", "Smad3", "Tgfbr2"),
  x_lab = "Liver HSC_Activated_ECM (Tumor vs Control)",
  color_up = "#8B4513",
  show_legend = FALSE, show_size_legend = FALSE
)

# Brain Microglia (clean, deliverables source)
mg <- load_deg("brain", "Tumor_vs_Control", "Microglia_Homeostatic",
               source = "deliverables_brain")
cat("[fig11] Brain Microglia panel (clean):\n")
p2 <- build_lollipop(
  mg, c("Smad3", "Apoe", "Lrp1", "Spred2", "Tgfbr2", "Cx3cr1"),
  x_lab = "Brain Microglia_Homeostatic clean (Tumor vs Control)",
  color_up = "#6A5ACD",
  show_legend = FALSE, show_size_legend = FALSE
)

p_combined <- p1 / p2 + plot_layout(heights = c(3, 6))

save_paper_fig(p_combined, "fig11_X1_TGFb_pair", canvas = "LARGE")
