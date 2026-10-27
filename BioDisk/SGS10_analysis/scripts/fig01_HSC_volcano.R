#!/usr/bin/env Rscript
# ============================================================================
# fig01_HSC_volcano.R  (v2: label repel hardened)
# HSC_Activated_ECM Tumor vs Control volcano (no IL-6 emphasis).
# Data: diff_stage1_results/deg_v4_refined/liver/
#       liver_tumor_vs_control_HSC_Activated_ECM.csv
# ============================================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
})
source("/mnt/BioDisk/wangmingyuan/SGS10_analysis/scripts/paper_theme.R")

csv_path <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis/diff_stage1_results/deg_v4_refined/liver/liver_tumor_vs_control_HSC_Activated_ECM.csv"
df <- read.csv(csv_path, stringsAsFactors = FALSE)

# ---- Thresholds (paper consistent) ----
fc_thr   <- 0.20
padj_thr <- 0.10

df$neglog10_padj <- -log10(df$p_val_adj)
df$sig <- "NS"
df$sig[df$avg_log2FC >  fc_thr & df$p_val_adj < padj_thr] <- "Up"
df$sig[df$avg_log2FC < -fc_thr & df$p_val_adj < padj_thr] <- "Down"
df$sig <- factor(df$sig, levels = c("Down", "NS", "Up"))

n_up   <- sum(df$sig == "Up")
n_down <- sum(df$sig == "Down")
cat(sprintf("[fig01] Up=%d  Down=%d  (expect 173/38)\n", n_up, n_down))

# ---- Labeled genes (IL-6 pathway already removed) ----
label_down <- c("Gm19951", "AY036118", "Gm48099",
                "Camk1d", "Il31ra", "Ntm", "Nrxn1", "Ank3",
                "Cmss1", "Egfem1")
label_up   <- c("Col1a1", "Col3a1", "Col4a1", "Sparc", "Slit3",
                "Ahnak", "Gas6", "C3", "H2-D1", "Hspg2")
df$label <- ifelse(df$gene %in% c(label_down, label_up), df$gene, NA)

# ---- Cap extreme padj ----
y_cap <- 95
df$neglog10_padj_capped <- pmin(df$neglog10_padj, y_cap)

cols <- c("Down" = PALETTE$down, "NS" = PALETTE$ns, "Up" = PALETTE$up)

# ---- Plot.  Repel hardened: high force, large box.padding, segments always on.
p <- ggplot(df, aes(x = avg_log2FC, y = neglog10_padj_capped, color = sig)) +
  geom_hline(yintercept = -log10(padj_thr), linetype = "dashed",
             color = "black", linewidth = 0.4) +
  geom_vline(xintercept = c(-fc_thr, fc_thr), linetype = "dashed",
             color = "black", linewidth = 0.4) +
  # NS points first
  geom_point(data = subset(df, sig == "NS"),
             size = 1.0, alpha = 0.30, shape = 16) +
  # Sig points on top
  geom_point(data = subset(df, sig != "NS"),
             size = 1.6, alpha = 0.80, shape = 16) +
  # Labels: hardened to prevent overlap with points
  geom_text_repel(
    data            = subset(df, !is.na(label)),
    aes(label = label),
    size            = 3.2,
    fontface        = "bold",
    max.overlaps    = Inf,
    box.padding     = 0.8,             # much larger pad
    point.padding   = 0.5,             # push away from point
    force           = 6,               # stronger repel force
    force_pull      = 0.2,
    min.segment.length = 0,            # always draw connector segment
    segment.size    = 0.30,
    segment.color   = "grey50",
    segment.alpha   = 0.7,
    seed            = 42,
    show.legend     = FALSE
  ) +
  scale_color_manual(values = cols, name = NULL) +
  labs(x = "avg_log2FC (Tumor vs Control)",
       y = expression(-log[10]~"(p_val_adj)")) +
  theme_paper() +
  theme(legend.position = "right") +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

save_paper_fig(p, "fig01_HSC_volcano", canvas = "WIDE")
