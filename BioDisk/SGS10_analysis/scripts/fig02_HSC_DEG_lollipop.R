#!/usr/bin/env Rscript
# ============================================================================
# fig02_HSC_DEG_lollipop.R
# HSC_Activated_ECM Top 10 DEG (myofibroblast transition genes).
# Dot size proportional to -log10(padj); red = sig (padj < 0.05), grey = ns.
# Same CSV as fig01.
# ============================================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})
source("/mnt/BioDisk/wangmingyuan/SGS10_analysis/scripts/paper_theme.R")

# ---- Load ----
csv_path <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis/diff_stage1_results/deg_v4_refined/liver/liver_tumor_vs_control_HSC_Activated_ECM.csv"
df <- read.csv(csv_path, stringsAsFactors = FALSE)

# ---- Pick the genes shown in original Fig 1A (myofibroblast / ECM panel) ----
top_genes <- c("Serpine1", "Acta2", "Col1a1", "F2r", "Tagln",
               "Pmepa1", "Lox", "Col1a2", "Loxl1", "Smad3")

sub <- df %>%
  filter(gene %in% top_genes) %>%
  mutate(
    gene = factor(gene, levels = rev(top_genes)),   # top of plot first
    neglog10_padj = -log10(p_val_adj),
    is_sig = p_val_adj < 0.05
  )

# Safety: report any gene missing in CSV
missing <- setdiff(top_genes, df$gene)
if (length(missing) > 0) {
  cat(sprintf("[fig02] WARNING: missing genes: %s\n",
              paste(missing, collapse = ", ")))
}

# ---- Plot ----
p <- ggplot(sub, aes(x = avg_log2FC, y = gene)) +
  geom_segment(aes(xend = 0, yend = gene, color = is_sig),
               linewidth = 0.7) +
  geom_point(aes(size = neglog10_padj, color = is_sig)) +
  geom_text(aes(label = sprintf("padj=%.1e", p_val_adj)),
            hjust = -0.15, size = 2.8, color = "grey25") +
  scale_color_manual(values = c("TRUE" = PALETTE$up, "FALSE" = PALETTE$ns),
                     labels = c("TRUE" = "sig (padj<0.05)", "FALSE" = "ns"),
                     name = NULL) +
  scale_size_continuous(range = c(2.5, 7),
                        name = expression(-log[10]~"(padj)")) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.30))) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.4) +
  labs(x = "avg_log2FC (Tumor vs Control)", y = NULL) +
  theme_paper() +
  theme(
    axis.text.y = element_text(size = FS$axis_text, face = "bold"),
    legend.position = "right"
  )

save_paper_fig(p, "fig02_HSC_DEG_lollipop", canvas = "WIDE",
               out_dir = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/figures/restyle_v1")  # [restyle] 新图目录
