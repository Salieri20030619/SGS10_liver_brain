#!/usr/bin/env Rscript
# ============================================================================
# fig04_HSC_DEG_lollipop_pb.R
# Alternative HSC DEG lollipop using pseudobulk edgeR (tumor_core vs baseline).
# README mentioned Acta2 log2FC 8.36 / Serpine1 8.83 / Col1a1 4.18, which are
# pseudobulk-scale fold changes (Seurat avg_log2FC cannot reach this range).
# Data: D3_DEG_output/pseudobulk_edgeR/A1_tumor_core_vs_baseline_host.csv
# Cols: logFC, logCPM, F, PValue, FDR, gene
# ============================================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})
source("/mnt/BioDisk/wangmingyuan/SGS10_analysis/scripts/paper_theme.R")

# ---- Load pseudobulk DEG ----
csv_path <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis/D3_DEG_output/pseudobulk_edgeR/A1_tumor_core_vs_baseline_host.csv"
df <- read.csv(csv_path, stringsAsFactors = FALSE)

# ---- Pick HSC myofibroblast / ECM transition genes ----
top_genes <- c("Serpine1", "Acta2", "Col1a1", "F2r", "Tagln",
               "Pmepa1", "Lox", "Col1a2", "Loxl1", "Smad3")

sub <- df %>%
  filter(gene %in% top_genes) %>%
  arrange(desc(logFC)) %>%             # biggest UP at top
  mutate(
    gene = factor(gene, levels = rev(gene)),
    neglog10_padj = -log10(FDR),
    is_sig = FDR < 0.05
  )

cat("[fig04] Pseudobulk values:\n")
print(sub %>% select(gene, logFC, FDR))

missing <- setdiff(top_genes, df$gene)
if (length(missing) > 0) {
  cat(sprintf("[fig04] WARNING: missing genes: %s\n",
              paste(missing, collapse = ", ")))
}

# ---- Plot ----
p <- ggplot(sub, aes(x = logFC, y = gene)) +
  geom_segment(aes(xend = 0, yend = gene, color = is_sig),
               linewidth = 0.7) +
  geom_point(aes(size = neglog10_padj, color = is_sig)) +
  geom_text(aes(label = sprintf("padj=%.1e", FDR)),
            hjust = -0.15, size = 2.8, color = "grey25") +
  scale_color_manual(values = c("TRUE" = PALETTE$up, "FALSE" = PALETTE$ns),
                     labels = c("TRUE" = "sig (FDR<0.05)", "FALSE" = "ns"),
                     name = NULL) +
  scale_size_continuous(range = c(2.5, 7),
                        name = expression(-log[10]~"(FDR)")) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.30))) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.4) +
  labs(x = "logFC (Tumor core vs Baseline host, pseudobulk)", y = NULL) +
  theme_paper() +
  theme(
    axis.text.y = element_text(size = FS$axis_text, face = "bold"),
    legend.position = "right"
  )

save_paper_fig(p, "fig04_HSC_DEG_lollipop_pb", canvas = "WIDE")
