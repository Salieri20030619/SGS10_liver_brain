#!/usr/bin/env Rscript
# ============================================================================
# fig06_BBB_receptor_lollipop.R
# Brain BBB endothelial thyroid receptor response (X2 cascade terminal).
# Genes: Thra, Thrb (receptor), Slco1c1, Slc16a2 (T4 transporters)
# Data: C4_brain_distant_vs_brain_baseline.csv (brain endothelial pseudobulk)
# Cols: logFC, logCPM, F, PValue, FDR, gene
# ============================================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})
source("/mnt/BioDisk/wangmingyuan/SGS10_analysis/scripts/paper_theme.R")

csv_path <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis/D3_DEG_output/pseudobulk_edgeR/C4_brain_distant_vs_brain_baseline.csv"
df <- read.csv(csv_path, stringsAsFactors = FALSE)

# ---- Genes of interest: thyroid hormone receptor + transporters ----
genes <- c("Thra", "Thrb", "Slco1c1", "Slc16a2", "Slc16a10",
           "Crym", "Dio2", "Hr")

sub <- df %>%
  filter(gene %in% genes) %>%
  arrange(desc(logFC)) %>%
  mutate(
    gene = factor(gene, levels = rev(gene)),
    is_sig = FDR < 0.05,
    neglog10_padj = -log10(FDR)
  )

missing <- setdiff(genes, df$gene)
if (length(missing) > 0) {
  cat(sprintf("[fig06] WARN missing genes: %s\n",
              paste(missing, collapse=", ")))
}
cat("--- fig06 values ---\n")
print(sub %>% select(gene, logFC, FDR))

p <- ggplot(sub, aes(x = logFC, y = gene)) +
  geom_segment(aes(xend = 0, yend = gene, color = is_sig),
               linewidth = 0.7) +
  geom_point(aes(size = neglog10_padj, color = is_sig)) +
  geom_text(aes(label = sprintf("padj=%.1e", FDR)),
            hjust = ifelse(sub$logFC >= 0, -0.15, 1.15),
            size = 2.8, color = "grey25") +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.4) +
  scale_color_manual(values = c("TRUE" = "#2A9D8F", "FALSE" = PALETTE$ns),
                     labels = c("TRUE" = "sig (FDR<0.05)", "FALSE" = "ns"),
                     name = NULL,
                     limits = c("TRUE", "FALSE")) +
  scale_size_continuous(range = c(2.5, 7),
                        name = expression(-log[10]~"(FDR)")) +
  scale_x_continuous(expand = expansion(mult = c(0.30, 0.30))) +
  labs(x = "logFC (Brain distant vs baseline, pseudobulk)", y = NULL) +
  theme_paper() +
  theme(
    axis.text.y = element_text(size = FS$axis_text, face = "bold"),
    legend.position = "right"
  )

save_paper_fig(p, "fig06_BBB_receptor_lollipop", canvas = "WIDE")
