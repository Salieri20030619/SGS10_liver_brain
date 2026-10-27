#!/usr/bin/env Rscript
# ============================================================================
# fig05_HSC_myeloid_hub.R
# HSC-myeloid signaling hub: ligand source -> receptor target (2-layer panel).
# Top: HSC_Activated_ECM ligand source (Gas6, C3, Csf1, Cxcl12)
# Bot: Kupffer receptor/marker target (Cd74, Csf1r, Axl, Mrc1, Mertk)
# Data: D3_DEG_output/pseudobulk_edgeR/A1_tumor_core_vs_baseline_host.csv
# Cols: logFC, logCPM, F, PValue, FDR, gene
# ============================================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})
source("/mnt/BioDisk/wangmingyuan/SGS10_analysis/scripts/paper_theme.R")

csv_path <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis/D3_DEG_output/pseudobulk_edgeR/A1_tumor_core_vs_baseline_host.csv"
df <- read.csv(csv_path, stringsAsFactors = FALSE)

# ---- Ligand panel (HSC -> outgoing signals) ----
ligand_genes <- c("Gas6", "C3", "Csf1", "Cxcl12")
ligand_color <- "#7E3F98"   # purple

# ---- Receptor / marker panel (Kupffer -> incoming) ----
receptor_genes <- c("Cd74", "Csf1r", "Axl", "Mrc1", "Mertk")
receptor_color <- "#3B82C4" # blue

build_panel <- function(genes, label_color) {
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
    cat(sprintf("[fig05] WARN missing: %s\n", paste(missing, collapse=", ")))
  }
  cat("--- panel values ---\n"); print(sub %>% select(gene, logFC, FDR))
  
  ggplot(sub, aes(x = logFC, y = gene)) +
    geom_segment(aes(xend = 0, yend = gene, color = is_sig),
                 linewidth = 0.7) +
    geom_point(aes(size = neglog10_padj, color = is_sig)) +
    geom_text(aes(label = sprintf("padj=%.1e", FDR)),
              hjust = ifelse(sub$logFC >= 0, -0.15, 1.15),
              size = 2.8, color = "grey25") +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.4) +
    scale_color_manual(values = c("TRUE" = label_color, "FALSE" = PALETTE$ns),
                       labels = c("TRUE" = "sig (FDR<0.05)", "FALSE" = "ns"),
                       name = NULL,
                       drop = FALSE,
                       limits = c("TRUE", "FALSE")) +
    scale_size_continuous(range = c(2.5, 7),
                          name = expression(-log[10]~"(FDR)")) +
    scale_x_continuous(expand = expansion(mult = c(0.20, 0.30))) +
    labs(x = "logFC (Tumor core vs Baseline host, pseudobulk)", y = NULL) +
    theme_paper() +
    theme(
      axis.text.y = element_text(size = FS$axis_text, face = "bold"),
      legend.position = "right"
    )
}

p_lig <- build_panel(ligand_genes,   ligand_color)
p_rec <- build_panel(receptor_genes, receptor_color)

# Stack two panels; collect legends
p_combined <- (p_lig / p_rec) +
  plot_layout(heights = c(4, 5), guides = "collect") &
  theme(legend.position = "right")

save_paper_fig(p_combined, "fig05_HSC_myeloid_hub", canvas = "LARGE")
