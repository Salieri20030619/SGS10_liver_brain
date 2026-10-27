#!/usr/bin/env Rscript
# ============================================================================
# figS1_qc_violin.R
# Supplementary Fig S1: SeekSpace QC indicators per sample.
# 2 organ panels x 3 QC metrics violin (nFeature_RNA, nCount_RNA, percent.mt).
# Data: umap_paper_grade/{liver,brain}_integrated.rds metadata.
# ============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})
source("/mnt/BioDisk/wangmingyuan/SGS10_analysis/scripts/paper_theme.R")

RDS_DIR <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis/umap_paper_grade"

extract_qc <- function(rds_path, organ) {
  cat(sprintf("[%s] loading %s...\n", organ, basename(rds_path)))
  obj <- readRDS(rds_path)
  md <- obj@meta.data[, c("orig.ident", "nCount_RNA",
                          "nFeature_RNA", "percent.mt")]
  md$organ <- organ
  cat(sprintf("  %d cells, %d samples\n", nrow(md),
              length(unique(md$orig.ident))))
  rm(obj); gc(verbose = FALSE)
  md
}

liver_md <- extract_qc(file.path(RDS_DIR, "liver_integrated.rds"), "Liver")
brain_md <- extract_qc(file.path(RDS_DIR, "brain_integrated.rds"), "Brain")
all_md <- rbind(liver_md, brain_md)

cat("\n[figS1] Sample distribution:\n")
print(table(all_md$orig.ident, all_md$organ))

# Long form for facet
all_long <- all_md %>%
  pivot_longer(cols = c(nCount_RNA, nFeature_RNA, percent.mt),
               names_to = "metric", values_to = "value")
all_long$metric <- factor(all_long$metric,
                          levels = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
                          labels = c("nCount_RNA", "nFeature_RNA", "percent.mt (%)"))

# Order samples within each organ by name
sample_order <- c(
  # liver
  sort(unique(liver_md$orig.ident)),
  # brain
  sort(unique(brain_md$orig.ident))
)
all_long$orig.ident <- factor(all_long$orig.ident, levels = sample_order)

# Color: liver = warm, brain = cool
sample_colors <- c(
  setNames(colorRampPalette(c("#FFD8A6","#E66D50","#7E3F98"))
           (length(unique(liver_md$orig.ident))),
           sort(unique(liver_md$orig.ident))),
  setNames(colorRampPalette(c("#A8D8E5","#3B82C4","#6A5ACD"))
           (length(unique(brain_md$orig.ident))),
           sort(unique(brain_md$orig.ident)))
)

build_panel <- function(organ_name, organ_colors) {
  d <- all_long[all_long$organ == organ_name, , drop = FALSE]
  ggplot(d, aes(x = orig.ident, y = value, fill = orig.ident)) +
    geom_violin(scale = "width", trim = TRUE, color = NA, alpha = 0.85) +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white",
                 color = "black", linewidth = 0.3) +
    scale_fill_manual(values = organ_colors, guide = "none") +
    scale_y_continuous(trans = "log10",
                       labels = function(x) format(x, scientific = FALSE, big.mark = ",")) +
    facet_wrap(~ metric, nrow = 1, scales = "free_y") +
    labs(x = NULL, y = NULL) +
    theme_paper() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9,
                                 face = "bold"),
      strip.text = element_text(face = "bold", size = 10),
      panel.spacing.x = unit(0.8, "lines")
    )
}

# percent.mt should not be log-scaled (linear small values)
build_panel_safe <- function(organ_name, organ_colors) {
  d <- all_long[all_long$organ == organ_name, , drop = FALSE]
  ggplot(d, aes(x = orig.ident, y = value, fill = orig.ident)) +
    geom_violin(scale = "width", trim = TRUE, color = NA, alpha = 0.85) +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white",
                 color = "black", linewidth = 0.3) +
    scale_fill_manual(values = organ_colors, guide = "none") +
    facet_wrap(~ metric, nrow = 1, scales = "free_y") +
    labs(x = NULL, y = NULL) +
    theme_paper() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9,
                                 face = "bold"),
      strip.text = element_text(face = "bold", size = 10),
      panel.spacing.x = unit(0.8, "lines")
    )
}

p_liver <- build_panel_safe("Liver", sample_colors) +
  labs(title = NULL) +
  theme(plot.title = element_blank())

p_brain <- build_panel_safe("Brain", sample_colors) +
  theme(plot.title = element_blank())

# Add organ label on the left
p_combined <- (p_liver / p_brain) +
  plot_annotation(theme = theme(plot.margin = margin(5,5,5,5)))

# Save
save_paper_fig(p_combined, "figS1_qc_violin",
               canvas = NULL, width = 11, height = 7)
