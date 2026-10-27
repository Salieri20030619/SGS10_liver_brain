#!/usr/bin/env Rscript
# ============================================================================
# fig12_liver_microenv_heatmap.R  (v3: no row_split, continuous matrix)
# 7 liver microenv x top 10 markers. Rows ordered by microenv, but no gaps.
# ============================================================================
suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
})
source("/mnt/BioDisk/wangmingyuan/SGS10_analysis/scripts/paper_theme.R")

csv_path <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis/D1_v2_RCTD_composition/microenv_markers_RCTD100/liver_microenv_markers_top20.csv"
df <- read.csv(csv_path, stringsAsFactors = FALSE)
cat(sprintf("[fig12] %d rows loaded\n", nrow(df)))

# Strip BOTH species prefixes
df$gene_clean <- gsub("^mm10---|^GRCh38-", "", df$gene)

top10 <- df %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 10) %>%
  ungroup()
cat(sprintf("[fig12] %d top markers across %d microenvs\n",
            nrow(top10), length(unique(top10$cluster))))

genes_ordered <- top10 %>%
  arrange(cluster, desc(avg_log2FC)) %>%
  pull(gene_clean) %>% unique()
microenvs <- sort(unique(top10$cluster))

mat <- matrix(NA, nrow = length(genes_ordered), ncol = length(microenvs),
              dimnames = list(genes_ordered, microenvs))
for (i in seq_len(nrow(top10))) {
  g <- top10$gene_clean[i]; c <- top10$cluster[i]
  mat[g, c] <- top10$avg_log2FC[i]
}

paper_key_genes <- c("Cd74", "H2-Aa", "Psap", "Inpp5d", "Ctsb", "Zeb2",
                     "Runx1", "Laptm5")
gene_face <- ifelse(genes_ordered %in% paper_key_genes, "bold.italic", "italic")
gene_color <- ifelse(genes_ordered %in% paper_key_genes, "#C0143B", "black")

col_fun <- colorRamp2(c(0, 1.5, 3, 5, 7),
                      c("#FFFFFF", "#FED8A6", "#F08C2A",
                        "#D44217", "#7E1A0F"))

col_label_colors <- ifelse(microenvs == "M05", "#C0143B", "black")

ht <- Heatmap(
  mat,
  name = "log2FC",
  col = col_fun,
  na_col = "#FAFAFA",
  cluster_rows = FALSE, cluster_columns = FALSE,
  show_row_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 8, fontface = gene_face, col = gene_color),
  show_column_names = TRUE,
  column_names_side = "top",
  column_names_rot = 0,
  column_names_centered = TRUE,
  column_names_gp = gpar(fontsize = 11, fontface = "bold",
                        col = col_label_colors),
  border = TRUE,
  rect_gp = gpar(col = "white", lwd = 0.3),
  heatmap_legend_param = list(
    title = "avg log2FC",
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 9),
    legend_height = unit(28, "mm")
  ),
  width = unit(55, "mm"),
  height = unit(150, "mm")
)

out_dir <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis/figures/paper_v2"
basename <- "fig12_liver_microenv_heatmap"
w <- 4.7; h <- 8

pdf(file.path(out_dir, "PDF", paste0(basename, ".pdf")), width = w, height = h)
draw(ht, padding = unit(c(4, 4, 4, 4), "mm"))
dev.off()

tiff(file.path(out_dir, "TIFF", paste0(basename, ".tif")),
     width = w, height = h, units = "in", res = 300, compression = "lzw")
draw(ht, padding = unit(c(4, 4, 4, 4), "mm"))
dev.off()

png(file.path(out_dir, "PNG", paste0(basename, ".png")),
    width = w, height = h, units = "in", res = 200, bg = "white")
draw(ht, padding = unit(c(4, 4, 4, 4), "mm"))
dev.off()

cat(sprintf("[OK] %s  %.1f x %.1f in  ->  PDF + TIFF + PNG\n",
            basename, w, h))

cat("\n[fig12] Microenv label key (for figure caption):\n")
key <- top10 %>% distinct(cluster, microenv_label) %>% arrange(cluster)
for (i in seq_len(nrow(key))) {
  cat(sprintf("  %s = %s\n", key$cluster[i], key$microenv_label[i]))
}
