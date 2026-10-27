options(cli.width=100L, cli.dynamic=FALSE, cli.ansi=FALSE)
Sys.setenv(R_CLI_DYNAMIC="false", NO_COLOR="1")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
})

project_root <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
source_root <- file.path(project_root, "liver_brain_v5_presto_results")
out_dir <- file.path(source_root, "08_client_figures", "shared_pathway_significance_dumbbell_v1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_organ_kegg <- function(organ) {
  organ_lc <- tolower(organ)
  input_dir <- file.path(
    source_root,
    if (organ == "Liver") "03_liver_enrichment" else "04_brain_enrichment"
  )
  pattern <- paste0("^", organ_lc, "_tumor_vs_control_.*_KEGG\\.csv$")
  files <- sort(list.files(input_dir, pattern = pattern, full.names = TRUE))
  if (length(files) == 0L) stop("No KEGG files found for ", organ)

  rows <- lapply(files, function(fp) {
    x <- read.csv(fp, check.names = FALSE, stringsAsFactors = FALSE)
    if (nrow(x) == 0L) return(NULL)
    x$source_file <- basename(fp)
    x$source_celltype <- sub(
      paste0("^", organ_lc, "_tumor_vs_control_(.*)_KEGG\\.csv$"),
      "\\1", basename(fp)
    )
    x$organ <- organ
    x
  })
  bind_rows(rows)
}

all_rows <- bind_rows(read_organ_kegg("Brain"), read_organ_kegg("Liver")) %>%
  filter(is.finite(pvalue), is.finite(p.adjust), p.adjust > 0, p.adjust <= 1)

significant_rows <- all_rows %>% filter(p.adjust < 0.05)
shared_ids <- Reduce(intersect, split(significant_rows$ID, significant_rows$organ))
if (length(shared_ids) == 0L) stop("No shared significant pathways")

best_rows <- significant_rows %>%
  filter(ID %in% shared_ids) %>%
  group_by(organ, ID, Description) %>%
  arrange(p.adjust, pvalue, desc(Count), source_celltype, .by_group = TRUE) %>%
  mutate(n_source_hits = n_distinct(source_celltype)) %>%
  slice(1L) %>%
  ungroup() %>%
  transmute(
    pathway_id = ID,
    pathway = Description,
    organ,
    source_celltype,
    source_file,
    p_value = pvalue,
    fdr = p.adjust,
    q_value = qvalue,
    neg_log10_fdr = -log10(p.adjust),
    gene_ratio = GeneRatio,
    background_ratio = BgRatio,
    rich_factor = RichFactor,
    fold_enrichment = FoldEnrichment,
    gene_count = Count,
    n_source_hits,
    fdr_family = "BH adjustment within the source cell-type KEGG result",
    aggregation_rule = "Minimum source-level BH-FDR within organ; ties by raw P, gene count, source name",
    evidence_level = "Exploratory: Tumor n=3 vs Control n=1; no organ-level replicate inference"
  )

if (anyDuplicated(best_rows[c("pathway_id", "organ")])) {
  stop("Duplicate pathway-organ rows after selection")
}

wide <- best_rows %>%
  select(pathway_id, pathway, organ, source_celltype, p_value, fdr,
         neg_log10_fdr, gene_ratio, rich_factor, fold_enrichment, gene_count,
         n_source_hits) %>%
  pivot_wider(
    names_from = organ,
    values_from = c(source_celltype, p_value, fdr, neg_log10_fdr,
                    gene_ratio, rich_factor, fold_enrichment, gene_count,
                    n_source_hits),
    names_glue = "{organ}_{.value}"
  ) %>%
  mutate(combined_score = Brain_neg_log10_fdr + Liver_neg_log10_fdr) %>%
  arrange(desc(combined_score), pathway_id)

write.table(
  best_rows %>% arrange(pathway, organ),
  file.path(out_dir, "shared_KEGG_pathway_significance_long_v1.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE, na = "NA"
)
write.table(
  wide,
  file.path(out_dir, "shared_KEGG_pathway_significance_wide_v1.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE, na = "NA"
)

wrap_label <- function(x, width = 38L) {
  vapply(x, function(z) paste(strwrap(z, width = width), collapse = "\n"), character(1))
}

order_df <- wide %>% arrange(combined_score)
path_levels <- order_df$pathway

plot_df <- best_rows %>%
  mutate(
    pathway_factor = factor(pathway, levels = path_levels),
    y_base = as.numeric(pathway_factor),
    y_plot = y_base + if_else(organ == "Liver", 0.13, -0.13),
    fdr_label = paste0("FDR ", formatC(fdr, format = "e", digits = 1L)),
    organ = factor(organ, levels = c("Brain", "Liver"))
  )

segment_df <- wide %>%
  mutate(
    pathway_factor = factor(pathway, levels = path_levels),
    y_base = as.numeric(pathway_factor)
  )

organ_colors <- c(Brain = "#2F6B8A", Liver = "#D95F02")
x_max <- max(plot_df$neg_log10_fdr) + 1.15

p <- ggplot() +
  geom_vline(
    xintercept = -log10(0.05), linetype = "22", linewidth = 0.45,
    colour = "#8A8A8A"
  ) +
  geom_segment(
    data = segment_df,
    aes(x = Brain_neg_log10_fdr, xend = Liver_neg_log10_fdr,
        y = y_base - 0.13, yend = y_base + 0.13),
    linewidth = 0.85, colour = "#B8BEC5", lineend = "round"
  ) +
  geom_point(
    data = plot_df,
    aes(x = neg_log10_fdr, y = y_plot, size = gene_count,
        fill = organ, colour = organ),
    shape = 21, stroke = 0.45, alpha = 0.96
  ) +
  geom_text(
    data = plot_df,
    aes(x = neg_log10_fdr + 0.08, y = y_plot, label = fdr_label,
        colour = organ),
    hjust = 0, size = 3.05, fontface = "bold", show.legend = FALSE
  ) +
  annotate(
    "text", x = -log10(0.05) + 0.05, y = length(path_levels) + 0.58,
    label = "FDR = 0.05", hjust = 0, vjust = 0.5,
    size = 3.05, colour = "#666666"
  ) +
  scale_x_continuous(
    limits = c(0, x_max),
    breaks = pretty_breaks(n = 7),
    expand = expansion(mult = c(0.005, 0.005))
  ) +
  scale_y_continuous(
    breaks = seq_along(path_levels),
    labels = wrap_label(path_levels, width = 38L),
    limits = c(0.45, length(path_levels) + 0.72),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_fill_manual(values = organ_colors, breaks = c("Brain", "Liver")) +
  scale_colour_manual(values = organ_colors, breaks = c("Brain", "Liver")) +
  scale_size_continuous(
    range = c(3.1, 8.2), breaks = c(5, 15, 30, 60, 100),
    name = "Enriched genes\n(min-FDR source)"
  ) +
  guides(
    fill = guide_legend(override.aes = list(size = 4.8), order = 1),
    colour = "none",
    size = guide_legend(order = 2)
  ) +
  labs(
    title = "Shared KEGG pathway enrichment in brain and liver",
    subtitle = "Tumor vs Control | paired significance profile across 13 shared pathways",
    x = expression(-log[10]("BH-FDR")),
    y = NULL,
    fill = "Organ",
    caption = paste0(
      "Point labels show exact source-level BH-FDR; point size shows enriched gene count. ",
      "Each organ value is the minimum FDR across its cell-type enrichment results.\n",
      "Exploratory evidence: Tumor n=3 vs Control n=1; values are not organ-level replicate-based inference."
    )
  ) +
  theme_minimal(base_size = 11.5, base_family = "Noto Sans CJK SC") +
  theme(
    plot.title = element_text(face = "bold", size = 17, colour = "#1E2A33", margin = margin(b = 5)),
    plot.subtitle = element_text(size = 10.8, colour = "#58636C", margin = margin(b = 13)),
    plot.caption = element_text(size = 8.2, colour = "#5F6870", hjust = 0, lineheight = 1.15,
                                margin = margin(t = 12)),
    axis.title.x = element_text(size = 11.5, colour = "#2A3137", margin = margin(t = 8)),
    axis.text.x = element_text(size = 9.5, colour = "#40484F"),
    axis.text.y = element_text(size = 9.1, colour = "#20272D", lineheight = 0.92,
                               margin = margin(r = 10)),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour = "#E5E8EB", linewidth = 0.42),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 9.5),
    legend.text = element_text(size = 9),
    plot.margin = margin(16, 18, 14, 18)
  )

pdf_file <- file.path(out_dir, "shared_KEGG_pathway_significance_dumbbell_v1.pdf")
png_file <- file.path(out_dir, "shared_KEGG_pathway_significance_dumbbell_v1.png")

ggsave(pdf_file, p, device = cairo_pdf, width = 14.5, height = 9.2, units = "in")
ggsave(png_file, p, width = 14.5, height = 9.2, units = "in", dpi = 240, bg = "white")

validation <- c(
  paste0("status\tPASS"),
  paste0("generated_at\t", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  paste0("n_shared_pathways\t", length(shared_ids)),
  paste0("n_long_rows\t", nrow(best_rows)),
  paste0("n_wide_rows\t", nrow(wide)),
  paste0("all_fdr_below_0.05\t", all(best_rows$fdr < 0.05)),
  paste0("unique_pathway_organ\t", !anyDuplicated(best_rows[c("pathway_id", "organ")])),
  paste0("pdf_exists\t", file.exists(pdf_file)),
  paste0("pdf_size_bytes\t", file.info(pdf_file)$size),
  paste0("png_exists\t", file.exists(png_file)),
  paste0("png_size_bytes\t", file.info(png_file)$size)
)
writeLines(validation, file.path(out_dir, "validation_v1.txt"))

cat(paste(validation, collapse = "\n"), "\n")
cat("OUTPUT_DIR\t", out_dir, "\n", sep = "")
