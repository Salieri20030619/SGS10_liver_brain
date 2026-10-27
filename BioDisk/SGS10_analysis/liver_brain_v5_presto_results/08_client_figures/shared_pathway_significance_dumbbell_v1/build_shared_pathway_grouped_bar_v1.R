options(cli.width=100L, cli.dynamic=FALSE, cli.ansi=FALSE)
Sys.setenv(R_CLI_DYNAMIC="false", NO_COLOR="1")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
})

project_root <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
out_dir <- file.path(
  project_root, "liver_brain_v5_presto_results", "08_client_figures",
  "shared_pathway_significance_dumbbell_v1"
)
input_tsv <- file.path(out_dir, "shared_KEGG_pathway_significance_long_v1.tsv")
if (!file.exists(input_tsv)) stop("Missing frozen significance table: ", input_tsv)

plot_df <- read.delim(input_tsv, check.names = FALSE, stringsAsFactors = FALSE)
if (nrow(plot_df) != 26L) stop("Expected 26 pathway-organ rows, got ", nrow(plot_df))

wrap_label <- function(x, width = 38L) {
  vapply(x, function(z) paste(strwrap(z, width = width), collapse = "\n"), character(1))
}

order_df <- plot_df %>%
  group_by(pathway_id, pathway) %>%
  summarise(combined_score = sum(neg_log10_fdr), .groups = "drop") %>%
  arrange(combined_score)
path_levels <- order_df$pathway

plot_df <- plot_df %>%
  mutate(
    pathway_factor = factor(pathway, levels = path_levels),
    organ = factor(organ, levels = c("Brain", "Liver")),
    fdr_label = paste0("FDR ", formatC(fdr, format = "e", digits = 1L))
  )

organ_colors <- c(Brain = "#2F6B8A", Liver = "#D95F02")
x_max <- max(plot_df$neg_log10_fdr) + 1.15

p <- ggplot(plot_df, aes(x = pathway_factor, y = neg_log10_fdr, fill = organ)) +
  geom_hline(
    yintercept = -log10(0.05), linetype = "22", linewidth = 0.45,
    colour = "#8A8A8A"
  ) +
  geom_col(
    position = position_dodge(width = 0.78), width = 0.66,
    colour = "white", linewidth = 0.22
  ) +
  geom_text(
    aes(label = fdr_label, colour = organ),
    position = position_dodge(width = 0.78), hjust = -0.08,
    size = 3.0, fontface = "bold", show.legend = FALSE
  ) +
  annotate(
    "text", x = length(path_levels) + 0.50, y = -log10(0.05) + 0.06,
    label = "FDR = 0.05", hjust = 0, vjust = 0.5,
    size = 3.0, colour = "#666666"
  ) +
  coord_flip(clip = "off") +
  scale_x_discrete(labels = wrap_label(path_levels, width = 38L)) +
  scale_y_continuous(
    limits = c(0, x_max), breaks = pretty_breaks(n = 7),
    expand = expansion(mult = c(0, 0.005))
  ) +
  scale_fill_manual(values = organ_colors, breaks = c("Brain", "Liver")) +
  scale_colour_manual(values = organ_colors, breaks = c("Brain", "Liver")) +
  labs(
    title = "Shared KEGG pathway enrichment in brain and liver",
    subtitle = "Tumor vs Control | classic grouped-bar presentation",
    x = NULL,
    y = expression(-log[10]("BH-FDR")),
    fill = "Organ",
    caption = paste0(
      "Bar-end labels show exact source-level BH-FDR. Each organ value is the minimum FDR ",
      "across its cell-type enrichment results.\n",
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
    panel.grid.major.y = element_line(colour = "#E5E8EB", linewidth = 0.42),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 9.5),
    legend.text = element_text(size = 9),
    plot.margin = margin(16, 18, 14, 18)
  )

pdf_file <- file.path(out_dir, "shared_KEGG_pathway_significance_grouped_bar_v1.pdf")
png_file <- file.path(out_dir, "shared_KEGG_pathway_significance_grouped_bar_v1.png")

ggsave(pdf_file, p, device = cairo_pdf, width = 14.5, height = 9.2, units = "in")
ggsave(png_file, p, width = 14.5, height = 9.2, units = "in", dpi = 240, bg = "white")

validation <- c(
  "status\tPASS",
  paste0("generated_at\t", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  paste0("n_rows\t", nrow(plot_df)),
  paste0("n_pathways\t", n_distinct(plot_df$pathway_id)),
  paste0("pdf_exists\t", file.exists(pdf_file)),
  paste0("pdf_size_bytes\t", file.info(pdf_file)$size),
  paste0("png_exists\t", file.exists(png_file)),
  paste0("png_size_bytes\t", file.info(png_file)$size)
)
writeLines(validation, file.path(out_dir, "grouped_bar_validation_v1.txt"))
cat(paste(validation, collapse = "\n"), "\n")
