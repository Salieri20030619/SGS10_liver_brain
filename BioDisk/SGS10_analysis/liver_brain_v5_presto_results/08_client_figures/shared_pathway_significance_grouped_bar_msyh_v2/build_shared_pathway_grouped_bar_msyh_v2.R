options(cli.width=100L, cli.dynamic=FALSE, cli.ansi=FALSE)
Sys.setenv(R_CLI_DYNAMIC="false", NO_COLOR="1")
Sys.setenv(
  FONTCONFIG_FILE = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/liver_brain_v5_presto_results/08_client_figures/shared_pathway_significance_grouped_bar_msyh_v2/sgs10_msyh_fontconfig_v1.conf"
)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
})

project_root <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
source_dir <- file.path(
  project_root, "liver_brain_v5_presto_results", "08_client_figures",
  "shared_pathway_significance_dumbbell_v1"
)
out_dir <- file.path(
  project_root, "liver_brain_v5_presto_results", "08_client_figures",
  "shared_pathway_significance_grouped_bar_msyh_v2"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

input_tsv <- file.path(source_dir, "shared_KEGG_pathway_significance_long_v1.tsv")
plot_df <- read.delim(input_tsv, check.names = FALSE, stringsAsFactors = FALSE)
if (nrow(plot_df) != 26L) stop("Expected 26 pathway-organ rows, got ", nrow(plot_df))

wrap_label <- function(x, width = 42L) {
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
    organ = factor(organ, levels = c("Brain", "Liver"), labels = c("脑", "肝")),
    fdr_label = paste0("FDR=", formatC(fdr, format = "e", digits = 1L))
  )

organ_colors <- c("脑" = "#3B82B4", "肝" = "#D94B42")
threshold <- -log10(0.05)
x_max <- max(plot_df$neg_log10_fdr) + 1.08

p <- ggplot(plot_df, aes(x = pathway_factor, y = neg_log10_fdr, fill = organ)) +
  geom_hline(
    yintercept = threshold,
    linetype = "22", linewidth = 0.55, colour = "#737373"
  ) +
  geom_col(
    position = position_dodge(width = 0.82), width = 0.68,
    colour = "white", linewidth = 0.28
  ) +
  geom_text(
    aes(label = fdr_label, colour = organ),
    position = position_dodge(width = 0.82),
    hjust = -0.10, size = 3.55, fontface = "bold",
    family = "Microsoft YaHei", show.legend = FALSE
  ) +
  annotate(
    "text",
    x = length(path_levels) + 0.53, y = threshold + 0.06,
    label = "FDR=0.05", hjust = 0, vjust = 0.5,
    size = 3.35, family = "Microsoft YaHei", colour = "#5E5E5E"
  ) +
  coord_flip(clip = "off") +
  scale_x_discrete(labels = wrap_label(path_levels, width = 42L)) +
  scale_y_continuous(
    limits = c(0, x_max), breaks = pretty_breaks(n = 8),
    expand = expansion(mult = c(0, 0.005))
  ) +
  scale_fill_manual(values = organ_colors, breaks = c("脑", "肝")) +
  scale_colour_manual(values = organ_colors, breaks = c("脑", "肝")) +
  labs(
    title = "肝脑共同富集 KEGG 通路显著性比较",
    subtitle = "肿瘤组 vs 对照组 | 13 条共同通路 | 柱端标注为 BH-FDR",
    x = NULL,
    y = expression(-log[10]("BH-FDR")),
    fill = "器官",
    caption = paste0(
      "每个器官显示各细胞类型富集结果中的最小 BH-FDR。",
      "Tumor n = 3; Control n = 1; 结果为探索性。"
    )
  ) +
  theme_minimal(base_size = 12.5, base_family = "Microsoft YaHei") +
  theme(
    plot.title = element_text(
      face = "bold", size = 19, colour = "#20262C",
      hjust = 0, margin = margin(b = 6)
    ),
    plot.subtitle = element_text(
      size = 11.5, colour = "#59636C", hjust = 0,
      margin = margin(b = 10)
    ),
    plot.caption = element_text(
      size = 9.2, colour = "#616A72", hjust = 0,
      lineheight = 1.12, margin = margin(t = 13)
    ),
    axis.title.x = element_text(
      size = 12.5, colour = "#30373D", margin = margin(t = 9)
    ),
    axis.text.x = element_text(size = 10.5, colour = "#3E474E"),
    axis.text.y = element_text(
      size = 10.8, colour = "#252C31", lineheight = 1.04,
      margin = margin(r = 12)
    ),
    panel.grid.major.y = element_line(colour = "#E2E6E9", linewidth = 0.48),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.justification = "left",
    legend.box.just = "left",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10.8),
    legend.key.height = unit(0.52, "cm"),
    legend.key.width = unit(0.72, "cm"),
    plot.margin = margin(18, 24, 16, 20)
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

pdf_file <- file.path(out_dir, "shared_KEGG_pathway_significance_grouped_bar_msyh_v2.pdf")
png_file <- file.path(out_dir, "shared_KEGG_pathway_significance_grouped_bar_msyh_v2.png")

ggsave(
  pdf_file, p, device = cairo_pdf,
  width = 16.0, height = 10.5, units = "in"
)
ggsave(
  png_file, p,
  width = 16.0, height = 10.5, units = "in",
  dpi = 240, bg = "white"
)

validation <- c(
  "status\tPASS",
  paste0("generated_at\t", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  paste0("font_family_requested\tMicrosoft YaHei"),
  paste0("n_rows\t", nrow(plot_df)),
  paste0("n_pathways\t", n_distinct(plot_df$pathway_id)),
  paste0("pdf_exists\t", file.exists(pdf_file)),
  paste0("pdf_size_bytes\t", file.info(pdf_file)$size),
  paste0("png_exists\t", file.exists(png_file)),
  paste0("png_size_bytes\t", file.info(png_file)$size)
)
writeLines(validation, file.path(out_dir, "generation_validation_v2.txt"))
cat(paste(validation, collapse = "\n"), "\n")
