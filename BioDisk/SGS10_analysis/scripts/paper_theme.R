# ============================================================================
# paper_theme.R - SGS10 Paper Figure Unified Theme
# ============================================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

CANVAS <- list(
  SMALL  = c(w = 3.54, h = 2.76),
  WIDE   = c(w = 7.09, h = 3.54),
  SQUARE = c(w = 3.54, h = 3.54),
  TALL   = c(w = 3.54, h = 4.50),
  LARGE  = c(w = 7.09, h = 5.50),
  HUGE   = c(w = 10.0, h = 11.0)
)

FS <- list(
  axis_text    = 13,   # [restyle] 10->13  (客户要求放大, 拼图后清晰)
  axis_title   = 14,   # [restyle] 11->14
  legend_text  = 12,   # [restyle]  9->12
  legend_title = 13,   # [restyle] 10->13
  annotate     = 11,   # [restyle]  9->11
  facet        = 13    # [restyle] 11->13
)

theme_paper <- function(base_size = FS$axis_text) {
  theme_bw(base_size = base_size, base_family = "Helvetica") +
    theme(
      text             = element_text(color = "black"),
      axis.text        = element_text(color = "black", size = FS$axis_text,
                                      face = "bold"),
      axis.title       = element_text(color = "black", size = FS$axis_title,
                                      face = "bold"),
      axis.title.x     = element_text(margin = margin(t = 5)),
      axis.title.y     = element_text(margin = margin(r = 5)),
      axis.line        = element_blank(),
      axis.ticks       = element_line(color = "black", linewidth = 0.4),
      axis.ticks.length = unit(2.5, "pt"),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.6),
      panel.grid       = element_blank(),
      plot.title       = element_blank(),
      plot.subtitle    = element_blank(),
      plot.caption     = element_blank(),
      legend.title     = element_text(size = FS$legend_title, face = "bold",
                                      color = "black"),
      legend.text      = element_text(size = FS$legend_text, color = "black"),
      legend.key       = element_blank(),
      legend.background = element_blank(),
      legend.key.size  = unit(10, "pt"),
      strip.background = element_rect(fill = "grey92", color = NA),
      strip.text       = element_text(size = FS$facet, face = "bold",
                                      color = "black"),
      plot.margin      = margin(4, 4, 4, 4)
    )
}

# ---- Triple-format export: PDF + TIFF + PNG --------------------------------
save_paper_fig <- function(plot, basename,
                           canvas = "SMALL",
                           out_dir = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/figures/paper_v2",
                           width = NULL, height = NULL) {
  pdf_dir <- file.path(out_dir, "PDF")
  tif_dir <- file.path(out_dir, "TIFF")
  png_dir <- file.path(out_dir, "PNG")
  dir.create(pdf_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(tif_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(png_dir, showWarnings = FALSE, recursive = TRUE)
  
  if (is.null(width) || is.null(height)) {
    dim <- CANVAS[[canvas]]
    width  <- as.numeric(dim["w"])
    height <- as.numeric(dim["h"])
  }
  
  pdf_path <- file.path(pdf_dir, paste0(basename, ".pdf"))
  tif_path <- file.path(tif_dir, paste0(basename, ".tif"))
  png_path <- file.path(png_dir, paste0(basename, ".png"))
  
  ggsave(pdf_path, plot, width = width, height = height,
         device = cairo_pdf, family = "Helvetica")
  ggsave(tif_path, plot, width = width, height = height,
         dpi = 300, compression = "lzw", device = "tiff")
  ggsave(png_path, plot, width = width, height = height,
         dpi = 300, device = "png", bg = "white")   # [restyle] 200->300, 与TIFF一致便于审图
  
  cat(sprintf("[OK] %-40s  %.2f x %.2f in  ->  PDF + TIFF + PNG\n",
              basename, width, height))
  invisible(list(pdf = pdf_path, tif = tif_path, png = png_path))
}

# ---- Data roots ------------------------------------------------------------
DATA_ROOT     <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
DEG_DIR_LIVER <- file.path(DATA_ROOT, "diff_stage1_results/deg_v4_refined/liver")
DEG_DIR_BRAIN <- file.path(DATA_ROOT, "diff_stage1_results/deg_v4_refined/brain")
DELIV_BRAIN_DIR <- file.path(DATA_ROOT,
  "SGS10_deliverables_2026-04-20/03_sc_DEG/brain_clean")
TF_DIR_BRAIN  <- file.path(DATA_ROOT,
  "dorothea_results_brain_clean/per_contrast")

# Load per-celltype DEG. Returns a normalized data frame with columns
# gene / avg_log2FC / pct1 / pct2 / p_val / p_val_adj (regardless of source).
#
# source: "diff_stage1" (default, liver + brain side) or "deliverables_brain"
# diff_stage1 columns:    gene, avg_log2FC, pct1, pct2, p_val, p_val_adj, sig
# deliverables columns:   cell_type, gene, avg_log2FC, pct.1, pct.2,
#                         p_val, p_val_adj, n_test, n_ctrl
load_deg <- function(organ, contrast, celltype, source = "diff_stage1") {
  if (source == "diff_stage1") {
    dir <- if (organ == "liver") DEG_DIR_LIVER else DEG_DIR_BRAIN
    fn  <- sprintf("%s_%s_%s.csv", organ, contrast, celltype)
    full <- file.path(dir, fn)
    if (!file.exists(full)) stop("DEG file not found: ", full)
    df <- read.csv(full, stringsAsFactors = FALSE)
    # already normalized: gene, avg_log2FC, pct1, pct2, p_val, p_val_adj
    return(df)
  } else if (source == "deliverables_brain") {
    # contrast like "Tumor_vs_Control" / "Tumor_vs_Blank" / "Blank_vs_Control"
    sub <- sprintf("brain_%s", contrast)
    fn  <- sprintf("%s_deg.csv", celltype)
    full <- file.path(DELIV_BRAIN_DIR, sub, fn)
    if (!file.exists(full)) stop("DEG file not found: ", full)
    df <- read.csv(full, stringsAsFactors = FALSE)
    # rename pct.1/pct.2 -> pct1/pct2, drop cell_type/n_test/n_ctrl
    if ("pct.1" %in% colnames(df)) names(df)[names(df) == "pct.1"] <- "pct1"
    if ("pct.2" %in% colnames(df)) names(df)[names(df) == "pct.2"] <- "pct2"
    keep <- c("gene", "avg_log2FC", "pct1", "pct2", "p_val", "p_val_adj")
    df <- df[, intersect(keep, colnames(df))]
    return(df)
  } else {
    stop("Unknown source: ", source)
  }
}

PALETTE <- list(
  up           = "#E04A4A",
  down         = "#3B82C4",
  ns           = "#BBBBBB",
  key_orange   = "#E66D50",
  liver        = "#D2691E",
  brain        = "#6A5ACD",
  baseline     = "#888888",
  tumor        = "#C0143B"
)

# Single-panel lollipop with hjust depending on sign of x.
build_lollipop <- function(df, genes,
                          x_col = "avg_log2FC", padj_col = "p_val_adj",
                          sig_thr = 0.05,
                          color_up = PALETTE$up, color_down = PALETTE$down,
                          color_ns = PALETTE$ns,
                          x_lab = NULL, show_legend = TRUE,
                          show_size_legend = TRUE,
                          padj_format = "%.1e") {
  if (!"gene" %in% colnames(df)) stop("df must have column 'gene'")
  sub <- df[df$gene %in% genes, , drop = FALSE]
  missing <- setdiff(genes, sub$gene)
  if (length(missing) > 0) {
    message("  [build_lollipop] missing: ", paste(missing, collapse=", "))
  }
  if (nrow(sub) == 0) stop("No genes matched.")
  
  sub$.x <- sub[[x_col]]
  sub$.padj <- sub[[padj_col]]
  sub <- sub[order(-sub$.x), , drop = FALSE]
  sub$gene <- factor(sub$gene, levels = rev(sub$gene))
  sub$is_sig <- sub$.padj < sig_thr
  sub$color_grp <- ifelse(!sub$is_sig, "ns",
                          ifelse(sub$.x >= 0, "up", "down"))
  sub$color_grp <- factor(sub$color_grp, levels = c("down", "ns", "up"))
  sub$neglog10_padj <- -log10(pmax(sub$.padj, 1e-300))
  
  hj <- ifelse(sub$.x >= 0, -0.15, 1.15)
  
  p <- ggplot(sub, aes(x = .x, y = gene)) +
    geom_segment(aes(xend = 0, yend = gene, color = color_grp),
                 linewidth = 0.7) +
    geom_point(aes(size = neglog10_padj, color = color_grp)) +
    geom_text(aes(label = sprintf(paste0("padj=", padj_format), .padj)),
              hjust = hj, size = 2.8, color = "grey25") +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.4) +
    scale_color_manual(
      values = c("down" = color_down, "ns" = color_ns, "up" = color_up),
      labels = c("down" = "Down (sig)", "ns" = "ns", "up" = "Up (sig)"),
      name = NULL, drop = FALSE,
      limits = c("up", "down", "ns")
    ) +
    scale_size_continuous(range = c(2.5, 7),
                          name = expression(-log[10]~"(padj)")) +
    scale_x_continuous(expand = expansion(mult = c(0.25, 0.30))) +
    labs(x = x_lab, y = NULL) +
    theme_paper() +
    theme(axis.text.y = element_text(size = FS$axis_text, face = "bold"))
  
  if (!show_legend)       p <- p + guides(color = "none")
  if (!show_size_legend)  p <- p + guides(size = "none")
  p
}

cat("paper_theme.R loaded.\n")
cat("  theme_paper()                                  - ggplot theme\n")
cat("  save_paper_fig(p, name)                        - PDF + TIFF + PNG\n")
cat("  load_deg(organ, contrast, celltype, source)    - load DEG\n")
cat("  build_lollipop(df, genes, ...)                 - standard lollipop\n")
cat("  CANVAS / FS / PALETTE                          - constants\n")
