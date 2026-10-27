#!/usr/bin/env Rscript
# ============================================================================
# fig20_21_umap_paper.R
# Liver + Brain UMAP cell atlas, paper-grade Batch 7.
# Based on paper_umap_grouped_v3_13_replot.R (v3.13 final, 10-row packed legend,
# UMAP arrow axes). Changes from v3.13:
#   1. No plot.title (already removed in v3.13)
#   2. Slightly larger fonts for paper:
#        idx label 2.6 -> 2.9
#        supergroup repel 4.5 -> 5.0
#        UMAP1/UMAP2 arrow text 3.2 -> 3.8 bold
#        legend supergroup 3.2 -> 3.6 bold
#        legend idx 2.8 -> 3.0
#        legend subtype text 2.85 -> 3.1
#   3. Output: paper_v2/{PDF,TIFF,PNG} with fig20/fig21 names
#   4. Triple-format (cairo_pdf + tiff 300dpi LZW + png 200dpi)
# Layout algorithm (build_legend_layout v3.13) preserved verbatim.
# ============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(grid)
})

source("/mnt/BioDisk/wangmingyuan/SGS10_analysis/scripts/paper_theme.R")

UMAP_RDS_DIR <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis/umap_paper_grade"
OUT_DIR      <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis/figures/paper_v2"
stopifnot(dir.exists(UMAP_RDS_DIR))
for (sub in c("PDF","TIFF","PNG")) {
  dir.create(file.path(OUT_DIR, sub), showWarnings = FALSE, recursive = TRUE)
}

cat("R version:", R.version.string, "\n\n")

# ============================================================================
# Group map / Color (verbatim from v3.13)
# ============================================================================
LIVER_GROUP_MAP <- list(
  Hep         = c("Hep_Periportal","Hep_Perivenous","Hep_Proliferating"),
  LSEC        = c("LSEC","LSEC_Pericentral","LSEC_Periportal","LVEC"),
  HSC         = c("HSC_Quiescent","HSC_Activated_ECM","HSC_Activated_Chemo",
                  "HSC_Activated_Remodel","HSC_Proliferating"),
  KC          = c("Kupffer_KC"),
  Mac         = c("Macrophage_M1","Macrophage_M2","Macrophage_MDM","Macrophage_NAM"),
  DC          = c("DC_cDC1","DC_cDC2","DC_pDC"),
  B           = c("B_Memory","B_Naive","B_Plasma"),
  T_NK        = c("T_CD4","T_CD8","T_NKT","T_Treg","NK"),
  Granulocyte = c("Neutrophil","Basophil"),
  Cholangio   = c("Cholangiocyte"),
  Vasc        = c("Pericyte"),
  Erythrocyte = c("Erythrocyte"),
  LCM         = c("LCM")
)
LIVER_GROUP_COLOR <- c(
  Hep="#FFD700", LSEC="#1E90FF", HSC="#8B4513", KC="#FF6347",
  Mac="#FF1493", DC="#FF8C00", B="#9400D3", T_NK="#00CED1",
  Granulocyte="#32CD32", Cholangio="#DC143C", Vasc="#87CEEB",
  Erythrocyte="#8B0000", LCM="#FFB6C1"
)
BRAIN_GROUP_MAP <- list(
  Ex_Cortex   = c("Excitatory_L23","Excitatory_L45","Excitatory_L6"),
  Inh_Cortex  = c("Inhibitory_Meis2","Inhibitory_Pvalb","Inhibitory_Sst","Inhibitory_Vip"),
  Hpc         = c("Hpc_CA1","Hpc_CA3","Hpc_DG"),
  Oligo       = c("Oligo_Mature","OPCs"),
  Astro       = c("Astrocyte_Resting","Astrocyte_Gfap"),
  Microglia   = c("Microglia_Homeostatic","Microglia_M1","Microglia_M2","Microglia_DAM_strict"),
  Endo        = c("Endothelial"),
  Pericyte    = c("Pericytes"),
  VLMC        = c("VLMC"),
  CPEC        = c("CPECs"),
  Ependymal   = c("Ependymal")
)
BRAIN_GROUP_COLOR <- c(
  Ex_Cortex="#228B22", Inh_Cortex="#9ACD32", Hpc="#DC143C",
  Oligo="#FFD700", Astro="#8B4513", Microglia="#FF1493",
  Endo="#1E90FF", Pericyte="#87CEEB", VLMC="#9400D3",
  CPEC="#FF8C00", Ependymal="#FFB6C1"
)

# ============================================================================
# Helpers (verbatim from v3.13)
# ============================================================================
shade_palette <- function(base_hex, n) {
  if (n == 1) return(base_hex)
  rgb_int <- col2rgb(base_hex)
  hsv_mat <- rgb2hsv(rgb_int)
  h <- hsv_mat["h", 1]; s <- hsv_mat["s", 1]; v <- hsv_mat["v", 1]
  v_seq <- seq(min(v * 1.10, 0.95), max(v * 0.45, 0.40), length.out = n)
  s_seq <- seq(max(s * 0.95, 0.85), 1.00, length.out = n)
  cols <- mapply(function(vv, ss) hsv(h, ss, vv), v_seq, s_seq)
  unname(cols)
}
build_subtype_palette <- function(group_map, group_color) {
  pal <- c()
  for (g in names(group_map)) {
    subs <- group_map[[g]]
    cols <- shade_palette(group_color[g], length(subs))
    names(cols) <- subs
    pal <- c(pal, cols)
  }
  pal
}

# ============================================================================
# v3.13 LEGEND PACKING (verbatim, do not modify)
# ============================================================================
build_legend_layout <- function(group_map, cols_per_row = 4L) {
  group_names <- names(group_map)
  group_sizes <- sapply(group_map, length)
  n_groups <- length(group_names)

  n_trailing_singles <- 0L
  for (i in n_groups:1) {
    if (group_sizes[i] == 1L) {
      n_trailing_singles <- n_trailing_singles + 1L
    } else {
      break
    }
  }
  reserve_for_full <- (n_trailing_singles %/% cols_per_row) * cols_per_row
  reserve_start_idx <- n_groups - reserve_for_full + 1L
  if (reserve_for_full == 0L) reserve_start_idx <- n_groups + 1L

  cat(sprintf("  [packing] trailing singles: %d, reserve %d for full row(s) from group #%d\n",
              n_trailing_singles, reserve_for_full, reserve_start_idx))

  rows_layout <- list(list())
  global_idx <- 0L
  group_first_used <- character(0)

  add_cell <- function(gname, sub) {
    global_idx <<- global_idx + 1L
    is_first <- !(gname %in% group_first_used)
    if (is_first) group_first_used <<- c(group_first_used, gname)
    cell <- list(gname = gname, subtype = sub, idx = global_idx,
                 is_first_of_group = is_first)
    cur_row <- rows_layout[[length(rows_layout)]]
    cur_row[[length(cur_row) + 1L]] <- cell
    rows_layout[[length(rows_layout)]] <<- cur_row
  }
  row_used <- function() length(rows_layout[[length(rows_layout)]])
  new_row <- function() rows_layout[[length(rows_layout) + 1L]] <<- list()

  for (i in seq_len(reserve_start_idx - 1L)) {
    gname <- group_names[i]
    subs <- group_map[[gname]]
    n <- length(subs)
    cur <- row_used()
    rem <- cols_per_row - cur

    if (cur > 0L && n <= rem) {
      for (s in subs) add_cell(gname, s)
    } else {
      if (cur > 0L) new_row()
      remaining <- n
      placed <- 0L
      while (remaining > 0L) {
        if (row_used() >= cols_per_row) new_row()
        take <- min(cols_per_row - row_used(), remaining)
        for (k in seq_len(take)) {
          placed <- placed + 1L
          add_cell(gname, subs[placed])
        }
        remaining <- remaining - take
      }
    }
  }

  if (reserve_for_full > 0L) {
    if (row_used() > 0L) new_row()
    for (i in reserve_start_idx:n_groups) {
      if (row_used() >= cols_per_row) new_row()
      gname <- group_names[i]
      add_cell(gname, group_map[[gname]][1])
    }
  }

  rows_df <- list()
  for (ri in seq_along(rows_layout)) {
    row <- rows_layout[[ri]]
    for (ci in seq_along(row)) {
      cell <- row[[ci]]
      rows_df[[length(rows_df) + 1L]] <- data.frame(
        subtype = cell$subtype,
        supergroup = cell$gname,
        idx = cell$idx,
        col_in_group = ci,
        phy_row = ri - 1L,
        is_first_of_group = cell$is_first_of_group,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows_df)
}

cat("[unit test] shade_palette: ")
test_red <- shade_palette("#DC143C", 3); cat(paste(test_red, collapse=" "), "\n")
if (any(grepl("^#0[0-2]0[0-2]0[0-2]$", test_red))) stop("[FATAL] near-black")
cat("[unit test] PASSED\n\n")

# ============================================================================
# PLOT FUNCTION — paper-grade (slightly bigger fonts, no title, arrow axes)
# ============================================================================
plot_paper_umap <- function(obj, group_map, group_color, organ_name,
                            cols_per_row = 4L) {
  pal <- build_subtype_palette(group_map, group_color)
  cat(sprintf("  [diag] palette: %d entries\n", length(pal)))

  obj$celltype_refined <- as.character(obj$celltype_refined)
  if (!"supergroup" %in% colnames(obj@meta.data)) {
    flat <- stack(group_map); colnames(flat) <- c("subtype","supergroup")
    flat$subtype <- as.character(flat$subtype)
    flat$supergroup <- as.character(flat$supergroup)
    obj$supergroup <- flat$supergroup[match(
      as.character(obj$celltype_refined), flat$subtype)]
  }

  emb <- as.data.frame(Embeddings(obj, "umap"))
  colnames(emb) <- c("UMAP1","UMAP2")
  emb$subtype    <- as.character(obj$celltype_refined)
  emb$supergroup <- as.character(obj$supergroup)
  emb <- emb[!is.na(emb$subtype) & emb$subtype %in% names(pal), ]
  emb$subtype <- factor(emb$subtype, levels = names(pal))

  centroids_g <- emb %>%
    group_by(supergroup) %>%
    summarise(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2),
              n = n(), .groups = "drop")

  centroids_s <- emb %>%
    mutate(subtype = as.character(subtype)) %>%
    group_by(subtype, supergroup) %>%
    summarise(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2),
              n = n(), .groups = "drop")

  emb_ordered <- emb %>% group_by(subtype) %>%
    mutate(n_in_group = n()) %>% ungroup() %>%
    arrange(desc(n_in_group))

  legend_layout <- build_legend_layout(group_map, cols_per_row)
  centroids_s$idx <- legend_layout$idx[match(centroids_s$subtype,
                                              legend_layout$subtype)]

  x_range <- range(emb$UMAP1)
  y_range <- range(emb$UMAP2)
  x_span <- diff(x_range)
  y_span <- diff(y_range)
  arrow_len_x <- x_span * 0.20
  arrow_len_y <- y_span * 0.20
  origin_x <- x_range[1] - x_span * 0.02
  origin_y <- y_range[1] - y_span * 0.02

  # ---- MAIN UMAP (paper-grade font bumps) ----
  p <- ggplot(emb_ordered, aes(UMAP1, UMAP2, color = subtype)) +
    geom_point(size = 0.25, alpha = 0.95, stroke = 0) +
    scale_color_manual(values = pal,
                       breaks = names(pal),
                       limits = names(pal),
                       drop = FALSE,
                       na.value = "#FF00FF",
                       guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.02))) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.02))) +
    geom_label(data = centroids_s,
               aes(x = UMAP1, y = UMAP2, label = idx, fill = supergroup),
               color = "white", fontface = "bold",
               size = 2.9, label.size = 0,
               label.padding = unit(0.10, "lines"),
               label.r = unit(0.3, "lines"),
               show.legend = FALSE,
               inherit.aes = FALSE) +
    scale_fill_manual(values = group_color) +
    ggrepel::geom_text_repel(
      data = centroids_g,
      aes(x = UMAP1, y = UMAP2, label = supergroup),
      color = "black", size = 5.0, fontface = "bold",
      bg.color = "white", bg.r = 0.18,
      box.padding = 0.9, point.padding = 0.6,
      force = 3.5, force_pull = 0.3,
      min.segment.length = 0,
      segment.color = "grey60",
      segment.size = 0.35,
      max.overlaps = Inf,
      seed = 42,
      inherit.aes = FALSE
    ) +
    # Arrow axes (v3.13 paper-grade signature)
    annotate("segment",
             x = origin_x, xend = origin_x + arrow_len_x,
             y = origin_y, yend = origin_y,
             arrow = arrow(length = unit(0.20, "cm"), type = "closed"),
             color = "black", linewidth = 0.6) +
    annotate("segment",
             x = origin_x, xend = origin_x,
             y = origin_y, yend = origin_y + arrow_len_y,
             arrow = arrow(length = unit(0.20, "cm"), type = "closed"),
             color = "black", linewidth = 0.6) +
    annotate("text",
             x = origin_x + arrow_len_x * 0.5,
             y = origin_y - y_span * 0.025,
             label = "UMAP1", hjust = 0.5, vjust = 1,
             size = 3.8, fontface = "bold") +
    annotate("text",
             x = origin_x - x_span * 0.020,
             y = origin_y + arrow_len_y * 0.5,
             label = "UMAP2", angle = 90, hjust = 0.5, vjust = 1,
             size = 3.8, fontface = "bold") +
    theme_classic(base_size = 11) +
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.margin = margin(5, 5, 2, 5)
    )

  # ---- LEGEND (slightly bigger fonts) ----
  y_step <- 1.3
  x_step <- 4.7
  idx_x_start <- 1.5
  txt_offset <- 0.45
  gname_offset <- 0.50

  leg_data <- legend_layout
  leg_data$x_idx <- idx_x_start + (leg_data$col_in_group - 1L) * x_step
  leg_data$x_txt <- leg_data$x_idx + txt_offset
  leg_data$x_gname <- leg_data$x_idx - gname_offset
  leg_data$y_pos <- -leg_data$phy_row * y_step

  group_label_data <- leg_data[leg_data$is_first_of_group, ]

  cat(sprintf("  [diag] legend rows: %d\n", max(leg_data$phy_row) + 1L))

  total_rows <- max(leg_data$phy_row) + 1L
  max_subtype_len <- max(nchar(leg_data$subtype))
  xlim_right <- idx_x_start + (cols_per_row - 1L) * x_step + txt_offset +
                max_subtype_len * 0.30 + 0.5
  xlim_left <- -1.5

  p_leg <- ggplot() +
    geom_text(data = group_label_data,
              aes(x = x_gname, y = y_pos, label = supergroup),
              hjust = 1, fontface = "bold", size = 3.6, color = "black") +
    geom_label(data = leg_data,
               aes(x = x_idx, y = y_pos, label = idx, fill = supergroup),
               color = "white", fontface = "bold", size = 3.0,
               label.size = 0, label.padding = unit(0.12, "lines")) +
    geom_text(data = leg_data,
              aes(x = x_txt, y = y_pos, label = subtype),
              hjust = 0, size = 3.1) +
    scale_fill_manual(values = group_color, guide = "none") +
    scale_x_continuous(limits = c(xlim_left, xlim_right),
                       expand = expansion(0)) +
    scale_y_continuous(limits = c(-(total_rows - 0.5) * y_step,
                                   0.5 * y_step),
                       expand = expansion(0)) +
    coord_cartesian(clip = "off") +
    theme_void() +
    theme(plot.margin = margin(8, 8, 8, 4))

  combined <- p / p_leg + plot_layout(heights = c(3.3, 1))
  return(list(main = p, legend = p_leg, combined = combined,
              centroids_s = centroids_s, leg_data = leg_data))
}

# ============================================================================
# Triple-format export
# ============================================================================
save_umap_triple <- function(plot_obj, basename, w = 10, h = 12) {
  pdf_path <- file.path(OUT_DIR, "PDF",  paste0(basename, ".pdf"))
  tif_path <- file.path(OUT_DIR, "TIFF", paste0(basename, ".tif"))
  png_path <- file.path(OUT_DIR, "PNG",  paste0(basename, ".png"))

  ggsave(pdf_path, plot_obj, width = w, height = h,
         device = cairo_pdf, family = "Helvetica")
  ggsave(tif_path, plot_obj, width = w, height = h,
         dpi = 300, compression = "lzw", device = "tiff")
  ggsave(png_path, plot_obj, width = w, height = h,
         dpi = 200, device = "png", bg = "white")

  cat(sprintf("[OK] %-30s  %.0f x %.0f in  ->  PDF + TIFF + PNG\n",
              basename, w, h))
}

# ============================================================================
# Run
# ============================================================================
cat("##############################################\n")
cat("# fig20 + fig21 UMAP paper-grade (based on v3.13)\n")
cat(sprintf("# Started: %s\n", Sys.time()))
cat("##############################################\n\n")

cat("===== Liver layout preview =====\n")
liver_lo <- build_legend_layout(LIVER_GROUP_MAP, 4L)
cat(sprintf("Liver: %d entries, %d rows\n",
            nrow(liver_lo), max(liver_lo$phy_row) + 1L))

cat("\n===== Brain layout preview =====\n")
brain_lo <- build_legend_layout(BRAIN_GROUP_MAP, 4L)
cat(sprintf("Brain: %d entries, %d rows\n\n",
            nrow(brain_lo), max(brain_lo$phy_row) + 1L))

# ---- LIVER (fig20) ----
cat("[fig20 Liver] loading rds...\n")
liver_obj <- readRDS(file.path(UMAP_RDS_DIR, "liver_integrated.rds"))
cat(sprintf("  %d cells, %d samples\n",
            ncol(liver_obj), length(unique(liver_obj$sample))))
liver_plot <- plot_paper_umap(liver_obj, LIVER_GROUP_MAP,
                              LIVER_GROUP_COLOR, "Liver", 4L)
save_umap_triple(liver_plot$combined, "fig20_liver_umap",
                 w = 10, h = 12)
write.csv(liver_plot$centroids_s,
          file.path(UMAP_RDS_DIR, "liver_subtype_centroids_paper.csv"),
          row.names = FALSE)
cat("[fig20 Liver] done.\n\n")
rm(liver_obj, liver_plot); gc(verbose = FALSE)

# ---- BRAIN (fig21) ----
cat("[fig21 Brain] loading rds...\n")
brain_obj <- readRDS(file.path(UMAP_RDS_DIR, "brain_integrated.rds"))
cat(sprintf("  %d cells, %d samples\n",
            ncol(brain_obj), length(unique(brain_obj$sample))))
brain_plot <- plot_paper_umap(brain_obj, BRAIN_GROUP_MAP,
                              BRAIN_GROUP_COLOR, "Brain", 4L)
save_umap_triple(brain_plot$combined, "fig21_brain_umap",
                 w = 10, h = 12)
write.csv(brain_plot$centroids_s,
          file.path(UMAP_RDS_DIR, "brain_subtype_centroids_paper.csv"),
          row.names = FALSE)
cat("[fig21 Brain] done.\n\n")

cat(sprintf("[FINISHED] %s\n", Sys.time()))
