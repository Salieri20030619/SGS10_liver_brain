#!/usr/bin/env Rscript
# =============================================================================
# taskB_spatial_feature_plots_brain_patch.R
# =============================================================================
#
# v2 patch - 只补画 brain 部分 (liver 部分已经成功)
# Bug fix: brain Seurat object 基因名没有 mm10--- 前缀 (pure_mouse),
#          liver 有前缀. 新版自动探测.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

PROJECT_ROOT <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
setwd(PROJECT_ROOT)

ENRICHED_DIR <- "spatial_results/rctd/enriched"
OUT_DIR <- "taskB/spatial_figures"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

BRAIN_SAMPLES <- list(
  list(name = "Con_brain_1", rds = "Con_brain_1_seekspace_rctd.rds", group = "Control"),
  list(name = "Bla_brain_1", rds = "Bla_brain_1_seekspace_rctd.rds", group = "Blank"),
  list(name = "Brain_3",     rds = "Brain_3_seekspace_rctd.rds",     group = "Tumor"),
  list(name = "Brain_8",     rds = "Brain_8_seekspace_rctd.rds",     group = "Tumor"),
  list(name = "nao56",       rds = "nao56_seekspace_rctd.rds",       group = "Tumor")
)

# Only brain figures
FIG_GROUPS <- list(
  Fig3a_N4_BBB = list(
    samples = BRAIN_SAMPLES,
    genes = c("Cldn5", "Pecam1"),
    title = "N4 Brain endothelial BBB tight junction",
    zones_keep = NULL
  ),
  Fig3b_X2_Thyroid = list(
    samples = BRAIN_SAMPLES,
    genes = c("Thra", "Ttr"),
    title = "X2 Thyroid axis brain (Thra + Ttr)",
    zones_keep = NULL
  ),
  Fig3c_N1_Microglia_CPECs_hub = list(
    samples = BRAIN_SAMPLES,
    genes = c("Smad3", "Apoe", "Psap", "Gas6"),
    title = "N1 Microglia state + CPECs broadcast hub",
    zones_keep = NULL
  )
)

.obj_cache <- new.env()

# =============================================================================
# Key fix: auto-detect gene name prefix
# =============================================================================
resolve_gene_names <- function(obj, genes) {
  feats <- rownames(obj@assays$RNA$data)
  
  # Try both with and without mm10--- prefix
  resolved <- sapply(genes, function(g) {
    # First try: pure name (brain case)
    if (g %in% feats) return(g)
    # Second try: mm10--- prefix (liver double-species case)
    prefixed <- paste0("mm10---", g)
    if (prefixed %in% feats) return(prefixed)
    # Not found
    return(NA_character_)
  }, USE.NAMES = TRUE)
  
  return(resolved)
}

extract_plot_df <- function(sample_info, genes, zones_keep = NULL) {
  key <- sample_info$name
  if (!exists(key, envir = .obj_cache)) {
    rds_path <- file.path(ENRICHED_DIR, sample_info$rds)
    cat(sprintf("  Loading %s ... ", basename(rds_path)))
    t0 <- Sys.time()
    obj <- readRDS(rds_path)
    .obj_cache[[key]] <- obj
    cat(sprintf("done (%.1fs, %d spots)\n",
                as.numeric(difftime(Sys.time(), t0, units = "secs")),
                ncol(obj)))
  }
  obj <- .obj_cache[[key]]
  
  # Resolve gene names (handle both prefix styles)
  resolved <- resolve_gene_names(obj, genes)
  missing <- genes[is.na(resolved)]
  if (length(missing) > 0) {
    cat(sprintf("    WARN: genes not found in %s: %s\n", 
                sample_info$name, paste(missing, collapse=", ")))
  }
  resolved <- resolved[!is.na(resolved)]
  if (length(resolved) == 0) {
    warning(sprintf("No genes found in %s", sample_info$name))
    return(NULL)
  }
  
  # Coords
  coords <- as.data.frame(Embeddings(obj, "spatial"))
  colnames(coords) <- c("x", "y")
  coords$barcode <- rownames(coords)
  
  # Fetch (using resolved names with or without prefix)
  fd <- FetchData(obj, vars = c(resolved, "rctd_dominant", "zone"))
  # Rename: strip mm10--- prefix to keep column name consistent across samples
  col_names <- colnames(fd)
  for (i in seq_along(resolved)) {
    clean_name <- sub("mm10---", "", resolved[i])
    col_names[col_names == resolved[i]] <- clean_name
  }
  colnames(fd) <- col_names
  fd$barcode <- rownames(fd)
  
  # Merge
  df <- merge(coords, fd, by = "barcode")
  df$sample <- sample_info$name
  df$group <- sample_info$group
  
  if (!is.null(zones_keep)) {
    df <- df[df$zone %in% zones_keep, ]
  }
  
  return(df)
}

# =============================================================================
# Plot functions (same as v1)
# =============================================================================
plot_gene_row <- function(combined_df, gene, sample_order) {
  combined_df$sample <- factor(combined_df$sample, levels = sample_order)
  combined_df$expr <- combined_df[[gene]]
  
  zero_df <- combined_df[combined_df$expr == 0, ]
  nz_df <- combined_df[combined_df$expr > 0, ]
  
  if (nrow(nz_df) > 0) {
    cap <- quantile(nz_df$expr, 0.99, na.rm = TRUE)
    nz_df$expr_capped <- pmin(nz_df$expr, cap)
  } else {
    cap <- 1
    nz_df$expr_capped <- nz_df$expr
  }
  
  p <- ggplot() +
    geom_point(data = zero_df, aes(x = x, y = y), 
               color = "grey90", size = 0.1, alpha = 0.4) +
    geom_point(data = nz_df, aes(x = x, y = y, color = expr_capped),
               size = 0.4, alpha = 0.9) +
    scale_color_gradient(low = "#FEE8C8", high = "#B30000",
                         name = gene, limits = c(0, cap)) +
    facet_wrap(~ sample, nrow = 1, scales = "free") +
    coord_fixed() +
    theme_void() +
    theme(
      strip.text = element_text(size = 10, face = "bold"),
      legend.position = "right",
      legend.title = element_text(size = 10, face = "italic"),
      legend.key.height = unit(0.8, "cm"),
      legend.key.width = unit(0.3, "cm"),
      plot.title = element_text(size = 11, face = "bold.italic", hjust = 0),
      plot.margin = margin(2, 2, 2, 2)
    ) +
    labs(title = gene)
  
  return(p)
}

plot_figure_group <- function(group_name, group_config) {
  cat(sprintf("\n========================================\n"))
  cat(sprintf("=== %s ===\n", group_name))
  cat(sprintf("=== %s ===\n", group_config$title))
  cat(sprintf("========================================\n"))
  
  samples <- group_config$samples
  genes <- group_config$genes
  zones <- group_config$zones_keep
  sample_order <- sapply(samples, function(x) x$name)
  
  cat("\nLoading samples + extracting data:\n")
  dfs <- list()
  for (s in samples) {
    df <- extract_plot_df(s, genes, zones_keep = zones)
    if (!is.null(df)) dfs[[s$name]] <- df
  }
  if (length(dfs) == 0) {
    cat("ERROR: no data for any sample!\n")
    return(NULL)
  }
  combined <- do.call(rbind, dfs)
  cat(sprintf("Total data: %d spots\n", nrow(combined)))
  
  plots <- lapply(genes, function(g) {
    if (!g %in% colnames(combined)) {
      cat(sprintf("  SKIP %s (not in data)\n", g))
      return(NULL)
    }
    cat(sprintf("  Plotting %s ... ", g))
    t0 <- Sys.time()
    p <- plot_gene_row(combined, g, sample_order)
    cat(sprintf("%.1fs\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    p
  })
  plots <- plots[!sapply(plots, is.null)]
  
  if (length(plots) == 0) {
    cat("ERROR: no plots generated\n")
    return(NULL)
  }
  
  final <- wrap_plots(plots, ncol = 1) +
    plot_annotation(
      title = group_config$title,
      subtitle = sprintf("Samples: %s | Zones: %s | %d genes x 5 samples",
                         paste(sample_order, collapse = ", "),
                         if (is.null(zones)) "all" else paste(zones, collapse=", "),
                         length(plots)),
      theme = theme(plot.title = element_text(size = 14, face = "bold"),
                    plot.subtitle = element_text(size = 10, color = "grey40"))
    )
  
  w <- 18
  h <- max(4, 3 * length(plots))
  
  pdf_out <- file.path(OUT_DIR, sprintf("%s.pdf", group_name))
  png_out <- file.path(OUT_DIR, sprintf("%s.png", group_name))
  
  cat(sprintf("Saving %s (%d x %d inch) ... ", pdf_out, w, h))
  ggsave(pdf_out, final, width = w, height = h, device = cairo_pdf)
  ggsave(png_out, final, width = w, height = h, dpi = 120)
  cat(sprintf("OK (PDF %.1f KB)\n", file.size(pdf_out) / 1024))
}

plot_celltype_map <- function(organ_name, samples, highlight_ct = NULL) {
  cat(sprintf("\n========================================\n"))
  cat(sprintf("=== Celltype map: %s ===\n", organ_name))
  cat(sprintf("========================================\n"))
  
  sample_order <- sapply(samples, function(x) x$name)
  
  dfs <- list()
  for (s in samples) {
    df <- extract_plot_df(s, genes = "Acta2", zones_keep = NULL)
    if (!is.null(df)) dfs[[s$name]] <- df
  }
  if (length(dfs) == 0) {
    cat("ERROR: no data\n")
    return(NULL)
  }
  combined <- do.call(rbind, dfs)
  combined$sample <- factor(combined$sample, levels = sample_order)
  
  if (!is.null(highlight_ct)) {
    combined$display <- ifelse(
      combined$rctd_dominant %in% highlight_ct,
      as.character(combined$rctd_dominant),
      "Other"
    )
    combined$display[is.na(combined$display) | combined$display == "NA"] <- "Other"
    combined$display <- factor(combined$display, levels = c(highlight_ct, "Other"))
    cols <- c(setNames(scales::hue_pal()(length(highlight_ct)), highlight_ct),
              "Other" = "grey90")
  } else {
    combined$display <- as.character(combined$rctd_dominant)
    combined$display[is.na(combined$display)] <- "unassigned"
    top_ct <- names(sort(table(combined$display), decreasing = TRUE))[1:15]
    combined$display <- ifelse(combined$display %in% top_ct, combined$display, "Other")
    combined$display <- factor(combined$display, levels = c(top_ct, "Other"))
    cols <- c(setNames(scales::hue_pal()(length(top_ct)), top_ct),
              "Other" = "grey90")
  }
  
  other_df <- combined[combined$display == "Other", ]
  main_df <- combined[combined$display != "Other", ]
  
  p <- ggplot() +
    geom_point(data = other_df, aes(x = x, y = y),
               color = "grey95", size = 0.1, alpha = 0.3) +
    geom_point(data = main_df, aes(x = x, y = y, color = display),
               size = 0.3, alpha = 0.9) +
    scale_color_manual(values = cols, name = "Cell type", drop = FALSE) +
    facet_wrap(~ sample, nrow = 1, scales = "free") +
    coord_fixed() +
    theme_void() +
    theme(
      strip.text = element_text(size = 10, face = "bold"),
      legend.position = "right",
      plot.title = element_text(size = 13, face = "bold")
    ) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    labs(title = sprintf("%s spatial cell type map", tools::toTitleCase(organ_name)),
         subtitle = sprintf("RCTD dominant celltype | %d spots | %d samples",
                            nrow(combined), length(samples)))
  
  out_base <- sprintf("Fig4_%s_celltype_map", organ_name)
  pdf_out <- file.path(OUT_DIR, sprintf("%s.pdf", out_base))
  png_out <- file.path(OUT_DIR, sprintf("%s.png", out_base))
  
  cat(sprintf("Saving %s ... ", pdf_out))
  ggsave(pdf_out, p, width = 22, height = 5, device = cairo_pdf)
  ggsave(png_out, p, width = 22, height = 5, dpi = 120)
  cat(sprintf("OK (PDF %.1f KB)\n", file.size(pdf_out) / 1024))
}

# =============================================================================
# MAIN - only brain
# =============================================================================
cat(sprintf("\n######################################################\n"))
cat(sprintf("### Task B Spatial Feature Plots v2 - BRAIN ONLY ###\n"))
cat(sprintf("### %s ###\n", Sys.time()))
cat(sprintf("######################################################\n"))

for (name in names(FIG_GROUPS)) {
  tryCatch({
    plot_figure_group(name, FIG_GROUPS[[name]])
  }, error = function(e) {
    cat(sprintf("\n!!! ERROR in %s: %s\n", name, e$message))
  })
}

tryCatch({
  plot_celltype_map("brain", BRAIN_SAMPLES,
                     highlight_ct = c("Microglia_Homeostatic", "Astrocyte_Resting",
                                      "Endothelial", "Pericytes", "CPECs",
                                      "Hpc_CA1"))
}, error = function(e) cat(sprintf("!!! %s\n", e$message)))

cat(sprintf("\n######################################################\n"))
cat(sprintf("### DONE ###\n"))
cat(sprintf("######################################################\n"))

list.files(OUT_DIR, pattern = "\\.pdf$", full.names = FALSE)
