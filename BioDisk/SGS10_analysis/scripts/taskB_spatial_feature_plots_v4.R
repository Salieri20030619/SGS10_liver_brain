#!/usr/bin/env Rscript
source("/mnt/BioDisk/wangmingyuan/SGS10_analysis/scripts/fig_save.R")  # [restyle] 三格式输出器
# =============================================================================
# taskB_spatial_feature_plots_v4.R
# =============================================================================
#
# v4 fixes:
#   1. FetchData mangles "mm10---Acta2" to "mm10.Acta2" etc.
#      Fix: use make.names() to predict the mangled name, then rename properly.
#      Better approach: GetAssayData(obj, layer="data")[gene_id, ] directly.
#   2. sprintf("%d", 17.5) fails. Use round() or "%g".
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
RCTD100_DIR  <- "spatial_results/rctd_full/per_sample"   # [restyle] Fig4 细胞身份改用 RCTD100
OUT_DIR <- "taskB/spatial_figures"
RESTYLE_DIR <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis/figures/restyle_v1"  # [restyle] 新图统一目录
dir.create(RESTYLE_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

LIVER_SAMPLES <- list(
  list(name = "Con_liver_1", rds = "Con_liver_1_seekspace_rctd.rds", group = "Control"),
  list(name = "Bla_liver_1", rds = "Bla_liver_1_seekspace_rctd.rds", group = "Blank"),
  list(name = "Liver_3",     rds = "Liver_3_seekspace_rctd.rds",     group = "Tumor"),
  list(name = "Liver_8",     rds = "Liver_8_seekspace_rctd.rds",     group = "Tumor"),
  list(name = "gan56",       rds = "gan56_seekspace_rctd.rds",       group = "Tumor")
)

BRAIN_SAMPLES <- list(
  list(name = "Con_brain_1", rds = "Con_brain_1_seekspace_rctd.rds", group = "Control"),
  list(name = "Bla_brain_1", rds = "Bla_brain_1_seekspace_rctd.rds", group = "Blank"),
  list(name = "Brain_3",     rds = "Brain_3_seekspace_rctd.rds",     group = "Tumor"),
  list(name = "Brain_8",     rds = "Brain_8_seekspace_rctd.rds",     group = "Tumor"),
  list(name = "nao56",       rds = "nao56_seekspace_rctd.rds",       group = "Tumor")
)

FIG_GROUPS <- list(
  Fig2a_L1_HSC_myofib = list(
    samples = LIVER_SAMPLES,
    genes = c("Acta2", "Serpine1", "Col1a1", "Lox", "Tagln"),
    title = "L1 HSC myofibroblast activation (liver, host zone)",
    zones_keep = c("host", "host_dom", "margin")
  ),
  Fig2b_L4_immune_vascular = list(
    samples = LIVER_SAMPLES,
    genes = c("Cd74", "Csf1r", "Mrc1", "Pecam1"),
    title = "L4 Liver immune-vascular hub (myeloid + endo)",
    zones_keep = c("host", "host_dom", "margin")
  ),
  Fig2c_L5_F2r = list(
    samples = LIVER_SAMPLES,
    genes = c("F2r"),
    title = "L5 F2r 20x receptor sensitization on HSC",
    zones_keep = c("host", "host_dom", "margin")
  ),
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

resolve_gene_names <- function(obj, genes) {
  feats <- rownames(obj@assays$RNA$data)
  sapply(genes, function(g) {
    if (g %in% feats) return(g)
    prefixed <- paste0("mm10---", g)
    if (prefixed %in% feats) return(prefixed)
    return(NA_character_)
  }, USE.NAMES = TRUE)
}

# =============================================================================
# Extract data frame: use GetAssayData directly (no name mangling)
# =============================================================================
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
  
  resolved <- resolve_gene_names(obj, genes)
  missing <- genes[is.na(resolved)]
  if (length(missing) > 0) {
    cat(sprintf("    WARN: not found in %s: %s\n", 
                sample_info$name, paste(missing, collapse=", ")))
  }
  resolved <- resolved[!is.na(resolved)]
  if (length(resolved) == 0) return(NULL)
  
  # Coords
  coords <- as.data.frame(Embeddings(obj, "spatial"))
  colnames(coords) <- c("x", "y")
  coords$barcode <- rownames(coords)
  
  # KEY FIX: Use GetAssayData directly, no FetchData name mangling
  # This gives us a dgCMatrix with EXACT rownames we asked for
  expr_mat <- GetAssayData(obj, assay = "RNA", layer = "data")[resolved, , drop = FALSE]
  expr_df <- as.data.frame(as.matrix(t(expr_mat)))
  # Strip mm10--- prefix from columns
  colnames(expr_df) <- sub("^mm10---", "", colnames(expr_df))
  expr_df$barcode <- rownames(expr_df)
  
  # Get metadata: rctd_dominant, zone
  meta_df <- obj@meta.data[, c("rctd_dominant", "zone"), drop = FALSE]
  meta_df$barcode <- rownames(meta_df)
  
  # Merge all
  df <- merge(coords, expr_df, by = "barcode")
  df <- merge(df, meta_df, by = "barcode")

  # [restyle] 叠加 RCTD100 主导细胞类型(Fig4 用), 按裸 barcode 对齐
  r100 <- file.path(RCTD100_DIR, paste0(sample_info$name, "_dominant_celltype.csv"))
  if (file.exists(r100)) {
    dom100 <- read.csv(r100, stringsAsFactors = FALSE)
    df <- merge(df, dom100[, c("barcode", "dominant_celltype")],
                by = "barcode", all.x = TRUE)
    names(df)[names(df) == "dominant_celltype"] <- "rctd100_dominant"
    .m <- sum(!is.na(df$rctd100_dominant))
    cat(sprintf("    [%s] RCTD100 主导类型匹配: %d / %d (%.1f%%)\n",
                sample_info$name, .m, nrow(df), 100 * .m / max(nrow(df), 1)))
  } else {
    df$rctd100_dominant <- NA_character_
    cat(sprintf("    [%s] 警告: 缺 RCTD100 dominant csv -> Fig4 该样本将全为 Other\n",
                sample_info$name))
  }

  df$sample <- sample_info$name
  df$group <- sample_info$group
  
  if (!is.null(zones_keep)) {
    df <- df[df$zone %in% zones_keep, ]
  }
  return(df)
}

# =============================================================================
# Single panel
# =============================================================================
plot_single_panel <- function(df, gene, show_title = TRUE, 
                               show_legend = FALSE, global_cap = NULL) {
  df$expr <- df[[gene]]
  zero_df <- df[!is.na(df$expr) & df$expr == 0, ]
  nz_df <- df[!is.na(df$expr) & df$expr > 0, ]
  
  cap <- if (!is.null(global_cap) && !is.na(global_cap)) global_cap else {
    if (nrow(nz_df) > 0) quantile(nz_df$expr, 0.99, na.rm = TRUE) else 1
  }
  if (nrow(nz_df) > 0) nz_df$expr_capped <- pmin(nz_df$expr, cap)
  
  sample_name <- unique(df$sample)[1]
  group_name <- unique(df$group)[1]
  
  p <- ggplot() +
    geom_point(data = zero_df, aes(x = x, y = y), 
               color = "grey90", size = 0.1, alpha = 0.4) +
    geom_point(data = nz_df, aes(x = x, y = y, color = expr_capped),
               size = 0.4, alpha = 0.9) +
    scale_color_gradient(low = "#FEE8C8", high = "#B30000",
                         name = gene, limits = c(0, cap)) +
    coord_fixed() +
    theme_void() +
    theme(
      legend.position = if (show_legend) "right" else "none",
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),  # [restyle] 9->12
      plot.margin = margin(1, 1, 1, 1)
    )
  
  if (show_title) {
    p <- p + labs(title = sprintf("%s (%s)", sample_name, group_name))
  }
  return(p)
}

# =============================================================================
# Figure group (grid layout)
# =============================================================================
plot_figure_group <- function(group_name, group_config) {
  cat(sprintf("\n========================================\n"))
  cat(sprintf("=== %s ===\n", group_name))
  cat(sprintf("========================================\n"))
  
  samples <- group_config$samples
  genes <- group_config$genes
  zones <- group_config$zones_keep
  sample_order <- sapply(samples, function(x) x$name)
  
  cat("\nLoading samples:\n")
  dfs <- list()
  for (s in samples) {
    df <- extract_plot_df(s, genes, zones_keep = zones)
    if (!is.null(df)) dfs[[s$name]] <- df
  }
  if (length(dfs) == 0) {
    cat("ERROR: no data\n"); return(NULL)
  }
  
  # DIAGNOSTIC: check what columns each df has
  cat("\nData check (columns per sample):\n")
  for (s_name in names(dfs)) {
    cat(sprintf("  %s: %d rows, genes present: %s\n",
                s_name, nrow(dfs[[s_name]]),
                paste(intersect(genes, colnames(dfs[[s_name]])), collapse=",")))
  }
  
  # Global cap per gene - check each sample independently, take max 0.99 quantile
  gene_caps <- sapply(genes, function(g) {
    all_nz <- c()
    for (s_name in names(dfs)) {
      df <- dfs[[s_name]]
      if (g %in% colnames(df)) {
        vals <- df[[g]]
        all_nz <- c(all_nz, vals[!is.na(vals) & vals > 0])
      }
    }
    if (length(all_nz) == 0) return(NA_real_)
    quantile(all_nz, 0.99, na.rm = TRUE)
  })
  
  for (g in genes) {
    cat(sprintf("  Gene cap %s: %.3f\n", g, gene_caps[g]))
  }
  
  panels <- list()
  valid_genes <- genes[!is.na(gene_caps)]
  if (length(valid_genes) == 0) {
    cat("ERROR: all gene caps are NA\n"); return(NULL)
  }
  
  for (gi in seq_along(valid_genes)) {
    g <- valid_genes[gi]
    for (si in seq_along(sample_order)) {
      s_name <- sample_order[si]
      if (is.null(dfs[[s_name]])) {
        panels[[length(panels) + 1]] <- patchwork::plot_spacer()
        next
      }
      df <- dfs[[s_name]]
      if (!g %in% colnames(df)) {
        panels[[length(panels) + 1]] <- patchwork::plot_spacer()
        next
      }
      
      show_title <- (gi == 1)
      show_legend <- (si == length(sample_order))
      
      p <- plot_single_panel(df, g,
                              show_title = show_title,
                              show_legend = show_legend,
                              global_cap = gene_caps[g])
      
      if (si == 1) {
        p <- p + labs(tag = g) +
          theme(plot.tag = element_text(size = 14, face = "bold.italic",  # [restyle] 11->14
                                         angle = 90, hjust = 0.5))
      }
      
      panels[[length(panels) + 1]] <- p
    }
    cat(sprintf("  Row %d: %s (cap=%.3f)\n", gi, g, gene_caps[g]))
  }
  
  if (length(panels) == 0) {
    cat("ERROR: no panels\n"); return(NULL)
  }
  
  final <- wrap_plots(panels, nrow = length(valid_genes), ncol = length(sample_order))  # [restyle] 大标题/副标题已移除
  
  # FIX sprintf %d bug: use round() for integer output
  w <- round(3.5 * length(sample_order))
  h <- round(max(4, 3 * length(valid_genes) + 0.5))
  
  pdf_out <- file.path(RESTYLE_DIR, sprintf("%s.pdf", group_name))
  png_out <- file.path(RESTYLE_DIR, sprintf("%s.png", group_name))
  
  cat(sprintf("Saving (%d x %d inch) ... ", w, h))
  tryCatch({
    save_fig_all(final, sub("\\.pdf$", "", pdf_out), width = w, height = h)  # [restyle] PDF+TIFF+PNG@300
    cat(sprintf("OK (PDF %.1f KB)\n", file.size(pdf_out) / 1024))
  }, error = function(e) {
    cat(sprintf("SAVE FAILED: %s\n", e$message))
  })
}

# =============================================================================
# Celltype map
# =============================================================================
plot_celltype_map <- function(organ_name, samples, highlight_ct = NULL) {
  cat(sprintf("\n=== Celltype map: %s ===\n", organ_name))
  sample_order <- sapply(samples, function(x) x$name)
  
  dfs <- list()
  for (s in samples) {
    df <- extract_plot_df(s, genes = "Acta2", zones_keep = NULL)
    if (!is.null(df)) dfs[[s$name]] <- df
  }
  if (length(dfs) == 0) return(NULL)
  
  if (!is.null(highlight_ct)) {
    used_cts <- c(highlight_ct, "Other")
    cols <- c(setNames(scales::hue_pal()(length(highlight_ct)), highlight_ct),
              "Other" = "grey90")
  } else {
    top_ct <- names(sort(table(unlist(lapply(dfs, function(d) as.character(d$rctd100_dominant)))),
                          decreasing = TRUE))[1:15]
    used_cts <- c(top_ct, "Other")
    cols <- c(setNames(scales::hue_pal()(length(top_ct)), top_ct),
              "Other" = "grey90")
  }
  
  panels <- list()
  for (si in seq_along(sample_order)) {
    s_name <- sample_order[si]
    df <- dfs[[s_name]]
    if (is.null(df)) next
    
    if (!is.null(highlight_ct)) {
      df$display <- ifelse(df$rctd100_dominant %in% highlight_ct,
                            as.character(df$rctd100_dominant), "Other")
    } else {
      df$display <- ifelse(df$rctd100_dominant %in% top_ct,
                            as.character(df$rctd100_dominant), "Other")
    }
    df$display[is.na(df$display)] <- "Other"
    df$display <- factor(df$display, levels = used_cts)
    
    other_df <- df[df$display == "Other", ]
    main_df <- df[df$display != "Other", ]
    
    group_name <- df$group[1]
    
    p <- ggplot() +
      geom_point(data = other_df, aes(x = x, y = y),
                 color = "grey95", size = 0.1, alpha = 0.3) +
      geom_point(data = main_df, aes(x = x, y = y, color = display),
                 size = 0.3, alpha = 0.9) +
      scale_color_manual(values = cols, name = "Cell type", drop = FALSE) +
      coord_fixed() +
      theme_void() +
      theme(
        legend.position = if (si == length(sample_order)) "right" else "none",
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),  # [restyle] 10->12
        plot.margin = margin(2, 2, 2, 2)
      ) +
      guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
      labs(title = sprintf("%s (%s)", s_name, group_name))
    
    panels[[length(panels) + 1]] <- p
  }
  
  final <- wrap_plots(panels, nrow = 1)  # [restyle] 大标题/副标题已移除
  
  w <- round(3.5 * length(sample_order) + 2)
  h <- 5
  
  out_base <- sprintf("Fig4_%s_celltype_map", organ_name)
  pdf_out <- file.path(RESTYLE_DIR, sprintf("%s.pdf", out_base))
  png_out <- file.path(RESTYLE_DIR, sprintf("%s.png", out_base))
  
  cat(sprintf("Saving (%d x %d inch) ... ", w, h))
  tryCatch({
    save_fig_all(final, sub("\\.pdf$", "", pdf_out), width = w, height = h)  # [restyle] PDF+TIFF+PNG@300
    cat(sprintf("OK (PDF %.1f KB)\n", file.size(pdf_out) / 1024))
  }, error = function(e) {
    cat(sprintf("SAVE FAILED: %s\n", e$message))
  })
}

# =============================================================================
# MAIN
# =============================================================================
cat(sprintf("\n######################################################\n"))
cat(sprintf("### Task B Spatial Feature Plots v4 ###\n"))
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
  plot_celltype_map("liver", LIVER_SAMPLES,
                     highlight_ct = c("HSC_Activated_ECM", "Macrophage_M1",
                                      "Kupffer", "LSEC_Mid", "cDC2"))
}, error = function(e) cat(sprintf("!!! liver celltype map: %s\n", e$message)))

tryCatch({
  plot_celltype_map("brain", BRAIN_SAMPLES,
                     highlight_ct = c("Microglia_Homeostatic", "Astrocyte_Resting",
                                      "Endothelial", "Pericytes", "CPECs",
                                      "Hpc_CA1"))
}, error = function(e) cat(sprintf("!!! brain celltype map: %s\n", e$message)))

cat(sprintf("\n######################################################\n"))
cat(sprintf("### DONE ###\n"))
cat(sprintf("######################################################\n"))
cat("\nAll outputs:\n")
list.files(OUT_DIR, pattern = "\\.pdf$", full.names = FALSE)

cat("\nFile sizes:\n")
pdfs <- list.files(OUT_DIR, pattern = "\\.pdf$", full.names = TRUE)
for (f in pdfs) {
  cat(sprintf("  %-40s  %8.1f KB\n", basename(f), file.size(f)/1024))
}
