#!/usr/bin/env Rscript
# ============================================================
# generate_niche_heatmap.R
# ============================================================
#
# 背景:
#   T1 niche 分析生成了每个 anchor cell type 的 _enrichment.csv
#   + _bar.pdf, 但缺一个汇总 heatmap (anchor × neighbor delta%).
#   本脚本从 all_niche_pairs.csv 或逐个 enrichment.csv 生成 heatmap.
#
# 输出:
#   taskB/niche_v2_10samples/liver_niche_heatmap_delta.pdf
#   taskB/niche_v2_10samples/brain_niche_heatmap_delta.pdf
#   taskB/niche_v2_10samples/<organ>_niche_delta_matrix.csv
#
# 用法:
#   cd ~/SGS10_analysis
#   Rscript scripts/generate_niche_heatmap.R
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
})

PROJECT_ROOT <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
setwd(PROJECT_ROOT)

NICHE_DIR <- "taskB/niche_v2_10samples"

# ============================================================
# Helper: auto-detect format and build anchor × neighbor matrix
# ============================================================
build_delta_matrix <- function(organ) {
  cat(sprintf("\n=== %s ===\n", toupper(organ)))
  
  # Try 1: unified all_niche_pairs.csv
  unified_path <- file.path(NICHE_DIR, sprintf("%s_all_niche_pairs.csv", organ))
  if (file.exists(unified_path)) {
    cat(sprintf("Using unified: %s\n", unified_path))
    df <- read.csv(unified_path, stringsAsFactors = FALSE)
    cat(sprintf("Columns: %s\n", paste(names(df), collapse = ", ")))
    cat(sprintf("Rows: %d\n", nrow(df)))
    cat(sprintf("Head:\n"))
    print(head(df, 3))
    
    # Auto-detect columns
    anchor_col <- intersect(names(df), c("anchor", "Anchor", "anchor_ct", "center"))[1]
    neighbor_col <- intersect(names(df), c("neighbor", "Neighbor", "neighbor_ct"))[1]
    delta_col <- intersect(names(df), c("delta_pct", "delta", "Delta", "delta_pct_tumor_vs_baseline"))[1]
    
    if (is.na(anchor_col) || is.na(neighbor_col) || is.na(delta_col)) {
      cat(sprintf("WARN: Unified schema mismatch. Need anchor/neighbor/delta columns. Falling back to per-anchor enrichment.csv\n"))
      df <- NULL
    }
  } else {
    df <- NULL
  }
  
  # Try 2: per-anchor enrichment.csv files
  if (is.null(df)) {
    patt <- sprintf("^%s_niche_.+_enrichment\\.csv$", organ)
    files <- list.files(NICHE_DIR, pattern = patt, full.names = TRUE)
    cat(sprintf("Found %d per-anchor enrichment CSVs\n", length(files)))
    
    if (length(files) == 0) {
      stop(sprintf("No enrichment data for %s", organ))
    }
    
    dfs <- lapply(files, function(f) {
      anchor <- sub(sprintf("^%s_niche_(.+)_enrichment\\.csv$", organ), "\\1", basename(f))
      x <- read.csv(f, stringsAsFactors = FALSE)
      cat(sprintf("  %s: cols=[%s], rows=%d\n", anchor, paste(names(x), collapse=","), nrow(x)))
      
      # Auto-detect columns
      neighbor_col <- intersect(names(x), c("neighbor", "Neighbor", "neighbor_ct", "celltype", "cell_type"))[1]
      # Delta column: try a few variants
      cand <- c("delta_pct","delta","Delta","delta_pct_tumor_vs_baseline","tumor_vs_baseline","diff","difference")
      delta_col <- intersect(names(x), cand)[1]
      
      if (is.na(neighbor_col) || is.na(delta_col)) {
        # Fallback: find first numeric column for delta
        num_cols <- names(x)[sapply(x, is.numeric)]
        if (length(num_cols) >= 2) {
          # Likely has tumor_pct and baseline_pct, delta = tumor - baseline
          tumor_col <- intersect(names(x), c("tumor_pct","tumor","tumor_mean"))[1]
          baseline_col <- intersect(names(x), c("baseline_pct","baseline","baseline_mean"))[1]
          if (!is.na(tumor_col) && !is.na(baseline_col)) {
            x$delta_pct <- x[[tumor_col]] - x[[baseline_col]]
            delta_col <- "delta_pct"
            cat(sprintf("    -> derived delta_pct = %s - %s\n", tumor_col, baseline_col))
          }
        }
        if (is.na(neighbor_col)) neighbor_col <- names(x)[1]
      }
      
      if (is.na(delta_col)) {
        cat(sprintf("    !! cannot determine delta column, skipping %s\n", anchor))
        return(NULL)
      }
      
      data.frame(
        anchor = anchor,
        neighbor = x[[neighbor_col]],
        delta = x[[delta_col]],
        stringsAsFactors = FALSE
      )
    })
    dfs <- dfs[!sapply(dfs, is.null)]
    df <- do.call(rbind, dfs)
    anchor_col <- "anchor"
    neighbor_col <- "neighbor"
    delta_col <- "delta"
  }
  
  if (is.null(df) || nrow(df) == 0) {
    stop(sprintf("Empty data for %s", organ))
  }
  
  cat(sprintf("\nFinal data: %d rows, %d unique anchors, %d unique neighbors\n",
              nrow(df), length(unique(df[[anchor_col]])), length(unique(df[[neighbor_col]]))))
  
  # Save long format CSV
  long_out <- file.path(NICHE_DIR, sprintf("%s_niche_delta_long.csv", organ))
  write.csv(df, long_out, row.names = FALSE)
  cat(sprintf("Wrote long CSV: %s\n", long_out))
  
  # Build matrix
  mat <- dcast(df, paste(anchor_col, "~", neighbor_col), value.var = delta_col, fun.aggregate = mean, fill = 0)
  rownames(mat) <- mat[[1]]
  mat[[1]] <- NULL
  mat <- as.matrix(mat)
  
  # Save matrix CSV
  mat_out <- file.path(NICHE_DIR, sprintf("%s_niche_delta_matrix.csv", organ))
  write.csv(mat, mat_out, row.names = TRUE)
  cat(sprintf("Wrote matrix CSV: %s  (%d × %d)\n", mat_out, nrow(mat), ncol(mat)))
  
  return(list(df = df, mat = mat, anchor_col = anchor_col, 
              neighbor_col = neighbor_col, delta_col = delta_col))
}

# ============================================================
# Heatmap plot function
# ============================================================
plot_heatmap <- function(data_obj, organ) {
  df <- data_obj$df
  anchor_col <- data_obj$anchor_col
  neighbor_col <- data_obj$neighbor_col
  delta_col <- data_obj$delta_col
  
  # Ensure factor order: put anchors with biggest variance first
  anchor_var <- tapply(abs(df[[delta_col]]), df[[anchor_col]], max, na.rm = TRUE)
  anchor_order <- names(sort(anchor_var, decreasing = TRUE))
  df[[anchor_col]] <- factor(df[[anchor_col]], levels = anchor_order)
  
  neighbor_max <- tapply(abs(df[[delta_col]]), df[[neighbor_col]], max, na.rm = TRUE)
  neighbor_order <- names(sort(neighbor_max, decreasing = TRUE))
  df[[neighbor_col]] <- factor(df[[neighbor_col]], levels = neighbor_order)
  
  # Cap extreme values for better color contrast
  vals <- df[[delta_col]]
  lim <- max(abs(quantile(vals, c(0.02, 0.98), na.rm = TRUE)))
  df$delta_capped <- pmin(pmax(vals, -lim), lim)
  
  n_anchor <- length(unique(df[[anchor_col]]))
  n_neighbor <- length(unique(df[[neighbor_col]]))
  
  # Dynamic size
  w <- max(7, min(16, 3 + n_neighbor * 0.35))
  h <- max(5, min(12, 2 + n_anchor * 0.35))
  
  p <- ggplot(df, aes_string(x = neighbor_col, y = anchor_col, fill = "delta_capped")) +
    geom_tile(color = "white", size = 0.3) +
    geom_text(aes_string(label = sprintf("sprintf('%%+.1f', %s)", delta_col)),
              size = 2.5, color = "black") +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, limit = c(-lim, lim),
      name = "Delta % (tumor - baseline)"
    ) +
    labs(
      title = sprintf("%s niche composition shift (tumor vs baseline)", 
                      tools::toTitleCase(organ)),
      subtitle = sprintf("Rows: anchor cell type | Columns: neighbor cell type | %d × %d",
                         n_anchor, n_neighbor),
      x = "Neighbor cell type",
      y = "Anchor cell type"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9),
      axis.text.y = element_text(size = 9),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", size = 13),
      legend.position = "right"
    )
  
  # Save PDF + PNG
  pdf_out <- file.path(NICHE_DIR, sprintf("%s_niche_heatmap_delta.pdf", organ))
  png_out <- file.path(NICHE_DIR, sprintf("%s_niche_heatmap_delta.png", organ))
  
  ggsave(pdf_out, p, width = w, height = h, device = cairo_pdf)
  ggsave(png_out, p, width = w, height = h, dpi = 150, device = "png")
  
  cat(sprintf("Wrote heatmap PDF: %s  (%.1f × %.1f inch)\n", pdf_out, w, h))
  cat(sprintf("Wrote heatmap PNG: %s\n", png_out))
  
  return(p)
}

# ============================================================
# Main
# ============================================================
for (organ in c("liver", "brain")) {
  cat(sprintf("\n##############################################\n"))
  cat(sprintf("### %s niche heatmap ###\n", toupper(organ)))
  cat(sprintf("##############################################\n"))
  
  tryCatch({
    data_obj <- build_delta_matrix(organ)
    plot_heatmap(data_obj, organ)
  }, error = function(e) {
    cat(sprintf("\n!!! ERROR for %s: %s\n", organ, e$message))
  })
}

cat(sprintf("\n\n==================================================\n"))
cat(sprintf("=== DONE ===\n"))
cat(sprintf("==================================================\n"))
cat(sprintf("Output files:\n"))
list.files(NICHE_DIR, pattern = "niche_heatmap_delta|niche_delta_matrix|niche_delta_long", full.names = TRUE)
