# =============================================================
# P0_split_human_mouse_genes_v2.R
# Purpose: 拆 SeekSpace 人鼠混合官方 RDS 为 human / mouse 两套 RDS
#
# 关键变化 vs v1：
#   - 坐标从 RCTD per_sample RDS 读（官方 RDS 没存坐标）
#   - 保留所有 spot（不按 human_frac 过滤）
#   - meta.data 构造用显式 rownames 赋值，避免 recycle bug
# =============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

OUT_DIR <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis/liver_species_split_v2"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Each sample: official RDS + RCTD RDS (for coords)
SAMPLES <- list(
  gan56 = list(
    official = "/mnt/BioData/bio/SGS10例/SGS24070264-gan56/data/rds/24070264_gan_56_expression/gan_56.rds",
    rctd     = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/spatial_results/rctd/per_sample/gan56_rctd_results.rds"
  ),
  Liver_3 = list(
    official = "/mnt/BioData/bio/SGS10例/SGS25090559-8例样本数据/data/rds/25090559_Liver_3_expression/Liver_3.rds",
    rctd     = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/spatial_results/rctd/per_sample/Liver_3_rctd_results.rds"
  ),
  Liver_8 = list(
    official = "/mnt/BioData/bio/SGS10例/SGS25090559-8例样本数据/data/rds/25090559_Liver_8_expression/Liver_8.rds",
    rctd     = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/spatial_results/rctd/per_sample/Liver_8_rctd_results.rds"
  )
)

log_msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

# -------------------------------------------------------------
# Extract coords from RCTD RDS
# -------------------------------------------------------------
get_coords_from_rctd <- function(rctd_path) {
  if (!file.exists(rctd_path)) {
    log_msg("  RCTD RDS not found: ", rctd_path)
    return(NULL)
  }
  rctd <- readRDS(rctd_path)
  # RCTD object structure: spot coordinates are in @spatialRNA@coords
  coords <- NULL
  if (!is.null(rctd@spatialRNA) && !is.null(rctd@spatialRNA@coords)) {
    coords <- rctd@spatialRNA@coords
    log_msg("  Coords from RCTD @spatialRNA@coords: ", nrow(coords), " spots")
  }
  rm(rctd); gc(verbose = FALSE)
  coords
}

# -------------------------------------------------------------
# Main split function
# -------------------------------------------------------------
split_one_sample <- function(sample_name, paths) {
  log_msg("==================================================")
  log_msg("Processing: ", sample_name)
  log_msg("==================================================")
  
  # ---- 1. Load official RDS ----
  log_msg("Loading official RDS: ", paths$official)
  obj <- readRDS(paths$official)
  
  counts <- tryCatch(
    GetAssayData(obj, assay = "RNA", layer = "counts"),
    error = function(e) GetAssayData(obj, assay = "RNA", slot = "counts")
  )
  genes <- rownames(counts)
  spot_bc <- colnames(counts)
  n_spots <- length(spot_bc)
  log_msg("Matrix: ", length(genes), " genes x ", n_spots, " spots")
  
  # ---- 2. Load image size from @misc ----
  size_x <- NA_real_; size_y <- NA_real_
  if (!is.null(obj@misc$info) && length(obj@misc$info) > 0) {
    info <- obj@misc$info[[1]]
    if (!is.null(info$size_x)) size_x <- as.numeric(info$size_x)
    if (!is.null(info$size_y)) size_y <- as.numeric(info$size_y)
  }
  log_msg("Image size: size_x=", size_x, ", size_y=", size_y)
  
  # ---- 3. Split genes by prefix ----
  human_idx <- grepl("^GRCh38-", genes)
  mouse_idx <- grepl("^mm10-", genes)
  log_msg("Human genes: ", sum(human_idx), ", Mouse genes: ", sum(mouse_idx))
  stopifnot(sum(human_idx) > 0, sum(mouse_idx) > 0)
  
  # ---- 4. Per-spot stats (before stripping prefix) ----
  h_per_spot <- Matrix::colSums(counts[human_idx, , drop = FALSE])
  m_per_spot <- Matrix::colSums(counts[mouse_idx, , drop = FALSE])
  h_frac     <- h_per_spot / pmax(h_per_spot + m_per_spot, 1)
  
  names(h_per_spot) <- spot_bc
  names(m_per_spot) <- spot_bc
  names(h_frac)     <- spot_bc
  
  log_msg("Per-spot human fraction: median=", round(median(h_frac), 3),
          ", mean=", round(mean(h_frac), 3))
  
  # ---- 5. Load coords from RCTD ----
  log_msg("Loading coords from RCTD RDS: ", paths$rctd)
  rctd_coords <- get_coords_from_rctd(paths$rctd)
  
  # Build per-spot coord table, aligned to official RDS barcodes
  x_vec <- rep(NA_real_, n_spots)
  y_vec <- rep(NA_real_, n_spots)
  names(x_vec) <- spot_bc
  names(y_vec) <- spot_bc
  
  if (!is.null(rctd_coords)) {
    rctd_bc <- rownames(rctd_coords)
    common <- intersect(spot_bc, rctd_bc)
    log_msg("  Barcodes in common (official & RCTD): ", length(common),
            " (", round(100*length(common)/n_spots, 1), "% of official)")
    if (length(common) > 0) {
      x_vec[common] <- rctd_coords[common, "x"]
      y_vec[common] <- rctd_coords[common, "y"]
    }
  }
  
  # ---- 6. Write spot_meta CSV ----
  spot_meta <- data.frame(
    barcode     = spot_bc,
    x           = x_vec,
    y           = y_vec,
    human_reads = as.integer(h_per_spot),
    mouse_reads = as.integer(m_per_spot),
    human_frac  = round(h_frac, 6),
    stringsAsFactors = FALSE
  )
  rownames(spot_meta) <- spot_bc
  
  meta_csv <- file.path(OUT_DIR, paste0(sample_name, "_spot_meta.csv"))
  write.csv(spot_meta, meta_csv, row.names = FALSE)
  log_msg("Wrote spot meta: ", meta_csv, " (", nrow(spot_meta), " rows)")
  log_msg("  spots with valid coords: ", sum(!is.na(x_vec)))
  
  # ---- 7. Build meta.data for Seurat (indexed by barcode, length = n_spots) ----
  md <- data.frame(
    sample      = rep(sample_name, n_spots),
    species     = character(n_spots),
    human_reads = spot_meta$human_reads,
    mouse_reads = spot_meta$mouse_reads,
    human_frac  = spot_meta$human_frac,
    x           = spot_meta$x,
    y           = spot_meta$y,
    stringsAsFactors = FALSE
  )
  rownames(md) <- spot_bc
  stopifnot(nrow(md) == n_spots)
  stopifnot(identical(rownames(md), spot_bc))
  
  # ---- 8. Build HUMAN subset RDS ----
  log_msg("Building human subset...")
  counts_human <- counts[human_idx, , drop = FALSE]
  rownames(counts_human) <- sub("^GRCh38-", "", rownames(counts_human))
  
  md_human <- md; md_human$species <- "human"
  
  human_obj <- CreateSeuratObject(
    counts    = counts_human,
    assay     = "RNA",
    project   = paste0(sample_name, "_human"),
    meta.data = md_human,
    min.cells = 0,
    min.features = 0
  )
  human_obj@misc$spatial <- list(
    coords = spot_meta[, c("barcode", "x", "y")],
    size_x = size_x,
    size_y = size_y,
    source = paths$official
  )
  
  human_out <- file.path(OUT_DIR, paste0(sample_name, "_human.rds"))
  saveRDS(human_obj, human_out, compress = "gzip")
  log_msg("Wrote: ", human_out,
          " (", nrow(counts_human), " genes x ", ncol(counts_human), " spots)")
  rm(human_obj, counts_human, md_human); gc(verbose = FALSE)
  
  # ---- 9. Build MOUSE subset RDS ----
  log_msg("Building mouse subset...")
  counts_mouse <- counts[mouse_idx, , drop = FALSE]
  rownames(counts_mouse) <- sub("^mm10-", "", rownames(counts_mouse))
  
  md_mouse <- md; md_mouse$species <- "mouse"
  
  mouse_obj <- CreateSeuratObject(
    counts    = counts_mouse,
    assay     = "RNA",
    project   = paste0(sample_name, "_mouse"),
    meta.data = md_mouse,
    min.cells = 0,
    min.features = 0
  )
  mouse_obj@misc$spatial <- list(
    coords = spot_meta[, c("barcode", "x", "y")],
    size_x = size_x,
    size_y = size_y,
    source = paths$official
  )
  
  mouse_out <- file.path(OUT_DIR, paste0(sample_name, "_mouse.rds"))
  saveRDS(mouse_obj, mouse_out, compress = "gzip")
  log_msg("Wrote: ", mouse_out,
          " (", nrow(counts_mouse), " genes x ", ncol(counts_mouse), " spots)")
  rm(mouse_obj, counts_mouse, md_mouse); gc(verbose = FALSE)
  
  # ---- 10. Sanity check output ----
  log_msg("--- Sanity: human_frac distribution (all ", n_spots, " spots retained) ---")
  for (thr in c(0.1, 0.3, 0.5, 0.7, 0.9)) {
    n <- sum(h_frac >= thr, na.rm = TRUE)
    log_msg(sprintf("  human_frac >= %.0f%%: %5d spots (%.1f%%)",
                    100*thr, n, 100*n/n_spots))
  }
  
  rm(obj, counts); gc(verbose = FALSE)
  log_msg("Done: ", sample_name, "\n")
  invisible(NULL)
}

# -------------------------------------------------------------
# Run
# -------------------------------------------------------------
log_msg("Output directory: ", OUT_DIR)
log_msg("Samples: ", paste(names(SAMPLES), collapse = ", "), "\n")

for (nm in names(SAMPLES)) {
  tryCatch(
    split_one_sample(nm, SAMPLES[[nm]]),
    error = function(e) {
      log_msg("!! ERROR on ", nm, ": ", conditionMessage(e))
    }
  )
}

log_msg("========== ALL DONE ==========")
for (f in list.files(OUT_DIR, full.names = FALSE)) {
  info <- file.info(file.path(OUT_DIR, f))
  log_msg(sprintf("  %-40s %8.1f MB", f, info$size / 1024 / 1024))
}
