#!/usr/bin/env Rscript
# ============================================================
# spatial_load_test.R --- SeekSpace Liver_8 单样本加载验证
# ============================================================
# 用法:
#   cd ~/SGS10_analysis
#   Rscript scripts/spatial_load_test.R 2>&1 | tee scripts/logs/spatial_load_test.log
#
# 目的:
#   验证 SeekSpace 数据能正确加载, 在跑全部 10 样本前先单样本试通流程
#
# 步骤:
#   1. 解压 Liver_8 的 zip
#   2. 读 Read10X (matrix + barcodes + features)
#   3. 读 cell_locations.tsv.gz, 作为 spatial reduction
#   4. 创建 Seurat 对象 + Normalize + 简单 PCA
#   5. 加载 H&E image
#   6. 画几个 marker (Acta2, Col1a1, F2r) 在 spatial reduction 上
#   7. 报告 spot 数量, 基因数量, marker 表达情况
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(Matrix)
  library(base64enc)
  library(png)
  library(grid)
})

PROJECT_ROOT <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
setwd(PROJECT_ROOT)

OUT_DIR <- file.path(PROJECT_ROOT, "spatial_test_results")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "plots"), recursive = TRUE, showWarnings = FALSE)

cat("========================================\n")
cat("SeekSpace Liver_8 load test\n")
cat("Start:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================\n\n")

t0 <- proc.time()

# ============================================================
# 1. 路径
# ============================================================
SPATIAL_BASE <- "/mnt/BioData/bio/SGS10例/SGS25090559-8例样本数据/data/matrix"
SAMPLE_DIR   <- file.path(SPATIAL_BASE, "25090559_Liver_8_expression")
ZIP_PATH     <- file.path(SAMPLE_DIR, "Liver_8_filtered_feature_bc_matrix.zip")
HE_IMG       <- file.path(SAMPLE_DIR, "Liver_8_aligned_HE_TIMG.png")
DAPI_IMG     <- file.path(SAMPLE_DIR, "Liver_8_aligned_DAPI.png")

cat("Sample dir:", SAMPLE_DIR, "\n")
cat("zip:       ", ZIP_PATH, "\n")
cat("HE image:  ", HE_IMG, "\n")
cat("DAPI image:", DAPI_IMG, "\n\n")

# ============================================================
# 2. 解压 zip 到临时目录
# ============================================================
TMP_DIR <- file.path(OUT_DIR, "tmp_extract")
dir.create(TMP_DIR, recursive = TRUE, showWarnings = FALSE)

cat("========== Unzipping ==========\n")
unzip(ZIP_PATH, exdir = TMP_DIR)
extracted <- list.files(TMP_DIR, recursive = TRUE, full.names = TRUE)
cat("Extracted files:\n")
for (f in extracted) cat(sprintf("  %s (%s)\n", basename(f), file.info(f)$size))

# 找解压后的 matrix 目录
matrix_dirs <- list.dirs(TMP_DIR, recursive = FALSE)
if (length(matrix_dirs) == 0) {
  # 文件可能直接在 TMP_DIR
  matrix_dir <- TMP_DIR
} else {
  matrix_dir <- matrix_dirs[1]
}
cat(sprintf("\nUsing matrix_dir: %s\n", matrix_dir))

# ============================================================
# 3. 读 expression matrix
# ============================================================
cat("\n========== Reading 10X matrix ==========\n")

# Read10X 需要 barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz
required_files <- c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")
for (f in required_files) {
  fp <- file.path(matrix_dir, f)
  if (!file.exists(fp)) {
    cat(sprintf("  ⚠ missing: %s\n", fp))
  } else {
    cat(sprintf("  ✓ %s (%.1f MB)\n", f, file.info(fp)$size / 1024^2))
  }
}

counts <- Read10X(matrix_dir)
cat(sprintf("\nMatrix loaded: %d genes × %d spots\n",
            nrow(counts), ncol(counts)))
cat(sprintf("  Total UMI: %.2e\n", sum(counts)))
cat(sprintf("  Mean UMI per spot: %.0f\n", mean(colSums(counts))))
cat(sprintf("  Median UMI per spot: %.0f\n", median(colSums(counts))))

# ============================================================
# 4. 读 cell_locations.tsv.gz
# ============================================================
cat("\n========== Reading cell_locations ==========\n")

loc_file <- file.path(matrix_dir, "cell_locations.tsv.gz")
if (!file.exists(loc_file)) {
  stop("cell_locations.tsv.gz not found in matrix dir")
}

cell_locs <- read.table(loc_file, sep = "\t", header = TRUE, row.names = 1)
cat(sprintf("Cell locations: %d rows × %d cols\n",
            nrow(cell_locs), ncol(cell_locs)))
cat("First 5 rows:\n")
print(head(cell_locs, 5))
cat("Column names:", paste(colnames(cell_locs), collapse = ", "), "\n")

# 取前两列作为 x, y (notebook 里是这么做的)
spatial_df <- cell_locs[, 1:2, drop = FALSE]
colnames(spatial_df) <- c("spatial_1", "spatial_2")
cat(sprintf("\nspatial coordinates range:\n"))
cat(sprintf("  spatial_1: [%.0f, %.0f]\n",
            min(spatial_df$spatial_1), max(spatial_df$spatial_1)))
cat(sprintf("  spatial_2: [%.0f, %.0f]\n",
            min(spatial_df$spatial_2), max(spatial_df$spatial_2)))

# barcodes 是否匹配?
common <- intersect(colnames(counts), rownames(spatial_df))
cat(sprintf("\nBarcodes overlap:\n"))
cat(sprintf("  matrix barcodes: %d\n", ncol(counts)))
cat(sprintf("  location barcodes: %d\n", nrow(spatial_df)))
cat(sprintf("  intersection: %d\n", length(common)))

if (length(common) == 0) {
  stop("No barcode overlap between matrix and cell_locations!")
}

# ============================================================
# 5. 创建 Seurat 对象
# ============================================================
cat("\n========== Creating Seurat object ==========\n")

obj <- CreateSeuratObject(counts = counts, project = "Liver_8")
cat(sprintf("Seurat object: %d genes × %d spots\n", nrow(obj), ncol(obj)))

# Add spatial as a reduction
spatial_matrix <- as.matrix(spatial_df[match(colnames(obj), rownames(spatial_df)), ])
obj@reductions$spatial <- CreateDimReducObject(
  embeddings = spatial_matrix,
  key = "spatial_",
  assay = "RNA"
)
cat("✓ spatial reduction added\n")

# ============================================================
# 6. QC + Normalize + PCA + UMAP
# ============================================================
cat("\n========== QC + Normalize ==========\n")

# Mt content
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
cat(sprintf("nFeature_RNA: median %.0f, IQR [%.0f, %.0f]\n",
            median(obj$nFeature_RNA),
            quantile(obj$nFeature_RNA, 0.25),
            quantile(obj$nFeature_RNA, 0.75)))
cat(sprintf("nCount_RNA: median %.0f\n", median(obj$nCount_RNA)))
cat(sprintf("percent.mt: median %.2f%%\n", median(obj$percent.mt)))

obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000,
                             verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:20, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:20, verbose = FALSE)

cat(sprintf("Clusters: %d\n", length(unique(Idents(obj)))))

# ============================================================
# 7. 加载 H&E image (按 SeekSpace 流程)
# ============================================================
cat("\n========== Loading H&E image ==========\n")

samplename <- "Liver_8"

if (file.exists(HE_IMG)) {
  cat(sprintf("HE image: %.1f MB\n", file.info(HE_IMG)$size / 1024^2))
  img_64 <- base64enc::dataURI(file = HE_IMG)
  img_gg <- png::readPNG(HE_IMG)
  img_grob <- grid::rasterGrob(img_gg, interpolate = FALSE,
                                width = grid::unit(1, "npc"),
                                height = grid::unit(1, "npc"))
  obj@misc$info[[samplename]]$img_he <- img_64
  obj@misc$info[[samplename]]$img_he_gg <- img_grob
  cat("✓ HE image loaded\n")
}

if (file.exists(DAPI_IMG)) {
  cat(sprintf("DAPI image: %.1f MB\n", file.info(DAPI_IMG)$size / 1024^2))
  img_64 <- base64enc::dataURI(file = DAPI_IMG)
  img_gg <- png::readPNG(DAPI_IMG)
  img_grob <- grid::rasterGrob(img_gg, interpolate = FALSE,
                                width = grid::unit(1, "npc"),
                                height = grid::unit(1, "npc"))
  obj@misc$info[[samplename]]$img_dapi <- img_64
  obj@misc$info[[samplename]]$img_dapi_gg <- img_grob
  cat("✓ DAPI image loaded\n")
}

# ============================================================
# 8. 检查 v2.3 关键 marker 的表达
# ============================================================
cat("\n========== v2.3 marker check ==========\n")

markers_v23 <- list(
  HSC_fibrosis = c("Acta2", "Col1a1", "Col1a2", "Tagln", "Lox", "Loxl1"),
  PARs_signature = c("F2", "Plg", "F2r", "Pard3"),
  HSC_myeloid_hub = c("Csf1", "Csf1r", "Gas6", "Mertk", "C3"),
  NAM = c("Cd74", "Cd68", "Csf1r", "Itgb2"),
  Cholangiocyte = c("Krt19", "Spp1", "Sox9"),
  Hep_zonation = c("Alb", "Glul", "Cyp2e1", "Pck1", "Cps1"),
  Tumor_proliferation = c("Mki67", "Top2a", "Pcna")
)

all_markers <- unique(unlist(markers_v23))
in_data <- intersect(all_markers, rownames(obj))
cat(sprintf("Markers in data: %d / %d\n", length(in_data), length(all_markers)))
cat(sprintf("Missing: %s\n",
            paste(setdiff(all_markers, in_data), collapse = ", ")))

# Marker pct + mean
data_layer <- GetAssayData(obj, assay = "RNA", layer = "data")
marker_df <- data.frame(
  gene = in_data,
  pct = round(rowMeans(data_layer[in_data, , drop = FALSE] > 0) * 100, 1),
  mean = round(rowMeans(data_layer[in_data, , drop = FALSE]), 3),
  stringsAsFactors = FALSE
)
marker_df <- marker_df[order(-marker_df$mean), ]
cat("\nMarker expression in Liver_8 (sorted by mean):\n")
print(marker_df)

# ============================================================
# 9. 画 spatial plots
# ============================================================
cat("\n========== Drawing spatial plots ==========\n")

# UMAP cluster
p_umap <- DimPlot(obj, reduction = "umap") +
  ggtitle("Liver_8 UMAP (clusters)")
ggsave(file.path(OUT_DIR, "plots", "01_umap.pdf"), p_umap,
       width = 8, height = 7, dpi = 200)

# Spatial cluster
p_spatial <- DimPlot(obj, reduction = "spatial", pt.size = 0.5) +
  ggtitle("Liver_8 spatial (clusters)") +
  coord_fixed()
ggsave(file.path(OUT_DIR, "plots", "02_spatial_clusters.pdf"), p_spatial,
       width = 10, height = 8, dpi = 200)

# Spatial featureplot for top markers
top_markers_for_plot <- c(
  "Acta2", "Col1a1",      # HSC fibrosis (v2.3 axis 1)
  "F2r", "Pard3",         # PARs signature (v2.3 signature finding)
  "Csf1", "Csf1r",        # myeloid hub
  "Cd74", "Cd68",         # NAM marker
  "Krt19", "Spp1",        # Cholangiocyte
  "Alb", "Mki67"          # baseline + proliferation
)
top_markers_for_plot <- intersect(top_markers_for_plot, rownames(obj))

for (g in top_markers_for_plot) {
  pct <- round(mean(data_layer[g, ] > 0) * 100, 1)
  if (pct < 0.5) {
    cat(sprintf("  skip %s (pct %.1f%% too low)\n", g, pct))
    next
  }
  tryCatch({
    p <- FeaturePlot(obj, features = g, reduction = "spatial",
                     pt.size = 0.5, order = TRUE) +
      ggtitle(sprintf("Liver_8 %s (pct %.1f%%)", g, pct)) +
      coord_fixed()
    ggsave(file.path(OUT_DIR, "plots", sprintf("03_spatial_%s.pdf", g)),
           p, width = 9, height = 8, dpi = 200)
    cat(sprintf("  ✓ %s\n", g))
  }, error = function(e) {
    cat(sprintf("  ⚠ %s: %s\n", g, e$message))
  })
}

# 也画一个 grid (multi-feature)
tryCatch({
  p_grid <- FeaturePlot(obj, features = top_markers_for_plot[1:6],
                        reduction = "spatial", pt.size = 0.3,
                        ncol = 3, order = TRUE)
  ggsave(file.path(OUT_DIR, "plots", "04_spatial_marker_grid.pdf"),
         p_grid, width = 18, height = 12, dpi = 200)
  cat("✓ marker grid\n")
}, error = function(e) {
  cat("⚠ grid:", e$message, "\n")
})

# ============================================================
# 10. 保存 RDS
# ============================================================
cat("\n========== Saving RDS ==========\n")

saveRDS(obj, file.path(OUT_DIR, "Liver_8_seekspace.rds"))
cat(sprintf("✓ saved to %s\n", file.path(OUT_DIR, "Liver_8_seekspace.rds")))

# ============================================================
# 11. 清理 + 完成
# ============================================================
cat("\n========== Cleanup ==========\n")
unlink(TMP_DIR, recursive = TRUE)
cat("✓ tmp_extract removed\n")

elapsed <- (proc.time() - t0)["elapsed"] / 60
cat(sprintf("\n========================================\n"))
cat(sprintf("Spatial test complete\n"))
cat(sprintf("  Sample: Liver_8\n"))
cat(sprintf("  Spots: %d\n", ncol(obj)))
cat(sprintf("  Genes: %d\n", nrow(obj)))
cat(sprintf("  Markers in data: %d / %d\n", length(in_data), length(all_markers)))
cat(sprintf("  Output: %s\n", OUT_DIR))
cat(sprintf("  Elapsed: %.1f min\n", elapsed))
cat(sprintf("========================================\n"))
