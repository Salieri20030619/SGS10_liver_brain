#!/usr/bin/env Rscript
# ============================================================
# spatial_batch_10samples.R --- SeekSpace 10 样本批量分析
# ============================================================
# 用法:
#   cd ~/SGS10_analysis
#   nohup Rscript scripts/spatial_batch_10samples.R \
#     > scripts/logs/spatial_batch.log 2>&1 &
#
# 设计:
#   每个样本独立处理 (不 merge, 因为 spatial 坐标不能跨样本合并)
#   1. 解压 zip
#   2. Read10X + cell_locations → Seurat + spatial reduction
#   3. QC + Normalize + PCA + Cluster
#   4. 物种分类 (pure_mouse / host_dom / mixed / tumor_dom / pure_human)
#   5. Zone stratified marker 表 (tumor / margin / host)
#   6. 画核心 plots (human_fraction, species_class, key markers,
#      TGFb module, PARs module, cytokine negative)
#   7. 保存 RDS (含 metadata 但不含 image 以节省空间)
#   8. 汇总 summary.csv
#
# 10 样本:
#   liver: Liver_3, Liver_8, gan56, Bla_liver_1, Con_liver_1
#   brain: Brain_3, Brain_8, nao56, Bla_brain_1, Con_brain_1
#
# 预期:
#   每样本 3-6 min (取决于 spot 数)
#   总耗时 30-60 min
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(Matrix)
  library(base64enc)
  library(png)
  library(grid)
  library(patchwork)
})

PROJECT_ROOT <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
setwd(PROJECT_ROOT)

# ============================================================
# 1. 样本配置
# ============================================================
SPATIAL_BASE_8 <- "/mnt/BioData/bio/SGS10例/SGS25090559-8例样本数据/data/matrix"
SPATIAL_BASE_GAN <- "/mnt/BioData/bio/SGS10例/SGS24070264-gan56/data/matrix"
SPATIAL_BASE_NAO <- "/mnt/BioData/bio/SGS10例/SGS24070263-nao56/data/matrix"

SAMPLES <- list(
  Liver_3     = list(
    folder = file.path(SPATIAL_BASE_8, "25090559_Liver_3_expression"),
    prefix = "Liver_3",
    organ  = "liver", group = "Tumor",  sc_sample = "liver3"
  ),
  Liver_8     = list(
    folder = file.path(SPATIAL_BASE_8, "25090559_Liver_8_expression"),
    prefix = "Liver_8",
    organ  = "liver", group = "Tumor",  sc_sample = "liver8"
  ),
  gan56       = list(
    folder = file.path(SPATIAL_BASE_GAN, "24070264_gan_56_expression"),
    prefix = NULL,  # 特殊处理, 文件名不同
    organ  = "liver", group = "Tumor",  sc_sample = "gan56"
  ),
  Bla_liver_1 = list(
    folder = file.path(SPATIAL_BASE_8, "25090559_Bla_liver_1_expression"),
    prefix = "Bla_liver_1",
    organ  = "liver", group = "Blank",  sc_sample = "bla_liver"
  ),
  Con_liver_1 = list(
    folder = file.path(SPATIAL_BASE_8, "25090559_Con_liver_1_expression"),
    prefix = "Con_liver_1",
    organ  = "liver", group = "Control",sc_sample = "con_liver"
  ),
  Brain_3     = list(
    folder = file.path(SPATIAL_BASE_8, "25090559_Brain_3_expression"),
    prefix = "Brain_3",
    organ  = "brain", group = "Tumor",  sc_sample = "brain3"
  ),
  Brain_8     = list(
    folder = file.path(SPATIAL_BASE_8, "25090559_Brain_8_expression"),
    prefix = "Brain_8",
    organ  = "brain", group = "Tumor",  sc_sample = "brain8"
  ),
  nao56       = list(
    folder = file.path(SPATIAL_BASE_NAO, "24070263_nao_56_expression"),
    prefix = NULL,
    organ  = "brain", group = "Tumor",  sc_sample = "nao56"
  ),
  Bla_brain_1 = list(
    folder = file.path(SPATIAL_BASE_8, "25090559_Bla_brain_1_expression"),
    prefix = "Bla_brain_1",
    organ  = "brain", group = "Blank",  sc_sample = "bla_brain"
  ),
  Con_brain_1 = list(
    folder = file.path(SPATIAL_BASE_8, "25090559_Con_brain_1_expression"),
    prefix = "Con_brain_1",
    organ  = "brain", group = "Control",sc_sample = "con_brain"
  )
)

OUT_DIR <- file.path(PROJECT_ROOT, "spatial_results")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "plots"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "summary"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "zone_tables"), showWarnings = FALSE, recursive = TRUE)

cat("========================================\n")
cat("SeekSpace 10 samples batch processing\n")
cat("Start:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Output:", OUT_DIR, "\n")
cat("========================================\n\n")

t0_global <- proc.time()

# ============================================================
# 2. 核心 marker 定义 (v2.4)
# ============================================================
MARKERS_HUMAN <- c(
  # HCC identity
  "GRCh38-AFP", "GRCh38-DLK1", "GRCh38-IGF2",
  "GRCh38-ALB", "GRCh38-APOB", "GRCh38-TF",
  # Tumor traits
  "GRCh38-MKI67", "GRCh38-EPCAM",
  # Signaling (PDX-derived candidates)
  "GRCh38-TGFB1", "GRCh38-TGFB2", "GRCh38-TGFB3",
  "GRCh38-VEGFA", "GRCh38-CSF1",
  # Coagulation
  "GRCh38-F2", "GRCh38-PLG",
  # Matrix
  "GRCh38-FN1", "GRCh38-COL1A1",
  # Cytokines (negative finding)
  "GRCh38-IL6", "GRCh38-IFNG", "GRCh38-TNF", "GRCh38-IL1B"
)

MARKERS_MOUSE <- c(
  # HSC / fibrosis
  "mm10---Acta2", "mm10---Col1a1", "mm10---Col1a2", "mm10---Tagln",
  "mm10---Lox", "mm10---Loxl1",
  # TGFb pathway
  "mm10---Tgfb1", "mm10---Smad3", "mm10---Smad2", "mm10---Smad7",
  "mm10---Pmepa1", "mm10---Serpine1",
  "mm10---Tgfbr1", "mm10---Tgfbr2", "mm10---Tgfbr3",
  # PARs
  "mm10---F2", "mm10---Plg", "mm10---F2r", "mm10---Pard3",
  # Myeloid hub
  "mm10---Csf1", "mm10---Csf1r", "mm10---Gas6", "mm10---Mertk", "mm10---C3",
  # NAM marker
  "mm10---Cd74", "mm10---Cd68", "mm10---Itgb2",
  # Cholangiocyte
  "mm10---Krt19", "mm10---Spp1",
  # Hepatocyte zonation
  "mm10---Alb", "mm10---Glul", "mm10---Cyp2e1", "mm10---Pck1", "mm10---Cps1",
  # VEGF
  "mm10---Vegfa",
  # Cytokines (negative finding)
  "mm10---Il6", "mm10---Ifng", "mm10---Tnf", "mm10---Il1b",
  # Brain markers (for brain samples)
  "mm10---Gfap", "mm10---Aif1", "mm10---Trem2", "mm10---P2ry12",
  "mm10---Apoe", "mm10---H2-D1", "mm10---Elmo1",
  "mm10---Mbp", "mm10---Olig2"
)

ALL_MARKERS <- c(MARKERS_HUMAN, MARKERS_MOUSE)

# 用于 module score
MODULE_GENES <- list(
  PDX_TGFb_lig    = c("GRCh38-TGFB1", "GRCh38-TGFB2", "GRCh38-TGFB3"),
  host_TGFb_rec   = c("mm10---Tgfbr1", "mm10---Tgfbr2", "mm10---Tgfbr3"),
  host_TGFb_resp  = c("mm10---Smad3", "mm10---Smad2", "mm10---Smad7",
                       "mm10---Pmepa1", "mm10---Serpine1"),
  PDX_angio       = c("GRCh38-VEGFA", "GRCh38-VEGFB", "GRCh38-VEGFC"),
  PDX_hcc_id      = c("GRCh38-AFP", "GRCh38-DLK1", "GRCh38-IGF2"),
  host_myofib     = c("mm10---Acta2", "mm10---Col1a1", "mm10---Col1a2",
                       "mm10---Tagln", "mm10---Lox", "mm10---Loxl1"),
  host_fibrosis_matrix = c("mm10---Col1a1", "mm10---Col1a2", "mm10---Col3a1",
                            "mm10---Col4a1", "mm10---Col4a2", "mm10---Fn1"),
  host_cytokines  = c("mm10---Il6", "mm10---Ifng", "mm10---Tnf", "mm10---Il1b"),
  PDX_cytokines   = c("GRCh38-IL6", "GRCh38-IFNG", "GRCh38-TNF", "GRCh38-IL1B")
)

# ============================================================
# 3. 辅助函数
# ============================================================

find_zip <- function(folder) {
  zips <- list.files(folder, pattern = "\\.zip$", full.names = TRUE)
  # 优先选 filtered_feature_bc_matrix
  filt <- zips[grepl("filtered", zips, ignore.case = TRUE)]
  if (length(filt) > 0) return(filt[1])
  if (length(zips) > 0) return(zips[1])
  return(NA)
}

find_image <- function(folder, type = "HE_TIMG") {
  imgs <- list.files(folder, pattern = sprintf("%s\\.png$", type),
                     full.names = TRUE)
  if (length(imgs) > 0) return(imgs[1])
  return(NA)
}

process_sample <- function(sample_name, cfg) {
  cat(sprintf("\n========================================\n"))
  cat(sprintf("Processing: %s (%s, %s)\n",
              sample_name, cfg$organ, cfg$group))
  cat(sprintf("  sc_sample: %s\n", cfg$sc_sample))
  cat(sprintf("  folder: %s\n", cfg$folder))
  cat(sprintf("========================================\n"))
  t0 <- proc.time()

  # 1. 找 zip
  zip_path <- find_zip(cfg$folder)
  if (is.na(zip_path)) {
    cat("  ✗ no zip found, skip\n")
    return(NULL)
  }
  cat(sprintf("  zip: %s (%.1f MB)\n", basename(zip_path),
              file.info(zip_path)$size / 1024^2))

  # 2. 解压到临时目录
  tmp_dir <- file.path(OUT_DIR, sprintf("tmp_%s", sample_name))
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  unzip(zip_path, exdir = tmp_dir)

  # 找解压后的 matrix 目录
  extracted <- list.dirs(tmp_dir, recursive = FALSE)
  matrix_dir <- if (length(extracted) > 0) extracted[1] else tmp_dir

  # 确认必要文件存在
  req_files <- c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz",
                 "cell_locations.tsv.gz")
  missing <- req_files[!file.exists(file.path(matrix_dir, req_files))]
  if (length(missing) > 0) {
    cat(sprintf("  ✗ missing files: %s\n", paste(missing, collapse = ", ")))
    unlink(tmp_dir, recursive = TRUE)
    return(NULL)
  }

  # 3. Read 10X
  cat("  Reading matrix...\n")
  counts <- tryCatch(
    Read10X(matrix_dir),
    error = function(e) {
      cat(sprintf("  ✗ Read10X failed: %s\n", e$message))
      NULL
    }
  )
  if (is.null(counts)) {
    unlink(tmp_dir, recursive = TRUE)
    return(NULL)
  }
  cat(sprintf("  matrix: %d genes × %d spots\n", nrow(counts), ncol(counts)))

  # 4. Read cell_locations
  cell_locs <- read.table(file.path(matrix_dir, "cell_locations.tsv.gz"),
                          sep = "\t", header = TRUE, row.names = 1)
  spatial_df <- cell_locs[, 1:2, drop = FALSE]
  colnames(spatial_df) <- c("spatial_1", "spatial_2")

  # 5. Create Seurat object
  obj <- CreateSeuratObject(counts = counts, project = sample_name)
  obj$sample <- sample_name
  obj$organ <- cfg$organ
  obj$group <- cfg$group
  obj$sc_sample <- cfg$sc_sample

  # Match spatial coords
  spatial_matrix <- as.matrix(
    spatial_df[match(colnames(obj), rownames(spatial_df)), ]
  )
  obj@reductions$spatial <- CreateDimReducObject(
    embeddings = spatial_matrix,
    key = "spatial_",
    assay = "RNA"
  )

  # 6. QC
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mm10---mt-")
  cat(sprintf("  median nFeature: %.0f, nCount: %.0f\n",
              median(obj$nFeature_RNA), median(obj$nCount_RNA)))

  # 7. Normalize + PCA + Cluster (light)
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:20, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)
  obj <- RunUMAP(obj, dims = 1:20, verbose = FALSE)

  # 8. 物种分类
  genes <- rownames(obj)
  counts_mat <- GetAssayData(obj, assay = "RNA", layer = "counts")
  h_genes <- grep("^GRCh38-", genes, value = TRUE)
  m_genes <- grep("^mm10---", genes, value = TRUE)

  if (length(h_genes) > 0 && length(m_genes) > 0) {
    h_total <- colSums(counts_mat[h_genes, , drop = FALSE])
    m_total <- colSums(counts_mat[m_genes, , drop = FALSE])
    total   <- h_total + m_total
    h_pct   <- h_total / pmax(total, 1) * 100

    obj$human_pct    <- h_pct
    obj$human_count  <- h_total
    obj$mouse_count  <- m_total
    obj$species_class <- cut(h_pct,
                              breaks = c(-1, 5, 30, 70, 95, 101),
                              labels = c("pure_mouse", "host_dom", "mixed",
                                         "tumor_dom", "pure_human"))

    zone_counts <- table(obj$species_class)
    cat("  species_class:\n")
    for (nm in names(zone_counts)) {
      cat(sprintf("    %-12s: %6d (%.1f%%)\n", nm,
                  zone_counts[nm], 100 * zone_counts[nm] / ncol(obj)))
    }

    # 定义 zone 分组用于后续分析
    obj$zone <- ifelse(h_pct > 70, "tumor",
                ifelse(h_pct >= 30, "margin",
                ifelse(h_pct < 5, "host", "host_dom")))
  } else {
    obj$human_pct <- 0
    obj$species_class <- "pure_mouse"
    obj$zone <- "host"
    cat("  (single species, no classification needed)\n")
  }

  # 9. Module scores
  data_layer <- GetAssayData(obj, assay = "RNA", layer = "data")
  for (mod_name in names(MODULE_GENES)) {
    gs <- intersect(MODULE_GENES[[mod_name]], rownames(obj))
    if (length(gs) >= 1) {
      tryCatch({
        obj <- AddModuleScore(obj, features = list(gs),
                              name = mod_name, nbin = 20)
      }, error = function(e) {
        cat(sprintf("    ⚠ module %s failed: %s\n", mod_name, e$message))
      })
    }
  }

  # 10. Zone-stratified marker 表
  zone_df <- data.frame()
  markers_in <- intersect(ALL_MARKERS, rownames(obj))

  for (zn in c("tumor", "margin", "host")) {
    sp <- which(obj$zone == zn)
    if (length(sp) < 5) next
    for (g in markers_in) {
      vals <- data_layer[g, sp]
      zone_df <- rbind(zone_df, data.frame(
        sample = sample_name,
        zone = zn,
        n_spots = length(sp),
        gene = g,
        pct = round(mean(vals > 0) * 100, 2),
        mean = round(mean(vals), 4),
        stringsAsFactors = FALSE
      ))
    }
  }

  if (nrow(zone_df) > 0) {
    write.csv(zone_df,
              file.path(OUT_DIR, "zone_tables",
                        sprintf("%s_zone_markers.csv", sample_name)),
              row.names = FALSE)
  }

  # 11. 画 plots
  plot_dir <- file.path(OUT_DIR, "plots", sample_name)
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

  # UMAP cluster
  tryCatch({
    p <- DimPlot(obj, reduction = "umap", label = TRUE) +
      ggtitle(sprintf("%s UMAP clusters", sample_name))
    ggsave(file.path(plot_dir, "01_umap.pdf"), p,
           width = 8, height = 7, dpi = 150)
  }, error = function(e) cat("    ⚠ umap:", e$message, "\n"))

  # Spatial cluster
  tryCatch({
    p <- DimPlot(obj, reduction = "spatial", pt.size = 0.2) +
      ggtitle(sprintf("%s spatial clusters", sample_name)) +
      coord_fixed()
    ggsave(file.path(plot_dir, "02_spatial_cluster.pdf"), p,
           width = 12, height = 8, dpi = 150)
  }, error = function(e) cat("    ⚠ spatial_cluster:", e$message, "\n"))

  # 如果有双物种 → 画 species / human_pct
  if ("human_pct" %in% colnames(obj@meta.data) && max(obj$human_pct) > 0) {
    tryCatch({
      p <- FeaturePlot(obj, features = "human_pct", reduction = "spatial",
                       pt.size = 0.2, order = TRUE) +
        scale_color_gradientn(colors = c("steelblue", "white", "red"),
                              limits = c(0, 100)) +
        coord_fixed() +
        ggtitle(sprintf("%s: human %% per spot", sample_name))
      ggsave(file.path(plot_dir, "03_human_fraction.pdf"), p,
             width = 12, height = 8, dpi = 150)
    }, error = function(e) cat("    ⚠ human_frac:", e$message, "\n"))

    tryCatch({
      p <- DimPlot(obj, group.by = "species_class", reduction = "spatial",
                   pt.size = 0.2) +
        ggtitle(sprintf("%s: species classification", sample_name)) +
        coord_fixed()
      ggsave(file.path(plot_dir, "04_species_class.pdf"), p,
             width = 12, height = 8, dpi = 150)
    }, error = function(e) cat("    ⚠ species:", e$message, "\n"))
  }

  # v2.4 signature plots: TGFB1 + Smad3 module
  if (cfg$organ == "liver" &&
      "GRCh38-TGFB1" %in% rownames(obj)) {
    tryCatch({
      key_plots <- list()
      for (g in c("GRCh38-TGFB1", "GRCh38-VEGFA", "GRCh38-AFP",
                  "mm10---Smad3", "mm10---Csf1r", "mm10---Csf1",
                  "mm10---Acta2", "mm10---F2r", "mm10---Alb")) {
        if (g %in% rownames(obj)) {
          key_plots[[g]] <- FeaturePlot(obj, features = g,
                                        reduction = "spatial",
                                        pt.size = 0.15, order = TRUE) +
            coord_fixed() + ggtitle(g)
        }
      }
      if (length(key_plots) >= 6) {
        combined <- wrap_plots(key_plots, ncol = 3)
        ggsave(file.path(plot_dir, "05_key_markers_grid.pdf"),
               combined, width = 18, height = 18, dpi = 150)
      }
    }, error = function(e) cat("    ⚠ key_markers:", e$message, "\n"))

    # Module score 4-panel
    tryCatch({
      mods_to_plot <- list()
      for (m in c("PDX_TGFb_lig1", "host_TGFb_rec1",
                  "host_TGFb_resp1", "host_myofib1")) {
        if (m %in% colnames(obj@meta.data)) {
          mods_to_plot[[m]] <- FeaturePlot(obj, features = m,
                                            reduction = "spatial",
                                            pt.size = 0.15, order = TRUE) +
            coord_fixed() + ggtitle(m)
        }
      }
      if (length(mods_to_plot) >= 4) {
        combined <- wrap_plots(mods_to_plot, ncol = 2)
        ggsave(file.path(plot_dir, "06_module_scores.pdf"),
               combined, width = 16, height = 14, dpi = 150)
      }
    }, error = function(e) cat("    ⚠ modules:", e$message, "\n"))
  }

  # Brain-specific: microglia + BBB
  if (cfg$organ == "brain") {
    tryCatch({
      bp <- list()
      for (g in c("mm10---Aif1", "mm10---Trem2", "mm10---P2ry12",
                  "mm10---Apoe", "mm10---H2-D1", "mm10---Elmo1",
                  "mm10---Gfap", "mm10---Mbp", "mm10---Olig2")) {
        if (g %in% rownames(obj)) {
          bp[[g]] <- FeaturePlot(obj, features = g,
                                 reduction = "spatial",
                                 pt.size = 0.15, order = TRUE) +
            coord_fixed() + ggtitle(g)
        }
      }
      if (length(bp) >= 4) {
        combined <- wrap_plots(bp, ncol = 3)
        ggsave(file.path(plot_dir, "05_brain_markers_grid.pdf"),
               combined, width = 18, height = 15, dpi = 150)
      }
    }, error = function(e) cat("    ⚠ brain_markers:", e$message, "\n"))
  }

  # 12. Spatial correlation (only if 双物种 + 有 PDX_TGFb_lig 模块)
  cor_result <- NULL
  if ("PDX_TGFb_lig1" %in% colnames(obj@meta.data) &&
      "host_TGFb_resp1" %in% colnames(obj@meta.data)) {
    cor_result <- list(
      sample = sample_name,
      n_tumor_spots  = sum(obj$zone == "tumor", na.rm = TRUE),
      n_margin_spots = sum(obj$zone == "margin", na.rm = TRUE),
      n_host_spots   = sum(obj$zone == "host", na.rm = TRUE),
      cor_hpct_PDXtgfb = cor(obj$human_pct, obj$PDX_TGFb_lig1,
                              use = "complete.obs"),
      cor_hpct_hostrec = cor(obj$human_pct, obj$host_TGFb_rec1,
                              use = "complete.obs"),
      cor_hpct_hostresp= cor(obj$human_pct, obj$host_TGFb_resp1,
                              use = "complete.obs"),
      cor_PDXtgfb_hostresp = cor(obj$PDX_TGFb_lig1, obj$host_TGFb_resp1,
                                   use = "complete.obs")
    )
    cat(sprintf("  PDX_TGFb_lig vs host_TGFb_resp: r = %.3f\n",
                cor_result$cor_PDXtgfb_hostresp))
  }

  # 13. 保存 RDS (不保存 image 以节省空间)
  obj@misc <- list()  # 清空 misc 省空间
  rds_path <- file.path(OUT_DIR, sprintf("%s_seekspace.rds", sample_name))
  saveRDS(obj, rds_path)
  cat(sprintf("  ✓ saved %s (%.0f MB)\n",
              basename(rds_path), file.info(rds_path)$size / 1024^2))

  # 14. 清理
  unlink(tmp_dir, recursive = TRUE)

  elapsed <- (proc.time() - t0)["elapsed"] / 60
  cat(sprintf("  ✓ completed in %.1f min\n", elapsed))

  # 返回 summary
  list(
    sample = sample_name,
    organ = cfg$organ,
    group = cfg$group,
    sc_sample = cfg$sc_sample,
    n_spots = ncol(obj),
    n_genes = nrow(obj),
    n_clusters = length(unique(Idents(obj))),
    n_human_genes = length(h_genes),
    n_mouse_genes = length(m_genes),
    median_nFeature = median(obj$nFeature_RNA),
    median_nCount = median(obj$nCount_RNA),
    n_tumor_spots = sum(obj$zone == "tumor", na.rm = TRUE),
    n_margin_spots = sum(obj$zone == "margin", na.rm = TRUE),
    n_host_spots = sum(obj$zone == "host", na.rm = TRUE),
    pct_tumor = round(100 * sum(obj$zone == "tumor", na.rm = TRUE) / ncol(obj), 1),
    pct_margin = round(100 * sum(obj$zone == "margin", na.rm = TRUE) / ncol(obj), 1),
    cor_PDXtgfb_hostresp = if (!is.null(cor_result))
      round(cor_result$cor_PDXtgfb_hostresp, 3) else NA,
    elapsed_min = round(elapsed, 1)
  )
}

# ============================================================
# 4. 主循环
# ============================================================
summary_list <- list()
error_log <- list()

for (sample_name in names(SAMPLES)) {
  result <- tryCatch(
    process_sample(sample_name, SAMPLES[[sample_name]]),
    error = function(e) {
      cat(sprintf("\n✗ %s FAILED: %s\n", sample_name, e$message))
      error_log[[sample_name]] <<- e$message
      NULL
    }
  )
  if (!is.null(result)) {
    summary_list[[sample_name]] <- result
  }

  # 强制 gc 节省内存
  gc(verbose = FALSE)
}

# ============================================================
# 5. 写 summary
# ============================================================
cat("\n========================================\n")
cat("Writing summary\n")
cat("========================================\n")

if (length(summary_list) > 0) {
  summary_df <- do.call(rbind, lapply(summary_list, function(x) {
    data.frame(
      sample = x$sample,
      organ = x$organ,
      group = x$group,
      sc_sample = x$sc_sample,
      n_spots = x$n_spots,
      n_genes = x$n_genes,
      n_human_genes = x$n_human_genes,
      n_mouse_genes = x$n_mouse_genes,
      median_nFeature = x$median_nFeature,
      median_nCount = x$median_nCount,
      n_tumor_spots = x$n_tumor_spots,
      n_margin_spots = x$n_margin_spots,
      n_host_spots = x$n_host_spots,
      pct_tumor = x$pct_tumor,
      pct_margin = x$pct_margin,
      cor_PDXtgfb_hostresp = x$cor_PDXtgfb_hostresp,
      elapsed_min = x$elapsed_min,
      stringsAsFactors = FALSE
    )
  }))
  write.csv(summary_df,
            file.path(OUT_DIR, "summary", "sample_summary.csv"),
            row.names = FALSE)
  cat("\nSample summary:\n")
  print(summary_df)
}

if (length(error_log) > 0) {
  err_df <- data.frame(
    sample = names(error_log),
    error = unlist(error_log),
    stringsAsFactors = FALSE
  )
  write.csv(err_df,
            file.path(OUT_DIR, "summary", "error_log.csv"),
            row.names = FALSE)
  cat("\nErrors:\n")
  print(err_df)
}

# ============================================================
# 6. 合并 zone tables
# ============================================================
zone_files <- list.files(file.path(OUT_DIR, "zone_tables"),
                          pattern = "_zone_markers\\.csv$",
                          full.names = TRUE)
if (length(zone_files) > 0) {
  all_zone <- do.call(rbind, lapply(zone_files, read.csv,
                                     stringsAsFactors = FALSE))
  write.csv(all_zone,
            file.path(OUT_DIR, "summary", "all_samples_zone_markers.csv"),
            row.names = FALSE)
  cat(sprintf("\nAll zone markers merged: %d rows\n", nrow(all_zone)))
}

elapsed_global <- (proc.time() - t0_global)["elapsed"] / 60
cat(sprintf("\n========================================\n"))
cat(sprintf("BATCH COMPLETE\n"))
cat(sprintf("  Processed: %d / %d samples\n",
            length(summary_list), length(SAMPLES)))
cat(sprintf("  Failed:    %d samples\n", length(error_log)))
cat(sprintf("  Total time: %.1f min\n", elapsed_global))
cat(sprintf("  Output: %s\n", OUT_DIR))
cat(sprintf("  Finish: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(sprintf("========================================\n"))
