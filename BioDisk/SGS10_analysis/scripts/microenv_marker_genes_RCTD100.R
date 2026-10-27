#!/usr/bin/env Rscript
source("/mnt/BioDisk/wangmingyuan/SGS10_analysis/scripts/fig_save.R")  # [restyle] 三格式输出器
# ============================================================
# microenv_marker_genes_RCTD100.R
# ============================================================
# Purpose:
#   Find marker genes per microenvironment (paper M01-M07 liver, M01-M06 brain).
#   Per-organ analysis with optional tumor-only subset for paper-key M05.
#
# Method:
#   1. Load 10 sample enriched seurat + RCTD weights
#   2. KNN smooth weights, cosine similarity assign domain → microenv
#   3. Merge samples, SCT normalize, PrepSCTFindMarkers
#   4. FindAllMarkers per microenv (one-vs-rest, only.pos)
#   5. Tumor-only subset analysis for M05 (liver) and all microenv (brain)
#
# Inputs:
#   - spatial_results/rctd/enriched/<s>_seekspace_rctd.rds
#   - spatial_results/rctd_full/per_sample/<s>_rctd_results.rds
#   - taskB/spatial_domain_RCTD100/<organ>_domain_ct_profile.csv
#   - taskB/spatial_domain_RCTD100/<organ>_microenv_naming.csv
#
# Outputs:
#   D1_v2_RCTD_composition/microenv_markers_RCTD100/
#     <organ>_combined_seurat.rds
#     <organ>_microenv_markers_all.csv
#     <organ>_microenv_markers_top20.csv
#     <organ>_microenv_markers_dotplot.pdf
#     <organ>_microenv_markers_heatmap.pdf
#     liver_M05_tumor_specific_markers.csv
#     brain_microenv_in_tumor_markers.csv
#     summary.md
#
# Time: ~75 min per organ, ~2.5 hr total
# ============================================================

options(future.globals.maxSize = 50 * 1024^3)

suppressPackageStartupMessages({
  library(Seurat)
  library(future)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(RANN)
})

plan("sequential")
cat(sprintf("future plan: sequential, globals.maxSize=%.0fGB\n",
            options("future.globals.maxSize")[[1]] / 1024^3))

# ---------- Config ----------
PR <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
RCTD_ENRICHED_DIR <- file.path(PR, "spatial_results/rctd/enriched")
RCTD_PER_DIR      <- file.path(PR, "spatial_results/rctd_full/per_sample")
SD_DIR            <- file.path(PR, "taskB/spatial_domain_RCTD100")
OUT_DIR           <- file.path(PR, "D1_v2_RCTD_composition/microenv_markers_RCTD100")
RESTYLE_DIR       <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis/figures/restyle_v1"  # [restyle] 新图统一目录
dir.create(RESTYLE_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

KNN_K <- 15
N_TOP_MARKERS <- 20
set.seed(42)

SAMPLES <- list(
  Bla_brain_1 = list(tissue = "brain", group = "baseline"),
  Brain_3     = list(tissue = "brain", group = "tumor"),
  Brain_8     = list(tissue = "brain", group = "tumor"),
  Con_brain_1 = list(tissue = "brain", group = "control"),
  nao56       = list(tissue = "brain", group = "tumor"),
  Bla_liver_1 = list(tissue = "liver", group = "baseline"),
  Con_liver_1 = list(tissue = "liver", group = "control"),
  Liver_3     = list(tissue = "liver", group = "tumor"),
  Liver_8     = list(tissue = "liver", group = "tumor"),
  gan56       = list(tissue = "liver", group = "tumor")
)

cat("=============================================================\n")
cat("Microenvironment marker genes\n")
cat("Time:", format(Sys.time()), "\n")
cat("Output:", OUT_DIR, "\n")
cat("=============================================================\n\n")

# ---------- Helpers ----------
smooth_weights <- function(coords, weights, k = KNN_K) {
  nn <- RANN::nn2(coords[, c("x", "y")], coords[, c("x", "y")], k = k + 1)
  nn_idx <- nn$nn.idx
  out <- matrix(0, nrow = nrow(weights), ncol = ncol(weights),
                dimnames = dimnames(weights))
  for (j in seq_len(ncol(weights))) {
    out[, j] <- rowMeans(matrix(weights[nn_idx, j], nrow = nrow(weights)))
  }
  rs <- rowSums(out); out / pmax(rs, 1e-9)
}

assign_domain <- function(weights_smoothed, profile_matrix) {
  common_ct <- intersect(colnames(weights_smoothed), colnames(profile_matrix))
  W <- weights_smoothed[, common_ct, drop = FALSE]
  P <- profile_matrix[, common_ct, drop = FALSE]
  W_norm <- sqrt(rowSums(W^2)); P_norm <- sqrt(rowSums(P^2))
  cos_sim <- (W %*% t(P)) / (pmax(W_norm, 1e-9) %o% pmax(P_norm, 1e-9))
  rownames(profile_matrix)[max.col(cos_sim, ties.method = "first")]
}

load_and_label_sample <- function(s, cfg, dom_profile, dom_to_me) {
  cat(sprintf("\n--- [%s] %s/%s ---\n", s, cfg$tissue, cfg$group))
  enriched_path <- file.path(RCTD_ENRICHED_DIR, paste0(s, "_seekspace_rctd.rds"))
  rctd_path <- file.path(RCTD_PER_DIR, paste0(s, "_rctd_results.rds"))
  if (!file.exists(enriched_path) || !file.exists(rctd_path)) {
    cat(sprintf("  SKIP: enriched=%s, rctd=%s\n",
                file.exists(enriched_path), file.exists(rctd_path))); return(NULL)
  }
  
  obj <- readRDS(enriched_path)
  rctd <- readRDS(rctd_path)
  coords  <- rctd@spatialRNA@coords
  weights <- as.matrix(rctd@results$weights)
  rs <- rowSums(weights); weights <- weights / pmax(rs, 1e-9)
  rm(rctd); gc(verbose = FALSE)
  
  common <- intersect(rownames(coords), rownames(weights))
  coords  <- coords[common, , drop = FALSE]
  weights <- weights[common, , drop = FALSE]
  
  weights_smoothed <- smooth_weights(coords, weights, k = KNN_K)
  domain <- assign_domain(weights_smoothed, dom_profile)
  microenv <- dom_to_me[domain]
  
  spot_df <- data.frame(barcode = rownames(weights),
                        domain = domain,
                        microenv_id = unname(microenv),
                        stringsAsFactors = FALSE)
  obj_spots <- intersect(colnames(obj), spot_df$barcode)
  obj <- subset(obj, cells = obj_spots)
  spot_df <- spot_df[match(colnames(obj), spot_df$barcode), ]
  obj$domain <- spot_df$domain
  obj$microenv_id <- spot_df$microenv_id
  obj$sample_id <- s
  obj$group <- cfg$group
  
  cat(sprintf("  spots=%d, microenv:\n", ncol(obj)))
  print(table(obj$microenv_id))
  obj
}

# ============================================================
process_organ <- function(organ) {
  cat(sprintf("\n##############################################\n"))
  cat(sprintf("###  PROCESSING %s\n", toupper(organ)))
  cat(sprintf("##############################################\n"))
  
  prof <- read.csv(file.path(SD_DIR, paste0(organ, "_domain_ct_profile.csv")),
                   row.names = 1, check.names = FALSE)
  prof_mat <- as.matrix(prof)
  me_naming <- read.csv(file.path(SD_DIR, paste0(organ, "_microenv_naming.csv")),
                        stringsAsFactors = FALSE)
  dom_to_me <- setNames(me_naming$microenv_id, me_naming$domain)
  me_to_label <- setNames(me_naming$microenv_label, me_naming$microenv_id)
  me_to_label <- me_to_label[!duplicated(names(me_to_label))]
  
  cat("\nMicroenv mapping:\n")
  for (id in sort(unique(names(me_to_label)))) {
    cat(sprintf("  %s -> %s\n", id, me_to_label[id]))
  }
  
  combined_path <- file.path(OUT_DIR, paste0(organ, "_combined_seurat.rds"))
  
  cat(sprintf("\n[Step 1] Build %s combined Seurat\n", organ))
  organ_samples <- names(SAMPLES)[sapply(SAMPLES, function(c) c$tissue == organ)]
  
  obj_list <- list()
  for (s in organ_samples) {
    obj_list[[s]] <- load_and_label_sample(s, SAMPLES[[s]], prof_mat, dom_to_me)
  }
  obj_list <- obj_list[!sapply(obj_list, is.null)]
  
  cat(sprintf("\nMerging %d samples...\n", length(obj_list)))
  combined <- merge(obj_list[[1]], y = obj_list[-1], add.cell.ids = names(obj_list))
  rm(obj_list); gc(verbose = FALSE)
  
  has_me <- !is.na(combined$microenv_id) & combined$microenv_id != ""
  combined <- subset(combined, cells = colnames(combined)[has_me])
  cat(sprintf("Filtered: %d spots × %d genes\n", ncol(combined), nrow(combined)))
  
  cat("\n[Step 2] SCTransform...\n")
  t0 <- Sys.time()
  combined <- SCTransform(combined, verbose = FALSE, variable.features.n = 3000)
  cat(sprintf("  done in %.1f min\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  
  cat("\n[Step 3] PrepSCTFindMarkers...\n")
  t0 <- Sys.time()
  combined <- PrepSCTFindMarkers(combined, assay = "SCT", verbose = FALSE)
  cat(sprintf("  done in %.1f min\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  saveRDS(combined, combined_path)
  cat(sprintf("Saved: %s\n", basename(combined_path)))
  
  cat("\n[Step 4] FindAllMarkers...\n")
  Idents(combined) <- "microenv_id"
  cat(sprintf("  Idents levels: %s\n",
              paste(levels(Idents(combined)), collapse=",")))
  t0 <- Sys.time()
  markers <- FindAllMarkers(combined, only.pos = TRUE,
                            min.pct = 0.25, logfc.threshold = 0.25,
                            verbose = FALSE)
  cat(sprintf("  done in %.1f min, %d markers\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins")), nrow(markers)))
  
  if (nrow(markers) == 0) { cat("⚠ No markers, abort.\n"); return(NULL) }
  
  markers$cluster <- as.character(markers$cluster)
  markers$microenv_label <- me_to_label[markers$cluster]
  markers$microenv_full <- paste0(markers$cluster, " ", markers$microenv_label)
  
  out_all <- file.path(OUT_DIR, paste0(organ, "_microenv_markers_all.csv"))
  write.csv(markers, out_all, row.names = FALSE)
  cat(sprintf("Saved: %s (%d rows)\n", out_all, nrow(markers)))
  
  cat("\nLabel mapping verification:\n")
  pairs <- unique(markers[, c("cluster", "microenv_label")])
  for (i in seq_len(nrow(pairs))) {
    cat(sprintf("  %s ↔ %s\n", pairs$cluster[i], pairs$microenv_label[i]))
  }
  
  top20 <- markers %>%
    group_by(cluster) %>%
    slice_max(avg_log2FC, n = N_TOP_MARKERS, with_ties = FALSE) %>%
    ungroup() %>% arrange(cluster, desc(avg_log2FC))
  out_top20 <- file.path(OUT_DIR, paste0(organ, "_microenv_markers_top20.csv"))
  write.csv(top20, out_top20, row.names = FALSE)
  cat(sprintf("Saved: %s\n", out_top20))
  
  cat("\n[Step 5] Visualization\n")
  top5_genes <- markers %>%
    group_by(cluster) %>% slice_max(avg_log2FC, n = 5, with_ties = FALSE) %>%
    pull(gene) %>% unique()
  
  p_dot <- DotPlot(combined, features = top5_genes,
                   group.by = "microenv_id", scale = TRUE) +
    RotatedAxis() +
    scale_color_gradient2(low = "navy", mid = "white", high = "firebrick3", midpoint = 0) +
    theme(axis.text.x = element_text(size = 10),   # [restyle] 7->10, 标题已移除
          axis.text.y = element_text(size = 13))   # [restyle] 11->13
  .dw <- max(10, length(top5_genes) * 0.18 + 2)
  .dh <- max(5, length(unique(markers$cluster)) * 0.5 + 2)
  save_fig_all(p_dot,                                        # [restyle] PDF+TIFF+PNG@300
               file.path(RESTYLE_DIR, paste0(organ, "_microenv_markers_dotplot")),
               width = .dw, height = .dh)
  cat("Saved: dotplot pdf\n")
  
  top10_genes <- markers %>%
    group_by(cluster) %>% slice_max(avg_log2FC, n = 5, with_ties = FALSE) %>%
    pull(gene) %>% unique()
  
  set.seed(42)
  cells_subset <- combined@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    group_by(microenv_id) %>%
    slice_sample(n = 200) %>%
    pull(cell)
  combined_sub <- subset(combined, cells = cells_subset)
  combined_sub <- ScaleData(combined_sub, features = top10_genes, verbose = FALSE)
  p_hm <- DoHeatmap(combined_sub, features = top10_genes,
                    group.by = "microenv_id", raster = FALSE, size = 5) +   # [restyle] 组标签 3->5
    NoLegend() +
    theme(axis.text.y = element_text(size = 14))   # [restyle] 基因名放大, 标题已移除
  .hmw <- 14; .hmh <- 9
  # [client] 重点基因: liver=fig12 key / brain=微环境M04 marker; 红名; 各微环境top1黑名; 其余隐藏
  .top1   <- markers %>% dplyr::group_by(cluster) %>%
               dplyr::slice_max(avg_log2FC, n = 1, with_ties = FALSE) %>% dplyr::pull(gene)
  .top1c  <- gsub("^mm10---|^GRCh38-", "", .top1)
  if (organ == "brain") {
    .is_m4 <- as.character(markers$cluster) %in% c("M04","4","M4","04")
    .keyc  <- gsub("^mm10---|^GRCh38-", "", markers$gene[.is_m4])
  } else {
    .keyc  <- c("Cd74","H2-Aa","Psap","Inpp5d","Ctsb","Zeb2","Runx1","Laptm5")
  }
  .ylc     <- gsub("^mm10---|^GRCh38-", "", levels(p_hm$data$Feature))
  .is_key  <- .ylc %in% .keyc
  .is_top1 <- .ylc %in% .top1c
  .lab_col <- ifelse(.is_key, "#C0143B", ifelse(.is_top1, "black", "transparent"))
  cat(sprintf("[%s] clusters = %s\n", organ, paste(unique(as.character(markers$cluster)), collapse=",")))
  cat(sprintf("[%s] red(key)=%d  black(top1)=%d  hidden=%d\n",
              organ, sum(.is_key), sum(.is_top1 & !.is_key), sum(!.is_key & !.is_top1)))
  p_hm <- p_hm + theme(axis.text.y = element_text(colour = .lab_col, size = 18))
  save_fig_all(p_hm,                                         # [restyle] PDF+TIFF+PNG@300
               file.path(RESTYLE_DIR, paste0(organ, "_microenv_markers_heatmap")),
               width = .hmw, height = .hmh)
  cat("Saved: heatmap png/pdf\n")
  rm(combined_sub); gc(verbose = FALSE)
  
  cat("\n[Step 6] Tumor-specific subset analysis\n")
  if (organ == "liver") {
    tumor_obj <- subset(combined, group == "tumor")
    cat(sprintf("Tumor liver: %d spots\n", ncol(tumor_obj)))
    
    if ("M05" %in% unique(tumor_obj$microenv_id)) {
      cat("Re-running PrepSCTFindMarkers on tumor subset...\n")
      t0 <- Sys.time()
      tumor_obj <- PrepSCTFindMarkers(tumor_obj, assay = "SCT", verbose = FALSE)
      cat(sprintf("  Prep done in %.1f min\n",
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
      
      Idents(tumor_obj) <- "microenv_id"
      m05_markers <- FindMarkers(tumor_obj, ident.1 = "M05",
                                  only.pos = TRUE, min.pct = 0.25,
                                  logfc.threshold = 0.25, verbose = FALSE)
      m05_markers$gene <- rownames(m05_markers)
      m05_markers$comparison <- "M05 vs other (in tumor only)"
      m05_markers <- m05_markers %>% arrange(desc(avg_log2FC))
      
      out_m05 <- file.path(OUT_DIR, "liver_M05_tumor_specific_markers.csv")
      write.csv(m05_markers, out_m05, row.names = FALSE)
      cat(sprintf("Saved: %s (%d markers)\n", out_m05, nrow(m05_markers)))
      
      cat("\n--- M05 top 10 markers (in tumor) ---\n")
      print(head(m05_markers[, c("gene","avg_log2FC","p_val_adj","pct.1","pct.2")], 10))
    }
    rm(tumor_obj); gc(verbose = FALSE)
  } else if (organ == "brain") {
    tumor_obj <- subset(combined, group == "tumor")
    cat(sprintf("Tumor brain: %d spots\n", ncol(tumor_obj)))
    
    cat("Re-running PrepSCTFindMarkers on brain tumor subset...\n")
    t0 <- Sys.time()
    tumor_obj <- PrepSCTFindMarkers(tumor_obj, assay = "SCT", verbose = FALSE)
    cat(sprintf("  Prep done in %.1f min\n",
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    
    Idents(tumor_obj) <- "microenv_id"
    tumor_markers <- FindAllMarkers(tumor_obj, only.pos = TRUE,
                                     min.pct = 0.25, logfc.threshold = 0.25,
                                     verbose = FALSE)
    if (nrow(tumor_markers) > 0) {
      tumor_markers$cluster <- as.character(tumor_markers$cluster)
      tumor_markers$microenv_label <- me_to_label[tumor_markers$cluster]
      out_bt <- file.path(OUT_DIR, "brain_microenv_in_tumor_markers.csv")
      write.csv(tumor_markers, out_bt, row.names = FALSE)
      cat(sprintf("Saved: %s (%d markers)\n", out_bt, nrow(tumor_markers)))
    }
    rm(tumor_obj); gc(verbose = FALSE)
  }
  
  rm(combined); gc(verbose = FALSE)
  
  list(organ = organ,
       n_microenv = length(unique(markers$cluster)),
       n_markers = nrow(markers),
       top_genes_per_me = top20 %>%
         group_by(cluster, microenv_label) %>%
         slice_head(n = 5) %>%
         summarise(top5 = paste(gene, collapse = ", "), .groups = "drop"))
}

# ============================================================
results <- list()
results[["liver"]] <- process_organ("liver")
results[["brain"]] <- process_organ("brain")

# Summary
md_file <- file.path(OUT_DIR, "summary.md")
con <- file(md_file, "w")
writeLines(c(
  "# Microenvironment Marker Genes Summary",
  "",
  paste0("**Generated:** ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## Method",
  "- 100% RCTD weights → KNN smooth → cosine assign domain → microenv",
  "- SCTransform per organ (10 samples merged)",
  "- FindAllMarkers with only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25",
  "- Tumor-only subset analysis (PrepSCTFindMarkers re-run)",
  "",
  "## Results"
), con)

for (organ in c("liver", "brain")) {
  r <- results[[organ]]
  if (is.null(r)) { writeLines(paste0("\n### ", toupper(organ), ": FAILED"), con); next }
  writeLines(c("",
               paste0("### ", toupper(organ), " (", r$n_microenv,
                      " microenvs, ", r$n_markers, " markers)"),
               "", "**Top 5 marker genes per microenv (by avg_log2FC):**", ""), con)
  for (i in seq_len(nrow(r$top_genes_per_me))) {
    writeLines(paste0("- **", r$top_genes_per_me$cluster[i], " ",
                      r$top_genes_per_me$microenv_label[i], "**: ",
                      r$top_genes_per_me$top5[i]), con)
  }
}

writeLines(c("",
             "## Files",
             "- `<organ>_microenv_markers_all.csv` / `_top20.csv` / `_dotplot.pdf` / `_heatmap.pdf`",
             "- `liver_M05_tumor_specific_markers.csv` (paper-key niche)",
             "- `brain_microenv_in_tumor_markers.csv` (brain microenv in tumor)",
             "- `<organ>_combined_seurat.rds` (intermediate)"), con)
close(con)
cat("\nSaved:", md_file, "\n")

cat("\n=============================================================\n")
cat("Done. Output:", OUT_DIR, "\n")
cat("=============================================================\n")
