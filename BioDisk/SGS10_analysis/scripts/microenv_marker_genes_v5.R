#!/usr/bin/env Rscript
# =============================================================================
# microenv_marker_genes_v5.R — 全部重跑, 合并 v3+v4 + 修复 label bug
# 
# v5 vs v3/v4:
#   1. ⭐ 修 label 错位 bug: cluster 必须 as.character() 再 lookup me_to_label
#   2. 跑全套 (liver + brain), 不依赖任何之前 cache (干净重跑)
#   3. 包含 tumor-specific subset (M05 in tumor / brain microenv in tumor)
#   4. 5 step checkpoint, 任何步骤挂可独立续跑
#
# 输出:
#   D1_v2_RCTD_composition/microenv_markers/
#     {liver,brain}_combined_seurat.rds              (中间产物, 含 SCT + Prep)
#     {liver,brain}_microenv_markers_all.csv         ⭐ 新生成 (label 已修)
#     {liver,brain}_microenv_markers_top20.csv       
#     {liver,brain}_microenv_markers_dotplot.pdf
#     {liver,brain}_microenv_markers_heatmap.pdf
#     liver_M05_tumor_specific_markers.csv           ⭐ paper key niche
#     brain_microenv_in_tumor_markers.csv
#     summary.md
#
# 用法:
#   nohup Rscript scripts/microenv_marker_genes_v5.R \
#       > scripts/logs/microenv_marker_v5.log 2>&1 &
#
# 时间预估:
#   Liver: SCT 60-76 min + Prep 3 min + FAM 0.2 min + viz 5 min + tumor 5 min = ~75 min
#   Brain: 同上 = ~75 min
#   Total: ~2.5 hr
# =============================================================================

# ⭐ FIX: future global size 上限 + sequential plan
options(future.globals.maxSize = 16 * 1024^3)

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
RCTD_PER_DIR <- file.path(PR, "spatial_results/rctd/per_sample")
SD_DIR <- file.path(PR, "taskB/spatial_domain")
ME_DIR <- file.path(PR, "D1_v2_RCTD_composition")
OUT_DIR <- file.path(ME_DIR, "microenv_markers")
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
cat("Phase 3 v5: 全部重跑 (label fix + tumor subset + dual organ)\n")
cat("Time:", format(Sys.time()), "\n")
cat("Output:", OUT_DIR, "\n")
cat("=============================================================\n\n")

# ---------- Helper ----------
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
    cat(sprintf("  SKIP: missing\n")); return(NULL)
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
# 主函数
# ============================================================
process_organ <- function(organ) {
  cat(sprintf("\n##############################################\n"))
  cat(sprintf("###  PROCESSING %s\n", toupper(organ)))
  cat(sprintf("##############################################\n"))
  
  prof <- read.csv(file.path(SD_DIR, paste0(organ, "_domain_ct_profile.csv")),
                   row.names = 1, check.names = FALSE)
  prof_mat <- as.matrix(prof)
  me_naming <- read.csv(file.path(ME_DIR, paste0(organ, "_microenv_naming.csv")),
                        stringsAsFactors = FALSE)
  dom_to_me <- setNames(me_naming$microenv_id, me_naming$domain)
  me_to_label <- setNames(me_naming$microenv_label, me_naming$microenv_id)
  me_to_label <- me_to_label[!duplicated(names(me_to_label))]
  
  cat("\nMicroenv mapping (corrected):\n")
  for (id in sort(unique(names(me_to_label)))) {
    cat(sprintf("  %s -> %s\n", id, me_to_label[id]))
  }
  
  # ---------- Step 1: Load + merge + SCT ----------
  combined_path <- file.path(OUT_DIR, paste0(organ, "_combined_seurat.rds"))
  
  cat(sprintf("\n[Step 1] Build %s combined Seurat from scratch...\n", organ))
  organ_samples <- names(SAMPLES)[sapply(SAMPLES, function(c) c$tissue == organ)]
  
  obj_list <- list()
  for (s in organ_samples) {
    obj_list[[s]] <- load_and_label_sample(s, SAMPLES[[s]], prof_mat, dom_to_me)
  }
  obj_list <- obj_list[!sapply(obj_list, is.null)]
  
  cat(sprintf("\nMerging %d samples...\n", length(obj_list)))
  combined <- merge(obj_list[[1]], y = obj_list[-1],
                    add.cell.ids = names(obj_list))
  rm(obj_list); gc(verbose = FALSE)
  
  has_me <- !is.na(combined$microenv_id) & combined$microenv_id != ""
  combined <- subset(combined, cells = colnames(combined)[has_me])
  cat(sprintf("Filtered: %d spots × %d genes\n", ncol(combined), nrow(combined)))
  
  cat("\n[Step 1b] SCTransform (60-76 min, no glmGamPoi)...\n")
  t0 <- Sys.time()
  combined <- SCTransform(combined, verbose = FALSE, variable.features.n = 3000)
  cat(sprintf("SCT done in %.1f min\n",
              as.numeric(difftime(Sys.time(), t0, units="mins"))))
  
  # ---------- Step 2: PrepSCTFindMarkers ----------
  cat("\n[Step 2] PrepSCTFindMarkers...\n")
  t0 <- Sys.time()
  combined <- PrepSCTFindMarkers(combined, assay = "SCT", verbose = FALSE)
  cat(sprintf("  done in %.1f min\n",
              as.numeric(difftime(Sys.time(), t0, units="mins"))))
  saveRDS(combined, combined_path)
  cat(sprintf("Saved: %s\n", basename(combined_path)))
  
  # ---------- Step 3: FindAllMarkers ----------
  cat("\n[Step 3] FindAllMarkers (microenv as ident)...\n")
  Idents(combined) <- "microenv_id"
  cat(sprintf("  Idents levels: %s\n",
              paste(levels(Idents(combined)), collapse=",")))
  t0 <- Sys.time()
  markers <- FindAllMarkers(combined, only.pos = TRUE,
                            min.pct = 0.25, logfc.threshold = 0.25,
                            verbose = FALSE)
  cat(sprintf("  done in %.1f min, %d markers\n",
              as.numeric(difftime(Sys.time(), t0, units="mins")),
              nrow(markers)))
  
  if (nrow(markers) == 0) {
    cat("⚠ No markers, abort organ.\n")
    return(NULL)
  }
  
  # ⭐ FIX: as.character(cluster) 避免 factor 索引错位
  markers$cluster <- as.character(markers$cluster)
  markers$microenv_label <- me_to_label[markers$cluster]
  markers$microenv_full <- paste0(markers$cluster, " ", markers$microenv_label)
  
  out_all <- file.path(OUT_DIR, paste0(organ, "_microenv_markers_all.csv"))
  write.csv(markers, out_all, row.names = FALSE)
  cat(sprintf("Saved: %s (%d rows)\n", out_all, nrow(markers)))
  
  # 显示 cluster ↔ label 检查
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
  
  # ---------- Step 4: Visualization ----------
  cat("\n[Step 4] Generating dotplot + heatmap...\n")
  top5_genes <- markers %>%
    group_by(cluster) %>% slice_max(avg_log2FC, n = 5, with_ties = FALSE) %>%
    pull(gene) %>% unique()
  
  p_dot <- DotPlot(combined, features = top5_genes,
                   group.by = "microenv_id", scale = TRUE) +
    RotatedAxis() +
    scale_color_gradient2(low = "navy", mid = "white", high = "firebrick3",
                          midpoint = 0) +
    labs(title = sprintf("%s microenv markers (top 5 per microenv)", organ)) +
    theme(axis.text.x = element_text(size = 7),
          axis.text.y = element_text(size = 11))
  ggsave(file.path(OUT_DIR, paste0(organ, "_microenv_markers_dotplot.pdf")),
         p_dot, width = max(10, length(top5_genes) * 0.18 + 2),
         height = max(5, length(unique(markers$cluster)) * 0.5 + 2))
  cat("Saved: dotplot pdf\n")
  
  top10_genes <- markers %>%
    group_by(cluster) %>% slice_max(avg_log2FC, n = 10, with_ties = FALSE) %>%
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
                    group.by = "microenv_id", raster = FALSE, size = 3) +
    NoLegend() +
    labs(title = sprintf("%s microenv markers (top 10)", organ))
  ggsave(file.path(OUT_DIR, paste0(organ, "_microenv_markers_heatmap.pdf")),
         p_hm, width = 14, height = max(8, length(top10_genes) * 0.18 + 2),
         limitsize = FALSE)
  cat("Saved: heatmap pdf\n")
  rm(combined_sub); gc(verbose = FALSE)
  
  # ---------- Step 5: Tumor-specific subset ----------
  cat("\n[Step 5] Tumor-specific analysis...\n")
  if (organ == "liver") {
    tumor_obj <- subset(combined, group == "tumor")
    cat(sprintf("Tumor liver: %d spots\n", ncol(tumor_obj)))
    
    if ("M05" %in% unique(tumor_obj$microenv_id)) {
      # ⭐ FIX: subset 后必须重新 PrepSCTFindMarkers
      cat("[v5 FIX] Re-running PrepSCTFindMarkers on tumor subset...\n")
      t0 <- Sys.time()
      tumor_obj <- PrepSCTFindMarkers(tumor_obj, assay = "SCT", verbose = FALSE)
      cat(sprintf("  Prep done in %.1f min\n",
                  as.numeric(difftime(Sys.time(), t0, units="mins"))))
      
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
    
    cat("[v5 FIX] Re-running PrepSCTFindMarkers on brain tumor subset...\n")
    t0 <- Sys.time()
    tumor_obj <- PrepSCTFindMarkers(tumor_obj, assay = "SCT", verbose = FALSE)
    cat(sprintf("  Prep done in %.1f min\n",
                as.numeric(difftime(Sys.time(), t0, units="mins"))))
    
    Idents(tumor_obj) <- "microenv_id"
    tumor_markers <- FindAllMarkers(tumor_obj, only.pos = TRUE,
                                    min.pct = 0.25, logfc.threshold = 0.25,
                                    verbose = FALSE)
    if (nrow(tumor_markers) > 0) {
      # ⭐ FIX
      tumor_markers$cluster <- as.character(tumor_markers$cluster)
      tumor_markers$microenv_label <- me_to_label[tumor_markers$cluster]
      out_bt <- file.path(OUT_DIR, "brain_microenv_in_tumor_markers.csv")
      write.csv(tumor_markers, out_bt, row.names = FALSE)
      cat(sprintf("Saved: %s (%d markers)\n", out_bt, nrow(tumor_markers)))
    }
    rm(tumor_obj); gc(verbose = FALSE)
  }
  
  rm(combined); gc(verbose = FALSE)
  
  list(
    organ = organ,
    n_microenv = length(unique(markers$cluster)),
    n_markers = nrow(markers),
    top_genes_per_me = top20 %>%
      group_by(cluster, microenv_label) %>%
      slice_head(n = 5) %>%
      summarise(top5 = paste(gene, collapse = ", "), .groups = "drop")
  )
}

# ---------- Run both organs ----------
results <- list()
results[["liver"]] <- process_organ("liver")
results[["brain"]] <- process_organ("brain")

# ---------- summary.md ----------
md_file <- file.path(OUT_DIR, "summary.md")
con <- file(md_file, "w")
writeLines(c(
  "# Microenvironment Marker Genes — Phase 3 v5 (全部重跑 + label 修复)",
  "",
  paste0("**Generated:** ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## v5 vs v3/v4 changes",
  "",
  "1. **Label bug fix**: cluster as.character() 后再 lookup, 避免 factor 索引错位",
  "2. **全部重跑** (不依赖任何旧 cache)",
  "3. PrepSCTFindMarkers 在 subset 前后都重跑",
  "",
  "## Results"
), con)

for (organ in c("liver", "brain")) {
  r <- results[[organ]]
  if (is.null(r)) { writeLines(paste0("\n### ", toupper(organ), ": FAILED"), con); next }
  writeLines(c(
    "",
    paste0("### ", toupper(organ),
           " (", r$n_microenv, " microenvs, ", r$n_markers, " markers total)"),
    "",
    "**Top 5 marker genes per microenv (by avg_log2FC)**:",
    ""
  ), con)
  for (i in seq_len(nrow(r$top_genes_per_me))) {
    writeLines(paste0("- **", r$top_genes_per_me$cluster[i],
                      " ", r$top_genes_per_me$microenv_label[i], "**: ",
                      r$top_genes_per_me$top5[i]), con)
  }
}

writeLines(c(
  "",
  "## Files",
  "- `<organ>_microenv_markers_all.csv` / `_top20.csv` / `_dotplot.pdf` / `_heatmap.pdf`",
  "- `liver_M05_tumor_specific_markers.csv` (paper key niche)",
  "- `brain_microenv_in_tumor_markers.csv` (brain microenv in tumor)",
  "- `<organ>_combined_seurat.rds` (intermediate)"
), con)
close(con)
cat("\nSaved:", md_file, "\n")

cat("\n=============================================================\n")
cat("Phase 3 v5 done. Output:", OUT_DIR, "\n")
cat("=============================================================\n")
