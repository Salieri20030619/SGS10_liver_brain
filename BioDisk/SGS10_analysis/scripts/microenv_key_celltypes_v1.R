#!/usr/bin/env Rscript
# =============================================================================
# microenv_key_celltypes_v1.R
# -----------------------------------------------------------------------------
# 计算每个微环境的关键细胞类型 (celltype mean weight ranking)
# 
# 算法:
#   1. 加载 10 sample 的 RCTD weights_matrix.csv
#   2. 加载 microenv_naming.csv 的 domain → microenv 映射
#   3. 复用 he_overlay_v5 的 cosine similarity 给每 spot 分配 domain
#   4. 把 domain 映射到 microenv
#   5. 每微环境取所有 spot 的 celltype mean weight, 按 weight 排序
#   6. 输出每微环境 top 5 关键细胞 + heatmap + dotplot
#
# 输入:
#   spatial_results/rctd/per_sample/<sample>_weights_matrix.csv (10个)
#   spatial_results/rctd/per_sample/<sample>_rctd_results.rds (10个, 取 coords)
#   taskB/spatial_domain/{liver,brain}_domain_ct_profile.csv  (D01-D10 profile)
#   D1_v2_RCTD_composition/{liver,brain}_microenv_naming.csv (domain→microenv)
#
# 输出:
#   D1_v2_RCTD_composition/microenv_key_celltypes/
#     {liver,brain}_microenv_celltypes_full.csv  (microenv × CT × mean_weight × n_spots)
#     {liver,brain}_microenv_celltypes_top5.csv  (每 microenv top 5)
#     {liver,brain}_microenv_celltype_heatmap.pdf
#     {liver,brain}_microenv_celltype_dotplot.pdf
#     {liver,brain}_microenv_specificity.csv  (M_x 的 CT 相对于 其他 M 的富集 log2FC)
#     summary.md
#
# 用法:
#   Rscript scripts/microenv_key_celltypes_v1.R
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(RANN)
})

# ---------- Config ----------
PR <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
RCTD_DIR <- file.path(PR, "spatial_results/rctd/per_sample")
SD_DIR   <- file.path(PR, "taskB/spatial_domain")
ME_DIR   <- file.path(PR, "D1_v2_RCTD_composition")
OUT_DIR  <- file.path(ME_DIR, "microenv_key_celltypes")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

KNN_K <- 15

set.seed(42)

cat("=============================================================\n")
cat("Phase 2: Microenvironment key celltype analysis\n")
cat("Output:", OUT_DIR, "\n")
cat("=============================================================\n\n")

# ---------- Sample list ----------
SAMPLES <- list(
  Bla_brain_1 = "brain", Brain_3 = "brain", Brain_8 = "brain",
  Con_brain_1 = "brain", nao56 = "brain",
  Bla_liver_1 = "liver", Con_liver_1 = "liver",
  Liver_3 = "liver", Liver_8 = "liver", gan56 = "liver"
)

# ---------- Helper: KNN smooth weights (复用 v4) ----------
smooth_weights <- function(coords, weights, k = KNN_K) {
  nn <- RANN::nn2(coords[, c("x", "y")], coords[, c("x", "y")], k = k + 1)
  nn_idx <- nn$nn.idx
  out <- matrix(0, nrow = nrow(weights), ncol = ncol(weights),
                dimnames = dimnames(weights))
  for (j in seq_len(ncol(weights))) {
    out[, j] <- rowMeans(matrix(weights[nn_idx, j], nrow = nrow(weights)))
  }
  rs <- rowSums(out)
  out <- out / pmax(rs, 1e-9)
  out
}

# ---------- Helper: cosine similarity assign domain (复用 v4) ----------
assign_domain <- function(weights_smoothed, profile_matrix) {
  common_ct <- intersect(colnames(weights_smoothed), colnames(profile_matrix))
  W <- weights_smoothed[, common_ct, drop = FALSE]
  P <- profile_matrix[, common_ct, drop = FALSE]
  W_norm <- sqrt(rowSums(W^2))
  P_norm <- sqrt(rowSums(P^2))
  dot    <- W %*% t(P)
  cos_sim <- dot / (pmax(W_norm, 1e-9) %o% pmax(P_norm, 1e-9))
  dom_idx <- max.col(cos_sim, ties.method = "first")
  rownames(profile_matrix)[dom_idx]
}

# ---------- Helper: 加载 + 标 microenv per sample ----------
load_sample_microenv <- function(s, organ, dom_profile, dom_to_me) {
  cat(sprintf("  [%s] %s\n", organ, s))
  
  rctd_path <- file.path(RCTD_DIR, paste0(s, "_rctd_results.rds"))
  if (!file.exists(rctd_path)) {
    cat(sprintf("    SKIP: %s not found\n", rctd_path))
    return(NULL)
  }
  
  rctd <- readRDS(rctd_path)
  coords  <- rctd@spatialRNA@coords
  weights <- as.matrix(rctd@results$weights)
  rs <- rowSums(weights); weights <- weights / pmax(rs, 1e-9)
  rm(rctd); gc(verbose = FALSE)
  
  common <- intersect(rownames(coords), rownames(weights))
  coords  <- coords[common, , drop = FALSE]
  weights <- weights[common, , drop = FALSE]
  
  # KNN smooth + cosine assign
  weights_smoothed <- smooth_weights(coords, weights, k = KNN_K)
  domain <- assign_domain(weights_smoothed, dom_profile)
  microenv <- dom_to_me[domain]
  
  # 输出 spot-level 表 (raw weights, NOT smoothed)
  data.frame(
    sample = s,
    barcode = rownames(weights),
    domain = domain,
    microenv_id = unname(microenv),
    weights,         # 22 (brain) / 36 (liver) celltype 列
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

# ============================================================
# Per-organ processing
# ============================================================

process_organ <- function(organ) {
  cat(sprintf("\n========== %s ==========\n", toupper(organ)))
  
  # Load domain profile + microenv naming
  prof <- read.csv(file.path(SD_DIR, paste0(organ, "_domain_ct_profile.csv")),
                   row.names = 1, check.names = FALSE)
  prof_mat <- as.matrix(prof)
  
  me_naming <- read.csv(file.path(ME_DIR, paste0(organ, "_microenv_naming.csv")),
                        stringsAsFactors = FALSE)
  dom_to_me <- setNames(me_naming$microenv_id, me_naming$domain)
  me_to_label <- setNames(me_naming$microenv_label, me_naming$microenv_id)
  me_to_label <- me_to_label[!duplicated(names(me_to_label))]
  
  cat(sprintf("Domain profile: %d × %d\n", nrow(prof_mat), ncol(prof_mat)))
  cat(sprintf("Microenv mapping: %d microenvs\n", length(unique(me_naming$microenv_id))))
  
  # Load all samples for this organ
  organ_samples <- names(SAMPLES)[unlist(SAMPLES) == organ]
  cat(sprintf("Loading %d samples...\n", length(organ_samples)))
  
  sample_dfs <- lapply(organ_samples, function(s) {
    load_sample_microenv(s, organ, prof_mat, dom_to_me)
  })
  sample_dfs <- sample_dfs[!sapply(sample_dfs, is.null)]
  
  # Combine
  combined <- bind_rows(sample_dfs)
  cat(sprintf("Combined: %d spots from %d samples\n",
              nrow(combined), length(unique(combined$sample))))
  
  # Add tumor/baseline group
  tumor_samples <- c("Liver_3", "Liver_8", "gan56", "Brain_3", "Brain_8", "nao56")
  combined$group <- ifelse(combined$sample %in% tumor_samples, "tumor", "baseline")
  
  # ---------- 1. microenv × celltype mean weight ----------
  ct_cols <- setdiff(colnames(combined),
                     c("sample","barcode","domain","microenv_id","group"))
  
  cat(sprintf("\nComputing microenv × celltype mean weights (%d celltypes)\n",
              length(ct_cols)))
  
  me_ct <- combined %>%
    select(microenv_id, all_of(ct_cols)) %>%
    group_by(microenv_id) %>%
    summarise(across(everything(), mean), .groups = "drop") %>%
    pivot_longer(cols = -microenv_id, names_to = "celltype",
                 values_to = "mean_weight")
  
  me_size <- combined %>%
    group_by(microenv_id) %>%
    summarise(n_spots = n(), .groups = "drop")
  
  me_ct <- me_ct %>%
    left_join(me_size, by = "microenv_id") %>%
    mutate(microenv_label = me_to_label[microenv_id]) %>%
    arrange(microenv_id, desc(mean_weight))
  
  out_full <- file.path(OUT_DIR, paste0(organ, "_microenv_celltypes_full.csv"))
  write.csv(me_ct, out_full, row.names = FALSE)
  cat(sprintf("Saved: %s\n", out_full))
  
  # ---------- 2. Top 5 per microenv ----------
  top5 <- me_ct %>%
    group_by(microenv_id) %>%
    slice_max(mean_weight, n = 5, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(rank = ave(mean_weight, microenv_id,
                      FUN = function(x) rank(-x, ties.method = "first")))
  
  out_top5 <- file.path(OUT_DIR, paste0(organ, "_microenv_celltypes_top5.csv"))
  write.csv(top5, out_top5, row.names = FALSE)
  cat(sprintf("Saved: %s\n", out_top5))
  
  cat("\n--- Top 3 celltypes per microenv ---\n")
  top3_print <- top5 %>%
    filter(rank <= 3) %>%
    group_by(microenv_id, microenv_label, n_spots) %>%
    summarise(top_celltypes = paste(sprintf("%s(%.2f)", celltype, mean_weight),
                                    collapse = " | "),
              .groups = "drop")
  print(as.data.frame(top3_print))
  
  # ---------- 3. Specificity (M_x 中 CT 相对于其他 M 的富集) ----------
  # 用每个 microenv 的 CT mean weight 除以全 sample mean
  cat("\nComputing microenv specificity (vs global mean)...\n")
  global_mean <- combined %>%
    summarise(across(all_of(ct_cols), mean)) %>%
    pivot_longer(cols = everything(), names_to = "celltype",
                 values_to = "global_weight")
  
  spec <- me_ct %>%
    select(microenv_id, microenv_label, celltype, mean_weight, n_spots) %>%
    left_join(global_mean, by = "celltype") %>%
    mutate(
      log2FC = log2((mean_weight + 1e-4) / (global_weight + 1e-4)),
      enriched = log2FC > 1   # 富集 2x 以上
    ) %>%
    arrange(microenv_id, desc(log2FC))
  
  out_spec <- file.path(OUT_DIR, paste0(organ, "_microenv_specificity.csv"))
  write.csv(spec, out_spec, row.names = FALSE)
  cat(sprintf("Saved: %s\n", out_spec))
  
  # 每微环境最 specific 的 top 3
  spec_top <- spec %>%
    filter(enriched) %>%
    group_by(microenv_id) %>%
    slice_max(log2FC, n = 3, with_ties = FALSE) %>%
    ungroup()
  
  cat("\n--- Top 3 enriched (vs global) celltypes per microenv ---\n")
  spec_top_print <- spec_top %>%
    group_by(microenv_id, microenv_label) %>%
    summarise(top_enriched = paste(sprintf("%s(log2FC=%.1f)",
                                            celltype, log2FC),
                                    collapse = " | "),
              .groups = "drop")
  print(as.data.frame(spec_top_print))
  
  # ---------- 4. Heatmap ----------
  cat("\nGenerating heatmap...\n")
  # 取每个 microenv top 3 enriched + present in all microenv
  top_cts_per_me <- top5 %>%
    filter(rank <= 3) %>%
    pull(celltype) %>%
    unique()
  
  hm_mat <- me_ct %>%
    filter(celltype %in% top_cts_per_me) %>%
    select(microenv_id, celltype, mean_weight) %>%
    pivot_wider(names_from = celltype, values_from = mean_weight, values_fill = 0) %>%
    as.data.frame()
  
  hm_labels <- me_to_label[hm_mat$microenv_id]
  rownames(hm_mat) <- paste0(hm_mat$microenv_id, " ", hm_labels)
  hm_mat <- as.matrix(hm_mat[, -1])
  
  pdf(file.path(OUT_DIR, paste0(organ, "_microenv_celltype_heatmap.pdf")),
      width = max(8, ncol(hm_mat) * 0.4 + 4), height = nrow(hm_mat) * 0.5 + 3)
  pheatmap(hm_mat,
           color = colorRampPalette(c("white", "lightblue", "steelblue", "darkblue"))(50),
           cluster_rows = FALSE, cluster_cols = TRUE,
           main = sprintf("%s: Microenvironment × top celltypes (mean RCTD weight)", organ),
           fontsize_row = 11, fontsize_col = 9,
           display_numbers = TRUE, number_format = "%.2f",
           fontsize_number = 7, angle_col = 45)
  dev.off()
  cat("Saved: heatmap pdf\n")
  
  # ---------- 5. Dotplot ----------
  # 每微环境 top 3 + 在所有微环境的 mean weight (横轴 microenv, 纵轴 celltype, 点大小=mean_weight)
  dot_df <- me_ct %>%
    filter(celltype %in% top_cts_per_me) %>%
    left_join(spec %>% select(microenv_id, celltype, log2FC),
              by = c("microenv_id", "celltype")) %>%
    mutate(microenv_full = paste0(microenv_id, "\n", microenv_label),
           celltype = factor(celltype, levels = rev(unique(celltype))))
  
  p_dot <- ggplot(dot_df, aes(x = microenv_full, y = celltype)) +
    geom_point(aes(size = mean_weight, color = log2FC)) +
    scale_size_continuous(range = c(0.5, 8), name = "Mean weight") +
    scale_color_gradient2(low = "blue", mid = "white", high = "red",
                          midpoint = 0, name = "log2FC vs global") +
    labs(title = sprintf("%s microenvironment celltype composition", organ),
         x = "Microenvironment", y = "Celltype") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9),
          panel.grid.minor = element_blank())
  ggsave(file.path(OUT_DIR, paste0(organ, "_microenv_celltype_dotplot.pdf")),
         p_dot, width = max(6, length(unique(dot_df$microenv_id)) * 1.4 + 2),
         height = max(4, length(unique(dot_df$celltype)) * 0.35 + 2))
  cat("Saved: dotplot pdf\n")
  
  # ---------- 6. tumor vs baseline 差异 (在每微环境内, top celltype 的 weight) ----------
  cat("\nComputing tumor vs baseline difference per microenv...\n")
  tb_diff <- combined %>%
    select(microenv_id, group, all_of(ct_cols)) %>%
    pivot_longer(cols = -c(microenv_id, group),
                 names_to = "celltype", values_to = "weight") %>%
    group_by(microenv_id, group, celltype) %>%
    summarise(mean_weight = mean(weight), .groups = "drop") %>%
    pivot_wider(names_from = group, values_from = mean_weight, values_fill = 0) %>%
    mutate(
      log2FC_tumor_vs_baseline = log2((tumor + 1e-4) / (baseline + 1e-4)),
      microenv_label = me_to_label[microenv_id]
    ) %>%
    arrange(microenv_id, desc(abs(log2FC_tumor_vs_baseline)))
  
  out_tb <- file.path(OUT_DIR, paste0(organ, "_microenv_tumor_vs_baseline.csv"))
  write.csv(tb_diff, out_tb, row.names = FALSE)
  cat(sprintf("Saved: %s\n", out_tb))
  
  return(list(
    organ = organ,
    n_microenv = length(unique(me_naming$microenv_id)),
    n_spots = nrow(combined),
    top5 = top5,
    spec_top = spec_top,
    tb_diff = tb_diff
  ))
}

# ---------- Run ----------
results <- list()
results[["liver"]] <- process_organ("liver")
results[["brain"]] <- process_organ("brain")

# ---------- Summary ----------
md_file <- file.path(OUT_DIR, "summary.md")
con <- file(md_file, "w")
writeLines(c(
  "# Microenvironment Key Celltypes — Phase 2",
  "",
  paste0("**Generated:** ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## Method",
  "",
  "1. RCTD weights (KNN smoothed k=15) → cosine similarity assign D01-D10",
  "2. D01-D10 mapped to microenvironments via biological merge",
  "3. Per microenv: mean RCTD weight per celltype across all spots from 5 samples",
  "4. Specificity: log2FC vs global mean (across all microenvs in same organ)",
  "5. Tumor vs baseline: separate group analysis within each microenv",
  "",
  "## Outputs"
), con)

for (organ in c("liver", "brain")) {
  r <- results[[organ]]
  writeLines(c(
    "",
    paste0("### ", toupper(organ),
           " (", r$n_microenv, " microenvs, ", r$n_spots, " spots)"),
    "",
    "**Key celltypes (top 3 enriched vs global)**:",
    ""
  ), con)
  
  spec_print <- r$spec_top %>%
    group_by(microenv_id, microenv_label) %>%
    summarise(top_enriched = paste(sprintf("%s (log2FC=%.1f)",
                                            celltype, log2FC),
                                    collapse = ", "),
              .groups = "drop")
  
  for (i in seq_len(nrow(spec_print))) {
    writeLines(paste0("- **", spec_print$microenv_id[i], " ",
                      spec_print$microenv_label[i], "**: ",
                      spec_print$top_enriched[i]), con)
  }
}

writeLines(c(
  "",
  "## Files per organ",
  "",
  "- `<organ>_microenv_celltypes_full.csv` — microenv × CT × mean_weight (full)",
  "- `<organ>_microenv_celltypes_top5.csv` — top 5 CT per microenv",
  "- `<organ>_microenv_specificity.csv` — log2FC vs global mean (which CT is microenv-specific)",
  "- `<organ>_microenv_tumor_vs_baseline.csv` — tumor vs baseline diff per microenv",
  "- `<organ>_microenv_celltype_heatmap.pdf` — visual",
  "- `<organ>_microenv_celltype_dotplot.pdf` — visual",
  "",
  "## Next: Phase 3 marker genes (FindMarkers per microenv)"
), con)
close(con)
cat("\nSaved:", md_file, "\n")

cat("\n=============================================================\n")
cat("Phase 2 done. Outputs in:", OUT_DIR, "\n")
cat("=============================================================\n")
