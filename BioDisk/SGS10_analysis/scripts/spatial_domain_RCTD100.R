#!/usr/bin/env Rscript
# ============================================================
# spatial_domain_RCTD100.R
# ============================================================
# 功能:
#   1. 加载 100% RCTD weights matrices (10 samples)
#   2. k-means 聚类 (k=10) 识别 10 个 spatial domain (D01-D10)
#   3. 按预设映射给 domain 命名为 microenv (M01-M07 liver, M01-M06 brain)
#   4. 输出 domain × cell type profile + per-sample distribution + heatmap
#
# 输入:
#   - spatial_results/rctd_full/per_sample/<sample>_weights_matrix.csv (10 个)
#
# 输出:
#   taskB/spatial_domain_RCTD100/
#     liver_domain_ct_profile.csv    domain × cell type mean weights
#     liver_domain_naming.csv        domain → top1 + auto_name
#     liver_microenv_naming.csv      domain → microenv_id + label (paper M01-M07)
#     liver_per_sample_domain.csv    sample × domain spot count
#     liver_tumor_vs_baseline.csv    domain proportion shift
#     liver_heatmap.pdf              domain × cell type heatmap
#     liver_per_sample_stacked.pdf   sample × domain stacked bar
#     brain_*  (同上 7 个文件)
#     summary.md
#
# 预计: 5-15 min
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
})

select    <- dplyr::select
filter    <- dplyr::filter
mutate    <- dplyr::mutate
summarise <- dplyr::summarise
group_by  <- dplyr::group_by
bind_rows <- dplyr::bind_rows

PROJECT_ROOT <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
RCTD_DIR     <- file.path(PROJECT_ROOT, "spatial_results", "rctd_full", "per_sample")
OUT_DIR      <- file.path(PROJECT_ROOT, "taskB", "spatial_domain_RCTD100")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

set.seed(42)

# ============================================================
# Sample groups
# ============================================================
SAMPLES <- list(
  liver = list(
    tumor    = c("Liver_3", "Liver_8", "gan56"),
    baseline = c("Bla_liver_1", "Con_liver_1")
  ),
  brain = list(
    tumor    = c("Brain_3", "Brain_8", "nao56"),
    baseline = c("Bla_brain_1", "Con_brain_1")
  )
)

K_CLUSTERS <- list(liver = 10, brain = 10)

# ============================================================
# Microenv mapping (domain new_id → microenv_id + label)
# 基于对齐分析: 100% RCTD k-means 出的 D 跟 v1 D 几乎完全对应
# 微调: D_new_05 (LSEC_Mid 60%) 划入 M03 (sinusoidal_endothelium)
#       而非 v1 的 M05 (因为它是过渡 LSEC, 非 tumor-induced)
# Paper-key M05 = D_new_01 (Mac_M1+HSC+Kupffer+cDC2+immune, 85.5% tumor)
# ============================================================
LIVER_MICROENV <- data.frame(
  domain         = c("D_new_03","D_new_07","D_new_02","D_new_06",
                     "D_new_09","D_new_05","D_new_08","D_new_01",
                     "D_new_04","D_new_10"),
  microenv_id    = c("M01","M01","M02","M02","M03","M03","M04","M05","M06","M07"),
  microenv_label = c("Hepatocyte_zonation","Hepatocyte_zonation",
                     "Vascular_perisinusoidal","Vascular_perisinusoidal",
                     "Sinusoidal_endothelium","Sinusoidal_endothelium",
                     "Stellate_quiescent",
                     "Tumor_immune_niche",
                     "B_plasma",
                     "Kupffer"),
  stringsAsFactors = FALSE
)

BRAIN_MICROENV <- data.frame(
  domain         = c("D_new_09","D_new_06","D_new_02","D_new_07","D_new_08",
                     "D_new_04","D_new_10","D_new_05","D_new_01","D_new_03"),
  microenv_id    = c("M01","M01","M01","M02","M02","M03","M04","M05","M05","M06"),
  microenv_label = c("Hippocampus","Hippocampus","Hippocampus",
                     "Cortex_excitatory","Cortex_excitatory",
                     "Cortex_inhibitory",
                     "Vascular_BBB",
                     "Glia","Glia",
                     "General_neurons"),
  stringsAsFactors = FALSE
)

cat("============================================================\n")
cat("Spatial Domain RCTD100\n")
cat(sprintf("Start: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("============================================================\n")

# ============================================================
# Run analysis for one organ
# ============================================================
run_organ <- function(organ, microenv_table) {
  cat(sprintf("\n\n### === %s ===\n", toupper(organ)))
  
  samples       <- c(SAMPLES[[organ]]$tumor, SAMPLES[[organ]]$baseline)
  tumor_samples <- SAMPLES[[organ]]$tumor
  base_samples  <- SAMPLES[[organ]]$baseline
  k             <- K_CLUSTERS[[organ]]
  
  # ---- 1. Load weights matrices ----
  cat(sprintf("\n[1/7] Loading %d samples\n", length(samples)))
  all_weights <- list()
  for (s in samples) {
    f <- file.path(RCTD_DIR, sprintf("%s_weights_matrix.csv", s))
    if (!file.exists(f)) { cat(sprintf("  miss: %s\n", s)); next }
    w <- read.csv(f, stringsAsFactors = FALSE, row.names = 1)
    rs <- rowSums(w, na.rm = TRUE)
    w <- w / pmax(rs, 1e-10)
    w[is.na(w)] <- 0
    w$spot_id <- paste0(s, "__", rownames(w))
    w$sample  <- s
    w$group   <- ifelse(s %in% tumor_samples, "tumor", "baseline")
    all_weights[[s]] <- w
    cat(sprintf("  %s: %d spots × %d ct\n", s, nrow(w), ncol(w) - 3))
  }
  
  # ---- 2. Align cell type columns ----
  cat(sprintf("\n[2/7] Aligning cell type columns\n"))
  meta_cols <- c("spot_id", "sample", "group")
  ct_cols <- unique(unlist(lapply(all_weights, function(w) setdiff(colnames(w), meta_cols))))
  cat(sprintf("  Union cell types: %d\n", length(ct_cols)))
  
  for (s in names(all_weights)) {
    miss <- setdiff(ct_cols, colnames(all_weights[[s]]))
    for (c in miss) all_weights[[s]][[c]] <- 0
    all_weights[[s]] <- all_weights[[s]][, c(meta_cols, ct_cols)]
  }
  combined <- bind_rows(all_weights)
  cat(sprintf("  Combined: %d spots\n", nrow(combined)))
  
  meta_df    <- combined[, meta_cols]
  weight_mat <- as.matrix(combined[, ct_cols])
  
  # ---- 3. k-means clustering ----
  cat(sprintf("\n[3/7] K-means clustering (k=%d)\n", k))
  t0 <- Sys.time()
  km <- kmeans(weight_mat, centers = k, nstart = 25, iter.max = 100)
  cat(sprintf("  done in %.1fs (tot.withinss = %.1f)\n",
              as.numeric(Sys.time() - t0), km$tot.withinss))
  
  meta_df$domain <- paste0("D_new_", sprintf("%02d", km$cluster))
  
  # ---- 4. Domain × cell type mean profile ----
  cat(sprintf("\n[4/7] Computing domain × cell type profile\n"))
  prof_long <- cbind(meta_df, weight_mat) %>%
    select(domain, all_of(ct_cols)) %>%
    group_by(domain) %>%
    summarise(across(everything(), mean), .groups = "drop")
  
  prof_mat <- as.matrix(prof_long[, -1])
  rownames(prof_mat) <- prof_long$domain
  
  # 写文件
  prof_df <- as.data.frame(prof_mat)
  write.csv(prof_df, file.path(OUT_DIR, sprintf("%s_domain_ct_profile.csv", organ)),
            row.names = TRUE)
  cat(sprintf("  Saved: %s_domain_ct_profile.csv\n", organ))
  
  # ---- 5. Domain naming (top 1-3 cell type) ----
  cat(sprintf("\n[5/7] Domain naming\n"))
  domain_size <- meta_df %>% count(domain, name = "n_spots")
  
  domain_naming <- data.frame(
    domain = rownames(prof_mat),
    top1 = colnames(prof_mat)[apply(prof_mat, 1, which.max)],
    top1_pct = round(apply(prof_mat, 1, max) * 100, 1),
    top2 = sapply(seq_len(nrow(prof_mat)), function(i) {
      x <- prof_mat[i, ]
      colnames(prof_mat)[order(x, decreasing = TRUE)[2]]
    }),
    top2_pct = sapply(seq_len(nrow(prof_mat)), function(i) {
      x <- prof_mat[i, ]
      round(sort(x, decreasing = TRUE)[2] * 100, 1)
    }),
    top3 = sapply(seq_len(nrow(prof_mat)), function(i) {
      x <- prof_mat[i, ]
      colnames(prof_mat)[order(x, decreasing = TRUE)[3]]
    }),
    top3_pct = sapply(seq_len(nrow(prof_mat)), function(i) {
      x <- prof_mat[i, ]
      round(sort(x, decreasing = TRUE)[3] * 100, 1)
    }),
    stringsAsFactors = FALSE
  ) %>% left_join(domain_size, by = "domain")
  
  domain_naming$auto_name <- sprintf("%s [%s]", domain_naming$domain, domain_naming$top1)
  
  write.csv(domain_naming, file.path(OUT_DIR, sprintf("%s_domain_naming.csv", organ)),
            row.names = FALSE)
  cat(sprintf("  Saved: %s_domain_naming.csv\n", organ))
  print(domain_naming)
  
  # ---- 6. Microenv mapping (domain → M_id) ----
  cat(sprintf("\n[6/7] Microenv mapping\n"))
  microenv_full <- microenv_table %>%
    left_join(domain_naming %>% select(domain, top1, n_spots), by = "domain") %>%
    mutate(auto_name = sprintf("%s [%s]", domain, top1),
           microenv_name = sprintf("%s %s", microenv_id, microenv_label))
  
  write.csv(microenv_full, file.path(OUT_DIR, sprintf("%s_microenv_naming.csv", organ)),
            row.names = FALSE)
  cat(sprintf("  Saved: %s_microenv_naming.csv\n", organ))
  print(microenv_full)
  
  # ---- 7. Per-sample × domain distribution ----
  cat(sprintf("\n[7/7] Per-sample distribution + heatmap\n"))
  per_sample <- meta_df %>%
    count(sample, domain) %>%
    pivot_wider(names_from = domain, values_from = n, values_fill = 0)
  write.csv(per_sample, file.path(OUT_DIR, sprintf("%s_per_sample_domain.csv", organ)),
            row.names = FALSE)
  
  # Tumor vs baseline
  tvb <- meta_df %>%
    count(group, domain) %>%
    group_by(group) %>%
    mutate(pct = n / sum(n) * 100) %>%
    select(group, domain, pct) %>%
    pivot_wider(names_from = group, values_from = pct, values_fill = 0)
  if (all(c("tumor","baseline") %in% colnames(tvb))) {
    tvb$diff <- tvb$tumor - tvb$baseline
    tvb$log2fc <- log2((tvb$tumor + 0.5) / (tvb$baseline + 0.5))
    tvb <- tvb %>% arrange(desc(log2fc))
  }
  
  # 加上 microenv label 帮助阅读
  tvb <- tvb %>% left_join(microenv_full %>% select(domain, microenv_name), by = "domain")
  
  write.csv(tvb, file.path(OUT_DIR, sprintf("%s_tumor_vs_baseline.csv", organ)),
            row.names = FALSE)
  cat(sprintf("  Saved: %s_tumor_vs_baseline.csv\n", organ))
  print(tvb)
  
  # Heatmap (domain × cell type)
  cat(sprintf("\n  Drawing heatmap...\n"))
  # 行序: 按 microenv_id 排序
  domain_order <- microenv_full %>% arrange(microenv_id, domain) %>% pull(domain)
  prof_ordered <- prof_mat[domain_order, ]
  
  # 列序: 按 max weight 排序 (top 25 cell types)
  ct_max <- apply(prof_ordered, 2, max)
  top_ct <- names(sort(ct_max, decreasing = TRUE))[1:min(25, length(ct_max))]
  prof_top <- prof_ordered[, top_ct]
  
  # 注释行
  row_anno <- data.frame(
    microenv = microenv_full$microenv_id[match(rownames(prof_top), microenv_full$domain)],
    row.names = rownames(prof_top)
  )
  
  pheatmap(prof_top,
           color = colorRampPalette(c("white","yellow","red","darkred"))(50),
           cluster_rows = FALSE, cluster_cols = FALSE,
           annotation_row = row_anno,
           main = sprintf("%s — Domain × Cell Type Profile (top 25 ct)", toupper(organ)),
           filename = file.path(OUT_DIR, sprintf("%s_heatmap.pdf", organ)),
           width = 12, height = 6,
           fontsize_row = 8, fontsize_col = 8)
  cat(sprintf("  Saved: %s_heatmap.pdf\n", organ))
  
  # Stacked bar (sample × domain)
  ps_long <- per_sample %>% pivot_longer(-sample, names_to = "domain", values_to = "n")
  ps_long$domain <- factor(ps_long$domain, levels = domain_order)
  ps_long <- ps_long %>% left_join(microenv_full %>% select(domain, microenv_id),
                                     by = "domain")
  
  p <- ggplot(ps_long, aes(x = sample, y = n, fill = domain)) +
    geom_col(position = "fill") +
    scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Set2"))(k)) +
    labs(title = sprintf("%s — Sample × Domain Composition", toupper(organ)),
         x = "Sample", y = "Proportion") +
    theme_classic(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(OUT_DIR, sprintf("%s_per_sample_stacked.pdf", organ)), p,
         width = 10, height = 6)
  ggsave(file.path(OUT_DIR, sprintf("%s_per_sample_stacked.png", organ)), p,
         width = 10, height = 6, dpi = 150)
  cat(sprintf("  Saved: %s_per_sample_stacked.pdf/.png\n", organ))
  
  # 返回 spot-level domain assignment (供下游 Phase 3 + H&E 用)
  list(
    meta_df       = cbind(meta_df, weight_mat),  # spot-level: spot_id, sample, group, domain, ct weights
    microenv_full = microenv_full,
    prof_mat      = prof_mat
  )
}

# ============================================================
# Run
# ============================================================
liver_res <- run_organ("liver", LIVER_MICROENV)
brain_res <- run_organ("brain", BRAIN_MICROENV)

# ============================================================
# Save spot-level assignment (for Phase 3 + H&E overlay)
# ============================================================
cat("\n\n[Saving spot-level domain assignment]\n")

spot_assign <- function(res, organ) {
  df <- res$meta_df %>%
    left_join(res$microenv_full %>% select(domain, microenv_id, microenv_label, microenv_name),
              by = "domain")
  
  # 只保留必要列 (spot_id, sample, group, domain, microenv_id, microenv_label)
  keep_cols <- c("spot_id", "sample", "group", "domain",
                 "microenv_id", "microenv_label", "microenv_name")
  df_meta_only <- df[, keep_cols]
  
  out_f <- file.path(OUT_DIR, sprintf("%s_spot_microenv_assignment.csv", organ))
  write.csv(df_meta_only, out_f, row.names = FALSE)
  cat(sprintf("  Saved: %s (%d spots)\n", basename(out_f), nrow(df_meta_only)))
}

spot_assign(liver_res, "liver")
spot_assign(brain_res, "brain")

# ============================================================
# Summary
# ============================================================
summary_md <- c(
  "# Spatial Domain RCTD100 Summary",
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## Liver",
  sprintf("- Total spots: %d", nrow(liver_res$meta_df)),
  sprintf("- Cell types: %d", ncol(liver_res$prof_mat)),
  sprintf("- Domains (k-means k=10): %s", paste(rownames(liver_res$prof_mat), collapse = ", ")),
  sprintf("- Microenvs: %s", 
          paste(unique(liver_res$microenv_full$microenv_id), collapse = ", ")),
  "",
  "## Brain",
  sprintf("- Total spots: %d", nrow(brain_res$meta_df)),
  sprintf("- Cell types: %d", ncol(brain_res$prof_mat)),
  sprintf("- Domains: %s", paste(rownames(brain_res$prof_mat), collapse = ", ")),
  sprintf("- Microenvs: %s",
          paste(unique(brain_res$microenv_full$microenv_id), collapse = ", ")),
  "",
  "## Output files (per organ)",
  "- `<organ>_domain_ct_profile.csv` — domain × cell type mean weights",
  "- `<organ>_domain_naming.csv` — auto domain name (top 1-3 cell types)",
  "- `<organ>_microenv_naming.csv` — manual D → M mapping",
  "- `<organ>_per_sample_domain.csv` — sample × domain count",
  "- `<organ>_tumor_vs_baseline.csv` — domain proportion shift",
  "- `<organ>_heatmap.pdf` — domain × cell type heatmap",
  "- `<organ>_per_sample_stacked.pdf` — sample × domain composition",
  "- `<organ>_spot_microenv_assignment.csv` — spot-level domain + microenv (for downstream)"
)
writeLines(summary_md, file.path(OUT_DIR, "summary.md"))

cat(sprintf("\n\nDone: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(sprintf("Output: %s\n", OUT_DIR))
system(sprintf("ls -la %s", OUT_DIR))
