#!/usr/bin/env Rscript
# ============================================================
# spatial_domain_by_composition.R  --  T4
# ============================================================
# 
# 目的: 用 RCTD weights matrix 的 k-means 聚类识别 spatial domain
#       (BayesSpace 替代方案, 基于 cell type composition 而非表达)
# 
# 逻辑:
#   1. 读 RCTD weights matrix (spot × cell type)
#   2. 对每个器官 (liver/brain), 合并 5 样本 weights
#   3. k-means (k=10 liver, k=10 brain), 得到每 spot 的 domain ID
#   4. 按 top 3 cell type 命名 domain
#   5. 算 tumor vs baseline 的 domain proportion 变化
#   6. 画: (a) domain × cell type heatmap, (b) per-sample domain map
# 
# 输出:
#   taskB/spatial_domain/
#     liver_domains_matrix.csv            domain × cell type mean weights
#     liver_domain_naming.csv             domain ID -> auto name
#     liver_per_sample_proportion.csv     sample × domain count
#     liver_tumor_vs_baseline.csv         domain shift (tumor - baseline)
#     liver_heatmap.pdf                   domain × ct heatmap
#     liver_per_sample_stacked.pdf        sample × domain stacked bar
#     (brain_ 同上 6 个文件)
#     summary.md
# 
# 用法:
#   cd ~/SGS10_analysis
#   screen -S t4 -dm bash -c '
#     Rscript scripts/spatial_domain_by_composition.R 2>&1 | \
#       tee scripts/logs/spatial_domain.log
#   '
# 预计: 10-15 min
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

PROJECT_ROOT  <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
RCTD_DIR      <- file.path(PROJECT_ROOT, "spatial_results", "rctd", "per_sample")
OUT_DIR       <- file.path(PROJECT_ROOT, "taskB", "spatial_domain")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

set.seed(42)

# ============================================================
# Samples per organ
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

cat("============================================================\n")
cat("T4: Spatial Domain by RCTD Composition (k-means)\n")
cat(sprintf("Start: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("============================================================\n")

# ============================================================
# Function: run domain analysis for one organ
# ============================================================
run_organ <- function(organ) {
  cat(sprintf("\n\n### === %s ===\n", toupper(organ)))
  
  samples      <- c(SAMPLES[[organ]]$tumor, SAMPLES[[organ]]$baseline)
  tumor_samples <- SAMPLES[[organ]]$tumor
  base_samples  <- SAMPLES[[organ]]$baseline
  k             <- K_CLUSTERS[[organ]]
  
  # --------------------------------------------------------
  # 1. Load weights matrices, concat
  # --------------------------------------------------------
  cat(sprintf("\n[1/6] Loading %d samples' weights matrices\n", length(samples)))
  all_weights <- list()
  
  for (s in samples) {
    w_path <- file.path(RCTD_DIR, paste0(s, "_weights_matrix.csv"))
    if (!file.exists(w_path)) {
      cat(sprintf("  ⚠ Missing: %s\n", w_path))
      next
    }
    w <- read.csv(w_path, stringsAsFactors = FALSE, row.names = 1)
    # Normalize per row (each spot's weights sum to 1)
    row_sums <- rowSums(w, na.rm = TRUE)
    w <- w / pmax(row_sums, 1e-10)
    w[is.na(w)] <- 0
    
    # tag sample
    w$spot_id <- paste0(s, "__", rownames(w))
    w$sample <- s
    w$group <- ifelse(s %in% tumor_samples, "tumor", "baseline")
    all_weights[[s]] <- w
    cat(sprintf("  %s: %d spots, %d cell types\n",
                s, nrow(w), ncol(w) - 3))
  }
  
  # --------------------------------------------------------
  # 2. Align cell type columns (union across samples)
  # --------------------------------------------------------
  cat(sprintf("\n[2/6] Aligning cell type columns\n"))
  meta_cols <- c("spot_id", "sample", "group")
  ct_cols_all <- unique(unlist(lapply(all_weights,
                                      function(w) setdiff(colnames(w), meta_cols))))
  cat(sprintf("  Union cell types: %d\n", length(ct_cols_all)))
  
  for (s in names(all_weights)) {
    w <- all_weights[[s]]
    missing_ct <- setdiff(ct_cols_all, colnames(w))
    for (ct in missing_ct) w[[ct]] <- 0
    all_weights[[s]] <- w[, c(meta_cols, ct_cols_all)]
  }
  
  combined <- bind_rows(all_weights)
  cat(sprintf("  Combined: %d spots total\n", nrow(combined)))
  
  meta_df <- combined[, meta_cols]
  weight_mat <- as.matrix(combined[, ct_cols_all])
  
  # --------------------------------------------------------
  # 3. k-means clustering
  # --------------------------------------------------------
  cat(sprintf("\n[3/6] K-means clustering (k=%d)\n", k))
  t0 <- Sys.time()
  km <- kmeans(weight_mat, centers = k, nstart = 25, iter.max = 100)
  cat(sprintf("  done in %.1fs, tot.withinss = %.1f\n",
              as.numeric(Sys.time()-t0), km$tot.withinss))
  
  meta_df$domain <- paste0("D", sprintf("%02d", km$cluster))
  
  # --------------------------------------------------------
  # 4. Domain × cell type matrix + auto-naming by top 3 CT
  # --------------------------------------------------------
  cat(sprintf("\n[4/6] Domain × cell type profile\n"))
  domain_profile <- combined %>%
    mutate(domain = meta_df$domain) %>%
    group_by(domain) %>%
    summarise(across(all_of(ct_cols_all), mean), .groups = "drop")
  
  domain_mat <- as.data.frame(domain_profile[, -1])
  rownames(domain_mat) <- domain_profile$domain
  
  # Auto-name: top 3 cell types per domain
  naming <- data.frame(
    domain = rownames(domain_mat),
    n_spots = as.integer(table(meta_df$domain)[rownames(domain_mat)]),
    top1 = NA_character_, top1_w = NA_real_,
    top2 = NA_character_, top2_w = NA_real_,
    top3 = NA_character_, top3_w = NA_real_,
    auto_name = NA_character_,
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(domain_mat))) {
    row <- domain_mat[i, ]
    ord <- order(as.numeric(row), decreasing = TRUE)
    top3_idx <- ord[1:min(3, length(ord))]
    naming$top1[i]   <- colnames(domain_mat)[top3_idx[1]]
    naming$top1_w[i] <- round(as.numeric(row[top3_idx[1]]), 3)
    if (length(top3_idx) >= 2) {
      naming$top2[i]   <- colnames(domain_mat)[top3_idx[2]]
      naming$top2_w[i] <- round(as.numeric(row[top3_idx[2]]), 3)
    }
    if (length(top3_idx) >= 3) {
      naming$top3[i]   <- colnames(domain_mat)[top3_idx[3]]
      naming$top3_w[i] <- round(as.numeric(row[top3_idx[3]]), 3)
    }
    # Compact auto name: use top1, add top2 if close
    nm <- naming$top1[i]
    if (!is.na(naming$top2_w[i]) && naming$top2_w[i] > 0.6 * naming$top1_w[i]) {
      nm <- paste(naming$top1[i], naming$top2[i], sep = "+")
    }
    naming$auto_name[i] <- sprintf("%s [%s]", naming$domain[i], nm)
  }
  
  cat("  Domain naming:\n")
  for (i in seq_len(nrow(naming))) {
    cat(sprintf("    %s: %-40s (n=%d)\n",
                naming$domain[i], naming$auto_name[i], naming$n_spots[i]))
  }
  
  write.csv(domain_mat, file.path(OUT_DIR, paste0(organ, "_domain_ct_profile.csv")))
  write.csv(naming, file.path(OUT_DIR, paste0(organ, "_domain_naming.csv")),
            row.names = FALSE)
  
  # --------------------------------------------------------
  # 5. Per-sample domain proportion + tumor vs baseline
  # --------------------------------------------------------
  cat(sprintf("\n[5/6] Per-sample domain proportion\n"))
  per_sample <- meta_df %>%
    group_by(sample, group, domain) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(sample) %>%
    mutate(pct = 100 * n / sum(n)) %>%
    ungroup()
  
  write.csv(per_sample, file.path(OUT_DIR, paste0(organ, "_per_sample_domain.csv")),
            row.names = FALSE)
  
  # Tumor vs baseline shift
  shift <- per_sample %>%
    group_by(group, domain) %>%
    summarise(mean_pct = mean(pct), .groups = "drop") %>%
    pivot_wider(names_from = group, values_from = mean_pct,
                values_fill = 0) %>%
    mutate(delta = tumor - baseline,
           fc = (tumor + 0.01) / (baseline + 0.01)) %>%
    arrange(desc(abs(delta))) %>%
    left_join(naming %>% select(domain, auto_name), by = "domain")
  
  write.csv(shift, file.path(OUT_DIR, paste0(organ, "_tumor_vs_baseline.csv")),
            row.names = FALSE)
  
  cat("  Tumor vs baseline domain shift (top 5 by |delta|):\n")
  print(head(shift, 5))
  
  # --------------------------------------------------------
  # 6. Plots
  # --------------------------------------------------------
  cat(sprintf("\n[6/6] Plotting\n"))
  
  # 6a. Domain × cell type heatmap
  # Restrict to cell types with max weight > 0.05 for readability
  ct_max <- apply(domain_mat, 2, max)
  ct_keep <- names(ct_max)[ct_max > 0.05]
  if (length(ct_keep) < 5) ct_keep <- names(sort(ct_max, decreasing = TRUE))[1:10]
  
  dm_sub <- as.matrix(domain_mat[, ct_keep, drop = FALSE])
  
  # Label rows by auto_name
  rownames(dm_sub) <- naming$auto_name[match(rownames(dm_sub), naming$domain)]
  
  pdf(file.path(OUT_DIR, paste0(organ, "_heatmap.pdf")), width = 10, height = 6.5)
  pheatmap(dm_sub,
           color = colorRampPalette(c("white", "orange", "red"))(100),
           cluster_rows = TRUE, cluster_cols = TRUE,
           display_numbers = TRUE, number_format = "%.2f",
           fontsize_number = 7,
           main = sprintf("%s spatial domains: mean cell type weights",
                         toupper(organ)))
  dev.off()
  cat(sprintf("  -> %s_heatmap.pdf\n", organ))
  
  # 6b. Per-sample stacked bar chart
  per_sample_wide <- per_sample %>%
    left_join(naming %>% select(domain, auto_name), by = "domain") %>%
    mutate(sample_label = paste0(sample, "\n[", group, "]"))
  
  p_stack <- ggplot(per_sample_wide,
                    aes(x = sample_label, y = pct, fill = auto_name)) +
    geom_col(position = "stack", color = "white", size = 0.2) +
    scale_fill_manual(values = colorRampPalette(brewer.pal(10, "Set3"))(
      length(unique(per_sample_wide$auto_name)))) +
    labs(title = sprintf("%s spatial domain composition per sample",
                         toupper(organ)),
         x = NULL, y = "Percent of spots (%)",
         fill = "Domain") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 0),
          legend.text = element_text(size = 8),
          legend.key.size = unit(0.4, "cm"))
  
  ggsave(file.path(OUT_DIR, paste0(organ, "_per_sample_stacked.pdf")),
         p_stack, width = 10, height = 6)
  ggsave(file.path(OUT_DIR, paste0(organ, "_per_sample_stacked.png")),
         p_stack, width = 10, height = 6, dpi = 200)
  cat(sprintf("  -> %s_per_sample_stacked.{pdf,png}\n", organ))
  
  # 6c. Tumor vs baseline bar plot
  shift_plot <- shift %>%
    mutate(direction = ifelse(delta > 0, "Expanded in tumor", "Reduced in tumor"))
  
  p_shift <- ggplot(shift_plot, aes(x = reorder(auto_name, delta),
                                     y = delta, fill = direction)) +
    geom_col(alpha = 0.85) +
    geom_hline(yintercept = 0, color = "gray40") +
    geom_text(aes(label = sprintf("%+.1f%%", delta)),
              hjust = ifelse(shift_plot$delta >= 0, -0.1, 1.1), size = 3) +
    coord_flip() +
    scale_fill_manual(values = c("Expanded in tumor" = "#d62728",
                                 "Reduced in tumor" = "#1f77b4")) +
    labs(title = sprintf("%s: domain proportion shift (tumor − baseline)",
                         toupper(organ)),
         subtitle = sprintf("n=%d tumor vs %d baseline samples",
                            length(tumor_samples), length(base_samples)),
         x = NULL, y = "Δ proportion (% points)") +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom")
  
  ggsave(file.path(OUT_DIR, paste0(organ, "_tumor_vs_baseline.pdf")),
         p_shift, width = 9, height = 6)
  cat(sprintf("  -> %s_tumor_vs_baseline.pdf\n", organ))
  
  # Return info for summary
  list(organ = organ, k = k, n_spots = nrow(meta_df),
       top_shift = shift[1:min(5, nrow(shift)), ])
}

# ============================================================
# Run for both organs
# ============================================================
results <- list()
for (organ in names(SAMPLES)) {
  tryCatch({
    results[[organ]] <- run_organ(organ)
  }, error = function(e) {
    cat(sprintf("\n⚠ ERROR in %s: %s\n", organ, conditionMessage(e)))
    results[[organ]] <<- NULL
  })
}

# ============================================================
# Summary.md
# ============================================================
cat("\n### Writing summary.md ###\n")

md <- c(
  "# T4 Spatial Domain by RCTD Composition",
  "",
  sprintf("**Generated:** %s", format(Sys.time(), "%Y-%m-%d %H:%M")),
  "",
  "## Method",
  "",
  "1. Load RCTD weights matrix (spot × cell type) for each sample",
  "2. Normalize per-spot to sum = 1",
  "3. Concat 5 samples per organ, k-means clustering (k=10)",
  "4. Auto-name each domain by top 3 cell types",
  "5. Compare tumor vs baseline domain proportions",
  "",
  "## Results"
)

for (organ in names(results)) {
  if (is.null(results[[organ]])) next
  r <- results[[organ]]
  md <- c(md, "",
    sprintf("### %s (k=%d, %d total spots)",
            toupper(r$organ), r$k, r$n_spots),
    "",
    "Top 5 tumor-shift domains:",
    "",
    "| Domain | Auto name | Tumor % | Baseline % | Δ |",
    "|--------|-----------|--------:|-----------:|--:|"
  )
  for (i in seq_len(nrow(r$top_shift))) {
    rs <- r$top_shift[i, ]
    md <- c(md, sprintf("| %s | %s | %.1f | %.1f | %+.1f |",
                        rs$domain, rs$auto_name, rs$tumor, rs$baseline, rs$delta))
  }
}

md <- c(md, "",
  "## Output files per organ",
  "",
  "- `<organ>_domain_ct_profile.csv`: domain × cell type mean weights",
  "- `<organ>_domain_naming.csv`: domain auto-names (top 3 CT)",
  "- `<organ>_per_sample_domain.csv`: sample × domain proportion",
  "- `<organ>_tumor_vs_baseline.csv`: tumor-shift per domain",
  "- `<organ>_heatmap.pdf`: domain × CT heatmap",
  "- `<organ>_per_sample_stacked.{pdf,png}`: sample domain composition",
  "- `<organ>_tumor_vs_baseline.pdf`: shift bar plot",
  "",
  "## Interpretation",
  "",
  "Domains with positive Δ are tumor-induced niches.",
  "Cross-reference domain auto-names with Task A v1.1 claims:",
  "- **Liver**: HSC_Activated_ECM+myeloid domain → L1/L4 validation",
  "- **Brain**: Microglia+Astrocyte / Endothelial+Pericyte → N1-N4 validation"
)

writeLines(md, file.path(OUT_DIR, "summary.md"))

cat("\n============================================================\n")
cat("T4 COMPLETE\n")
cat(sprintf("Output: %s\n", OUT_DIR))
cat(sprintf("Files:  %d\n", length(list.files(OUT_DIR))))
cat("============================================================\n")
