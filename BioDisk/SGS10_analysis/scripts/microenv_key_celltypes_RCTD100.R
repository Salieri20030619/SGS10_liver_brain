#!/usr/bin/env Rscript
source("/mnt/BioDisk/wangmingyuan/SGS10_analysis/scripts/fig_save.R")  # [restyle] 三格式输出器
# ============================================================
# microenv_key_celltypes_RCTD100.R
# ============================================================
# Purpose:
#   For each microenvironment (M01-M07 liver, M01-M06 brain),
#   identify key cell types by mean RCTD weight ranking.
#
# Method:
#   1. Load 10 sample RCTD weights + spot-level microenv assignment
#   2. Per microenv, compute mean cell type weight (across all spots)
#   3. Rank cell types per microenv (top 5)
#   4. Compute specificity: log2FC vs global mean weight
#   5. Compute tumor vs baseline shift per (microenv, celltype)
#
# Inputs:
#   - spatial_results/rctd_full/per_sample/<s>_weights_matrix.csv
#   - taskB/spatial_domain_RCTD100/<organ>_microenv_naming.csv
#   - taskB/spatial_domain_RCTD100/<organ>_spot_microenv_assignment.csv
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
})

PR <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
RCTD_DIR <- file.path(PR, "spatial_results/rctd_full/per_sample")
SD_DIR   <- file.path(PR, "taskB/spatial_domain_RCTD100")
OUT_DIR  <- file.path(PR, "D1_v2_RCTD_composition/microenv_key_celltypes_RCTD100")
RESTYLE_DIR <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis/figures/restyle_v1"  # [restyle] 新图统一目录
dir.create(RESTYLE_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

set.seed(42)

LIVER_SAMPLES <- c("Liver_3","Liver_8","gan56","Bla_liver_1","Con_liver_1")
BRAIN_SAMPLES <- c("Brain_3","Brain_8","nao56","Bla_brain_1","Con_brain_1")
TUMOR_SAMPLES <- c("Liver_3","Liver_8","gan56","Brain_3","Brain_8","nao56")

cat("=============================================================\n")
cat("Microenvironment key cell types\n")
cat("Time:", format(Sys.time()), "\n")
cat("Output:", OUT_DIR, "\n")
cat("=============================================================\n\n")

process_organ <- function(organ) {
  samples <- if (organ == "liver") LIVER_SAMPLES else BRAIN_SAMPLES
  cat(sprintf("\n##### %s #####\n", toupper(organ)))
  
  # ---- Load microenv naming + spot assignment ----
  me_naming <- read.csv(file.path(SD_DIR, sprintf("%s_microenv_naming.csv", organ)),
                         stringsAsFactors = FALSE)
  me_to_label <- setNames(me_naming$microenv_label, me_naming$microenv_id)
  me_to_label <- me_to_label[!duplicated(names(me_to_label))]
  
  spot_assign <- read.csv(file.path(SD_DIR, sprintf("%s_spot_microenv_assignment.csv", organ)),
                           stringsAsFactors = FALSE)
  cat(sprintf("Spot assignment: %d spots\n", nrow(spot_assign)))
  cat(sprintf("Spot assignment columns: %s\n", paste(colnames(spot_assign), collapse=",")))
  cat("Sample distribution in assignment:\n")
  print(table(spot_assign$sample))
  cat("Group distribution in assignment:\n")
  print(table(spot_assign$group))
  
  # spot_assign 已经有 sample, group, microenv_id 列, 不需要再加
  # 但 group 可能是 "tumor"/"baseline" 也可能是 "Tumor"/"Blank"/"Control"
  # → 我们重新规范化: tumor vs baseline (control 也算 baseline)
  spot_assign$group_norm <- ifelse(spot_assign$sample %in% TUMOR_SAMPLES,
                                    "tumor", "baseline")
  cat("Normalized group distribution:\n")
  print(table(spot_assign$group_norm))
  
  # ---- Load weights matrices, 用 sample__barcode 作为 spot_id ----
  cat("\nLoading RCTD weights...\n")
  weights_long <- list()
  for (s in samples) {
    f <- file.path(RCTD_DIR, sprintf("%s_weights_matrix.csv", s))
    w <- read.csv(f, stringsAsFactors = FALSE, row.names = 1)
    rs <- rowSums(w, na.rm = TRUE); w <- w / pmax(rs, 1e-10); w[is.na(w)] <- 0
    
    # spot_id 跟 spot_assign 一样: sample__barcode
    w$spot_id <- paste0(s, "__", rownames(w))
    weights_long[[s]] <- w
    cat(sprintf("  %s: %d spots × %d ct\n", s, nrow(w), ncol(w) - 1))
  }
  
  # 列对齐 (除 spot_id)
  all_ct <- unique(unlist(lapply(weights_long, function(w) setdiff(colnames(w), "spot_id"))))
  for (s in names(weights_long)) {
    miss <- setdiff(all_ct, colnames(weights_long[[s]]))
    for (c in miss) weights_long[[s]][[c]] <- 0
    weights_long[[s]] <- weights_long[[s]][, c("spot_id", all_ct)]
  }
  combined_w <- do.call(rbind, weights_long)
  cat(sprintf("Combined weights: %d spots × %d ct\n",
              nrow(combined_w), length(all_ct)))
  
  # 看 spot_id 例子
  cat(sprintf("\nspot_id format:\n"))
  cat(sprintf("  weights:    %s\n", combined_w$spot_id[1]))
  cat(sprintf("  assignment: %s\n", spot_assign$spot_id[1]))
  
  # Match
  common_spots <- intersect(combined_w$spot_id, spot_assign$spot_id)
  cat(sprintf("\nMatched spots (weights ∩ assignment): %d / %d\n",
              length(common_spots), nrow(spot_assign)))
  
  if (length(common_spots) == 0) {
    cat("\n⚠ NO MATCH! Check spot_id format consistency.\n")
    cat("First 5 weights spot_ids:\n")
    print(head(combined_w$spot_id, 5))
    cat("First 5 assignment spot_ids:\n")
    print(head(spot_assign$spot_id, 5))
    stop("No matching spot_ids")
  }
  
  # 子集 + 对齐顺序
  combined_w <- combined_w %>% filter(spot_id %in% common_spots)
  spot_assign <- spot_assign %>% filter(spot_id %in% common_spots)
  
  # 排序保证一致 (用 spot_id 排序)
  combined_w <- combined_w[order(combined_w$spot_id), ]
  spot_assign <- spot_assign[order(spot_assign$spot_id), ]
  
  # 验证一致性
  if (!identical(combined_w$spot_id, spot_assign$spot_id)) {
    stop("spot_id 顺序不一致 after merge sort")
  }
  
  # 提取 weight matrix (不含 spot_id 列)
  weight_mat <- as.matrix(combined_w[, all_ct])
  rownames(weight_mat) <- combined_w$spot_id
  
  cat(sprintf("Final: %d spots × %d cell types\n", nrow(weight_mat), ncol(weight_mat)))
  
  # ============ 1. Microenv × cell type mean weight ============
  cat("\n[1] Microenv × cell type mean weight\n")
  long_df <- data.frame(
    microenv_id = spot_assign$microenv_id,
    group_norm = spot_assign$group_norm,
    weight_mat,
    check.names = FALSE
  )
  
  me_ct <- long_df %>%
    select(microenv_id, all_of(all_ct)) %>%
    group_by(microenv_id) %>%
    summarise(across(everything(), mean), .groups = "drop") %>%
    pivot_longer(cols = -microenv_id, names_to = "celltype", values_to = "mean_weight")
  
  me_size <- long_df %>% count(microenv_id, name = "n_spots")
  
  me_ct <- me_ct %>% left_join(me_size, by = "microenv_id") %>%
    mutate(microenv_label = me_to_label[microenv_id]) %>%
    arrange(microenv_id, desc(mean_weight))
  
  out_full <- file.path(OUT_DIR, paste0(organ, "_microenv_celltypes_full.csv"))
  write.csv(me_ct, out_full, row.names = FALSE)
  cat(sprintf("Saved: %s\n", basename(out_full)))
  
  # ============ 2. Top 5 per microenv ============
  cat("\n[2] Top 5 per microenv\n")
  top5 <- me_ct %>%
    group_by(microenv_id) %>%
    slice_max(mean_weight, n = 5, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(microenv_id, desc(mean_weight))
  out_top5 <- file.path(OUT_DIR, paste0(organ, "_microenv_celltypes_top5.csv"))
  write.csv(top5, out_top5, row.names = FALSE)
  cat(sprintf("Saved: %s\n", basename(out_top5)))
  
  cat("\n--- Top 5 cell types per microenv ---\n")
  for (m in unique(top5$microenv_id)) {
    cat(sprintf("\n%s (%s):\n", m, me_to_label[m]))
    sub <- top5[top5$microenv_id == m, ]
    for (i in seq_len(nrow(sub))) {
      cat(sprintf("  %-25s %.3f (%.1f%%)\n",
                  sub$celltype[i], sub$mean_weight[i], sub$mean_weight[i] * 100))
    }
  }
  
  # ============ 3. Specificity (log2FC vs global) ============
  cat("\n[3] Specificity (log2FC vs global mean weight)\n")
  global_w <- colMeans(weight_mat)
  
  spec_df <- me_ct %>%
    rowwise() %>%
    mutate(global_weight = global_w[celltype],
           log2FC = log2((mean_weight + 1e-5) / (global_weight + 1e-5)),
           enriched = log2FC > 1) %>%
    ungroup() %>%
    arrange(microenv_id, desc(log2FC))
  
  out_spec <- file.path(OUT_DIR, paste0(organ, "_microenv_specificity.csv"))
  write.csv(spec_df, out_spec, row.names = FALSE)
  cat(sprintf("Saved: %s\n", basename(out_spec)))
  
  # ============ 4. Tumor vs baseline shift ============
  cat("\n[4] Tumor vs baseline weight shift per microenv\n")
  
  # 看下 group_norm 的实际值
  cat(sprintf("group_norm values: %s\n",
              paste(unique(long_df$group_norm), collapse = ",")))
  
  tvb_long <- long_df %>%
    select(microenv_id, group_norm, all_of(all_ct)) %>%
    pivot_longer(-c(microenv_id, group_norm), names_to = "celltype", values_to = "weight") %>%
    group_by(microenv_id, group_norm, celltype) %>%
    summarise(mean_weight = mean(weight), .groups = "drop")
  
  cat(sprintf("tvb_long rows: %d, group_norm distinct: %s\n",
              nrow(tvb_long),
              paste(unique(tvb_long$group_norm), collapse = ",")))
  
  tvb_wide <- tvb_long %>%
    pivot_wider(names_from = group_norm, values_from = mean_weight, values_fill = 0)
  
  cat(sprintf("tvb_wide cols: %s\n", paste(colnames(tvb_wide), collapse = ",")))
  
  # 安全检查: 必须有 tumor 和 baseline 两列才能算 log2FC
  if (all(c("tumor", "baseline") %in% colnames(tvb_wide))) {
    tvb_wide <- tvb_wide %>%
      mutate(log2FC_tumor_vs_baseline = log2((tumor + 1e-5) / (baseline + 1e-5)),
             microenv_label = me_to_label[microenv_id]) %>%
      arrange(microenv_id, desc(log2FC_tumor_vs_baseline))
    
    out_tvb <- file.path(OUT_DIR, paste0(organ, "_microenv_tumor_vs_baseline.csv"))
    write.csv(tvb_wide, out_tvb, row.names = FALSE)
    cat(sprintf("Saved: %s\n", basename(out_tvb)))
    
    # 打印 paper-key M05 的 tumor vs baseline (liver 才有)
    if (organ == "liver" && "M05" %in% tvb_wide$microenv_id) {
      cat("\n--- M05 paper-key niche: top 10 tumor-enriched cell types ---\n")
      m05_top <- tvb_wide %>%
        filter(microenv_id == "M05") %>%
        slice_max(log2FC_tumor_vs_baseline, n = 10, with_ties = FALSE)
      print(m05_top[, c("microenv_id","celltype","baseline","tumor","log2FC_tumor_vs_baseline")])
    }
  } else {
    cat(sprintf("⚠ Cannot compute tumor vs baseline: missing column.\n"))
    cat(sprintf("  Available: %s\n", paste(colnames(tvb_wide), collapse=",")))
    
    # 退而求其次: 写出 wide 表
    out_tvb <- file.path(OUT_DIR, paste0(organ, "_microenv_group_weight.csv"))
    write.csv(tvb_wide, out_tvb, row.names = FALSE)
    cat(sprintf("Saved (without log2FC): %s\n", basename(out_tvb)))
  }
  
  # ============ 5. Heatmap (microenv × cell type) ============
  cat("\n[5] Heatmap\n")
  heat_mat <- me_ct %>%
    pivot_wider(id_cols = microenv_id, names_from = celltype, values_from = mean_weight) %>%
    tibble::column_to_rownames("microenv_id") %>% as.matrix()
  
  row_order <- sort(rownames(heat_mat))
  heat_mat <- heat_mat[row_order, , drop = FALSE]
  
  col_max <- apply(heat_mat, 2, max)
  top_cts <- names(sort(col_max, decreasing = TRUE))[1:min(25, length(col_max))]
  heat_top <- heat_mat[, top_cts, drop = FALSE]
  
  rownames(heat_top) <- paste0(rownames(heat_top), " ",
                                me_to_label[rownames(heat_top)])
  
  .hm_h <- max(4, nrow(heat_top) * 0.45 + 2)
  .ph <- pheatmap(heat_top,
           color = colorRampPalette(c("white","yellow","red","darkred"))(50),
           cluster_rows = FALSE, cluster_cols = FALSE,
           main = NA,                                   # [restyle] 标题已移除
           silent = TRUE,                               # [restyle] 不直接出图,手动控分辨率
           fontsize_row = 12, fontsize_col = 11,        # [restyle] 9/8 -> 12/11
           display_numbers = TRUE, number_format = "%.2f", number_color = "black",
           fontsize_number = 9)                         # [restyle] 6 -> 9
  save_fig_all(.ph$gtable,                                   # [restyle] PDF+TIFF+PNG@300
               file.path(RESTYLE_DIR, sprintf("%s_microenv_celltype_heatmap", organ)),
               width = 12, height = .hm_h)
  cat(sprintf("Saved: %s_microenv_celltype_heatmap.png/.pdf\n", organ))
  
  # ============ 6. Dotplot (top 5 per microenv) ============
  cat("\n[6] Dotplot\n")
  top5_for_dot <- top5 %>% mutate(microenv_full = paste0(microenv_id, " ", microenv_label))
  
  p_dot <- ggplot(top5_for_dot,
                  aes(x = celltype, y = microenv_full,
                      size = n_spots, color = mean_weight)) +
    geom_point() +
    scale_color_gradient(low = "lightyellow", high = "red", name = "Mean weight") +
    scale_size_continuous(range = c(3, 12), name = "n spots") +
    labs(x = "Cell type", y = "Microenvironment") +   # [restyle] 标题已移除
    theme_classic(base_size = 14) +                    # [restyle] 11->14
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # [restyle]
          axis.text.y = element_text(size = 12),
          legend.title = element_text(size = 13),
          legend.text  = element_text(size = 11))

  .dot_w <- max(8, length(unique(top5_for_dot$celltype)) * 0.3 + 4)
  .dot_h <- max(5, length(unique(top5_for_dot$microenv_full)) * 0.5 + 2)
  save_fig_all(p_dot,                                        # [restyle] PDF+TIFF+PNG@300
               file.path(RESTYLE_DIR, sprintf("%s_microenv_celltype_dotplot", organ)),
               width = .dot_w, height = .dot_h)
  cat(sprintf("Saved: %s_microenv_celltype_dotplot.pdf\n", organ))
  
  # ============ Return ============
  list(organ = organ,
       n_microenv = length(unique(spot_assign$microenv_id)),
       n_spots = nrow(spot_assign),
       top5 = top5)
}

# ============================================================
results <- list()
results[["liver"]] <- process_organ("liver")
results[["brain"]] <- process_organ("brain")

# Summary
md_file <- file.path(OUT_DIR, "summary.md")
con <- file(md_file, "w")
writeLines(c(
  "# Microenvironment Key Cell Types Summary",
  "",
  paste0("**Generated:** ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## Method",
  "- For each microenvironment, compute mean RCTD weight per cell type",
  "- Top 5 cell types per microenv (by mean weight)",
  "- Specificity: log2FC vs global mean weight",
  "- Tumor vs baseline weight shift per (microenv x cell type)",
  "",
  "## Results"
), con)

for (organ in c("liver", "brain")) {
  r <- results[[organ]]
  if (is.null(r)) next
  writeLines(c("",
               paste0("### ", toupper(organ),
                      " (", r$n_microenv, " microenvs, ", r$n_spots, " spots)"),
               "", "**Top 5 cell types per microenv:**", ""), con)
  for (m in unique(r$top5$microenv_id)) {
    sub <- r$top5[r$top5$microenv_id == m, ]
    label <- sub$microenv_label[1]
    cts <- paste(sprintf("%s (%.1f%%)", sub$celltype, sub$mean_weight * 100),
                  collapse = ", ")
    writeLines(paste0("- **", m, " ", label, "**: ", cts), con)
  }
}

writeLines(c("",
             "## Files",
             "- `<organ>_microenv_celltypes_full.csv` - microenv x ct x mean_weight x n_spots",
             "- `<organ>_microenv_celltypes_top5.csv` - top 5 per microenv",
             "- `<organ>_microenv_specificity.csv` - log2FC vs global",
             "- `<organ>_microenv_tumor_vs_baseline.csv` - tumor vs baseline shift",
             "- `<organ>_microenv_celltype_heatmap.pdf`",
             "- `<organ>_microenv_celltype_dotplot.pdf`"), con)
close(con)
cat("\nSaved:", md_file, "\n")

cat("\n=============================================================\n")
cat("Done. Output:", OUT_DIR, "\n")
cat("=============================================================\n")
system(sprintf("ls -la %s", OUT_DIR))
