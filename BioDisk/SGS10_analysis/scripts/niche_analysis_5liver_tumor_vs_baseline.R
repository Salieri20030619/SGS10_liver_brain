#!/usr/bin/env Rscript
# ============================================================
# niche_analysis_5liver_tumor_vs_baseline.R
# Niche 分析 (正式版): Tumor vs Baseline 对比
# ============================================================
#
# 样本:
#   Tumor   (3): Liver_3, Liver_8, gan56
#   Blank   (1): Bla_liver_1
#   Control (1): Con_liver_1
#
# 目的: 识别 *tumor-induced* niche 变化 (不是生理 niche)
#
# 核心方法:
#   1. 每个样本独立做 spatial neighborhood k=6 NN + permutation
#   2. 对每个 cell type pair, 计算在 tumor 组和 baseline 组的 log2 enrichment
#   3. 差值 Δlog2 = log2(tumor) - log2(baseline) = tumor-induced niche change
#
# 关键输出:
#   "delta_enrichment_matrix.csv" ← tumor vs baseline 的 niche 变化
#   Top tumor-induced pairs (Δ > 1.5, baseline 近 0)
#   Top tumor-lost pairs (Δ < -1.5)
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(Matrix)
  library(Seurat)
  library(FNN)
  library(pheatmap)
  library(RColorBrewer)
})

PROJECT_ROOT <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
setwd(PROJECT_ROOT)

OUT_DIR <- file.path(PROJECT_ROOT, "taskB", "niche_analysis_v1")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Sample groups
TUMOR_SAMPLES    <- c("Liver_3", "Liver_8", "gan56")
BASELINE_SAMPLES <- c("Bla_liver_1", "Con_liver_1")
ALL_SAMPLES      <- c(TUMOR_SAMPLES, BASELINE_SAMPLES)

K_NN <- 6
PERMUTATION_N <- 100

cat("============================================\n")
cat("Niche 分析 v1: Tumor vs Baseline 对比\n")
cat(sprintf("  Tumor:    %s\n", paste(TUMOR_SAMPLES, collapse=", ")))
cat(sprintf("  Baseline: %s\n", paste(BASELINE_SAMPLES, collapse=", ")))
cat("============================================\n\n")

# ============================================================
# 1. 加载 5 样本数据
# ============================================================
cat("### 1. 加载数据 ###\n")

all_data <- list()
for (s in ALL_SAMPLES) {
  cat(sprintf("[%s] ", s))
  
  dom_file <- file.path("spatial_results/rctd/per_sample",
                         sprintf("%s_dominant_celltype.csv", s))
  if (!file.exists(dom_file)) {
    cat("✗ RCTD not found, skip\n"); next
  }
  dom <- read.csv(dom_file, stringsAsFactors = FALSE)
  
  obj <- readRDS(sprintf("spatial_results/%s_seekspace.rds", s))
  meta <- obj@meta.data[, c("zone"), drop = FALSE]
  meta$barcode <- rownames(meta)
  coords <- as.data.frame(obj@reductions$spatial@cell.embeddings)
  coords$barcode <- rownames(coords)
  colnames(coords)[1:2] <- c("x", "y")
  rm(obj); gc(verbose = FALSE)
  
  df <- dom %>%
    left_join(meta, by = "barcode") %>%
    left_join(coords, by = "barcode") %>%
    mutate(sample = s,
           group  = ifelse(s %in% TUMOR_SAMPLES, "tumor", "baseline"))
  
  all_data[[s]] <- df
  cat(sprintf("%d spots, group=%s\n", nrow(df), df$group[1]))
}

cat(sprintf("\nLoaded %d samples, total spots: %d\n\n",
            length(all_data), sum(sapply(all_data, nrow))))

# ============================================================
# 2. 每样本计算 spatial neighborhood enrichment
# ============================================================
cat("### 2. 单样本 spatial NN enrichment ###\n")

per_sample_enrich <- list()

for (s in names(all_data)) {
  df <- all_data[[s]] %>% filter(!is.na(x), !is.na(y), !is.na(dominant_celltype))
  if (nrow(df) < 100) next
  
  cat(sprintf("[%s] k=%d NN, perm N=%d... ", s, K_NN, PERMUTATION_N))
  
  coords_mat <- as.matrix(df[, c("x", "y")])
  nn_result <- get.knn(coords_mat, k = K_NN)
  ct_vec <- df$dominant_celltype
  
  # Observed
  nn_edges <- data.frame(
    center_ct   = rep(ct_vec, each = K_NN),
    neighbor_ct = ct_vec[as.vector(t(nn_result$nn.index))]
  )
  obs_count <- nn_edges %>% count(center_ct, neighbor_ct, name = "obs")
  
  # Permutation expected
  rand_accum <- NULL
  for (p in seq_len(PERMUTATION_N)) {
    ct_shuf <- sample(ct_vec)
    rand_edges <- data.frame(
      center_ct   = rep(ct_shuf, each = K_NN),
      neighbor_ct = ct_shuf[as.vector(t(nn_result$nn.index))]
    )
    rc <- rand_edges %>% count(center_ct, neighbor_ct, name = "rand")
    if (is.null(rand_accum)) {
      rand_accum <- rc %>% mutate(rand_sum = rand) %>% select(-rand)
    } else {
      rand_accum <- full_join(rand_accum, rc,
                                by = c("center_ct", "neighbor_ct")) %>%
        mutate(rand = ifelse(is.na(rand), 0, rand),
               rand_sum = rand_sum + rand) %>%
        select(-rand)
    }
  }
  rand_accum <- rand_accum %>% mutate(rand_mean = rand_sum / PERMUTATION_N)
  
  enrich_sample <- obs_count %>%
    full_join(rand_accum, by = c("center_ct", "neighbor_ct")) %>%
    mutate(
      sample = s,
      group = df$group[1],
      obs = ifelse(is.na(obs), 0, obs),
      rand_mean = ifelse(is.na(rand_mean) | rand_mean < 1, 1, rand_mean),
      log2_enrichment = log2((obs + 1) / (rand_mean + 1)),  # pseudocount 防 -Inf
      z_score = (obs - rand_mean) / sqrt(rand_mean + 1)
    )
  
  per_sample_enrich[[s]] <- enrich_sample
  cat(sprintf("✓ %d pairs\n", nrow(enrich_sample)))
}

# ============================================================
# 3. 组合统计: Tumor vs Baseline 对比
# ============================================================
cat("\n### 3. Tumor vs Baseline 对比 ###\n")

all_enrich <- bind_rows(per_sample_enrich)

# Per-pair, per-group 平均 log2 enrichment
group_stats <- all_enrich %>%
  group_by(center_ct, neighbor_ct, group) %>%
  summarise(
    n_samples    = n(),
    mean_log2    = mean(log2_enrichment, na.rm = TRUE),
    mean_obs     = mean(obs, na.rm = TRUE),
    mean_rand    = mean(rand_mean, na.rm = TRUE),
    .groups = "drop"
  )

# Pivot 成 tumor / baseline 两列
delta_df <- group_stats %>%
  select(center_ct, neighbor_ct, group, mean_log2) %>%
  pivot_wider(names_from = group, values_from = mean_log2,
              names_prefix = "log2_") %>%
  mutate(
    log2_tumor    = ifelse(is.na(log2_tumor), 0, log2_tumor),
    log2_baseline = ifelse(is.na(log2_baseline), 0, log2_baseline),
    delta_log2    = log2_tumor - log2_baseline,
    # tumor_induced: 肿瘤里富集, baseline 里不富集 (或负)
    tumor_induced = (log2_tumor > 0.5) & (log2_baseline < 0.3) & (delta_log2 > 0.5),
    tumor_lost    = (log2_baseline > 0.5) & (log2_tumor < 0.3) & (delta_log2 < -0.5)
  )

# 加 obs count (两组合计, 用 tumor 的)
obs_stats <- group_stats %>%
  select(center_ct, neighbor_ct, group, mean_obs) %>%
  pivot_wider(names_from = group, values_from = mean_obs,
              names_prefix = "obs_")

delta_df <- delta_df %>%
  left_join(obs_stats, by = c("center_ct", "neighbor_ct"))

delta_df <- delta_df %>%
  arrange(desc(delta_log2))

cat(sprintf("Unique cell type pairs: %d\n", nrow(delta_df)))
cat(sprintf("  tumor_induced (Δ>0.5, baseline<0.3, tumor>0.5): %d\n",
            sum(delta_df$tumor_induced, na.rm = TRUE)))
cat(sprintf("  tumor_lost (Δ<-0.5, baseline>0.5, tumor<0.3): %d\n",
            sum(delta_df$tumor_lost, na.rm = TRUE)))

write.csv(delta_df,
          file.path(OUT_DIR, "02_delta_enrichment_tumor_vs_baseline.csv"),
          row.names = FALSE)

# ============================================================
# 4. Top tumor-induced niche pairs (最有意义的!)
# ============================================================
cat("\n### 4. Top TUMOR-INDUCED niche pairs (⭐核心) ###\n")

tumor_induced_top <- delta_df %>%
  filter(tumor_induced, center_ct != neighbor_ct, obs_tumor >= 3) %>%
  arrange(desc(delta_log2))

cat(sprintf("\nTop 20 tumor-induced niche (tumor ↑, baseline 平):\n"))
cat("| Center → Neighbor | log2 Tumor | log2 Baseline | Δ log2 | mean obs (tumor) |\n")
cat("|-------------------|-----------:|--------------:|-------:|------------------:|\n")
for (i in seq_len(min(20, nrow(tumor_induced_top)))) {
  r <- tumor_induced_top[i, ]
  cat(sprintf("| %s → %s | %+.2f | %+.2f | %+.2f | %.1f |\n",
              r$center_ct, r$neighbor_ct, r$log2_tumor, r$log2_baseline,
              r$delta_log2, r$obs_tumor))
}

write.csv(tumor_induced_top,
          file.path(OUT_DIR, "03_top_tumor_induced_niche.csv"),
          row.names = FALSE)

# ============================================================
# 5. Top tumor-LOST niche pairs (正常生理 niche 丢失)
# ============================================================
cat("\n### 5. Top TUMOR-LOST niche pairs (正常 niche 消失) ###\n")

tumor_lost_top <- delta_df %>%
  filter(tumor_lost, center_ct != neighbor_ct, obs_baseline >= 3) %>%
  arrange(delta_log2)

cat(sprintf("\nTop 15 tumor-lost niche (baseline 有, tumor 无):\n"))
cat("| Center → Neighbor | log2 Tumor | log2 Baseline | Δ log2 |\n")
cat("|-------------------|-----------:|--------------:|-------:|\n")
for (i in seq_len(min(15, nrow(tumor_lost_top)))) {
  r <- tumor_lost_top[i, ]
  cat(sprintf("| %s → %s | %+.2f | %+.2f | %+.2f |\n",
              r$center_ct, r$neighbor_ct, r$log2_tumor, r$log2_baseline,
              r$delta_log2))
}

write.csv(tumor_lost_top,
          file.path(OUT_DIR, "04_top_tumor_lost_niche.csv"),
          row.names = FALSE)

# ============================================================
# 6. Focus triad: HSC_Activated_ECM / Kupffer (Task A L1/L4)
# ============================================================
cat("\n### 6. Focus triad: Task A 论点的 cell type ###\n")

for (focus_ct in c("HSC_Activated_ECM", "Kupffer")) {
  cat(sprintf("\n--- %s niche (tumor vs baseline) ---\n", focus_ct))
  focus_df <- delta_df %>%
    filter(center_ct == focus_ct, neighbor_ct != focus_ct,
            (obs_tumor >= 3 | obs_baseline >= 3)) %>%
    arrange(desc(delta_log2))
  
  cat(sprintf("Top 10 %s neighbors (by delta log2):\n", focus_ct))
  head(focus_df, 10) %>% print(n = 10)
  
  write.csv(focus_df,
            file.path(OUT_DIR, sprintf("05_%s_niche_tumor_vs_baseline.csv", focus_ct)),
            row.names = FALSE)
}

# ============================================================
# 7. Immune × Endothelial (triangle 肝端)
# ============================================================
cat("\n### 7. Immune × Endothelial (triangle 肝端) ###\n")

immune_cells <- c("Kupffer", "KC_Homeostatic", "Macrophage_M1", "Macrophage_M2",
                   "Macrophage_NAM", "Macrophage_MDM", "cDC2", "pDC",
                   "CD4_T", "CD8_T", "NK_cells", "B_Plasma")
endo_cells <- c("Endothelial", "LSEC", "LSEC_Mid", "LSEC_Periportal",
                 "LSEC_Pericentral", "LVEC")

immune_endo_delta <- delta_df %>%
  filter(center_ct %in% immune_cells, neighbor_ct %in% endo_cells) %>%
  arrange(desc(delta_log2))

cat(sprintf("\nTop 15 immune-endothelial niche changes:\n"))
head(immune_endo_delta, 15) %>% print(n = 15)

write.csv(immune_endo_delta,
          file.path(OUT_DIR, "06_immune_endothelial_delta.csv"),
          row.names = FALSE)

# ============================================================
# 8. 可视化
# ============================================================
cat("\n### 8. 可视化 ###\n")

# 取 top 20 cell type
combined_for_top <- bind_rows(all_data) %>%
  filter(!is.na(dominant_celltype)) %>%
  count(dominant_celltype, sort = TRUE)
top20_cts <- head(combined_for_top$dominant_celltype, 20)

# Delta heatmap (tumor-induced niche 核心可视化)
delta_mat <- delta_df %>%
  filter(center_ct %in% top20_cts, neighbor_ct %in% top20_cts) %>%
  select(center_ct, neighbor_ct, delta_log2) %>%
  pivot_wider(names_from = neighbor_ct, values_from = delta_log2,
              values_fill = 0)
rn <- delta_mat$center_ct
delta_mat <- as.matrix(delta_mat[, -1])
rownames(delta_mat) <- rn
delta_mat[is.infinite(delta_mat)] <- NA
delta_mat[delta_mat > 3] <- 3
delta_mat[delta_mat < -3] <- -3

pdf(file.path(OUT_DIR, "07_delta_enrichment_heatmap.pdf"), width = 11, height = 10)
pheatmap(delta_mat,
         main = "Δ log2 enrichment (tumor - baseline): tumor-induced niche",
         color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
         breaks = seq(-3, 3, length.out = 101),
         fontsize = 8,
         cluster_rows = TRUE, cluster_cols = TRUE,
         display_numbers = TRUE, number_format = "%.1f",
         fontsize_number = 5.5,
         border_color = "gray80")
dev.off()
cat("  → 07_delta_enrichment_heatmap.pdf\n")

# Bar chart: Top 20 tumor-induced pairs
plot_df <- tumor_induced_top %>%
  head(20) %>%
  mutate(
    pair = paste(center_ct, "→", neighbor_ct),
    pair = factor(pair, levels = rev(pair))
  )

p <- ggplot(plot_df) +
  geom_col(aes(x = delta_log2, y = pair, fill = log2_tumor), width = 0.7) +
  scale_fill_gradient(low = "#FEE0D2", high = "#A50F15",
                       name = "log2 in\ntumor") +
  geom_text(aes(x = delta_log2, y = pair,
                 label = sprintf("tumor=%+.1f, base=%+.1f", log2_tumor, log2_baseline)),
             hjust = -0.05, size = 2.5) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "gray60") +
  xlim(0, max(plot_df$delta_log2, na.rm = TRUE) * 1.8) +
  theme_minimal(base_size = 10) +
  labs(title = "Top 20 tumor-induced niche pairs",
        subtitle = "tumor-baseline Δ log2 > 0.5, 基于 3 tumor + 2 baseline liver samples",
        x = "Δ log2 (tumor - baseline)", y = "")

ggsave(file.path(OUT_DIR, "08_top_tumor_induced_barplot.pdf"),
       p, width = 11, height = 7, dpi = 150)
cat("  → 08_top_tumor_induced_barplot.pdf\n")

# Focus bar charts (HSC_ECM, Kupffer)
for (focus_ct in c("HSC_Activated_ECM", "Kupffer")) {
  f_df <- delta_df %>%
    filter(center_ct == focus_ct, neighbor_ct != focus_ct,
            (obs_tumor >= 3 | obs_baseline >= 3)) %>%
    arrange(desc(delta_log2)) %>%
    head(12) %>%
    mutate(neighbor_ct = factor(neighbor_ct, levels = rev(neighbor_ct)))
  
  if (nrow(f_df) == 0) next
  
  f_long <- f_df %>%
    select(neighbor_ct, log2_tumor, log2_baseline) %>%
    pivot_longer(c(log2_tumor, log2_baseline),
                  names_to = "group", values_to = "log2_enrich") %>%
    mutate(group = gsub("log2_", "", group))
  
  p <- ggplot(f_long, aes(x = log2_enrich, y = neighbor_ct, fill = group)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c(tumor = "#E41A1C", baseline = "#377EB8"),
                       name = "Group") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme_minimal(base_size = 10) +
    labs(title = sprintf("%s 的 niche 邻居: tumor vs baseline", focus_ct),
          subtitle = "Top 12 by delta log2 enrichment",
          x = "log2 enrichment", y = "Neighbor cell type")
  
  ggsave(file.path(OUT_DIR, sprintf("09_%s_niche_compare.pdf", focus_ct)),
         p, width = 10, height = 6, dpi = 150)
}
cat("  → 09_HSC_Activated_ECM_niche_compare.pdf\n")
cat("  → 09_Kupffer_niche_compare.pdf\n")

# ============================================================
# 9. Summary markdown
# ============================================================
cat("\n### 9. 生成 summary.md ###\n")

lines <- c(
  "# Liver Niche 分析 v1: Tumor vs Baseline 对比",
  "",
  sprintf("生成时间: %s", format(Sys.time(), "%Y-%m-%d %H:%M")),
  "",
  "## 样本",
  "",
  sprintf("- **Tumor** (n=3): %s", paste(TUMOR_SAMPLES, collapse=", ")),
  sprintf("- **Baseline** (n=2): %s (%d blank + %d control)",
          paste(BASELINE_SAMPLES, collapse=", "), 1, 1),
  "",
  "## 方法",
  "",
  sprintf("1. 每样本 spatial neighborhood k=%d NN + permutation N=%d", K_NN, PERMUTATION_N),
  "2. log2 enrichment = log2((obs+1) / (rand+1))",
  "3. **Tumor-induced niche**: log2_tumor > 0.5 AND log2_baseline < 0.3 AND Δ > 0.5",
  "4. **Tumor-lost niche**: log2_baseline > 0.5 AND log2_tumor < 0.3 AND Δ < -0.5",
  "",
  "## 核心发现",
  "",
  sprintf("- Tumor-induced niche pairs: **%d**", sum(delta_df$tumor_induced, na.rm = TRUE)),
  sprintf("- Tumor-lost niche pairs: %d", sum(delta_df$tumor_lost, na.rm = TRUE)),
  "",
  "### Top 15 Tumor-Induced Niche (⭐ 核心)",
  "",
  "这些是 **tumor 特异出现** 的空间共定位 pattern, 正常肝中不存在:",
  "",
  "| Center | Neighbor | log2 tumor | log2 baseline | Δ log2 |",
  "|--------|----------|-----------:|--------------:|-------:|"
)

for (i in seq_len(min(15, nrow(tumor_induced_top)))) {
  r <- tumor_induced_top[i, ]
  lines <- c(lines, sprintf("| %s | %s | %+.2f | %+.2f | **%+.2f** |",
                              r$center_ct, r$neighbor_ct,
                              r$log2_tumor, r$log2_baseline, r$delta_log2))
}

lines <- c(lines, "",
            "### Task A 论点 niche 验证",
            "",
            "#### HSC_Activated_ECM niche (L1 论点)",
            "",
            "Top 10 邻居 (按 Δ log2):",
            "",
            "| Neighbor | log2 tumor | log2 baseline | Δ log2 |",
            "|----------|-----------:|--------------:|-------:|")

hsc_df <- delta_df %>%
  filter(center_ct == "HSC_Activated_ECM", neighbor_ct != "HSC_Activated_ECM",
          (obs_tumor >= 3 | obs_baseline >= 3)) %>%
  arrange(desc(delta_log2)) %>%
  head(10)
for (i in seq_len(nrow(hsc_df))) {
  r <- hsc_df[i, ]
  lines <- c(lines, sprintf("| %s | %+.2f | %+.2f | %+.2f |",
                              r$neighbor_ct, r$log2_tumor,
                              r$log2_baseline, r$delta_log2))
}

lines <- c(lines, "",
            "#### Kupffer niche (L4 论点)",
            "",
            "| Neighbor | log2 tumor | log2 baseline | Δ log2 |",
            "|----------|-----------:|--------------:|-------:|")

kupf_df <- delta_df %>%
  filter(center_ct == "Kupffer", neighbor_ct != "Kupffer",
          (obs_tumor >= 3 | obs_baseline >= 3)) %>%
  arrange(desc(delta_log2)) %>%
  head(10)
for (i in seq_len(nrow(kupf_df))) {
  r <- kupf_df[i, ]
  lines <- c(lines, sprintf("| %s | %+.2f | %+.2f | %+.2f |",
                              r$neighbor_ct, r$log2_tumor,
                              r$log2_baseline, r$delta_log2))
}

lines <- c(lines, "",
            "## 输出",
            "",
            "- `02_delta_enrichment_tumor_vs_baseline.csv` — 所有 ct pair 的 Δ log2",
            "- `03_top_tumor_induced_niche.csv` — ⭐ 核心: tumor-induced pair",
            "- `04_top_tumor_lost_niche.csv` — tumor-lost pair",
            "- `05_HSC_Activated_ECM_niche_tumor_vs_baseline.csv`",
            "- `05_Kupffer_niche_tumor_vs_baseline.csv`",
            "- `06_immune_endothelial_delta.csv`",
            "- `07_delta_enrichment_heatmap.pdf` — 20×20 热图",
            "- `08_top_tumor_induced_barplot.pdf` — Top 20 bar 图",
            "- `09_HSC_Activated_ECM_niche_compare.pdf`, `09_Kupffer_niche_compare.pdf`",
            "",
            "## 下一步",
            "",
            "- 等 brain 样本全完成后, 跑 brain niche (真 immune-endo-neuron triangle)",
            "- Spatial CellChat 基于 tumor-induced niche, 识别 L-R 通讯 pair",
            "- 与 Task A L1/L4 论点整合, 加入 Task B v1 报告"
)

writeLines(lines, file.path(OUT_DIR, "summary.md"))
cat(sprintf("  → summary.md (%d lines)\n", length(lines)))

cat("\n============================================\n")
cat("Liver niche 分析 v1 完成!\n")
cat(sprintf("  输出: %s\n", OUT_DIR))
cat(sprintf("  文件数: %d\n", length(list.files(OUT_DIR))))
cat("============================================\n")
