#!/usr/bin/env Rscript
# ============================================================
# zbtb7a_liver_vs_brain.R
# Task 1.2: Liver Zbtb7a vs Brain Zbtb7a 对比分析
# ============================================================
# 
# 背景:
#   Task A v1.1 发现 Zbtb7a 在 12/22 brain cell type 活性一致下调 (pan-brain signature).
#   Liu 2014 Cell Metab 报告 Zbtb7a 参与肝代谢 (Pck1 / gluconeogenesis).
#   由于 Thyroid hormone 也调控 Pck1, Zbtb7a 可能是一个**新的肝脑轴候选论点**.
#
# 本分析:
#   对比 liver 8 cell types 和 brain 22 cell types 的 Zbtb7a activity,
#   验证 Zbtb7a 是否在 liver 也呈广泛下调 (pan-organ signature).
#   如果成立, 写入 Task A v1.2 作为新的肝脑轴论点候选.
#
# 输入:
#   dorothea_results/per_contrast/liver__tumor_vs_control__*_tf_activity.csv  (8 个)
#   dorothea_results_brain_clean/per_contrast/brain__tumor_vs_control__*_tf_activity.csv  (22 个)
#
# 输出:
#   taskA_v1_2/zbtb7a_analysis/
#     05_zbtb7a_liver_vs_brain.csv   — 所有 cell type 的 Zbtb7a activity
#     06_zbtb7a_organ_comparison.md  — 对比报告
#     07_zbtb7a_organ_scatter.pdf    — 可视化
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

PROJECT_ROOT <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
setwd(PROJECT_ROOT)

OUT_DIR <- file.path(PROJECT_ROOT, "taskA_v1_2", "zbtb7a_analysis")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("============================================\n")
cat("Zbtb7a: Liver vs Brain 对比分析\n")
cat("============================================\n\n")

# ============================================================
# 1. 加载 Liver 和 Brain 的 Zbtb7a activity
# ============================================================

# Helper: 从 per_contrast 文件里拿 Zbtb7a 的 score + padj
extract_zbtb7a <- function(file) {
  if (!file.exists(file)) {
    return(list(score = NA, p_value = NA, padj = NA, sig = NA, 
                sig_strict = NA, sig_relaxed = NA))
  }
  df <- read.csv(file, stringsAsFactors = FALSE)
  row <- df[df$TF == "Zbtb7a", ]
  if (nrow(row) == 0) {
    return(list(score = NA, p_value = NA, padj = NA, sig = NA,
                sig_strict = NA, sig_relaxed = NA))
  }
  list(
    score       = row$score[1],
    p_value     = row$p_value[1],
    padj        = row$padj[1],
    sig         = row$sig[1],
    sig_strict  = row$sig_strict[1],
    sig_relaxed = row$sig_relaxed[1]
  )
}

# Liver contrasts (dorothea_results/, 8 cell types × 2 contrasts)
cat("### 加载 Liver Zbtb7a activity ###\n")
LIVER_DIR <- "dorothea_results/per_contrast"
liver_files <- list.files(LIVER_DIR, pattern = "^liver__tumor_vs.*_tf_activity\\.csv$",
                          full.names = TRUE)
cat(sprintf("Liver contrast files: %d\n", length(liver_files)))

liver_df <- data.frame()
for (f in liver_files) {
  fname <- basename(f)
  # parse: liver__tumor_vs_control__HSC_Activated_ECM_tf_activity.csv
  parts <- strsplit(gsub("_tf_activity\\.csv$", "", fname), "__")[[1]]
  # parts = c("liver", "tumor_vs_blank" or "tumor_vs_control", "<celltype>")
  organ    <- parts[1]
  contrast <- parts[2]
  celltype <- parts[3]
  
  z <- extract_zbtb7a(f)
  liver_df <- bind_rows(liver_df, data.frame(
    organ = organ, contrast = contrast, celltype = celltype,
    score = z$score, p_value = z$p_value, padj = z$padj,
    sig_strict = z$sig_strict, sig_relaxed = z$sig_relaxed,
    stringsAsFactors = FALSE
  ))
}
cat(sprintf("  Liver rows: %d\n", nrow(liver_df)))

# Brain clean contrasts (dorothea_results_brain_clean/, 22 cell types × 2 contrasts)
cat("\n### 加载 Brain Zbtb7a activity ###\n")
BRAIN_DIR <- "dorothea_results_brain_clean/per_contrast"
brain_files <- list.files(BRAIN_DIR, pattern = "^brain__tumor_vs.*_tf_activity\\.csv$",
                          full.names = TRUE)
cat(sprintf("Brain contrast files: %d\n", length(brain_files)))

brain_df <- data.frame()
for (f in brain_files) {
  fname <- basename(f)
  parts <- strsplit(gsub("_tf_activity\\.csv$", "", fname), "__")[[1]]
  organ    <- parts[1]
  contrast <- parts[2]
  celltype <- parts[3]
  
  z <- extract_zbtb7a(f)
  brain_df <- bind_rows(brain_df, data.frame(
    organ = organ, contrast = contrast, celltype = celltype,
    score = z$score, p_value = z$p_value, padj = z$padj,
    sig_strict = z$sig_strict, sig_relaxed = z$sig_relaxed,
    stringsAsFactors = FALSE
  ))
}
cat(sprintf("  Brain rows: %d\n", nrow(brain_df)))

# ============================================================
# 2. 合并
# ============================================================
combined <- bind_rows(liver_df, brain_df) %>%
  arrange(organ, contrast, score)

write.csv(combined,
          file.path(OUT_DIR, "05_zbtb7a_liver_vs_brain.csv"),
          row.names = FALSE)
cat(sprintf("\n✓ Saved: 05_zbtb7a_liver_vs_brain.csv (%d rows)\n", nrow(combined)))

# ============================================================
# 3. 关键对比: Tumor_vs_Control contrast 下 Zbtb7a 的 organ pattern
# ============================================================
cat("\n### 3. Tumor_vs_Control 下 Zbtb7a activity ###\n")

tvc <- combined %>% filter(contrast == "tumor_vs_control")
cat("\n--- Liver (Tumor_vs_Control) ---\n")
tvc_liver <- tvc %>% filter(organ == "liver") %>% arrange(score)
print(tvc_liver, row.names = FALSE)

cat("\n--- Brain (Tumor_vs_Control) ---\n")
tvc_brain <- tvc %>% filter(organ == "brain") %>% arrange(score)
print(tvc_brain, row.names = FALSE)

# ============================================================
# 4. 统计 summary
# ============================================================
cat("\n### 4. 统计 summary ###\n")

summary_stats <- tvc %>%
  filter(!is.na(score)) %>%
  group_by(organ) %>%
  summarise(
    n_celltypes       = n(),
    n_down            = sum(score < 0),
    n_down_strict     = sum(sig_strict == "Down", na.rm = TRUE),
    n_down_relaxed    = sum(sig_relaxed == "Down", na.rm = TRUE),
    pct_down          = round(100 * n_down / n_celltypes, 1),
    median_score      = round(median(score), 2),
    min_score         = round(min(score), 2),
    max_score         = round(max(score), 2),
    .groups = "drop"
  )

cat("\n跨器官 summary (Tumor_vs_Control):\n")
print(summary_stats)

# 关键判断: 两个器官都呈 consistent down?
liver_down_pct <- summary_stats %>% filter(organ == "liver") %>% pull(pct_down)
brain_down_pct <- summary_stats %>% filter(organ == "brain") %>% pull(pct_down)
liver_down_strict <- summary_stats %>% filter(organ == "liver") %>% pull(n_down_strict)
brain_down_strict <- summary_stats %>% filter(organ == "brain") %>% pull(n_down_strict)

# ============================================================
# 5. 可视化: liver vs brain 的 Zbtb7a activity 横向对比
# ============================================================
cat("\n### 5. Visualization ###\n")

# Bar plot: 每个 cell type 的 Zbtb7a score, organ 着色, 按 score 排序
tvc_plot <- tvc %>%
  filter(!is.na(score)) %>%
  arrange(organ, score) %>%
  mutate(
    label = paste0(organ, ":", celltype),
    label = factor(label, levels = label)
  )

p1 <- ggplot(tvc_plot, aes(x = score, y = label, fill = organ)) +
  geom_col() +
  scale_fill_manual(values = c(liver = "#E41A1C", brain = "#377EB8")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = c(-3, 3), linetype = "dotted", color = "gray70") +
  labs(
    title = "Zbtb7a TF activity in Tumor_vs_Control\n(Liver vs Brain 跨器官对比)",
    subtitle = sprintf("Liver: %d/%d celltypes down (%s%% strict sig)\nBrain: %d/%d celltypes down (%s%% strict sig)",
                         sum(tvc_liver$score < 0, na.rm = TRUE), nrow(tvc_liver), 
                         liver_down_strict,
                         sum(tvc_brain$score < 0, na.rm = TRUE), nrow(tvc_brain),
                         brain_down_strict),
    x = "Zbtb7a activity score (negative = down-regulated)",
    y = "Cell type",
    fill = "Organ"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "gray30"),
    panel.grid.major.y = element_blank()
  )

ggsave(file.path(OUT_DIR, "07_zbtb7a_organ_scatter.pdf"),
       p1, width = 10, height = 10, dpi = 150)
cat("  → 07_zbtb7a_organ_scatter.pdf\n")

# ============================================================
# 6. 报告 summary.md
# ============================================================
cat("\n### 6. 生成 summary report ###\n")

summary_lines <- c(
  "# Zbtb7a: Liver vs Brain 跨器官对比 (Task A v1.2 新论点候选)",
  "",
  sprintf("生成时间: %s", format(Sys.time(), "%Y-%m-%d %H:%M")),
  "",
  "## 假设",
  "",
  "Task A v1.1 Section 4.4 发现 Zbtb7a 在 **12/22 brain cell type** 活性一致下调.",
  "已知 Zbtb7a (Liu 2014 Cell Metab) 调控肝 gluconeogenesis (Pck1) 代谢.",
  "Thyroid hormone 也调控 Pck1, 因此 **Zbtb7a 可能是新的肝脑轴候选论点**.",
  "",
  "本分析检验: **Zbtb7a 是否在 liver 也呈系统性下调** (类似 brain 的 pan-cell-type pattern)?",
  "",
  "## 关键数据",
  "",
  sprintf("**Tumor_vs_Control 对照下的 Zbtb7a activity**:"),
  ""
)

# Liver table
summary_lines <- c(summary_lines,
                    "### Liver (8 cell types)",
                    "",
                    "| Cell type | Zbtb7a score | padj | sig_strict | sig_relaxed |",
                    "|-----------|-------------:|-----:|-----------|-------------|")
for (i in seq_len(nrow(tvc_liver))) {
  r <- tvc_liver[i, ]
  if (is.na(r$score)) next
  summary_lines <- c(summary_lines,
                      sprintf("| %s | %.2f | %.3g | %s | %s |",
                              r$celltype, r$score, r$padj,
                              ifelse(is.na(r$sig_strict), "-", r$sig_strict),
                              ifelse(is.na(r$sig_relaxed), "-", r$sig_relaxed)))
}

summary_lines <- c(summary_lines, "",
                    "### Brain (22 cell types, clean)",
                    "",
                    "| Cell type | Zbtb7a score | padj | sig_strict | sig_relaxed |",
                    "|-----------|-------------:|-----:|-----------|-------------|")
for (i in seq_len(nrow(tvc_brain))) {
  r <- tvc_brain[i, ]
  if (is.na(r$score)) next
  summary_lines <- c(summary_lines,
                      sprintf("| %s | %.2f | %.3g | %s | %s |",
                              r$celltype, r$score, r$padj,
                              ifelse(is.na(r$sig_strict), "-", r$sig_strict),
                              ifelse(is.na(r$sig_relaxed), "-", r$sig_relaxed)))
}

# 关键结论
conclusion <- if (liver_down_pct >= 70 && liver_down_strict >= 2) {
  "✅ **假设成立**: Zbtb7a 在 liver 也呈系统性下调, 与 brain 类似. 这是**新的肝脑轴论点候选**, 值得写入 Task A v1.2."
} else if (liver_down_pct >= 50) {
  "⚠️ **假设部分成立**: Liver 有一定比例 cell type 呈下调, 但不如 brain 显著. 可作为辅助论点."
} else {
  "❌ **假设不成立**: Zbtb7a 在 liver 未呈广泛下调. Pan-brain signature 可能是 brain 特异性."
}

summary_lines <- c(summary_lines, "",
                    "## 结果解读",
                    "",
                    "### 跨器官统计",
                    "",
                    "| Organ | # cell types | # down | # down strict | # down relaxed | % down | median score |",
                    "|-------|-------------:|-------:|--------------:|---------------:|-------:|-------------:|")
for (i in seq_len(nrow(summary_stats))) {
  r <- summary_stats[i, ]
  summary_lines <- c(summary_lines,
                      sprintf("| %s | %d | %d | %d | %d | %.1f%% | %.2f |",
                              r$organ, r$n_celltypes, r$n_down, 
                              r$n_down_strict, r$n_down_relaxed,
                              r$pct_down, r$median_score))
}

summary_lines <- c(summary_lines, "",
                    "### 结论",
                    "",
                    conclusion,
                    "",
                    "## 如果肝脑共同下调, 潜在的机制解释",
                    "",
                    "1. **HCC-PDX 的全身性代谢应激**: 肿瘤消耗能量, 肝端和脑端都启动代谢重编程",
                    "   - Zbtb7a 作为肝 gluconeogenesis repressor, 下调 → Pck1/G6pc ↑",
                    "   - Zbtb7a 作为 brain 多细胞类型 transcription repressor, 下调 → 广泛 target 释放",
                    "",
                    "2. **与 Task A 肝脑轴 2 (Thyroid) 的协同**:",
                    "   - Thyroid hormone 通过 THR 调控 Pck1 (肝) + 神经发育 (脑)",
                    "   - Zbtb7a 同时在这两个器官下调, 可能是 Thyroid 信号的上游或下游 effector",
                    "",
                    "3. **潜在 Task A v1.2 新论点 X4 'Zbtb7a 肝脑共下调轴'**:",
                    "   - Liver 8/N 个 cell type Zbtb7a activity ↓",
                    "   - Brain 12/22 cell type Zbtb7a activity ↓ (已有)",
                    "   - 涉及 Thyroid-responsive gene Pck1",
                    "   - 可能是 metabolic-neural cross-talk 的 new marker",
                    "",
                    "## 下一步",
                    "",
                    "- 如果成立, 写入 Task A v1.2 Section 5.4 (肝脑轴 4)",
                    "- Wet-lab 验证: serum Zbtb7a target protein (e.g. Pck1) 测量",
                    "- 跨器官 DEG 比较: liver Pck1 vs brain Pck1 表达变化",
                    "",
                    "## 输出",
                    "",
                    "- `05_zbtb7a_liver_vs_brain.csv` - 所有 cell type 的 Zbtb7a activity",
                    "- `06_zbtb7a_organ_comparison.md` - 本报告",
                    "- `07_zbtb7a_organ_scatter.pdf` - 可视化对比图"
)

writeLines(summary_lines, file.path(OUT_DIR, "06_zbtb7a_organ_comparison.md"))
cat(sprintf("  → 06_zbtb7a_organ_comparison.md (%d lines)\n", length(summary_lines)))

# ============================================================
cat("\n============================================\n")
cat("Task 1.2 完成!\n")
cat(sprintf("Output: %s\n", OUT_DIR))
cat("\n核心数字:\n")
cat(sprintf("  Liver Zbtb7a (Tumor_vs_Control):  %d/%d down (%.0f%% negative, strict sig: %d)\n",
            sum(tvc_liver$score < 0, na.rm = TRUE), 
            nrow(tvc_liver),
            liver_down_pct,
            liver_down_strict))
cat(sprintf("  Brain Zbtb7a (Tumor_vs_Control):  %d/%d down (%.0f%% negative, strict sig: %d)\n",
            sum(tvc_brain$score < 0, na.rm = TRUE),
            nrow(tvc_brain),
            brain_down_pct,
            brain_down_strict))
cat(sprintf("\n结论: %s\n", conclusion))
cat("============================================\n")

