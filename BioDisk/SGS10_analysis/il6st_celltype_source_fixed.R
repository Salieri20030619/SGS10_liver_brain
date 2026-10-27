#!/usr/bin/env Rscript
# =============================================================================
# 文件名: il6st_celltype_source_fixed.R
# 用途:   查看各组细胞类型组成，确定 Il6st 的表达来源（修复版）
# =============================================================================

cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║              Il6st 细胞来源分析 - 各细胞类型表达谱                    ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n\n")

# =============================================================================
# 1. 环境配置
# =============================================================================
.libPaths("/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/lib/R/library")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(reshape2)
  library(ggrepel)
})

cat("✓ 包加载完成\n\n")

# =============================================================================
# 2. 参数设置
# =============================================================================
setwd("/mnt/BioDisk/wangmingyuan/SGS10_analysis")

output_dir <- "il6st_source_analysis"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive=TRUE)

cat(sprintf("输出目录: %s\n\n", output_dir))

# =============================================================================
# 3. 加载样本数据
# =============================================================================
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("步骤1: 加载样本数据\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

load_sample <- function(sample_name) {
  fp <- file.path(sample_name, paste0(sample_name, "_manual_annotated.rds"))
  
  if (!file.exists(fp)) {
    cat(sprintf("  ✗ 文件不存在: %s\n", fp))
    return(NULL)
  }
  
  obj <- tryCatch(readRDS(fp), error = function(e) {
    cat(sprintf("  ✗ 读取失败: %s - %s\n", sample_name, e$message))
    return(NULL)
  })
  
  if (is.null(obj)) return(NULL)
  
  obj <- JoinLayers(obj)
  obj <- RenameCells(obj, add.cell.id = sample_name)
  cat(sprintf("  ✓ %-12s 总细胞数: %d\n", sample_name, ncol(obj)))
  
  return(obj)
}

# 样本列表
tumor_samples <- c("nao56", "brain3", "brain8")
control_sample <- "con_brain"

# 加载所有样本
tumor_objs <- list()
for (s in tumor_samples) {
  tumor_objs[[s]] <- load_sample(s)
}
tumor_objs <- Filter(Negate(is.null), tumor_objs)

control_obj <- load_sample(control_sample)

if (is.null(control_obj)) {
  stop("对照组加载失败")
}

if (length(tumor_objs) == 0) {
  stop("肿瘤组加载失败")
}

# 合并肿瘤组
tumor_merged <- Reduce(merge, tumor_objs)
tumor_merged <- JoinLayers(tumor_merged)
tumor_merged$group <- "Tumor"

control_obj$group <- "Control"

cat(sprintf("\n✓ 肿瘤组合并: %d 个细胞\n", ncol(tumor_merged)))
cat(sprintf("✓ 对照组: %d 个细胞\n\n", ncol(control_obj)))

# =============================================================================
# 4. 检查细胞类型注释
# =============================================================================
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("步骤2: 检查细胞类型注释\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

# 确定细胞类型注释列
celltype_col <- NULL
if ("final_celltype" %in% colnames(tumor_merged@meta.data)) {
  celltype_col <- "final_celltype"
} else if ("celltype" %in% colnames(tumor_merged@meta.data)) {
  celltype_col <- "celltype"
} else if ("CellType" %in% colnames(tumor_merged@meta.data)) {
  celltype_col <- "CellType"
}

if (is.null(celltype_col)) {
  cat("⚠️ 未找到细胞类型注释列！\n")
  cat("可用的列: ", paste(colnames(tumor_merged@meta.data), collapse=", "), "\n")
  stop("无法继续分析")
}

cat(sprintf("✓ 使用细胞类型注释列: %s\n\n", celltype_col))

# =============================================================================
# 5. 细胞类型组成
# =============================================================================
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("步骤3: 细胞类型组成\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

# 肿瘤组细胞类型分布
tumor_celltypes <- table(tumor_merged@meta.data[[celltype_col]])
tumor_celltypes <- sort(tumor_celltypes, decreasing = TRUE)

cat("【肿瘤组】细胞类型分布:\n")
cat(sprintf("  %-30s %10s %10s\n", "CellType", "Count", "Percentage"))
cat("  ─────────────────────────────────────────────────\n")
for (ct in names(tumor_celltypes)) {
  pct <- tumor_celltypes[ct] / ncol(tumor_merged) * 100
  cat(sprintf("  %-30s %10d %9.2f%%\n", ct, tumor_celltypes[ct], pct))
}

# 对照组细胞类型分布
control_celltypes <- table(control_obj@meta.data[[celltype_col]])
control_celltypes <- sort(control_celltypes, decreasing = TRUE)

cat("\n【对照组】细胞类型分布:\n")
cat(sprintf("  %-30s %10s %10s\n", "CellType", "Count", "Percentage"))
cat("  ─────────────────────────────────────────────────\n")
for (ct in names(control_celltypes)) {
  pct <- control_celltypes[ct] / ncol(control_obj) * 100
  cat(sprintf("  %-30s %10d %9.2f%%\n", ct, control_celltypes[ct], pct))
}

# =============================================================================
# 6. 合并数据并提取 Il6st 表达
# =============================================================================
cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("步骤4: 合并数据并提取 Il6st 表达\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

# 合并肿瘤和对照组
combined_obj <- merge(tumor_merged, control_obj)
combined_obj <- JoinLayers(combined_obj)

# 检查 Il6st
if (!"Il6st" %in% rownames(combined_obj)) {
  stop("Il6st 不在基因列表中")
}

# 提取表达数据
il6st_exp <- FetchData(combined_obj, vars = "Il6st")
colnames(il6st_exp) <- "Il6st"
il6st_exp$CellType <- combined_obj@meta.data[[celltype_col]]
il6st_exp$Group <- combined_obj$group

# =============================================================================
# 7. 各细胞类型 Il6st 表达统计
# =============================================================================
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("步骤5: 各细胞类型 Il6st 表达统计\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

celltype_stats <- il6st_exp %>%
  group_by(CellType, Group) %>%
  summarise(
    N_cells = n(),
    Pct_positive = mean(Il6st > 0) * 100,
    Mean_expr = mean(Il6st),
    Median_expr = median(Il6st),
    .groups = 'drop'
  ) %>%
  arrange(desc(Pct_positive))

cat("各细胞类型 Il6st 表达统计:\n\n")
cat(sprintf("  %-30s %10s %12s %12s\n", 
            "CellType", "Group", "N_cells", "Pct_pos%"))
cat("  ────────────────────────────────────────────────────────────────\n")

for (i in 1:min(30, nrow(celltype_stats))) {
  cat(sprintf("  %-30s %10s %10d %11.1f%%\n",
              substr(celltype_stats$CellType[i], 1, 30),
              celltype_stats$Group[i],
              celltype_stats$N_cells[i],
              celltype_stats$Pct_positive[i]))
}

# =============================================================================
# 8. 差异分析
# =============================================================================
cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("步骤6: 各细胞类型中 Tumor vs Control 的 Il6st 差异\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

all_celltypes <- unique(il6st_exp$CellType)
wilcox_results <- data.frame()

for (ct in all_celltypes) {
  ct_data <- il6st_exp[il6st_exp$CellType == ct, ]
  n_tumor <- sum(ct_data$Group == "Tumor")
  n_control <- sum(ct_data$Group == "Control")
  
  if (n_tumor >= 3 && n_control >= 3) {
    tumor_exp <- ct_data$Il6st[ct_data$Group == "Tumor"]
    control_exp <- ct_data$Il6st[ct_data$Group == "Control"]
    
    wt <- wilcox.test(tumor_exp, control_exp)
    
    pct_tumor <- mean(tumor_exp > 0) * 100
    pct_control <- mean(control_exp > 0) * 100
    mean_tumor <- mean(tumor_exp)
    mean_control <- mean(control_exp)
    log2fc <- log2((mean_tumor + 0.001) / (mean_control + 0.001))
    
    wilcox_results <- rbind(wilcox_results, data.frame(
      CellType = ct,
      N_Tumor = n_tumor,
      N_Control = n_control,
      Pct_Tumor = pct_tumor,
      Pct_Control = pct_control,
      log2FC = log2fc,
      p_value = wt$p.value,
      stringsAsFactors = FALSE
    ))
  }
}

wilcox_results$p_adj <- p.adjust(wilcox_results$p_value, method = "bonferroni")
wilcox_results$Significant <- ifelse(
  wilcox_results$p_adj < 0.05 & abs(wilcox_results$log2FC) > 0.2,
  ifelse(wilcox_results$log2FC > 0, "Up", "Down"), 
  "NS"
)
wilcox_results <- wilcox_results %>% arrange(p_value)

cat("各细胞类型 Il6st 差异分析结果 (Tumor vs Control):\n\n")
cat(sprintf("  %-30s %8s %8s %10s %10s %8s %12s %8s\n",
            "CellType", "N_Tumor", "N_Control", "Pct_T%", "Pct_C%", "log2FC", "p_adj", "Sig"))
cat("  ───────────────────────────────────────────────────────────────────────────────────────\n")

for (i in 1:nrow(wilcox_results)) {
  sig_symbol <- ifelse(wilcox_results$Significant[i] == "Up", "▲",
                      ifelse(wilcox_results$Significant[i] == "Down", "▼", " "))
  cat(sprintf("  %-30s %8d %8d %9.1f%% %9.1f%% %8.3f %11.2e %s %s\n",
              substr(wilcox_results$CellType[i], 1, 30),
              wilcox_results$N_Tumor[i],
              wilcox_results$N_Control[i],
              wilcox_results$Pct_Tumor[i],
              wilcox_results$Pct_Control[i],
              wilcox_results$log2FC[i],
              wilcox_results$p_adj[i],
              wilcox_results$Significant[i],
              sig_symbol))
}

# =============================================================================
# 9. 可视化
# =============================================================================
cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("步骤7: 生成可视化图表\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

# 图1: 检测率条形图
p1_data <- celltype_stats %>%
  filter(Group %in% c("Tumor", "Control")) %>%
  tidyr::pivot_wider(id_cols = CellType, names_from = Group, values_from = Pct_positive) %>%
  mutate(Tumor = ifelse(is.na(Tumor), 0, Tumor),
         Control = ifelse(is.na(Control), 0, Control)) %>%
  arrange(desc(Tumor)) %>%
  head(15)

p1_data_long <- p1_data %>%
  tidyr::pivot_longer(cols = c(Tumor, Control), names_to = "Group", values_to = "Pct_positive")

p1 <- ggplot(p1_data_long, aes(x = reorder(CellType, -Pct_positive), y = Pct_positive, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  coord_flip() +
  scale_fill_manual(values = c("Tumor" = "#E41A1C", "Control" = "#377EB8")) +
  labs(title = "Il6st Detection Rate by Cell Type",
       x = NULL, y = "Il6st Positive Cells (%)") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(output_dir, "Il6st_detection_rate.pdf"), p1, width = 10, height = 8)
ggsave(file.path(output_dir, "Il6st_detection_rate.png"), p1, width = 10, height = 8, dpi = 200)
cat("  ✓ 图1: Il6st_detection_rate.pdf/png\n")

# 图2: 火山图
sig_cells <- wilcox_results %>% filter(Significant != "NS")
p2 <- ggplot(wilcox_results, aes(x = log2FC, y = -log10(p_adj), color = Significant)) +
  geom_point(alpha = 0.7, size = 3) +
  geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed", linewidth = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.3) +
  scale_color_manual(values = c("Up" = "#E41A1C", "Down" = "#377EB8", "NS" = "grey70")) +
  labs(title = "Il6st Differential Expression by Cell Type",
       subtitle = "Tumor vs Control",
       x = "log2FC (Tumor/Control)", y = "-log10(adjusted p-value)") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

if (nrow(sig_cells) > 0) {
  p2 <- p2 + geom_text_repel(data = sig_cells, aes(label = CellType), 
                             size = 3, max.overlaps = 20)
}

ggsave(file.path(output_dir, "Il6st_DE_volcano.pdf"), p2, width = 8, height = 6)
ggsave(file.path(output_dir, "Il6st_DE_volcano.png"), p2, width = 8, height = 6, dpi = 200)
cat("  ✓ 图2: Il6st_DE_volcano.pdf/png\n")

# =============================================================================
# 10. 保存结果
# =============================================================================
write.csv(data.frame(
  CellType = names(tumor_celltypes),
  Tumor_Count = as.numeric(tumor_celltypes),
  Tumor_Percent = as.numeric(tumor_celltypes) / ncol(tumor_merged) * 100
), file.path(output_dir, "tumor_celltype_composition.csv"), row.names = FALSE)

write.csv(data.frame(
  CellType = names(control_celltypes),
  Control_Count = as.numeric(control_celltypes),
  Control_Percent = as.numeric(control_celltypes) / ncol(control_obj) * 100
), file.path(output_dir, "control_celltype_composition.csv"), row.names = FALSE)

write.csv(celltype_stats, file.path(output_dir, "il6st_expression_by_celltype.csv"), row.names = FALSE)
write.csv(wilcox_results, file.path(output_dir, "il6st_de_by_celltype.csv"), row.names = FALSE)

cat("\n✓ 所有数据已保存到 CSV 文件\n")

# =============================================================================
# 11. 总结报告
# =============================================================================
cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("步骤8: 总结报告\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

sig_up <- wilcox_results %>% filter(Significant == "Up") %>% arrange(desc(log2FC))
sig_down <- wilcox_results %>% filter(Significant == "Down") %>% arrange(log2FC)

cat("【Il6st 显著上调的细胞类型】\n")
if (nrow(sig_up) > 0) {
  for (i in 1:nrow(sig_up)) {
    cat(sprintf("  • %s: log2FC = %.3f, p_adj = %.2e, Tumor检测率 = %.1f%%, Control检测率 = %.1f%%\n", 
                sig_up$CellType[i], sig_up$log2FC[i], sig_up$p_adj[i],
                sig_up$Pct_Tumor[i], sig_up$Pct_Control[i]))
  }
} else {
  cat("  无\n")
}

cat("\n【Il6st 显著下调的细胞类型】\n")
if (nrow(sig_down) > 0) {
  for (i in 1:nrow(sig_down)) {
    cat(sprintf("  • %s: log2FC = %.3f, p_adj = %.2e\n", 
                sig_down$CellType[i], sig_down$log2FC[i], sig_down$p_adj[i]))
  }
} else {
  cat("  无\n")
}

cat("\n【Il6st 检测率最高的细胞类型 (Tumor)】\n")
top_pct <- celltype_stats %>% 
  filter(Group == "Tumor") %>% 
  arrange(desc(Pct_positive)) %>% 
  head(5)
for (i in 1:nrow(top_pct)) {
  cat(sprintf("  • %s: %.1f%% (%d cells)\n", 
              top_pct$CellType[i], top_pct$Pct_positive[i], top_pct$N_cells[i]))
}

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║                        分析完成！                                    ║\n")
cat("╠══════════════════════════════════════════════════════════════════════╣\n")
cat(sprintf("║  输出目录: %-50s║\n", output_dir))
cat("╚══════════════════════════════════════════════════════════════════════╝\n")
