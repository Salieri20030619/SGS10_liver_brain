#!/usr/bin/env Rscript
# =============================================================================
# 文件名: il6st_celltype_source.R
# 用途:   查看各组细胞类型组成，确定 Il6st 的表达来源
# 用法:   conda activate R_4.5.2 && Rscript il6st_celltype_source.R
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

# 检查可用的注释列
cat("肿瘤组 metadata 列名:\n")
cat(paste("  ", colnames(tumor_merged@meta.data), collapse="\n"), "\n\n")

cat("对照组 metadata 列名:\n")
cat(paste("  ", colnames(control_obj@meta.data), collapse="\n"), "\n\n")

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

# 查看细胞类型分布
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("步骤3: 细胞类型组成\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

# 肿瘤组细胞类型分布
tumor_celltypes <- table(tumor_merged@meta.data[[celltype_col]])
tumor_celltypes <- sort(tumor_celltypes, decreasing = TRUE)

cat("【肿瘤组】细胞类型分布:\n")
cat(sprintf("  %-25s %10s %10s\n", "CellType", "Count", "Percentage"))
cat("  ─────────────────────────────────────────────\n")
for (ct in names(tumor_celltypes)) {
  pct <- tumor_celltypes[ct] / ncol(tumor_merged) * 100
  cat(sprintf("  %-25s %10d %9.2f%%\n", ct, tumor_celltypes[ct], pct))
}

# 对照组细胞类型分布
control_celltypes <- table(control_obj@meta.data[[celltype_col]])
control_celltypes <- sort(control_celltypes, decreasing = TRUE)

cat("\n【对照组】细胞类型分布:\n")
cat(sprintf("  %-25s %10s %10s\n", "CellType", "Count", "Percentage"))
cat("  ─────────────────────────────────────────────\n")
for (ct in names(control_celltypes)) {
  pct <- control_celltypes[ct] / ncol(control_obj) * 100
  cat(sprintf("  %-25s %10d %9.2f%%\n", ct, control_celltypes[ct], pct))
}

# =============================================================================
# 5. Il6st 在各细胞类型中的表达
# =============================================================================
cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("步骤4: Il6st 在各细胞类型中的表达\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

# 合并肿瘤和对照组以便统一分析
combined_obj <- merge(tumor_merged, control_obj)
combined_obj <- JoinLayers(combined_obj)

# 确保 Il6st 存在
if (!"Il6st" %in% rownames(combined_obj)) {
  stop("Il6st 不在基因列表中")
}

# 获取 Il6st 表达数据
il6st_exp <- FetchData(combined_obj, vars = "Il6st")
colnames(il6st_exp) <- "Il6st"
il6st_exp$CellType <- combined_obj@meta.data[[celltype_col]]
il6st_exp$Group <- combined_obj$group

# 计算每个细胞类型的 Il6st 表达统计
celltype_stats <- il6st_exp %>%
  group_by(CellType, Group) %>%
  summarise(
    N_cells = n(),
    Pct_positive = mean(Il6st > 0) * 100,
    Mean_expr = mean(Il6st),
    Median_expr = median(Il6st),
    Max_expr = max(Il6st),
    .groups = 'drop'
  ) %>%
  arrange(desc(Pct_positive))

# 输出详细表格
cat("各细胞类型 Il6st 表达统计:\n")
cat("\n")
cat(sprintf("  %-25s %10s %12s %12s %12s\n", 
            "CellType", "Group", "N_cells", "Pct_pos%", "Mean_expr"))
cat("  ─────────────────────────────────────────────────────────────────\n")

for (i in 1:nrow(celltype_stats)) {
  cat(sprintf("  %-25s %10s %10d %11.1f%% %11.4f\n",
              celltype_stats$CellType[i],
              celltype_stats$Group[i],
              celltype_stats$N_cells[i],
              celltype_stats$Pct_positive[i],
              celltype_stats$Mean_expr[i]))
}

# =============================================================================
# 6. 计算各细胞类型中 Tumor vs Control 的差异
# =============================================================================
cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("步骤5: 各细胞类型中 Tumor vs Control 的 Il6st 差异\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

# 获取所有细胞类型
all_celltypes <- unique(il6st_exp$CellType)

# 存储结果
wilcox_results <- data.frame()

for (ct in all_celltypes) {
  # 提取该细胞类型的数据
  ct_data <- il6st_exp[il6st_exp$CellType == ct, ]
  
  # 确保有足够的细胞数
  n_tumor <- sum(ct_data$Group == "Tumor")
  n_control <- sum(ct_data$Group == "Control")
  
  if (n_tumor >= 3 && n_control >= 3) {
    # 提取表达值
    tumor_exp <- ct_data$Il6st[ct_data$Group == "Tumor"]
    control_exp <- ct_data$Il6st[ct_data$Group == "Control"]
    
    # Wilcoxon 检验
    wt <- wilcox.test(tumor_exp, control_exp)
    
    # 计算检测率
    pct_tumor <- mean(tumor_exp > 0) * 100
    pct_control <- mean(control_exp > 0) * 100
    
    # 计算平均表达量的 log2FC
    mean_tumor <- mean(tumor_exp)
    mean_control <- mean(control_exp)
    log2fc <- log2((mean_tumor + 0.001) / (mean_control + 0.001))
    
    # 计算效应量 (Cliff's Delta)
    u <- wt$statistic
    cliff_delta <- (2 * u) / (n_tumor * n_control) - 1
    
    wilcox_results <- rbind(wilcox_results, data.frame(
      CellType = ct,
      N_Tumor = n_tumor,
      N_Control = n_control,
      Pct_Tumor = pct_tumor,
      Pct_Control = pct_control,
      Mean_Tumor = mean_tumor,
      Mean_Control = mean_control,
      log2FC = log2fc,
      p_value = wt$p.value,
      Cliff_Delta = cliff_delta,
      stringsAsFactors = FALSE
    ))
  }
}

# 添加校正后的 p 值
wilcox_results$p_adj <- p.adjust(wilcox_results$p_value, method = "bonferroni")
wilcox_results$Significant <- ifelse(wilcox_results$p_adj < 0.05 & abs(wilcox_results$log2FC) > 0.2, 
                                      ifelse(wilcox_results$log2FC > 0, "Up", "Down"), "NS")

# 按 p 值排序
wilcox_results <- wilcox_results %>% arrange(p_value)

cat("各细胞类型 Il6st 差异分析结果 (Tumor vs Control):\n")
cat("\n")
cat(sprintf("  %-22s %8s %8s %10s %10s %8s %12s %8s\n",
            "CellType", "N_Tumor", "N_Control", "Pct_T%", "Pct_C%", "log2FC", "p_adj", "Sig"))
cat("  ────────────────────────────────────────────────────────────────────────────────\n")

for (i in 1:nrow(wilcox_results)) {
  sig_symbol <- ""
  if (wilcox_results$Significant[i] == "Up") sig_symbol <- "▲"
  else if (wilcox_results$Significant[i] == "Down") sig_symbol <- "▼"
  else sig_symbol <- " "
  
  cat(sprintf("  %-22s %8d %8d %9.1f%% %9.1f%% %8.3f %11.2e %s%s\n",
              substr(wilcox_results$CellType[i], 1, 22),
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
# 7. 可视化
# =============================================================================
cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("步骤6: 生成可视化图表\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

# 图1: 各细胞类型 Il6st 检测率对比
p1_data <- celltype_stats %>%
  filter(Group %in% c("Tumor", "Control")) %>%
  dcast(CellType ~ Group, value.var = "Pct_positive") %>%
  arrange(desc(Tumor))

p1 <- ggplot(p1_data, aes(x = reorder(CellType, Tumor), y = Tumor, fill = "Tumor")) +
  geom_bar(stat = "identity", alpha = 0.7) +
  geom_point(aes(y = Control, color = "Control"), size = 3) +
  coord_flip() +
  scale_fill_manual(values = c("Tumor" = "#E41A1C")) +
  scale_color_manual(values = c("Control" = "#377EB8")) +
  labs(title = "Il6st Detection Rate by Cell Type",
       subtitle = paste0("Tumor (bars) vs Control (points)"),
       x = NULL, y = "Il6st Positive Cells (%)") +
  theme_bw(base_size = 10) +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(output_dir, "Il6st_detection_rate_by_celltype.pdf"), 
       p1, width = 10, height = 8)
ggsave(file.path(output_dir, "Il6st_detection_rate_by_celltype.png"), 
       p1, width = 10, height = 8, dpi = 200)
cat("  ✓ 图1: Il6st_detection_rate_by_celltype.pdf/png\n")

# 图2: 显著差异的细胞类型火山图
sig_cells <- wilcox_results %>% filter(Significant != "NS")
if (nrow(sig_cells) > 0) {
  p2 <- ggplot(wilcox_results, aes(x = log2FC, y = -log10(p_adj), 
                                     color = Significant, size = N_Tumor + N_Control)) +
    geom_point(alpha = 0.7) +
    geom_text_repel(data = sig_cells, aes(label = CellType), size = 3, max.overlaps = 20) +
    geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed", linewidth = 0.3) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.3) +
    scale_color_manual(values = c("Up" = "#E41A1C", "Down" = "#377EB8", "NS" = "grey70")) +
    labs(title = "Il6st Differential Expression by Cell Type",
         subtitle = "Tumor vs Control | Points size = total cell count",
         x = "log2FC (Tumor/Control)", y = "-log10(adjusted p-value)") +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave(file.path(output_dir, "Il6st_DE_by_celltype.pdf"), 
         p2, width = 8, height = 6)
  ggsave(file.path(output_dir, "Il6st_DE_by_celltype.png"), 
         p2, width = 8, height = 6, dpi = 200)
  cat("  ✓ 图2: Il6st_DE_by_celltype.pdf/png\n")
}

# 图3: 主要细胞类型的 Il6st 表达小提琴图
# 选择细胞数最多的前8种细胞类型
top_celltypes <- names(sort(table(combined_obj@meta.data[[celltype_col]]), decreasing = TRUE))[1:8]
plot_data <- il6st_exp[il6st_exp$CellType %in% top_celltypes, ]

p3 <- ggplot(plot_data, aes(x = Group, y = Il6st, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.2, outlier.size = 0.3) +
  facet_wrap(~CellType, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = c("Tumor" = "#E41A1C", "Control" = "#377EB8")) +
  labs(title = "Il6st Expression in Major Cell Types",
       subtitle = "Tumor vs Control",
       y = "Il6st Expression Level", x = NULL) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_dir, "Il6st_violin_by_celltype.pdf"), 
       p3, width = 12, height = 8)
ggsave(file.path(output_dir, "Il6st_violin_by_celltype.png"), 
       p3, width = 12, height = 8, dpi = 200)
cat("  ✓ 图3: Il6st_violin_by_celltype.pdf/png\n")

# =============================================================================
# 8. 保存结果
# =============================================================================
# 保存细胞类型组成
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

# 保存 Il6st 表达统计
write.csv(celltype_stats, file.path(output_dir, "il6st_expression_by_celltype.csv"), row.names = FALSE)

# 保存差异分析结果
write.csv(wilcox_results, file.path(output_dir, "il6st_de_by_celltype.csv"), row.names = FALSE)

cat("\n✓ 所有数据已保存到 CSV 文件\n")

# =============================================================================
# 9. 总结报告
# =============================================================================
cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("步骤7: 总结报告\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

# 找出显著上调的细胞类型
sig_up <- wilcox_results %>% filter(Significant == "Up") %>% arrange(desc(log2FC))
sig_down <- wilcox_results %>% filter(Significant == "Down") %>% arrange(log2FC)

cat("【Il6st 显著上调的细胞类型】\n")
if (nrow(sig_up) > 0) {
  for (i in 1:nrow(sig_up)) {
    cat(sprintf("  • %s: log2FC = %.3f, p_adj = %.2e\n", 
                sig_up$CellType[i], sig_up$log2FC[i], sig_up$p_adj[i]))
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

# 找出检测率最高的细胞类型
cat("\n【Il6st 检测率最高的细胞类型 (Tumor)】\n")
top_pct <- celltype_stats %>% 
  filter(Group == "Tumor") %>% 
  arrange(desc(Pct_positive)) %>% 
  head(5)
for (i in 1:nrow(top_pct)) {
  cat(sprintf("  • %s: %.1f%% (%d cells)\n", 
              top_pct$CellType[i], top_pct$Pct_positive[i], top_pct$N_cells[i]))
}

# =============================================================================
# 10. 完成
# =============================================================================
cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║                        分析完成！                                    ║\n")
cat("╠══════════════════════════════════════════════════════════════════════╣\n")
cat(sprintf("║  输出目录: %-50s║\n", output_dir))
cat("╠══════════════════════════════════════════════════════════════════════╣\n")
cat("║  生成文件:                                                          ║\n")
cat("║    - tumor_celltype_composition.csv        肿瘤组细胞类型组成       ║\n")
cat("║    - control_celltype_composition.csv      对照组细胞类型组成       ║\n")
cat("║    - il6st_expression_by_celltype.csv      各细胞类型Il6st表达      ║\n")
cat("║    - il6st_de_by_celltype.csv              各细胞类型差异分析       ║\n")
cat("║    - Il6st_*.pdf/png                       可视化图表               ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n")
