# ==================== 判断样本中人鼠基因比例 ====================
# 文件名: check_human_mouse_ratio.R
# 用法: Rscript check_human_mouse_ratio.R

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

# ==================== 工作目录 ====================
setwd("/mnt/BioDisk/wangmingyuan/SGS10_analysis")

# 创建输出目录
output_dir <- "human_mouse_check"
if (!dir.exists(output_dir)) dir.create(output_dir)

# 定义需要检查的样本
samples_to_check <- c(
  "gan56", "liver3", "liver8",
  "nao56", "brain3", "brain8",
  "con_brain", "bla_brain",
  "bla_liver", "con_liver"
)

# ==================== 定义分类函数 ====================

# 方法1: 基于基因命名规则（全大写=人源，首字母大写+小写=鼠源）
classify_by_name <- function(obj, sample_name) {
  cat("\n方法1: 基于基因命名规则\n")
  all_genes <- rownames(obj)
  
  # 人源基因：全大写（如 ALB, AFP, HNF4A）
  human_genes <- all_genes[grepl("^[A-Z][A-Z0-9-]+$", all_genes)]
  
  # 鼠源基因：首字母大写，其余小写（如 Alb, Afp, Hnf4a）
  mouse_genes <- all_genes[grepl("^[A-Z][a-z]", all_genes)]
  
  # 其他基因（如线粒体基因、非编码RNA等）
  other_genes <- setdiff(all_genes, c(human_genes, mouse_genes))
  
  cat(sprintf("  人源基因数: %d (%.2f%%)\n", 
              length(human_genes), length(human_genes)/length(all_genes)*100))
  cat(sprintf("  鼠源基因数: %d (%.2f%%)\n", 
              length(mouse_genes), length(mouse_genes)/length(all_genes)*100))
  cat(sprintf("  其他基因数: %d (%.2f%%)\n", 
              length(other_genes), length(other_genes)/length(all_genes)*100))
  
  return(list(
    human = length(human_genes),
    mouse = length(mouse_genes),
    other = length(other_genes),
    total = length(all_genes),
    human_genes = human_genes,
    mouse_genes = mouse_genes
  ))
}

# 方法2: 基于表达量（计算每个细胞中人源和小鼠源基因的表达比例）
classify_by_expression <- function(obj, sample_name) {
  cat("\n方法2: 基于表达量（每个细胞）\n")
  
  all_genes <- rownames(obj)
  human_genes <- all_genes[grepl("^[A-Z][A-Z0-9-]+$", all_genes)]
  mouse_genes <- all_genes[grepl("^[A-Z][a-z]", all_genes)]
  
  # 取交集（确保基因存在于表达矩阵中）
  human_genes_present <- intersect(human_genes, rownames(obj))
  mouse_genes_present <- intersect(mouse_genes, rownames(obj))
  
  expr <- GetAssayData(obj, layer = "data")
  
  # 计算每个细胞中人源基因的总表达量
  if (length(human_genes_present) > 0) {
    human_expr <- colSums(expr[human_genes_present, , drop = FALSE])
  } else {
    human_expr <- rep(0, ncol(obj))
  }
  
  # 计算每个细胞中鼠源基因的总表达量
  if (length(mouse_genes_present) > 0) {
    mouse_expr <- colSums(expr[mouse_genes_present, , drop = FALSE])
  } else {
    mouse_expr <- rep(0, ncol(obj))
  }
  
  total_expr <- human_expr + mouse_expr
  human_ratio <- ifelse(total_expr > 0, human_expr / total_expr, 0)
  
  # 分类细胞
  cell_type <- ifelse(human_ratio >= 0.8, "Human",
               ifelse(human_ratio <= 0.2, "Mouse", "Mixed"))
  
  cat(sprintf("  人源细胞 (Human): %d (%.2f%%)\n", 
              sum(cell_type == "Human"), mean(cell_type == "Human")*100))
  cat(sprintf("  鼠源细胞 (Mouse): %d (%.2f%%)\n", 
              sum(cell_type == "Mouse"), mean(cell_type == "Mouse")*100))
  cat(sprintf("  混合细胞 (Mixed): %d (%.2f%%)\n", 
              sum(cell_type == "Mixed"), mean(cell_type == "Mixed")*100))
  
  return(list(
    human_cells = sum(cell_type == "Human"),
    mouse_cells = sum(cell_type == "Mouse"),
    mixed_cells = sum(cell_type == "Mixed"),
    total_cells = ncol(obj),
    human_ratio_mean = mean(human_ratio),
    human_ratio_median = median(human_ratio),
    cell_type = cell_type,
    human_ratio = human_ratio
  ))
}

# 方法3: 基于已知标记基因
classify_by_markers <- function(obj, sample_name) {
  cat("\n方法3: 基于已知标记基因\n")
  
  # 人源标记基因
  human_markers <- c("ALB", "AFP", "TTR", "HNF4A", "HNF1A", "GPC3", "EPCAM")
  # 鼠源标记基因
  mouse_markers <- c("Alb", "Afp", "Ttr", "Hnf4a", "Hnf1a", "Gpc3", "Epcam")
  
  human_markers_present <- intersect(human_markers, rownames(obj))
  mouse_markers_present <- intersect(mouse_markers, rownames(obj))
  
  expr <- GetAssayData(obj, layer = "data")
  
  # 计算标记基因表达
  if (length(human_markers_present) > 0) {
    human_marker_expr <- colMeans(expr[human_markers_present, , drop = FALSE])
  } else {
    human_marker_expr <- rep(0, ncol(obj))
  }
  
  if (length(mouse_markers_present) > 0) {
    mouse_marker_expr <- colMeans(expr[mouse_markers_present, , drop = FALSE])
  } else {
    mouse_marker_expr <- rep(0, ncol(obj))
  }
  
  # 分类
  cell_type <- ifelse(human_marker_expr > mouse_marker_expr, "Human",
               ifelse(mouse_marker_expr > human_marker_expr, "Mouse", "Unknown"))
  
  cat(sprintf("  人源细胞 (Human): %d (%.2f%%)\n", 
              sum(cell_type == "Human"), mean(cell_type == "Human")*100))
  cat(sprintf("  鼠源细胞 (Mouse): %d (%.2f%%)\n", 
              sum(cell_type == "Mouse"), mean(cell_type == "Mouse")*100))
  cat(sprintf("  未知细胞 (Unknown): %d (%.2f%%)\n", 
              sum(cell_type == "Unknown"), mean(cell_type == "Unknown")*100))
  
  return(list(
    human_cells = sum(cell_type == "Human"),
    mouse_cells = sum(cell_type == "Mouse"),
    unknown_cells = sum(cell_type == "Unknown"),
    total_cells = ncol(obj),
    human_marker_expr_mean = mean(human_marker_expr),
    mouse_marker_expr_mean = mean(mouse_marker_expr)
  ))
}

# ==================== 主分析 ====================
cat("========================================\n")
cat("判断样本中人鼠基因比例\n")
cat("========================================\n\n")

all_results <- list()

for (s in samples_to_check) {
  cat("\n========================================\n")
  cat("样本:", s)
  cat("\n========================================\n")
  
  fp <- file.path(s, paste0(s, "_manual_annotated.rds"))
  if (!file.exists(fp)) {
    cat("文件不存在:", fp, "\n")
    next
  }
  
  obj <- readRDS(fp)
  cat("细胞数:", ncol(obj), "\n")
  cat("基因数:", nrow(obj), "\n")
  
  # 方法1: 基因命名规则
  result1 <- classify_by_name(obj, s)
  
  # 方法2: 基于表达量
  result2 <- classify_by_expression(obj, s)
  
  # 方法3: 基于标记基因
  result3 <- classify_by_markers(obj, s)
  
  # 保存结果
  all_results[[s]] <- list(
    method1 = result1,
    method2 = result2,
    method3 = result3
  )
  
  # 可视化：每个细胞的人源比例分布
  pdf(file.path(output_dir, paste0(s, "_human_ratio_distribution.pdf")), width = 10, height = 6)
  hist(result2$human_ratio, breaks = 50, 
       main = paste0(s, " - Human Ratio Distribution"),
       xlab = "Human Gene Expression Ratio",
       ylab = "Cell Count",
       col = "steelblue")
  abline(v = 0.5, col = "red", lty = 2, lwd = 2)
  legend("topright", legend = c("Human (≥0.8)", "Mixed (0.2-0.8)", "Mouse (≤0.2)"),
         fill = c("red", "yellow", "blue"), cex = 0.8)
  dev.off()
  
  # 可视化：UMAP按细胞来源着色
  obj$cell_origin <- result2$cell_type
  obj$human_ratio <- result2$human_ratio
  
  pdf(file.path(output_dir, paste0(s, "_umap_by_origin.pdf")), width = 12, height = 5)
  p1 <- DimPlot(obj, reduction = "umap", group.by = "cell_origin", 
                label = TRUE, repel = TRUE) +
    ggtitle(paste0(s, " - Cell Origin (Expression-based)")) +
    scale_color_manual(values = c("Human" = "red", "Mouse" = "blue", "Mixed" = "purple"))
  
  p2 <- FeaturePlot(obj, features = "human_ratio", reduction = "umap") +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.5) +
    ggtitle(paste0(s, " - Human Ratio"))
  
  print(p1 | p2)
  dev.off()
}

# ==================== 生成汇总报告 ====================
cat("\n\n========================================\n")
cat("汇总报告\n")
cat("========================================\n")

# 创建汇总表
summary_df <- data.frame()

for (s in names(all_results)) {
  res <- all_results[[s]]
  
  # 方法1结果
  summary_df <- rbind(summary_df, data.frame(
    Sample = s,
    Method = "GeneName",
    Human_Genes = res$method1$human,
    Mouse_Genes = res$method1$mouse,
    Other_Genes = res$method1$other,
    Total_Genes = res$method1$total,
    Human_Pct_Genes = res$method1$human / res$method1$total * 100,
    Mouse_Pct_Genes = res$method1$mouse / res$method1$total * 100
  ))
  
  # 方法2结果
  summary_df <- rbind(summary_df, data.frame(
    Sample = s,
    Method = "Expression",
    Human_Genes = res$method2$human_cells,
    Mouse_Genes = res$method2$mouse_cells,
    Other_Genes = res$method2$mixed_cells,
    Total_Genes = res$method2$total_cells,
    Human_Pct_Genes = res$method2$human_cells / res$method2$total_cells * 100,
    Mouse_Pct_Genes = res$method2$mouse_cells / res$method2$total_cells * 100
  ))
  
  # 方法3结果
  summary_df <- rbind(summary_df, data.frame(
    Sample = s,
    Method = "Markers",
    Human_Genes = res$method3$human_cells,
    Mouse_Genes = res$method3$mouse_cells,
    Other_Genes = res$method3$unknown_cells,
    Total_Genes = res$method3$total_cells,
    Human_Pct_Genes = res$method3$human_cells / res$method3$total_cells * 100,
    Mouse_Pct_Genes = res$method3$mouse_cells / res$method3$total_cells * 100
  ))
}

# 保存汇总表
write.csv(summary_df, file.path(output_dir, "human_mouse_summary.csv"), row.names = FALSE)

# 打印汇总
print(summary_df)

# 绘制汇总图
library(ggplot2)

p_summary <- ggplot(summary_df, aes(x = Sample, y = Human_Pct_Genes, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Human Gene/Cell Percentage by Different Methods",
       x = "Sample", y = "Human Percentage (%)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_dir, "human_percentage_summary.pdf"), p_summary, width = 12, height = 6, dpi = 300)

cat("\n\n========================================\n")
cat("分析完成！\n")
cat("结果保存在:", output_dir, "\n")
cat("========================================\n")
