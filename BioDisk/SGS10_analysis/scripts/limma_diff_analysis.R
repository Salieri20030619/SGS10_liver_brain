# ==================== 使用limma进行差异分析 ====================
# 文件名: ~/SGS10_analysis/scripts/limma_diff_analysis.R
# 用法: Rscript limma_diff_analysis.R

library(Seurat)
library(limma)
library(edgeR)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(RColorBrewer)

# 设置工作目录
setwd("/mnt/BioDisk/wangmingyuan/SGS10_analysis")

# 创建输出目录
output_dir <- "limma_diff_results"
if(!dir.exists(output_dir)) dir.create(output_dir)

plot_dir <- file.path(output_dir, "plots")
if(!dir.exists(plot_dir)) dir.create(plot_dir)

# ==================== 参数设置 ====================
MIN_CELLS <- 10          # 最小细胞数阈值
LOG2FC_CUTOFF <- 1       # log2 fold change阈值（可调整）
P_VALUE_CUTOFF <- 0.05   # p值阈值
ADJ_P_CUTOFF <- 0.05     # 校正后p值阈值
MIN_PCT <- 0.1           # 最小表达百分比

# 备用阈值（如果筛选结果不好）
LOG2FC_CUTOFF_ALT <- 0.5
P_VALUE_CUTOFF_ALT <- 0.05

cat("========================================\n")
cat("limma差异分析\n")
cat("========================================\n")
cat("主阈值:\n")
cat("  - |log2FC| >", LOG2FC_CUTOFF, "\n")
cat("  - adj.P.Val <", ADJ_P_CUTOFF, "\n")
cat("备用阈值:\n")
cat("  - |log2FC| >", LOG2FC_CUTOFF_ALT, "\n")
cat("  - P.Value <", P_VALUE_CUTOFF_ALT, "\n")
cat("  - min.pct:", MIN_PCT, "\n")
cat("  - 最小细胞数:", MIN_CELLS, "\n\n")

# ==================== 定义样本对应关系 ====================

# 1. 肿瘤组 vs 对照组的比较
tumor_vs_control <- list(
  # gan56 vs 各对照
  list(name = "gan56_vs_nao56", group1 = "gan56", group2 = "nao56", type = "tumor_vs_brain"),
  list(name = "gan56_vs_bla_liver", group1 = "gan56", group2 = "bla_liver", type = "tumor_vs_liver"),
  list(name = "gan56_vs_con_liver", group1 = "gan56", group2 = "con_liver", type = "tumor_vs_liver"),
  
  # liver3 vs 各对照
  list(name = "liver3_vs_brain3", group1 = "liver3", group2 = "brain3", type = "tumor_vs_brain"),
  list(name = "liver3_vs_bla_liver", group1 = "liver3", group2 = "bla_liver", type = "tumor_vs_liver"),
  list(name = "liver3_vs_con_liver", group1 = "liver3", group2 = "con_liver", type = "tumor_vs_liver"),
  
  # liver8 vs 各对照
  list(name = "liver8_vs_brain8", group1 = "liver8", group2 = "brain8", type = "tumor_vs_brain"),
  list(name = "liver8_vs_bla_liver", group1 = "liver8", group2 = "bla_liver", type = "tumor_vs_liver"),
  list(name = "liver8_vs_con_liver", group1 = "liver8", group2 = "con_liver", type = "tumor_vs_liver")
)

# 2. 空白组之间的比较（4种组合）
control_vs_control <- list(
  list(name = "bla_liver_vs_bla_brain", group1 = "bla_liver", group2 = "bla_brain", type = "blanco_liver_vs_blanco_brain"),
  list(name = "con_liver_vs_con_brain", group1 = "con_liver", group2 = "con_brain", type = "control_liver_vs_control_brain"),
  list(name = "bla_liver_vs_con_liver", group1 = "bla_liver", group2 = "con_liver", type = "liver_blanco_vs_control"),
  list(name = "bla_brain_vs_con_brain", group1 = "bla_brain", group2 = "con_brain", type = "brain_blanco_vs_control")
)

# ==================== 加载所有样本 ====================
cat("\n========== 加载所有样本 ==========\n")

samples_to_load <- unique(c(
  "gan56", "liver3", "liver8",
  "nao56", "brain3", "brain8", "con_brain", "bla_brain",
  "bla_liver", "con_liver"
))

seurat_objs <- list()
for(s in samples_to_load) {
  cat("加载", s, "...\n")
  file_path <- file.path(s, paste0(s, "_manual_annotated.rds"))
  if(file.exists(file_path)) {
    obj <- readRDS(file_path)
    obj <- JoinLayers(obj)
    obj <- RenameCells(obj, add.cell.id = s)
    seurat_objs[[s]] <- obj
    cat("  ✓ 细胞数:", ncol(seurat_objs[[s]]), "\n")
  } else {
    cat("  ✗ 文件不存在:", file_path, "\n")
  }
}

# ==================== 区分人源和小鼠源细胞 ====================
human_markers <- c("ALB", "AFP", "TTR", "HNF4A", "HNF1A", "GPC3", "EPCAM")
mouse_markers <- c("Alb", "Afp", "Ttr", "Hnf4a", "Hnf1a", "Gpc3", "Epcam")

classify_cells <- function(obj, sample_name) {
  cat("\n分类", sample_name, "细胞...\n")
  
  expr <- GetAssayData(obj, layer = "data")
  
  human_genes <- intersect(human_markers, rownames(expr))
  if(length(human_genes) > 0) {
    human_score <- colMeans(expr[human_genes, , drop = FALSE])
  } else {
    human_score <- rep(0, ncol(obj))
  }
  
  mouse_genes <- intersect(mouse_markers, rownames(expr))
  if(length(mouse_genes) > 0) {
    mouse_score <- colMeans(expr[mouse_genes, , drop = FALSE])
  } else {
    mouse_score <- rep(0, ncol(obj))
  }
  
  cell_origin <- ifelse(human_score > mouse_score, "Human", 
                         ifelse(mouse_score > human_score, "Mouse", "Unknown"))
  
  obj$cell_origin <- cell_origin
  obj$human_score <- human_score
  obj$mouse_score <- mouse_score
  
  cat("  人源细胞:", sum(cell_origin == "Human"), "\n")
  cat("  小鼠细胞:", sum(cell_origin == "Mouse"), "\n")
  cat("  未知细胞:", sum(cell_origin == "Unknown"), "\n")
  
  return(obj)
}

# 对肿瘤样本进行分类
tumor_samples <- c("gan56", "liver3", "liver8")
for(s in tumor_samples) {
  if(s %in% names(seurat_objs)) {
    seurat_objs[[s]] <- classify_cells(seurat_objs[[s]], s)
  }
}

# ==================== 辅助函数：提取和合并 ====================
extract_and_merge <- function(obj1, obj2, cell_type = NULL, origin = NULL, group1_name = NULL, group2_name = NULL) {
  # 获取细胞
  if(!is.null(cell_type)) {
    cells1 <- WhichCells(obj1, expression = final_celltype == cell_type)
    cells2 <- WhichCells(obj2, expression = final_celltype == cell_type)
  } else {
    cells1 <- colnames(obj1)
    cells2 <- colnames(obj2)
  }
  
  # 筛选细胞来源（仅对肿瘤样本）
  if(!is.null(origin) && "cell_origin" %in% colnames(obj1@meta.data)) {
    cells1 <- intersect(cells1, WhichCells(obj1, expression = cell_origin == origin))
  }
  
  if(length(cells1) < MIN_CELLS || length(cells2) < MIN_CELLS) {
    return(NULL)
  }
  
  # 创建子集
  subset1 <- subset(obj1, cells = cells1)
  subset2 <- subset(obj2, cells = cells2)
  
  # 添加分组信息
  if(is.null(group1_name)) group1_name <- deparse(substitute(obj1))
  if(is.null(group2_name)) group2_name <- deparse(substitute(obj2))
  
  subset1$group <- group1_name
  subset2$group <- group2_name
  
  # 合并
  merged <- merge(subset1, subset2)
  merged <- JoinLayers(merged)
  
  return(merged)
}

# ==================== 使用limma进行差异分析 ====================
run_limma_diff <- function(merged_obj, comparison_name, output_dir, plot_dir, use_alt_threshold = FALSE) {
  
  # 设置log2FC和p值阈值
  if(use_alt_threshold) {
    log2fc_cutoff <- LOG2FC_CUTOFF_ALT
    p_cutoff <- P_VALUE_CUTOFF_ALT
    use_adj_p <- FALSE
    cat("    使用备用阈值: |log2FC| >", log2fc_cutoff, ", P.Value <", p_cutoff, "\n")
  } else {
    log2fc_cutoff <- LOG2FC_CUTOFF
    p_cutoff <- ADJ_P_CUTOFF
    use_adj_p <- TRUE
    cat("    使用主阈值: |log2FC| >", log2fc_cutoff, ", adj.P.Val <", p_cutoff, "\n")
  }
  
  # 设置身份
  Idents(merged_obj) <- merged_obj$group
  groups <- unique(merged_obj$group)
  if(length(groups) < 2) return(NULL)
  
  group1 <- groups[1]
  group2 <- groups[2]
  
  # 检查细胞数
  count1 <- sum(merged_obj$group == group1)
  count2 <- sum(merged_obj$group == group2)
  
  if(count1 < MIN_CELLS || count2 < MIN_CELLS) {
    cat(sprintf("  ⚠️ 细胞数不足 (G1=%d, G2=%d)，跳过\n", count1, count2))
    return(NULL)
  }
  
  # 提取表达矩阵
  expr_matrix <- GetAssayData(merged_obj, layer = "data")
  
  # 过滤低表达基因
  keep_genes <- rowSums(expr_matrix > 0) >= (MIN_PCT * ncol(expr_matrix))
  expr_matrix_filtered <- expr_matrix[keep_genes, , drop = FALSE]
  
  if(nrow(expr_matrix_filtered) < 100) {
    cat("  ⚠️ 过滤后基因数不足100，跳过\n")
    return(NULL)
  }
  
  # 创建分组向量
  group_vector <- factor(merged_obj$group, levels = c(group1, group2))
  
  # 构建设计矩阵
  design <- model.matrix(~ 0 + group_vector)
  colnames(design) <- levels(group_vector)
  contrast_matrix <- makeContrasts(contrasts = paste(group1, "-", group2), levels = design)
  
  # 使用voom进行转换
  dge <- DGEList(counts = expr_matrix_filtered)
  dge <- calcNormFactors(dge)
  v <- voom(dge, design, plot = FALSE)
  
  # 线性拟合
  fit <- lmFit(v, design)
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)
  
  # 提取结果
  results <- topTable(fit2, number = Inf, adjust.method = "BH")
  results$gene <- rownames(results)
  
  # 添加显著性标记
  if(use_adj_p) {
    results$sig <- ifelse(abs(results$logFC) > log2fc_cutoff & results$adj.P.Val < p_cutoff,
                          ifelse(results$logFC > 0, "Up", "Down"), "NS")
  } else {
    results$sig <- ifelse(abs(results$logFC) > log2fc_cutoff & results$P.Value < p_cutoff,
                          ifelse(results$logFC > 0, "Up", "Down"), "NS")
  }
  
  # 统计
  up_count <- sum(results$sig == "Up")
  down_count <- sum(results$sig == "Down")
  
  cat(sprintf("  ✓ %s: G1=%d, G2=%d, Up=%d, Down=%d\n", 
              comparison_name, count1, count2, up_count, down_count))
  
  # 如果使用主阈值且结果太少，尝试备用阈值
  if(!use_alt_threshold && (up_count + down_count) < 10) {
    cat("    主阈值结果较少，尝试使用备用阈值...\n")
    return(run_limma_diff(merged_obj, comparison_name, output_dir, plot_dir, use_alt_threshold = TRUE))
  }
  
  # 保存结果
  write.csv(results, file.path(output_dir, paste0(comparison_name, "_all_genes.csv")), row.names = FALSE)
  
  # 保存显著基因
  sig_genes <- results %>% filter(sig != "NS")
  if(nrow(sig_genes) > 0) {
    write.csv(sig_genes, file.path(output_dir, paste0(comparison_name, "_significant_genes.csv")), row.names = FALSE)
  }
  
  # ==================== 生成火山图 ====================
  if(nrow(results) > 0) {
    cat("    生成火山图...\n")
    
    results$log10_p <- -log10(ifelse(use_adj_p, results$adj.P.Val, results$P.Value))
    results$log10_p[is.infinite(results$log10_p)] <- max(results$log10_p[is.finite(results$log10_p)], na.rm = TRUE)
    
    # 标记top基因
    top_up <- results %>% filter(sig == "Up") %>% arrange(adj.P.Val) %>% head(10)
    top_down <- results %>% filter(sig == "Down") %>% arrange(adj.P.Val) %>% head(10)
    top_genes <- rbind(top_up, top_down)
    
    # 火山图
    p_volcano <- ggplot(results, aes(x = logFC, y = log10_p, color = sig)) +
      geom_point(alpha = 0.6, size = 0.8) +
      scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
      geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff), linetype = "dashed", color = "black", alpha = 0.5) +
      geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", color = "black", alpha = 0.5) +
      labs(title = comparison_name,
           x = expression(Log[2] ~ "Fold Change"),
           y = expression(-Log[10] ~ P - value)) +
      theme_minimal() +
      theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5, size = 10),
            axis.title = element_text(size = 10))
    
    # 添加基因标签
    if(nrow(top_genes) > 0) {
      p_volcano <- p_volcano + 
        geom_text(data = top_genes, 
                  aes(x = logFC, y = log10_p, label = gene),
                  size = 2.5, vjust = -0.5, hjust = 0.5, check_overlap = TRUE)
    }
    
    ggsave(file.path(plot_dir, paste0(comparison_name, "_volcano.pdf")), p_volcano, width = 10, height = 8, dpi = 300)
  }
  
  # ==================== 生成热图 ====================
  sig_genes_for_heatmap <- results %>% filter(sig != "NS") %>% arrange(desc(abs(logFC))) %>% head(30)
  
  if(nrow(sig_genes_for_heatmap) >= 2) {
    cat("    生成热图...\n")
    
    # 提取表达矩阵
    expr_heatmap <- expr_matrix_filtered[rownames(expr_matrix_filtered) %in% sig_genes_for_heatmap$gene, , drop = FALSE]
    
    if(nrow(expr_heatmap) >= 2) {
      # 计算z-score
      expr_z <- t(scale(t(as.matrix(expr_heatmap))))
      expr_z[is.na(expr_z)] <- 0
      
      # 获取分组信息
      annotation_df <- data.frame(Group = merged_obj$group)
      rownames(annotation_df) <- colnames(expr_z)
      
      # 设置颜色
      group_colors <- setNames(c("#E41A1C", "#377EB8"), unique(merged_obj$group))
      annotation_colors <- list(Group = group_colors)
      
      # 如果细胞太多，随机采样
      if(ncol(expr_z) > 200) {
        set.seed(42)
        sample_cols <- sample(colnames(expr_z), 200)
        expr_z <- expr_z[, sample_cols]
        annotation_df <- annotation_df[sample_cols, , drop = FALSE]
      }
      
      # 绘制热图
      pdf(file.path(plot_dir, paste0(comparison_name, "_heatmap.pdf")), width = 12, height = 10)
      pheatmap(expr_z,
               annotation_col = annotation_df,
               annotation_colors = annotation_colors,
               show_rownames = TRUE,
               show_colnames = FALSE,
               cluster_rows = TRUE,
               cluster_cols = TRUE,
               fontsize_row = 8,
               main = paste0(comparison_name, " - Top ", nrow(sig_genes_for_heatmap), " DEGs"),
               color = colorRampPalette(c("blue", "white", "red"))(100))
      dev.off()
    }
  }
  
  return(results)
}

# ==================== 分析特定细胞类型 ====================
analyze_cell_type_comparison <- function(comp, cell_type, seurat_objs, output_dir, plot_dir) {
  obj1 <- seurat_objs[[comp$group1]]
  obj2 <- seurat_objs[[comp$group2]]
  
  if(is.null(obj1) || is.null(obj2)) return(NULL)
  
  # 检查细胞类型是否存在
  if(!cell_type %in% unique(obj1$final_celltype)) return(NULL)
  if(!cell_type %in% unique(obj2$final_celltype)) return(NULL)
  
  # 合并
  merged <- tryCatch({
    extract_and_merge(obj1, obj2, cell_type = cell_type,
                       group1_name = comp$group1, group2_name = comp$group2)
  }, error = function(e) { NULL })
  
  if(is.null(merged)) return(NULL)
  
  # 差异分析
  result_name <- paste0(comp$name, "_", cell_type)
  markers <- run_limma_diff(merged, result_name, output_dir, plot_dir)
  
  return(markers)
}

# ==================== 主分析 ====================
cat("\n========== 开始差异分析 ==========\n\n")

# 定义要分析的细胞类型
cell_types <- c("Hepatocytes", "Kupffer", "Endothelial", "T_cells", "B_cells", "Stellate")

# 存储结果
all_results <- list()
successful_comparisons <- list()

# ==================== 1. 肿瘤组 vs 对照组 ====================
cat("\n========== 1. 肿瘤组 vs 对照组分析 ==========\n")

for(comp in tumor_vs_control) {
  cat("\n---", comp$name, "---\n")
  
  obj1 <- seurat_objs[[comp$group1]]
  obj2 <- seurat_objs[[comp$group2]]
  
  if(is.null(obj1) || is.null(obj2)) {
    cat("  样本不存在\n")
    next
  }
  
  # 分析小鼠源细胞（肿瘤微环境）
  if("cell_origin" %in% colnames(obj1@meta.data)) {
    merged_mouse <- tryCatch({
      extract_and_merge(obj1, obj2, origin = "Mouse",
                         group1_name = comp$group1, group2_name = comp$group2)
    }, error = function(e) { NULL })
    
    if(!is.null(merged_mouse)) {
      result_name <- paste0(comp$name, "_mouse_origin")
      markers <- run_limma_diff(merged_mouse, result_name, output_dir, plot_dir)
      if(!is.null(markers)) {
        all_results[[result_name]] <- markers
        successful_comparisons <- c(successful_comparisons, result_name)
      }
    }
    
    # 分析人源细胞（肿瘤细胞）
    merged_human <- tryCatch({
      extract_and_merge(obj1, obj2, origin = "Human",
                         group1_name = comp$group1, group2_name = comp$group2)
    }, error = function(e) { NULL })
    
    if(!is.null(merged_human)) {
      result_name <- paste0(comp$name, "_human_origin")
      markers <- run_limma_diff(merged_human, result_name, output_dir, plot_dir)
      if(!is.null(markers)) {
        all_results[[result_name]] <- markers
        successful_comparisons <- c(successful_comparisons, result_name)
      }
    }
  }
  
  # 按细胞类型分析
  for(ct in cell_types) {
    markers <- analyze_cell_type_comparison(comp, ct, seurat_objs, output_dir, plot_dir)
    if(!is.null(markers)) {
      result_name <- paste0(comp$name, "_", ct)
      all_results[[result_name]] <- markers
      successful_comparisons <- c(successful_comparisons, result_name)
    }
  }
}

# ==================== 2. 空白组之间的比较 ====================
cat("\n\n========== 2. 空白组之间比较 ==========\n")
cat("  1. bla_liver vs bla_brain (空白组肝 vs 空白组脑)\n")
cat("  2. con_liver vs con_brain (对照组肝 vs 对照组脑)\n")
cat("  3. bla_liver vs con_liver (空白组肝 vs 对照组肝)\n")
cat("  4. bla_brain vs con_brain (空白组脑 vs 对照组脑)\n\n")

for(comp in control_vs_control) {
  cat("\n---", comp$name, "(", comp$type, ")---\n")
  
  obj1 <- seurat_objs[[comp$group1]]
  obj2 <- seurat_objs[[comp$group2]]
  
  if(is.null(obj1) || is.null(obj2)) {
    cat("  样本不存在\n")
    next
  }
  
  # 全细胞分析
  merged_all <- tryCatch({
    extract_and_merge(obj1, obj2, group1_name = comp$group1, group2_name = comp$group2)
  }, error = function(e) { NULL })
  
  if(!is.null(merged_all)) {
    markers <- run_limma_diff(merged_all, comp$name, output_dir, plot_dir)
    if(!is.null(markers)) {
      all_results[[comp$name]] <- markers
      successful_comparisons <- c(successful_comparisons, comp$name)
    }
  }
  
  # 按细胞类型分析
  for(ct in cell_types) {
    markers <- analyze_cell_type_comparison(comp, ct, seurat_objs, output_dir, plot_dir)
    if(!is.null(markers)) {
      result_name <- paste0(comp$name, "_", ct)
      all_results[[result_name]] <- markers
      successful_comparisons <- c(successful_comparisons, result_name)
    }
  }
}

# ==================== 3. 生成汇总统计 ====================
cat("\n\n========== 3. 生成汇总统计 ==========\n")

# 创建汇总表
summary_df <- data.frame()

for(name in names(all_results)) {
  data <- all_results[[name]]
  
  # 提取显著基因
  sig_up <- sum(data$sig == "Up", na.rm = TRUE)
  sig_down <- sum(data$sig == "Down", na.rm = TRUE)
  
  # 提取比较类型
  if(grepl("mouse_origin", name)) {
    comp_type <- "Tumor_Microenvironment"
  } else if(grepl("human_origin", name)) {
    comp_type <- "Tumor_Cells"
  } else if(grepl("bla_liver_vs_bla_brain", name)) {
    comp_type <- "Blanco_Liver_vs_Blanco_Brain"
  } else if(grepl("con_liver_vs_con_brain", name)) {
    comp_type <- "Control_Liver_vs_Control_Brain"
  } else if(grepl("bla_liver_vs_con_liver", name)) {
    comp_type <- "Liver_Blanco_vs_Control"
  } else if(grepl("bla_brain_vs_con_brain", name)) {
    comp_type <- "Brain_Blanco_vs_Control"
  } else {
    comp_type <- "Cell_Type"
  }
  
  summary_df <- rbind(summary_df,
                       data.frame(
                         Comparison = name,
                         Type = comp_type,
                         Up = sig_up,
                         Down = sig_down,
                         Total = sig_up + sig_down
                       ))
}

# 保存汇总
write.csv(summary_df, file.path(output_dir, "diff_analysis_summary.csv"), row.names = FALSE)

# 生成汇总柱状图
if(nrow(summary_df) > 0) {
  p_summary <- ggplot(summary_df, aes(x = reorder(Comparison, Total), y = Total, fill = Type)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = "limma Differential Expression Analysis Summary",
         x = "Comparison",
         y = "Number of Significant DEGs") +
    theme_minimal() +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5))
  
  ggsave(file.path(plot_dir, "summary_barplot.pdf"), p_summary, width = 14, height = max(8, nrow(summary_df)*0.3), dpi = 300)
}

# ==================== 4. 生成报告 ====================
cat("\n\n========== 4. 生成报告 ==========\n")

report_file <- file.path(plot_dir, "limma_diff_report.txt")
sink(report_file)

cat("========================================\n")
cat("limma差异分析报告\n")
cat("========================================\n\n")
cat("分析时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("分析参数:\n")
cat("  - 主阈值: |log2FC| >", LOG2FC_CUTOFF, ", adj.P.Val <", ADJ_P_CUTOFF, "\n")
cat("  - 备用阈值: |log2FC| >", LOG2FC_CUTOFF_ALT, ", P.Value <", P_VALUE_CUTOFF_ALT, "\n")
cat("  - min.pct:", MIN_PCT, "\n")
cat("  - 最小细胞数:", MIN_CELLS, "\n\n")

cat("成功完成的分析:\n")
cat("----------------------------------------\n")
for(name in successful_comparisons) {
  cat("  ✓", name, "\n")
}
cat("\n")

cat("差异基因统计:\n")
cat("----------------------------------------\n")
for(i in 1:nrow(summary_df)) {
  cat(sprintf("%-60s: Up=%d, Down=%d, Total=%d\n", 
              summary_df$Comparison[i], 
              summary_df$Up[i], 
              summary_df$Down[i], 
              summary_df$Total[i]))
}

sink()
cat("  ✓ 报告已保存\n")

# ==================== 5. 最终输出 ====================
cat("\n\n========================================\n")
cat("分析完成！\n")
cat("========================================\n")
cat("结果保存在:\n")
cat("  - 差异基因CSV: ", output_dir, "\n")
cat("  - 火山图和热图: ", plot_dir, "\n")
cat("\n成功完成的比较数:", length(successful_comparisons), "\n")
cat("========================================\n")
