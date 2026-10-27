# ==================== 修正版limma差异分析（正确对应关系） ====================
# 文件名: limma_diff_analysis_corrected.R
# 用法: Rscript limma_diff_analysis_corrected.R

suppressPackageStartupMessages({
  library(Seurat)
  library(limma)
  library(edgeR)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(pheatmap)
  library(RColorBrewer)
  library(patchwork)
})

# ==================== 工作目录 ====================
setwd("/mnt/BioDisk/wangmingyuan/SGS10_analysis")

# 创建输出目录
output_dir <- "limma_diff_corrected"
dirs <- c(
  output_dir,
  file.path(output_dir, "plots/volcano"),
  file.path(output_dir, "plots/heatmap"),
  file.path(output_dir, "deg_tables")
)
for (d in dirs) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

plot_dir <- file.path(output_dir, "plots")
deg_dir  <- file.path(output_dir, "deg_tables")

# ==================== 参数设置 ====================
LOG2FC_CUTOFF <- 1
ADJ_P_CUTOFF  <- 0.05
LOG2FC_ALT    <- 0.5
P_ALT         <- 0.05
MIN_CELLS     <- 10
MIN_PCT       <- 0.1

cat("========================================\n")
cat("修正版limma差异分析\n")
cat("========================================\n")
cat("主阈值: |log2FC| >", LOG2FC_CUTOFF, " & adj.P.Val <", ADJ_P_CUTOFF, "\n")
cat("备用阈值: |log2FC| >", LOG2FC_ALT, " & P.Value <", P_ALT, "\n\n")

# ==================== 定义比较组 ====================

# 1. 肿瘤 vs 对应脑对照（配对）
tumor_vs_brain <- list(
  list(name = "gan56_vs_nao56",   tumor = "gan56", control = "nao56", type = "tumor_vs_brain"),
  list(name = "liver3_vs_brain3", tumor = "liver3", control = "brain3", type = "tumor_vs_brain"),
  list(name = "liver8_vs_brain8", tumor = "liver8", control = "brain8", type = "tumor_vs_brain")
)

# 2. 肿瘤 vs 肝空白
tumor_vs_bla_liver <- list(
  list(name = "gan56_vs_bla_liver",   tumor = "gan56", control = "bla_liver", type = "tumor_vs_blank_liver"),
  list(name = "liver3_vs_bla_liver",  tumor = "liver3", control = "bla_liver", type = "tumor_vs_blank_liver"),
  list(name = "liver8_vs_bla_liver",  tumor = "liver8", control = "bla_liver", type = "tumor_vs_blank_liver")
)

# 3. 肿瘤 vs 肝对照
tumor_vs_con_liver <- list(
  list(name = "gan56_vs_con_liver",   tumor = "gan56", control = "con_liver", type = "tumor_vs_control_liver"),
  list(name = "liver3_vs_con_liver",  tumor = "liver3", control = "con_liver", type = "tumor_vs_control_liver"),
  list(name = "liver8_vs_con_liver",  tumor = "liver8", control = "con_liver", type = "tumor_vs_control_liver")
)

# 4. 脑样本 vs 脑空白
brain_vs_bla_brain <- list(
  list(name = "nao56_vs_bla_brain",   brain = "nao56", control = "bla_brain", type = "brain_vs_blank"),
  list(name = "brain3_vs_bla_brain",  brain = "brain3", control = "bla_brain", type = "brain_vs_blank"),
  list(name = "brain8_vs_bla_brain",  brain = "brain8", control = "bla_brain", type = "brain_vs_blank")
)

# 5. 脑样本 vs 脑对照
brain_vs_con_brain <- list(
  list(name = "nao56_vs_con_brain",   brain = "nao56", control = "con_brain", type = "brain_vs_control"),
  list(name = "brain3_vs_con_brain",  brain = "brain3", control = "con_brain", type = "brain_vs_control"),
  list(name = "brain8_vs_con_brain",  brain = "brain8", control = "con_brain", type = "brain_vs_control")
)

# 6. 空白与对照之间的比较
blank_vs_control <- list(
  list(name = "bla_liver_vs_con_liver", group1 = "bla_liver", group2 = "con_liver", type = "liver_blank_vs_control"),
  list(name = "bla_brain_vs_con_brain", group1 = "bla_brain", group2 = "con_brain", type = "brain_blank_vs_control"),
  list(name = "bla_liver_vs_bla_brain", group1 = "bla_liver", group2 = "bla_brain", type = "liver_blank_vs_brain_blank"),
  list(name = "con_liver_vs_con_brain", group1 = "con_liver", group2 = "con_brain", type = "liver_control_vs_brain_control")
)

# 合并所有比较
all_comparisons <- c(tumor_vs_brain, tumor_vs_bla_liver, tumor_vs_con_liver,
                     brain_vs_bla_brain, brain_vs_con_brain, blank_vs_control)

# 分析的细胞类型
cell_types <- c("Hepatocytes", "Kupffer", "Endothelial", "T_cells", "B_cells", "Stellate")

# ==================== 1. 加载所有样本 ====================
cat("\n========== 1. 加载所有样本 ==========\n")

samples_to_load <- c(
  "gan56", "liver3", "liver8",
  "nao56", "brain3", "brain8",
  "con_brain", "bla_brain",
  "bla_liver", "con_liver"
)

seurat_objs <- list()
for (s in samples_to_load) {
  cat("加载", s, "...\n")
  fp <- file.path(s, paste0(s, "_manual_annotated.rds"))
  if (file.exists(fp)) {
    obj <- readRDS(fp)
    obj <- JoinLayers(obj)
    obj <- RenameCells(obj, add.cell.id = s)
    seurat_objs[[s]] <- obj
    cat("  ✓ 细胞数:", ncol(obj), "\n")
  } else {
    cat("  ✗ 文件不存在:", fp, "\n")
  }
}

# ==================== 2. 人鼠细胞分类（仅肿瘤样本） ====================
classify_cells_by_name <- function(obj, sample_name) {
  cat("\n分类", sample_name, "（基因命名规则法）...\n")
  all_genes   <- rownames(obj)
  human_genes <- all_genes[grepl("^[A-Z][A-Z0-9-]+$", all_genes)]
  mouse_genes <- all_genes[grepl("^[A-Z][a-z]", all_genes)]
  
  expr         <- GetAssayData(obj, layer = "data")
  human_score  <- colSums(expr[intersect(human_genes, rownames(expr)), , drop = FALSE] > 0)
  mouse_score  <- colSums(expr[intersect(mouse_genes,  rownames(expr)), , drop = FALSE] > 0)
  total_score  <- human_score + mouse_score
  human_ratio  <- ifelse(total_score > 0, human_score / total_score, 0)

  cell_origin  <- ifelse(human_ratio >= 0.6, "Human",
                   ifelse(human_ratio <= 0.4, "Mouse", "Unknown"))

  obj$cell_origin <- cell_origin
  obj$human_ratio <- human_ratio

  cat(sprintf("  Human: %d (%.1f%%)  Mouse: %d (%.1f%%)  Unknown: %d (%.1f%%)\n",
    sum(cell_origin=="Human"),   100*mean(cell_origin=="Human"),
    sum(cell_origin=="Mouse"),   100*mean(cell_origin=="Mouse"),
    sum(cell_origin=="Unknown"), 100*mean(cell_origin=="Unknown")))
  return(obj)
}

for (s in c("gan56", "liver3", "liver8")) {
  if (s %in% names(seurat_objs))
    seurat_objs[[s]] <- classify_cells_by_name(seurat_objs[[s]], s)
}

# ==================== 3. 辅助函数 ====================
extract_and_merge <- function(obj1, obj2, cell_type = NULL, origin = NULL,
                               group1_name = NULL, group2_name = NULL) {
  get_cells <- function(obj, ct) {
    if (!is.null(ct) && "final_celltype" %in% colnames(obj@meta.data))
      WhichCells(obj, expression = final_celltype == ct)
    else colnames(obj)
  }
  cells1 <- get_cells(obj1, cell_type)
  cells2 <- get_cells(obj2, cell_type)

  # 肿瘤样本只取小鼠源细胞（微环境）
  if (!is.null(origin) && "cell_origin" %in% colnames(obj1@meta.data))
    cells1 <- intersect(cells1, WhichCells(obj1, expression = cell_origin == origin))

  if (length(cells1) < MIN_CELLS || length(cells2) < MIN_CELLS) return(NULL)

  s1 <- subset(obj1, cells = cells1)
  s2 <- subset(obj2, cells = cells2)
  
  s1$group <- if(is.null(group1_name)) deparse(substitute(obj1)) else group1_name
  s2$group <- if(is.null(group2_name)) deparse(substitute(obj2)) else group2_name
  
  merged <- merge(s1, s2)
  merged <- JoinLayers(merged)
  
  # 过滤表达量为0的细胞
  expr <- GetAssayData(merged, layer = "data")
  lib_sizes <- colSums(expr)
  keep_cells <- lib_sizes > 0
  if (sum(!keep_cells) > 0) {
    merged <- subset(merged, cells = colnames(expr)[keep_cells])
  }
  
  return(merged)
}

# ==================== 4. 差异分析函数 ====================
run_limma_diff <- function(merged_obj, comp_name, use_alt = FALSE) {
  if (use_alt) {
    lfc_cut <- LOG2FC_ALT; p_cut <- P_ALT; use_adj <- FALSE
    cat("    [备用阈值] |log2FC|>", lfc_cut, " P.Value<", p_cut, "\n")
  } else {
    lfc_cut <- LOG2FC_CUTOFF; p_cut <- ADJ_P_CUTOFF; use_adj <- TRUE
    cat("    [主阈值]   |log2FC|>", lfc_cut, " adj.P.Val<", p_cut, "\n")
  }

  Idents(merged_obj) <- merged_obj$group
  groups <- unique(merged_obj$group)
  if (length(groups) < 2) return(NULL)
  g1 <- groups[1]; g2 <- groups[2]
  n1 <- sum(merged_obj$group == g1); n2 <- sum(merged_obj$group == g2)

  if (n1 < MIN_CELLS || n2 < MIN_CELLS) {
    cat(sprintf("  ⚠️ 细胞数不足 (G1=%d, G2=%d)，跳过\n", n1, n2))
    return(NULL)
  }

  expr <- GetAssayData(merged_obj, layer = "data")
  keep <- rowSums(expr > 0) >= (MIN_PCT * ncol(expr))
  expr_f <- expr[keep, , drop = FALSE]
  
  if (nrow(expr_f) < 100) {
    cat("    ⚠️ 过滤后基因数 < 100，跳过\n")
    return(NULL)
  }
  
  # 检查并过滤表达量为0的细胞
  lib_sizes <- colSums(expr_f)
  if (any(lib_sizes == 0)) {
    keep_cells <- lib_sizes > 0
    expr_f <- expr_f[, keep_cells, drop = FALSE]
    merged_obj <- subset(merged_obj, cells = colnames(expr_f))
    if (ncol(expr_f) < MIN_CELLS * 2) return(NULL)
  }

  # 构建设计矩阵
  gv <- factor(merged_obj$group, levels = c(g1, g2))
  design <- model.matrix(~ 0 + gv); colnames(design) <- levels(gv)
  cont <- makeContrasts(contrasts = paste(g1, "-", g2), levels = design)

  dge <- DGEList(counts = expr_f)
  dge <- calcNormFactors(dge)
  
  v <- tryCatch({
    voom(dge, design, plot = FALSE)
  }, error = function(e) {
    cat("    voom错误:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(v)) return(NULL)

  fit <- lmFit(v, design)
  fit2 <- contrasts.fit(fit, cont)
  fit2 <- eBayes(fit2)

  res <- topTable(fit2, number = Inf, adjust.method = "BH")
  if (is.null(res) || nrow(res) == 0) return(NULL)
  
  res$gene <- rownames(res)

  p_col <- if (use_adj) "adj.P.Val" else "P.Value"
  res$sig <- ifelse(abs(res$logFC) > lfc_cut & res[[p_col]] < p_cut,
               ifelse(res$logFC > 0, "Up", "Down"), "NS")

  n_up <- sum(res$sig == "Up"); n_down <- sum(res$sig == "Down"); n_total <- n_up + n_down
  cat(sprintf("  ✓ %s | G1=%d G2=%d | Up=%d Down=%d Total=%d\n",
              comp_name, n1, n2, n_up, n_down, n_total))

  if (!use_alt && n_total < 10) {
    cat("    ⚠️ DEG<10，自动切换备用阈值...\n")
    return(run_limma_diff(merged_obj, comp_name, use_alt = TRUE))
  }
  
  return(res)
}

# ==================== 5. 火山图 ====================
plot_volcano <- function(res, comp_name, lfc_cut, p_cut, use_adj, save_dir) {
  if (is.null(res) || nrow(res) == 0) return()
  
  p_col <- if (use_adj) "adj.P.Val" else "P.Value"
  y_lab <- if (use_adj) expression(-Log[10]~"adj.P.Val") else expression(-Log[10]~"P.Value")

  res$log10p <- -log10(res[[p_col]] + 1e-300)
  res$log10p[is.infinite(res$log10p)] <- max(res$log10p[is.finite(res$log10p)], na.rm = TRUE)

  top_genes <- bind_rows(
    res %>% filter(sig == "Up") %>% arrange(adj.P.Val) %>% head(10),
    res %>% filter(sig == "Down") %>% arrange(adj.P.Val) %>% head(10)
  )

  p <- ggplot(res, aes(x = logFC, y = log10p, color = sig)) +
    geom_point(alpha = 0.5, size = 0.7) +
    scale_color_manual(values = c("Up"="#E41A1C","Down"="#377EB8","NS"="grey70")) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", color = "black", linewidth = 0.4) +
    geom_hline(yintercept = -log10(p_cut), linetype = "dashed", color = "black", linewidth = 0.4) +
    geom_text_repel(data = top_genes, aes(label = gene), size = 2.5, max.overlaps = 25) +
    labs(title = comp_name, x = expression(Log[2]~"Fold Change"), y = y_lab) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, size = 10))

  ggsave(file.path(save_dir, "volcano", paste0(comp_name, "_volcano.pdf")), p, width = 8, height = 6, dpi = 300)
  cat("    ✓ 火山图已保存\n")
}

# ==================== 6. 热图 ====================
plot_heatmap <- function(res, merged_obj, comp_name, save_dir) {
  if (is.null(res) || nrow(res) == 0) return()
  
  top_genes <- res %>% filter(sig != "NS") %>% arrange(desc(abs(logFC))) %>% head(30)
  if (nrow(top_genes) < 2) return()

  expr <- GetAssayData(merged_obj, layer = "data")
  keep <- rowSums(expr > 0) >= (MIN_PCT * ncol(expr))
  expr_f <- as.matrix(expr[keep, ])
  genes_use <- intersect(top_genes$gene, rownames(expr_f))
  if (length(genes_use) < 2) return()

  expr_sub <- expr_f[genes_use, , drop = FALSE]
  expr_z <- t(scale(t(expr_sub))); expr_z[is.na(expr_z)] <- 0

  anno_df <- data.frame(Group = merged_obj$group, row.names = colnames(expr_z))
  groups <- unique(merged_obj$group)
  anno_col <- list(Group = setNames(c("#E41A1C","#377EB8"), groups))

  if (ncol(expr_z) > 200) {
    set.seed(42)
    samp <- unlist(lapply(groups, function(g) {
      cs <- colnames(expr_z)[merged_obj$group[colnames(expr_z)] == g]
      sample(cs, min(100, length(cs)))
    }))
    expr_z <- expr_z[, samp]; anno_df <- anno_df[samp, , drop = FALSE]
  }

  pdf(file.path(save_dir, "heatmap", paste0(comp_name, "_heatmap.pdf")), width = 12, height = 10)
  pheatmap(expr_z, annotation_col = anno_df, annotation_colors = anno_col,
           show_rownames = TRUE, show_colnames = FALSE, cluster_rows = TRUE, cluster_cols = TRUE,
           fontsize_row = 8, main = paste0(comp_name, " Top", nrow(top_genes), " DEGs"),
           color = colorRampPalette(c("#4575B4","white","#D73027"))(100))
  dev.off()
  cat("    ✓ 热图已保存\n")
}

# ==================== 7. 主分析流程 ====================
cat("\n========== 开始差异分析 ==========\n\n")

all_results <- list()
successful_comps <- character(0)

# 分析函数
analyze_comparison <- function(comp, is_tumor_comparison = FALSE, origin = NULL) {
  if (is_tumor_comparison) {
    obj1 <- seurat_objs[[comp$tumor]]
    obj2 <- seurat_objs[[comp$control]]
    group1_name <- comp$tumor
    group2_name <- comp$control
  } else if ("brain" %in% names(comp)) {
    obj1 <- seurat_objs[[comp$brain]]
    obj2 <- seurat_objs[[comp$control]]
    group1_name <- comp$brain
    group2_name <- comp$control
  } else {
    obj1 <- seurat_objs[[comp$group1]]
    obj2 <- seurat_objs[[comp$group2]]
    group1_name <- comp$group1
    group2_name <- comp$group2
  }
  
  if (is.null(obj1) || is.null(obj2)) return(NULL)
  
  for (ct in cell_types) {
    # 检查细胞类型是否存在
    if (!"final_celltype" %in% colnames(obj1@meta.data) || !ct %in% unique(obj1$final_celltype) ||
        !"final_celltype" %in% colnames(obj2@meta.data) || !ct %in% unique(obj2$final_celltype)) next
    
    comp_name <- paste0(comp$name, "_", ct)
    if (is_tumor_comparison) comp_name <- paste0(comp$name, "_mouse_", ct)
    
    cat("\n  >>", comp_name, "\n")
    
    # 合并细胞
    merged <- extract_and_merge(obj1, obj2, cell_type = ct, origin = origin,
                                 group1_name = group1_name, group2_name = group2_name)
    if (is.null(merged)) next
    
    # 差异分析
    res <- run_limma_diff(merged, comp_name)
    if (is.null(res)) next
    
    # 保存结果
    write.csv(res, file.path(deg_dir, paste0(comp_name, "_all_genes.csv")), row.names = FALSE)
    
    # 判断使用的阈值
    used_adj <- any(res$sig != "NS" & res$adj.P.Val < ADJ_P_CUTOFF)
    lfc_used <- if (any(abs(res$logFC[res$sig!="NS"]) >= LOG2FC_CUTOFF)) LOG2FC_CUTOFF else LOG2FC_ALT
    p_used <- if (used_adj) ADJ_P_CUTOFF else P_ALT
    
    # 绘图
    plot_volcano(res, comp_name, lfc_used, p_used, used_adj, plot_dir)
    plot_heatmap(res, merged, comp_name, plot_dir)
    
    all_results[[comp_name]] <- res
    successful_comps <- c(successful_comps, comp_name)
  }
  return(NULL)
}

# 执行所有比较
cat("\n===== 1. 肿瘤 vs 对应脑对照 =====\n")
for (comp in tumor_vs_brain) {
  cat("\n---", comp$name, "---\n")
  analyze_comparison(comp, is_tumor_comparison = TRUE, origin = "Mouse")
}

cat("\n\n===== 2. 肿瘤 vs 肝空白 =====\n")
for (comp in tumor_vs_bla_liver) {
  cat("\n---", comp$name, "---\n")
  analyze_comparison(comp, is_tumor_comparison = TRUE, origin = "Mouse")
}

cat("\n\n===== 3. 肿瘤 vs 肝对照 =====\n")
for (comp in tumor_vs_con_liver) {
  cat("\n---", comp$name, "---\n")
  analyze_comparison(comp, is_tumor_comparison = TRUE, origin = "Mouse")
}

cat("\n\n===== 4. 脑样本 vs 脑空白 =====\n")
for (comp in brain_vs_bla_brain) {
  cat("\n---", comp$name, "---\n")
  analyze_comparison(comp, is_tumor_comparison = FALSE)
}

cat("\n\n===== 5. 脑样本 vs 脑对照 =====\n")
for (comp in brain_vs_con_brain) {
  cat("\n---", comp$name, "---\n")
  analyze_comparison(comp, is_tumor_comparison = FALSE)
}

cat("\n\n===== 6. 空白与对照比较 =====\n")
for (comp in blank_vs_control) {
  cat("\n---", comp$name, "---\n")
  analyze_comparison(comp, is_tumor_comparison = FALSE)
}

# ==================== 8. 汇总统计 ====================
cat("\n\n========== 汇总统计 ==========\n")

summary_df <- data.frame()
for (nm in names(all_results)) {
  d <- all_results[[nm]]
  if (is.null(d)) next
  summary_df <- rbind(summary_df, data.frame(
    Comparison = nm,
    Up = sum(d$sig == "Up", na.rm = TRUE),
    Down = sum(d$sig == "Down", na.rm = TRUE),
    Total = sum(d$sig != "NS", na.rm = TRUE)
  ))
}

write.csv(summary_df, file.path(output_dir, "diff_analysis_summary.csv"), row.names = FALSE)

cat("\n成功完成的比较数:", length(successful_comps), "\n")
if (nrow(summary_df) > 0) {
  cat("\nDEG数量统计:\n")
  print(summary_df)
}

# ==================== 9. 报告 ====================
report_file <- file.path(output_dir, "analysis_report.txt")
sink(report_file)
cat("========================================\n")
cat("修正版limma差异分析报告\n")
cat("========================================\n")
cat("分析时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("比较类型:\n")
cat("  1. 肿瘤 vs 对应脑对照 (gan56→nao56, liver3→brain3, liver8→brain8)\n")
cat("  2. 肿瘤 vs 肝空白 (bla_liver)\n")
cat("  3. 肿瘤 vs 肝对照 (con_liver)\n")
cat("  4. 脑样本 vs 脑空白 (bla_brain)\n")
cat("  5. 脑样本 vs 脑对照 (con_brain)\n")
cat("  6. 空白与对照比较\n\n")
cat("主阈值: |log2FC| >", LOG2FC_CUTOFF, " & adj.P.Val <", ADJ_P_CUTOFF, "\n")
cat("备用阈值: |log2FC| >", LOG2FC_ALT, " & P.Value <", P_ALT, "\n\n")
cat("成功完成的比较数:", length(successful_comps), "\n\n")
cat("DEG统计:\n")
for (i in seq_len(nrow(summary_df))) {
  cat(sprintf("  %s: Up=%d, Down=%d, Total=%d\n", 
              summary_df$Comparison[i],
              summary_df$Up[i], summary_df$Down[i], summary_df$Total[i]))
}
sink()

cat("\n========================================\n")
cat("分析完成！\n")
cat("成功比较数:", length(successful_comps), "\n")
cat("输出目录: ", output_dir, "\n")
cat("========================================\n")
