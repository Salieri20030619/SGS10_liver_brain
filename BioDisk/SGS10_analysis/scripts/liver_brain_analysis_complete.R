# ==================== 肝脑肿瘤完整分析脚本（逐级宽松阈值版） ====================
# 文件名: liver_brain_analysis_complete.R
# 用法: Rscript liver_brain_analysis_complete.R

suppressPackageStartupMessages({
  library(Seurat)
  library(limma)
  library(edgeR)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(pheatmap)
  library(RColorBrewer)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
  library(patchwork)
})

# ==================== 工作目录 ====================
setwd("/mnt/BioDisk/wangmingyuan/SGS10_analysis")

output_dir <- "liver_brain_complete_results"
dirs <- c(
  output_dir,
  file.path(output_dir, "liver"),
  file.path(output_dir, "brain"),
  file.path(output_dir, "crosstalk"),
  file.path(output_dir, "plots/volcano"),
  file.path(output_dir, "plots/heatmap"),
  file.path(output_dir, "deg_tables"),
  file.path(output_dir, "threshold_summary"),
  file.path(output_dir, "enrichment")
)
for (d in dirs) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

# ==================== 参数设置 ====================
MIN_CELLS <- 10
MIN_PCT <- 0.1

# 阈值层级（从严格到宽松）
threshold_levels <- list(
  level1 = list(logFC = 1.0, p_cut = 0.05, use_adj = TRUE, name = "Very_Strict"),
  level2 = list(logFC = 0.8, p_cut = 0.05, use_adj = TRUE, name = "Strict"),
  level3 = list(logFC = 0.6, p_cut = 0.05, use_adj = TRUE, name = "Moderate"),
  level4 = list(logFC = 0.5, p_cut = 0.05, use_adj = FALSE, name = "Relaxed"),
  level5 = list(logFC = 0.4, p_cut = 0.05, use_adj = FALSE, name = "Very_Relaxed")
)

# 定义分析的细胞类型层级
# 层级1: 全部细胞（All_Cells）
# 层级2: 最主要细胞（Hepatocytes - 肝细胞）
# 层级3: 两种最主要亚型（Endothelial - 内皮细胞, Kupffer - 枯否细胞）
analysis_cell_types <- list(
  list(name = "All_Cells", type = "All_Cells", description = "全部细胞"),
  list(name = "Hepatocytes", type = "Cell_Type_Specific", description = "肝细胞（最主要细胞）"),
  list(name = "Endothelial", type = "Cell_Type_Specific", description = "内皮细胞（最主要亚型）"),
  list(name = "Kupffer", type = "Cell_Type_Specific", description = "枯否细胞（最主要亚型）")
)

cat("========================================\n")
cat("肝脑肿瘤完整分析（逐级递减阈值）\n")
cat("========================================\n")
cat("R版本:", version$version.string, "\n")
cat("clusterProfiler版本:", as.character(packageVersion("clusterProfiler")), "\n")
cat("阈值层级（自动选择第一个产生≥10个DEG的阈值）:\n")
for (i in 1:length(threshold_levels)) {
  lvl <- threshold_levels[[i]]
  cat(sprintf("  Level %d (%s): |log2FC| > %.1f, %s < 0.05\n", 
              i, lvl$name, lvl$logFC, ifelse(lvl$use_adj, "adj.P.Val", "P.Value")))
}
cat("\n分析的细胞类型:\n")
for (ct in analysis_cell_types) {
  cat(sprintf("  - %s: %s\n", ct$name, ct$description))
}
cat("最小细胞数:", MIN_CELLS, "\n\n")

# ==================== 定义比较组 ====================

# 1. 肝脏内部比较
liver_comparisons <- list(
  list(name = "tumor_vs_blank", tumor = c("gan56", "liver3", "liver8"), control = "bla_liver", type = "liver"),
  list(name = "tumor_vs_control", tumor = c("gan56", "liver3", "liver8"), control = "con_liver", type = "liver"),
  list(name = "blank_vs_control", group1 = "bla_liver", group2 = "con_liver", type = "liver")
)

# 2. 脑内部比较
brain_comparisons <- list(
  list(name = "brain_tumor_vs_blank", tumor_brain = c("nao56", "brain3", "brain8"), control = "bla_brain", type = "brain"),
  list(name = "brain_tumor_vs_control", tumor_brain = c("nao56", "brain3", "brain8"), control = "con_brain", type = "brain"),
  list(name = "brain_blank_vs_control", group1 = "bla_brain", group2 = "con_brain", type = "brain")
)

# 3. 跨器官比较（配对）
crosstalk_comparisons <- list(
  list(name = "gan56_vs_nao56", liver = "gan56", brain = "nao56"),
  list(name = "liver3_vs_brain3", liver = "liver3", brain = "brain3"),
  list(name = "liver8_vs_brain8", liver = "liver8", brain = "brain8")
)

# ==================== 1. 加载所有样本 ====================
cat("\n========== 1. 加载所有样本 ==========\n")

samples_to_load <- c("gan56", "liver3", "liver8", "nao56", "brain3", "brain8",
                     "con_brain", "bla_brain", "bla_liver", "con_liver")

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

# ==================== 2. 辅助函数 ====================
extract_cells <- function(obj, cell_type = NULL) {
  if (!is.null(cell_type) && cell_type != "All_Cells" && "final_celltype" %in% colnames(obj@meta.data)) {
    return(WhichCells(obj, expression = final_celltype == cell_type))
  }
  return(colnames(obj))
}

count_cells <- function(obj, cell_type = NULL) {
  return(length(extract_cells(obj, cell_type)))
}

# 合并多个对象（用于肿瘤组多个样本合并）
merge_objects <- function(obj_list, group_names, cell_type = NULL) {
  merged_list <- list()
  for (i in seq_along(obj_list)) {
    cells <- extract_cells(obj_list[[i]], cell_type)
    if (length(cells) < MIN_CELLS) return(NULL)
    subset_obj <- subset(obj_list[[i]], cells = cells)
    subset_obj$group <- group_names[i]
    merged_list[[i]] <- subset_obj
  }
  
  merged <- Reduce(merge, merged_list)
  merged <- JoinLayers(merged)
  
  expr <- GetAssayData(merged, layer = "data")
  lib_sizes <- colSums(expr)
  keep_cells <- lib_sizes > 0
  if (sum(!keep_cells) > 0) {
    merged <- subset(merged, cells = colnames(expr)[keep_cells])
  }
  
  return(merged)
}

# 合并两个对象
merge_two_objects <- function(obj1, obj2, cell_type = NULL, group1_name, group2_name) {
  cells1 <- extract_cells(obj1, cell_type)
  cells2 <- extract_cells(obj2, cell_type)
  
  if (length(cells1) < MIN_CELLS || length(cells2) < MIN_CELLS) {
    return(NULL)
  }
  
  subset1 <- subset(obj1, cells = cells1)
  subset1$group <- group1_name
  subset2 <- subset(obj2, cells = cells2)
  subset2$group <- group2_name
  
  merged <- merge(subset1, subset2)
  merged <- JoinLayers(merged)
  
  expr <- GetAssayData(merged, layer = "data")
  lib_sizes <- colSums(expr)
  keep_cells <- lib_sizes > 0
  if (sum(!keep_cells) > 0) {
    merged <- subset(merged, cells = colnames(expr)[keep_cells])
  }
  
  return(merged)
}

# ==================== 3. 差异分析函数（逐级阈值） ====================
run_limma_diff_adaptive <- function(merged_obj, comp_name, group1, group2, cell_type_name) {
  Idents(merged_obj) <- merged_obj$group
  n1 <- sum(merged_obj$group == group1)
  n2 <- sum(merged_obj$group == group2)
  
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
  
  lib_sizes <- colSums(expr_f)
  if (any(lib_sizes == 0)) {
    keep_cells <- lib_sizes > 0
    expr_f <- expr_f[, keep_cells, drop = FALSE]
    merged_obj <- subset(merged_obj, cells = colnames(expr_f))
    if (ncol(expr_f) < MIN_CELLS * 2) return(NULL)
  }

  gv <- factor(merged_obj$group, levels = c(group1, group2))
  design <- model.matrix(~ 0 + gv)
  colnames(design) <- levels(gv)
  
  if (qr(design)$rank < ncol(design)) {
    cat("    ⚠️ 设计矩阵不满秩，跳过\n")
    return(NULL)
  }
  
  cont <- makeContrasts(contrasts = paste(group1, "-", group2), levels = design)

  dge <- tryCatch({
    DGEList(counts = expr_f)
  }, error = function(e) { return(NULL) })
  
  if (is.null(dge)) return(NULL)
  
  dge <- tryCatch({
    calcNormFactors(dge)
  }, error = function(e) { return(NULL) })
  
  if (is.null(dge)) return(NULL)
  
  v <- tryCatch({
    voom(dge, design, plot = FALSE)
  }, error = function(e) { return(NULL) })
  
  if (is.null(v)) return(NULL)

  fit <- lmFit(v, design)
  fit2 <- contrasts.fit(fit, cont)
  fit2 <- eBayes(fit2)

  res <- topTable(fit2, number = Inf, adjust.method = "BH")
  if (is.null(res) || nrow(res) == 0) return(NULL)
  
  res$gene <- rownames(res)
  
  # 逐级尝试阈值
  result <- NULL
  used_threshold <- NULL
  
  for (lvl_idx in 1:length(threshold_levels)) {
    lvl <- threshold_levels[[lvl_idx]]
    p_col <- if (lvl$use_adj) "adj.P.Val" else "P.Value"
    
    res$sig <- ifelse(abs(res$logFC) > lvl$logFC & res[[p_col]] < lvl$p_cut,
                      ifelse(res$logFC > 0, "Up", "Down"), "NS")
    
    n_up <- sum(res$sig == "Up")
    n_down <- sum(res$sig == "Down")
    
    if (n_up + n_down >= 10) {
      result <- res
      used_threshold <- lvl
      cat(sprintf("  ✓ %s | %s | 阈值Level %d (%s) | Up=%d Down=%d Total=%d\n", 
                  comp_name, cell_type_name, lvl_idx, lvl$name, n_up, n_down, n_up + n_down))
      break
    }
  }
  
  if (is.null(result)) {
    # 使用最宽松的阈值
    lvl <- threshold_levels[[length(threshold_levels)]]
    p_col <- if (lvl$use_adj) "adj.P.Val" else "P.Value"
    res$sig <- ifelse(abs(res$logFC) > lvl$logFC & res[[p_col]] < lvl$p_cut,
                      ifelse(res$logFC > 0, "Up", "Down"), "NS")
    result <- res
    used_threshold <- lvl
    n_up <- sum(result$sig == "Up")
    n_down <- sum(result$sig == "Down")
    cat(sprintf("  ✓ %s | %s | 阈值Level %d (%s) | Up=%d Down=%d Total=%d\n", 
                comp_name, cell_type_name, length(threshold_levels), lvl$name, n_up, n_down, n_up + n_down))
  }
  
  # 记录使用的阈值
  attr(result, "threshold") <- used_threshold
  
  return(result)
}

# ==================== 4. 富集分析函数 ====================
run_enrichment <- function(deg_df, comp_name, output_dir, cell_type_name) {
  if (is.null(deg_df) || nrow(deg_df) == 0) return(NULL)
  
  sig_genes <- deg_df %>% filter(sig != "NS") %>% pull(gene)
  if (length(sig_genes) < 10) {
    cat("    ⚠️ 显著基因数 < 10，跳过富集分析\n")
    return(NULL)
  }
  
  gene_df <- tryCatch({
    bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  }, error = function(e) { 
    cat("    ⚠️ bitr转换失败:", e$message, "\n")
    return(NULL) 
  })
  
  if (is.null(gene_df) || nrow(gene_df) < 10) return(NULL)
  
  entrez_ids <- unique(gene_df$ENTREZID)
  
  # KEGG富集
  kegg_res <- tryCatch({
    enrichKEGG(gene = entrez_ids, organism = "mmu", pvalueCutoff = 0.05)
  }, error = function(e) return(NULL))
  
  if (!is.null(kegg_res) && nrow(kegg_res@result) > 0) {
    write.csv(as.data.frame(kegg_res), 
              file.path(output_dir, paste0(comp_name, "_", cell_type_name, "_KEGG.csv")), 
              row.names = FALSE)
    cat("    ✓ KEGG富集完成，通路数:", nrow(kegg_res@result), "\n")
  }
  
  # GO富集
  go_res <- tryCatch({
    enrichGO(gene = entrez_ids, OrgDb = org.Mm.eg.db, ont = "BP", 
             pvalueCutoff = 0.05, readable = TRUE)
  }, error = function(e) return(NULL))
  
  if (!is.null(go_res) && nrow(go_res@result) > 0) {
    write.csv(as.data.frame(go_res), 
              file.path(output_dir, paste0(comp_name, "_", cell_type_name, "_GO_BP.csv")), 
              row.names = FALSE)
    cat("    ✓ GO富集完成，条目数:", nrow(go_res@result), "\n")
  }
  
  return(list(kegg = kegg_res, go = go_res))
}

# ==================== 5. 可视化函数 ====================
plot_volcano <- function(res, comp_name, cell_type_name, save_dir, threshold_info) {
  if (is.null(res) || nrow(res) == 0) return()
  
  res$log10p <- -log10(res$adj.P.Val + 1e-300)
  top_genes <- bind_rows(
    res %>% filter(sig == "Up") %>% arrange(adj.P.Val) %>% head(10),
    res %>% filter(sig == "Down") %>% arrange(adj.P.Val) %>% head(10)
  )

  p <- ggplot(res, aes(x = logFC, y = log10p, color = sig)) +
    geom_point(alpha = 0.5, size = 0.7) +
    scale_color_manual(values = c("Up"="#E41A1C","Down"="#377EB8","NS"="grey70")) +
    geom_vline(xintercept = c(-threshold_info$logFC, threshold_info$logFC), 
               linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(threshold_info$p_cut), 
               linetype = "dashed", color = "black") +
    geom_text_repel(data = top_genes, aes(label = gene), size = 2.5, max.overlaps = 25) +
    labs(title = paste0(comp_name, " (", cell_type_name, ")\n阈值: |log2FC|>", threshold_info$logFC,
                        ifelse(threshold_info$use_adj, ", adj.P.Val<0.05", ", P.Value<0.05")), 
         x = "Log2 Fold Change", y = "-Log10 adj.P.Val") +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5, size = 9))

  ggsave(file.path(save_dir, paste0(comp_name, "_", cell_type_name, "_volcano.pdf")), 
         p, width = 8, height = 6, dpi = 300)
  cat("    ✓ 火山图已保存\n")
}

plot_heatmap <- function(res, merged_obj, comp_name, cell_type_name, save_dir) {
  top_genes <- res %>% filter(sig != "NS") %>% arrange(desc(abs(logFC))) %>% head(30)
  if (nrow(top_genes) < 2) return()

  expr <- GetAssayData(merged_obj, layer = "data")
  expr_f <- as.matrix(expr[top_genes$gene, ])
  expr_z <- t(scale(t(expr_f))); expr_z[is.na(expr_z)] <- 0

  anno_df <- data.frame(Group = merged_obj$group, row.names = colnames(expr_z))
  if (ncol(expr_z) > 200) {
    set.seed(42)
    samp <- sample(colnames(expr_z), 200)
    expr_z <- expr_z[, samp]; anno_df <- anno_df[samp, , drop = FALSE]
  }

  pdf(file.path(save_dir, paste0(comp_name, "_", cell_type_name, "_heatmap.pdf")), 
      width = 12, height = 10)
  pheatmap(expr_z, annotation_col = anno_df, show_rownames = TRUE, show_colnames = FALSE,
           cluster_rows = TRUE, cluster_cols = TRUE, fontsize_row = 8,
           main = paste0(comp_name, " (", cell_type_name, ") Top DEGs"),
           color = colorRampPalette(c("#4575B4","white","#D73027"))(100))
  dev.off()
  cat("    ✓ 热图已保存\n")
}

# ==================== 6. 统计细胞数量 ====================
cat("\n========== 2. 统计各样本细胞数量 ==========\n")

cell_count_summary <- data.frame()

for (s in names(seurat_objs)) {
  for (ct in analysis_cell_types) {
    cnt <- count_cells(seurat_objs[[s]], ct$name)
    cell_count_summary <- rbind(cell_count_summary, 
                                 data.frame(Sample = s, CellType = ct$name, Count = cnt))
  }
}

write.csv(cell_count_summary, file.path(output_dir, "threshold_summary", "cell_counts.csv"), row.names = FALSE)
cat("✓ 细胞数量统计已保存\n")

# 打印统计
cat("\n细胞数量统计:\n")
print(cell_count_summary)

# ==================== 7. 肝脏内部比较 ====================
cat("\n\n========== 3. 肝脏内部比较 ==========\n")

liver_deg <- list()
liver_enrich <- list()

# 7.1 肿瘤 vs 空白
cat("\n--- 3.1 肿瘤 vs 空白 ---\n")

for (ct in analysis_cell_types) {
  comp_name <- paste0("liver_tumor_vs_blank_", ct$name)
  cat("\n  >>", comp_name, "\n")
  
  tumor_objs <- list()
  for (t in c("gan56", "liver3", "liver8")) {
    cnt <- count_cells(seurat_objs[[t]], ct$name)
    if (cnt >= MIN_CELLS) {
      tumor_objs[[t]] <- seurat_objs[[t]]
    }
  }
  
  blank_obj <- seurat_objs[["bla_liver"]]
  blank_cnt <- count_cells(blank_obj, ct$name)
  
  if (length(tumor_objs) > 0 && blank_cnt >= MIN_CELLS) {
    merged_tumor <- merge_objects(tumor_objs, rep("Tumor", length(tumor_objs)), ct$name)
    if (!is.null(merged_tumor)) {
      blank_subset <- subset(blank_obj, cells = extract_cells(blank_obj, ct$name))
      blank_subset$group <- "Blank"
      merged <- merge(merged_tumor, blank_subset)
      merged <- JoinLayers(merged)
      
      if (!is.null(merged)) {
        res <- run_limma_diff_adaptive(merged, comp_name, "Tumor", "Blank", ct$name)
        if (!is.null(res)) {
          liver_deg[[comp_name]] <- res
          write.csv(res, file.path(output_dir, "deg_tables", paste0(comp_name, ".csv")), row.names = FALSE)
          plot_volcano(res, comp_name, ct$name, file.path(output_dir, "plots/volcano"), attr(res, "threshold"))
          plot_heatmap(res, merged, comp_name, ct$name, file.path(output_dir, "plots/heatmap"))
          enrich <- run_enrichment(res, comp_name, file.path(output_dir, "liver"), ct$name)
          if (!is.null(enrich)) liver_enrich[[comp_name]] <- enrich
        }
      }
    }
  } else {
    cat("    跳过：肿瘤样本数=", length(tumor_objs), ", 空白细胞数=", blank_cnt, "\n")
  }
}

# 7.2 肿瘤 vs 对照
cat("\n\n--- 3.2 肿瘤 vs 对照 ---\n")

for (ct in analysis_cell_types) {
  comp_name <- paste0("liver_tumor_vs_control_", ct$name)
  cat("\n  >>", comp_name, "\n")
  
  tumor_objs <- list()
  for (t in c("gan56", "liver3", "liver8")) {
    cnt <- count_cells(seurat_objs[[t]], ct$name)
    if (cnt >= MIN_CELLS) {
      tumor_objs[[t]] <- seurat_objs[[t]]
    }
  }
  
  control_obj <- seurat_objs[["con_liver"]]
  control_cnt <- count_cells(control_obj, ct$name)
  
  if (length(tumor_objs) > 0 && control_cnt >= MIN_CELLS) {
    merged_tumor <- merge_objects(tumor_objs, rep("Tumor", length(tumor_objs)), ct$name)
    if (!is.null(merged_tumor)) {
      control_subset <- subset(control_obj, cells = extract_cells(control_obj, ct$name))
      control_subset$group <- "Control"
      merged <- merge(merged_tumor, control_subset)
      merged <- JoinLayers(merged)
      
      if (!is.null(merged)) {
        res <- run_limma_diff_adaptive(merged, comp_name, "Tumor", "Control", ct$name)
        if (!is.null(res)) {
          liver_deg[[comp_name]] <- res
          write.csv(res, file.path(output_dir, "deg_tables", paste0(comp_name, ".csv")), row.names = FALSE)
          plot_volcano(res, comp_name, ct$name, file.path(output_dir, "plots/volcano"), attr(res, "threshold"))
          plot_heatmap(res, merged, comp_name, ct$name, file.path(output_dir, "plots/heatmap"))
          enrich <- run_enrichment(res, comp_name, file.path(output_dir, "liver"), ct$name)
          if (!is.null(enrich)) liver_enrich[[comp_name]] <- enrich
        }
      }
    }
  } else {
    cat("    跳过：肿瘤样本数=", length(tumor_objs), ", 对照细胞数=", control_cnt, "\n")
  }
}

# 7.3 空白 vs 对照
cat("\n\n--- 3.3 空白 vs 对照 ---\n")

for (ct in analysis_cell_types) {
  comp_name <- paste0("liver_blank_vs_control_", ct$name)
  cat("\n  >>", comp_name, "\n")
  
  blank_obj <- seurat_objs[["bla_liver"]]
  control_obj <- seurat_objs[["con_liver"]]
  
  blank_cnt <- count_cells(blank_obj, ct$name)
  control_cnt <- count_cells(control_obj, ct$name)
  
  if (blank_cnt >= MIN_CELLS && control_cnt >= MIN_CELLS) {
    merged <- merge_two_objects(blank_obj, control_obj, ct$name, "Blank", "Control")
    if (!is.null(merged)) {
      res <- run_limma_diff_adaptive(merged, comp_name, "Blank", "Control", ct$name)
      if (!is.null(res)) {
        liver_deg[[comp_name]] <- res
        write.csv(res, file.path(output_dir, "deg_tables", paste0(comp_name, ".csv")), row.names = FALSE)
        plot_volcano(res, comp_name, ct$name, file.path(output_dir, "plots/volcano"), attr(res, "threshold"))
        plot_heatmap(res, merged, comp_name, ct$name, file.path(output_dir, "plots/heatmap"))
        enrich <- run_enrichment(res, comp_name, file.path(output_dir, "liver"), ct$name)
        if (!is.null(enrich)) liver_enrich[[comp_name]] <- enrich
      }
    }
  } else {
    cat("    跳过：空白细胞数=", blank_cnt, ", 对照细胞数=", control_cnt, "\n")
  }
}

# ==================== 8. 脑内部比较 ====================
cat("\n\n========== 4. 脑内部比较 ==========\n")

brain_deg <- list()
brain_enrich <- list()

# 8.1 肿瘤对应脑 vs 脑空白
cat("\n--- 4.1 肿瘤对应脑 vs 脑空白 ---\n")

for (ct in analysis_cell_types) {
  comp_name <- paste0("brain_tumor_vs_blank_", ct$name)
  cat("\n  >>", comp_name, "\n")
  
  brain_objs <- list()
  for (b in c("nao56", "brain3", "brain8")) {
    cnt <- count_cells(seurat_objs[[b]], ct$name)
    if (cnt >= MIN_CELLS) {
      brain_objs[[b]] <- seurat_objs[[b]]
    }
  }
  
  blank_brain_obj <- seurat_objs[["bla_brain"]]
  blank_cnt <- count_cells(blank_brain_obj, ct$name)
  
  if (length(brain_objs) > 0 && blank_cnt >= MIN_CELLS) {
    merged_brain <- merge_objects(brain_objs, rep("TumorBrain", length(brain_objs)), ct$name)
    if (!is.null(merged_brain)) {
      blank_subset <- subset(blank_brain_obj, cells = extract_cells(blank_brain_obj, ct$name))
      blank_subset$group <- "Blank"
      merged <- merge(merged_brain, blank_subset)
      merged <- JoinLayers(merged)
      
      if (!is.null(merged)) {
        res <- run_limma_diff_adaptive(merged, comp_name, "TumorBrain", "Blank", ct$name)
        if (!is.null(res)) {
          brain_deg[[comp_name]] <- res
          write.csv(res, file.path(output_dir, "deg_tables", paste0(comp_name, ".csv")), row.names = FALSE)
          plot_volcano(res, comp_name, ct$name, file.path(output_dir, "plots/volcano"), attr(res, "threshold"))
          plot_heatmap(res, merged, comp_name, ct$name, file.path(output_dir, "plots/heatmap"))
          enrich <- run_enrichment(res, comp_name, file.path(output_dir, "brain"), ct$name)
          if (!is.null(enrich)) brain_enrich[[comp_name]] <- enrich
        }
      }
    }
  } else {
    cat("    跳过：脑肿瘤样本数=", length(brain_objs), ", 脑空白细胞数=", blank_cnt, "\n")
  }
}

# 8.2 肿瘤对应脑 vs 脑对照
cat("\n\n--- 4.2 肿瘤对应脑 vs 脑对照 ---\n")

for (ct in analysis_cell_types) {
  comp_name <- paste0("brain_tumor_vs_control_", ct$name)
  cat("\n  >>", comp_name, "\n")
  
  brain_objs <- list()
  for (b in c("nao56", "brain3", "brain8")) {
    cnt <- count_cells(seurat_objs[[b]], ct$name)
    if (cnt >= MIN_CELLS) {
      brain_objs[[b]] <- seurat_objs[[b]]
    }
  }
  
  control_brain_obj <- seurat_objs[["con_brain"]]
  control_cnt <- count_cells(control_brain_obj, ct$name)
  
  if (length(brain_objs) > 0 && control_cnt >= MIN_CELLS) {
    merged_brain <- merge_objects(brain_objs, rep("TumorBrain", length(brain_objs)), ct$name)
    if (!is.null(merged_brain)) {
      control_subset <- subset(control_brain_obj, cells = extract_cells(control_brain_obj, ct$name))
      control_subset$group <- "Control"
      merged <- merge(merged_brain, control_subset)
      merged <- JoinLayers(merged)
      
      if (!is.null(merged)) {
        res <- run_limma_diff_adaptive(merged, comp_name, "TumorBrain", "Control", ct$name)
        if (!is.null(res)) {
          brain_deg[[comp_name]] <- res
          write.csv(res, file.path(output_dir, "deg_tables", paste0(comp_name, ".csv")), row.names = FALSE)
          plot_volcano(res, comp_name, ct$name, file.path(output_dir, "plots/volcano"), attr(res, "threshold"))
          plot_heatmap(res, merged, comp_name, ct$name, file.path(output_dir, "plots/heatmap"))
          enrich <- run_enrichment(res, comp_name, file.path(output_dir, "brain"), ct$name)
          if (!is.null(enrich)) brain_enrich[[comp_name]] <- enrich
        }
      }
    }
  } else {
    cat("    跳过：脑肿瘤样本数=", length(brain_objs), ", 脑对照细胞数=", control_cnt, "\n")
  }
}

# 8.3 脑空白 vs 脑对照
cat("\n\n--- 4.3 脑空白 vs 脑对照 ---\n")

for (ct in analysis_cell_types) {
  comp_name <- paste0("brain_blank_vs_control_", ct$name)
  cat("\n  >>", comp_name, "\n")
  
  blank_brain_obj <- seurat_objs[["bla_brain"]]
  control_brain_obj <- seurat_objs[["con_brain"]]
  
  blank_cnt <- count_cells(blank_brain_obj, ct$name)
  control_cnt <- count_cells(control_brain_obj, ct$name)
  
  if (blank_cnt >= MIN_CELLS && control_cnt >= MIN_CELLS) {
    merged <- merge_two_objects(blank_brain_obj, control_brain_obj, ct$name, "Blank", "Control")
    if (!is.null(merged)) {
      res <- run_limma_diff_adaptive(merged, comp_name, "Blank", "Control", ct$name)
      if (!is.null(res)) {
        brain_deg[[comp_name]] <- res
        write.csv(res, file.path(output_dir, "deg_tables", paste0(comp_name, ".csv")), row.names = FALSE)
        plot_volcano(res, comp_name, ct$name, file.path(output_dir, "plots/volcano"), attr(res, "threshold"))
        plot_heatmap(res, merged, comp_name, ct$name, file.path(output_dir, "plots/heatmap"))
        enrich <- run_enrichment(res, comp_name, file.path(output_dir, "brain"), ct$name)
        if (!is.null(enrich)) brain_enrich[[comp_name]] <- enrich
      }
    }
  } else {
    cat("    跳过：脑空白细胞数=", blank_cnt, ", 脑对照细胞数=", control_cnt, "\n")
  }
}

# ==================== 9. 跨器官比较 ====================
cat("\n\n========== 5. 跨器官比较 ====================\n")

crosstalk_deg <- list()

for (comp in crosstalk_comparisons) {
  cat("\n---", comp$name, "---\n")
  liver_obj <- seurat_objs[[comp$liver]]
  brain_obj <- seurat_objs[[comp$brain]]
  
  if (is.null(liver_obj) || is.null(brain_obj)) next
  
  for (ct in analysis_cell_types) {
    comp_name <- paste0(comp$name, "_", ct$name)
    cat("\n  >>", comp_name, "\n")
    
    liver_cnt <- count_cells(liver_obj, ct$name)
    brain_cnt <- count_cells(brain_obj, ct$name)
    
    if (liver_cnt >= MIN_CELLS && brain_cnt >= MIN_CELLS) {
      cat("    肝脏细胞数:", liver_cnt, ", 脑细胞数:", brain_cnt, "\n")
      
      merged <- merge_two_objects(liver_obj, brain_obj, ct$name, "Liver", "Brain")
      
      if (!is.null(merged)) {
        res <- run_limma_diff_adaptive(merged, comp_name, "Liver", "Brain", ct$name)
        if (!is.null(res)) {
          crosstalk_deg[[comp_name]] <- res
          write.csv(res, file.path(output_dir, "deg_tables", paste0(comp_name, ".csv")), row.names = FALSE)
          plot_volcano(res, comp_name, ct$name, file.path(output_dir, "plots/volcano"), attr(res, "threshold"))
          plot_heatmap(res, merged, comp_name, ct$name, file.path(output_dir, "plots/heatmap"))
          enrich <- run_enrichment(res, comp_name, file.path(output_dir, "crosstalk"), ct$name)
        }
      }
    } else {
      cat("    跳过：肝脏细胞数=", liver_cnt, ", 脑细胞数=", brain_cnt, "\n")
    }
  }
}

# ==================== 10. 肝脑交互分析 ====================
cat("\n\n========== 6. 肝脑交互分析 ====================\n")

# 收集跨器官比较的KEGG富集结果
crosstalk_kegg_files <- list.files(file.path(output_dir, "crosstalk"), pattern = "*_KEGG.csv", full.names = TRUE)

if (length(crosstalk_kegg_files) > 0) {
  all_crosstalk_pathways <- list()
  
  for (file in crosstalk_kegg_files) {
    df <- read.csv(file)
    comp_name <- gsub("_KEGG.csv", "", basename(file))
    all_crosstalk_pathways[[comp_name]] <- df
  }
  
  # 合并所有跨器官比较的通路
  crosstalk_all <- bind_rows(all_crosstalk_pathways, .id = "Comparison") %>%
    group_by(ID, Description) %>%
    summarise(
      n_comparisons = n(),
      mean_p = mean(p.adjust, na.rm = TRUE),
      min_p = min(p.adjust, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(n_comparisons), mean_p)
  
  cat(sprintf("  总通路数: %d\n", nrow(crosstalk_all)))
  cat(sprintf("  出现在所有3个比较中的通路: %d\n", sum(crosstalk_all$n_comparisons == 3)))
  cat(sprintf("  出现在至少2个比较中的通路: %d\n", sum(crosstalk_all$n_comparisons >= 2)))
  
  if (nrow(crosstalk_all) > 0) {
    write.csv(crosstalk_all, 
              file.path(output_dir, "crosstalk", "crosstalk_shared_pathways.csv"), 
              row.names = FALSE)
    
    # 绘制气泡图
    if (nrow(crosstalk_all) >= 5) {
      plot_df <- crosstalk_all %>% 
        head(30) %>%
        mutate(logp = -log10(mean_p),
               Description = ifelse(nchar(Description) > 50, 
                                    paste0(substr(Description, 1, 47), "..."), 
                                    Description))
      
      p <- ggplot(plot_df, aes(x = reorder(Description, -logp), y = logp, size = n_comparisons)) +
        geom_point(color = "steelblue", alpha = 0.7) +
        coord_flip() +
        labs(title = "肝癌对脑影响 - 共同富集通路",
             x = "", 
             y = "-Log10 P-value (average)", 
             size = "出现次数") +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5, size = 14),
              axis.text.y = element_text(size = 8))
      
      ggsave(file.path(output_dir, "crosstalk", "crosstalk_pathways_bubble.pdf"), 
             p, width = 12, height = max(6, nrow(plot_df)*0.3), dpi = 300)
      cat("  ✓ 共同通路气泡图已保存\n")
    }
    
    # 打印前20条共同通路
    cat("\nTop 20 共同通路:\n")
    print(head(crosstalk_all, 20))
  }
} else {
  cat("  未找到跨器官比较的KEGG富集结果\n")
}

# ==================== 11. 生成报告 ====================
cat("\n\n========== 7. 生成报告 ==========\n")

report_file <- file.path(output_dir, "analysis_report.txt")
sink(report_file)

cat("========================================\n")
cat("肝脑肿瘤完整分析报告\n")
cat("========================================\n")
cat("分析时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("R版本:", version$version.string, "\n")
cat("clusterProfiler版本:", as.character(packageVersion("clusterProfiler")), "\n\n")

cat("【分析参数】\n")
cat("阈值层级:\n")
for (i in 1:length(threshold_levels)) {
  lvl <- threshold_levels[[i]]
  cat(sprintf("  Level %d (%s): |log2FC| > %.1f, %s < 0.05\n", 
              i, lvl$name, lvl$logFC, ifelse(lvl$use_adj, "adj.P.Val", "P.Value")))
}
cat("最小细胞数:", MIN_CELLS, "\n\n")

cat("【分析的细胞类型】\n")
for (ct in analysis_cell_types) {
  cat(sprintf("  - %s: %s\n", ct$name, ct$description))
}
cat("\n")

cat("【成功完成的比较】\n\n")
cat("肝脏内部比较:\n")
for (nm in names(liver_deg)) cat("  ✓", nm, "\n")
cat("\n脑内部比较:\n")
for (nm in names(brain_deg)) cat("  ✓", nm, "\n")
cat("\n跨器官比较:\n")
for (nm in names(crosstalk_deg)) cat("  ✓", nm, "\n")

cat("\n【DEG数量统计】\n\n")
all_deg <- c(liver_deg, brain_deg, crosstalk_deg)
for (nm in names(all_deg)) {
  res <- all_deg[[nm]]
  up <- sum(res$sig == "Up")
  down <- sum(res$sig == "Down")
  thr <- attr(res, "threshold")
  cat(sprintf("  %s: Up=%d, Down=%d, Total=%d (阈值: %s)\n", 
              nm, up, down, up+down, thr$name))
}

sink()

cat("\n========================================\n")
cat("分析完成！\n")
cat("成功比较数:\n")
cat("  肝脏:", length(names(liver_deg)), "\n")
cat("  脑:", length(names(brain_deg)), "\n")
cat("  跨器官:", length(names(crosstalk_deg)), "\n")
cat("输出目录: ", output_dir, "\n")
cat("========================================\n")
