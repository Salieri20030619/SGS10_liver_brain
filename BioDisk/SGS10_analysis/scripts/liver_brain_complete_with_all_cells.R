# ==================== 肝脑肿瘤完整分析脚本（包含所有细胞混合分析） ====================
# 文件名: liver_brain_complete_with_all_cells.R
# 用法: Rscript liver_brain_complete_with_all_cells.R

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
  file.path(output_dir, "cell_count_summary")
)
for (d in dirs) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

# ==================== 参数设置 ====================
LOG2FC_CUTOFF <- 1
ADJ_P_CUTOFF  <- 0.05
MIN_CELLS     <- 10
MIN_PCT       <- 0.1

cat("========================================\n")
cat("肝脑肿瘤完整分析（包含所有细胞混合分析）\n")
cat("========================================\n")
cat("R版本:", version$version.string, "\n")
cat("clusterProfiler版本:", as.character(packageVersion("clusterProfiler")), "\n")
cat("阈值: |log2FC| >", LOG2FC_CUTOFF, " & adj.P.Val <", ADJ_P_CUTOFF, "\n")
cat("最小细胞数:", MIN_CELLS, "\n\n")

# ==================== 定义比较组 ====================

# 1. 肝脏内部比较
liver_comparisons <- list(
  # 肿瘤 vs 空白
  list(name = "tumor_vs_blank", tumor = c("gan56", "liver3", "liver8"), control = "bla_liver", type = "liver"),
  # 肿瘤 vs 对照
  list(name = "tumor_vs_control", tumor = c("gan56", "liver3", "liver8"), control = "con_liver", type = "liver"),
  # 空白 vs 对照
  list(name = "blank_vs_control", group1 = "bla_liver", group2 = "con_liver", type = "liver")
)

# 2. 脑内部比较
brain_comparisons <- list(
  # 肿瘤对应脑 vs 脑空白
  list(name = "brain_tumor_vs_blank", tumor_brain = c("nao56", "brain3", "brain8"), control = "bla_brain", type = "brain"),
  # 肿瘤对应脑 vs 脑对照
  list(name = "brain_tumor_vs_control", tumor_brain = c("nao56", "brain3", "brain8"), control = "con_brain", type = "brain"),
  # 脑空白 vs 脑对照
  list(name = "brain_blank_vs_control", group1 = "bla_brain", group2 = "con_brain", type = "brain")
)

# 3. 跨器官比较（配对）
crosstalk_comparisons <- list(
  list(name = "gan56_vs_nao56", liver = "gan56", brain = "nao56", pair = "pair1"),
  list(name = "liver3_vs_brain3", liver = "liver3", brain = "brain3", pair = "pair2"),
  list(name = "liver8_vs_brain8", liver = "liver8", brain = "brain8", pair = "pair3")
)

# 细胞类型（用于细胞类型特异性分析）
cell_types <- c("Hepatocytes", "Endothelial", "Kupffer", "T_cells", "B_cells", "Stellate")

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
  if (!is.null(cell_type) && "final_celltype" %in% colnames(obj@meta.data)) {
    return(WhichCells(obj, expression = final_celltype == cell_type))
  }
  return(colnames(obj))
}

count_cells <- function(obj, cell_type = NULL) {
  cells <- extract_cells(obj, cell_type)
  return(length(cells))
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

# ==================== 3. 差异分析函数 ====================
run_limma_diff <- function(merged_obj, comp_name, group1, group2, analysis_type) {
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
  res$sig <- ifelse(abs(res$logFC) > LOG2FC_CUTOFF & res$adj.P.Val < ADJ_P_CUTOFF,
               ifelse(res$logFC > 0, "Up", "Down"), "NS")

  n_up <- sum(res$sig == "Up"); n_down <- sum(res$sig == "Down")
  cat(sprintf("  ✓ %s | %s | Up=%d Down=%d Total=%d\n", 
              comp_name, analysis_type, n_up, n_down, n_up + n_down))
  
  return(res)
}

# ==================== 4. 富集分析函数 ====================
run_enrichment <- function(deg_df, comp_name, output_dir, analysis_type) {
  if (is.null(deg_df) || nrow(deg_df) == 0) return(NULL)
  
  sig_genes <- deg_df %>% filter(sig != "NS") %>% pull(gene)
  if (length(sig_genes) < 10) return(NULL)
  
  gene_df <- tryCatch({
    bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  }, error = function(e) { return(NULL) })
  
  if (is.null(gene_df) || nrow(gene_df) < 10) return(NULL)
  
  entrez_ids <- unique(gene_df$ENTREZID)
  
  kegg_res <- tryCatch({
    enrichKEGG(gene = entrez_ids, organism = "mmu", pvalueCutoff = 0.05)
  }, error = function(e) return(NULL))
  
  if (!is.null(kegg_res) && nrow(kegg_res@result) > 0) {
    write.csv(as.data.frame(kegg_res), 
              file.path(output_dir, paste0(comp_name, "_", analysis_type, "_KEGG.csv")), 
              row.names = FALSE)
    cat("    ✓ KEGG富集完成，通路数:", nrow(kegg_res@result), "\n")
  }
  
  go_res <- tryCatch({
    enrichGO(gene = entrez_ids, OrgDb = org.Mm.eg.db, ont = "BP", 
             pvalueCutoff = 0.05, readable = TRUE)
  }, error = function(e) return(NULL))
  
  if (!is.null(go_res) && nrow(go_res@result) > 0) {
    write.csv(as.data.frame(go_res), 
              file.path(output_dir, paste0(comp_name, "_", analysis_type, "_GO_BP.csv")), 
              row.names = FALSE)
    cat("    ✓ GO富集完成，条目数:", nrow(go_res@result), "\n")
  }
  
  return(list(kegg = kegg_res, go = go_res))
}

# ==================== 5. 可视化函数 ====================
plot_volcano <- function(res, comp_name, analysis_type, save_dir) {
  if (is.null(res) || nrow(res) == 0) return()
  
  res$log10p <- -log10(res$adj.P.Val + 1e-300)
  top_genes <- bind_rows(
    res %>% filter(sig == "Up") %>% arrange(adj.P.Val) %>% head(10),
    res %>% filter(sig == "Down") %>% arrange(adj.P.Val) %>% head(10)
  )

  p <- ggplot(res, aes(x = logFC, y = log10p, color = sig)) +
    geom_point(alpha = 0.5, size = 0.7) +
    scale_color_manual(values = c("Up"="#E41A1C","Down"="#377EB8","NS"="grey70")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    geom_text_repel(data = top_genes, aes(label = gene), size = 2.5, max.overlaps = 25) +
    labs(title = paste0(comp_name, " (", analysis_type, ")"), 
         x = "Log2 Fold Change", y = "-Log10 adj.P.Val") +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5, size = 10))

  ggsave(file.path(save_dir, paste0(comp_name, "_", analysis_type, "_volcano.pdf")), 
         p, width = 8, height = 6, dpi = 300)
  cat("    ✓ 火山图已保存\n")
}

plot_heatmap <- function(res, merged_obj, comp_name, analysis_type, save_dir) {
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

  pdf(file.path(save_dir, paste0(comp_name, "_", analysis_type, "_heatmap.pdf")), 
      width = 12, height = 10)
  pheatmap(expr_z, annotation_col = anno_df, show_rownames = TRUE, show_colnames = FALSE,
           cluster_rows = TRUE, cluster_cols = TRUE, fontsize_row = 8,
           main = paste0(comp_name, " (", analysis_type, ") Top DEGs"),
           color = colorRampPalette(c("#4575B4","white","#D73027"))(100))
  dev.off()
  cat("    ✓ 热图已保存\n")
}

# ==================== 6. 统计细胞数量 ====================
cat("\n========== 2. 统计各样本细胞数量 ==========\n")

cell_count_summary <- data.frame()

for (s in names(seurat_objs)) {
  # 所有细胞
  total_cnt <- count_cells(seurat_objs[[s]], NULL)
  cell_count_summary <- rbind(cell_count_summary, 
                               data.frame(Sample = s, CellType = "All_Cells", Count = total_cnt))
  
  # 各细胞类型
  for (ct in cell_types) {
    cnt <- count_cells(seurat_objs[[s]], ct)
    cell_count_summary <- rbind(cell_count_summary, 
                                 data.frame(Sample = s, CellType = ct, Count = cnt))
  }
}

write.csv(cell_count_summary, file.path(output_dir, "cell_count_summary", "cell_counts.csv"), row.names = FALSE)
cat("✓ 细胞数量统计已保存\n")

# 打印统计
cat("\n细胞数量统计（细胞数 >= 10）:\n")
print(cell_count_summary %>% filter(Count >= MIN_CELLS))

# ==================== 7. 肝脏内部比较 ====================
cat("\n\n========== 3. 肝脏内部比较 ==========\n")
cat("分析类型说明:\n")
cat("  - All_Cells: 所有细胞混合分析（整体水平）\n")
cat("  - Cell_Type_Specific: 细胞类型特异性分析（特定细胞类型）\n\n")

liver_deg <- list()
liver_enrich <- list()

# ==================== 3.1 肿瘤 vs 空白 ====================
cat("\n--- 3.1 肿瘤 vs 空白 ---\n")

# 所有细胞混合分析
cat("\n  【所有细胞混合分析】\n")
comp_name_all <- "liver_tumor_vs_blank_All_Cells"
cat("\n  >>", comp_name_all, "\n")

tumor_objs <- list()
for (t in c("gan56", "liver3", "liver8")) {
  tumor_objs[[t]] <- seurat_objs[[t]]
}

blank_cnt <- count_cells(seurat_objs[["bla_liver"]], NULL)
cat("    肿瘤样本细胞数:", paste(sapply(tumor_objs, ncol), collapse=", "), "\n")
cat("    空白样本细胞数:", blank_cnt, "\n")

merged_tumor <- merge_objects(tumor_objs, rep("Tumor", length(tumor_objs)), NULL)

if (!is.null(merged_tumor)) {
  blank_subset <- subset(seurat_objs[["bla_liver"]], cells = colnames(seurat_objs[["bla_liver"]]))
  blank_subset$group <- "Blank"
  merged <- merge(merged_tumor, blank_subset)
  merged <- JoinLayers(merged)
  
  if (!is.null(merged)) {
    res <- run_limma_diff(merged, comp_name_all, "Tumor", "Blank", "All_Cells")
    if (!is.null(res)) {
      liver_deg[[comp_name_all]] <- res
      write.csv(res, file.path(output_dir, "deg_tables", paste0(comp_name_all, ".csv")), row.names = FALSE)
      plot_volcano(res, comp_name_all, "All_Cells", file.path(output_dir, "plots/volcano"))
      plot_heatmap(res, merged, comp_name_all, "All_Cells", file.path(output_dir, "plots/heatmap"))
      enrich <- run_enrichment(res, comp_name_all, file.path(output_dir, "liver"), "All_Cells")
      if (!is.null(enrich)) liver_enrich[[comp_name_all]] <- enrich
    }
  }
}

# 细胞类型特异性分析
cat("\n  【细胞类型特异性分析】\n")
for (ct in cell_types) {
  comp_name <- paste0("liver_tumor_vs_blank_", ct)
  cat("\n  >>", comp_name, "\n")
  
  tumor_objs <- list()
  tumor_counts <- c()
  for (t in c("gan56", "liver3", "liver8")) {
    cnt <- count_cells(seurat_objs[[t]], ct)
    if (cnt >= MIN_CELLS) {
      tumor_objs[[t]] <- seurat_objs[[t]]
      tumor_counts <- c(tumor_counts, cnt)
    } else {
      cat("    ", t, "细胞数不足 (", cnt, ")，跳过\n")
    }
  }
  
  blank_cnt <- count_cells(seurat_objs[["bla_liver"]], ct)
  
  if (length(tumor_objs) > 0 && blank_cnt >= MIN_CELLS) {
    cat("    肿瘤样本细胞数:", paste(tumor_counts, collapse=", "), "\n")
    cat("    空白样本细胞数:", blank_cnt, "\n")
    
    merged_tumor <- merge_objects(tumor_objs, rep("Tumor", length(tumor_objs)), ct)
    
    if (!is.null(merged_tumor)) {
      blank_subset <- subset(seurat_objs[["bla_liver"]], cells = extract_cells(seurat_objs[["bla_liver"]], ct))
      blank_subset$group <- "Blank"
      merged <- merge(merged_tumor, blank_subset)
      merged <- JoinLayers(merged)
      
      if (!is.null(merged)) {
        res <- run_limma_diff(merged, comp_name, "Tumor", "Blank", "Cell_Type_Specific")
        if (!is.null(res)) {
          liver_deg[[comp_name]] <- res
          write.csv(res, file.path(output_dir, "deg_tables", paste0(comp_name, ".csv")), row.names = FALSE)
          plot_volcano(res, comp_name, "Cell_Type_Specific", file.path(output_dir, "plots/volcano"))
          plot_heatmap(res, merged, comp_name, "Cell_Type_Specific", file.path(output_dir, "plots/heatmap"))
          enrich <- run_enrichment(res, comp_name, file.path(output_dir, "liver"), "Cell_Type_Specific")
          if (!is.null(enrich)) liver_enrich[[comp_name]] <- enrich
        }
      }
    }
  } else {
    cat("    跳过：肿瘤样本数=", length(tumor_objs), ", 空白细胞数=", blank_cnt, "\n")
  }
}

# ==================== 3.2 肿瘤 vs 对照 ====================
cat("\n\n--- 3.2 肿瘤 vs 对照 ---\n")

# 所有细胞混合分析
cat("\n  【所有细胞混合分析】\n")
comp_name_all <- "liver_tumor_vs_control_All_Cells"
cat("\n  >>", comp_name_all, "\n")

tumor_objs <- list()
for (t in c("gan56", "liver3", "liver8")) {
  tumor_objs[[t]] <- seurat_objs[[t]]
}

control_cnt <- count_cells(seurat_objs[["con_liver"]], NULL)
cat("    肿瘤样本细胞数:", paste(sapply(tumor_objs, ncol), collapse=", "), "\n")
cat("    对照样本细胞数:", control_cnt, "\n")

merged_tumor <- merge_objects(tumor_objs, rep("Tumor", length(tumor_objs)), NULL)

if (!is.null(merged_tumor)) {
  control_subset <- subset(seurat_objs[["con_liver"]], cells = colnames(seurat_objs[["con_liver"]]))
  control_subset$group <- "Control"
  merged <- merge(merged_tumor, control_subset)
  merged <- JoinLayers(merged)
  
  if (!is.null(merged)) {
    res <- run_limma_diff(merged, comp_name_all, "Tumor", "Control", "All_Cells")
    if (!is.null(res)) {
      liver_deg[[comp_name_all]] <- res
      write.csv(res, file.path(output_dir, "deg_tables", paste0(comp_name_all, ".csv")), row.names = FALSE)
      plot_volcano(res, comp_name_all, "All_Cells", file.path(output_dir, "plots/volcano"))
      plot_heatmap(res, merged, comp_name_all, "All_Cells", file.path(output_dir, "plots/heatmap"))
      enrich <- run_enrichment(res, comp_name_all, file.path(output_dir, "liver"), "All_Cells")
      if (!is.null(enrich)) liver_enrich[[comp_name_all]] <- enrich
    }
  }
}

# 细胞类型特异性分析
cat("\n  【细胞类型特异性分析】\n")
for (ct in cell_types) {
  comp_name <- paste0("liver_tumor_vs_control_", ct)
  cat("\n  >>", comp_name, "\n")
  
  tumor_objs <- list()
  tumor_counts <- c()
  for (t in c("gan56", "liver3", "liver8")) {
    cnt <- count_cells(seurat_objs[[t]], ct)
    if (cnt >= MIN_CELLS) {
      tumor_objs[[t]] <- seurat_objs[[t]]
      tumor_counts <- c(tumor_counts, cnt)
    } else {
      cat("    ", t, "细胞数不足 (", cnt, ")，跳过\n")
    }
  }
  
  control_cnt <- count_cells(seurat_objs[["con_liver"]], ct)
  
  if (length(tumor_objs) > 0 && control_cnt >= MIN_CELLS) {
    cat("    肿瘤样本细胞数:", paste(tumor_counts, collapse=", "), "\n")
    cat("    对照样本细胞数:", control_cnt, "\n")
    
    merged_tumor <- merge_objects(tumor_objs, rep("Tumor", length(tumor_objs)), ct)
    
    if (!is.null(merged_tumor)) {
      control_subset <- subset(seurat_objs[["con_liver"]], cells = extract_cells(seurat_objs[["con_liver"]], ct))
      control_subset$group <- "Control"
      merged <- merge(merged_tumor, control_subset)
      merged <- JoinLayers(merged)
      
      if (!is.null(merged)) {
        res <- run_limma_diff(merged, comp_name, "Tumor", "Control", "Cell_Type_Specific")
        if (!is.null(res)) {
          liver_deg[[comp_name]] <- res
          write.csv(res, file.path(output_dir, "deg_tables", paste0(comp_name, ".csv")), row.names = FALSE)
          plot_volcano(res, comp_name, "Cell_Type_Specific", file.path(output_dir, "plots/volcano"))
          plot_heatmap(res, merged, comp_name, "Cell_Type_Specific", file.path(output_dir, "plots/heatmap"))
          enrich <- run_enrichment(res, comp_name, file.path(output_dir, "liver"), "Cell_Type_Specific")
          if (!is.null(enrich)) liver_enrich[[comp_name]] <- enrich
        }
      }
    }
  } else {
    cat("    跳过：肿瘤样本数=", length(tumor_objs), ", 对照细胞数=", control_cnt, "\n")
  }
}

# ==================== 3.3 空白 vs 对照 ====================
cat("\n\n--- 3.3 空白 vs 对照 ---\n")

# 所有细胞混合分析
cat("\n  【所有细胞混合分析】\n")
comp_name_all <- "liver_blank_vs_control_All_Cells"
cat("\n  >>", comp_name_all, "\n")

blank_cnt <- count_cells(seurat_objs[["bla_liver"]], NULL)
control_cnt <- count_cells(seurat_objs[["con_liver"]], NULL)
cat("    空白样本细胞数:", blank_cnt, "\n")
cat("    对照样本细胞数:", control_cnt, "\n")

if (blank_cnt >= MIN_CELLS && control_cnt >= MIN_CELLS) {
  blank_subset <- subset(seurat_objs[["bla_liver"]], cells = colnames(seurat_objs[["bla_liver"]]))
  blank_subset$group <- "Blank"
  control_subset <- subset(seurat_objs[["con_liver"]], cells = colnames(seurat_objs[["con_liver"]]))
  control_subset$group <- "Control"
  merged <- merge(blank_subset, control_subset)
  merged <- JoinLayers(merged)
  
  if (!is.null(merged)) {
    res <- run_limma_diff(merged, comp_name_all, "Blank", "Control", "All_Cells")
    if (!is.null(res)) {
      liver_deg[[comp_name_all]] <- res
      write.csv(res, file.path(output_dir, "deg_tables", paste0(comp_name_all, ".csv")), row.names = FALSE)
      plot_volcano(res, comp_name_all, "All_Cells", file.path(output_dir, "plots/volcano"))
      plot_heatmap(res, merged, comp_name_all, "All_Cells", file.path(output_dir, "plots/heatmap"))
      enrich <- run_enrichment(res, comp_name_all, file.path(output_dir, "liver"), "All_Cells")
      if (!is.null(enrich)) liver_enrich[[comp_name_all]] <- enrich
    }
  }
}

# 细胞类型特异性分析
cat("\n  【细胞类型特异性分析】\n")
for (ct in cell_types) {
  comp_name <- paste0("liver_blank_vs_control_", ct)
  cat("\n  >>", comp_name, "\n")
  
  blank_cnt <- count_cells(seurat_objs[["bla_liver"]], ct)
  control_cnt <- count_cells(seurat_objs[["con_liver"]], ct)
  
  if (blank_cnt >= MIN_CELLS && control_cnt >= MIN_CELLS) {
    cat("    空白样本细胞数:", blank_cnt, "\n")
    cat("    对照样本细胞数:", control_cnt, "\n")
    
    merged <- merge_two_objects(seurat_objs[["bla_liver"]], seurat_objs[["con_liver"]], 
                                 ct, "Blank", "Control")
    
    if (!is.null(merged)) {
      res <- run_limma_diff(merged, comp_name, "Blank", "Control", "Cell_Type_Specific")
      if (!is.null(res)) {
        liver_deg[[comp_name]] <- res
        write.csv(res, file.path(output_dir, "deg_tables", paste0(comp_name, ".csv")), row.names = FALSE)
        plot_volcano(res, comp_name, "Cell_Type_Specific", file.path(output_dir, "plots/volcano"))
        plot_heatmap(res, merged, comp_name, "Cell_Type_Specific", file.path(output_dir, "plots/heatmap"))
        enrich <- run_enrichment(res, comp_name, file.path(output_dir, "liver"), "Cell_Type_Specific")
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

# ==================== 4.1 肿瘤对应脑 vs 脑空白 ====================
cat("\n--- 4.1 肿瘤对应脑 vs 脑空白 ---\n")

# 所有细胞混合分析
cat("\n  【所有细胞混合分析】\n")
comp_name_all <- "brain_tumor_vs_blank_All_Cells"
cat("\n  >>", comp_name_all, "\n")

brain_objs <- list()
for (b in c("nao56", "brain3", "brain8")) {
  brain_objs[[b]] <- seurat_objs[[b]]
}

blank_cnt <- count_cells(seurat_objs[["bla_brain"]], NULL)
cat("    脑肿瘤样本细胞数:", paste(sapply(brain_objs, ncol), collapse=", "), "\n")
cat("    脑空白样本细胞数:", blank_cnt, "\n")

merged_brain <- merge_objects(brain_objs, rep("TumorBrain", length(brain_objs)), NULL)

if (!is.null(merged_brain)) {
  blank_subset <- subset(seurat_objs[["bla_brain"]], cells = colnames(seurat_objs[["bla_brain"]]))
  blank_subset$group <- "Blank"
  merged <- merge(merged_brain, blank_subset)
  merged <- JoinLayers(merged)
  
  if (!is.null(merged)) {
    res <- run_limma_diff(merged, comp_name_all, "TumorBrain", "Blank", "All_Cells")
    if (!is.null(res)) {
      brain_deg[[comp_name_all]] <- res
      write.csv(res, file.path(output_dir, "deg_tables", paste0(comp_name_all, ".csv")), row.names = FALSE)
      plot_volcano(res, comp_name_all, "All_Cells", file.path(output_dir, "plots/volcano"))
      plot_heatmap(res, merged, comp_name_all, "All_Cells", file.path(output_dir, "plots/heatmap"))
      enrich <- run_enrichment(res, comp_name_all, file.path(output_dir, "brain"), "All_Cells")
      if (!is.null(enrich)) brain_enrich[[comp_name_all]] <- enrich
    }
  }
}

# 细胞类型特异性分析
cat("\n  【细胞类型特异性分析】\n")
for (ct in cell_types) {
  comp_name <- paste0("brain_tumor_vs_blank_", ct)
  cat("\n  >>", comp_name, "\n")
  
  brain_objs <- list()
  brain_counts <- c()
  for (b in c("nao56", "brain3", "brain8")) {
    cnt <- count_cells(seurat_objs[[b]], ct)
    if (cnt >= MIN_CELLS) {
      brain_objs[[b]] <- seurat_objs[[b]]
      brain_counts <- c(brain_counts, cnt)
    } else {
      cat("    ", b, "细胞数不足 (", cnt, ")，跳过\n")
    }
  }
  
  blank_cnt <- count_cells(seurat_objs[["bla_brain"]], ct)
  
  if (length(brain_objs) > 0 && blank_cnt >= MIN_CELLS) {
    cat("    脑肿瘤样本细胞数:", paste(brain_counts, collapse=", "), "\n")
    cat("    脑空白样本细胞数:", blank_cnt, "\n")
    
    merged_brain <- merge_objects(brain_objs, rep("TumorBrain", length(brain_objs)), ct)
    
    if (!is.null(merged_brain)) {
      blank_subset <- subset(seurat_objs[["bla_brain"]], cells = extract_cells(seurat_objs[["bla_brain"]], ct))
      blank_subset$group <- "Blank"
      merged <- merge(merged_brain, blank_subset)
      merged <- JoinLayers(merged)
      
      if (!is.null(merged)) {
        res <- run_limma_diff(merged, comp_name, "TumorBrain", "Blank", "Cell_Type_Specific")
        if (!is.null(res)) {
          brain_deg[[comp_name]] <- res
          write.csv(res, file.path(output_dir, "deg_tables", paste0(comp_name, ".csv")), row.names = FALSE)
          plot_volcano(res, comp_name, "Cell_Type_Specific", file.path(output_dir, "plots/volcano"))
          plot_heatmap(res, merged, comp_name, "Cell_Type_Specific", file.path(output_dir, "plots/heatmap"))
          enrich <- run_enrichment(res, comp_name, file.path(output_dir, "brain"), "Cell_Type_Specific")
          if (!is.null(enrich)) brain_enrich[[comp_name]] <- enrich
        }
      }
    }
  } else {
    cat("    跳过：脑肿瘤样本数=", length(brain_objs), ", 脑空白细胞数=", blank_cnt, "\n")
  }
}

# ==================== 4.2 肿瘤对应脑 vs 脑对照 ====================
cat("\n\n--- 4.2 肿瘤对应脑 vs 脑对照 ---\n")

# 所有细胞混合分析
cat("\n  【所有细胞混合分析】\n")
comp_name_all <- "brain_tumor_vs_control_All_Cells"
cat("\n  >>", comp_name_all, "\n")

brain_objs <- list()
for (b in c("nao56", "brain3", "brain8")) {
  brain_objs[[b]] <- seurat_objs[[b]]
}

control_cnt <- count_cells(seurat_objs[["con_brain"]], NULL)
cat("    脑肿瘤样本细胞数:", paste(sapply(brain_objs, ncol), collapse=", "), "\n")
cat("    脑对照样本细胞数:", control_cnt, "\n")

merged_brain <- merge_objects(brain_objs, rep("TumorBrain", length(brain_objs)), NULL)

if (!is.null(merged_brain)) {
  control_subset <- subset(seurat_objs[["con_brain"]], cells = colnames(seurat_objs[["con_brain"]]))
  control_subset$group <- "Control"
  merged <- merge(merged_brain, control_subset)
  merged <- JoinLayers(merged)
  
  if (!is.null(merged)) {
    res <- run_limma_diff(merged, comp_name_all, "TumorBrain", "Control", "All_Cells")
    if (!is.null(res)) {
      brain_deg[[comp_name_all]] <- res
      write.csv(res, file.path(output_dir, "deg_tables", paste0(comp_name_all, ".csv")), row.names = FALSE)
      plot_volcano(res, comp_name_all, "All_Cells", file.path(output_dir, "plots/volcano"))
      plot_heatmap(res, merged, comp_name_all, "All_Cells", file.path(output_dir, "plots/heatmap"))
      enrich <- run_enrichment(res, comp_name_all, file.path(output_dir, "brain"), "All_Cells")
      if (!is.null(enrich)) brain_enrich[[comp_name_all]] <- enrich
    }
  }
}

# 细胞类型特异性分析
cat("\n  【细胞类型特异性分析】\n")
for (ct in cell_types) {
  comp_name <- paste0("brain_tumor_vs_control_", ct)
  cat("\n  >>", comp_name, "\n")
  
  brain_objs <- list()
  brain_counts <- c()
  for (b in c("nao56", "brain3", "brain8")) {
    cnt <- count_cells(seurat_objs[[b]], ct)
    if (cnt >= MIN_CELLS) {
      brain_objs[[b]] <- seurat_objs[[b]]
      brain_counts <- c(brain_counts, cnt)
    } else {
      cat("    ", b, "细胞数不足 (", cnt, ")，跳过\n")
    }
  }
  
  control_cnt <- count_cells(seurat_objs[["con_brain"]], ct)
  
  if (length(brain_objs) > 0 && control_cnt >= MIN_CELLS) {
    cat("    脑肿瘤样本细胞数:", paste(brain_counts, collapse=", "), "\n")
    cat("    脑对照样本细胞数:", control_cnt, "\n")
    
    merged_brain <- merge_objects(brain_objs, rep("TumorBrain", length(brain_objs)), ct)
    
    if (!is.null(merged_brain)) {
      control_subset <- subset(seurat_objs[["con_brain"]], cells = extract_cells(seurat_objs[["con_brain"]], ct))
      control_subset$group <- "Control"
      merged <- merge(merged_brain, control_subset)
      merged <- JoinLayers(merged)
      
      if (!is.null(merged)) {
        res <- run_limma_diff(merged, comp_name, "TumorBrain", "Control", "Cell_Type_Specific")
        if (!is.null(res)) {
          brain_deg[[comp_name]] <- res
          write.csv(res, file.path(output_dir, "deg_tables", paste0(comp_name, ".csv")), row.names = FALSE)
          plot_volcano(res, comp_name, "Cell_Type_Specific", file.path(output_dir, "plots/volcano"))
          plot_heatmap(res, merged, comp_name, "Cell_Type_Specific", file.path(output_dir, "plots/heatmap"))
          enrich <- run_enrichment(res, comp_name, file.path(output_dir, "brain"), "Cell_Type_Specific")
          if (!is.null(enrich)) brain_enrich[[comp_name]] <- enrich
        }
      }
    }
  } else {
    cat("    跳过：脑肿瘤样本数=", length(brain_objs), ", 脑对照细胞数=", control_cnt, "\n")
  }
}

# ==================== 4.3 脑空白 vs 脑对照 ====================
cat("\n\n--- 4.3 脑空白 vs 脑对照 ---\n")

# 所有细胞混合分析
cat("\n  【所有细胞混合分析】\n")
comp_name_all <- "brain_blank_vs_control_All_Cells"
cat("\n  >>", comp_name_all, "\n")

blank_cnt <- count_cells(seurat_objs[["bla_brain"]], NULL)
control_cnt <- count_cells(seurat_objs[["con_brain"]], NULL)
cat("    脑空白样本细胞数:", blank_cnt, "\n")
cat("    脑对照样本细胞数:", control_cnt, "\n")

if (blank_cnt >= MIN_CELLS && control_cnt >= MIN_CELLS) {
  blank_subset <- subset(seurat_objs[["bla_brain"]], cells = colnames(seurat_objs[["bla_brain"]]))
  blank_subset$group <- "Blank"
  control_subset <- subset(seurat_objs[["con_brain"]], cells = colnames(seurat_objs[["con_brain"]]))
  control_subset$group <- "Control"
  merged <- merge(blank_subset, control_subset)
  merged <- JoinLayers(merged)
  
  if (!is.null(merged)) {
    res <- run_limma_diff(merged, comp_name_all, "Blank", "Control", "All_Cells")
    if (!is.null(res)) {
      brain_deg[[comp_name_all]] <- res
      write.csv(res, file.path(output_dir, "deg_tables", paste0(comp_name_all, ".csv")), row.names = FALSE)
      plot_volcano(res, comp_name_all, "All_Cells", file.path(output_dir, "plots/volcano"))
      plot_heatmap(res, merged, comp_name_all, "All_Cells", file.path(output_dir, "plots/heatmap"))
      enrich <- run_enrichment(res, comp_name_all, file.path(output_dir, "brain"), "All_Cells")
      if (!is.null(enrich)) brain_enrich[[comp_name_all]] <- enrich
    }
  }
}

# 细胞类型特异性分析
cat("\n  【细胞类型特异性分析】\n")
for (ct in cell_types) {
  comp_name <- paste0("brain_blank_vs_control_", ct)
  cat("\n  >>", comp_name, "\n")
  
  blank_cnt <- count_cells(seurat_objs[["bla_brain"]], ct)
  control_cnt <- count_cells(seurat_objs[["con_brain"]], ct)
  
  if (blank_cnt >= MIN_CELLS && control_cnt >= MIN_CELLS) {
    cat("    脑空白样本细胞数:", blank_cnt, "\n")
    cat("    脑对照样本细胞数:", control_cnt, "\n")
    
    merged <- merge_two_objects(seurat_objs[["bla_brain"]], seurat_objs[["con_brain"]], 
                                 ct, "Blank", "Control")
    
    if (!is.null(merged)) {
      res <- run_limma_diff(merged, comp_name, "Blank", "Control", "Cell_Type_Specific")
      if (!is.null(res)) {
        brain_deg[[comp_name]] <- res
        write.csv(res, file.path(output_dir, "deg_tables", paste0(comp_name, ".csv")), row.names = FALSE)
        plot_volcano(res, comp_name, "Cell_Type_Specific", file.path(output_dir, "plots/volcano"))
        plot_heatmap(res, merged, comp_name, "Cell_Type_Specific", file.path(output_dir, "plots/heatmap"))
        enrich <- run_enrichment(res, comp_name, file.path(output_dir, "brain"), "Cell_Type_Specific")
        if (!is.null(enrich)) brain_enrich[[comp_name]] <- enrich
      }
    }
  } else {
    cat("    跳过：脑空白细胞数=", blank_cnt, ", 脑对照细胞数=", control_cnt, "\n")
  }
}

# ==================== 9. 跨器官比较 ====================
cat("\n\n========== 5. 跨器官比较 ==========\n")

crosstalk_deg <- list()

for (comp in crosstalk_comparisons) {
  cat("\n---", comp$name, "---\n")
  liver_obj <- seurat_objs[[comp$liver]]
  brain_obj <- seurat_objs[[comp$brain]]
  
  if (is.null(liver_obj) || is.null(brain_obj)) next
  
  # 所有细胞混合分析
  cat("\n  【所有细胞混合分析】\n")
  comp_name_all <- paste0(comp$name, "_All_Cells")
  cat("\n  >>", comp_name_all, "\n")
  
  liver_cnt <- count_cells(liver_obj, NULL)
  brain_cnt <- count_cells(brain_obj, NULL)
  cat("    肝脏细胞数:", liver_cnt, "\n")
  cat("    脑细胞数:", brain_cnt, "\n")
  
  if (liver_cnt >= MIN_CELLS && brain_cnt >= MIN_CELLS) {
    liver_subset <- subset(liver_obj, cells = colnames(liver_obj))
    liver_subset$group <- "Liver"
    brain_subset <- subset(brain_obj, cells = colnames(brain_obj))
    brain_subset$group <- "Brain"
    merged <- merge(liver_subset, brain_subset)
    merged <- JoinLayers(merged)
    
    if (!is.null(merged)) {
      res <- run_limma_diff(merged, comp_name_all, "Liver", "Brain", "All_Cells")
      if (!is.null(res)) {
        crosstalk_deg[[comp_name_all]] <- res
        write.csv(res, file.path(output_dir, "deg_tables", paste0(comp_name_all, ".csv")), row.names = FALSE)
        plot_volcano(res, comp_name_all, "All_Cells", file.path(output_dir, "plots/volcano"))
        plot_heatmap(res, merged, comp_name_all, "All_Cells", file.path(output_dir, "plots/heatmap"))
        enrich <- run_enrichment(res, comp_name_all, file.path(output_dir, "crosstalk"), "All_Cells")
      }
    }
  }
  
  # 细胞类型特异性分析
  cat("\n  【细胞类型特异性分析】\n")
  for (ct in cell_types) {
    comp_name <- paste0(comp$name, "_", ct)
    cat("\n  >>", comp_name, "\n")
    
    liver_cnt <- count_cells(liver_obj, ct)
    brain_cnt <- count_cells(brain_obj, ct)
    
    if (liver_cnt >= MIN_CELLS && brain_cnt >= MIN_CELLS) {
      cat("    肝脏细胞数:", liver_cnt, "\n")
      cat("    脑细胞数:", brain_cnt, "\n")
      
      merged <- merge_two_objects(liver_obj, brain_obj, ct, "Liver", "Brain")
      
      if (!is.null(merged)) {
        res <- run_limma_diff(merged, comp_name, "Liver", "Brain", "Cell_Type_Specific")
        if (!is.null(res)) {
          crosstalk_deg[[comp_name]] <- res
          write.csv(res, file.path(output_dir, "deg_tables", paste0(comp_name, ".csv")), row.names = FALSE)
          plot_volcano(res, comp_name, "Cell_Type_Specific", file.path(output_dir, "plots/volcano"))
          plot_heatmap(res, merged, comp_name, "Cell_Type_Specific", file.path(output_dir, "plots/heatmap"))
          enrich <- run_enrichment(res, comp_name, file.path(output_dir, "crosstalk"), "Cell_Type_Specific")
        }
      }
    } else {
      cat("    跳过：肝脏细胞数=", liver_cnt, ", 脑细胞数=", brain_cnt, "\n")
    }
  }
}

# ==================== 10. 肝脑交互分析 ====================
cat("\n\n========== 6. 肝脑交互分析 ==========\n")

# 收集肝脏和脑的KEGG通路（包括所有细胞混合分析和细胞类型特异性分析）
liver_pathways <- list()
brain_pathways <- list()

for (name in names(liver_enrich)) {
  if (!is.null(liver_enrich[[name]]$kegg)) {
    df <- as.data.frame(liver_enrich[[name]]$kegg)
    liver_pathways[[name]] <- df
  }
}

for (name in names(brain_enrich)) {
  if (!is.null(brain_enrich[[name]]$kegg)) {
    df <- as.data.frame(brain_enrich[[name]]$kegg)
    brain_pathways[[name]] <- df
  }
}

if (length(liver_pathways) > 0 && length(brain_pathways) > 0) {
  liver_all <- bind_rows(liver_pathways, .id = "Comparison") %>%
    group_by(ID, Description) %>%
    summarise(n_liver = n(), mean_p_liver = mean(p.adjust), .groups = "drop")
  
  brain_all <- bind_rows(brain_pathways, .id = "Comparison") %>%
    group_by(ID, Description) %>%
    summarise(n_brain = n(), mean_p_brain = mean(p.adjust), .groups = "drop")
  
  shared <- inner_join(liver_all, brain_all, by = c("ID", "Description")) %>%
    arrange(desc(n_liver + n_brain))
  
  cat(sprintf("  肝脏显著通路数: %d\n", nrow(liver_all)))
  cat(sprintf("  脑显著通路数: %d\n", nrow(brain_all)))
  cat(sprintf("  共同通路数: %d\n", nrow(shared)))
  
  if (nrow(shared) > 0) {
    write.csv(shared, file.path(output_dir, "crosstalk", "shared_pathways.csv"), row.names = FALSE)
    
    # 绘制共同通路散点图
    if (nrow(shared) >= 3) {
      plot_df <- shared %>% head(30) %>%
        mutate(logp_liver = -log10(mean_p_liver),
               logp_brain = -log10(mean_p_brain))
      
      p <- ggplot(plot_df, aes(x = logp_liver, y = logp_brain, 
                                size = n_liver + n_brain, label = Description)) +
        geom_point(alpha = 0.7, color = "steelblue") +
        geom_text_repel(size = 2.5, max.overlaps = 20) +
        labs(title = "肝脑共同通路", x = "-Log10 P-value (Liver)", 
             y = "-Log10 P-value (Brain)", size = "出现次数") +
        theme_bw()
      
      ggsave(file.path(output_dir, "crosstalk", "shared_pathways_scatter.pdf"), p, width = 10, height = 8, dpi = 300)
      cat("  ✓ 共同通路散点图已保存\n")
    }
  }
} else {
  cat("  富集结果不足，跳过交互分析\n")
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
cat("  差异分析阈值: |log2FC| >", LOG2FC_CUTOFF, " & adj.P.Val <", ADJ_P_CUTOFF, "\n")
cat("  最小细胞数:", MIN_CELLS, "\n\n")

cat("【分析类型说明】\n")
cat("  - All_Cells: 所有细胞混合分析（整体水平）\n")
cat("  - Cell_Type_Specific: 细胞类型特异性分析（特定细胞类型）\n\n")

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
  cat(sprintf("  %s: Up=%d, Down=%d, Total=%d\n", nm, up, down, up+down))
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
