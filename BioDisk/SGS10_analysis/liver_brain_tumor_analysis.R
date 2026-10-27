# ==================== 肝脑肿瘤差异分析及交互作用 ====================
# 文件名: liver_brain_tumor_analysis.R
# 用法: Rscript liver_brain_tumor_analysis.R

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

output_dir <- "liver_brain_analysis"
dirs <- c(
  output_dir,
  file.path(output_dir, "liver"),
  file.path(output_dir, "brain"),
  file.path(output_dir, "crosstalk"),
  file.path(output_dir, "plots/volcano"),
  file.path(output_dir, "plots/heatmap"),
  file.path(output_dir, "plots/enrichment"),
  file.path(output_dir, "deg_tables")
)
for (d in dirs) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

# ==================== 参数设置 ====================
LOG2FC_CUTOFF <- 1
ADJ_P_CUTOFF  <- 0.05
MIN_CELLS     <- 10
MIN_PCT       <- 0.1

cat("========================================\n")
cat("肝脑肿瘤差异分析及交互作用\n")
cat("========================================\n")
cat("阈值: |log2FC| >", LOG2FC_CUTOFF, " & adj.P.Val <", ADJ_P_CUTOFF, "\n\n")

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

# ==================== 2. 辅助函数 ====================
extract_cells <- function(obj, cell_type = NULL) {
  if (!is.null(cell_type) && "final_celltype" %in% colnames(obj@meta.data)) {
    return(WhichCells(obj, expression = final_celltype == cell_type))
  }
  return(colnames(obj))
}

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
  
  # 过滤低质量细胞
  expr <- GetAssayData(merged, layer = "data")
  lib_sizes <- colSums(expr)
  keep_cells <- lib_sizes > 0
  if (sum(!keep_cells) > 0) {
    merged <- subset(merged, cells = colnames(expr)[keep_cells])
  }
  
  return(merged)
}

# ==================== 3. 差异分析函数 ====================
run_limma_diff <- function(merged_obj, comp_name, group1, group2) {
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

  gv <- factor(merged_obj$group, levels = c(group1, group2))
  design <- model.matrix(~ 0 + gv); colnames(design) <- levels(gv)
  cont <- makeContrasts(contrasts = paste(group1, "-", group2), levels = design)

  dge <- DGEList(counts = expr_f)
  dge <- calcNormFactors(dge)
  
  v <- tryCatch({
    voom(dge, design, plot = FALSE)
  }, error = function(e) {
    return(NULL)
  })
  
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
  cat(sprintf("  ✓ %s | Up=%d Down=%d Total=%d\n", comp_name, n_up, n_down, n_up + n_down))
  
  return(res)
}

# ==================== 4. 富集分析 ====================
run_enrichment <- function(deg_df, comp_name, output_dir) {
  if (is.null(deg_df) || nrow(deg_df) == 0) return(NULL)
  
  sig_genes <- deg_df %>% filter(sig != "NS") %>% pull(gene)
  if (length(sig_genes) < 10) return(NULL)
  
  # 转换Entrez ID
  gene_df <- tryCatch({
    bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  }, error = function(e) return(NULL))
  
  if (is.null(gene_df) || nrow(gene_df) < 10) return(NULL)
  
  entrez_ids <- unique(gene_df$ENTREZID)
  
  # KEGG富集
  kegg_res <- tryCatch({
    enrichKEGG(gene = entrez_ids, organism = "mmu", pvalueCutoff = 0.05)
  }, error = function(e) return(NULL))
  
  if (!is.null(kegg_res) && nrow(kegg_res@result) > 0) {
    write.csv(as.data.frame(kegg_res), 
              file.path(output_dir, paste0(comp_name, "_KEGG.csv")), 
              row.names = FALSE)
    
    # 绘制气泡图
    if (nrow(kegg_res@result) >= 3) {
      pdf(file.path(output_dir, paste0(comp_name, "_KEGG_dotplot.pdf")), width = 10, height = 8)
      p <- dotplot(kegg_res, showCategory = 15, title = comp_name)
      print(p)
      dev.off()
    }
  }
  
  # GO富集
  go_res <- tryCatch({
    enrichGO(gene = entrez_ids, OrgDb = org.Mm.eg.db, ont = "BP", 
             pvalueCutoff = 0.05, readable = TRUE)
  }, error = function(e) return(NULL))
  
  if (!is.null(go_res) && nrow(go_res@result) > 0) {
    write.csv(as.data.frame(go_res), 
              file.path(output_dir, paste0(comp_name, "_GO_BP.csv")), 
              row.names = FALSE)
  }
  
  return(list(kegg = kegg_res, go = go_res))
}

# ==================== 5. 火山图和热图 ====================
plot_volcano <- function(res, comp_name, save_dir) {
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
    labs(title = comp_name, x = "Log2 Fold Change", y = "-Log10 adj.P.Val") +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5, size = 10))

  ggsave(file.path(save_dir, paste0(comp_name, "_volcano.pdf")), p, width = 8, height = 6, dpi = 300)
}

plot_heatmap <- function(res, merged_obj, comp_name, save_dir) {
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

  pdf(file.path(save_dir, paste0(comp_name, "_heatmap.pdf")), width = 12, height = 10)
  pheatmap(expr_z, annotation_col = anno_df, show_rownames = TRUE, show_colnames = FALSE,
           cluster_rows = TRUE, cluster_cols = TRUE, fontsize_row = 8,
           main = paste0(comp_name, " Top DEGs"),
           color = colorRampPalette(c("#4575B4","white","#D73027"))(100))
  dev.off()
}

# ==================== 6. 肝脑交互分析 ====================
analyze_crosstalk <- function(liver_results, brain_results, crosstalk_results, output_dir) {
  cat("\n========== 肝脑交互分析 ==========\n")
  
  # 提取所有通路的KEGG结果
  liver_pathways <- list()
  brain_pathways <- list()
  
  for (name in names(liver_results)) {
    if (!is.null(liver_results[[name]]$kegg)) {
      df <- as.data.frame(liver_results[[name]]$kegg)
      liver_pathways[[name]] <- df
    }
  }
  
  for (name in names(brain_results)) {
    if (!is.null(brain_results[[name]]$kegg)) {
      df <- as.data.frame(brain_results[[name]]$kegg)
      brain_pathways[[name]] <- df
    }
  }
  
  if (length(liver_pathways) == 0 || length(brain_pathways) == 0) {
    cat("  ⚠️ 富集结果不足，跳过交互分析\n")
    return(NULL)
  }
  
  # 合并所有肝脏通路的ID
  liver_all <- bind_rows(liver_pathways, .id = "Comparison") %>%
    group_by(ID, Description) %>%
    summarise(n_liver = n(), mean_p_liver = mean(p.adjust), .groups = "drop")
  
  brain_all <- bind_rows(brain_pathways, .id = "Comparison") %>%
    group_by(ID, Description) %>%
    summarise(n_brain = n(), mean_p_brain = mean(p.adjust), .groups = "drop")
  
  # 找出共同通路
  shared <- inner_join(liver_all, brain_all, by = c("ID", "Description")) %>%
    arrange(desc(n_liver + n_brain))
  
  cat(sprintf("  肝脏显著通路数: %d\n", nrow(liver_all)))
  cat(sprintf("  脑显著通路数: %d\n", nrow(brain_all)))
  cat(sprintf("  共同通路数: %d\n", nrow(shared)))
  
  if (nrow(shared) > 0) {
    write.csv(shared, file.path(output_dir, "shared_pathways.csv"), row.names = FALSE)
    
    # 绘制共同通路热图
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
      
      ggsave(file.path(output_dir, "shared_pathways_scatter.pdf"), p, width = 10, height = 8, dpi = 300)
    }
  }
  
  return(shared)
}

# ==================== 7. 主分析流程 ====================
cat("\n========== 开始分析 ==========\n\n")

# 定义分析的细胞类型
cell_types <- c("Hepatocytes", "Endothelial", "Kupffer", "T_cells", "B_cells", "Stellate")

# 存储结果
liver_deg_results <- list()
brain_deg_results <- list()
crosstalk_deg_results <- list()
liver_enrich_results <- list()
brain_enrich_results <- list()

# ==================== 7.1 肝脏内部比较 ====================
cat("\n===== 1. 肝脏内部比较 =====\n")

# 肿瘤 vs 空白（合并多个肿瘤样本）
cat("\n--- 肿瘤 vs 空白 ---\n")
for (ct in cell_types) {
  comp_name <- paste0("liver_tumor_vs_blank_", ct)
  cat("\n  >>", comp_name, "\n")
  
  # 获取肿瘤样本
  tumor_objs <- list()
  for (t in c("gan56", "liver3", "liver8")) {
    if (t %in% names(seurat_objs)) tumor_objs[[t]] <- seurat_objs[[t]]
  }
  
  # 获取空白样本
  blank_obj <- seurat_objs[["bla_liver"]]
  
  if (length(tumor_objs) > 0 && !is.null(blank_obj)) {
    # 合并肿瘤样本
    tumor_group <- rep("Tumor", length(tumor_objs))
    merged_tumor <- merge_objects(tumor_objs, tumor_group, cell_type = ct)
    
    # 提取空白细胞
    blank_cells <- extract_cells(blank_obj, ct)
    if (length(blank_cells) >= MIN_CELLS) {
      blank_subset <- subset(blank_obj, cells = blank_cells)
      blank_subset$group <- "Blank"
      
      # 合并肿瘤和空白
      merged <- merge(merged_tumor, blank_subset)
      merged <- JoinLayers(merged)
      
      if (!is.null(merged)) {
        res <- run_limma_diff(merged, comp_name, "Tumor", "Blank")
        if (!is.null(res)) {
          liver_deg_results[[comp_name]] <- res
          
          # 保存结果
          write.csv(res, file.path(output_dir, "deg_tables", paste0(comp_name, ".csv")), row.names = FALSE)
          plot_volcano(res, comp_name, file.path(output_dir, "plots/volcano"))
          plot_heatmap(res, merged, comp_name, file.path(output_dir, "plots/heatmap"))
          
          # 富集分析
          enrich <- run_enrichment(res, comp_name, file.path(output_dir, "liver"))
          if (!is.null(enrich)) liver_enrich_results[[comp_name]] <- enrich
        }
      }
    }
  }
}

# 肿瘤 vs 对照
cat("\n\n--- 肿瘤 vs 对照 ---\n")
for (ct in cell_types) {
  comp_name <- paste0("liver_tumor_vs_control_", ct)
  cat("\n  >>", comp_name, "\n")
  
  tumor_objs <- list()
  for (t in c("gan56", "liver3", "liver8")) {
    if (t %in% names(seurat_objs)) tumor_objs[[t]] <- seurat_objs[[t]]
  }
  
  control_obj <- seurat_objs[["con_liver"]]
  
  if (length(tumor_objs) > 0 && !is.null(control_obj)) {
    tumor_group <- rep("Tumor", length(tumor_objs))
    merged_tumor <- merge_objects(tumor_objs, tumor_group, cell_type = ct)
    
    control_cells <- extract_cells(control_obj, ct)
    if (length(control_cells) >= MIN_CELLS) {
      control_subset <- subset(control_obj, cells = control_cells)
      control_subset$group <- "Control"
      
      merged <- merge(merged_tumor, control_subset)
      merged <- JoinLayers(merged)
      
      if (!is.null(merged)) {
        res <- run_limma_diff(merged, comp_name, "Tumor", "Control")
        if (!is.null(res)) {
          liver_deg_results[[comp_name]] <- res
          write.csv(res, file.path(output_dir, "deg_tables", paste0(comp_name, ".csv")), row.names = FALSE)
          plot_volcano(res, comp_name, file.path(output_dir, "plots/volcano"))
          plot_heatmap(res, merged, comp_name, file.path(output_dir, "plots/heatmap"))
          enrich <- run_enrichment(res, comp_name, file.path(output_dir, "liver"))
          if (!is.null(enrich)) liver_enrich_results[[comp_name]] <- enrich
        }
      }
    }
  }
}

# 空白 vs 对照
cat("\n\n--- 空白 vs 对照 ---\n")
for (ct in cell_types) {
  comp_name <- paste0("liver_blank_vs_control_", ct)
  cat("\n  >>", comp_name, "\n")
  
  blank_obj <- seurat_objs[["bla_liver"]]
  control_obj <- seurat_objs[["con_liver"]]
  
  if (!is.null(blank_obj) && !is.null(control_obj)) {
    blank_cells <- extract_cells(blank_obj, ct)
    control_cells <- extract_cells(control_obj, ct)
    
    if (length(blank_cells) >= MIN_CELLS && length(control_cells) >= MIN_CELLS) {
      blank_subset <- subset(blank_obj, cells = blank_cells)
      blank_subset$group <- "Blank"
      control_subset <- subset(control_obj, cells = control_cells)
      control_subset$group <- "Control"
      
      merged <- merge(blank_subset, control_subset)
      merged <- JoinLayers(merged)
      
      if (!is.null(merged)) {
        res <- run_limma_diff(merged, comp_name, "Blank", "Control")
        if (!is.null(res)) {
          liver_deg_results[[comp_name]] <- res
          write.csv(res, file.path(output_dir, "deg_tables", paste0(comp_name, ".csv")), row.names = FALSE)
          plot_volcano(res, comp_name, file.path(output_dir, "plots/volcano"))
          plot_heatmap(res, merged, comp_name, file.path(output_dir, "plots/heatmap"))
          enrich <- run_enrichment(res, comp_name, file.path(output_dir, "liver"))
          if (!is.null(enrich)) liver_enrich_results[[comp_name]] <- enrich
        }
      }
    }
  }
}

# ==================== 7.2 脑内部比较 ====================
cat("\n\n===== 2. 脑内部比较 =====\n")

# 肿瘤对应脑 vs 脑空白
cat("\n--- 肿瘤对应脑 vs 脑空白 ---\n")
for (ct in cell_types) {
  comp_name <- paste0("brain_tumor_vs_blank_", ct)
  cat("\n  >>", comp_name, "\n")
  
  brain_objs <- list()
  for (b in c("nao56", "brain3", "brain8")) {
    if (b %in% names(seurat_objs)) brain_objs[[b]] <- seurat_objs[[b]]
  }
  
  blank_brain_obj <- seurat_objs[["bla_brain"]]
  
  if (length(brain_objs) > 0 && !is.null(blank_brain_obj)) {
    brain_group <- rep("TumorBrain", length(brain_objs))
    merged_brain <- merge_objects(brain_objs, brain_group, cell_type = ct)
    
    blank_cells <- extract_cells(blank_brain_obj, ct)
    if (length(blank_cells) >= MIN_CELLS) {
      blank_subset <- subset(blank_brain_obj, cells = blank_cells)
      blank_subset$group <- "Blank"
      
      merged <- merge(merged_brain, blank_subset)
      merged <- JoinLayers(merged)
      
      if (!is.null(merged)) {
        res <- run_limma_diff(merged, comp_name, "TumorBrain", "Blank")
        if (!is.null(res)) {
          brain_deg_results[[comp_name]] <- res
          write.csv(res, file.path(output_dir, "deg_tables", paste0(comp_name, ".csv")), row.names = FALSE)
          plot_volcano(res, comp_name, file.path(output_dir, "plots/volcano"))
          plot_heatmap(res, merged, comp_name, file.path(output_dir, "plots/heatmap"))
          enrich <- run_enrichment(res, comp_name, file.path(output_dir, "brain"))
          if (!is.null(enrich)) brain_enrich_results[[comp_name]] <- enrich
        }
      }
    }
  }
}

# 肿瘤对应脑 vs 脑对照
cat("\n\n--- 肿瘤对应脑 vs 脑对照 ---\n")
for (ct in cell_types) {
  comp_name <- paste0("brain_tumor_vs_control_", ct)
  cat("\n  >>", comp_name, "\n")
  
  brain_objs <- list()
  for (b in c("nao56", "brain3", "brain8")) {
    if (b %in% names(seurat_objs)) brain_objs[[b]] <- seurat_objs[[b]]
  }
  
  control_brain_obj <- seurat_objs[["con_brain"]]
  
  if (length(brain_objs) > 0 && !is.null(control_brain_obj)) {
    brain_group <- rep("TumorBrain", length(brain_objs))
    merged_brain <- merge_objects(brain_objs, brain_group, cell_type = ct)
    
    control_cells <- extract_cells(control_brain_obj, ct)
    if (length(control_cells) >= MIN_CELLS) {
      control_subset <- subset(control_brain_obj, cells = control_cells)
      control_subset$group <- "Control"
      
      merged <- merge(merged_brain, control_subset)
      merged <- JoinLayers(merged)
      
      if (!is.null(merged)) {
        res <- run_limma_diff(merged, comp_name, "TumorBrain", "Control")
        if (!is.null(res)) {
          brain_deg_results[[comp_name]] <- res
          write.csv(res, file.path(output_dir, "deg_tables", paste0(comp_name, ".csv")), row.names = FALSE)
          plot_volcano(res, comp_name, file.path(output_dir, "plots/volcano"))
          plot_heatmap(res, merged, comp_name, file.path(output_dir, "plots/heatmap"))
          enrich <- run_enrichment(res, comp_name, file.path(output_dir, "brain"))
          if (!is.null(enrich)) brain_enrich_results[[comp_name]] <- enrich
        }
      }
    }
  }
}

# 脑空白 vs 脑对照
cat("\n\n--- 脑空白 vs 脑对照 ---\n")
for (ct in cell_types) {
  comp_name <- paste0("brain_blank_vs_control_", ct)
  cat("\n  >>", comp_name, "\n")
  
  blank_obj <- seurat_objs[["bla_brain"]]
  control_obj <- seurat_objs[["con_brain"]]
  
  if (!is.null(blank_obj) && !is.null(control_obj)) {
    blank_cells <- extract_cells(blank_obj, ct)
    control_cells <- extract_cells(control_obj, ct)
    
    if (length(blank_cells) >= MIN_CELLS && length(control_cells) >= MIN_CELLS) {
      blank_subset <- subset(blank_obj, cells = blank_cells)
      blank_subset$group <- "Blank"
      control_subset <- subset(control_obj, cells = control_cells)
      control_subset$group <- "Control"
      
      merged <- merge(blank_subset, control_subset)
      merged <- JoinLayers(merged)
      
      if (!is.null(merged)) {
        res <- run_limma_diff(merged, comp_name, "Blank", "Control")
        if (!is.null(res)) {
          brain_deg_results[[comp_name]] <- res
          write.csv(res, file.path(output_dir, "deg_tables", paste0(comp_name, ".csv")), row.names = FALSE)
          plot_volcano(res, comp_name, file.path(output_dir, "plots/volcano"))
          plot_heatmap(res, merged, comp_name, file.path(output_dir, "plots/heatmap"))
          enrich <- run_enrichment(res, comp_name, file.path(output_dir, "brain"))
          if (!is.null(enrich)) brain_enrich_results[[comp_name]] <- enrich
        }
      }
    }
  }
}

# ==================== 7.3 跨器官比较（配对） ====================
cat("\n\n===== 3. 跨器官比较（肝癌对脑的影响） =====\n")

for (comp in crosstalk_comparisons) {
  cat("\n---", comp$name, "---\n")
  liver_obj <- seurat_objs[[comp$liver]]
  brain_obj <- seurat_objs[[comp$brain]]
  
  if (is.null(liver_obj) || is.null(brain_obj)) next
  
  for (ct in cell_types) {
    comp_name <- paste0(comp$name, "_", ct)
    cat("\n  >>", comp_name, "\n")
    
    # 提取细胞
    liver_cells <- extract_cells(liver_obj, ct)
    brain_cells <- extract_cells(brain_obj, ct)
    
    if (length(liver_cells) >= MIN_CELLS && length(brain_cells) >= MIN_CELLS) {
      liver_subset <- subset(liver_obj, cells = liver_cells)
      liver_subset$group <- "Liver"
      brain_subset <- subset(brain_obj, cells = brain_cells)
      brain_subset$group <- "Brain"
      
      merged <- merge(liver_subset, brain_subset)
      merged <- JoinLayers(merged)
      
      if (!is.null(merged)) {
        res <- run_limma_diff(merged, comp_name, "Liver", "Brain")
        if (!is.null(res)) {
          crosstalk_deg_results[[comp_name]] <- res
          write.csv(res, file.path(output_dir, "deg_tables", paste0(comp_name, ".csv")), row.names = FALSE)
          plot_volcano(res, comp_name, file.path(output_dir, "plots/volcano"))
          plot_heatmap(res, merged, comp_name, file.path(output_dir, "plots/heatmap"))
          
          # 富集分析
          enrich <- run_enrichment(res, comp_name, file.path(output_dir, "crosstalk"))
          if (!is.null(enrich)) crosstalk_deg_results[[comp_name]] <- enrich
        }
      }
    }
  }
}

# ==================== 8. 肝脑交互分析 ====================
crosstalk_shared <- analyze_crosstalk(liver_enrich_results, brain_enrich_results, 
                                       crosstalk_deg_results, file.path(output_dir, "crosstalk"))

# ==================== 9. 汇总报告 ====================
cat("\n\n========== 生成汇总报告 ==========\n")

report_file <- file.path(output_dir, "analysis_report.txt")
sink(report_file)

cat("========================================\n")
cat("肝脑肿瘤差异分析及交互作用报告\n")
cat("========================================\n")
cat("分析时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("【分析策略】\n")
cat("1. 同器官组间差异分析\n")
cat("   - 肝脏: 肿瘤组 vs 空白组, 肿瘤组 vs 对照组, 空白组 vs 对照组\n")
cat("   - 脑: 肿瘤对应脑 vs 脑空白, 肿瘤对应脑 vs 脑对照, 脑空白 vs 脑对照\n")
cat("2. 跨器官分析（肝癌对脑的影响）\n")
cat("   - gan56 vs nao56, liver3 vs brain3, liver8 vs brain8\n")
cat("3. 肝脑交互分析\n")
cat("   - 找出肝脏和脑中共同富集的通路\n\n")

cat("【分析参数】\n")
cat("  差异分析阈值: |log2FC| >", LOG2FC_CUTOFF, " & adj.P.Val <", ADJ_P_CUTOFF, "\n")
cat("  富集分析阈值: p.adjust < 0.05\n\n")

cat("【成功完成的比较】\n")
cat("\n肝脏内部比较:\n")
for (nm in names(liver_deg_results)) cat("  ✓", nm, "\n")
cat("\n脑内部比较:\n")
for (nm in names(brain_deg_results)) cat("  ✓", nm, "\n")
cat("\n跨器官比较:\n")
for (nm in names(crosstalk_deg_results)) cat("  ✓", nm, "\n")

sink()

cat("\n========================================\n")
cat("分析完成！\n")
cat("输出目录: ", output_dir, "\n")
cat("========================================\n")
