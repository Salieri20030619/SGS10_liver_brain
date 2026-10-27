# ==================== 肝脑肿瘤完整分析脚本（修复版） ====================
# 文件名: liver_brain_complete_analysis_fixed.R
# 用法: Rscript liver_brain_complete_analysis_fixed.R

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
cat("肝脑肿瘤完整分析（修复版）\n")
cat("========================================\n")
cat("R版本:", version$version.string, "\n")
cat("clusterProfiler版本:", as.character(packageVersion("clusterProfiler")), "\n")
cat("阈值: |log2FC| >", LOG2FC_CUTOFF, " & adj.P.Val <", ADJ_P_CUTOFF, "\n\n")

# ==================== 定义比较组 ====================

# 肝脏比较
liver_comparisons <- list(
  list(name = "tumor_vs_blank", tumor = c("gan56", "liver3", "liver8"), control = "bla_liver", type = "liver"),
  list(name = "tumor_vs_control", tumor = c("gan56", "liver3", "liver8"), control = "con_liver", type = "liver"),
  list(name = "blank_vs_control", group1 = "bla_liver", group2 = "con_liver", type = "liver")
)

# 脑比较
brain_comparisons <- list(
  list(name = "brain_tumor_vs_blank", tumor_brain = c("nao56", "brain3", "brain8"), control = "bla_brain", type = "brain"),
  list(name = "brain_tumor_vs_control", tumor_brain = c("nao56", "brain3", "brain8"), control = "con_brain", type = "brain"),
  list(name = "brain_blank_vs_control", group1 = "bla_brain", group2 = "con_brain", type = "brain")
)

# 跨器官比较
crosstalk_comparisons <- list(
  list(name = "gan56_vs_nao56", liver = "gan56", brain = "nao56"),
  list(name = "liver3_vs_brain3", liver = "liver3", brain = "brain3"),
  list(name = "liver8_vs_brain8", liver = "liver8", brain = "brain8")
)

# 细胞类型
cell_types <- c("Hepatocytes", "Endothelial", "Kupffer", "T_cells", "B_cells", "Stellate")

# ==================== 1. 加载样本 ====================
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
  
  # 检查并过滤
  expr <- GetAssayData(merged, layer = "data")
  lib_sizes <- colSums(expr)
  keep_cells <- lib_sizes > 0
  if (sum(!keep_cells) > 0) {
    merged <- subset(merged, cells = colnames(expr)[keep_cells])
  }
  
  return(merged)
}

# ==================== 3. 差异分析函数（修复版） ====================
run_limma_diff <- function(merged_obj, comp_name, group1, group2) {
  Idents(merged_obj) <- merged_obj$group
  n1 <- sum(merged_obj$group == group1)
  n2 <- sum(merged_obj$group == group2)
  
  if (n1 < MIN_CELLS || n2 < MIN_CELLS) {
    cat(sprintf("  ⚠️ 细胞数不足 (G1=%d, G2=%d)，跳过\n", n1, n2))
    return(NULL)
  }

  expr <- GetAssayData(merged_obj, layer = "data")
  
  # 过滤低表达基因
  keep <- rowSums(expr > 0) >= (MIN_PCT * ncol(expr))
  expr_f <- expr[keep, , drop = FALSE]
  
  if (nrow(expr_f) < 100) {
    cat("    ⚠️ 过滤后基因数 < 100，跳过\n")
    return(NULL)
  }
  
  # 检查并过滤库大小为0的细胞
  lib_sizes <- colSums(expr_f)
  if (any(lib_sizes == 0)) {
    keep_cells <- lib_sizes > 0
    expr_f <- expr_f[, keep_cells, drop = FALSE]
    merged_obj <- subset(merged_obj, cells = colnames(expr_f))
    if (ncol(expr_f) < MIN_CELLS * 2) return(NULL)
  }

  # 构建设计矩阵
  gv <- factor(merged_obj$group, levels = c(group1, group2))
  design <- model.matrix(~ 0 + gv)
  colnames(design) <- levels(gv)
  
  # 检查设计矩阵是否满秩
  if (qr(design)$rank < ncol(design)) {
    cat("    ⚠️ 设计矩阵不满秩，跳过\n")
    return(NULL)
  }
  
  cont <- makeContrasts(contrasts = paste(group1, "-", group2), levels = design)

  dge <- tryCatch({
    DGEList(counts = expr_f)
  }, error = function(e) {
    cat("    DGEList错误:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(dge)) return(NULL)
  
  # 计算归一化因子
  dge <- tryCatch({
    calcNormFactors(dge)
  }, error = function(e) {
    cat("    calcNormFactors错误:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(dge)) return(NULL)
  
  # voom转换
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
  res$sig <- ifelse(abs(res$logFC) > LOG2FC_CUTOFF & res$adj.P.Val < ADJ_P_CUTOFF,
               ifelse(res$logFC > 0, "Up", "Down"), "NS")

  n_up <- sum(res$sig == "Up"); n_down <- sum(res$sig == "Down")
  cat(sprintf("  ✓ %s | Up=%d Down=%d Total=%d\n", comp_name, n_up, n_down, n_up + n_down))
  
  return(res)
}

# ==================== 4. 富集分析函数 ====================
run_enrichment <- function(deg_df, comp_name, output_dir) {
  if (is.null(deg_df) || nrow(deg_df) == 0) return(NULL)
  
  sig_genes <- deg_df %>% filter(sig != "NS") %>% pull(gene)
  if (length(sig_genes) < 10) {
    cat("    ⚠️ 显著基因数 < 10，跳过富集分析\n")
    return(NULL)
  }
  
  gene_df <- tryCatch({
    bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  }, error = function(e) {
    cat("    bitr转换失败:", e$message, "\n")
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
              file.path(output_dir, paste0(comp_name, "_KEGG.csv")), 
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
              file.path(output_dir, paste0(comp_name, "_GO_BP.csv")), 
              row.names = FALSE)
    cat("    ✓ GO富集完成，条目数:", nrow(go_res@result), "\n")
  }
  
  return(list(kegg = kegg_res, go = go_res))
}

# ==================== 5. 可视化函数 ====================
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
  cat("    ✓ 火山图已保存\n")
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
  cat("    ✓ 热图已保存\n")
}

# ==================== 6. 肝脑交互分析 ====================
analyze_crosstalk <- function(liver_enrich, brain_enrich, output_dir) {
  cat("\n========== 肝脑交互分析 ==========\n")
  
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
  
  if (length(liver_pathways) == 0 || length(brain_pathways) == 0) {
    cat("  ⚠️ 富集结果不足，跳过交互分析\n")
    return(NULL)
  }
  
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
    write.csv(shared, file.path(output_dir, "shared_pathways.csv"), row.names = FALSE)
  }
  
  return(shared)
}

# ==================== 7. 主分析 ====================
cat("\n========== 开始分析 ==========\n\n")

liver_deg <- list()
brain_deg <- list()
crosstalk_deg <- list()
liver_enrich <- list()
brain_enrich <- list()

# 肝脏分析
cat("\n===== 1. 肝脏内部比较 =====\n")

# 肿瘤 vs 空白
cat("\n--- 肿瘤 vs 空白 ---\n")
for (ct in cell_types) {
  comp_name <- paste0("liver_tumor_vs_blank_", ct)
  cat("\n  >>", comp_name, "\n")
  
  tumor_objs <- list()
  for (t in c("gan56", "liver3", "liver8")) {
    if (t %in% names(seurat_objs)) tumor_objs[[t]] <- seurat_objs[[t]]
  }
  blank_obj <- seurat_objs[["bla_liver"]]
  
  if (length(tumor_objs) > 0 && !is.null(blank_obj)) {
    # 检查每个肿瘤样本中该细胞类型的细胞数
    has_cells <- TRUE
    for (t_obj in tumor_objs) {
      cells <- extract_cells(t_obj, ct)
      if (length(cells) < MIN_CELLS) {
        cat("    肿瘤样本细胞数不足，跳过\n")
        has_cells <- FALSE
        break
      }
    }
    
    blank_cells <- extract_cells(blank_obj, ct)
    if (length(blank_cells) < MIN_CELLS) {
      cat("    空白样本细胞数不足，跳过\n")
      has_cells <- FALSE
    }
    
    if (has_cells) {
      merged_tumor <- merge_objects(tumor_objs, rep("Tumor", length(tumor_objs)), ct)
      
      if (!is.null(merged_tumor)) {
        blank_subset <- subset(blank_obj, cells = blank_cells)
        blank_subset$group <- "Blank"
        merged <- merge(merged_tumor, blank_subset)
        merged <- JoinLayers(merged)
        
        if (!is.null(merged)) {
          res <- run_limma_diff(merged, comp_name, "Tumor", "Blank")
          if (!is.null(res)) {
            liver_deg[[comp_name]] <- res
            write.csv(res, file.path(output_dir, "deg_tables", paste0(comp_name, ".csv")), row.names = FALSE)
            plot_volcano(res, comp_name, file.path(output_dir, "plots/volcano"))
            plot_heatmap(res, merged, comp_name, file.path(output_dir, "plots/heatmap"))
            enrich <- run_enrichment(res, comp_name, file.path(output_dir, "liver"))
            if (!is.null(enrich)) liver_enrich[[comp_name]] <- enrich
          }
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
    has_cells <- TRUE
    for (t_obj in tumor_objs) {
      cells <- extract_cells(t_obj, ct)
      if (length(cells) < MIN_CELLS) {
        has_cells <- FALSE
        break
      }
    }
    
    control_cells <- extract_cells(control_obj, ct)
    if (length(control_cells) < MIN_CELLS) has_cells <- FALSE
    
    if (has_cells) {
      merged_tumor <- merge_objects(tumor_objs, rep("Tumor", length(tumor_objs)), ct)
      
      if (!is.null(merged_tumor)) {
        control_subset <- subset(control_obj, cells = control_cells)
        control_subset$group <- "Control"
        merged <- merge(merged_tumor, control_subset)
        merged <- JoinLayers(merged)
        
        if (!is.null(merged)) {
          res <- run_limma_diff(merged, comp_name, "Tumor", "Control")
          if (!is.null(res)) {
            liver_deg[[comp_name]] <- res
            write.csv(res, file.path(output_dir, "deg_tables", paste0(comp_name, ".csv")), row.names = FALSE)
            plot_volcano(res, comp_name, file.path(output_dir, "plots/volcano"))
            plot_heatmap(res, merged, comp_name, file.path(output_dir, "plots/heatmap"))
            enrich <- run_enrichment(res, comp_name, file.path(output_dir, "liver"))
            if (!is.null(enrich)) liver_enrich[[comp_name]] <- enrich
          }
        }
      }
    }
  }
}

cat("\n\n========================================\n")
cat("分析完成！\n")
cat("输出目录: ", output_dir, "\n")
cat("========================================\n")
