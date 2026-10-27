# ==================== 肝脑肿瘤分析脚本 v2（添加小胶质细胞） ====================
# 文件名: liver_brain_analysis_v2_microglia.R
# 用法: Rscript liver_brain_analysis_v2_microglia.R
#
# 【分析逻辑说明】
# 1. 同器官差异分析（DEG + 富集）
#    - 肝：肿瘤组 vs 空白组 / 肿瘤组 vs 对照组 / 空白组 vs 对照组
#    - 脑：同上，**新增小胶质细胞单独分析**
# 2. 肝脑联系分析（通路层面）
#    - 寻找肝和脑共同富集的KEGG/GO通路 → 肝癌对远端脑的调控机制
#    - 配对样本（同一小鼠）通路相关性分析
# 3. 不做肝vs脑的直接limma差异（器官间基因表达本来差异极大，无意义）
# =========================================================================

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
  library(VennDiagram)
  library(grid)
  library(ComplexHeatmap)
  library(circlize)
})

# ==================== 工作目录与输出 ====================
setwd("/mnt/BioDisk/wangmingyuan/SGS10_analysis")

output_dir <- "liver_brain_v2_microglia_results"
dirs <- c(
  output_dir,
  file.path(output_dir, "01_liver_DEG"),
  file.path(output_dir, "02_brain_DEG"),
  file.path(output_dir, "03_liver_enrichment"),
  file.path(output_dir, "04_brain_enrichment"),
  file.path(output_dir, "05_crosstalk"),
  file.path(output_dir, "06_plots/volcano"),
  file.path(output_dir, "06_plots/heatmap"),
  file.path(output_dir, "06_plots/venn"),
  file.path(output_dir, "06_plots/bubble"),
  file.path(output_dir, "07_summary")
)
for (d in dirs) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

# ==================== 全局参数 ====================
MIN_CELLS <- 5  # 小胶质细胞可能较少，降低到5
MIN_PCT   <- 0.1

# 阈值层级（从严格到宽松，自动选取第一个产生≥10个DEG的层级）
threshold_levels <- list(
  level1 = list(logFC = 1.0, p_cut = 0.05, use_adj = TRUE,  name = "Very_Strict"),
  level2 = list(logFC = 0.8, p_cut = 0.05, use_adj = TRUE,  name = "Strict"),
  level3 = list(logFC = 0.6, p_cut = 0.05, use_adj = TRUE,  name = "Moderate"),
  level4 = list(logFC = 0.5, p_cut = 0.05, use_adj = FALSE, name = "Relaxed"),
  level5 = list(logFC = 0.4, p_cut = 0.05, use_adj = FALSE, name = "Very_Relaxed")
)

# 分析的细胞类型 - **新增小胶质细胞**
analysis_cell_types <- list(
  list(name = "All_Cells",   description = "全部细胞"),
  list(name = "Hepatocytes", description = "肝细胞"),
  list(name = "Endothelial", description = "内皮细胞"),
  list(name = "Kupffer",     description = "枯否细胞"),
  list(name = "Microglia",   description = "小胶质细胞（脑）")  # 新增
)

# 配对样本（同一只小鼠）
paired_samples <- list(
  list(liver = "gan56",  brain = "nao56",  id = "mouse56"),
  list(liver = "liver3", brain = "brain3", id = "mouse3"),
  list(liver = "liver8", brain = "brain8", id = "mouse8")
)

cat("========================================\n")
cat("肝脑肿瘤分析 v2（添加小胶质细胞）\n")
cat("========================================\n\n")

# ==================== 1. 加载所有样本 ====================
cat("========== 1. 加载样本 ==========\n")

samples_to_load <- c("gan56","liver3","liver8","nao56","brain3","brain8",
                     "con_brain","bla_brain","bla_liver","con_liver")

seurat_objs <- list()
for (s in samples_to_load) {
  fp <- file.path(s, paste0(s, "_manual_annotated.rds"))
  if (file.exists(fp)) {
    obj <- readRDS(fp)
    obj <- JoinLayers(obj)
    obj <- RenameCells(obj, add.cell.id = s)
    seurat_objs[[s]] <- obj
    cat(sprintf("  ✓ %-15s 细胞数: %d\n", s, ncol(obj)))
  } else {
    cat(sprintf("  ✗ %-15s 文件不存在\n", s))
  }
}

# ==================== 2. 辅助函数 ====================

# 提取特定细胞类型的细胞ID
extract_cells <- function(obj, cell_type = "All_Cells") {
  if (cell_type != "All_Cells" && "final_celltype" %in% colnames(obj@meta.data)) {
    return(WhichCells(obj, expression = final_celltype == cell_type))
  }
  return(colnames(obj))
}

# 统计细胞数量
count_cells <- function(obj, cell_type = "All_Cells") {
  return(length(extract_cells(obj, cell_type)))
}

# 合并多个对象的特定细胞类型
merge_to_group <- function(obj_list, group_label, cell_type) {
  sub_list <- lapply(obj_list, function(obj) {
    cells <- extract_cells(obj, cell_type)
    if (length(cells) < MIN_CELLS) return(NULL)
    s <- subset(obj, cells = cells)
    s$group <- group_label
    return(s)
  })
  sub_list <- Filter(Negate(is.null), sub_list)
  if (length(sub_list) == 0) return(NULL)
  merged <- Reduce(merge, sub_list)
  merged <- JoinLayers(merged)
  # 移除空库大小的细胞
  expr <- GetAssayData(merged, layer = "data")
  merged <- subset(merged, cells = colnames(expr)[colSums(expr) > 0])
  return(merged)
}

# 合并两组对象
merge_two_groups <- function(obj1, obj2, cell_type, label1, label2) {
  cells1 <- extract_cells(obj1, cell_type)
  cells2 <- extract_cells(obj2, cell_type)
  if (length(cells1) < MIN_CELLS || length(cells2) < MIN_CELLS) return(NULL)
  s1 <- subset(obj1, cells = cells1); s1$group <- label1
  s2 <- subset(obj2, cells = cells2); s2$group <- label2
  merged <- merge(s1, s2)
  merged <- JoinLayers(merged)
  expr <- GetAssayData(merged, layer = "data")
  merged <- subset(merged, cells = colnames(expr)[colSums(expr) > 0])
  return(merged)
}

# ==================== 3. limma差异分析（自适应阈值） ====================
run_limma_adaptive <- function(merged_obj, comp_name, group1, group2, ct_name) {
  n1 <- sum(merged_obj$group == group1)
  n2 <- sum(merged_obj$group == group2)
  if (n1 < MIN_CELLS || n2 < MIN_CELLS) {
    cat(sprintf("    ⚠ 细胞数不足 G1=%d G2=%d，跳过\n", n1, n2))
    return(NULL)
  }

  expr <- GetAssayData(merged_obj, layer = "data")
  keep <- rowSums(expr > 0) >= (MIN_PCT * ncol(expr))
  expr_f <- expr[keep, , drop = FALSE]
  if (nrow(expr_f) < 100) { cat("    ⚠ 基因数<100，跳过\n"); return(NULL) }

  gv <- factor(merged_obj$group, levels = c(group1, group2))
  design <- model.matrix(~ 0 + gv)
  colnames(design) <- levels(gv)

  dge <- tryCatch(calcNormFactors(DGEList(counts = expr_f)), error = function(e) NULL)
  if (is.null(dge)) return(NULL)
  v <- tryCatch(voom(dge, design, plot = FALSE), error = function(e) NULL)
  if (is.null(v)) return(NULL)

  cont <- makeContrasts(contrasts = paste(group1, "-", group2), levels = design)
  fit  <- lmFit(v, design)
  fit2 <- contrasts.fit(fit, cont)
  fit2 <- eBayes(fit2)
  res  <- topTable(fit2, number = Inf, adjust.method = "BH")
  res$gene <- rownames(res)

  # 自适应阈值选择
  used_thr <- threshold_levels[[length(threshold_levels)]]  # 兜底最宽松
  for (lvl in threshold_levels) {
    p_col <- if (lvl$use_adj) "adj.P.Val" else "P.Value"
    sig <- abs(res$logFC) > lvl$logFC & res[[p_col]] < lvl$p_cut
    if (sum(sig) >= 10) { used_thr <- lvl; break }
  }

  p_col <- if (used_thr$use_adj) "adj.P.Val" else "P.Value"
  res$sig <- ifelse(abs(res$logFC) > used_thr$logFC & res[[p_col]] < used_thr$p_cut,
                    ifelse(res$logFC > 0, "Up", "Down"), "NS")
  attr(res, "threshold") <- used_thr

  n_up   <- sum(res$sig == "Up")
  n_down <- sum(res$sig == "Down")
  cat(sprintf("  ✓ %-50s [%s] Up=%-4d Down=%-4d 阈值=%s\n",
              comp_name, ct_name, n_up, n_down, used_thr$name))
  return(res)
}

# ==================== 4. 可视化函数 ====================

# 4.1 火山图（同时输出PDF+PNG）
save_volcano <- function(res, comp_name, ct_name, save_dir) {
  if (is.null(res) || nrow(res) == 0) return()
  thr <- attr(res, "threshold")
  res$log10p <- -log10(res$adj.P.Val + 1e-300)

  top_genes <- bind_rows(
    res %>% filter(sig == "Up")   %>% arrange(adj.P.Val) %>% head(10),
    res %>% filter(sig == "Down") %>% arrange(adj.P.Val) %>% head(10)
  )

  p <- ggplot(res, aes(x = logFC, y = log10p, color = sig)) +
    geom_point(alpha = 0.5, size = 0.8) +
    scale_color_manual(values = c("Up" = "#E41A1C", "Down" = "#377EB8", "NS" = "grey75")) +
    geom_vline(xintercept = c(-thr$logFC, thr$logFC), linetype = "dashed", linewidth = 0.4) +
    geom_hline(yintercept = -log10(thr$p_cut), linetype = "dashed", linewidth = 0.4) +
    geom_text_repel(data = top_genes, aes(label = gene), size = 2.5,
                    max.overlaps = 30, segment.size = 0.2) +
    labs(
      title = paste0(comp_name, "\n", ct_name),
      subtitle = sprintf("|log2FC|>%.1f, %s<%.2f | Up=%d Down=%d",
                         thr$logFC,
                         if (thr$use_adj) "adj.P" else "P",
                         thr$p_cut,
                         sum(res$sig == "Up"),
                         sum(res$sig == "Down")),
      x = "Log2 Fold Change", y = "-Log10 adj.P.Val", color = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(plot.title    = element_text(hjust = 0.5, size = 10, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 8, color = "grey40"),
          legend.position = "right")

  base_path <- file.path(save_dir, paste0(comp_name, "_", ct_name))
  ggsave(paste0(base_path, "_volcano.pdf"), p, width = 8, height = 6, dpi = 300)
  ggsave(paste0(base_path, "_volcano.png"), p, width = 8, height = 6, dpi = 200)
  cat("    ✓ 火山图已保存\n")
}

# 4.2 热图（同时输出PDF+PNG）
save_heatmap <- function(res, merged_obj, comp_name, ct_name, save_dir, n_genes = 40) {
  sig_genes <- res %>% filter(sig != "NS") %>% arrange(desc(abs(logFC))) %>% head(n_genes) %>% pull(gene)
  if (length(sig_genes) < 3) return()

  expr <- GetAssayData(merged_obj, layer = "data")
  expr_f <- as.matrix(expr[intersect(sig_genes, rownames(expr)), ])
  expr_z <- t(scale(t(expr_f))); expr_z[is.na(expr_z)] <- 0

  anno_df <- data.frame(Group = merged_obj$group, row.names = colnames(expr_z))
  
  # 动态创建颜色映射
  unique_groups <- unique(anno_df$Group)
  group_colors_list <- list()
  if (any(grepl("Tumor", unique_groups))) group_colors_list$Tumor <- "#E41A1C"
  if (any(grepl("Blank", unique_groups))) group_colors_list$Blank <- "#984EA3"
  if (any(grepl("Control", unique_groups))) group_colors_list$Control <- "#4DAF4A"
  if (any(grepl("Liver", unique_groups))) group_colors_list$Liver <- "#E41A1C"
  if (any(grepl("Brain", unique_groups))) group_colors_list$Brain <- "#377EB8"
  
  group_colors <- list(Group = group_colors_list)

  # 随机抽样防止热图过宽
  if (ncol(expr_z) > 300) {
    set.seed(42)
    samp <- sample(colnames(expr_z), 300)
    expr_z <- expr_z[, samp]; anno_df <- anno_df[samp, , drop = FALSE]
  }

  col_fun <- colorRamp2(c(-2, 0, 2), c("#4575B4", "white", "#D73027"))

  base_path <- file.path(save_dir, paste0(comp_name, "_", ct_name, "_heatmap"))

  # PDF
  pdf(paste0(base_path, ".pdf"), width = 12, height = max(6, length(sig_genes) * 0.3 + 2))
  ha <- HeatmapAnnotation(df = anno_df)
  ht <- Heatmap(expr_z, name = "Z-score", top_annotation = ha,
                col = col_fun, show_column_names = FALSE, row_names_gp = gpar(fontsize = 8),
                column_title = paste0(comp_name, " (", ct_name, ")"),
                cluster_rows = TRUE, cluster_columns = TRUE)
  draw(ht)
  dev.off()

  # PNG
  png(paste0(base_path, ".png"), width = 2400, height = max(1200, length(sig_genes) * 60 + 400), res = 200)
  draw(ht)
  dev.off()
  cat("    ✓ 热图已保存\n")
}

# ==================== 5. 富集分析函数 ====================
run_enrichment <- function(res, comp_name, ct_name, enrich_dir) {
  if (is.null(res) || nrow(res) == 0) return(NULL)
  
  sig_genes <- res %>% filter(sig != "NS") %>% pull(gene)
  if (length(sig_genes) < 10) {
    cat("    ⚠ DEG<10，跳过富集\n")
    return(NULL)
  }

  gene_df <- tryCatch(
    bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db),
    error = function(e) NULL
  )
  if (is.null(gene_df) || nrow(gene_df) < 5) return(NULL)
  entrez <- unique(gene_df$ENTREZID)

  results <- list()

  # KEGG
  kegg <- tryCatch(enrichKEGG(gene = entrez, organism = "mmu", pvalueCutoff = 0.05),
                   error = function(e) NULL)
  if (!is.null(kegg) && nrow(kegg@result) > 0) {
    write.csv(as.data.frame(kegg),
              file.path(enrich_dir, paste0(comp_name, "_", ct_name, "_KEGG.csv")),
              row.names = FALSE)
    results$kegg <- kegg
    cat(sprintf("    ✓ KEGG: %d 通路\n", nrow(kegg@result[kegg@result$p.adjust < 0.05, ])))
    
    # 富集点图
    if (nrow(kegg@result) >= 5) {
      base_path <- file.path(enrich_dir, paste0(comp_name, "_", ct_name, "_KEGG_dotplot"))
      p_dot <- dotplot(kegg, showCategory = 20, title = paste0("KEGG: ", comp_name, " (", ct_name, ")"))
      ggsave(paste0(base_path, ".pdf"), p_dot, width = 10, height = 8, dpi = 300)
      ggsave(paste0(base_path, ".png"), p_dot, width = 10, height = 8, dpi = 200)
    }
  }

  # GO BP
  go <- tryCatch(
    enrichGO(gene = entrez, OrgDb = org.Mm.eg.db, ont = "BP",
             pvalueCutoff = 0.05, readable = TRUE),
    error = function(e) NULL
  )
  if (!is.null(go) && nrow(go@result) > 0) {
    write.csv(as.data.frame(go),
              file.path(enrich_dir, paste0(comp_name, "_", ct_name, "_GO_BP.csv")),
              row.names = FALSE)
    results$go <- go
    cat(sprintf("    ✓ GO BP: %d 条目\n", nrow(go@result[go@result$p.adjust < 0.05, ])))
  }

  return(results)
}

# ==================== 6. 执行单器官比较的通用函数 ====================
run_organ_comparisons <- function(comparisons, cell_types, seurat_objs,
                                  deg_dir, enrich_dir, plot_dir) {
  all_deg <- list()

  for (comp in comparisons) {
    cat("\n--- ", comp$comp, " ---\n")

    for (ct in cell_types) {
      ct_name <- ct$name
      
      # 跳过肝脏中不存在的小胶质细胞
      if (ct_name == "Microglia" && grepl("liver", comp$comp)) {
        cat(sprintf("  ⚠ %-50s [%s] 肝脏无小胶质细胞，跳过\n", comp$comp, ct_name))
        next
      }

      # 构建合并对象
      merged <- NULL

      if (!is.null(comp$g1_samples)) {
        # 多样本合并 vs 单样本
        g1_objs <- Filter(Negate(is.null), seurat_objs[comp$g1_samples])
        g1_merged <- merge_to_group(g1_objs, comp$g1_label, ct_name)
        g2_cells <- extract_cells(seurat_objs[[comp$g2_obj]], ct_name)

        if (!is.null(g1_merged) && length(g2_cells) >= MIN_CELLS) {
          g2_sub <- subset(seurat_objs[[comp$g2_obj]], cells = g2_cells)
          g2_sub$group <- comp$g2_label
          merged <- merge(g1_merged, g2_sub)
          merged <- JoinLayers(merged)
        }
      } else {
        # 单样本 vs 单样本
        merged <- merge_two_groups(seurat_objs[[comp$g1_obj]],
                                   seurat_objs[[comp$g2_obj]],
                                   ct_name, comp$g1_label, comp$g2_label)
      }

      if (is.null(merged)) {
        cat(sprintf("  ⚠ %-50s [%s] 细胞数不足，跳过\n", comp$comp, ct_name))
        next
      }

      # limma差异分析
      res <- run_limma_adaptive(merged, comp$comp, comp$g1_label, comp$g2_label, ct_name)
      if (is.null(res)) next

      key <- paste0(comp$comp, "_", ct_name)
      all_deg[[key]] <- res

      # 保存DEG表
      write.csv(res, file.path(deg_dir, paste0(key, ".csv")), row.names = FALSE)

      # 可视化
      save_volcano(res, comp$comp, ct_name, file.path(plot_dir, "volcano"))
      save_heatmap(res, merged, comp$comp, ct_name, file.path(plot_dir, "heatmap"))

      # 富集分析
      run_enrichment(res, comp$comp, ct_name, enrich_dir)
    }
  }
  return(all_deg)
}

# ==================== 6.1 肝脏同器官分析 ====================
cat("\n\n========== 3. 肝脏同器官分析 ==========\n")

# 定义肝脏比较组
liver_comparisons <- list(
  list(comp = "liver_tumor_vs_blank",
       g1_samples = c("gan56","liver3","liver8"), g1_label = "Tumor",
       g2_obj = "bla_liver",                     g2_label = "Blank"),
  list(comp = "liver_tumor_vs_control",
       g1_samples = c("gan56","liver3","liver8"), g1_label = "Tumor",
       g2_obj = "con_liver",                     g2_label = "Control"),
  list(comp = "liver_blank_vs_control",
       g1_obj = "bla_liver",  g1_label = "Blank",
       g2_obj = "con_liver",  g2_label = "Control")
)

# 肝脏分析（排除Microglia）
liver_cell_types <- analysis_cell_types[!sapply(analysis_cell_types, function(x) x$name == "Microglia")]

liver_deg <- run_organ_comparisons(
  liver_comparisons, liver_cell_types, seurat_objs,
  file.path(output_dir, "01_liver_DEG"),
  file.path(output_dir, "03_liver_enrichment"),
  file.path(output_dir, "06_plots")
)

# ==================== 6.2 脑同器官分析（包含小胶质细胞） ====================
cat("\n\n========== 4. 脑同器官分析（包含小胶质细胞） ==========\n")

# 定义脑比较组
brain_comparisons <- list(
  list(comp = "brain_tumor_vs_blank",
       g1_samples = c("nao56","brain3","brain8"), g1_label = "TumorBrain",
       g2_obj = "bla_brain",                      g2_label = "Blank"),
  list(comp = "brain_tumor_vs_control",
       g1_samples = c("nao56","brain3","brain8"), g1_label = "TumorBrain",
       g2_obj = "con_brain",                      g2_label = "Control"),
  list(comp = "brain_blank_vs_control",
       g1_obj = "bla_brain",  g1_label = "Blank",
       g2_obj = "con_brain",  g2_label = "Control")
)

# 脑分析（包含所有细胞类型，包括Microglia）
brain_deg <- run_organ_comparisons(
  brain_comparisons, analysis_cell_types, seurat_objs,
  file.path(output_dir, "02_brain_DEG"),
  file.path(output_dir, "04_brain_enrichment"),
  file.path(output_dir, "06_plots")
)

# ==================== 7. 肝脑通路联系分析 ====================
cat("\n\n========== 5. 肝脑通路联系分析 ==========\n")

# 读取肝和脑的KEGG富集结果（肿瘤vs对照，最具生物学意义）
read_kegg_results <- function(enrich_dir, pattern) {
  files <- list.files(enrich_dir, pattern = pattern, full.names = TRUE)
  results <- list()
  for (f in files) {
    key <- gsub("_KEGG.csv$", "", basename(f))
    df <- tryCatch(read.csv(f), error = function(e) NULL)
    if (!is.null(df) && nrow(df) > 0) results[[key]] <- df
  }
  return(results)
}

liver_kegg <- read_kegg_results(file.path(output_dir, "03_liver_enrichment"), "*_KEGG.csv")
brain_kegg  <- read_kegg_results(file.path(output_dir, "04_brain_enrichment"), "*_KEGG.csv")

cat(sprintf("  肝脏KEGG富集结果集合: %d 个\n", length(liver_kegg)))
cat(sprintf("  脑部KEGG富集结果集合: %d 个\n", length(brain_kegg)))

if (length(liver_kegg) > 0 && length(brain_kegg) > 0) {

  # 7.1 汇总所有肝/脑通路（肿瘤vs对照）
  summarize_pathways <- function(kegg_list, tag) {
    bind_rows(lapply(names(kegg_list), function(nm) {
      df <- kegg_list[[nm]]
      df$source <- nm
      df$organ  <- tag
      return(df)
    }))
  }

  liver_path_all <- summarize_pathways(liver_kegg, "Liver")
  brain_path_all <- summarize_pathways(brain_kegg,  "Brain")

  # 7.2 肿瘤vs对照的通路（最核心比较）
  liver_tvsc <- liver_path_all %>% filter(grepl("tumor_vs_control", source))
  brain_tvsc  <- brain_path_all  %>% filter(grepl("tumor_vs_control", source))

  liver_path_ids <- unique(liver_tvsc$ID)
  brain_path_ids  <- unique(brain_tvsc$ID)
  shared_ids      <- intersect(liver_path_ids, brain_path_ids)

  cat(sprintf("\n  肝脏显著通路数: %d\n", length(liver_path_ids)))
  cat(sprintf("  脑部显著通路数: %d\n", length(brain_path_ids)))
  cat(sprintf("  共同显著通路数: %d\n", length(shared_ids)))

  # 7.3 Venn图（PDF+PNG）
  if (length(liver_path_ids) > 0 || length(brain_path_ids) > 0) {
    venn_base <- file.path(output_dir, "06_plots/venn", "liver_brain_pathway_venn")
    venn_pdf <- paste0(venn_base, ".pdf")
    venn_png <- paste0(venn_base, ".png")

    venn_list <- list(Liver = liver_path_ids, Brain = brain_path_ids)

    pdf(venn_pdf, width = 7, height = 6)
    vd <- venn.diagram(venn_list, filename = NULL,
                       fill = c("#E41A1C", "#377EB8"), alpha = 0.5,
                       cex = 1.5, cat.cex = 1.5, cat.pos = c(-20, 20),
                       main = "肿瘤组 vs 对照组：肝脑共同富集KEGG通路",
                       main.cex = 1.2)
    grid.draw(vd)
    dev.off()

    png(venn_png, width = 1400, height = 1200, res = 200)
    grid.draw(vd)
    dev.off()
    cat("  ✓ Venn图已保存\n")
  }

  # 7.4 共同通路汇总表
  if (length(shared_ids) > 0) {
    liver_shared <- liver_tvsc %>%
      filter(ID %in% shared_ids) %>%
      group_by(ID, Description) %>%
      summarise(Liver_minP = min(p.adjust), Liver_meanP = mean(p.adjust),
                Liver_n = n(), .groups = "drop")

    brain_shared <- brain_tvsc %>%
      filter(ID %in% shared_ids) %>%
      group_by(ID, Description) %>%
      summarise(Brain_minP = min(p.adjust), Brain_meanP = mean(p.adjust),
                Brain_n = n(), .groups = "drop")

    shared_df <- inner_join(liver_shared, brain_shared, by = c("ID", "Description")) %>%
      mutate(combined_score = -log10(Liver_minP) + -log10(Brain_minP)) %>%
      arrange(desc(combined_score))

    write.csv(shared_df,
              file.path(output_dir, "05_crosstalk", "shared_KEGG_pathways.csv"),
              row.names = FALSE)

    cat("\n  【肝脑共同富集通路 Top20】\n")
    print(head(shared_df %>% select(ID, Description, Liver_minP, Brain_minP, combined_score), 20))

    # 7.5 联合气泡图
    if (nrow(shared_df) >= 3) {
      plot_df <- shared_df %>%
        head(25) %>%
        mutate(
          Liver_logP = -log10(Liver_minP),
          Brain_logP = -log10(Brain_minP),
          Desc_short = ifelse(nchar(Description) > 45,
                              paste0(substr(Description, 1, 42), "..."),
                              Description)
        )

      # 散点图
      p_scatter <- ggplot(plot_df, aes(x = Liver_logP, y = Brain_logP, label = Desc_short)) +
        geom_point(aes(size = combined_score, color = combined_score), alpha = 0.8) +
        scale_color_gradientn(colors = c("#4575B4","#91BFDB","#FFFFBF","#FC8D59","#D73027")) +
        geom_text_repel(size = 2.8, max.overlaps = 20, segment.size = 0.2) +
        labs(title = "肝脑共同KEGG通路",
             subtitle = "X轴：肝脏显著性；Y轴：脑部显著性",
             x = "-Log10 P (Liver)", y = "-Log10 P (Brain)",
             color = "Combined\nScore", size = "Combined\nScore") +
        theme_bw(base_size = 12) +
        theme(plot.title    = element_text(hjust = 0.5, face = "bold"),
              plot.subtitle = element_text(hjust = 0.5, color = "grey40"))

      ggsave(file.path(output_dir, "06_plots/bubble", "shared_pathway_scatter.pdf"),
             p_scatter, width = 10, height = 8, dpi = 300)
      ggsave(file.path(output_dir, "06_plots/bubble", "shared_pathway_scatter.png"),
             p_scatter, width = 10, height = 8, dpi = 200)

      # 条形图
      plot_long <- plot_df %>%
        select(Desc_short, Liver_logP, Brain_logP) %>%
        tidyr::pivot_longer(cols = c(Liver_logP, Brain_logP),
                            names_to = "Organ", values_to = "logP") %>%
        mutate(Organ = recode(Organ, "Liver_logP" = "Liver", "Brain_logP" = "Brain"))

      p_bar <- ggplot(plot_long, aes(x = reorder(Desc_short, logP), y = logP, fill = Organ)) +
        geom_bar(stat = "identity", position = "dodge") +
        scale_fill_manual(values = c("Liver" = "#E41A1C", "Brain" = "#377EB8")) +
        coord_flip() +
        labs(title = "肝脑共同富集通路显著性比较",
             x = "", y = "-Log10 P-value", fill = "器官") +
        theme_bw(base_size = 11) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              axis.text.y = element_text(size = 8))

      h <- max(6, nrow(plot_df) * 0.4 + 2)
      ggsave(file.path(output_dir, "06_plots/bubble", "shared_pathway_barplot.pdf"),
             p_bar, width = 12, height = h, dpi = 300)
      ggsave(file.path(output_dir, "06_plots/bubble", "shared_pathway_barplot.png"),
             p_bar, width = 12, height = h, dpi = 200)
      cat("  ✓ 共同通路可视化已保存\n")
    }
  }

  # 7.6 配对样本通路分析（重点分析小胶质细胞）
  cat("\n\n--- 5.2 配对样本通路分析（包含小胶质细胞）---\n")

  for (pair in paired_samples) {
    cat(sprintf("\n  配对: %s (肝: %s, 脑: %s)\n", pair$id, pair$liver, pair$brain))

    for (ct in analysis_cell_types) {
      ct_name <- ct$name
      
      # 肝脏配对中跳过小胶质细胞
      if (ct_name == "Microglia") {
        cat(sprintf("    [%s] 肝脏无小胶质细胞，跳过\n", ct_name))
        next
      }

      comp_name <- paste0("paired_", pair$id, "_", ct_name)

      liver_obj <- seurat_objs[[pair$liver]]
      brain_obj  <- seurat_objs[[pair$brain]]

      liver_cnt <- length(extract_cells(liver_obj, ct_name))
      brain_cnt  <- length(extract_cells(brain_obj, ct_name))

      cat(sprintf("    [%s] 肝=%d 脑=%d 细胞 ", ct_name, liver_cnt, brain_cnt))

      if (liver_cnt < MIN_CELLS || brain_cnt < MIN_CELLS) {
        cat("→ 跳过\n"); next
      }

      cat("→ 分析中\n")

      # 肝脏：肿瘤 vs 对照
      liver_vs_ctrl <- merge_two_groups(
        liver_obj, seurat_objs[["con_liver"]], ct_name, "TumorLiver", "Control"
      )

      # 脑：肿瘤对应脑 vs 脑对照
      brain_vs_ctrl <- merge_two_groups(
        brain_obj, seurat_objs[["con_brain"]], ct_name, "TumorBrain", "Control"
      )

      if (!is.null(liver_vs_ctrl)) {
        res_l <- run_limma_adaptive(liver_vs_ctrl, paste0("paired_", pair$id, "_liver"), 
                                    "TumorLiver", "Control", ct_name)
        if (!is.null(res_l)) {
          run_enrichment(res_l, paste0("paired_", pair$id, "_liver"), ct_name,
                         file.path(output_dir, "05_crosstalk"))
        }
      }

      if (!is.null(brain_vs_ctrl)) {
        res_b <- run_limma_adaptive(brain_vs_ctrl, paste0("paired_", pair$id, "_brain"), 
                                    "TumorBrain", "Control", ct_name)
        if (!is.null(res_b)) {
          run_enrichment(res_b, paste0("paired_", pair$id, "_brain"), ct_name,
                         file.path(output_dir, "05_crosstalk"))
        }
      }
    }
  }
} else {
  cat("  ⚠ 肝或脑的KEGG富集结果不足，跳过肝脑交互分析\n")
}

# ==================== 8. 小胶质细胞特异性分析报告 ====================
cat("\n\n========== 6. 小胶质细胞特异性分析报告 ==========\n")

# 提取小胶质细胞相关结果
microglia_results <- list()

for (nm in names(brain_deg)) {
  if (grepl("Microglia", nm)) {
    microglia_results[[nm]] <- brain_deg[[nm]]
  }
}

if (length(microglia_results) > 0) {
  cat(sprintf("  成功分析的小胶质细胞比较: %d 个\n", length(microglia_results)))
  
  # 小胶质细胞DEG汇总
  microglia_summary <- do.call(rbind, lapply(names(microglia_results), function(nm) {
    res <- microglia_results[[nm]]
    data.frame(
      Comparison = nm,
      Up = sum(res$sig == "Up", na.rm = TRUE),
      Down = sum(res$sig == "Down", na.rm = TRUE),
      Total = sum(res$sig != "NS", na.rm = TRUE),
      Threshold = attr(res, "threshold")$name,
      stringsAsFactors = FALSE
    )
  }))
  
  write.csv(microglia_summary, 
            file.path(output_dir, "07_summary", "microglia_DEG_summary.csv"), 
            row.names = FALSE)
  
  cat("\n  小胶质细胞DEG统计:\n")
  print(microglia_summary)
} else {
  cat("  ⚠ 未找到小胶质细胞差异分析结果\n")
}

# ==================== 9. 汇总报告 ====================
cat("\n\n========== 7. 生成汇总报告 ==========\n")

all_deg <- c(liver_deg, brain_deg)

# DEG数量汇总表
if (length(all_deg) > 0) {
  deg_summary <- do.call(rbind, lapply(names(all_deg), function(nm) {
    res <- all_deg[[nm]]
    data.frame(
      Comparison = nm,
      Up         = sum(res$sig == "Up", na.rm = TRUE),
      Down       = sum(res$sig == "Down", na.rm = TRUE),
      Total      = sum(res$sig != "NS", na.rm = TRUE),
      Threshold  = attr(res, "threshold")$name,
      stringsAsFactors = FALSE
    )
  }))
  
  write.csv(deg_summary, file.path(output_dir, "07_summary", "DEG_count_summary.csv"), row.names = FALSE)
} else {
  deg_summary <- data.frame()
  cat("  ⚠ 没有成功完成的差异分析\n")
}

# 文字报告
report_path <- file.path(output_dir, "07_summary", "analysis_report.txt")
sink(report_path)
cat("========================================\n")
cat("肝脑肿瘤分析报告 v2（添加小胶质细胞）\n")
cat("========================================\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("【分析框架】\n")
cat("1. 同器官差异分析：肝（3组对比）+ 脑（3组对比，包含小胶质细胞）\n")
cat("2. 肝脑通路联系：通过共同富集KEGG通路揭示肝癌对脑的远端影响\n")
cat("3. 配对样本分析：同一小鼠肝脑通路对比\n")
cat("4. 小胶质细胞特异性分析：重点关注脑内小胶质细胞的变化\n\n")

cat("【分析细胞类型】\n")
for (ct in analysis_cell_types) {
  cat(sprintf("  - %s: %s\n", ct$name, ct$description))
}
cat("\n")

cat("【阈值策略】\n")
for (i in 1:length(threshold_levels)) {
  lvl <- threshold_levels[[i]]
  cat(sprintf("  Level %d (%s): |log2FC| > %.1f, %s < 0.05\n", 
              i, lvl$name, lvl$logFC, ifelse(lvl$use_adj, "adj.P.Val", "P.Value")))
}
cat("  自动选择第一个产生≥10个DEG的阈值\n")
cat("  小胶质细胞使用最小细胞数: 5\n\n")

cat("【DEG数量统计】\n")
if (nrow(deg_summary) > 0) {
  print(deg_summary)
} else {
  cat("  无显著DEG\n")
}

cat("\n【小胶质细胞特异性结果】\n")
if (exists("microglia_summary") && nrow(microglia_summary) > 0) {
  print(microglia_summary)
} else {
  cat("  小胶质细胞分析无显著结果\n")
}

sink()

cat("✓ 报告已保存:", report_path, "\n")

# 打印最终统计
cat("\n========================================\n")
cat("分析完成！\n")
cat(sprintf("  肝脏比较: %d 个\n", length(liver_deg)))
cat(sprintf("  脑部比较: %d 个\n", length(brain_deg)))
cat(sprintf("  小胶质细胞比较: %d 个\n", length(microglia_results)))
cat(sprintf("  输出目录: %s\n", output_dir))
cat("========================================\n")
