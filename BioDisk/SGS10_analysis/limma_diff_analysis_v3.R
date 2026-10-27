# ==================== 使用limma进行单细胞差异分析 ====================
# 文件名: limma_diff_analysis_v3.R
# 版本说明: 符合数据分析部门项目分析标准 v5.0
# 分析内容:
#   1. 质控可视化（PCA图 + 箱线图）
#   2. 差异分析（limma/voom，单细胞数据）
#      - 主阈值: |log2FC|>1 & adj.P.Val<0.05
#      - 备用阈值: |log2FC|>0.5 & P.Value<0.05（主阈值DEG<100时自动启用）
#   3. 肿瘤特异性DEG筛选（排除对照组基础差异）
#   4. 火山图 + 热图可视化
#   5. KEGG/GO富集分析（ORA，优先padj，基因数≥10）
#   6. GSEA（|NES|>1 & p<0.05，背景基因集优先KEGG）
#   7. 肝脑通路交集分析与可视化
# 用法: Rscript limma_diff_analysis_v3.R

# ==================== 加载依赖包 ====================
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
  library(org.Mm.eg.db)      # 鼠源基因注释（鼠源细胞用Mm）
  library(enrichplot)
  library(ggpubr)
  library(patchwork)
})

# ==================== 工作目录与输出目录 ====================
setwd("/mnt/BioDisk/wangmingyuan/SGS10_analysis")

dirs <- c(
  "limma_diff_results",
  "limma_diff_results/plots",
  "limma_diff_results/plots/qc",
  "limma_diff_results/plots/volcano",
  "limma_diff_results/plots/heatmap",
  "limma_diff_results/plots/enrichment",
  "limma_diff_results/plots/gsea",
  "limma_diff_results/plots/liver_brain_crosstalk",
  "limma_diff_results/deg_tables",
  "limma_diff_results/tumor_specific_deg"
)
for (d in dirs) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

output_dir   <- "limma_diff_results"
plot_dir     <- file.path(output_dir, "plots")
deg_dir      <- file.path(output_dir, "deg_tables")
specific_dir <- file.path(output_dir, "tumor_specific_deg")

# ==================== 全局参数（符合标准v5.0） ====================
# 差异分析阈值
LOG2FC_CUTOFF    <- 1      # 主log2FC阈值
ADJ_P_CUTOFF     <- 0.05  # 主校正p值阈值（优先使用padj）
LOG2FC_ALT       <- 0.5   # 备用log2FC阈值（最低标准）
P_ALT            <- 0.05  # 备用p值阈值（最低标准）
MIN_DEG_COUNT    <- 100   # DEG数量下限（低于此需特别处理）
MIN_CELLS        <- 10    # 最小细胞数
MIN_PCT          <- 0.1   # 最小表达比例（基因过滤）

# 富集分析阈值（符合标准v5.0）
ORA_PADJ_CUTOFF  <- 0.05  # ORA校正p值阈值（优先padj）
ORA_MIN_GENES    <- 10    # 富集最少基因数（低于10转GSEA）
GSEA_NES_CUTOFF  <- 1     # GSEA NES绝对值阈值
GSEA_P_CUTOFF    <- 0.05  # GSEA p值阈值（adjust.P<0.05）

# 相关性分析（标准v5.0：优先spearman，最低阈值0.3）
COR_METHOD       <- "spearman"
COR_MIN          <- 0.3

cat("========================================\n")
cat("limma 差异分析 v3（符合标准v5.0）\n")
cat("========================================\n")
cat("主阈值:   |log2FC| >", LOG2FC_CUTOFF, " & adj.P.Val <", ADJ_P_CUTOFF, "\n")
cat("备用阈值: |log2FC| >", LOG2FC_ALT,   " & P.Value <",   P_ALT, "\n")
cat("DEG数量下限:", MIN_DEG_COUNT, "个\n")
cat("富集最少基因数:", ORA_MIN_GENES, "个\n\n")

# ==================== 比较组定义 ====================
# 肿瘤组 vs 空白/对照组（核心比较）
tumor_vs_blank <- list(
  # 肝脏轴：肿瘤 vs 空白肝（排除对照影响的主要比较）
  list(name = "gan56_vs_bla_liver",  group1 = "gan56",  group2 = "bla_liver",
       organ = "liver", type = "tumor_vs_blank"),
  list(name = "liver3_vs_bla_liver", group1 = "liver3", group2 = "bla_liver",
       organ = "liver", type = "tumor_vs_blank"),
  list(name = "liver8_vs_bla_liver", group1 = "liver8", group2 = "bla_liver",
       organ = "liver", type = "tumor_vs_blank"),

  # 脑轴：肿瘤 vs 空白脑
  list(name = "gan56_vs_nao56",      group1 = "gan56",  group2 = "nao56",
       organ = "brain", type = "tumor_vs_blank"),
  list(name = "liver3_vs_brain3",    group1 = "liver3", group2 = "brain3",
       organ = "brain", type = "tumor_vs_blank"),
  list(name = "liver8_vs_brain8",    group1 = "liver8", group2 = "brain8",
       organ = "brain", type = "tumor_vs_blank")
)

# 对照组 vs 空白组（用于排除基础差异）
ctrl_vs_blank <- list(
  list(name = "con_liver_vs_bla_liver", group1 = "con_liver", group2 = "bla_liver",
       organ = "liver", type = "ctrl_vs_blank"),
  list(name = "con_brain_vs_bla_brain", group1 = "con_brain", group2 = "bla_brain",
       organ = "brain", type = "ctrl_vs_blank")
)

# 其余组间比较（可选扩展）
other_comps <- list(
  list(name = "bla_liver_vs_bla_brain", group1 = "bla_liver", group2 = "bla_brain",
       organ = "both", type = "blank_liver_vs_brain"),
  list(name = "con_liver_vs_con_brain", group1 = "con_liver", group2 = "con_brain",
       organ = "both", type = "ctrl_liver_vs_brain"),
  list(name = "gan56_vs_con_liver",     group1 = "gan56",     group2 = "con_liver",
       organ = "liver", type = "tumor_vs_ctrl")
)

all_comparisons <- c(tumor_vs_blank, ctrl_vs_blank, other_comps)

# 分析的细胞类型
cell_types <- c("Hepatocytes", "Kupffer", "Endothelial", "T_cells", "B_cells", "Stellate")

# ==================== 1. 加载样本 ====================
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

# ==================== 2. 【修复】人鼠细胞分类（基因命名规则法） ====================
# 原版用7个肝细胞标记导致大量Unknown（liver8达52%）
# 修复：全大写=人源，首字母大写+小写=鼠源，使用比例阈值（>=0.6 Human，<=0.4 Mouse）
classify_cells_by_name <- function(obj, sample_name) {
  cat("\n分类", sample_name, "（基因命名规则法）...\n")
  all_genes   <- rownames(obj)
  human_genes <- all_genes[grepl("^[A-Z][A-Z0-9-]+$", all_genes)]
  mouse_genes <- all_genes[grepl("^[A-Z][a-z]", all_genes)]
  cat("  人源基因数:", length(human_genes), " 鼠源基因数:", length(mouse_genes), "\n")

  expr         <- GetAssayData(obj, layer = "data")
  human_score  <- colSums(expr[intersect(human_genes, rownames(expr)), , drop = FALSE] > 0)
  mouse_score  <- colSums(expr[intersect(mouse_genes,  rownames(expr)), , drop = FALSE] > 0)
  total_score  <- human_score + mouse_score
  human_ratio  <- ifelse(total_score > 0, human_score / total_score, 0)

  cell_origin  <- ifelse(human_ratio >= 0.6, "Human",
                   ifelse(human_ratio <= 0.4, "Mouse", "Unknown"))

  obj$cell_origin <- cell_origin
  obj$human_score <- human_score
  obj$mouse_score <- mouse_score
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

# ==================== 3. 质控可视化（标准v5.0：需提前做PCA和箱线图） ====================
cat("\n========== 3. 质控可视化（PCA + 箱线图）==========\n")

plot_qc_for_comparison <- function(obj1, obj2, comp_name, cell_type = NULL) {
  # 提取细胞
  get_cells <- function(obj, ct) {
    if (!is.null(ct) && "final_celltype" %in% colnames(obj@meta.data))
      WhichCells(obj, expression = final_celltype == ct)
    else colnames(obj)
  }
  c1 <- get_cells(obj1, cell_type)
  c2 <- get_cells(obj2, cell_type)
  if (length(c1) < MIN_CELLS || length(c2) < MIN_CELLS) return(invisible(NULL))

  s1 <- subset(obj1, cells = c1); s1$plot_group <- unique(obj1$orig.ident)[1]
  s2 <- subset(obj2, cells = c2); s2$plot_group <- unique(obj2$orig.ident)[1]
  merged <- merge(s1, s2); merged <- JoinLayers(merged)

  expr <- GetAssayData(merged, layer = "data")
  keep <- rowSums(expr > 0) >= (MIN_PCT * ncol(expr))
  expr_f <- as.matrix(expr[keep, ])

  # --- 箱线图：各样本基因表达中位数分布 ---
  # 用于检查组内是否有批次效应或离群样本
  expr_df <- data.frame(
    sample     = merged$plot_group,
    median_exp = apply(expr_f, 2, function(x) median(x[x > 0]))
  )
  p_box <- ggplot(expr_df, aes(x = sample, y = median_exp, fill = sample)) +
    geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
    geom_jitter(width = 0.2, size = 0.3, alpha = 0.3) +
    labs(title  = paste0(comp_name, if(!is.null(cell_type)) paste0("_", cell_type)),
         x      = "样本组", y     = "基因表达中位数（仅表达基因）",
         caption= "箱线图用于检查组内样本是否存在批次或离群") +
    theme_bw(base_size = 11) +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

  # --- PCA图：检查样本聚类是否符合预期 ---
  # 取top2000高变基因做PCA
  n_top <- min(2000, nrow(expr_f))
  gene_var  <- apply(expr_f, 1, var)
  top_genes <- names(sort(gene_var, decreasing = TRUE))[1:n_top]
  pca_res   <- prcomp(t(expr_f[top_genes, ]), scale. = TRUE, center = TRUE)
  pca_df    <- data.frame(
    PC1    = pca_res$x[, 1],
    PC2    = pca_res$x[, 2],
    sample = merged$plot_group
  )
  var_exp <- round(summary(pca_res)$importance[2, 1:2] * 100, 1)

  p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = sample)) +
    geom_point(size = 0.8, alpha = 0.7) +
    stat_ellipse(level = 0.95, linetype = "dashed") +
    labs(title   = paste0("PCA: ", comp_name, if(!is.null(cell_type)) paste0("_", cell_type)),
         x       = paste0("PC1 (", var_exp[1], "%)"),
         y       = paste0("PC2 (", var_exp[2], "%)"),
         caption = "PCA用于检查样本分离情况，椭圆为95%置信区间") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5))

  # 保存
  fname <- paste0(comp_name, if(!is.null(cell_type)) paste0("_", cell_type))
  ggsave(file.path(plot_dir, "qc", paste0(fname, "_boxplot.pdf")),
         p_box, width = 6, height = 5, dpi = 300)
  ggsave(file.path(plot_dir, "qc", paste0(fname, "_PCA.pdf")),
         p_pca, width = 7, height = 6, dpi = 300)

  cat("    ✓ QC图已保存:", fname, "\n")
}

# ==================== 4. 辅助：提取合并两组细胞 ====================
extract_and_merge <- function(obj1, obj2,
                               cell_type   = NULL,
                               origin      = NULL,
                               group1_name = NULL,
                               group2_name = NULL) {
  get_cells <- function(obj, ct) {
    if (!is.null(ct) && "final_celltype" %in% colnames(obj@meta.data))
      WhichCells(obj, expression = final_celltype == ct)
    else colnames(obj)
  }
  cells1 <- get_cells(obj1, cell_type)
  cells2 <- get_cells(obj2, cell_type)

  # 筛选人鼠来源（仅对肿瘤样本）
  if (!is.null(origin) && "cell_origin" %in% colnames(obj1@meta.data))
    cells1 <- intersect(cells1, WhichCells(obj1, expression = cell_origin == origin))

  if (length(cells1) < MIN_CELLS || length(cells2) < MIN_CELLS) return(NULL)

  s1 <- subset(obj1, cells = cells1); s1$group <- group1_name
  s2 <- subset(obj2, cells = cells2); s2$group <- group2_name
  merged <- merge(s1, s2)
  merged <- JoinLayers(merged)
  return(merged)
}

# ==================== 5. 核心差异分析（limma/voom） ====================
# 符合标准v5.0：
#   - 优先adj.P.Val；备用P.Value（比最低标准更低需申请，此处自动降级）
#   - DEG < 100时记录警告（标准要求反馈主管确认）
#   - 分组顺序可通过makeContrasts结果验证（+1=Case，-1=Control）
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
    cat(sprintf("  ⚠️  细胞数不足 (G1=%d, G2=%d)，跳过\n", n1, n2)); return(NULL)
  }

  expr    <- GetAssayData(merged_obj, layer = "data")
  keep    <- rowSums(expr > 0) >= (MIN_PCT * ncol(expr))
  expr_f  <- expr[keep, , drop = FALSE]
  if (nrow(expr_f) < 100) { cat("  ⚠️  过滤后基因<100，跳过\n"); return(NULL) }

  # 构建设计矩阵（标准v5.0：makeContrasts +1=Case, -1=Control）
  gv     <- factor(merged_obj$group, levels = c(g1, g2))
  design <- model.matrix(~ 0 + gv); colnames(design) <- levels(gv)
  cont   <- makeContrasts(contrasts = paste(g1, "-", g2), levels = design)
  # 打印分组方向确认（标准v5.0要求）
  cat("    分组方向（+1=Case,-1=Control）:", g1, "vs", g2, "\n")

  dge  <- DGEList(counts = expr_f)
  dge  <- calcNormFactors(dge)
  v    <- voom(dge, design, plot = FALSE)
  fit  <- lmFit(v, design)
  fit2 <- contrasts.fit(fit, cont)
  fit2 <- eBayes(fit2)

  res       <- topTable(fit2, number = Inf, adjust.method = "BH")
  res$gene  <- rownames(res)

  # 标记显著性（优先adj.P.Val）
  p_col   <- if (use_adj) "adj.P.Val" else "P.Value"
  res$sig <- ifelse(abs(res$logFC) > lfc_cut & res[[p_col]] < p_cut,
               ifelse(res$logFC > 0, "Up", "Down"), "NS")

  # 同时计算另一列便于报告
  res$pct.1 <- rowMeans(expr_f[res$gene, merged_obj$group == g1, drop = FALSE] > 0)
  res$pct.2 <- rowMeans(expr_f[res$gene, merged_obj$group == g2, drop = FALSE] > 0)

  n_up   <- sum(res$sig == "Up")
  n_down <- sum(res$sig == "Down")
  n_total<- n_up + n_down

  cat(sprintf("  ✓ %s | G1=%d G2=%d | Up=%d Down=%d Total=%d\n",
              comp_name, n1, n2, n_up, n_down, n_total))

  # 标准v5.0：DEG<100需要记录警告
  if (!use_alt && n_total < MIN_DEG_COUNT) {
    if (n_total < 10) {
      cat("    ⚠️  DEG<10，自动切换备用阈值（|log2FC|>0.5, P<0.05）...\n")
      return(run_limma_diff(merged_obj, comp_name, use_alt = TRUE))
    } else {
      cat(sprintf("    ⚠️  【注意】DEG=%d<100，按标准v5.0需反馈主管确认是否可行\n", n_total))
    }
  }

  return(res)
}

# ==================== 6. 火山图 ====================
plot_volcano <- function(res, comp_name, lfc_cut, p_cut, use_adj, save_dir) {
  p_col  <- if (use_adj) "adj.P.Val" else "P.Value"
  y_lab  <- if (use_adj) expression(-Log[10]~"adj.P.Val")
             else        expression(-Log[10]~"P.Value")

  res$log10p <- -log10(res[[p_col]] + 1e-300)
  res$log10p[is.infinite(res$log10p)] <-
    max(res$log10p[is.finite(res$log10p)], na.rm = TRUE)

  # top10 Up + top10 Down标签（按adj.P.Val排序）
  top_genes <- bind_rows(
    res %>% filter(sig == "Up")   %>% arrange(adj.P.Val) %>% head(10),
    res %>% filter(sig == "Down") %>% arrange(adj.P.Val) %>% head(10)
  )

  p <- ggplot(res, aes(x = logFC, y = log10p, color = sig)) +
    geom_point(alpha = 0.5, size = 0.7) +
    scale_color_manual(values = c("Up"="#E41A1C","Down"="#377EB8","NS"="grey70"),
                       name   = "DEG") +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut),
               linetype = "dashed", color = "black", linewidth = 0.4) +
    geom_hline(yintercept = -log10(p_cut),
               linetype = "dashed", color = "black", linewidth = 0.4) +
    geom_text_repel(data = top_genes, aes(label = gene),
                    size = 2.5, max.overlaps = 25,
                    box.padding = 0.3, show.legend = FALSE) +
    labs(title   = comp_name, x = expression(Log[2]~"Fold Change"), y = y_lab,
         caption = paste0("阈值: |log2FC|>", lfc_cut, ", ",
                          if(use_adj) paste0("adj.P.Val<",p_cut)
                          else        paste0("P.Value<",p_cut))) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, size = 11))

  ggsave(file.path(save_dir, "volcano", paste0(comp_name, "_volcano.pdf")),
         p, width = 8, height = 6, dpi = 300)
}

# ==================== 7. 热图 ====================
plot_heatmap <- function(res, merged_obj, comp_name, save_dir) {
  top_genes <- res %>%
    filter(sig != "NS") %>%
    arrange(desc(abs(logFC))) %>%
    head(30)
  if (nrow(top_genes) < 2) return(invisible(NULL))

  expr    <- GetAssayData(merged_obj, layer = "data")
  keep    <- rowSums(expr > 0) >= (MIN_PCT * ncol(expr))
  expr_f  <- as.matrix(expr[keep, ])
  genes_use <- intersect(top_genes$gene, rownames(expr_f))
  if (length(genes_use) < 2) return(invisible(NULL))

  expr_sub <- expr_f[genes_use, , drop = FALSE]
  expr_z   <- t(scale(t(expr_sub))); expr_z[is.na(expr_z)] <- 0

  # 标准v5.0：热图观察组内表现（检查离群样本）
  anno_df  <- data.frame(Group = merged_obj$group,
                          row.names = colnames(expr_z))
  groups   <- unique(merged_obj$group)
  anno_col <- list(Group = setNames(c("#E41A1C","#377EB8"), groups))

  # 细胞过多时按组等比例采样（每组最多100）
  if (ncol(expr_z) > 200) {
    set.seed(42)
    samp <- unlist(lapply(groups, function(g) {
      cs <- colnames(expr_z)[merged_obj$group[colnames(expr_z)] == g]
      sample(cs, min(100, length(cs)))
    }))
    expr_z  <- expr_z[, samp]
    anno_df <- anno_df[samp, , drop = FALSE]
  }

  pdf(file.path(save_dir, "heatmap", paste0(comp_name, "_heatmap.pdf")),
      width = 12, height = 10)
  pheatmap(expr_z,
           annotation_col    = anno_df,
           annotation_colors = anno_col,
           show_rownames      = TRUE, show_colnames = FALSE,
           cluster_rows       = TRUE, cluster_cols  = TRUE,
           fontsize_row       = 8,
           main               = paste0(comp_name, " Top", nrow(top_genes), " DEGs"),
           color              = colorRampPalette(c("#4575B4","white","#D73027"))(100))
  dev.off()
  cat("    ✓ 热图已保存\n")
}

# ==================== 8. ORA富集分析（KEGG + GO） ====================
# 符合标准v5.0：优先padj，基因数<10转GSEA
run_ora_enrichment <- function(res, comp_name, save_dir) {
  sig_genes <- res %>% filter(sig != "NS") %>% pull(gene)
  cat("    ORA输入基因数:", length(sig_genes), "\n")

  if (length(sig_genes) < ORA_MIN_GENES) {
    cat("    ⚠️  DEG<10，跳过ORA（标准v5.0：基因数<10转单基因GSEA）\n")
    return(invisible(NULL))
  }

  # 转换Entrez ID（鼠源用Mm）
  gene_df <- bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID",
                   OrgDb = org.Mm.eg.db)
  entrez_ids <- unique(gene_df$ENTREZID)
  if (length(entrez_ids) < ORA_MIN_GENES) {
    cat("    ⚠️  Entrez转换后基因<10，跳过ORA\n"); return(invisible(NULL))
  }

  # --- KEGG富集 ---
  kegg_res <- tryCatch(
    enrichKEGG(gene         = entrez_ids,
               organism     = "mmu",
               pAdjustMethod= "BH",
               pvalueCutoff = ORA_PADJ_CUTOFF,   # 优先padj（标准v5.0）
               qvalueCutoff = ORA_PADJ_CUTOFF),
    error = function(e) NULL)

  if (!is.null(kegg_res) && nrow(kegg_res@result) > 0) {
    write.csv(as.data.frame(kegg_res),
              file.path(deg_dir, paste0(comp_name, "_KEGG_ORA.csv")), row.names = FALSE)
    # 条形图
    p_kegg <- barplot(kegg_res, showCategory = 20,
                       title = paste0(comp_name, " KEGG ORA")) +
      theme(axis.text.y = element_text(size = 8))
    ggsave(file.path(save_dir, "enrichment", paste0(comp_name, "_KEGG_barplot.pdf")),
           p_kegg, width = 10, height = 8, dpi = 300)
    cat("    ✓ KEGG ORA完成，通路数:", nrow(kegg_res@result), "\n")
  } else {
    cat("    ⚠️  KEGG ORA无显著通路\n")
  }

  # --- GO富集（BP）---
  go_res <- tryCatch(
    enrichGO(gene          = entrez_ids,
              OrgDb         = org.Mm.eg.db,
              ont           = "BP",
              pAdjustMethod = "BH",
              pvalueCutoff  = ORA_PADJ_CUTOFF,
              qvalueCutoff  = ORA_PADJ_CUTOFF,
              readable      = TRUE),
    error = function(e) NULL)

  if (!is.null(go_res) && nrow(go_res@result) > 0) {
    write.csv(as.data.frame(go_res),
              file.path(deg_dir, paste0(comp_name, "_GO_BP_ORA.csv")), row.names = FALSE)
    p_go <- dotplot(go_res, showCategory = 20,
                     title = paste0(comp_name, " GO-BP ORA"))
    ggsave(file.path(save_dir, "enrichment", paste0(comp_name, "_GO_BP_dotplot.pdf")),
           p_go, width = 10, height = 8, dpi = 300)
    cat("    ✓ GO ORA完成，条目数:", nrow(go_res@result), "\n")
  } else {
    cat("    ⚠️  GO ORA无显著条目\n")
  }

  return(list(kegg = kegg_res, go = go_res))
}

# ==================== 9. GSEA分析 ====================
# 符合标准v5.0：
#   - 背景基因集优先KEGG
#   - 阈值：|NES|>1 & p<0.05（adjust.P<0.05）
#   - 使用seed参数避免随机性（clusterProfiler）
run_gsea <- function(res, comp_name, save_dir) {
  # 构建排序基因列表（按log2FC排序）
  gene_rank <- res %>%
    filter(!is.na(logFC)) %>%
    arrange(desc(logFC))

  gene_df <- bitr(gene_rank$gene, fromType = "SYMBOL", toType = "ENTREZID",
                   OrgDb = org.Mm.eg.db)
  gene_rank2 <- gene_rank %>%
    left_join(gene_df, by = c("gene" = "SYMBOL")) %>%
    filter(!is.na(ENTREZID)) %>%
    distinct(ENTREZID, .keep_all = TRUE)

  rank_vec <- setNames(gene_rank2$logFC, gene_rank2$ENTREZID)
  rank_vec <- sort(rank_vec, decreasing = TRUE)

  if (length(rank_vec) < 50) {
    cat("    ⚠️  GSEA排序基因<50，跳过\n"); return(invisible(NULL))
  }

  # GSEA-KEGG（背景基因集优先KEGG，标准v5.0）
  # seed参数避免随机性（标准v5.0要求）
  set.seed(2025)
  gsea_res <- tryCatch(
    gseKEGG(geneList      = rank_vec,
             organism      = "mmu",
             minGSSize     = 15,
             maxGSSize     = 500,
             pvalueCutoff  = GSEA_P_CUTOFF,
             pAdjustMethod = "BH",
             verbose       = FALSE,
             seed          = 2025),
    error = function(e) NULL)

  if (is.null(gsea_res) || nrow(gsea_res@result) == 0) {
    cat("    ⚠️  GSEA-KEGG无显著通路\n"); return(invisible(NULL))
  }

  # 筛选标准v5.0阈值：|NES|>1 & p<0.05
  gsea_df <- as.data.frame(gsea_res) %>%
    filter(abs(NES) > GSEA_NES_CUTOFF & pvalue < GSEA_P_CUTOFF)
  cat("    ✓ GSEA完成，|NES|>1通路数:", nrow(gsea_df), "\n")

  write.csv(gsea_df,
            file.path(deg_dir, paste0(comp_name, "_GSEA_KEGG.csv")), row.names = FALSE)

  # GSEA气泡图
  if (nrow(gsea_df) > 0) {
    gsea_plot_df <- gsea_df %>%
      arrange(desc(NES)) %>%
      head(30) %>%
      mutate(direction = ifelse(NES > 0, "Activated", "Suppressed"))

    p_gsea <- ggplot(gsea_plot_df,
                      aes(x = NES, y = reorder(Description, NES),
                          size = setSize, color = p.adjust)) +
      geom_point() +
      scale_color_gradient(low = "#E41A1C", high = "#377EB8",
                            name = "adj.P.Val") +
      scale_size_continuous(name = "Gene Set Size") +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.5) +
      labs(title = paste0(comp_name, " GSEA-KEGG"),
           x     = "NES", y = "",
           caption = "阈值: |NES|>1 & adj.P.Val<0.05") +
      theme_bw(base_size = 10) +
      theme(plot.title = element_text(hjust = 0.5))

    ggsave(file.path(save_dir, "gsea", paste0(comp_name, "_GSEA_bubble.pdf")),
           p_gsea, width = 10, height = max(6, nrow(gsea_plot_df)*0.3), dpi = 300)
  }

  return(gsea_res)
}

# ==================== 10. 肿瘤特异性DEG筛选（排除对照组基础差异） ====================
# 用户核心需求：用肿瘤 vs 空白 的DEG，减去 对照 vs 空白 的基础差异
# 肿瘤特异性DEG = 在肿瘤组显著 AND 在对照组基础比较中不显著
extract_tumor_specific_deg <- function(tumor_deg, ctrl_deg, comp_name) {
  if (is.null(tumor_deg) || is.null(ctrl_deg)) return(NULL)

  # 对照组基础差异基因（需要排除的）
  baseline_genes <- ctrl_deg %>% filter(sig != "NS") %>% pull(gene)
  cat("    基线差异基因数（对照 vs 空白）:", length(baseline_genes), "\n")

  # 肿瘤特异性 = 肿瘤中显著 & 基线中不显著
  specific_deg <- tumor_deg %>%
    filter(sig != "NS") %>%
    filter(!gene %in% baseline_genes)

  cat("    肿瘤特异性DEG数:", nrow(specific_deg),
      "(Up:", sum(specific_deg$sig=="Up"),
      " Down:", sum(specific_deg$sig=="Down"), ")\n")

  write.csv(specific_deg,
            file.path(specific_dir, paste0(comp_name, "_tumor_specific_DEG.csv")),
            row.names = FALSE)
  return(specific_deg)
}

# ==================== 11. 肝脑通路交集分析 ====================
analyze_liver_brain_crosstalk <- function(all_gsea_results, all_results, save_dir) {
  cat("\n========== 肝脑通路交集分析 ==========\n")

  # 收集肝脏轴和脑轴的显著GSEA通路
  liver_pathways <- list(); brain_pathways <- list()

  for (name in names(all_gsea_results)) {
    gsea_obj <- all_gsea_results[[name]]
    if (is.null(gsea_obj)) next
    df <- as.data.frame(gsea_obj) %>%
      filter(abs(NES) > GSEA_NES_CUTOFF & pvalue < GSEA_P_CUTOFF)
    if (nrow(df) == 0) next

    # 根据比较名中的器官标签归类
    if (grepl("liver|bla_liver|con_liver", name, ignore.case = TRUE) &&
        !grepl("brain|nao|bla_brain|con_brain", name, ignore.case = TRUE)) {
      liver_pathways[[name]] <- df
    } else if (grepl("brain|nao|bla_brain|con_brain", name, ignore.case = TRUE)) {
      brain_pathways[[name]] <- df
    }
  }

  if (length(liver_pathways) == 0 || length(brain_pathways) == 0) {
    cat("  ⚠️  肝脏或脑轴GSEA结果不足，跳过交集分析\n"); return(invisible(NULL))
  }

  # 合并同一器官的通路（取并集后统计出现频次）
  liver_all <- bind_rows(liver_pathways) %>%
    group_by(ID, Description) %>%
    summarise(n_liver = n(), mean_NES_liver = mean(NES),
              min_padj_liver = min(p.adjust), .groups = "drop")

  brain_all <- bind_rows(brain_pathways) %>%
    group_by(ID, Description) %>%
    summarise(n_brain = n(), mean_NES_brain = mean(NES),
              min_padj_brain = min(p.adjust), .groups = "drop")

  # 交集通路
  shared <- inner_join(liver_all, brain_all, by = c("ID","Description")) %>%
    mutate(NES_concordant = sign(mean_NES_liver) == sign(mean_NES_brain)) %>%
    arrange(desc(n_liver + n_brain))

  cat("  肝脏轴显著通路数:", nrow(liver_all), "\n")
  cat("  脑轴显著通路数:  ", nrow(brain_all), "\n")
  cat("  共同通路数:      ", nrow(shared), "\n")
  cat("  方向一致通路数:  ", sum(shared$NES_concordant), "\n")

  write.csv(shared, file.path(save_dir, "liver_brain_crosstalk",
                               "liver_brain_shared_pathways.csv"), row.names = FALSE)
  write.csv(liver_all, file.path(save_dir, "liver_brain_crosstalk",
                                  "liver_axis_pathways.csv"), row.names = FALSE)
  write.csv(brain_all, file.path(save_dir, "liver_brain_crosstalk",
                                  "brain_axis_pathways.csv"), row.names = FALSE)

  # --- 可视化1：韦恩图式气泡图（肝 vs 脑 NES对比） ---
  if (nrow(shared) > 0) {
    plot_df <- shared %>% head(30) %>%
      mutate(direction = case_when(
        mean_NES_liver > 0 & mean_NES_brain > 0 ~ "Both_Up",
        mean_NES_liver < 0 & mean_NES_brain < 0 ~ "Both_Down",
        TRUE                                     ~ "Discordant"
      ))

    p_cross <- ggplot(plot_df,
                       aes(x = mean_NES_liver, y = mean_NES_brain,
                           size = n_liver + n_brain,
                           color = direction,
                           label = Description)) +
      geom_point(alpha = 0.8) +
      geom_text_repel(size = 2.5, max.overlaps = 20,
                       box.padding = 0.3, show.legend = FALSE) +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
      geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
      scale_color_manual(values = c("Both_Up"="#E41A1C","Both_Down"="#377EB8",
                                     "Discordant"="grey50")) +
      scale_size_continuous(range = c(3, 10), name = "频次和") +
      labs(title   = "肝脑共同富集通路（GSEA-KEGG）",
           x       = "肝脏轴 平均NES",
           y       = "脑轴 平均NES",
           caption = "气泡大小=在肝脑比较中出现的总频次；颜色=富集方向") +
      theme_bw(base_size = 11) +
      theme(plot.title = element_text(hjust = 0.5, size = 12))

    ggsave(file.path(save_dir, "liver_brain_crosstalk",
                      "liver_brain_NES_scatter.pdf"),
           p_cross, width = 10, height = 8, dpi = 300)

    # --- 可视化2：共同通路热图（按器官分列NES） ---
    if (nrow(shared) >= 2) {
      hm_data <- shared %>%
        head(30) %>%
        select(Description, mean_NES_liver, mean_NES_brain) %>%
        tibble::column_to_rownames("Description")

      pdf(file.path(save_dir, "liver_brain_crosstalk",
                     "liver_brain_shared_heatmap.pdf"), width = 8, height = 10)
      pheatmap(hm_data,
               color           = colorRampPalette(c("#4575B4","white","#D73027"))(100),
               cluster_cols    = FALSE,
               fontsize_row    = 8,
               main            = "肝脑共同通路 NES热图",
               labels_col      = c("肝脏轴 NES","脑轴 NES"))
      dev.off()
    }
  }

  cat("  ✓ 肝脑交集分析完成\n")
  return(shared)
}

# ==================== 12. 主分析流程 ====================
cat("\n========== 开始主分析流程 ==========\n\n")

all_results        <- list()
all_gsea_results   <- list()
all_enrichment     <- list()
successful_comps   <- character(0)

# 分析所有比较组（按细胞类型分层）
run_full_analysis <- function(comp, seurat_objs, origin = NULL, origin_tag = "") {
  obj1 <- seurat_objs[[comp$group1]]
  obj2 <- seurat_objs[[comp$group2]]
  if (is.null(obj1) || is.null(obj2)) return(list())

  results_list <- list()

  for (ct in cell_types) {
    # 检查细胞类型是否存在于两组
    has_ct1 <- "final_celltype" %in% colnames(obj1@meta.data) &&
               ct %in% unique(obj1$final_celltype)
    has_ct2 <- "final_celltype" %in% colnames(obj2@meta.data) &&
               ct %in% unique(obj2$final_celltype)
    if (!has_ct1 || !has_ct2) next

    comp_name <- paste0(comp$name, origin_tag, "_", ct)
    cat("\n  >>", comp_name, "\n")

    # QC图
    plot_qc_for_comparison(obj1, obj2, comp_name = comp$name,
                            cell_type = ct)

    # 提取合并
    merged <- tryCatch(
      extract_and_merge(obj1, obj2,
                         cell_type   = ct,
                         origin      = origin,
                         group1_name = comp$group1,
                         group2_name = comp$group2),
      error = function(e) { cat("  ✗ 合并出错:", e$message, "\n"); NULL })
    if (is.null(merged)) next

    # 差异分析
    res <- run_limma_diff(merged, comp_name)
    if (is.null(res)) next

    # 保存全量DEG
    write.csv(res, file.path(deg_dir, paste0(comp_name, "_all_genes.csv")),
              row.names = FALSE)
    sig_df <- res %>% filter(sig != "NS")
    if (nrow(sig_df) > 0)
      write.csv(sig_df, file.path(deg_dir, paste0(comp_name, "_significant_genes.csv")),
                row.names = FALSE)

    # 判断实际使用的阈值（查看sig列是否用备用阈值）
    used_adj <- any(res$sig != "NS" & res$adj.P.Val < ADJ_P_CUTOFF)
    lfc_used <- if (any(abs(res$logFC[res$sig!="NS"]) >= LOG2FC_CUTOFF)) LOG2FC_CUTOFF else LOG2FC_ALT
    p_used   <- if (used_adj) ADJ_P_CUTOFF else P_ALT

    # 火山图
    plot_volcano(res, comp_name, lfc_used, p_used, used_adj, plot_dir)
    # 热图
    plot_heatmap(res, merged, comp_name, plot_dir)
    # ORA富集
    enrich_res <- run_ora_enrichment(res, comp_name, plot_dir)
    # GSEA
    gsea_res   <- run_gsea(res, comp_name, plot_dir)

    results_list[[comp_name]] <- list(
      deg     = res,
      enrich  = enrich_res,
      gsea    = gsea_res,
      organ   = comp$organ,
      type    = comp$type
    )
  }
  return(results_list)
}

# ---------- A. 肿瘤 vs 空白（鼠源细胞）----------
cat("\n===== A. 肿瘤 vs 空白 =====\n")
for (comp in tumor_vs_blank) {
  cat("\n--- ", comp$name, "---\n")
  obj1 <- seurat_objs[[comp$group1]]
  obj2 <- seurat_objs[[comp$group2]]
  if (is.null(obj1) || is.null(obj2)) { cat("样本未加载\n"); next }

  # 肿瘤样本只用鼠源细胞（与纯鼠对照可比）
  origin_tag <- if ("cell_origin" %in% colnames(obj1@meta.data)) "_mouse" else ""
  origin_val <- if ("cell_origin" %in% colnames(obj1@meta.data)) "Mouse"  else NULL

  res_list <- run_full_analysis(comp, seurat_objs,
                                 origin = origin_val, origin_tag = origin_tag)
  for (nm in names(res_list)) {
    all_results[[nm]]      <- res_list[[nm]]$deg
    all_gsea_results[[nm]] <- res_list[[nm]]$gsea
    all_enrichment[[nm]]   <- res_list[[nm]]$enrich
    successful_comps       <- c(successful_comps, nm)
  }
}

# ---------- B. 对照 vs 空白（基础差异，用于后续排除） ----------
cat("\n\n===== B. 对照 vs 空白（基础差异） =====\n")
ctrl_blank_results <- list()
for (comp in ctrl_vs_blank) {
  cat("\n---", comp$name, "---\n")
  res_list <- run_full_analysis(comp, seurat_objs)
  for (nm in names(res_list)) {
    ctrl_blank_results[[nm]] <- res_list[[nm]]$deg
    all_results[[nm]]        <- res_list[[nm]]$deg
    successful_comps         <- c(successful_comps, nm)
  }
}

# ---------- C. 肿瘤特异性DEG筛选 ----------
cat("\n\n===== C. 肿瘤特异性DEG（排除对照基础差异）=====\n")
# 肝脏：gan56_vs_bla_liver 减去 con_liver_vs_bla_liver
# 脑轴：gan56_vs_nao56 减去 con_brain_vs_bla_brain（同类型比较）

organ_pairs <- list(
  list(
    tumor_pattern = "gan56_vs_bla_liver_mouse_",
    ctrl_pattern  = "con_liver_vs_bla_liver_",
    organ         = "liver"
  ),
  list(
    tumor_pattern = "gan56_vs_nao56_mouse_",
    ctrl_pattern  = "con_brain_vs_bla_brain_",
    organ         = "brain"
  )
)

tumor_specific_all <- list()
for (pair in organ_pairs) {
  tumor_names <- grep(pair$tumor_pattern, names(all_results), value = TRUE, fixed = TRUE)
  ctrl_names  <- grep(pair$ctrl_pattern,  names(all_results), value = TRUE, fixed = TRUE)

  for (tn in tumor_names) {
    ct <- sub(pair$tumor_pattern, "", tn)
    cn <- paste0(pair$ctrl_pattern, ct)
    if (!cn %in% names(all_results)) next

    cat("\n  肿瘤特异性:", tn, "\n  基线比较:  ", cn, "\n")
    specific <- extract_tumor_specific_deg(
      tumor_deg = all_results[[tn]],
      ctrl_deg  = all_results[[cn]],
      comp_name = paste0(pair$organ, "_", ct, "_tumor_specific")
    )
    if (!is.null(specific))
      tumor_specific_all[[paste0(pair$organ, "_", ct)]] <- specific
  }
}

# ---------- D. 其他比较（扩展） ----------
cat("\n\n===== D. 其他组间比较 =====\n")
for (comp in other_comps) {
  cat("\n---", comp$name, "---\n")
  res_list <- run_full_analysis(comp, seurat_objs)
  for (nm in names(res_list)) {
    all_results[[nm]]      <- res_list[[nm]]$deg
    all_gsea_results[[nm]] <- res_list[[nm]]$gsea
    successful_comps       <- c(successful_comps, nm)
  }
}

# ---------- E. 肝脑通路交集分析 ----------
crosstalk_res <- analyze_liver_brain_crosstalk(all_gsea_results, all_results, plot_dir)

# ==================== 13. 汇总统计 ====================
cat("\n\n========== 汇总统计 ==========\n")

summary_df <- data.frame()
for (nm in names(all_results)) {
  d <- all_results[[nm]]
  if (is.null(d)) next
  summary_df <- rbind(summary_df, data.frame(
    Comparison = nm,
    Up         = sum(d$sig == "Up",   na.rm = TRUE),
    Down       = sum(d$sig == "Down", na.rm = TRUE),
    Total      = sum(d$sig != "NS",   na.rm = TRUE)
  ))
}
write.csv(summary_df, file.path(output_dir, "diff_analysis_summary.csv"),
          row.names = FALSE)

# 汇总柱状图
if (nrow(summary_df) > 0) {
  p_sum <- ggplot(summary_df %>% arrange(Total) %>%
                    tail(40) %>%
                    tidyr::pivot_longer(cols = c("Up","Down"),
                                        names_to = "direction", values_to = "count"),
                   aes(x = reorder(Comparison, count), y = count, fill = direction)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = c("Up"="#E41A1C","Down"="#377EB8")) +
    coord_flip() +
    labs(title = "差异基因数量汇总", x = "", y = "DEG数") +
    theme_bw(base_size = 9) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(file.path(plot_dir, "summary_barplot.pdf"),
         p_sum, width = 14, height = max(8, nrow(summary_df)*0.25), dpi = 300)
}

# ==================== 14. 分析报告 ====================
report_file <- file.path(output_dir, "limma_analysis_report.txt")
sink(report_file)
cat("========================================\n")
cat("limma 差异分析报告 v3（标准v5.0）\n")
cat("========================================\n")
cat("分析时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("【主要修改点】\n")
cat("  1. 人鼠分类：基因命名规则法（全大写=Human，首字母大写+小写=Mouse）\n")
cat("     比例阈值：human_ratio>=0.6→Human，<=0.4→Mouse，中间→Unknown\n")
cat("  2. 新增QC图：PCA图 + 箱线图（检查批次和离群样本，符合标准v5.0）\n")
cat("  3. 新增肿瘤特异性DEG筛选（排除对照组基础差异）\n")
cat("  4. 新增GSEA分析（|NES|>1 & adj.P.Val<0.05，背景KEGG，set.seed固定）\n")
cat("  5. 新增肝脑交集通路分析与可视化\n\n")

cat("【分析参数（符合标准v5.0）】\n")
cat("  主阈值:   |log2FC|>", LOG2FC_CUTOFF, " & adj.P.Val<", ADJ_P_CUTOFF, "\n")
cat("  备用阈值: |log2FC|>", LOG2FC_ALT, " & P.Value<", P_ALT, "\n")
cat("  DEG数量下限:", MIN_DEG_COUNT, "（低于需反馈主管）\n")
cat("  富集分析: 优先padj，基因数≥10，KEGG背景\n")
cat("  GSEA阈值: |NES|>1 & adj.P.Val<0.05\n")
cat("  相关性: spearman，最低阈值0.3\n\n")

cat("【成功完成的比较】\n")
for (nm in successful_comps) cat("  ✓", nm, "\n")

cat("\n【DEG数量统计】\n")
for (i in seq_len(nrow(summary_df))) {
  flag <- if (summary_df$Total[i] < MIN_DEG_COUNT) " ⚠️ <100需确认" else ""
  cat(sprintf("  %-70s Up=%-4d Down=%-4d Total=%-4d%s\n",
              summary_df$Comparison[i],
              summary_df$Up[i], summary_df$Down[i], summary_df$Total[i], flag))
}
sink()

# ==================== 完成 ====================
cat("\n========================================\n")
cat("分析完成！\n")
cat("成功比较数:", length(successful_comps), "\n")
cat("输出目录:\n")
cat("  DEG表格:    ", deg_dir, "\n")
cat("  肿瘤特异性: ", specific_dir, "\n")
cat("  QC图:       ", file.path(plot_dir, "qc"), "\n")
cat("  火山图:     ", file.path(plot_dir, "volcano"), "\n")
cat("  热图:       ", file.path(plot_dir, "heatmap"), "\n")
cat("  富集分析:   ", file.path(plot_dir, "enrichment"), "\n")
cat("  GSEA:       ", file.path(plot_dir, "gsea"), "\n")
cat("  肝脑交集:   ", file.path(plot_dir, "liver_brain_crosstalk"), "\n")
cat("  分析报告:   ", report_file, "\n")
cat("========================================\n")
