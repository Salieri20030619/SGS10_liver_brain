#!/usr/bin/env Rscript
# ============================================================
# diff_v4_viz_enrich.R  --  阶段 2: 可视化 + 富集 + 机制验证
# ============================================================
# 用法:
#   cd ~/SGS10_analysis
#   Rscript scripts/diff_v4_viz_enrich.R 2>&1 | tee scripts/logs/diff_viz.log
#
# 前置: 阶段 1 (diff_v4_core.R) 已完成, diff_stage1_results/ 下有 DEG CSV
#
# 7 个模块 (均 try-catch 包裹, 单个失败不影响其他):
#   Module 1: 火山图
#   Module 2: Top DEG 热图
#   Module 3: KEGG 富集 + dotplot
#   Module 4: GO BP 富集
#   Module 5: 9 机制模块跨对比热图 (核心产出)
#   Module 6: Il6st 专项可视化 (Violin + pct 条形图)
#   Module 7: 肝脑 KEGG Venn + shared pathway
#
# 主分析版本: v4_refined
# v32 和 v4_broad 只在 Module 6 机制热图里作为参考列对比
#
# 输出目录结构:
#   diff_stage1_results/
#     viz/
#       01_volcano/{liver,brain}/<key>_volcano.{pdf,png}
#       02_heatmap/{liver,brain}/<key>_heatmap.{pdf,png}
#       03_enrichment_kegg/{liver,brain}/<key>_KEGG.csv + dotplot
#       04_enrichment_go/{liver,brain}/<key>_GO_BP.csv
#       05_mechanism/mechanism_heatmap.{pdf,png} + mechanism_FC_matrix.csv
#       06_il6st_panel/Il6st_violin.{pdf,png}, Il6st_pct_bar.{pdf,png}
#       07_kegg_venn/liver_brain_venn.{pdf,png}, shared_pathways.csv
# ============================================================

suppressPackageStartupMessages({
  library(igraph)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(ComplexHeatmap)
  library(circlize)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
  library(tidyr)
  library(VennDiagram)
  library(grid)
  library(parallel)
})

# 锁定 dplyr (AnnotationDbi::select 冲突)
select    <- dplyr::select
rename    <- dplyr::rename
filter    <- dplyr::filter
mutate    <- dplyr::mutate
arrange   <- dplyr::arrange
summarise <- dplyr::summarise
group_by  <- dplyr::group_by
bind_rows <- dplyr::bind_rows
distinct  <- dplyr::distinct
left_join <- dplyr::left_join
inner_join<- dplyr::inner_join
pull      <- dplyr::pull

# ============================================================
# 0. 路径 & 配置
# ============================================================
PROJECT_ROOT <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
SCRIPT_DIR   <- file.path(PROJECT_ROOT, "scripts")
setwd(PROJECT_ROOT)

source(file.path(SCRIPT_DIR, "diff_common.R"))

STAGE1_DIR <- file.path(PROJECT_ROOT, "diff_stage1_results")
VIZ_ROOT   <- file.path(STAGE1_DIR, "viz")

VIZ_DIRS <- c(
  "01_volcano/liver", "01_volcano/brain",
  "02_heatmap/liver", "02_heatmap/brain",
  "03_enrichment_kegg/liver", "03_enrichment_kegg/brain",
  "04_enrichment_go/liver",   "04_enrichment_go/brain",
  "05_mechanism",
  "06_il6st_panel",
  "07_kegg_venn"
)
for (d in VIZ_DIRS) {
  dir.create(file.path(VIZ_ROOT, d), recursive = TRUE, showWarnings = FALSE)
}

VIZ_LOG <- file.path(STAGE1_DIR, "summary", "viz_error_log.txt")
writeLines(sprintf("== viz error log - %s ==", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
           VIZ_LOG)

log_viz_error <- function(module, context, msg) {
  line <- sprintf("[%s] [%s] %s => %s",
                  format(Sys.time(), "%H:%M:%S"), module, context, msg)
  cat("  ⚠", line, "\n")
  cat(line, "\n", file = VIZ_LOG, append = TRUE)
}

FC_THRESHOLD <- DEFAULT_PARAMS$FC_THRESHOLD
PADJ_CUT     <- DEFAULT_PARAMS$PADJ_CUT

cat("========================================\n")
cat("Stage 2: Visualization + Enrichment + Mechanism\n")
cat("Start:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================\n\n")

# ============================================================
# 1. 读取阶段 1 的所有 DEG CSV
# ============================================================
cat("========== Loading stage 1 DEG CSVs ==========\n")

load_deg_csvs <- function(version_dir, organ) {
  csv_dir <- file.path(STAGE1_DIR, version_dir, organ)
  if (!dir.exists(csv_dir)) return(list())
  files <- list.files(csv_dir, pattern = "\\.csv$", full.names = TRUE)
  out <- list()
  for (f in files) {
    key <- tools::file_path_sans_ext(basename(f))
    df <- tryCatch(read.csv(f, stringsAsFactors = FALSE),
                   error = function(e) NULL)
    if (!is.null(df) && nrow(df) > 0) out[[key]] <- df
  }
  out
}

deg_v4r_liver <- load_deg_csvs("deg_v4_refined", "liver")
deg_v4r_brain <- load_deg_csvs("deg_v4_refined", "brain")
deg_v4b_liver <- load_deg_csvs("deg_v4_broad", "liver")
deg_v4b_brain <- load_deg_csvs("deg_v4_broad", "brain")
deg_v32_liver <- load_deg_csvs("deg_v32", "liver")
deg_v32_brain <- load_deg_csvs("deg_v32", "brain")

cat(sprintf("  v4_refined: liver=%d brain=%d\n",
            length(deg_v4r_liver), length(deg_v4r_brain)))
cat(sprintf("  v4_broad  : liver=%d brain=%d\n",
            length(deg_v4b_liver), length(deg_v4b_brain)))
cat(sprintf("  v32       : liver=%d brain=%d\n\n",
            length(deg_v32_liver), length(deg_v32_brain)))

# 主分析用 v4_refined
deg_main_liver <- deg_v4r_liver
deg_main_brain <- deg_v4r_brain

# ============================================================
# 2. Module 1: 火山图
# ============================================================
cat("========== Module 1: Volcano plots ==========\n")

IL6_PATHWAY_FOR_VOLCANO <- c("Il6st", "Il6ra", "Il6", "Jak1", "Jak2",
                             "Stat3", "Socs3")

save_volcano <- function(res, key, organ) {
  tryCatch({
    rv <- res
    rv$log10p <- -log10(rv$p_val_adj + 1e-300)
    tg <- dplyr::bind_rows(
      rv %>% dplyr::filter(sig == "Up")   %>% dplyr::arrange(p_val_adj) %>% head(10),
      rv %>% dplyr::filter(sig == "Down") %>% dplyr::arrange(p_val_adj) %>% head(10)
    )
    rv$is_il6 <- rv$gene %in% IL6_PATHWAY_FOR_VOLCANO

    pv <- ggplot(rv, aes(x = avg_log2FC, y = log10p, color = sig)) +
      geom_point(alpha = 0.5, size = 0.8) +
      scale_color_manual(values = c("Up" = "#E41A1C", "Down" = "#377EB8",
                                    "NS" = "grey75")) +
      geom_vline(xintercept = c(-FC_THRESHOLD, FC_THRESHOLD),
                 linetype = "dashed", linewidth = 0.4) +
      geom_hline(yintercept = -log10(PADJ_CUT),
                 linetype = "dashed", linewidth = 0.4) +
      geom_text_repel(data = tg, aes(label = gene),
                      size = 2.5, max.overlaps = 30, segment.size = 0.2) +
      geom_point(data = rv %>% dplyr::filter(is_il6 & sig != "NS"),
                 aes(x = avg_log2FC, y = log10p),
                 color = "#FF7F00", size = 2.5, shape = 17) +
      geom_text_repel(data = rv %>% dplyr::filter(is_il6),
                      aes(label = gene),
                      color = "#FF7F00", size = 3, fontface = "bold",
                      max.overlaps = Inf, segment.size = 0.3) +
      labs(title = key,
           subtitle = sprintf("FC>%.2f | padj<%.2f | Up=%d Down=%d | 橙色 ▲ = IL-6 通路",
                              FC_THRESHOLD, PADJ_CUT,
                              sum(rv$sig == "Up"),
                              sum(rv$sig == "Down")),
           x = "avg_log2FC", y = "-log10(p_val_adj)", color = NULL) +
      theme_bw(base_size = 11) +
      theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, size = 8, color = "grey40"))

    base <- file.path(VIZ_ROOT, "01_volcano", organ, paste0(key, "_volcano"))
    ggsave(paste0(base, ".pdf"), pv, width = 9, height = 7, dpi = 300)
    ggsave(paste0(base, ".png"), pv, width = 9, height = 7, dpi = 200)
  }, error = function(e) log_viz_error("volcano", key, e$message))
}

for (organ in c("liver", "brain")) {
  deg_list <- if (organ == "liver") deg_main_liver else deg_main_brain
  for (key in names(deg_list)) save_volcano(deg_list[[key]], key, organ)
}
cat(sprintf("  liver: %d brain: %d volcanoes\n\n",
            length(deg_main_liver), length(deg_main_brain)))

# ============================================================
# 3. Module 2: Top DEG 热图 (基于 FC matrix, 不需要 Seurat 对象)
# ============================================================
cat("========== Module 2: Top DEG heatmaps ==========\n")

# 用 FC matrix 画 top DEG 热图: 每个 key 的 top40 基因在当前 key 里的 FC
save_fc_heatmap <- function(res, key, organ) {
  tryCatch({
    sig_genes <- res %>%
      dplyr::filter(sig != "NS") %>%
      dplyr::arrange(dplyr::desc(abs(avg_log2FC))) %>%
      head(40)
    if (nrow(sig_genes) < 3) return(invisible(NULL))

    mat <- matrix(sig_genes$avg_log2FC, ncol = 1,
                  dimnames = list(sig_genes$gene, key))

    cf <- colorRamp2(c(-max(abs(mat)), 0, max(abs(mat))),
                     c("#4575B4", "white", "#D73027"))

    base <- file.path(VIZ_ROOT, "02_heatmap", organ, paste0(key, "_heatmap"))
    for (fmt in c("pdf", "png")) {
      h <- max(6, nrow(sig_genes) * 0.25 + 2)
      if (fmt == "pdf") pdf(paste0(base, ".pdf"), width = 5, height = h)
      else              png(paste0(base, ".png"), width = 1000,
                            height = max(1200, nrow(sig_genes) * 50 + 400),
                            res = 200)
      draw(Heatmap(mat,
                   name = "log2FC",
                   col = cf,
                   row_names_gp = gpar(fontsize = 8),
                   column_names_gp = gpar(fontsize = 9),
                   cluster_rows = FALSE,
                   cluster_columns = FALSE,
                   show_row_names = TRUE,
                   column_title = paste0("Top DEG: ", key),
                   column_title_gp = gpar(fontsize = 10, fontface = "bold")))
      dev.off()
    }
  }, error = function(e) log_viz_error("heatmap", key, e$message))
}

for (organ in c("liver", "brain")) {
  deg_list <- if (organ == "liver") deg_main_liver else deg_main_brain
  for (key in names(deg_list)) save_fc_heatmap(deg_list[[key]], key, organ)
}
cat("  done\n\n")

# ============================================================
# 4. Module 3 + 4: KEGG + GO BP 富集
# ============================================================
cat("========== Module 3+4: KEGG + GO enrichment ==========\n")

run_enrichment <- function(res, key, organ) {
  tryCatch({
    sig_genes <- res %>% dplyr::filter(sig != "NS") %>% dplyr::pull(gene)
    if (length(sig_genes) < 10) return(invisible(NULL))

    gene_df <- tryCatch(
      bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID",
           OrgDb = org.Mm.eg.db),
      error = function(e) NULL
    )
    if (is.null(gene_df) || nrow(gene_df) < 5) return(invisible(NULL))
    entrez <- unique(gene_df$ENTREZID)

    # KEGG
    kegg <- tryCatch(
      enrichKEGG(gene = entrez, organism = "mmu", pvalueCutoff = 0.05),
      error = function(e) NULL
    )
    if (!is.null(kegg) && nrow(kegg@result) > 0) {
      write.csv(as.data.frame(kegg),
                file.path(VIZ_ROOT, "03_enrichment_kegg", organ,
                          paste0(key, "_KEGG.csv")),
                row.names = FALSE)
      n_sig <- sum(kegg@result$p.adjust < 0.05)
      if (n_sig >= 5) {
        p_dot <- tryCatch(
          dotplot(kegg, showCategory = 20, title = paste0("KEGG: ", key)),
          error = function(e) NULL
        )
        if (!is.null(p_dot)) {
          bk <- file.path(VIZ_ROOT, "03_enrichment_kegg", organ,
                          paste0(key, "_KEGG_dotplot"))
          ggsave(paste0(bk, ".pdf"), p_dot, width = 10, height = 8, dpi = 300)
          ggsave(paste0(bk, ".png"), p_dot, width = 10, height = 8, dpi = 200)
        }
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
                file.path(VIZ_ROOT, "04_enrichment_go", organ,
                          paste0(key, "_GO_BP.csv")),
                row.names = FALSE)
    }

  }, error = function(e) log_viz_error("enrichment", key, e$message))
}

cat("  running enrichment for liver...\n")
for (key in names(deg_main_liver)) {
  run_enrichment(deg_main_liver[[key]], key, "liver")
}
cat("  running enrichment for brain...\n")
for (key in names(deg_main_brain)) {
  run_enrichment(deg_main_brain[[key]], key, "brain")
}
cat("  done\n\n")

# ============================================================
# 5. Module 5: 9 机制模块跨对比热图 (核心产出)
# ============================================================
cat("========== Module 5: Mechanism heatmap (9 modules) ==========\n")

# 构建 FC 矩阵: 行 = 机制基因, 列 = 对比_细类
build_mechanism_fc_matrix <- function(deg_list, genes) {
  mat_df <- data.frame(gene = genes, stringsAsFactors = FALSE)
  for (nm in names(deg_list)) {
    res <- deg_list[[nm]]
    if (!"gene" %in% colnames(res)) next
    tmp <- res %>%
      dplyr::filter(gene %in% genes) %>%
      dplyr::select(gene, avg_log2FC)
    if (nrow(tmp) == 0) next
    tmp <- tmp[!duplicated(tmp$gene), ]
    new_col <- setNames(tmp$avg_log2FC, tmp$gene)
    mat_df[[nm]] <- new_col[mat_df$gene]
  }
  rownames(mat_df) <- mat_df$gene
  m <- as.matrix(mat_df[, -1, drop = FALSE])
  m[is.na(m)] <- 0
  m
}

tryCatch({
  fc_liver <- build_mechanism_fc_matrix(deg_main_liver, ALL_MECH_GENES)
  fc_brain <- build_mechanism_fc_matrix(deg_main_brain, ALL_MECH_GENES)

  # 保存 raw FC 矩阵
  write.csv(
    data.frame(gene = rownames(fc_liver), fc_liver),
    file.path(VIZ_ROOT, "05_mechanism", "mechanism_FC_liver.csv"),
    row.names = FALSE
  )
  write.csv(
    data.frame(gene = rownames(fc_brain), fc_brain),
    file.path(VIZ_ROOT, "05_mechanism", "mechanism_FC_brain.csv"),
    row.names = FALSE
  )

  # 只保留至少一列 |FC| > 0.25 的基因 (过滤弱信号, 避免稀释强机制信号)
  # 注意: 这个阈值只影响基因级别热图的展示过滤, 不影响模块强度热图
  keep_liver <- apply(abs(fc_liver), 1, max) > 0.25
  keep_brain <- apply(abs(fc_brain), 1, max) > 0.25
  gene_show  <- union(rownames(fc_liver)[keep_liver],
                      rownames(fc_brain)[keep_brain])
  gene_show  <- gene_show[gene_show %in% ALL_MECH_GENES]

  if (length(gene_show) < 3) {
    stop("too few mechanism genes pass filter")
  }

  # 对齐行
  fc_liver_show <- fc_liver[gene_show, , drop = FALSE]
  fc_brain_show <- fc_brain[gene_show, , drop = FALSE]

  # 过滤掉 All_Cells 列 (解读价值低)
  fc_liver_show <- fc_liver_show[, !grepl("_All_Cells$", colnames(fc_liver_show)),
                                 drop = FALSE]
  fc_brain_show <- fc_brain_show[, !grepl("_All_Cells$", colnames(fc_brain_show)),
                                 drop = FALSE]

  fc_comb <- cbind(fc_liver_show, fc_brain_show)

  # 列注释: 器官
  col_anno <- data.frame(
    Organ = c(rep("Liver", ncol(fc_liver_show)),
              rep("Brain", ncol(fc_brain_show))),
    row.names = colnames(fc_comb)
  )

  # 行注释: 机制模块
  pw_map <- dplyr::bind_rows(lapply(names(MECHANISM_GENES), function(nm) {
    data.frame(gene = MECHANISM_GENES[[nm]], Pathway = nm,
               stringsAsFactors = FALSE)
  }))
  pw_map <- pw_map[!duplicated(pw_map$gene), ]
  row_pw <- pw_map$Pathway[match(gene_show, pw_map$gene)]
  row_anno <- data.frame(Pathway = row_pw, row.names = gene_show)

  # 配色
  pw_col <- c(
    IL6_JAK_STAT3 = "#E41A1C", TGFb_SMAD = "#FF7F00",
    NFkB          = "#984EA3", PI3K_AKT  = "#4DAF4A",
    MAPK          = "#377EB8", Complement_C3 = "#A65628",
    CXCL_CXCR     = "#F781BF", S100_RAGE = "#FFFF33",
    AcutePhase    = "#999999", Microglia_Activation = "#66C2A5"
  )
  cf_fc <- colorRamp2(c(-2, 0, 2), c("#377EB8", "white", "#E41A1C"))

  base <- file.path(VIZ_ROOT, "05_mechanism", "mechanism_heatmap")
  for (fmt in c("pdf", "png")) {
    if (fmt == "pdf") pdf(paste0(base, ".pdf"),
                          width = max(12, ncol(fc_comb) * 0.35 + 5),
                          height = max(8, length(gene_show) * 0.25 + 3))
    else              png(paste0(base, ".png"),
                          width = max(2400, ncol(fc_comb) * 90 + 800),
                          height = max(1400, length(gene_show) * 50 + 600),
                          res = 200)
    draw(Heatmap(
      fc_comb, name = "avg_log2FC",
      col = cf_fc,
      top_annotation = HeatmapAnnotation(
        df = col_anno,
        col = list(Organ = c("Liver" = "#E41A1C", "Brain" = "#377EB8")),
        annotation_name_side = "left"
      ),
      left_annotation = rowAnnotation(
        df = row_anno,
        col = list(Pathway = pw_col),
        width = unit(4, "mm")
      ),
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      show_row_names = TRUE,
      show_column_names = TRUE,
      row_names_side = "left",
      row_names_gp = gpar(fontsize = 8),
      column_names_gp = gpar(fontsize = 7),
      column_names_rot = 45,
      column_title = "9 机制模块: 肝脑跨器官 FC 面板 (v4_refined)",
      column_title_gp = gpar(fontsize = 12, fontface = "bold"),
      cell_fun = function(j, i, x, y, width, height, fill) {
        v <- fc_comb[i, j]
        if (abs(v) > FC_THRESHOLD)
          grid.text(sprintf("%.2f", v), x, y,
                    gp = gpar(fontsize = 6, col = "black"))
      }
    ))
    dev.off()
  }

  # 按机制模块聚合: 每个对比_细类下, 每个模块的平均 |FC| 作为"模块强度"
  module_strength <- function(mat, gene_list_byname) {
    # mat: gene x comparison_ct
    comps <- colnames(mat)
    out <- matrix(0, nrow = length(gene_list_byname), ncol = length(comps),
                  dimnames = list(names(gene_list_byname), comps))
    for (mod in names(gene_list_byname)) {
      genes <- intersect(gene_list_byname[[mod]], rownames(mat))
      if (length(genes) == 0) next
      sub <- mat[genes, , drop = FALSE]
      # 模块强度 = 该模块所有基因 |FC| 的均值
      out[mod, ] <- colMeans(abs(sub), na.rm = TRUE)
    }
    out
  }

  strength_liver <- module_strength(fc_liver_show, MECHANISM_GENES)
  strength_brain <- module_strength(fc_brain_show, MECHANISM_GENES)
  strength_comb  <- cbind(strength_liver, strength_brain)

  write.csv(
    data.frame(module = rownames(strength_comb), strength_comb),
    file.path(VIZ_ROOT, "05_mechanism", "mechanism_strength_matrix.csv"),
    row.names = FALSE
  )

  col_anno_s <- data.frame(
    Organ = c(rep("Liver", ncol(strength_liver)),
              rep("Brain", ncol(strength_brain))),
    row.names = colnames(strength_comb)
  )
  cf_s <- colorRamp2(c(0, 0.5, 1.5), c("white", "#FEB24C", "#B10026"))

  base2 <- file.path(VIZ_ROOT, "05_mechanism", "mechanism_strength_heatmap")
  for (fmt in c("pdf", "png")) {
    if (fmt == "pdf") pdf(paste0(base2, ".pdf"),
                          width = max(12, ncol(strength_comb) * 0.35 + 4),
                          height = 6)
    else              png(paste0(base2, ".png"),
                          width = max(2400, ncol(strength_comb) * 90 + 700),
                          height = 1200, res = 200)
    draw(Heatmap(
      strength_comb, name = "mean|FC|",
      col = cf_s,
      top_annotation = HeatmapAnnotation(
        df = col_anno_s,
        col = list(Organ = c("Liver" = "#E41A1C", "Brain" = "#377EB8"))
      ),
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      row_names_gp = gpar(fontsize = 10, fontface = "bold"),
      column_names_gp = gpar(fontsize = 7),
      column_names_rot = 45,
      column_title = "9 机制模块强度 (每模块基因 |FC| 均值)",
      column_title_gp = gpar(fontsize = 12, fontface = "bold"),
      cell_fun = function(j, i, x, y, width, height, fill) {
        v <- strength_comb[i, j]
        if (v > 0.1)
          grid.text(sprintf("%.2f", v), x, y,
                    gp = gpar(fontsize = 7, col = "black"))
      }
    ))
    dev.off()
  }

  cat(sprintf("  ✓ mechanism heatmap: %d genes x %d comparisons\n",
              length(gene_show), ncol(fc_comb)))
  cat(sprintf("  ✓ mechanism strength heatmap: %d modules\n",
              nrow(strength_comb)))
}, error = function(e) log_viz_error("mechanism", "main", e$message))

# ============================================================
# 6. Module 6: Il6st 专项可视化 (绕过 padj 阈值)
# ============================================================
cat("\n========== Module 6: Il6st dedicated visualization ==========\n")

tryCatch({
  # 从 v4 RDS 读取所有脑样本, 提取 Microglia_Homeostatic 和 Microglia_DAM 细胞的 Il6st 表达
  cat("  loading brain v4 RDS for Il6st expression extraction...\n")
  brain_samples <- BRAIN_SAMPLES

  il6st_long <- list()
  for (s in brain_samples) {
    fp <- resolve_rds_path(s, "v4")
    if (!file.exists(fp)) next
    obj <- tryCatch(readRDS(fp), error = function(e) NULL)
    if (is.null(obj)) next
    obj <- tryCatch(JoinLayers(obj), error = function(e) obj)

    if (!"Il6st" %in% rownames(obj)) {
      cat(sprintf("  [%s] Il6st not in rownames, skip\n", s))
      next
    }
    if (!"celltype_refined" %in% colnames(obj@meta.data)) {
      cat(sprintf("  [%s] celltype_refined not in metadata, skip\n", s))
      next
    }

    # 只保留感兴趣的细类
    target_cts <- c("Microglia_Homeostatic", "Microglia_DAM", "Endothelial")
    cell_mask <- obj$celltype_refined %in% target_cts
    if (sum(cell_mask) == 0) next

    expr <- GetAssayData(obj, layer = "data")["Il6st", cell_mask]
    meta <- obj@meta.data[cell_mask, , drop = FALSE]

    group_label <- if (s == "con_brain") "Control" else
                   if (s == "bla_brain") "Blank" else "Tumor"

    il6st_long[[s]] <- data.frame(
      sample   = s,
      group    = group_label,
      celltype = as.character(meta$celltype_refined),
      il6st    = as.numeric(expr),
      stringsAsFactors = FALSE
    )
    rm(obj); gc(verbose = FALSE)
    cat(sprintf("  [%s] extracted %d cells\n", s, nrow(il6st_long[[s]])))
  }

  if (length(il6st_long) > 0) {
    il6st_df <- dplyr::bind_rows(il6st_long)
    il6st_df$group <- factor(il6st_df$group, levels = c("Control", "Blank", "Tumor"))
    il6st_df$celltype <- factor(il6st_df$celltype,
                                 levels = c("Microglia_Homeostatic",
                                            "Microglia_DAM", "Endothelial"))

    write.csv(il6st_df,
              file.path(VIZ_ROOT, "06_il6st_panel", "Il6st_cell_level_expr.csv"),
              row.names = FALSE)

    # Violin plot
    p_violin <- ggplot(il6st_df, aes(x = group, y = il6st, fill = group)) +
      geom_violin(scale = "width", trim = FALSE) +
      geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
      facet_wrap(~celltype, scales = "free_y", ncol = 3) +
      scale_fill_manual(values = c(Control = "#4DAF4A",
                                    Blank = "#377EB8",
                                    Tumor = "#E41A1C")) +
      labs(title = "Il6st expression: Brain cells across groups",
           subtitle = "v4_refined annotation, all cells (no downsample)",
           y = "Normalized expression (log)", x = "") +
      theme_bw(base_size = 11) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
            strip.text = element_text(face = "bold"),
            legend.position = "bottom")

    ggsave(file.path(VIZ_ROOT, "06_il6st_panel", "Il6st_violin.pdf"),
           p_violin, width = 12, height = 6, dpi = 300)
    ggsave(file.path(VIZ_ROOT, "06_il6st_panel", "Il6st_violin.png"),
           p_violin, width = 12, height = 6, dpi = 200)

    # pct (expressor) 条形图
    pct_df <- il6st_df %>%
      dplyr::group_by(sample, group, celltype) %>%
      dplyr::summarise(
        n_total   = dplyr::n(),
        n_expr    = sum(il6st > 0),
        pct_expr  = 100 * n_expr / n_total,
        mean_expr = mean(il6st),
        .groups   = "drop"
      )

    write.csv(pct_df,
              file.path(VIZ_ROOT, "06_il6st_panel", "Il6st_pct_by_sample.csv"),
              row.names = FALSE)

    p_pct <- ggplot(pct_df, aes(x = sample, y = pct_expr, fill = group)) +
      geom_col() +
      facet_wrap(~celltype, scales = "free_y", ncol = 3) +
      scale_fill_manual(values = c(Control = "#4DAF4A",
                                    Blank = "#377EB8",
                                    Tumor = "#E41A1C")) +
      labs(title = "Il6st % expressor cells: Brain samples",
           subtitle = "pct = fraction of cells with Il6st > 0 (v4_refined)",
           y = "% cells expressing Il6st", x = "") +
      theme_bw(base_size = 11) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
            strip.text = element_text(face = "bold"),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            legend.position = "bottom")

    ggsave(file.path(VIZ_ROOT, "06_il6st_panel", "Il6st_pct_bar.pdf"),
           p_pct, width = 12, height = 6, dpi = 300)
    ggsave(file.path(VIZ_ROOT, "06_il6st_panel", "Il6st_pct_bar.png"),
           p_pct, width = 12, height = 6, dpi = 200)

    # 组间聚合 pct (PDX 3 个合并)
    pct_group_df <- il6st_df %>%
      dplyr::group_by(group, celltype) %>%
      dplyr::summarise(
        n_total   = dplyr::n(),
        n_expr    = sum(il6st > 0),
        pct_expr  = 100 * n_expr / n_total,
        mean_expr = mean(il6st),
        .groups   = "drop"
      )

    write.csv(pct_group_df,
              file.path(VIZ_ROOT, "06_il6st_panel", "Il6st_pct_by_group.csv"),
              row.names = FALSE)

    p_pct_g <- ggplot(pct_group_df,
                      aes(x = group, y = pct_expr, fill = group)) +
      geom_col(width = 0.6) +
      geom_text(aes(label = sprintf("%.1f%%", pct_expr)),
                vjust = -0.3, size = 3) +
      facet_wrap(~celltype, scales = "free_y", ncol = 3) +
      scale_fill_manual(values = c(Control = "#4DAF4A",
                                    Blank = "#377EB8",
                                    Tumor = "#E41A1C")) +
      labs(title = "Il6st % expressor: group-level aggregation",
           subtitle = "Tumor = 3 PDX samples pooled",
           y = "% cells expressing Il6st", x = "") +
      theme_bw(base_size = 12) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
            strip.text = element_text(face = "bold"),
            legend.position = "none")

    ggsave(file.path(VIZ_ROOT, "06_il6st_panel", "Il6st_pct_bar_grouped.pdf"),
           p_pct_g, width = 12, height = 5, dpi = 300)
    ggsave(file.path(VIZ_ROOT, "06_il6st_panel", "Il6st_pct_bar_grouped.png"),
           p_pct_g, width = 12, height = 5, dpi = 200)

    cat("  ✓ Il6st panel saved\n")
    cat("\n  【Il6st % expressor by group】\n")
    print(pct_group_df, row.names = FALSE)
  } else {
    cat("  ⚠ No Il6st data extracted\n")
  }
}, error = function(e) log_viz_error("il6st", "main", e$message))

# ============================================================
# 7. Module 7: 肝脑 KEGG Venn + shared pathways
# ============================================================
cat("\n========== Module 7: Liver-Brain shared KEGG ==========\n")

tryCatch({
  read_kegg_dir <- function(dir) {
    if (!dir.exists(dir)) return(list())
    files <- list.files(dir, pattern = "_KEGG\\.csv$", full.names = TRUE)
    out <- list()
    for (f in files) {
      key <- gsub("_KEGG\\.csv$", "", basename(f))
      df <- tryCatch(read.csv(f, stringsAsFactors = FALSE),
                     error = function(e) NULL)
      if (!is.null(df) && nrow(df) > 0) out[[key]] <- df
    }
    out
  }

  lk <- read_kegg_dir(file.path(VIZ_ROOT, "03_enrichment_kegg", "liver"))
  bk <- read_kegg_dir(file.path(VIZ_ROOT, "03_enrichment_kegg", "brain"))
  cat(sprintf("  liver KEGG: %d sets, brain KEGG: %d sets\n",
              length(lk), length(bk)))

  if (length(lk) > 0 && length(bk) > 0) {
    # 只看 tumor_vs_control 的
    lk_tc <- lk[grepl("tumor_vs_control", names(lk))]
    bk_tc <- bk[grepl("tumor_vs_control", names(bk))]

    if (length(lk_tc) > 0 && length(bk_tc) > 0) {
      pool <- function(kl, tag) dplyr::bind_rows(
        lapply(names(kl), function(nm) {
          df <- kl[[nm]]
          df$source <- nm
          df$organ  <- tag
          df
        })
      )
      lp <- pool(lk_tc, "Liver")
      bp <- pool(bk_tc, "Brain")
      li <- unique(lp$ID)
      bi <- unique(bp$ID)
      si <- intersect(li, bi)
      cat(sprintf("  liver pathways: %d, brain pathways: %d, shared: %d\n",
                  length(li), length(bi), length(si)))

      # Venn
      if (length(li) > 0 || length(bi) > 0) {
        vd <- venn.diagram(
          list(Liver = li, Brain = bi),
          filename = NULL,
          fill = c("#E41A1C", "#377EB8"),
          alpha = 0.5, cex = 1.5, cat.cex = 1.5, cat.pos = c(-20, 20),
          main = "Shared KEGG pathways: tumor vs control",
          main.cex = 1.2
        )
        vb <- file.path(VIZ_ROOT, "07_kegg_venn", "liver_brain_venn")
        pdf(paste0(vb, ".pdf"), width = 7, height = 6)
        grid.draw(vd); dev.off()
        png(paste0(vb, ".png"), width = 1400, height = 1200, res = 200)
        grid.draw(vd); dev.off()
      }

      # Shared pathway 排名
      if (length(si) > 0) {
        ls_s <- lp %>%
          dplyr::filter(ID %in% si) %>%
          dplyr::group_by(ID, Description) %>%
          dplyr::summarise(Liver_minP = min(p.adjust), .groups = "drop")
        bs_s <- bp %>%
          dplyr::filter(ID %in% si) %>%
          dplyr::group_by(ID, Description) %>%
          dplyr::summarise(Brain_minP = min(p.adjust), .groups = "drop")

        shared_df <- dplyr::inner_join(ls_s, bs_s, by = c("ID", "Description")) %>%
          dplyr::mutate(combined_score = -log10(Liver_minP) - log10(Brain_minP)) %>%
          dplyr::arrange(dplyr::desc(combined_score))

        write.csv(shared_df,
                  file.path(VIZ_ROOT, "07_kegg_venn", "shared_pathways.csv"),
                  row.names = FALSE)

        # Scatter
        if (nrow(shared_df) >= 3) {
          pd <- shared_df %>% head(25) %>%
            dplyr::mutate(
              Liver_logP = -log10(Liver_minP),
              Brain_logP = -log10(Brain_minP),
              Desc_short = ifelse(nchar(Description) > 45,
                                  paste0(substr(Description, 1, 42), "..."),
                                  Description)
            )

          p_sc <- ggplot(pd, aes(x = Liver_logP, y = Brain_logP,
                                  label = Desc_short)) +
            geom_point(aes(size = combined_score, color = combined_score),
                       alpha = 0.8) +
            scale_color_gradientn(
              colors = c("#4575B4", "#91BFDB", "#FFFFBF",
                         "#FC8D59", "#D73027")
            ) +
            geom_text_repel(size = 2.8, max.overlaps = 20,
                            segment.size = 0.2) +
            labs(title = "Liver-Brain shared KEGG pathways",
                 subtitle = "Tumor vs Control (both organs)",
                 x = "-log10 padj (Liver)",
                 y = "-log10 padj (Brain)") +
            theme_bw(base_size = 12) +
            theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                  plot.subtitle = element_text(hjust = 0.5, color = "grey40"))

          ggsave(file.path(VIZ_ROOT, "07_kegg_venn", "shared_scatter.pdf"),
                 p_sc, width = 10, height = 8, dpi = 300)
          ggsave(file.path(VIZ_ROOT, "07_kegg_venn", "shared_scatter.png"),
                 p_sc, width = 10, height = 8, dpi = 200)
        }

        cat("\n  【Shared pathway top 15】\n")
        top15 <- shared_df %>% head(15) %>%
          dplyr::select(ID, Description, Liver_minP, Brain_minP, combined_score)
        print(top15, row.names = FALSE)
      }
    } else {
      cat("  (no tumor_vs_control KEGG data to compare)\n")
    }
  }
}, error = function(e) log_viz_error("kegg_venn", "main", e$message))

# ============================================================
# 8. 完成
# ============================================================
cat("\n========================================\n")
cat("Stage 2 complete\n")
cat(sprintf("  Output dir: %s\n", VIZ_ROOT))
cat(sprintf("  Error log : %s\n", VIZ_LOG))
cat("========================================\n")
