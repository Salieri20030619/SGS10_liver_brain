#!/usr/bin/env Rscript
# Stage 1 富集分析: GO + KEGG for brain clean DEG
# 输出: Top 15 柱形图 + 点图 + 富集表 CSV

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(dplyr)
})

PROJECT_ROOT <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
setwd(PROJECT_ROOT)

DEG_DIR <- file.path(PROJECT_ROOT, "diff_stage1_brain_clean", "deg_v4_refined")
OUT_DIR <- file.path(PROJECT_ROOT, "diff_stage1_brain_clean", "enrichment")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# 配置
PADJ_CUT <- 0.05
LOGFC_CUT <- 0.25
TOP_N_TERM <- 15

# 排除的 cell types
EXCLUDE_CT <- c("low_quality","Unknown","Other","Microglia_DAM_strict")

# 遍历 3 个 contrasts
contrasts <- list.dirs(DEG_DIR, recursive = FALSE)
if (length(contrasts) == 0) {
  stop(sprintf("No contrast dirs found in %s", DEG_DIR))
}

cat("========================================\n")
cat("Brain Stage 1 enrichment (GO + KEGG)\n")
cat("========================================\n\n")

summary_rows <- list()

for (ct_dir in contrasts) {
  contrast_name <- basename(ct_dir)
  cat(sprintf("\n\n### Contrast: %s ###\n", contrast_name))
  
  # per-contrast output dirs
  go_dir <- file.path(OUT_DIR, contrast_name, "GO")
  kegg_dir <- file.path(OUT_DIR, contrast_name, "KEGG")
  dir.create(go_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(kegg_dir, showWarnings = FALSE, recursive = TRUE)
  
  deg_files <- list.files(ct_dir, pattern = "_deg\\.csv$", full.names = TRUE)
  cat(sprintf("  found %d DEG files\n", length(deg_files)))
  
  for (fp in deg_files) {
    cell_type <- sub("_deg\\.csv$", "", basename(fp))
    if (cell_type %in% EXCLUDE_CT) next
    
    deg <- tryCatch(read.csv(fp, stringsAsFactors = FALSE),
                    error = function(e) NULL)
    if (is.null(deg) || nrow(deg) == 0) next
    
    # 取显著 + up-regulated gene
    sig_up <- deg %>%
      filter(p_val_adj < PADJ_CUT, avg_log2FC > LOGFC_CUT) %>%
      pull(gene)
    sig_down <- deg %>%
      filter(p_val_adj < PADJ_CUT, avg_log2FC < -LOGFC_CUT) %>%
      pull(gene)
    
    cat(sprintf("  %s: %d up, %d down genes\n",
                cell_type, length(sig_up), length(sig_down)))
    
    if (length(sig_up) < 10) next  # 太少没富集意义
    
    # Map symbol -> ENTREZ
    ez_up <- bitr(sig_up, fromType = "SYMBOL", toType = "ENTREZID",
                  OrgDb = org.Mm.eg.db) %>% suppressMessages()
    if (nrow(ez_up) < 5) next
    
    # === GO enrichment (BP) ===
    go_res <- tryCatch({
      enrichGO(gene = ez_up$ENTREZID,
               OrgDb = org.Mm.eg.db,
               keyType = "ENTREZID",
               ont = "BP",
               pAdjustMethod = "BH",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.2,
               readable = TRUE)
    }, error = function(e) NULL)
    
    if (!is.null(go_res) && nrow(as.data.frame(go_res)) > 0) {
      # Simplify 去冗余
      go_res <- tryCatch(simplify(go_res, cutoff = 0.7),
                         error = function(e) go_res)
      
      go_df <- as.data.frame(go_res)
      write.csv(go_df,
                file.path(go_dir, sprintf("%s_GO_BP.csv", cell_type)),
                row.names = FALSE)
      
      if (nrow(go_df) > 0) {
        n_top <- min(TOP_N_TERM, nrow(go_df))
        
        # 柱形图 ⭐ 客户要求
        p_bar <- tryCatch({
          barplot(go_res, showCategory = n_top, font.size = 9) +
            ggtitle(sprintf("%s - %s\nGO BP (UP)", contrast_name, cell_type)) +
            theme(plot.title = element_text(size = 10))
        }, error = function(e) NULL)
        if (!is.null(p_bar)) {
          ggsave(file.path(go_dir, sprintf("%s_GO_BP_barplot.pdf", cell_type)),
                 p_bar, width = 9, height = 7)
        }
        
        # 点图
        p_dot <- tryCatch({
          dotplot(go_res, showCategory = n_top, font.size = 9) +
            ggtitle(sprintf("%s - %s\nGO BP (UP)", contrast_name, cell_type)) +
            theme(plot.title = element_text(size = 10))
        }, error = function(e) NULL)
        if (!is.null(p_dot)) {
          ggsave(file.path(go_dir, sprintf("%s_GO_BP_dotplot.pdf", cell_type)),
                 p_dot, width = 9, height = 7)
        }
        
        summary_rows[[length(summary_rows) + 1]] <- data.frame(
          contrast = contrast_name,
          cell_type = cell_type,
          source = "GO_BP",
          n_terms = nrow(go_df),
          top_term = go_df$Description[1]
        )
      }
    }
    
    # === KEGG enrichment ===
    kegg_res <- tryCatch({
      enrichKEGG(gene = ez_up$ENTREZID,
                 organism = "mmu",
                 keyType = "kegg",
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.2)
    }, error = function(e) NULL)
    
    if (!is.null(kegg_res) && nrow(as.data.frame(kegg_res)) > 0) {
      # Readable (symbol 替换 entrez)
      kegg_res <- tryCatch(setReadable(kegg_res, OrgDb = org.Mm.eg.db,
                                        keyType = "ENTREZID"),
                           error = function(e) kegg_res)
      
      kegg_df <- as.data.frame(kegg_res)
      write.csv(kegg_df,
                file.path(kegg_dir, sprintf("%s_KEGG.csv", cell_type)),
                row.names = FALSE)
      
      if (nrow(kegg_df) > 0) {
        n_top <- min(TOP_N_TERM, nrow(kegg_df))
        
        # 柱形图 ⭐ 客户要求
        p_bar <- tryCatch({
          barplot(kegg_res, showCategory = n_top, font.size = 9) +
            ggtitle(sprintf("%s - %s\nKEGG pathways (UP)",
                            contrast_name, cell_type)) +
            theme(plot.title = element_text(size = 10))
        }, error = function(e) NULL)
        if (!is.null(p_bar)) {
          ggsave(file.path(kegg_dir, sprintf("%s_KEGG_barplot.pdf", cell_type)),
                 p_bar, width = 9, height = 7)
        }
        
        # 点图
        p_dot <- tryCatch({
          dotplot(kegg_res, showCategory = n_top, font.size = 9) +
            ggtitle(sprintf("%s - %s\nKEGG pathways (UP)",
                            contrast_name, cell_type)) +
            theme(plot.title = element_text(size = 10))
        }, error = function(e) NULL)
        if (!is.null(p_dot)) {
          ggsave(file.path(kegg_dir, sprintf("%s_KEGG_dotplot.pdf", cell_type)),
                 p_dot, width = 9, height = 7)
        }
        
        summary_rows[[length(summary_rows) + 1]] <- data.frame(
          contrast = contrast_name,
          cell_type = cell_type,
          source = "KEGG",
          n_terms = nrow(kegg_df),
          top_term = kegg_df$Description[1]
        )
      }
    }
  }
}

# 合并 summary
if (length(summary_rows) > 0) {
  summary_df <- do.call(rbind, summary_rows)
  summary_csv <- file.path(OUT_DIR, "enrichment_summary.csv")
  write.csv(summary_df, summary_csv, row.names = FALSE)
  
  cat(sprintf("\n\n=== Summary ===\n"))
  cat(sprintf("Total enrichment blocks: %d\n", nrow(summary_df)))
  cat(sprintf("Summary saved: %s\n", summary_csv))
  print(table(summary_df$contrast, summary_df$source))
}

cat(sprintf("\n\nOutput dir: %s\n", OUT_DIR))
cat("Per contrast: GO/ + KEGG/ 各含 CSV + 柱形图 + 点图\n")
