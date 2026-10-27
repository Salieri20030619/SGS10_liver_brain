#!/usr/bin/env Rscript
# ==============================================================================
# SGS10 - D3 富集分析补充: GO + KEGG + ReactomePA
# 范围: 5 个 contrast × {UP, DOWN} = 10 个 task
#   A1: liver tumor_core vs baseline
#   A2: liver invaded_host vs baseline  ⭐
#   A3: liver distant_host vs baseline
#   B3: liver distant_host late vs early (sample-driven, 但客户要)
#   C4: brain distant_host vs brain baseline  ⭐⭐⭐
# 输入: _D3_FINAL/signature_output/filtered_DEG/{contrast}_{UP,DOWN}_filtered.csv
# 输出: D4_GO_output_v2/
#   - enrichment/  per-contrast 三件套表 (GO_BP/GO_MF/GO_CC/KEGG/Reactome)
#   - figures/     per-contrast dotplot + emapplot + cnetplot
#   - cross/       跨 contrast compareCluster 比较 (UP only)
#   - logs/        运行日志
#   - summary.csv  各 task 的成功/失败状态 + 显著 term 数
#
# 运行方式 (服务器):
#   cd /mnt/BioDisk/wangmingyuan/SGS10_analysis
#   nohup Rscript scripts/D3_enrichment_GO_KEGG_Reactome.R > logs/D3_enrich.log 2>&1 &
# 预计耗时: ~30-60 min (clusterProfiler/Reactome 各 9 task)
# 预计内存: < 8 GB
# ==============================================================================

# ---------- 0. 配置 ----------
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(ReactomePA)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
  library(enrichplot)
  library(DOSE)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(patchwork)
})

CFG <- list(
  PROJECT_ROOT     = "/mnt/BioDisk/wangmingyuan/SGS10_analysis",
  IN_DIR           = "_D3_FINAL/signature_output/filtered_DEG",
  COMBINED_DIR     = "_D3_FINAL/DEG_output/combined",
  OUT_DIR          = "D4_GO_output_v2",
  CONTRASTS        = c("A1", "A2", "A3", "B3", "C4"),
  DIRECTIONS       = c("UP", "DOWN"),
  ORG_DB           = "org.Mm.eg.db",
  ORG_KEGG         = "mmu",          # KEGG 物种代码 mouse
  REACTOME_ORG     = "mouse",
  MIN_GENES        = 10,             # 富集最少基因数, 不够则跳过
  PVAL_CUT         = 0.05,
  QVAL_CUT         = 0.10,
  PADJUST          = "BH",
  TOP_N_DOTPLOT    = 20,
  GO_ONTOLOGIES    = c("BP", "MF", "CC"),
  # B3 是 sample-driven (gan56 单样本 pseudogene), 标注但仍跑
  KNOWN_LIMITATIONS = list(
    A1_UP    = "Only 2 robust UP genes (Celf2, Gm42418); enrichment skipped.",
    B3       = "Single-sample (gan56-only) signature, dominated by Mup/Gm pseudogenes."
  )
)

# null-coalesce 兜底 (用于 summary 行的 0 fallback; 必须在主循环前定义)
`%||%` <- function(a, b) if (is.null(a)) b else a

setwd(CFG$PROJECT_ROOT)

OUT_BASE <- file.path(CFG$PROJECT_ROOT, CFG$OUT_DIR)
DIR_ENRICH <- file.path(OUT_BASE, "enrichment")
DIR_FIG    <- file.path(OUT_BASE, "figures")
DIR_CROSS  <- file.path(OUT_BASE, "cross")
DIR_LOG    <- file.path(OUT_BASE, "logs")
for (d in c(DIR_ENRICH, DIR_FIG, DIR_CROSS, DIR_LOG)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ---------- 1. logging ----------
LOG_FILE <- file.path(DIR_LOG, sprintf("D3_enrich_%s.log",
                                       format(Sys.time(), "%Y%m%d_%H%M%S")))
log_msg <- function(...) {
  msg <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = LOG_FILE, append = TRUE)
}

log_msg("=== D3 enrichment 启动 ===")
log_msg("R version: ", R.version.string)
log_msg("clusterProfiler: v", as.character(packageVersion("clusterProfiler")))
log_msg("ReactomePA:      v", as.character(packageVersion("ReactomePA")))
log_msg("org.Mm.eg.db:    v", as.character(packageVersion("org.Mm.eg.db")))
log_msg("OUT_DIR:         ", OUT_BASE)

# ---------- 2. 工具函数 ----------

# Symbol → Entrez (mouse), 返回 list(entrez, n_input, n_mapped, missing)
symbols_to_entrez <- function(symbols) {
  symbols <- unique(symbols[!is.na(symbols) & symbols != ""])
  if (length(symbols) == 0) {
    return(list(entrez = character(0), n_input = 0, n_mapped = 0,
                missing = character(0)))
  }
  map <- tryCatch(
    AnnotationDbi::select(
      org.Mm.eg.db,
      keys = symbols,
      keytype = "SYMBOL",
      columns = c("SYMBOL", "ENTREZID")
    ),
    error = function(e) NULL
  )
  if (is.null(map)) {
    return(list(entrez = character(0), n_input = length(symbols),
                n_mapped = 0, missing = symbols))
  }
  map <- map[!is.na(map$ENTREZID), ]
  # 处理一对多: 每 SYMBOL 取第一个 ENTREZID
  map <- map[!duplicated(map$SYMBOL), ]
  mapped_symbols <- map$SYMBOL
  list(
    entrez   = map$ENTREZID,
    n_input  = length(symbols),
    n_mapped = nrow(map),
    missing  = setdiff(symbols, mapped_symbols)
  )
}

# 读取 filtered DEG 文件 → 返回 gene symbol 向量
read_filtered <- function(contrast, direction) {
  fp <- file.path(CFG$PROJECT_ROOT, CFG$IN_DIR,
                  sprintf("%s_%s_filtered.csv", contrast, direction))
  if (!file.exists(fp)) {
    log_msg("  [WARN] 文件不存在: ", fp); return(NULL)
  }
  df <- tryCatch(read_csv(fp, show_col_types = FALSE),
                 error = function(e) {
                   log_msg("  [ERR] 读 csv 失败: ", conditionMessage(e))
                   NULL
                 })
  if (is.null(df) || nrow(df) == 0) return(NULL)
  if (!"gene" %in% colnames(df)) {
    log_msg("  [ERR] 没有 gene 列, contrast=", contrast, " ", direction)
    return(NULL)
  }
  df
}

# 读 combined 文件作为 expressed background (per contrast)
get_background <- function(contrast) {
  # B3 没 combined 文件 (D3 报告: pseudobulk 跳过), 用 wilcoxon 文件代替
  fp_combined <- file.path(CFG$PROJECT_ROOT, CFG$COMBINED_DIR,
                           sprintf("%s_%s_combined.csv", contrast,
                                   c(A1 = "tumor_core_vs_baseline_host",
                                     A2 = "invaded_host_vs_baseline_host",
                                     A3 = "distant_host_vs_baseline_host",
                                     C4 = "brain_distant_vs_brain_baseline")[contrast]))
  fp_wilcox <- file.path(CFG$PROJECT_ROOT, "_D3_FINAL/DEG_output/wilcoxon",
                         sprintf("%s_%s.csv", contrast,
                                 c(A1 = "tumor_core_vs_baseline_host",
                                   A2 = "invaded_host_vs_baseline_host",
                                   A3 = "distant_host_vs_baseline_host",
                                   B3 = "distant_late_vs_distant_early",
                                   C4 = "brain_distant_vs_brain_baseline")[contrast]))
  bg_file <- if (file.exists(fp_combined)) fp_combined else fp_wilcox
  if (!file.exists(bg_file)) {
    log_msg("  [WARN] 背景文件不存在 (", contrast, "), 用 NULL = 全转录组")
    return(NULL)
  }
  df <- read_csv(bg_file, show_col_types = FALSE)
  if (!"gene" %in% colnames(df)) {
    log_msg("  [WARN] 背景文件没 gene 列 (", contrast, ")")
    return(NULL)
  }
  bg <- symbols_to_entrez(df$gene)
  log_msg("  背景 (contrast=", contrast, "): ", bg$n_mapped,
          "/", bg$n_input, " gene mapped to entrez (from ", basename(bg_file), ")")
  bg$entrez
}

# 单 task 富集主流程 (返回 list 含 GO_BP / GO_MF / GO_CC / KEGG / Reactome)
run_enrichment <- function(gene_entrez, bg_entrez, label) {
  res <- list()
  for (ont in CFG$GO_ONTOLOGIES) {
    log_msg("    GO_", ont, " ...")
    go <- tryCatch(
      enrichGO(gene = gene_entrez,
               universe = bg_entrez,
               OrgDb    = org.Mm.eg.db,
               keyType  = "ENTREZID",
               ont      = ont,
               pvalueCutoff  = CFG$PVAL_CUT,
               qvalueCutoff  = CFG$QVAL_CUT,
               pAdjustMethod = CFG$PADJUST,
               readable = TRUE),
      error = function(e) {
        log_msg("    [ERR] GO_", ont, ": ", conditionMessage(e))
        NULL
      })
    res[[paste0("GO_", ont)]] <- go
  }
  log_msg("    KEGG ...")
  res[["KEGG"]] <- tryCatch(
    enrichKEGG(gene = gene_entrez,
               universe = bg_entrez,
               organism = CFG$ORG_KEGG,
               keyType  = "kegg",
               pvalueCutoff  = CFG$PVAL_CUT,
               qvalueCutoff  = CFG$QVAL_CUT,
               pAdjustMethod = CFG$PADJUST),
    error = function(e) {
      log_msg("    [ERR] KEGG: ", conditionMessage(e)); NULL
    })
  if (!is.null(res[["KEGG"]]) && nrow(as.data.frame(res[["KEGG"]])) > 0) {
    res[["KEGG"]] <- tryCatch(setReadable(res[["KEGG"]], OrgDb = org.Mm.eg.db,
                                          keyType = "ENTREZID"),
                              error = function(e) res[["KEGG"]])
  }
  log_msg("    Reactome ...")
  res[["Reactome"]] <- tryCatch(
    enrichPathway(gene = gene_entrez,
                  universe = bg_entrez,
                  organism = CFG$REACTOME_ORG,
                  pvalueCutoff  = CFG$PVAL_CUT,
                  qvalueCutoff  = CFG$QVAL_CUT,
                  pAdjustMethod = CFG$PADJUST,
                  readable = TRUE),
    error = function(e) {
      log_msg("    [ERR] Reactome: ", conditionMessage(e)); NULL
    })
  res
}

# 保存 enrichResult 到 csv (添加 task 元信息列)
save_result <- function(er, label, type) {
  if (is.null(er)) return(0)
  df <- tryCatch(as.data.frame(er), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) {
    log_msg("    [", type, "] 无显著 term")
    return(0)
  }
  df$task <- label
  df$enrich_type <- type
  fp <- file.path(DIR_ENRICH, sprintf("%s_%s.csv", label, type))
  write_csv(df, fp)
  log_msg("    [", type, "] ", nrow(df), " term, 写入 ", basename(fp))
  nrow(df)
}

# 出图: dotplot (top N) + emapplot + cnetplot
make_plots <- function(er_list, label, gene_logfc) {
  for (type in names(er_list)) {
    er <- er_list[[type]]
    if (is.null(er)) next
    df <- as.data.frame(er)
    if (nrow(df) == 0) next
    n_show <- min(CFG$TOP_N_DOTPLOT, nrow(df))
    
    # dotplot
    p_dot <- tryCatch(
      dotplot(er, showCategory = n_show) +
        ggtitle(sprintf("%s | %s | top %d", label, type, n_show)) +
        theme(plot.title = element_text(size = 11)),
      error = function(e) NULL
    )
    if (!is.null(p_dot)) {
      ggsave(file.path(DIR_FIG, sprintf("%s_%s_dotplot.pdf", label, type)),
             p_dot, width = 9, height = max(6, 0.3 * n_show + 2),
             limitsize = FALSE)
    }
    
    # emapplot (term 间相似性网络)
    if (nrow(df) >= 5) {
      p_emap <- tryCatch({
        er2 <- pairwise_termsim(er)
        emapplot(er2, showCategory = min(30, nrow(df))) +
          ggtitle(sprintf("%s | %s | term similarity", label, type))
      }, error = function(e) {
        log_msg("    [WARN] emapplot 失败 (", label, " ", type, "): ",
                conditionMessage(e)); NULL
      })
      if (!is.null(p_emap)) {
        ggsave(file.path(DIR_FIG, sprintf("%s_%s_emapplot.pdf", label, type)),
               p_emap, width = 11, height = 9)
      }
    }
    
    # cnetplot (gene-term 网络, 仅 top 5 term 避免太密)
    if (!is.null(gene_logfc) && length(gene_logfc) > 0 && nrow(df) >= 1) {
      p_cnet <- tryCatch(
        cnetplot(er, showCategory = min(5, nrow(df)),
                 foldChange = gene_logfc, circular = FALSE,
                 colorEdge = TRUE) +
          ggtitle(sprintf("%s | %s | gene-term network (top 5)", label, type)),
        error = function(e) {
          log_msg("    [WARN] cnetplot 失败 (", label, " ", type, "): ",
                  conditionMessage(e)); NULL
        }
      )
      if (!is.null(p_cnet)) {
        ggsave(file.path(DIR_FIG, sprintf("%s_%s_cnetplot.pdf", label, type)),
               p_cnet, width = 11, height = 9)
      }
    }
  }
}

# ---------- 3. 主循环: per (contrast, direction) task ----------
summary_rows <- list()
gene_lists_for_compare <- list()  # 给 compareCluster 用 (UP only)

for (contrast in CFG$CONTRASTS) {
  for (direction in CFG$DIRECTIONS) {
    label <- sprintf("%s_%s", contrast, direction)
    log_msg("=== Task: ", label, " ===")
    
    df <- read_filtered(contrast, direction)
    if (is.null(df)) {
      summary_rows[[label]] <- data.frame(
        task = label, status = "no_input_file", n_genes = 0,
        n_mapped = 0, GO_BP = 0, GO_MF = 0, GO_CC = 0,
        KEGG = 0, Reactome = 0)
      next
    }
    
    n_genes <- nrow(df)
    log_msg("  原始 gene 数: ", n_genes)
    
    if (n_genes < CFG$MIN_GENES) {
      log_msg("  [SKIP] gene 数 < ", CFG$MIN_GENES, ", 跳过富集")
      summary_rows[[label]] <- data.frame(
        task = label, status = "skipped_too_few_genes",
        n_genes = n_genes, n_mapped = 0,
        GO_BP = 0, GO_MF = 0, GO_CC = 0, KEGG = 0, Reactome = 0)
      next
    }
    
    # Symbol → Entrez
    mapped <- symbols_to_entrez(df$gene)
    log_msg("  Entrez 映射: ", mapped$n_mapped, "/", mapped$n_input,
            " (", round(100*mapped$n_mapped/mapped$n_input, 1), "%)")
    if (length(mapped$missing) > 0) {
      miss_show <- head(mapped$missing, 10)
      log_msg("  未映射示例 (", length(mapped$missing), "): ",
              paste(miss_show, collapse = ", "),
              if (length(mapped$missing) > 10) " ..." else "")
    }
    
    if (mapped$n_mapped < CFG$MIN_GENES) {
      log_msg("  [SKIP] 映射后 gene 数 < ", CFG$MIN_GENES)
      summary_rows[[label]] <- data.frame(
        task = label, status = "skipped_too_few_mapped",
        n_genes = n_genes, n_mapped = mapped$n_mapped,
        GO_BP = 0, GO_MF = 0, GO_CC = 0, KEGG = 0, Reactome = 0)
      next
    }
    
    # 背景: per contrast 的 expressed gene set
    bg_entrez <- get_background(contrast)
    
    # 准备 cnetplot 用的 fold-change vector (entrez-named)
    gene_lookup <- AnnotationDbi::select(
      org.Mm.eg.db, keys = df$gene, keytype = "SYMBOL",
      columns = c("SYMBOL", "ENTREZID")
    )
    gene_lookup <- gene_lookup[!is.na(gene_lookup$ENTREZID), ]
    gene_lookup <- gene_lookup[!duplicated(gene_lookup$SYMBOL), ]
    df_with_ez <- merge(df, gene_lookup, by.x = "gene", by.y = "SYMBOL")
    fc_vec <- df_with_ez$avg_log2FC
    names(fc_vec) <- df_with_ez$ENTREZID
    fc_vec <- sort(fc_vec, decreasing = TRUE)
    
    # 富集
    er_list <- run_enrichment(mapped$entrez, bg_entrez, label)
    
    # 保存 csv
    n_terms <- list()
    for (type in names(er_list)) {
      n_terms[[type]] <- save_result(er_list[[type]], label, type)
    }
    
    # 出图
    make_plots(er_list, label, fc_vec)
    
    # 给 compareCluster 准备 (UP, 且 > MIN_GENES)
    if (direction == "UP") {
      gene_lists_for_compare[[contrast]] <- mapped$entrez
    }
    
    summary_rows[[label]] <- data.frame(
      task = label, status = "ok",
      n_genes = n_genes, n_mapped = mapped$n_mapped,
      GO_BP = n_terms[["GO_BP"]] %||% 0,
      GO_MF = n_terms[["GO_MF"]] %||% 0,
      GO_CC = n_terms[["GO_CC"]] %||% 0,
      KEGG  = n_terms[["KEGG"]] %||% 0,
      Reactome = n_terms[["Reactome"]] %||% 0)
  }
}

# ---------- 4. 跨 contrast compareCluster ----------
summary_df <- do.call(rbind, summary_rows)
write_csv(summary_df, file.path(OUT_BASE, "task_summary.csv"))
log_msg("\n=== Task summary ===")
print(summary_df)

log_msg("\n=== compareCluster (UP genes, A2/A3/B3/C4) ===")

# 只比较有足够 gene 的 contrast
cc_input <- gene_lists_for_compare[
  sapply(gene_lists_for_compare, length) >= CFG$MIN_GENES
]
log_msg("compareCluster input contrasts: ",
        paste(names(cc_input), collapse = ", "),
        " (sizes: ", paste(sapply(cc_input, length), collapse = ", "), ")")

if (length(cc_input) >= 2) {
  # GO BP cross-contrast
  cc_go <- tryCatch(
    compareCluster(geneClusters = cc_input,
                   fun = "enrichGO",
                   OrgDb = org.Mm.eg.db,
                   keyType = "ENTREZID",
                   ont = "BP",
                   pvalueCutoff  = CFG$PVAL_CUT,
                   qvalueCutoff  = CFG$QVAL_CUT,
                   pAdjustMethod = CFG$PADJUST,
                   readable = TRUE),
    error = function(e) {
      log_msg("[ERR] compareCluster GO_BP: ", conditionMessage(e)); NULL
    })
  if (!is.null(cc_go)) {
    write_csv(as.data.frame(cc_go),
              file.path(DIR_CROSS, "compareCluster_GO_BP.csv"))
    p_cc_go <- dotplot(cc_go, showCategory = 10) +
      ggtitle("Cross-contrast GO_BP (UP)") +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))
    ggsave(file.path(DIR_CROSS, "compareCluster_GO_BP_dotplot.pdf"),
           p_cc_go, width = 11, height = 10)
    log_msg("compareCluster GO_BP: ",
            length(unique(as.data.frame(cc_go)$ID)), " unique terms")
  }
  
  # KEGG cross-contrast
  cc_kegg <- tryCatch(
    compareCluster(geneClusters = cc_input,
                   fun = "enrichKEGG",
                   organism = CFG$ORG_KEGG,
                   pvalueCutoff  = CFG$PVAL_CUT,
                   qvalueCutoff  = CFG$QVAL_CUT,
                   pAdjustMethod = CFG$PADJUST),
    error = function(e) {
      log_msg("[ERR] compareCluster KEGG: ", conditionMessage(e)); NULL
    })
  if (!is.null(cc_kegg)) {
    cc_kegg <- tryCatch(setReadable(cc_kegg, OrgDb = org.Mm.eg.db,
                                    keyType = "ENTREZID"),
                        error = function(e) cc_kegg)
    write_csv(as.data.frame(cc_kegg),
              file.path(DIR_CROSS, "compareCluster_KEGG.csv"))
    if (nrow(as.data.frame(cc_kegg)) > 0) {
      p_cc_kegg <- dotplot(cc_kegg, showCategory = 10) +
        ggtitle("Cross-contrast KEGG (UP)") +
        theme(axis.text.x = element_text(angle = 30, hjust = 1))
      ggsave(file.path(DIR_CROSS, "compareCluster_KEGG_dotplot.pdf"),
             p_cc_kegg, width = 11, height = 10)
      log_msg("compareCluster KEGG: ",
              length(unique(as.data.frame(cc_kegg)$ID)), " unique terms")
    }
  }
  
  # Reactome cross-contrast
  cc_react <- tryCatch(
    compareCluster(geneClusters = cc_input,
                   fun = "enrichPathway",
                   organism = CFG$REACTOME_ORG,
                   pvalueCutoff  = CFG$PVAL_CUT,
                   qvalueCutoff  = CFG$QVAL_CUT,
                   pAdjustMethod = CFG$PADJUST,
                   readable = TRUE),
    error = function(e) {
      log_msg("[ERR] compareCluster Reactome: ", conditionMessage(e)); NULL
    })
  if (!is.null(cc_react) && nrow(as.data.frame(cc_react)) > 0) {
    write_csv(as.data.frame(cc_react),
              file.path(DIR_CROSS, "compareCluster_Reactome.csv"))
    p_cc_re <- dotplot(cc_react, showCategory = 10) +
      ggtitle("Cross-contrast Reactome (UP)") +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))
    ggsave(file.path(DIR_CROSS, "compareCluster_Reactome_dotplot.pdf"),
           p_cc_re, width = 11, height = 10)
    log_msg("compareCluster Reactome: ",
            length(unique(as.data.frame(cc_react)$ID)), " unique terms")
  }
}

# ---------- 5. 收尾 ----------
log_msg("\n=== 全部完成 ===")
log_msg("输出目录: ", OUT_BASE)
log_msg("  enrichment/: per-contrast 富集表 (csv)")
log_msg("  figures/:    dotplot + emapplot + cnetplot per task")
log_msg("  cross/:      compareCluster 跨 contrast 比较")
log_msg("  task_summary.csv: 各 task 状态")
log_msg("  logs/D3_enrich_*.log: 本次运行日志")

sessionInfo()
