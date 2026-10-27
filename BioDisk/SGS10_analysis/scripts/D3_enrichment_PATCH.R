#!/usr/bin/env Rscript
# ==============================================================================
# SGS10 - D3 富集脚本 PATCH
# 修复点:
#   (1) cnetplot 全失败 -> enrichplot 1.30.5 API 变更, 移除 circular/colorEdge 参数
#   (2) A1_DOWN_KEGG 因网络瞬时断开而失败 -> 重跑这一个
# 不重做: 已经成功的 csv / dotplot / emapplot (~25 min 工作不浪费)
# 仅生成: 缺失的 cnetplot pdf + A1_DOWN_KEGG.csv (+ 该 KEGG 的 plots)
# ==============================================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(ReactomePA)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
  library(enrichplot)
  library(DOSE)
  library(ggplot2)
  library(dplyr)
  library(readr)
})

CFG <- list(
  PROJECT_ROOT = "/mnt/BioDisk/wangmingyuan/SGS10_analysis",
  IN_DIR       = "_D3_FINAL/signature_output/filtered_DEG",
  COMBINED_DIR = "_D3_FINAL/DEG_output/combined",
  WILCOX_DIR   = "_D3_FINAL/DEG_output/wilcoxon",
  OUT_DIR      = "D4_GO_output_v2",
  ORG_KEGG     = "mmu",
  REACTOME_ORG = "mouse",
  PVAL_CUT     = 0.05,
  QVAL_CUT     = 0.10,
  PADJUST      = "BH",
  TOP_N_DOTPLOT = 20
)
setwd(CFG$PROJECT_ROOT)

OUT_BASE   <- file.path(CFG$PROJECT_ROOT, CFG$OUT_DIR)
DIR_ENRICH <- file.path(OUT_BASE, "enrichment")
DIR_FIG    <- file.path(OUT_BASE, "figures")
DIR_LOG    <- file.path(OUT_BASE, "logs")
LOG_FILE   <- file.path(DIR_LOG, sprintf("D3_enrich_PATCH_%s.log",
                                         format(Sys.time(), "%Y%m%d_%H%M%S")))
log_msg <- function(...) {
  msg <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"),
                 paste0(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = LOG_FILE, append = TRUE)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

log_msg("=== PATCH 启动 ===")
log_msg("enrichplot version: ", as.character(packageVersion("enrichplot")))

# ---------- 工具: Symbol -> Entrez ----------
symbols_to_entrez <- function(symbols) {
  symbols <- unique(symbols[!is.na(symbols) & symbols != ""])
  if (length(symbols) == 0) return(list(entrez = character(0), missing = character(0)))
  map <- AnnotationDbi::select(org.Mm.eg.db, keys = symbols,
                               keytype = "SYMBOL",
                               columns = c("SYMBOL", "ENTREZID"))
  map <- map[!is.na(map$ENTREZID), ]
  map <- map[!duplicated(map$SYMBOL), ]
  list(entrez = map$ENTREZID, mapped_df = map,
       missing = setdiff(symbols, map$SYMBOL))
}

# ---------- 工具: 读 filtered DEG ----------
read_filtered <- function(contrast, direction) {
  fp <- file.path(CFG$PROJECT_ROOT, CFG$IN_DIR,
                  sprintf("%s_%s_filtered.csv", contrast, direction))
  if (!file.exists(fp)) return(NULL)
  read_csv(fp, show_col_types = FALSE)
}

# ---------- 工具: 准备 cnetplot 用的 fold-change vector ----------
make_fc_vec <- function(df) {
  if (is.null(df) || !"gene" %in% colnames(df)) return(NULL)
  lk <- AnnotationDbi::select(org.Mm.eg.db, keys = df$gene,
                              keytype = "SYMBOL",
                              columns = c("SYMBOL", "ENTREZID"))
  lk <- lk[!is.na(lk$ENTREZID), ]
  lk <- lk[!duplicated(lk$SYMBOL), ]
  m <- merge(df, lk, by.x = "gene", by.y = "SYMBOL")
  v <- m$avg_log2FC; names(v) <- m$ENTREZID
  sort(v, decreasing = TRUE)
}

# ============================================================================
# 部分 A: cnetplot 修复 - 用新 API 重画
# 策略: 不重跑富集 (csv 都还在), 而是从 csv 反向构造 enrichResult-like 对象
#       → 但 cnetplot 需要 enrichResult 对象, 不能直接从 csv 来
#       → 替代方案: 重跑富集本身 (内存里, 不写 csv), 然后只画 cnetplot
#       → 富集本身只读 ~10s/task, 不慢
# ============================================================================

# 跑富集 (in-memory, 不写盘) 用于支撑 cnetplot
run_one_enrich <- function(gene_entrez, bg_entrez, type) {
  if (length(gene_entrez) < 3) return(NULL)
  if (type == "GO_BP") {
    enrichGO(gene = gene_entrez, universe = bg_entrez,
             OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "BP",
             pvalueCutoff = CFG$PVAL_CUT, qvalueCutoff = CFG$QVAL_CUT,
             pAdjustMethod = CFG$PADJUST, readable = TRUE)
  } else if (type == "GO_MF") {
    enrichGO(gene = gene_entrez, universe = bg_entrez,
             OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "MF",
             pvalueCutoff = CFG$PVAL_CUT, qvalueCutoff = CFG$QVAL_CUT,
             pAdjustMethod = CFG$PADJUST, readable = TRUE)
  } else if (type == "GO_CC") {
    enrichGO(gene = gene_entrez, universe = bg_entrez,
             OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "CC",
             pvalueCutoff = CFG$PVAL_CUT, qvalueCutoff = CFG$QVAL_CUT,
             pAdjustMethod = CFG$PADJUST, readable = TRUE)
  } else if (type == "KEGG") {
    er <- enrichKEGG(gene = gene_entrez, universe = bg_entrez,
                     organism = CFG$ORG_KEGG, keyType = "kegg",
                     pvalueCutoff = CFG$PVAL_CUT, qvalueCutoff = CFG$QVAL_CUT,
                     pAdjustMethod = CFG$PADJUST)
    if (!is.null(er) && nrow(as.data.frame(er)) > 0) {
      er <- tryCatch(setReadable(er, OrgDb = org.Mm.eg.db, keyType = "ENTREZID"),
                     error = function(e) er)
    }
    er
  } else if (type == "Reactome") {
    enrichPathway(gene = gene_entrez, universe = bg_entrez,
                  organism = CFG$REACTOME_ORG,
                  pvalueCutoff = CFG$PVAL_CUT, qvalueCutoff = CFG$QVAL_CUT,
                  pAdjustMethod = CFG$PADJUST, readable = TRUE)
  }
}

# 新 API cnetplot (移除 circular / colorEdge)
make_cnetplot <- function(er, label, type, fc_vec) {
  if (is.null(er) || nrow(as.data.frame(er)) == 0) return(invisible(NULL))
  n_show <- min(5, nrow(as.data.frame(er)))
  
  # enrichplot 1.30+ API: 用 color.params=list(foldChange=fc_vec)
  # 多版本兼容: 先试新签名, 失败则降级
  p <- tryCatch(
    cnetplot(er, showCategory = n_show,
             color.params = list(foldChange = fc_vec)) +
      ggtitle(sprintf("%s | %s | gene-term network (top %d)",
                      label, type, n_show)),
    error = function(e1) {
      log_msg("    [INFO] cnetplot 新 API 失败 (", e1$message,
              "), 降级到无 foldChange 版本")
      tryCatch(
        cnetplot(er, showCategory = n_show) +
          ggtitle(sprintf("%s | %s | gene-term network (top %d, no FC)",
                          label, type, n_show)),
        error = function(e2) {
          log_msg("    [ERR] cnetplot 完全失败: ", e2$message); NULL
        })
    })
  
  if (!is.null(p)) {
    fp <- file.path(DIR_FIG, sprintf("%s_%s_cnetplot.pdf", label, type))
    ggsave(fp, p, width = 11, height = 9)
    log_msg("    [OK] ", basename(fp))
  }
}

# 已经有 csv 的 task 列表 (扫 enrichment/ 目录)
existing_csvs <- list.files(DIR_ENRICH, pattern = "\\.csv$", full.names = FALSE)
log_msg("\n现有 enrichment csv: ", length(existing_csvs), " 个")

# 解析 csv 文件名 → (label, type)
parse_csv_name <- function(fname) {
  # 例: A2_UP_GO_BP.csv → label="A2_UP", type="GO_BP"
  base <- sub("\\.csv$", "", fname)
  # 找出 type (GO_BP/GO_MF/GO_CC/KEGG/Reactome)
  for (t in c("GO_BP", "GO_MF", "GO_CC", "KEGG", "Reactome")) {
    if (endsWith(base, t)) {
      lab <- sub(paste0("_", t, "$"), "", base)
      return(list(label = lab, type = t))
    }
  }
  NULL
}

# Get expressed background per contrast (与主脚本同逻辑)
get_background <- function(contrast) {
  fp_combined_map <- c(A1 = "A1_tumor_core_vs_baseline_host_combined.csv",
                       A2 = "A2_invaded_host_vs_baseline_host_combined.csv",
                       A3 = "A3_distant_host_vs_baseline_host_combined.csv",
                       C4 = "C4_brain_distant_vs_brain_baseline_combined.csv")
  fp_wilcox_map <- c(A1 = "A1_tumor_core_vs_baseline_host.csv",
                     A2 = "A2_invaded_host_vs_baseline_host.csv",
                     A3 = "A3_distant_host_vs_baseline_host.csv",
                     B3 = "B3_distant_late_vs_distant_early.csv",
                     C4 = "C4_brain_distant_vs_brain_baseline.csv")
  fp_combined <- if (contrast %in% names(fp_combined_map))
    file.path(CFG$PROJECT_ROOT, CFG$COMBINED_DIR, fp_combined_map[contrast]) else ""
  fp_wilcox <- if (contrast %in% names(fp_wilcox_map))
    file.path(CFG$PROJECT_ROOT, CFG$WILCOX_DIR, fp_wilcox_map[contrast]) else ""
  bg_file <- if (file.exists(fp_combined)) fp_combined else fp_wilcox
  if (!file.exists(bg_file)) return(NULL)
  df <- read_csv(bg_file, show_col_types = FALSE)
  if (!"gene" %in% colnames(df)) return(NULL)
  symbols_to_entrez(df$gene)$entrez
}

# 按 (label, type) 分组, 同 label 只重跑一次富集
csv_meta <- do.call(rbind, lapply(existing_csvs, function(f) {
  p <- parse_csv_name(f)
  if (is.null(p)) return(NULL)
  data.frame(file = f, label = p$label, type = p$type, stringsAsFactors = FALSE)
}))
csv_meta <- csv_meta[!is.na(csv_meta$label), ]
log_msg("解析后 task: ", nrow(csv_meta), " 个 type-level (覆盖 ",
        length(unique(csv_meta$label)), " label)")
print(table(csv_meta$label))

# ---------- 修补 cnetplot ----------
log_msg("\n=== 部分 A: 重画 cnetplot (用新 API) ===")
unique_labels <- unique(csv_meta$label)

for (label in unique_labels) {
  parts <- strsplit(label, "_")[[1]]
  contrast <- parts[1]
  direction <- parts[2]
  log_msg("\n[", label, "]")
  
  df <- read_filtered(contrast, direction)
  if (is.null(df)) {
    log_msg("  [SKIP] 没有 filtered DEG"); next
  }
  mapped <- symbols_to_entrez(df$gene)
  if (length(mapped$entrez) < 10) {
    log_msg("  [SKIP] mapped < 10"); next
  }
  bg <- get_background(contrast)
  fc <- make_fc_vec(df)
  
  types_for_label <- csv_meta$type[csv_meta$label == label]
  for (type in types_for_label) {
    cnet_fp <- file.path(DIR_FIG, sprintf("%s_%s_cnetplot.pdf", label, type))
    if (file.exists(cnet_fp) && file.size(cnet_fp) > 1024) {
      log_msg("  [", type, "] cnetplot 已存在, 跳过")
      next
    }
    log_msg("  [", type, "] 重跑富集 + 画 cnetplot ...")
    er <- tryCatch(run_one_enrich(mapped$entrez, bg, type),
                   error = function(e) {
                     log_msg("    [ERR] 富集失败: ", e$message); NULL
                   })
    make_cnetplot(er, label, type, fc)
  }
}

# ---------- 部分 B: 重跑 KEGG (网络断了的几个 task) ----------
# 主 log 显示: A1_DOWN_KEGG 在 16:16:35 网络断, C4_DOWN_KEGG 在 16:36 网络再次抽
# 检查: enrichment/ 下哪些 contrast 的 KEGG.csv 缺失或 size 极小, 都重跑
log_msg("\n=== 部分 B: 重跑可能受网络影响的 KEGG task ===")

# 网络重试包装器: 自动重试最多 3 次, 间隔 30 秒
run_kegg_retry <- function(gene_entrez, bg_entrez, label, max_tries = 3) {
  for (i in seq_len(max_tries)) {
    log_msg("  [", label, "] KEGG 第 ", i, "/", max_tries, " 次尝试 ...")
    er <- tryCatch(
      enrichKEGG(gene = gene_entrez, universe = bg_entrez,
                 organism = CFG$ORG_KEGG, keyType = "kegg",
                 pvalueCutoff = CFG$PVAL_CUT, qvalueCutoff = CFG$QVAL_CUT,
                 pAdjustMethod = CFG$PADJUST),
      error = function(e) {
        log_msg("  [WARN] ", e$message)
        NULL
      })
    if (!is.null(er)) {
      log_msg("  [OK] KEGG 成功, 显著 term = ",
              nrow(as.data.frame(er)))
      if (nrow(as.data.frame(er)) > 0) {
        er <- tryCatch(setReadable(er, OrgDb = org.Mm.eg.db, keyType = "ENTREZID"),
                       error = function(e) er)
      }
      return(er)
    }
    if (i < max_tries) {
      log_msg("  [INFO] 等 30 秒后重试 ...")
      Sys.sleep(30)
    }
  }
  log_msg("  [FAIL] ", label, " KEGG 所有重试均失败")
  NULL
}

# 处理一个 KEGG task: 重跑 + 写 csv + 画 dotplot/emap/cnet
do_kegg_task <- function(contrast, direction) {
  label <- sprintf("%s_%s", contrast, direction)
  csv_fp <- file.path(DIR_ENRICH, sprintf("%s_KEGG.csv", label))
  
  # 已存在且不是网络失败标记 → 跳过
  if (file.exists(csv_fp)) {
    sz <- file.size(csv_fp)
    df_check <- tryCatch(read_csv(csv_fp, show_col_types = FALSE),
                         error = function(e) NULL)
    if (!is.null(df_check) && nrow(df_check) > 0 &&
        !"note" %in% colnames(df_check)) {
      log_msg("[", label, " KEGG] csv 已存在 (n=", nrow(df_check), "), 跳过")
      return(invisible())
    }
  }
  
  log_msg("[", label, " KEGG] 重跑 ...")
  df <- read_filtered(contrast, direction)
  if (is.null(df)) { log_msg("  [SKIP] 没有 filtered DEG"); return(invisible()) }
  mapped <- symbols_to_entrez(df$gene)
  if (length(mapped$entrez) < 10) {
    log_msg("  [SKIP] mapped < 10"); return(invisible())
  }
  bg <- get_background(contrast)
  fc <- make_fc_vec(df)
  
  er <- run_kegg_retry(mapped$entrez, bg, label)
  if (is.null(er)) {
    log_msg("  [FAIL] 写空 csv 标记")
    write_csv(data.frame(task = label, enrich_type = "KEGG",
                         note = "network_failure"), csv_fp)
    return(invisible())
  }
  
  ok_df <- as.data.frame(er)
  if (nrow(ok_df) == 0) {
    log_msg("  [INFO] 无显著 term, 写空标记")
    write_csv(data.frame(task = label, enrich_type = "KEGG",
                         note = "no_significant_terms"), csv_fp)
    return(invisible())
  }
  
  ok_df$task <- label
  ok_df$enrich_type <- "KEGG"
  write_csv(ok_df, csv_fp)
  log_msg("  [OK] ", basename(csv_fp), ": ", nrow(ok_df), " term")
  
  # plots
  n_show <- min(20, nrow(ok_df))
  p_dot <- tryCatch(
    dotplot(er, showCategory = n_show) +
      ggtitle(sprintf("%s | KEGG | top %d", label, n_show)),
    error = function(e) NULL)
  if (!is.null(p_dot)) {
    ggsave(file.path(DIR_FIG, sprintf("%s_KEGG_dotplot.pdf", label)),
           p_dot, width = 9, height = max(6, 0.3 * n_show + 2),
           limitsize = FALSE)
  }
  if (nrow(ok_df) >= 5) {
    p_emap <- tryCatch({
      er2 <- pairwise_termsim(er)
      emapplot(er2, showCategory = min(30, nrow(ok_df))) +
        ggtitle(sprintf("%s | KEGG | term similarity", label))
    }, error = function(e) NULL)
    if (!is.null(p_emap)) {
      ggsave(file.path(DIR_FIG, sprintf("%s_KEGG_emapplot.pdf", label)),
             p_emap, width = 11, height = 9)
    }
  }
  make_cnetplot(er, label, "KEGG", fc)
}

# 检查所有 contrast (除 A1_UP) 的 KEGG csv 是否健康
kegg_checklist <- list(
  c("A1", "DOWN"),  # 第一次网络断
  c("C4", "DOWN")   # log 末尾网络再次抽 - 但可能是真的无显著
)
for (cd in kegg_checklist) {
  do_kegg_task(cd[1], cd[2])
}

# ---------- 收尾: 更新 task_summary.csv 中的 KEGG 列 ----------
log_msg("\n=== 部分 C: 更新 task_summary.csv ===")
ts_fp <- file.path(OUT_BASE, "task_summary.csv")
if (file.exists(ts_fp)) {
  ts <- read_csv(ts_fp, show_col_types = FALSE)
  for (lbl in c("A1_DOWN", "C4_DOWN")) {
    if (lbl %in% ts$task) {
      csv_path <- file.path(DIR_ENRICH, sprintf("%s_KEGG.csv", lbl))
      n_terms <- 0
      if (file.exists(csv_path)) {
        d <- tryCatch(read_csv(csv_path, show_col_types = FALSE),
                      error = function(e) NULL)
        if (!is.null(d) && !"note" %in% colnames(d)) n_terms <- nrow(d)
      }
      ts$KEGG[ts$task == lbl] <- n_terms
      log_msg("task_summary.csv 已更新 ", lbl, ".KEGG = ", n_terms)
    }
  }
  write_csv(ts, ts_fp)
}

log_msg("\n=== PATCH 完成 ===")
log_msg("现在 figures/ 下的 cnetplot 数:")
n_cnet <- length(list.files(DIR_FIG, pattern = "cnetplot\\.pdf$"))
log_msg("  cnetplot 数: ", n_cnet, " (期望 ~9-15)")
log_msg("查看输出: ls ", DIR_FIG, "/*cnetplot*")
