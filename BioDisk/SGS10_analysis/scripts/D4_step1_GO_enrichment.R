#!/usr/bin/env Rscript
# =============================================================================
# SGS10 D4 step1: GO/KEGG/Reactome enrichment
# -----------------------------------------------------------------------------
# Input:  D3_signature_output/filtered_DEG/{A1,A2,A3,C4}_{UP,DOWN}_filtered.csv
#         (B3 skipped - genes are mostly pseudogenes/predicted, no enrichment)
#
# Outputs:
#   D4_GO_output/
#     enrichment/<contrast>_<dir>_GO.csv        (BP+MF+CC merged)
#     enrichment/<contrast>_<dir>_KEGG.csv
#     enrichment/<contrast>_<dir>_Reactome.csv  (optional, if package available)
#     figures/<contrast>_<dir>_dotplot.pdf      (top 10 BP + KEGG + Reactome)
#     figures/combined_top_pathways.pdf         (cross-contrast comparison)
#     pathway_summary.csv
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
})

# ReactomePA is optional
HAS_REACTOME <- tryCatch({
  suppressPackageStartupMessages(library(ReactomePA))
  TRUE
}, error = function(e) FALSE)

option_list <- list(
  make_option("--d3_dir", type="character",
              default="/mnt/BioDisk/wangmingyuan/SGS10_analysis/D3_signature_output"),
  make_option("--out_dir", type="character",
              default="/mnt/BioDisk/wangmingyuan/SGS10_analysis/D4_GO_output"),
  make_option("--pvalue_cutoff", type="double", default=0.05),
  make_option("--qvalue_cutoff", type="double", default=0.10),
  make_option("--top_n_show", type="integer", default=10,
              help="top N pathways to show in dotplot")
)
opt <- parse_args(OptionParser(option_list=option_list))

CFG <- list(
  d3_dir = opt$d3_dir,
  out_dir = opt$out_dir,
  pvalue_cutoff = opt$pvalue_cutoff,
  qvalue_cutoff = opt$qvalue_cutoff,
  top_n_show = opt$top_n_show,
  contrasts = list(
    A2 = list(label = "Liver invaded vs baseline",  organ = "liver", color = "#d7191c"),
    A3 = list(label = "Liver distant vs baseline",  organ = "liver", color = "#fdae61"),
    C4 = list(label = "Brain distant vs baseline",  organ = "brain", color = "#2c7bb6"),
    A1 = list(label = "Liver tumor_core vs baseline", organ = "liver", color = "#762a83")
  ),
  directions = c("UP", "DOWN")
)

dir.create(CFG$out_dir, showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(CFG$out_dir, "enrichment"), showWarnings=FALSE)
dir.create(file.path(CFG$out_dir, "figures"),    showWarnings=FALSE)

cat_log <- function(...) cat(..., sep="")
sec <- function(name) {
  cat_log("\n", strrep("-", 70), "\n  ", name, "\n",
          strrep("-", 70), "\n")
}

cat_log("=====================================================================\n")
cat_log("D4 Step 1: GO/KEGG/Reactome enrichment\n")
cat_log("Start: ", format(Sys.time()), "\n")
cat_log("Reactome available: ", HAS_REACTOME, "\n")
cat_log("=====================================================================\n")

# ============ Helper: convert symbols to ENTREZID ============

symbols_to_entrez <- function(symbols) {
  res <- bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID",
              OrgDb = org.Mm.eg.db, drop = TRUE)
  res$ENTREZID
}

# ============ Helper: run all enrichments for one gene set ============

run_enrichment_panel <- function(gene_symbols, label, out_prefix) {
  cat_log("\n  [", label, "] gene set size = ", length(gene_symbols), "\n")
  
  if (length(gene_symbols) < 5) {
    cat_log("    too few genes, skip\n")
    return(NULL)
  }
  
  # Symbol -> ENTREZID
  entrez <- symbols_to_entrez(gene_symbols)
  cat_log("    SYMBOL->ENTREZID: ", length(gene_symbols), " -> ",
          length(entrez), " mapped\n")
  if (length(entrez) < 5) return(NULL)
  
  results <- list()
  
  # === GO BP / MF / CC ===
  go_combined <- list()
  for (ont in c("BP", "MF", "CC")) {
    cat_log("    enrichGO ", ont, " ...\n")
    res <- tryCatch({
      enrichGO(gene = entrez,
               OrgDb = org.Mm.eg.db,
               keyType = "ENTREZID",
               ont = ont,
               pAdjustMethod = "BH",
               pvalueCutoff = CFG$pvalue_cutoff,
               qvalueCutoff = CFG$qvalue_cutoff,
               readable = TRUE)
    }, error = function(e) {
      cat_log("      GO ", ont, " failed: ", conditionMessage(e), "\n")
      NULL
    })
    if (!is.null(res) && nrow(as.data.frame(res)) > 0) {
      df <- as.data.frame(res)
      df$ontology <- ont
      go_combined[[ont]] <- df
      cat_log("      ", ont, ": ", nrow(df), " enriched\n")
    } else {
      cat_log("      ", ont, ": 0 enriched\n")
    }
  }
  if (length(go_combined) > 0) {
    go_all <- bind_rows(go_combined) %>% arrange(p.adjust)
    write_csv(go_all, paste0(out_prefix, "_GO.csv"))
    results$GO <- go_all
    cat_log("    saved GO: ", paste0(out_prefix, "_GO.csv"),
            " (", nrow(go_all), " rows)\n")
  }
  
  # === KEGG ===
  cat_log("    enrichKEGG ...\n")
  kegg_res <- tryCatch({
    enrichKEGG(gene = entrez,
               organism = "mmu",
               pvalueCutoff = CFG$pvalue_cutoff,
               qvalueCutoff = CFG$qvalue_cutoff)
  }, error = function(e) {
    cat_log("      KEGG failed: ", conditionMessage(e), "\n")
    NULL
  })
  if (!is.null(kegg_res) && nrow(as.data.frame(kegg_res)) > 0) {
    # readable conversion (manual, since readable=TRUE not supported for KEGG)
    kegg_df <- as.data.frame(kegg_res)
    # convert IDs in geneID column
    kegg_df$geneSymbols <- sapply(strsplit(kegg_df$geneID, "/"), function(ids) {
      sym <- bitr(ids, fromType="ENTREZID", toType="SYMBOL",
                  OrgDb=org.Mm.eg.db, drop=TRUE)$SYMBOL
      paste(sym, collapse="/")
    })
    write_csv(kegg_df, paste0(out_prefix, "_KEGG.csv"))
    results$KEGG <- kegg_df
    cat_log("    saved KEGG: ", paste0(out_prefix, "_KEGG.csv"),
            " (", nrow(kegg_df), " rows)\n")
  } else {
    cat_log("    KEGG: 0 enriched\n")
  }
  
  # === Reactome (optional) ===
  if (HAS_REACTOME) {
    cat_log("    enrichPathway (Reactome) ...\n")
    re_res <- tryCatch({
      enrichPathway(gene = entrez,
                    organism = "mouse",
                    pvalueCutoff = CFG$pvalue_cutoff,
                    qvalueCutoff = CFG$qvalue_cutoff,
                    readable = TRUE)
    }, error = function(e) {
      cat_log("      Reactome failed: ", conditionMessage(e), "\n")
      NULL
    })
    if (!is.null(re_res) && nrow(as.data.frame(re_res)) > 0) {
      re_df <- as.data.frame(re_res)
      write_csv(re_df, paste0(out_prefix, "_Reactome.csv"))
      results$Reactome <- re_df
      cat_log("    saved Reactome: ", paste0(out_prefix, "_Reactome.csv"),
              " (", nrow(re_df), " rows)\n")
    } else {
      cat_log("    Reactome: 0 enriched\n")
    }
  }
  
  results
}

# ============ Helper: build dotplot from enrichment df ============

build_dotplot <- function(enrich_df, title, top_n = 10, color_high = "red") {
  if (is.null(enrich_df) || nrow(enrich_df) == 0) {
    return(ggplot() +
             annotate("text", x=0.5, y=0.5,
                      label=paste0(title, "\n(no enrichment)"),
                      size=4) +
             theme_void())
  }
  df <- enrich_df %>%
    arrange(p.adjust) %>%
    head(top_n) %>%
    mutate(
      Description = factor(Description, levels = rev(Description)),
      neg_log10_padj = -log10(p.adjust),
      GeneRatioNum = sapply(strsplit(as.character(GeneRatio), "/"),
                             function(x) as.numeric(x[1]) / as.numeric(x[2]))
    )
  ggplot(df, aes(x = neg_log10_padj, y = Description)) +
    geom_point(aes(size = Count, color = GeneRatioNum)) +
    scale_color_gradient(low = "lightgray", high = color_high,
                         name = "GeneRatio") +
    scale_size(range = c(2, 7), name = "Count") +
    labs(title = title, x = "-log10(p.adjust)", y = NULL) +
    theme_minimal(base_size = 9) +
    theme(plot.title = element_text(size = 10, hjust = 0.5),
          axis.text.y = element_text(size = 8))
}

# ============ Main loop ============

all_results <- list()
summary_rows <- list()

for (ct_id in names(CFG$contrasts)) {
  ci <- CFG$contrasts[[ct_id]]
  for (dir_label in CFG$directions) {
    sec(paste0("Contrast ", ct_id, " ", dir_label, ": ", ci$label))
    
    fp <- file.path(CFG$d3_dir, "filtered_DEG",
                    paste0(ct_id, "_", dir_label, "_filtered.csv"))
    if (!file.exists(fp)) {
      cat_log("  filter file missing: ", fp, "\n")
      next
    }
    
    df <- read_csv(fp, show_col_types = FALSE)
    cat_log("  loaded ", nrow(df), " genes from ", fp, "\n")
    if (nrow(df) < 5) {
      cat_log("  too few genes, skip\n")
      next
    }
    
    out_prefix <- file.path(CFG$out_dir, "enrichment",
                            paste0(ct_id, "_", dir_label))
    res <- run_enrichment_panel(df$gene, paste0(ct_id, "_", dir_label),
                                 out_prefix)
    
    if (is.null(res)) next
    
    all_results[[paste0(ct_id, "_", dir_label)]] <- res
    
    # Build dotplot panel (BP + KEGG + Reactome if available)
    plots <- list()
    
    if (!is.null(res$GO)) {
      bp_df <- res$GO %>% filter(ontology == "BP")
      plots$BP <- build_dotplot(bp_df,
                                 paste0(ct_id, " ", dir_label, " - GO:BP"),
                                 CFG$top_n_show, ci$color)
    }
    if (!is.null(res$KEGG)) {
      plots$KEGG <- build_dotplot(res$KEGG,
                                   paste0(ct_id, " ", dir_label, " - KEGG"),
                                   CFG$top_n_show, ci$color)
    }
    if (!is.null(res$Reactome)) {
      plots$Reactome <- build_dotplot(res$Reactome,
                                       paste0(ct_id, " ", dir_label, " - Reactome"),
                                       CFG$top_n_show, ci$color)
    }
    
    if (length(plots) > 0) {
      combined <- wrap_plots(plots, ncol = 1)
      out_pdf <- file.path(CFG$out_dir, "figures",
                           paste0(ct_id, "_", dir_label, "_dotplot.pdf"))
      ggsave(out_pdf, combined,
             width = 9, height = 3 + length(plots) * 3.5,
             limitsize = FALSE)
      cat_log("  saved dotplot: ", out_pdf, "\n")
    }
    
    # Summary row
    summary_rows[[paste0(ct_id, "_", dir_label)]] <- data.frame(
      contrast = ct_id,
      direction = dir_label,
      n_genes_input = nrow(df),
      n_GO_BP_enriched = if (!is.null(res$GO)) sum(res$GO$ontology == "BP") else 0,
      n_GO_MF_enriched = if (!is.null(res$GO)) sum(res$GO$ontology == "MF") else 0,
      n_GO_CC_enriched = if (!is.null(res$GO)) sum(res$GO$ontology == "CC") else 0,
      n_KEGG_enriched  = if (!is.null(res$KEGG)) nrow(res$KEGG) else 0,
      n_Reactome_enriched = if (!is.null(res$Reactome)) nrow(res$Reactome) else 0,
      top_BP_term = if (!is.null(res$GO) && any(res$GO$ontology == "BP"))
                    res$GO$Description[res$GO$ontology == "BP"][1] else NA,
      top_KEGG_term = if (!is.null(res$KEGG) && nrow(res$KEGG) > 0)
                      res$KEGG$Description[1] else NA,
      stringsAsFactors = FALSE
    )
  }
}

# ============ Summary ============
sec("Summary")

if (length(summary_rows) > 0) {
  summary_df <- bind_rows(summary_rows)
  write_csv(summary_df, file.path(CFG$out_dir, "pathway_summary.csv"))
  cat_log("\nPathway enrichment summary:\n")
  print(as.data.frame(summary_df))
}

# ============ Combined cross-contrast pathway plot ============
sec("Cross-contrast comparison")

# For UP only: pull top 5 BP per contrast, build comparison plot
combined_top <- list()
for (key in names(all_results)) {
  parts <- strsplit(key, "_")[[1]]
  ct_id <- parts[1]; dir_label <- parts[2]
  if (dir_label != "UP") next
  res <- all_results[[key]]
  if (is.null(res$GO)) next
  bp <- res$GO %>%
    filter(ontology == "BP") %>%
    arrange(p.adjust) %>%
    head(5) %>%
    select(Description, p.adjust, Count) %>%
    mutate(contrast = ct_id,
           organ = CFG$contrasts[[ct_id]]$organ,
           neg_log10_padj = -log10(p.adjust))
  combined_top[[ct_id]] <- bp
}

if (length(combined_top) > 0) {
  comb_df <- bind_rows(combined_top)
  comb_df$contrast_label <- paste0(comb_df$contrast, " (", comb_df$organ, ")")
  
  p_combined <- ggplot(comb_df,
                       aes(x = contrast_label, y = Description,
                           size = Count, color = neg_log10_padj)) +
    geom_point() +
    scale_color_gradient(low = "lightgray", high = "darkred",
                         name = "-log10(p.adj)") +
    scale_size(range = c(2, 8), name = "Count") +
    labs(title = "Top 5 GO:BP enriched per contrast (UP genes only)",
         x = "contrast", y = NULL) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 8),
          plot.title = element_text(hjust = 0.5))
  
  ggsave(file.path(CFG$out_dir, "figures", "combined_top_pathways.pdf"),
         p_combined, width = 12, height = 10, limitsize = FALSE)
  cat_log("  saved combined plot\n")
}

cat_log("\n=== DONE: ", format(Sys.time()), " ===\n")
cat_log("Output: ", CFG$out_dir, "\n\n")
