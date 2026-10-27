#!/usr/bin/env Rscript
source("/mnt/BioDisk/wangmingyuan/SGS10_analysis/scripts/fig_save.R")  # [restyle] 三格式输出器
# ============================================================
# regenerate_cellchat_pathway_plots.R  (v3)
# ============================================================
#
# v3 fixes vs v2:
#   - Remove title.name= from netVisual_aggregate (Circle & Chord)
#     — Conflict: CellChat internally passes title.name to plot.igraph
#     — Solution: use title() AFTER the plot call (works for base graphics)
#   - Drop Chord's title.name too (same issue)
#
# Expected: 4/4 panels per pathway (Circle, Chord, Heatmap, Contribution)
# ============================================================

suppressPackageStartupMessages({
  library(CellChat)
  library(ggplot2)
})

PROJECT_ROOT <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
setwd(PROJECT_ROOT)

# [restyle] 可选: 只跑指定 organ / pathway, 省时.
#   全跑:           Rscript scripts/regenerate_cellchat_pathway_plots.R
#   只跑脑某通路:   Rscript scripts/regenerate_cellchat_pathway_plots.R brain SPP1
.args        <- commandArgs(trailingOnly = TRUE)
ONLY_ORGAN   <- if (length(.args) >= 1) .args[1] else ""
ONLY_PATHWAY <- if (length(.args) >= 2) .args[2] else ""

CELLCHAT_DIR <- "taskB/cellchat"
OUT_DIR <- file.path(CELLCHAT_DIR, "pathway_plots")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

sanitize_fname <- function(s) gsub("[/\\\\ :]", "_", s)

safe_plot <- function(plot_expr, panel_name) {
  tryCatch({ plot_expr(); TRUE },
    error = function(e) {
      message(sprintf("      [%s FAIL] %s", panel_name, substr(e$message, 1, 60)))
      FALSE
    })
}

# [restyle] circle 节点名缩写表(脑); 未列出的保持原名. 用于整圈消叠字
ABBREV <- c(
  Astrocyte_Gfap        = "Ast_Gfap",
  Astrocyte_Resting     = "Ast_Rest",
  Endothelial           = "Endo",
  Ependymal             = "Epen",
  Excitatory_L23        = "Ex_L23",
  Excitatory_L45        = "Ex_L45",
  Excitatory_L6         = "Ex_L6",
  Hpc_CA1               = "CA1",
  Hpc_CA3               = "CA3",
  Hpc_DG                = "DG",
  Inhibitory_Meis2      = "In_Meis2",
  Inhibitory_Pvalb      = "In_Pvalb",
  Inhibitory_Sst        = "In_Sst",
  Inhibitory_Vip        = "In_Vip",
  Microglia_Homeostatic = "MG_Homeo",
  Microglia_other       = "MG_other",
  Oligo_Mature          = "Oligo",
  Pericytes             = "Peri"
)
.ren <- function(x) { i <- match(x, names(ABBREV)); ifelse(is.na(i), x, ABBREV[i]) }

# 在 cellchat 对象各处把全名替换成缩写, 返回副本(只动维名/levels, 不动顺序->颜色不变)
abbrev_cellchat <- function(cc) {
  try(levels(cc@idents) <- .ren(levels(cc@idents)), silent = TRUE)
  for (nm in names(cc@net)) {
    a <- cc@net[[nm]]
    if (!is.null(dim(a)) && !is.null(dimnames(a)) && length(dimnames(a)) >= 2) {
      dn <- dimnames(a); dn[[1]] <- .ren(dn[[1]]); dn[[2]] <- .ren(dn[[2]])
      dimnames(a) <- dn; cc@net[[nm]] <- a
    }
  }
  if (!is.null(cc@netP$prob)) {
    dn <- dimnames(cc@netP$prob); dn[[1]] <- .ren(dn[[1]]); dn[[2]] <- .ren(dn[[2]])
    dimnames(cc@netP$prob) <- dn
  }
  cc
}

plot_pathway <- function(cellchat, pw_name, outfile) {
  cat(sprintf("    [%s] ", pw_name))
  
  pdf(outfile, width = 10, height = 8)
  n_panels_ok <- 0
  
  # Panel 1: Circle (no title.name; netVisual_aggregate will auto-title)
  if (safe_plot(function() {
    par(mfrow = c(1, 1), mar = c(2, 2, 3, 2))
    netVisual_aggregate(
      abbrev_cellchat(cellchat), signaling = pw_name,
      layout = "circle",
      vertex.label.cex = 1.5,   # [restyle] 0.6->1.5
      edge.width.max = 5
    )
  }, "Circle")) n_panels_ok <- n_panels_ok + 1
  
  # Panel 2: Chord (用 netVisual_chord_cell, 默认无标题)
  if (safe_plot(function() {
    par(mfrow = c(1, 1), mar = c(2, 2, 3, 2))
    netVisual_chord_cell(
      cellchat, signaling = pw_name,
      lab.cex = 1.5,       # [restyle] 扇区标签字号
      title.name = ""      # [restyle] 彻底去标题
    )
  }, "Chord")) n_panels_ok <- n_panels_ok + 1
  
  # Panel 3: Heatmap
  if (safe_plot(function() {
    ht <- netVisual_heatmap(
      cellchat, signaling = pw_name,
      color.heatmap = "Reds",
      title.name = sprintf("%s - Communication Probability", pw_name)
    )
    ComplexHeatmap::draw(ht)
  }, "Heatmap")) n_panels_ok <- n_panels_ok + 1
  
  # Panel 4: Contribution (ggplot - supports ggtitle safely)
  if (safe_plot(function() {
    p <- netAnalysis_contribution(cellchat, signaling = pw_name) +
      ggtitle(sprintf("%s - LR Pair Contribution", pw_name))
    print(p)
  }, "Contribution")) n_panels_ok <- n_panels_ok + 1
  
  dev.off()

  # [restyle] 额外: Circle / Chord 各自出 PDF+TIFF+PNG@300 (拼图用, 节点标签已放大)
  .cc_base <- sub("\\.pdf$", "", outfile)
  try(save_fig_all(drawfun = function() {
    par(mar = c(2, 2, 2, 2))
    netVisual_aggregate(abbrev_cellchat(cellchat), signaling = pw_name, layout = "circle",
                        vertex.label.cex = 1.5, edge.width.max = 5)
  }, file = paste0(.cc_base, "_circle"), width = 12, height = 10), silent = TRUE)
  try(save_fig_all(drawfun = function() {
    par(mar = c(2, 2, 2, 2))
    netVisual_chord_cell(cellchat, signaling = pw_name,
                         lab.cex = 1.5, title.name = "")
  }, file = paste0(.cc_base, "_chord"), width = 12, height = 10), silent = TRUE)

  sz <- file.size(outfile)
  if (n_panels_ok == 0 || sz < 5000) {
    cat(sprintf("FAIL (%d panels, size=%d)\n", n_panels_ok, sz))
    if (file.exists(outfile)) file.remove(outfile)
    return(FALSE)
  }
  cat(sprintf("OK (%d/4 panels, %.1f KB)\n", n_panels_ok, sz/1024))
  return(TRUE)
}

# =====================================================
organs <- list(
  list(name = "liver", rds = file.path(CELLCHAT_DIR, "liver_tumor_cellchat.rds")),
  list(name = "brain", rds = file.path(CELLCHAT_DIR, "brain_tumor_cellchat.rds"))
)

summary_rows <- list()

if (nzchar(ONLY_ORGAN)) {
  organs <- Filter(function(o) o$name == ONLY_ORGAN, organs)
  cat(sprintf(">>> 只跑 organ = %s\n", ONLY_ORGAN))
}

for (org in organs) {
  cat(sprintf("\n############################################################\n"))
  cat(sprintf("### === %s ===\n", toupper(org$name)))
  cat(sprintf("############################################################\n"))
  
  t0 <- Sys.time()
  cat(sprintf("Loading %s ... ", org$rds))
  cellchat <- readRDS(org$rds)
  cat(sprintf("done (%.1fs, %d cell types)\n", 
              as.numeric(difftime(Sys.time(), t0, units = "secs")),
              length(levels(cellchat@idents))))
  
  pathways <- cellchat@netP$pathways
  if (nzchar(ONLY_PATHWAY)) {
    pathways <- intersect(pathways, ONLY_PATHWAY)
    if (length(pathways) == 0)
      cat(sprintf("!!! %s 中无通路 '%s', 可用: %s\n",
                  org$name, ONLY_PATHWAY,
                  paste(head(cellchat@netP$pathways, 20), collapse = ", ")))
  }
  cat(sprintf("Active pathways: %d\n", length(pathways)))
  
  org_out <- file.path(OUT_DIR, org$name)
  dir.create(org_out, showWarnings = FALSE, recursive = TRUE)
  
  n_ok <- 0; n_fail <- 0
  failed_pathways <- c()
  for (i in seq_along(pathways)) {
    pw <- pathways[i]
    fname <- sanitize_fname(pw)
    outfile <- file.path(org_out, sprintf("%s.pdf", fname))
    
    cat(sprintf("  [%d/%d] ", i, length(pathways)))
    
    while (!is.null(dev.list())) try(dev.off(), silent = TRUE)
    
    ok <- tryCatch(
      plot_pathway(cellchat, pw, outfile),
      error = function(e) {
        while (!is.null(dev.list())) try(dev.off(), silent = TRUE)
        cat(sprintf("ERR %s\n", substr(e$message, 1, 50)))
        if (file.exists(outfile)) file.remove(outfile)
        FALSE
      }
    )
    if (ok) n_ok <- n_ok + 1 else { 
      n_fail <- n_fail + 1
      failed_pathways <- c(failed_pathways, pw)
    }
  }
  
  cat(sprintf("\n=== %s: %d OK, %d FAILED ===\n", toupper(org$name), n_ok, n_fail))
  if (length(failed_pathways) > 0) {
    cat(sprintf("Failed: %s\n", paste(failed_pathways, collapse = ", ")))
  }
  
  summary_rows[[org$name]] <- data.frame(
    organ = org$name, n_pathway = length(pathways),
    n_ok = n_ok, n_fail = n_fail, stringsAsFactors = FALSE
  )
  
  rm(cellchat); gc(verbose = FALSE)
}

final_summary <- do.call(rbind, summary_rows)
write.csv(final_summary, file.path(OUT_DIR, "regenerate_summary.csv"), row.names = FALSE)

cat(sprintf("\n\n############################################################\n"))
cat(sprintf("### PATHWAY PDF REGENERATION COMPLETE\n"))
cat(sprintf("############################################################\n"))
print(final_summary)
cat(sprintf("\nOutput: %s\n", OUT_DIR))

for (org in organs) {
  org_out <- file.path(OUT_DIR, org$name)
  pdfs <- list.files(org_out, pattern = "\\.pdf$", full.names = TRUE)
  if (length(pdfs) > 0) {
    sz <- file.info(pdfs)$size
    df <- data.frame(file = basename(pdfs), size_KB = round(sz / 1024, 1))
    df <- df[order(-df$size_KB), ]
    cat(sprintf("\n=== %s TOP 10 ===\n", toupper(org$name)))
    print(head(df, 10), row.names = FALSE)
  }
}
