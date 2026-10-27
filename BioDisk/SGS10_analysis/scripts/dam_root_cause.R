#!/usr/bin/env Rscript
# 诊断 Microglia_DAM 分类错误的根因

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

PROJECT_ROOT <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
setwd(PROJECT_ROOT)
source("scripts/diff_common.R")

# === Part 1: 看 brain marker library 里 Microglia_DAM 的定义 ===
cat("========================================\n")
cat("Part 1: v4 brain marker library Microglia_DAM definition\n")
cat("========================================\n\n")

if (file.exists("scripts/marker_lib_brain_v4.R")) {
  source("scripts/marker_lib_brain_v4.R")

  dam_markers <- NULL
  # 试几个常见变量名
  for (var_name in c("BRAIN_MARKERS","brain_markers","BRAIN_MARKER_LIB",
                     "MARKERS_BRAIN","markers_brain")) {
    if (exists(var_name)) {
      ml <- get(var_name)
      if (is.list(ml) && "Microglia_DAM" %in% names(ml)) {
        dam_markers <- ml[["Microglia_DAM"]]
        cat(sprintf("Found in variable: %s\n", var_name))
        cat(sprintf("Microglia_DAM marker panel (%d genes):\n", length(dam_markers)))
        print(dam_markers)
        break
      }
    }
  }

  if (is.null(dam_markers)) {
    cat("Could not auto-find DAM markers, listing all brain marker variables:\n")
    vars <- ls()
    for (v in vars) {
      val <- tryCatch(get(v), error = function(e) NULL)
      if (is.list(val) && length(val) > 0 && any(grepl("Microglia", names(val)))) {
        cat(sprintf("\n--- Variable: %s (list) ---\n", v))
        cat(sprintf("  keys: %s\n", paste(names(val), collapse=", ")))
        if ("Microglia_DAM" %in% names(val)) {
          dam_markers <- val[["Microglia_DAM"]]
          cat(sprintf("  Microglia_DAM genes: %s\n", paste(dam_markers, collapse=", ")))
        }
      }
    }
  }
}

# === Part 2: 看 v3.2 把这些 cell 叫什么 ===
cat("\n\n========================================\n")
cat("Part 2: v3.2 original labels of v4 Microglia_DAM cells\n")
cat("========================================\n")

samples_brain <- c("brain3","brain8","nao56","bla_brain","con_brain")
all_trace <- list()

for (s in samples_brain) {
  fp <- tryCatch(resolve_rds_path(s, "v4"), error = function(e) NULL)
  if (is.null(fp)) next

  obj <- readRDS(fp)
  md <- obj@meta.data

  # filter v4 Microglia_DAM
  if ("celltype_refined" %in% colnames(md)) {
    dam_mask <- md$celltype_refined == "Microglia_DAM"
  } else {
    next
  }
  n_dam <- sum(dam_mask, na.rm = TRUE)
  if (n_dam == 0) next

  cat(sprintf("\n--- %s (%d DAM cells) ---\n", s, n_dam))

  # v3.2 original labels
  if ("celltype_v32" %in% colnames(md)) {
    v32_tab <- table(md$celltype_v32[dam_mask], useNA="ifany")
    v32_tab <- sort(v32_tab, decreasing = TRUE)
    cat("  v3.2 original labels:\n")
    for (i in seq_len(min(10, length(v32_tab)))) {
      cat(sprintf("    %-30s %5d  (%.1f%%)\n",
                  names(v32_tab)[i], v32_tab[i],
                  v32_tab[i] / n_dam * 100))
    }
  }

  # v4 refinement source
  if ("celltype_refinement_source" %in% colnames(md)) {
    src_tab <- table(md$celltype_refinement_source[dam_mask], useNA="ifany")
    cat("\n  v4 refinement source:\n")
    for (i in seq_along(src_tab)) {
      cat(sprintf("    %-30s %5d  (%.1f%%)\n",
                  names(src_tab)[i], src_tab[i],
                  src_tab[i] / n_dam * 100))
    }
  }

  # v4 rank (rank = where was v3.2 name in v4 score? 1=top agreement)
  if ("celltype_v32_rank" %in% colnames(md)) {
    rank_tab <- table(md$celltype_v32_rank[dam_mask], useNA="ifany")
    cat("\n  v4 rank of v3.2 label (in score):\n")
    for (i in seq_along(rank_tab)) {
      cat(sprintf("    rank %-5s %5d  (%.1f%%)\n",
                  names(rank_tab)[i], rank_tab[i],
                  rank_tab[i] / n_dam * 100))
    }
  }

  # Seurat cluster (原始聚类)
  if ("seurat_clusters" %in% colnames(md)) {
    cl_tab <- table(md$seurat_clusters[dam_mask])
    cl_tab <- sort(cl_tab, decreasing=TRUE)
    cat("\n  Seurat cluster (top 5):\n")
    for (i in seq_len(min(5, length(cl_tab)))) {
      cat(sprintf("    cluster %-5s %5d  (%.1f%%)\n",
                  names(cl_tab)[i], cl_tab[i],
                  cl_tab[i] / n_dam * 100))
    }
  }

  # QC metrics
  cat("\n  QC (DAM vs Homeostatic as reference):\n")
  qc_cols <- intersect(c("nFeature_RNA","nCount_RNA","percent.mt","percent.mito"),
                       colnames(md))
  for (q in qc_cols) {
    dam_med <- median(md[[q]][dam_mask], na.rm=TRUE)
    hmo_mask <- md$celltype_refined == "Microglia_Homeostatic"
    hmo_med <- median(md[[q]][hmo_mask], na.rm=TRUE)
    cat(sprintf("    %-15s DAM=%.1f  Homeostatic=%.1f  (ratio %.2f)\n",
                q, dam_med, hmo_med, dam_med / hmo_med))
  }

  rm(obj); gc(verbose = FALSE)
}

cat("\n\n========================================\n")
cat("Part 3: check which DAM markers are driving score\n")
cat("========================================\n")

# 对 brain3 Microglia_DAM, 测每个 DAM marker gene 的表达率
if (!is.null(dam_markers) && length(dam_markers) > 0) {
  obj <- readRDS(resolve_rds_path("brain3","v4"))
  obj <- tryCatch(JoinLayers(obj), error = function(e) obj)
  if (!"data" %in% Layers(obj[["RNA"]])) {
    obj <- NormalizeData(obj, verbose = FALSE)
  }

  md <- obj@meta.data
  dam_cells <- colnames(obj)[md$celltype_refined == "Microglia_DAM"]
  hmo_cells <- colnames(obj)[md$celltype_refined == "Microglia_Homeostatic"]

  dam_sub <- subset(obj, cells = dam_cells)
  hmo_sub <- subset(obj, cells = hmo_cells)

  dam_exp <- GetAssayData(dam_sub, assay="RNA", layer="data")
  hmo_exp <- GetAssayData(hmo_sub, assay="RNA", layer="data")

  genes_in <- intersect(dam_markers, rownames(dam_exp))
  cat(sprintf("\nDAM marker panel gene expression (DAM subset vs Homeostatic):\n"))
  cat(sprintf("%-20s %10s %10s %10s %10s\n",
              "gene", "DAM_pct", "DAM_mean", "Hmo_pct", "Hmo_mean"))
  for (g in genes_in) {
    d_vals <- dam_exp[g, ]
    h_vals <- hmo_exp[g, ]
    cat(sprintf("%-20s %9.1f%% %10.3f %9.1f%% %10.3f\n",
                g, mean(d_vals > 0) * 100, mean(d_vals),
                mean(h_vals > 0) * 100, mean(h_vals)))
  }

  # 哪些 DAM panel gene 在 DAM > Homeostatic (= pro-DAM, 这些驱动了 score)
  # 哪些 gene 在 DAM < Homeostatic (= false signals, 这些应该更高反而没有)
  cat("\n\nWhich genes drive the 'DAM' score classification:\n")
  for (g in genes_in) {
    d_mean <- mean(dam_exp[g, ])
    h_mean <- mean(hmo_exp[g, ])
    ratio <- d_mean / (h_mean + 1e-6)
    flag <- if (ratio > 1.5) "[DRIVE score]"
            else if (ratio < 0.5) "[OPPOSITE - real microglia marker missing]"
            else "[~same]"
    cat(sprintf("  %-20s DAM/Hmo = %.2f  %s\n", g, ratio, flag))
  }

  rm(obj, dam_sub, hmo_sub); gc(verbose = FALSE)
}
