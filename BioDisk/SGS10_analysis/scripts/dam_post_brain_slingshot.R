#!/usr/bin/env Rscript
# Slingshot + tradeSeq pseudotime (v4: 彻底重构)
# 关键改动:
#   1. Harmony 真的能 work (多重 fallback)
#   2. Slingshot 限制最多 2 lineages (end.clus 参数)
#   3. tradeSeq 默认 without conditions (避免 predict NULL), 单独尝试 conditionTest

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(slingshot)
  library(tradeSeq)
  library(SingleCellExperiment)
  library(BiocParallel)
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  library(gridExtra)
})

PROJECT_ROOT <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
setwd(PROJECT_ROOT)

OUT_DIR <- file.path(PROJECT_ROOT, "pseudotime_brain_slingshot")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
PLOT_DIR <- file.path(OUT_DIR, "plots")
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

get_clean_brain_rds <- function(sample) {
  fp <- file.path(PROJECT_ROOT, sample, "mouse",
                  sprintf("%s_manual_annotated_v4_clean.rds", sample))
  if (!file.exists(fp)) stop(sprintf("Not found: %s", fp))
  fp
}

samples_brain <- c("brain3","brain8","nao56","bla_brain","con_brain")
SAMPLE_GROUP <- c(brain3="Tumor", brain8="Tumor", nao56="Tumor",
                  bla_brain="Blank", con_brain="Control")

CELLTYPE_GROUPS <- list(
  Microglia = c("Microglia_Homeostatic","Microglia_DAM_strict"),
  Astrocyte = c("Astrocyte_Resting","Astrocyte_Gfap"),
  Pericytes = c("Pericytes")
)

KNOWN_MARKERS <- list(
  Microglia = c("P2ry12","Tmem119","Cx3cr1","Csf1r","Hexb","Apoe","Trem2",
                "Spp1","Smad3","Elmo1","Mertk","H2-D1","B2m","Cst7","Lpl",
                "Il6st","Lrp1","Cd180","Ppard","Tspan18","Nav2","Basp1"),
  Astrocyte = c("Sox9","Gfap","Slc1a2","Slc1a3","Aqp4","Cdh4","Grm3","Atp1a2",
                "Kcnj10","Fasn","Aldh1l1","Plpp3","Sorcs2","Bach2","Mtss2",
                "Ubc","Xylt1","Neat1"),
  Pericytes = c("Cldn5","Cdh5","Klf4","Notch3","Cspg4","Pdgfrb","Vtn","H2-D1",
                "H2-K1","Gbp4","Egr1","Foxf1","Cd248","Bsg","Sik1","Ubc",
                "Jund","Nr4a1")
)

# ============================================================
# 提取 subset data (避免 Seurat 5 subset 的 metadata sync 问题)
# ============================================================
extract_subset_data <- function(o, cells_keep_idx, sample_name, condition_str) {
  cells_keep <- colnames(o)[cells_keep_idx]
  counts_mat <- GetAssayData(o, assay = "RNA", layer = "counts")[, cells_keep]
  md_keep <- o@meta.data[cells_keep, c("celltype_refined","celltype_broad",
                                        "nCount_RNA","nFeature_RNA",
                                        "percent.mt","seurat_clusters"),
                         drop = FALSE]
  md_keep$sample_id <- sample_name
  md_keep$condition <- condition_str
  md_keep$orig_barcode <- cells_keep
  new_names <- paste0(sample_name, "_", cells_keep)
  colnames(counts_mat) <- new_names
  rownames(md_keep) <- new_names
  list(counts = counts_mat, meta = md_keep)
}

subset_celltype_across_samples <- function(target_celltypes, samples) {
  cat(sprintf("\n  Subsetting %s across %d samples...\n",
              paste(target_celltypes, collapse="+"), length(samples)))
  all_data <- list()
  for (s in samples) {
    fp <- get_clean_brain_rds(s)
    o <- readRDS(fp)
    o <- tryCatch(JoinLayers(o), error = function(e) o)
    cells_keep_idx <- which(o$celltype_refined %in% target_celltypes)
    cat(sprintf("    %s: %d cells\n", s, length(cells_keep_idx)))
    if (length(cells_keep_idx) < 10) { rm(o); gc(verbose=FALSE); next }
    cond_str <- unname(SAMPLE_GROUP[s])
    all_data[[s]] <- extract_subset_data(o, cells_keep_idx, s, cond_str)
    rm(o); gc(verbose = FALSE)
  }
  if (length(all_data) == 0) return(NULL)
  all_counts <- do.call(cbind, lapply(all_data, function(x) x$counts))
  all_meta <- do.call(rbind, lapply(all_data, function(x) x$meta))
  rm(all_data); gc(verbose = FALSE)
  cat(sprintf("  Creating Seurat object: %d cells x %d genes\n",
              ncol(all_counts), nrow(all_counts)))
  merged <- CreateSeuratObject(counts = all_counts, meta.data = all_meta,
                               project = "brain_pseudotime")
  rm(all_counts, all_meta); gc(verbose = FALSE)
  cat("  Conditions:\n"); print(table(merged$condition))
  cat("  Samples:\n"); print(table(merged$sample_id))
  merged
}

# ============================================================
# Harmony with verification + fallback
# ============================================================
run_harmony_with_verification <- function(seu) {
  cat("  [Harmony] Attempting harmony integration...\n")
  
  # 尝试 1: 最简单调用
  seu_h <- tryCatch({
    res <- harmony::RunHarmony(seu, group.by.vars = "sample_id", verbose = TRUE)
    res
  }, error = function(e) {
    cat(sprintf("  [Harmony] Attempt 1 failed: %s\n", e$message))
    NULL
  })
  
  # 尝试 2: 带具体参数
  if (is.null(seu_h) || !"harmony" %in% Reductions(seu_h)) {
    cat("  [Harmony] Attempt 2 with explicit params...\n")
    seu_h <- tryCatch({
      harmony::RunHarmony(seu, 
                          group.by.vars = "sample_id",
                          reduction.use = "pca",
                          dims.use = 1:30,
                          reduction.save = "harmony",
                          verbose = TRUE)
    }, error = function(e) {
      cat(sprintf("  [Harmony] Attempt 2 failed: %s\n", e$message))
      NULL
    })
  }
  
  # 验证: harmony reduction 必须存在 且 embedding 维度正确
  if (!is.null(seu_h) && "harmony" %in% Reductions(seu_h)) {
    harm_dim <- dim(Embeddings(seu_h, "harmony"))
    cat(sprintf("  [Harmony] ✓ SUCCESS: harmony embedding %d x %d\n",
                harm_dim[1], harm_dim[2]))
    return(list(seu = seu_h, red = "harmony"))
  } else {
    cat("  [Harmony] ✗ FAILED - falling back to PCA (no batch correction!)\n")
    return(list(seu = seu, red = "pca"))
  }
}

# ============================================================
# Run Slingshot (限制最多 2 lineages)
# ============================================================
run_slingshot_limited <- function(sce, cluster_col, start_cluster, 
                                    pick_end_cluster = TRUE) {
  cat("  [Slingshot] Running...\n")
  
  # 基础 slingshot
  sce_base <- slingshot(sce,
                        clusterLabels = cluster_col,
                        reducedDim = "UMAP",
                        start.clus = start_cluster,
                        approx_points = 150,
                        extend = "n",
                        stretch = 0)
  
  n_lin <- length(slingLineages(sce_base))
  cat(sprintf("  [Slingshot] Base run: %d lineages\n", n_lin))
  
  # 如果 lineage 太多 (> 2), 自动限制
  if (n_lin > 2 && pick_end_cluster) {
    # 选出 2 个最可能是 "activated endpoint" 的 cluster
    # 策略: 最远离 start cluster (根据 cluster centroid 在 UMAP 的距离)
    cluster_labels <- as.character(sce[[cluster_col]])
    unique_clusters <- setdiff(unique(cluster_labels), start_cluster)
    
    # 算每个 cluster 在 UMAP 空间中与 start cluster centroid 的距离
    umap_emb <- reducedDim(sce, "UMAP")
    start_mask <- cluster_labels == start_cluster
    start_centroid <- colMeans(umap_emb[start_mask, , drop=FALSE])
    
    dist_to_start <- sapply(unique_clusters, function(cl) {
      cl_mask <- cluster_labels == cl
      if (sum(cl_mask) < 5) return(0)
      cl_centroid <- colMeans(umap_emb[cl_mask, , drop=FALSE])
      sqrt(sum((cl_centroid - start_centroid)^2))
    })
    
    # 取距离最远的 2 个
    top2_end <- names(sort(dist_to_start, decreasing = TRUE)[1:min(2, length(dist_to_start))])
    cat(sprintf("  [Slingshot] Limiting to end.clus: %s\n",
                paste(top2_end, collapse=", ")))
    
    sce_limited <- tryCatch({
      slingshot(sce,
                clusterLabels = cluster_col,
                reducedDim = "UMAP",
                start.clus = start_cluster,
                end.clus = top2_end,
                approx_points = 150,
                extend = "n",
                stretch = 0)
    }, error = function(e) {
      cat(sprintf("  [Slingshot] Limited run failed: %s, using base\n", e$message))
      sce_base
    })
    
    n_lin_new <- length(slingLineages(sce_limited))
    cat(sprintf("  [Slingshot] After limiting: %d lineages\n", n_lin_new))
    return(sce_limited)
  }
  
  sce_base
}

# ============================================================
# Main run function
# ============================================================
run_pseudotime <- function(seu, group_name, n_knots = 6) {
  cat(sprintf("\n=== %s: pipeline start ===\n", group_name))
  t0 <- proc.time()
  
  # === Step 1: Preprocessing ===
  seu <- NormalizeData(seu, verbose = FALSE)
  seu <- FindVariableFeatures(seu, nfeatures = 2000, verbose = FALSE)
  seu <- ScaleData(seu, verbose = FALSE)
  seu <- RunPCA(seu, npcs = 30, verbose = FALSE)
  
  # === Step 2: Harmony (verified) ===
  harm_result <- run_harmony_with_verification(seu)
  seu <- harm_result$seu
  red_use <- harm_result$red
  
  # === Step 3: UMAP + clustering ===
  seu <- RunUMAP(seu, reduction = red_use, dims = 1:20, verbose = FALSE)
  seu <- FindNeighbors(seu, reduction = red_use, dims = 1:20, verbose = FALSE)
  
  cluster_res <- if (group_name == "Pericytes") 0.2 else 0.3
  seu <- FindClusters(seu, resolution = cluster_res, verbose = FALSE)
  
  n_clusters <- length(unique(seu$seurat_clusters))
  cat(sprintf("\n  Clusters: %d\n", n_clusters))
  cat("  Cluster x Condition:\n")
  print(table(seu$seurat_clusters, seu$condition))
  cat("  Cluster x Sample (batch check):\n")
  print(table(seu$seurat_clusters, seu$sample_id))
  
  # === Step 4: Convert to SCE ===
  sce <- as.SingleCellExperiment(seu)
  reducedDim(sce, "UMAP") <- Embeddings(seu, "umap")
  
  # === Step 5: Slingshot (limited to ≤ 2 lineages) ===
  cluster_col <- "seurat_clusters"
  cluster_cond_tab <- table(sce[[cluster_col]], sce$condition)
  
  if ("Control" %in% colnames(cluster_cond_tab)) {
    ctrl_pct <- cluster_cond_tab[, "Control"] / rowSums(cluster_cond_tab)
    start_cluster <- names(ctrl_pct)[which.max(ctrl_pct)]
  } else {
    start_cluster <- names(table(sce[[cluster_col]]))[1]
  }
  cat(sprintf("\n  Start cluster (Control-enriched): %s\n", start_cluster))
  
  sce <- run_slingshot_limited(sce, cluster_col, start_cluster, 
                                pick_end_cluster = TRUE)
  
  n_lineages <- length(slingLineages(sce))
  cat(sprintf("\n  === Final: %d lineages ===\n", n_lineages))
  for (i in seq_len(n_lineages)) {
    cat(sprintf("    Lineage %d: %s\n", i,
                paste(slingLineages(sce)[[i]], collapse=" -> ")))
  }
  
  pseudotime_mat <- slingPseudotime(sce, na = FALSE)
  cell_weights <- slingCurveWeights(sce)
  
  # Check condition balance in each lineage
  cat("\n  Condition balance per lineage (valid cells with weight > 0.5):\n")
  cond_balance_ok <- TRUE
  for (i in seq_len(ncol(pseudotime_mat))) {
    valid <- !is.na(pseudotime_mat[, i]) & cell_weights[, i] > 0.5
    cond_tab <- table(sce$condition[valid])
    cat(sprintf("    Lineage %d: ", i))
    print(cond_tab)
    if (any(cond_tab < 10)) {
      cat(sprintf("    ⚠️ Lineage %d has very few cells in some condition\n", i))
      cond_balance_ok <- FALSE
    }
  }
  
  # === Step 6: Save pseudotime CSV ===
  pt_df <- data.frame(
    cell = colnames(sce),
    sample = sce$sample_id,
    condition = sce$condition,
    celltype = sce$celltype_refined,
    cluster = sce[[cluster_col]]
  )
  for (i in seq_len(ncol(pseudotime_mat))) {
    pt_df[[sprintf("pseudotime_lineage_%d", i)]] <- pseudotime_mat[, i]
  }
  write.csv(pt_df, file.path(OUT_DIR, sprintf("%s_pseudotime.csv", group_name)),
            row.names = FALSE)
  
  # === Step 7: Plots (trajectory + distribution) ===
  cond_colors <- c("Tumor"="red","Blank"="orange","Control"="blue")
  
  pdf(file.path(PLOT_DIR, sprintf("%s_01_umap_trajectory.pdf", group_name)),
      width = 8, height = 7)
  plot(reducedDim(sce, "UMAP"),
       col = cond_colors[sce$condition],
       pch = 16, cex = 0.5, asp = 1,
       main = sprintf("%s: trajectory (%s)", group_name, 
                       ifelse(red_use == "harmony", "Harmony", "PCA-only")))
  lines(SlingshotDataSet(sce), lwd = 2, col = "black")
  legend("topright", legend = names(cond_colors),
         col = cond_colors, pch = 16)
  dev.off()
  
  if (ncol(pseudotime_mat) > 0) {
    pdf(file.path(PLOT_DIR, sprintf("%s_02_umap_pseudotime.pdf", group_name)),
        width = 8, height = 7)
    pt1 <- pseudotime_mat[, 1]
    pt_valid <- !is.na(pt1)
    if (sum(pt_valid) > 10) {
      pt_colors <- colorRampPalette(c("blue","white","red"))(100)
      pt_norm <- (pt1 - min(pt1, na.rm=TRUE)) /
                 (max(pt1, na.rm=TRUE) - min(pt1, na.rm=TRUE))
      pt_idx <- pmin(pmax(round(pt_norm * 99) + 1, 1), 100)
      pt_idx[!pt_valid] <- 50
      plot(reducedDim(sce, "UMAP"),
           col = pt_colors[pt_idx], pch = 16, cex = 0.5, asp = 1,
           main = sprintf("%s: lineage 1 pseudotime", group_name))
      lines(SlingshotDataSet(sce), lwd = 2)
    }
    dev.off()
  }
  
  # Pseudotime violin
  pt_long <- data.frame()
  for (i in seq_len(ncol(pseudotime_mat))) {
    valid <- !is.na(pseudotime_mat[, i]) & cell_weights[, i] > 0.5
    if (sum(valid) < 10) next
    pt_long <- rbind(pt_long, data.frame(
      pt = pseudotime_mat[valid, i],
      condition = sce$condition[valid],
      lineage = sprintf("Lineage %d", i)
    ))
  }
  if (nrow(pt_long) > 0) {
    p3 <- ggplot(pt_long, aes(x = condition, y = pt, fill = condition)) +
      geom_violin() +
      geom_boxplot(width = 0.1, fill = "white") +
      facet_wrap(~lineage) +
      scale_fill_manual(values = cond_colors) +
      theme_bw() +
      ggtitle(sprintf("%s: pseudotime by condition", group_name))
    ggsave(file.path(PLOT_DIR, sprintf("%s_03_pt_by_condition.pdf", group_name)),
           p3, width = 7, height = 5)
  }
  
  # Sample QC plot
  p_sample <- DimPlot(seu, group.by = "sample_id", reduction = "umap",
                       pt.size = 0.3) +
    ggtitle(sprintf("%s: UMAP by sample (red=%s)", group_name, red_use))
  ggsave(file.path(PLOT_DIR, sprintf("%s_06_umap_by_sample.pdf", group_name)),
         p_sample, width = 8, height = 6)
  
  # Save Seurat + SCE before tradeSeq (以防 tradeSeq 崩)
  saveRDS(seu, file.path(OUT_DIR, sprintf("%s_seurat.rds", group_name)))
  saveRDS(sce, file.path(OUT_DIR, sprintf("%s_sce_slingshot.rds", group_name)))
  
  # === Step 8: tradeSeq (ROBUST - associationTest primary) ===
  cat("\n  === tradeSeq ===\n")
  cat("  Fitting GAM (without conditions, safer)...\n")
  
  top_genes <- VariableFeatures(seu)[1:min(800, length(VariableFeatures(seu)))]
  counts <- as.matrix(GetAssayData(seu, assay = "RNA", layer = "counts")[top_genes, ])
  
  set.seed(42)
  sce_ts <- tryCatch({
    fitGAM(counts = counts,
           pseudotime = pseudotime_mat,
           cellWeights = cell_weights,
           nknots = n_knots,
           parallel = TRUE,
           BPPARAM = MulticoreParam(workers = 4))
  }, error = function(e) {
    cat(sprintf("    fitGAM failed: %s\n", e$message))
    NULL
  })
  
  if (is.null(sce_ts)) {
    cat("  tradeSeq fitGAM failed completely. Saving slingshot-only results.\n")
    elapsed <- (proc.time() - t0)["elapsed"] / 60
    cat(sprintf("  ✓ %s done in %.1f min (slingshot only)\n", group_name, elapsed))
    rm(seu, sce, counts); gc(verbose = FALSE)
    return(invisible(NULL))
  }
  
  # associationTest (没 condition, 稳定)
  cat("  Running associationTest...\n")
  assoc_res <- tryCatch(associationTest(sce_ts), error = function(e) NULL)
  if (!is.null(assoc_res)) {
    assoc_res$gene <- rownames(assoc_res)
    assoc_res$padj <- p.adjust(assoc_res$pvalue, method = "BH")
    assoc_res <- assoc_res %>% arrange(pvalue)
    write.csv(assoc_res,
              file.path(OUT_DIR, sprintf("%s_associationTest.csv", group_name)),
              row.names = FALSE)
    cat(sprintf("  ✓ associationTest: %d dynamic genes (padj<0.05)\n",
                sum(assoc_res$padj < 0.05, na.rm = TRUE)))
    
    # Known marker sanity check
    known <- KNOWN_MARKERS[[group_name]]
    if (!is.null(known)) {
      cat(sprintf("\n  === Known %s markers in trajectory ===\n", group_name))
      hits <- assoc_res[assoc_res$gene %in% known, ]
      if (nrow(hits) > 0) {
        for (i in seq_len(min(15, nrow(hits)))) {
          sig <- if (!is.na(hits$padj[i]) && hits$padj[i] < 0.05) "*" else ""
          cat(sprintf("    %-12s waldStat=%7.1f padj=%.2e %s\n",
                      hits$gene[i], hits$waldStat[i],
                      hits$padj[i], sig))
        }
      }
    }
    
    # Top 12 dynamic gene plots
    top_assoc <- head(assoc_res$gene[assoc_res$padj < 0.05], 12)
    if (length(top_assoc) > 0) {
      plots <- lapply(top_assoc, function(g) {
        tryCatch({
          plotSmoothers(sce_ts, counts = counts, gene = g, alpha = 0.3) +
            ggtitle(g)
        }, error = function(e) NULL)
      })
      plots <- plots[!sapply(plots, is.null)]
      if (length(plots) > 0) {
        pdf(file.path(PLOT_DIR, sprintf("%s_04_top_dynamic_genes.pdf", group_name)),
            width = 16, height = 10)
        gridExtra::grid.arrange(grobs = plots, ncol = 4)
        dev.off()
      }
    }
  }
  
  # 可选: 尝试 conditionTest (可能失败, 不影响主流程)
  if (cond_balance_ok) {
    cat("\n  Attempting conditionTest (optional, may fail)...\n")
    cond_factor <- factor(sce$condition, levels = c("Control","Blank","Tumor"))
    
    sce_ts_cond <- tryCatch({
      fitGAM(counts = counts,
             pseudotime = pseudotime_mat,
             cellWeights = cell_weights,
             conditions = cond_factor,
             nknots = n_knots,
             parallel = TRUE,
             BPPARAM = MulticoreParam(workers = 4))
    }, error = function(e) {
      cat(sprintf("    fitGAM w/ conditions failed: %s\n", e$message))
      NULL
    })
    
    if (!is.null(sce_ts_cond)) {
      cond_res <- tryCatch({
        conditionTest(sce_ts_cond, l2fc = log2(1.5))
      }, error = function(e) NULL)
      
      if (!is.null(cond_res)) {
        cond_res$gene <- rownames(cond_res)
        cond_res$padj <- p.adjust(cond_res$pvalue, method = "BH")
        cond_res <- cond_res %>% arrange(pvalue)
        write.csv(cond_res,
                  file.path(OUT_DIR, sprintf("%s_conditionTest.csv", group_name)),
                  row.names = FALSE)
        cat(sprintf("  ✓ conditionTest: %d genes differ between conditions\n",
                    sum(cond_res$padj < 0.05, na.rm = TRUE)))
      }
    }
  } else {
    cat("  Skipping conditionTest (condition imbalance in lineages)\n")
  }
  
  saveRDS(sce_ts, file.path(OUT_DIR, sprintf("%s_sce_tradeSeq.rds", group_name)))
  
  elapsed <- (proc.time() - t0)["elapsed"] / 60
  cat(sprintf("\n  ✓ %s done in %.1f min\n", group_name, elapsed))
  
  rm(seu, sce, sce_ts, counts); gc(verbose = FALSE)
}

# ============================================================
cat("========================================\n")
cat("Brain Slingshot + tradeSeq pseudotime (v4 rewrite)\n")
cat("========================================\n")

for (gn in names(CELLTYPE_GROUPS)) {
  cat(sprintf("\n\n############# %s #############\n", gn))
  
  merged <- subset_celltype_across_samples(
    target_celltypes = CELLTYPE_GROUPS[[gn]],
    samples = samples_brain
  )
  
  if (is.null(merged) || ncol(merged) < 50) {
    cat(sprintf("  Skip %s: too few cells\n", gn))
    next
  }
  
  run_pseudotime(merged, group_name = gn)
}

cat("\n\n========================================\n")
cat("ALL DONE\n")
cat(sprintf("Output: %s\n", OUT_DIR))
