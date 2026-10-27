#!/usr/bin/env Rscript
# ============================================================
# cellchat_niche_communication.R  --  T6
# ============================================================
# 
# 目的: 在空间组层面跑 CellChat, 识别 tumor 微环境的 L-R 网络
# 
# 输入:
#   - spatial_results/<sample>_seekspace.rds (gene expression)
#   - spatial_results/rctd/per_sample/<sample>_dominant_celltype.csv
# 
# 策略:
#   按客户确认的计划版: Liver+Brain 各 merge 3 tumor mice
#   liver_tumor:  Liver_3 + Liver_8 + gan56
#   brain_tumor:  Brain_3 + Brain_8 + nao56
#   每组独立跑 CellChat, 输出 L-R network
# 
# 处理 dual-species liver:
#   所有 mm10--- prefix 基因 rename 为 bare name
#   丢弃 GRCh38- 基因 (我们只关心 mouse host 响应)
# 
# 输出:
#   taskB/cellchat/
#     liver_tumor_cellchat.rds            CellChat object
#     liver_tumor_lr_pairs.csv            all significant LR pairs
#     liver_tumor_pathway_network.pdf     signaling pathway overview
#     liver_tumor_top_interactions.pdf    top 20 LR visualization
#     brain_tumor_* (同上)
#     summary.md
# 
# 预计: 1-2 hr (主要耗时在 CellChat 的 computeCommunProb)
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(dplyr)
  library(patchwork)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  # NMF 被 CellChat 内部用
})

options(stringsAsFactors = FALSE)
future::plan("multisession", workers = 4)
options(future.globals.maxSize = 64 * 1024^3)

PROJECT_ROOT  <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
SPATIAL_DIR   <- file.path(PROJECT_ROOT, "spatial_results")
RCTD_DIR      <- file.path(SPATIAL_DIR, "rctd", "per_sample")
OUT_DIR       <- file.path(PROJECT_ROOT, "taskB", "cellchat")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("============================================================\n")
cat("T6: CellChat on Spatial (Liver tumor + Brain tumor merged)\n")
cat(sprintf("Start: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("============================================================\n")

# ============================================================
# Tumor groups (client-confirmed)
# ============================================================
GROUPS <- list(
  liver_tumor = c("Liver_3", "Liver_8", "gan56"),
  brain_tumor = c("Brain_3", "Brain_8", "nao56")
)

# Min cells per cell type (CellChat requirement)
MIN_CELLS_PER_CT <- 10

# ============================================================
# Function: load + prep one sample for merge
# ============================================================
load_for_cellchat <- function(sample) {
  rds_path  <- file.path(SPATIAL_DIR, paste0(sample, "_seekspace.rds"))
  rctd_path <- file.path(RCTD_DIR, paste0(sample, "_dominant_celltype.csv"))
  
  if (!file.exists(rds_path) || !file.exists(rctd_path)) return(NULL)
  
  cat(sprintf("  Loading %s ... ", sample))
  t0 <- Sys.time()
  obj <- readRDS(rds_path)
  
  # Rename rownames: mm10---Gene -> Gene, drop GRCh38- (human)
  old_names <- rownames(obj)
  has_mm10  <- grepl("^mm10---", old_names)
  has_grch  <- grepl("^GRCh38-", old_names)
  
  if (any(has_mm10) || any(has_grch)) {
    # Keep only mouse features, rename to bare
    # Rebuild RNA assay with only mouse features
    mouse_features <- old_names[has_mm10 | (!has_mm10 & !has_grch)]
    new_names <- sub("^mm10---", "", mouse_features)
    
    counts <- GetAssayData(obj, layer = "counts")[mouse_features, ]
    rownames(counts) <- new_names
    
    # Handle duplicated new_names (if any)
    if (any(duplicated(new_names))) {
      cat(sprintf("\n    ⚠ Duplicated bare names after rename, keeping first\n"))
      counts <- counts[!duplicated(new_names), ]
    }
    
    # Rebuild seurat object with mouse-only features
    obj2 <- CreateSeuratObject(counts = counts, meta.data = obj@meta.data,
                               project = sample)
    obj2 <- NormalizeData(obj2, verbose = FALSE)
    rm(obj)
    obj <- obj2
    cat(sprintf("rebuilt mouse-only (%d features) ... ", nrow(obj)))
  }
  
  # Merge RCTD dominant_celltype
  rctd_df <- read.csv(rctd_path, stringsAsFactors = FALSE)
  rctd_df$barcode <- gsub('"', '', rctd_df$barcode)
  
  md <- obj@meta.data
  md$barcode <- rownames(md)
  md <- md %>%
    left_join(rctd_df[, c("barcode", "dominant_celltype")], by = "barcode")
  rownames(md) <- md$barcode
  obj@meta.data <- md
  
  # Filter: keep only spots with RCTD label
  keep <- !is.na(obj$dominant_celltype)
  obj <- obj[, keep]
  
  # Tag sample
  obj$sample_id <- sample
  
  cat(sprintf("%d spots kept (%.1fs)\n",
              ncol(obj), as.numeric(Sys.time()-t0)))
  obj
}

# ============================================================
# Run CellChat for one merged group
# ============================================================
run_cellchat_group <- function(group_name, samples) {
  cat(sprintf("\n\n#### === %s (samples: %s) ===\n",
              group_name, paste(samples, collapse = ", ")))
  
  # 1. Load + merge samples
  cat(sprintf("\n[1/6] Loading %d samples\n", length(samples)))
  obj_list <- list()
  for (s in samples) {
    o <- load_for_cellchat(s)
    if (!is.null(o)) obj_list[[s]] <- o
  }
  if (length(obj_list) == 0) return(NULL)
  
  if (length(obj_list) == 1) {
    merged <- obj_list[[1]]
  } else {
    # Merge: may need to handle different gene sets
    common_genes <- Reduce(intersect, lapply(obj_list, rownames))
    cat(sprintf("  Common genes across samples: %d\n", length(common_genes)))
    obj_list <- lapply(obj_list, function(o) o[common_genes, ])
    merged <- merge(obj_list[[1]], obj_list[-1],
                    add.cell.ids = names(obj_list))
    merged <- JoinLayers(merged)
    merged <- NormalizeData(merged, verbose = FALSE)
  }
  
  cat(sprintf("  Merged: %d spots, %d features\n", ncol(merged), nrow(merged)))
  
  # 2. Filter cell types with enough cells
  ct_counts <- table(merged$dominant_celltype)
  keep_ct <- names(ct_counts)[ct_counts >= MIN_CELLS_PER_CT]
  cat(sprintf("\n[2/6] Filtering: keep %d / %d cell types (n >= %d)\n",
              length(keep_ct), length(ct_counts), MIN_CELLS_PER_CT))
  keep <- merged$dominant_celltype %in% keep_ct
  merged <- merged[, keep]
  cat(sprintf("  After filter: %d spots\n", ncol(merged)))
  
  # 3. Build CellChat object
  cat(sprintf("\n[3/6] Building CellChat object\n"))
  data_input <- GetAssayData(merged, layer = "data")
  meta <- data.frame(
    labels = merged$dominant_celltype,
    row.names = colnames(merged)
  )
  
  cellchat <- createCellChat(object = data_input, meta = meta, group.by = "labels")
  
  # Set DB
  cellchat@DB <- CellChatDB.mouse
  
  # 4. Pre-processing
  cat(sprintf("\n[4/6] CellChat pre-processing\n"))
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  # Not all versions have projectData
  cellchat <- tryCatch(projectData(cellchat, PPI.mouse),
                       error = function(e) cellchat)
  
  # 5. Compute communication probability
  cat(sprintf("\n[5/6] Computing communication probability (slowest step, 20-60 min)\n"))
  t0 <- Sys.time()
  cellchat <- computeCommunProb(cellchat,
                                type = "triMean",
                                raw.use = TRUE,
                                population.size = FALSE)
  cat(sprintf("  computeCommunProb done: %.1f min\n",
              as.numeric(Sys.time()-t0, units = "mins")))
  
  cellchat <- filterCommunication(cellchat, min.cells = MIN_CELLS_PER_CT)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  # 6. Save + export
  cat(sprintf("\n[6/6] Saving outputs\n"))
  saveRDS(cellchat, file.path(OUT_DIR, sprintf("%s_cellchat.rds", group_name)))
  
  # LR pairs table
  df_net <- subsetCommunication(cellchat)
  write.csv(df_net, file.path(OUT_DIR, sprintf("%s_lr_pairs.csv", group_name)),
            row.names = FALSE)
  cat(sprintf("  -> %s_lr_pairs.csv (%d pairs)\n", group_name, nrow(df_net)))
  
  # Pathway table
  df_path <- subsetCommunication(cellchat, slot.name = "netP")
  write.csv(df_path, file.path(OUT_DIR, sprintf("%s_pathway_pairs.csv", group_name)),
            row.names = FALSE)
  cat(sprintf("  -> %s_pathway_pairs.csv (%d pairs)\n",
              group_name, nrow(df_path)))
  
  # Pathway list
  pathways <- cellchat@netP$pathways
  cat(sprintf("  Active pathways (%d): %s\n",
              length(pathways), paste(head(pathways, 20), collapse = ", ")))
  
  # Aggregate plot
  pdf(file.path(OUT_DIR, sprintf("%s_aggregate_network.pdf", group_name)),
      width = 10, height = 8)
  par(mfrow = c(1, 2), xpd = TRUE)
  tryCatch({
    netVisual_circle(cellchat@net$count, vertex.weight = as.numeric(table(cellchat@idents)),
                     weight.scale = TRUE, label.edge = FALSE,
                     title.name = "Number of interactions")
    netVisual_circle(cellchat@net$weight, vertex.weight = as.numeric(table(cellchat@idents)),
                     weight.scale = TRUE, label.edge = FALSE,
                     title.name = "Interaction weight/strength")
  }, error = function(e) cat(sprintf("  plot error: %s\n", conditionMessage(e))))
  dev.off()
  
  # Pathway network heatmap (top 5 pathways)
  top_pathways <- head(pathways, 5)
  for (pw in top_pathways) {
    tryCatch({
      pdf(file.path(OUT_DIR, sprintf("%s_pathway_%s.pdf", group_name, pw)),
          width = 7, height = 5)
      netVisual_heatmap(cellchat, signaling = pw, color.heatmap = "Reds")
      dev.off()
    }, error = function(e) cat(sprintf("  pathway %s error: %s\n",
                                       pw, conditionMessage(e))))
  }
  
  rm(merged, obj_list); gc(verbose = FALSE)
  
  list(
    group = group_name,
    n_pairs = nrow(df_net),
    n_pathways = length(pathways),
    pathways = pathways
  )
}

# ============================================================
# Run for liver_tumor + brain_tumor
# ============================================================
results <- list()
for (grp in names(GROUPS)) {
  tryCatch({
    results[[grp]] <- run_cellchat_group(grp, GROUPS[[grp]])
  }, error = function(e) {
    cat(sprintf("\n⚠ ERROR in %s: %s\n", grp, conditionMessage(e)))
  })
}

# ============================================================
# Summary.md
# ============================================================
cat("\n### Writing summary.md ###\n")

md <- c(
  "# T6 Spatial CellChat (Liver tumor + Brain tumor merged)",
  "",
  sprintf("**Generated:** %s", format(Sys.time(), "%Y-%m-%d %H:%M")),
  "",
  "## Groups",
  "",
  "- **liver_tumor**: Liver_3 + Liver_8 + gan56 (merged)",
  "- **brain_tumor**: Brain_3 + Brain_8 + nao56 (merged)",
  "",
  "## Results"
)
for (grp in names(results)) {
  if (is.null(results[[grp]])) next
  r <- results[[grp]]
  md <- c(md, "", sprintf("### %s", r$group),
    sprintf("- Total LR pairs: %d", r$n_pairs),
    sprintf("- Active pathways: %d", r$n_pathways),
    sprintf("- Top pathways: %s",
            paste(head(r$pathways, 10), collapse = ", ")))
}

md <- c(md, "",
  "## Output files per group",
  "",
  "- `<group>_cellchat.rds`: full CellChat object",
  "- `<group>_lr_pairs.csv`: all significant LR",
  "- `<group>_pathway_pairs.csv`: pathway-level",
  "- `<group>_aggregate_network.pdf`: cell-cell interaction overview",
  "- `<group>_pathway_<name>.pdf`: per-pathway heatmap (top 5)",
  "",
  "## Next steps",
  "- Inspect LR pairs for X1 TGFb axis (TGFB1 / TGFBR2 / SMAD3)",
  "- Inspect for X2 Thyroid (TTR / THRA, note: CellChatDB coverage may limit)",
  "- Inspect HSC_ECM ↔ myeloid (CSF1, GAS6-MERTK)"
)

writeLines(md, file.path(OUT_DIR, "summary.md"))

cat("\n============================================================\n")
cat("T6 COMPLETE\n")
cat(sprintf("Output: %s\n", OUT_DIR))
cat(sprintf("Files:  %d\n", length(list.files(OUT_DIR))))
cat("============================================================\n")
