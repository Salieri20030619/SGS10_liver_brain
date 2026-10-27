#!/usr/bin/env Rscript
# =============================================================================
# SGS10 D4 step1: Brain niche-level DEG (Wilcoxon + edgeR pseudobulk + null filter)
# -----------------------------------------------------------------------------
# 设计沿用 D3 step1 v2 的框架, 但 contrast 单位从 zone 切到 niche:
#   - 9 个 target niche: B-N01/02/03/04/05/06/09/11/13
#   - per niche × {tumor 3 sample vs baseline 2 sample}: 9 个主 contrast
#   - per niche × {Bla vs Con baseline 2 vs 2}: 9 个 null contrast
#   - 每个主 contrast 的 robust UP/DOWN = wilcoxon padj<0.05 + sample_consistency=1
#                                          - null contrast 中也显著的基因
#
# 输入:
#   D2_zone_output/spot_zone_stage_brain.csv  (293K spot × 18 col, 已含 niche+zone+sample)
#   spatial_results/rctd/enriched/{Bla_brain_1,Con_brain_1,Brain_3,Brain_8,nao56}_seekspace_rctd.rds
#
# 输出 (类似 D3 结构):
#   D4_brain_niche_DEG_output/
#     ├── wilcoxon/
#     │   ├── B-N01_tumor_vs_baseline.csv
#     │   ├── B-N01_NULL_Bla_vs_Con.csv
#     │   └── ... (×9 niche × 2 contrast = 18)
#     ├── pseudobulk_edgeR/  (only main contrasts × 9)
#     ├── combined/          (wilc × pb merge × 9)
#     ├── filtered_DEG/      (null-filtered robust UP/DOWN × 9)
#     └── logs/
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(edgeR)
})

option_list <- list(
  make_option("--base_dir", type="character",
              default="/mnt/BioDisk/wangmingyuan/SGS10_analysis"),
  make_option("--zsd_csv", type="character",
              default="/mnt/BioDisk/wangmingyuan/SGS10_analysis/D2_zone_output/spot_zone_stage_brain.csv"),
  make_option("--out_dir",  type="character",
              default="/mnt/BioDisk/wangmingyuan/SGS10_analysis/D4_brain_niche_DEG_output"),
  make_option("--target_niches", type="character",
              default="B-N01,B-N02,B-N03,B-N04,B-N05,B-N06,B-N09,B-N11,B-N13",
              help="comma-separated niche IDs"),
  make_option("--niche", type="character", default="all",
              help="run single niche or 'all'"),
  make_option("--min_pct", type="double", default=0.10),
  make_option("--logfc_thresh", type="double", default=0.25),
  make_option("--max_spots_per_group", type="integer", default=20000),
  make_option("--padj_cutoff", type="double", default=0.05),
  make_option("--skip_wilcoxon",   action="store_true", default=FALSE),
  make_option("--skip_pseudobulk", action="store_true", default=FALSE),
  make_option("--skip_null",       action="store_true", default=FALSE),
  make_option("--dry_run",         action="store_true", default=FALSE)
)
opt <- parse_args(OptionParser(option_list=option_list))

CFG <- list(
  base_dir = opt$base_dir,
  enriched_dir = file.path(opt$base_dir, "spatial_results/rctd/enriched"),
  zsd_csv = opt$zsd_csv,
  out_dir = opt$out_dir,
  target_niches = strsplit(opt$target_niches, ",")[[1]],
  niche_filter = opt$niche,
  min_pct = opt$min_pct,
  logfc_thresh = opt$logfc_thresh,
  max_spots_per_group = opt$max_spots_per_group,
  padj_cutoff = opt$padj_cutoff,
  skip_wilcoxon = opt$skip_wilcoxon,
  skip_pseudobulk = opt$skip_pseudobulk,
  skip_null = opt$skip_null,
  dry_run = opt$dry_run,
  brain_samples = c("Bla_brain_1","Con_brain_1","Brain_3","Brain_8","nao56"),
  brain_baseline = c("Bla_brain_1","Con_brain_1"),
  brain_tumor    = c("Brain_3","Brain_8","nao56")
)

# null-coalesce 兜底 (在使用前定义)
`%||%` <- function(a, b) if (is.null(a)) b else a

dir.create(CFG$out_dir, showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(CFG$out_dir, "wilcoxon"),         showWarnings=FALSE)
dir.create(file.path(CFG$out_dir, "pseudobulk_edgeR"), showWarnings=FALSE)
dir.create(file.path(CFG$out_dir, "combined"),         showWarnings=FALSE)
dir.create(file.path(CFG$out_dir, "filtered_DEG"),     showWarnings=FALSE)
dir.create(file.path(CFG$out_dir, "logs"),             showWarnings=FALSE)

LOG_FILE <- file.path(CFG$out_dir, "logs",
                      sprintf("D4_step1_%s.log",
                              format(Sys.time(), "%Y%m%d_%H%M%S")))
cat_log <- function(...) {
  msg <- paste0(..., collapse="")
  cat(msg)
  cat(msg, file=LOG_FILE, append=TRUE)
}
sec <- function(name) {
  cat_log("\n", strrep("=", 70), "\n  ", name, "\n",
          strrep("=", 70), "\n")
}

cat_log("=====================================================================\n")
cat_log("D4 Step 1: Brain niche-level DEG\n")
cat_log("Start: ", format(Sys.time()), "\n")
cat_log("Target niches: ", paste(CFG$target_niches, collapse=", "), "\n")
cat_log("Single niche filter: ", CFG$niche_filter, "\n")
cat_log("Seurat version: ", as.character(packageVersion("Seurat")), "\n")
cat_log("=====================================================================\n")

# ============ 1. Load metadata ============
sec("1. Load spot-level metadata (D2 output)")

zsd <- read_csv(CFG$zsd_csv, show_col_types=FALSE)
cat_log("zsd dim: ", paste(dim(zsd), collapse=" × "), "\n")
cat_log("zsd colnames: ", paste(colnames(zsd), collapse=", "), "\n")
cat_log("\nNiche × sample 分布:\n")
print(table(zsd$niche, zsd$sample))

# 仅保留 target niche 的 spot
zsd_target <- zsd %>% filter(niche %in% CFG$target_niches)
cat_log("\n保留 target niche 后: ", nrow(zsd_target), " spots\n")

# ============ 2. Define 18 contrasts (9 main + 9 null) ============
sec("2. Define contrasts (9 main + 9 null)")

contrasts_def <- list()
for (n in CFG$target_niches) {
  # Main: niche 内 tumor 3 sample vs baseline 2 sample
  contrasts_def[[paste0(n, "_tumor_vs_baseline")]] <- list(
    organ = "brain",
    niche = n,
    type = "main",
    group1_filter = local({nn <- n; function(z) z$niche == nn & z$sample %in% CFG$brain_tumor}),
    group2_filter = local({nn <- n; function(z) z$niche == nn & z$sample %in% CFG$brain_baseline}),
    group1_label = "tumor",
    group2_label = "baseline",
    can_pseudobulk = TRUE,
    can_consistency = TRUE,  # 3 vs 2 sample, 可以算 sample_consistency
    description = paste0("Brain ", n, " tumor (3 sample) vs baseline (2 sample)")
  )
  # Null: niche 内 Bla vs Con baseline batch noise
  contrasts_def[[paste0(n, "_NULL_Bla_vs_Con")]] <- list(
    organ = "brain",
    niche = n,
    type = "null",
    group1_filter = local({nn <- n; function(z) z$niche == nn & z$sample == "Bla_brain_1"}),
    group2_filter = local({nn <- n; function(z) z$niche == nn & z$sample == "Con_brain_1"}),
    group1_label = "Bla",
    group2_label = "Con",
    can_pseudobulk = FALSE,  # 1 vs 1 sample
    can_consistency = FALSE,
    description = paste0("Brain ", n, " Bla vs Con (null reference)")
  )
}

# 单 niche 过滤
if (CFG$niche_filter != "all") {
  keep <- grepl(paste0("^", CFG$niche_filter, "_"), names(contrasts_def))
  contrasts_def <- contrasts_def[keep]
}
if (CFG$skip_null) {
  contrasts_def <- contrasts_def[!grepl("_NULL_", names(contrasts_def))]
}

cat_log("\n共 ", length(contrasts_def), " 个 contrast 待跑:\n")
for (cn in names(contrasts_def)) {
  cat_log("  ", cn, " (", contrasts_def[[cn]]$type, ")\n")
}

if (CFG$dry_run) {
  cat_log("\n--dry_run 模式, 仅打印 contrast, 不实际跑\n")
  quit(status=0)
}

# ============ 3. Load Seurat objects + merge metadata ============
sec("3. Load Seurat objects")

load_brain_seurat <- function() {
  obj_list <- list()
  for (s in CFG$brain_samples) {
    fp <- file.path(CFG$enriched_dir, paste0(s, "_seekspace_rctd.rds"))
    cat_log("  reading ", basename(fp), " ... ")
    t0 <- Sys.time()
    obj <- readRDS(fp)
    cat_log(round(as.numeric(Sys.time()-t0, units="secs"), 1), "s, ",
            ncol(obj), " spots\n")
    # 添加 spot_key 列 (与 zsd 对应)
    obj$spot_key <- paste0(s, "_", colnames(obj))
    obj_list[[s]] <- obj
  }
  obj_list
}

# Merge 5 samples into 1 Seurat object
merge_brain <- function(obj_list) {
  cat_log("\nMerging 5 brain samples ...\n")
  merged <- merge(obj_list[[1]], y = obj_list[2:5],
                  add.cell.ids = names(obj_list)[1:5],
                  project = "brain_merged")
  # JoinLayers (Seurat v5 必须做; v4 不需要, try-catch)
  cat_log("\nJoinLayers (Seurat v5) ...\n")
  merged <- tryCatch(
    JoinLayers(merged),
    error = function(e) {
      cat_log("  JoinLayers 报错: ", conditionMessage(e),
              "  -- 假设 Seurat v4, 跳过\n")
      merged
    }
  )
  cat_log("  layers after JoinLayers: ",
          paste(Layers(merged), collapse=","), "\n")
  cat_log("  merged: ", ncol(merged), " spots\n")
  merged
}

# Map zsd metadata onto merged Seurat object
merge_metadata <- function(merged, zsd) {
  # merged barcode 格式: <sample>_<原 barcode>
  # 但 merge() 加的是 cell.ids prefix, 不是 spot_key (原 sample 加 "_")
  # 实际 colnames 是: "Bla_brain_1_TGAGACTG..."
  # 等于 zsd$spot_key
  
  spot_key <- colnames(merged)
  cat_log("  merged colnames example: ", head(spot_key, 2), "\n")
  cat_log("  zsd spot_key example: ", head(zsd$spot_key, 2), "\n")
  
  # zsd 索引化
  zsd_indexed <- zsd %>% select(spot_key, niche, zone, sample, group, stage_label, animal_id) %>%
                         distinct(spot_key, .keep_all = TRUE)
  
  zsd_ord <- zsd_indexed[match(spot_key, zsd_indexed$spot_key), ]
  
  n_match <- sum(!is.na(zsd_ord$niche))
  cat_log("  match: ", n_match, "/", length(spot_key), " (", round(100*n_match/length(spot_key),1), "%)\n")
  
  merged$niche       <- zsd_ord$niche
  merged$zsd_zone    <- zsd_ord$zone
  merged$stage_label <- zsd_ord$stage_label
  merged$animal_id   <- zsd_ord$animal_id
  
  cat_log("\n  niche distribution after merge:\n")
  print(table(merged$niche, useNA="ifany"))
  
  merged
}

cat_log("Loading 5 brain Seurat objects ...\n")
obj_list <- load_brain_seurat()
brain_merged <- merge_brain(obj_list)
brain_merged <- merge_metadata(brain_merged, zsd)
rm(obj_list); gc()

# ============ 4. Wilcoxon function ============
run_wilcoxon <- function(seurat_obj, ct, group1_idx, group2_idx) {
  cat_log("    [Wilcoxon] g1=", length(group1_idx), " g2=", length(group2_idx), " spots\n")
  
  if (length(group1_idx) < 3 || length(group2_idx) < 3) {
    cat_log("    [Wilcoxon] SKIP: too few spots\n")
    return(NULL)
  }
  
  # Subsample
  if (CFG$max_spots_per_group > 0) {
    if (length(group1_idx) > CFG$max_spots_per_group) {
      set.seed(42)
      group1_idx <- sample(group1_idx, CFG$max_spots_per_group)
    }
    if (length(group2_idx) > CFG$max_spots_per_group) {
      set.seed(42)
      group2_idx <- sample(group2_idx, CFG$max_spots_per_group)
    }
  }
  
  # Subset
  cells_use <- c(group1_idx, group2_idx)
  ident_vec <- c(rep(ct$group1_label, length(group1_idx)),
                 rep(ct$group2_label, length(group2_idx)))
  names(ident_vec) <- cells_use
  
  seurat_sub <- subset(seurat_obj, cells = cells_use)
  Idents(seurat_sub) <- ident_vec[colnames(seurat_sub)]
  
  cat_log("    [Wilcoxon] running FindMarkers ...\n")
  t0 <- Sys.time()
  markers <- FindMarkers(seurat_sub,
                         ident.1 = ct$group1_label,
                         ident.2 = ct$group2_label,
                         test.use = "wilcox",
                         min.pct = CFG$min_pct,
                         logfc.threshold = CFG$logfc_thresh,
                         only.pos = FALSE,
                         verbose = FALSE)
  dt <- as.numeric(Sys.time() - t0, units="secs")
  cat_log("    [Wilcoxon] FindMarkers done. ", round(dt/60, 1),
          " min, ", nrow(markers), " markers\n")
  
  markers$gene <- rownames(markers)
  rownames(markers) <- NULL
  
  # Sample consistency (only if can_consistency)
  if (ct$can_consistency) {
    # 算每个基因在 group1 各 sample 内 vs group2 各 sample 内的 log2FC sign 一致性
    g1_samples <- unique(seurat_obj$sample[group1_idx])
    g2_samples <- unique(seurat_obj$sample[group2_idx])
    n_pairs <- length(g1_samples) * length(g2_samples)
    
    if (length(g1_samples) >= 2 && length(g2_samples) >= 1 && nrow(markers) > 0) {
      cat_log("    [Wilcoxon] computing sample_consistency over ", n_pairs, " pairs ...\n")
      
      # 选择并行实现 (pbmcapply 没装就用 mclapply, 都没就串行)
      apply_fn <- if (requireNamespace("pbmcapply", quietly=TRUE)) {
        function(idx, fn) pbmcapply::pbmclapply(idx, fn, mc.cores=4)
      } else if (.Platform$OS.type == "unix") {
        function(idx, fn) parallel::mclapply(idx, fn, mc.cores=4)
      } else {
        function(idx, fn) lapply(idx, fn)
      }
      
      cnt <- apply_fn(seq_len(nrow(markers)), function(i) {
        g <- markers$gene[i]
        if (!(g %in% rownames(seurat_obj))) return(NA)
        same_dir <- 0; total <- 0
        for (s1 in g1_samples) {
          for (s2 in g2_samples) {
            v1 <- seurat_obj@assays$RNA$data[g, intersect(group1_idx,
                          colnames(seurat_obj)[seurat_obj$sample == s1])]
            v2 <- seurat_obj@assays$RNA$data[g, intersect(group2_idx,
                          colnames(seurat_obj)[seurat_obj$sample == s2])]
            if (length(v1) < 3 || length(v2) < 3) next
            d <- mean(v1) - mean(v2)
            same_dir <- same_dir + (sign(d) == sign(markers$avg_log2FC[i]))
            total <- total + 1
          }
        }
        if (total == 0) NA else same_dir / total
      })
      markers$sample_consistency <- unlist(cnt)
    } else {
      markers$sample_consistency <- NA
    }
  } else {
    markers$sample_consistency <- NA
  }
  
  markers
}

# ============ 5. Pseudobulk edgeR ============
run_pseudobulk_edgeR <- function(seurat_obj, ct, group1_idx, group2_idx) {
  cat_log("    [Pseudobulk] g1=", length(group1_idx), " g2=", length(group2_idx), " spots\n")
  
  if (!ct$can_pseudobulk) {
    cat_log("    [Pseudobulk] SKIP (per ct$can_pseudobulk=FALSE)\n")
    return(NULL)
  }
  
  cells_use <- c(group1_idx, group2_idx)
  samples <- seurat_obj$sample[cells_use]
  groups  <- c(rep(ct$group1_label, length(group1_idx)),
               rep(ct$group2_label, length(group2_idx)))
  
  # aggregate counts by sample
  counts <- GetAssayData(seurat_obj, assay="RNA", layer="counts")[, cells_use]
  
  agg <- list()
  meta <- data.frame(sample=character(), group=character())
  for (s in unique(samples)) {
    idx <- which(samples == s)
    g_s <- unique(groups[idx])
    if (length(g_s) > 1) next  # 该 sample 跨组, 不是有效 pseudobulk
    if (length(idx) < 5) next  # spots 太少
    agg[[s]] <- Matrix::rowSums(counts[, idx, drop=FALSE])
    meta <- rbind(meta, data.frame(sample=s, group=g_s))
  }
  
  if (length(agg) < 4) {
    cat_log("    [Pseudobulk] SKIP: only ", length(agg), " valid samples\n")
    return(NULL)
  }
  
  agg_mat <- do.call(cbind, agg)
  colnames(agg_mat) <- names(agg)
  meta <- meta[match(colnames(agg_mat), meta$sample), ]
  cat_log("    [Pseudobulk] N samples per group: ",
          paste(table(meta$group), collapse=" / "), "\n")
  
  if (length(unique(meta$group)) < 2) {
    cat_log("    [Pseudobulk] SKIP: only 1 group has samples\n")
    return(NULL)
  }
  
  # edgeR
  y <- DGEList(counts=agg_mat, group=meta$group)
  keep <- filterByExpr(y, group=meta$group)
  y <- y[keep, , keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  
  design <- model.matrix(~ 0 + meta$group)
  colnames(design) <- levels(factor(meta$group))
  
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  
  contrast <- makeContrasts(
    contrasts = paste0(ct$group1_label, " - ", ct$group2_label),
    levels = design
  )
  qlf <- glmQLFTest(fit, contrast=contrast)
  
  tt <- topTags(qlf, n=Inf, sort.by="PValue")$table
  tt$gene <- rownames(tt)
  rownames(tt) <- NULL
  cat_log("    [Pseudobulk] done. ", nrow(tt), " genes after filterByExpr\n")
  
  tt
}

# ============ 6. Per-contrast main loop ============
sec("6. Run contrasts")

# 给 brain_merged 加 cell index
all_cells <- colnames(brain_merged)
zsd_brain <- zsd %>% filter(spot_key %in% all_cells)
cat_log("zsd_brain rows: ", nrow(zsd_brain), "\n")

results_summary <- list()

for (cn in names(contrasts_def)) {
  ct <- contrasts_def[[cn]]
  cat_log("\n", strrep("-", 70), "\n",
          "Contrast: ", cn, " [", ct$type, "]\n",
          ct$description, "\n",
          strrep("-", 70), "\n")
  
  group1_keys <- zsd_brain$spot_key[ct$group1_filter(zsd_brain)]
  group2_keys <- zsd_brain$spot_key[ct$group2_filter(zsd_brain)]
  group1_idx <- intersect(group1_keys, all_cells)
  group2_idx <- intersect(group2_keys, all_cells)
  
  cat_log("  matched g1=", length(group1_idx), " g2=", length(group2_idx), "\n")
  
  res_summary <- list(contrast=cn, type=ct$type,
                      n_g1=length(group1_idx), n_g2=length(group2_idx),
                      wilc_n=NA, pb_n=NA)
  
  if (length(group1_idx) < 3 || length(group2_idx) < 3) {
    cat_log("  SKIP: too few spots\n")
    res_summary$status <- "skipped_too_few"
    results_summary[[cn]] <- res_summary
    next
  }
  
  # Wilcoxon
  if (!CFG$skip_wilcoxon) {
    wilc <- tryCatch(
      run_wilcoxon(brain_merged, ct, group1_idx, group2_idx),
      error = function(e) {
        cat_log("    [Wilcoxon] ERROR: ", conditionMessage(e), "\n"); NULL
      })
    if (!is.null(wilc) && nrow(wilc) > 0) {
      out <- file.path(CFG$out_dir, "wilcoxon", paste0(cn, ".csv"))
      write_csv(wilc, out)
      cat_log("  [Wilcoxon] saved: ", basename(out), "\n")
      res_summary$wilc_n <- nrow(wilc)
    }
  }
  
  # Pseudobulk (main only)
  if (!CFG$skip_pseudobulk && ct$can_pseudobulk) {
    pb <- tryCatch(
      run_pseudobulk_edgeR(brain_merged, ct, group1_idx, group2_idx),
      error = function(e) {
        cat_log("    [Pseudobulk] ERROR: ", conditionMessage(e), "\n"); NULL
      })
    if (!is.null(pb) && nrow(pb) > 0) {
      out <- file.path(CFG$out_dir, "pseudobulk_edgeR", paste0(cn, ".csv"))
      write_csv(pb, out)
      cat_log("  [Pseudobulk] saved: ", basename(out), "\n")
      res_summary$pb_n <- nrow(pb)
    }
  }
  
  res_summary$status <- "ok"
  results_summary[[cn]] <- res_summary
}

# ============ 7. Combined: wilc × pb merge ============
sec("7. Combined tables (wilcoxon × pseudobulk_edgeR)")

for (n in CFG$target_niches) {
  cn <- paste0(n, "_tumor_vs_baseline")
  wf <- file.path(CFG$out_dir, "wilcoxon", paste0(cn, ".csv"))
  pf <- file.path(CFG$out_dir, "pseudobulk_edgeR", paste0(cn, ".csv"))
  if (!file.exists(wf) || !file.exists(pf)) {
    cat_log("  ", n, ": skip combined (missing input)\n")
    next
  }
  wilc <- read_csv(wf, show_col_types=FALSE)
  pb <- read_csv(pf, show_col_types=FALSE)
  
  comb <- wilc %>%
    select(gene, wilc_log2FC=avg_log2FC, wilc_p_adj=p_val_adj,
           wilc_pct1=pct.1, wilc_pct2=pct.2, sample_consistency) %>%
    full_join(
      pb %>% select(gene, pb_log2FC=logFC, pb_p_adj=FDR, pb_logCPM=logCPM),
      by="gene"
    ) %>%
    mutate(
      both_up = !is.na(wilc_p_adj) & wilc_p_adj < CFG$padj_cutoff & wilc_log2FC > CFG$logfc_thresh &
                !is.na(pb_p_adj) & pb_p_adj < CFG$padj_cutoff & pb_log2FC > 0,
      both_down = !is.na(wilc_p_adj) & wilc_p_adj < CFG$padj_cutoff & wilc_log2FC < -CFG$logfc_thresh &
                  !is.na(pb_p_adj) & pb_p_adj < CFG$padj_cutoff & pb_log2FC < 0,
      direction = case_when(both_up ~ "UP", both_down ~ "DOWN", TRUE ~ "discordant_or_ns")
    )
  
  out <- file.path(CFG$out_dir, "combined", paste0(cn, "_combined.csv"))
  write_csv(comb, out)
  cat_log("  ", n, " combined saved (", nrow(comb), " genes)\n")
}

# ============ 8. Filtered DEG (null filter) ============
sec("8. Filtered DEG (apply null filter per niche)")

for (n in CFG$target_niches) {
  main_cn <- paste0(n, "_tumor_vs_baseline")
  null_cn <- paste0(n, "_NULL_Bla_vs_Con")
  
  main_wf <- file.path(CFG$out_dir, "wilcoxon", paste0(main_cn, ".csv"))
  null_wf <- file.path(CFG$out_dir, "wilcoxon", paste0(null_cn, ".csv"))
  
  if (!file.exists(main_wf)) {
    cat_log("  ", n, ": skip filtered (no main wilcoxon)\n"); next
  }
  
  main_w <- read_csv(main_wf, show_col_types=FALSE)
  
  # Null filter: 找出 null 中也显著的基因, 从 main 中减去
  null_sig_genes <- character(0)
  if (file.exists(null_wf) && !CFG$skip_null) {
    null_w <- read_csv(null_wf, show_col_types=FALSE)
    null_sig <- null_w %>% filter(p_val_adj < CFG$padj_cutoff & abs(avg_log2FC) > CFG$logfc_thresh)
    null_sig_genes <- null_sig$gene
    cat_log("  ", n, " null-significant genes: ", length(null_sig_genes), "\n")
  }
  
  # Robust UP: padj<0.05, log2FC>thresh, sample_consistency=1 (or NA), NOT in null
  up_filt <- main_w %>%
    filter(p_val_adj < CFG$padj_cutoff,
           avg_log2FC > CFG$logfc_thresh,
           is.na(sample_consistency) | sample_consistency == 1,
           !(gene %in% null_sig_genes)) %>%
    arrange(desc(avg_log2FC))
  
  down_filt <- main_w %>%
    filter(p_val_adj < CFG$padj_cutoff,
           avg_log2FC < -CFG$logfc_thresh,
           is.na(sample_consistency) | sample_consistency == 1,
           !(gene %in% null_sig_genes)) %>%
    arrange(avg_log2FC)
  
  out_up <- file.path(CFG$out_dir, "filtered_DEG", paste0(n, "_UP_filtered.csv"))
  out_dn <- file.path(CFG$out_dir, "filtered_DEG", paste0(n, "_DOWN_filtered.csv"))
  write_csv(up_filt, out_up)
  write_csv(down_filt, out_dn)
  cat_log("  ", n, " UP=", nrow(up_filt), " DOWN=", nrow(down_filt),
          " (after null filter)\n")
}

# ============ 9. Summary ============
sec("9. Summary")

summary_df <- do.call(rbind, lapply(names(results_summary), function(cn) {
  rs <- results_summary[[cn]]
  data.frame(
    contrast = cn,
    type = rs$type %||% NA,
    n_g1 = rs$n_g1 %||% NA,
    n_g2 = rs$n_g2 %||% NA,
    wilc_n = rs$wilc_n %||% NA,
    pb_n = rs$pb_n %||% NA,
    status = rs$status %||% NA
  )
}))
write_csv(summary_df, file.path(CFG$out_dir, "task_summary.csv"))
cat_log("\nFull summary:\n")
print(summary_df)

# Also: filtered_DEG summary
filt_summary <- do.call(rbind, lapply(CFG$target_niches, function(n) {
  up_fp <- file.path(CFG$out_dir, "filtered_DEG", paste0(n, "_UP_filtered.csv"))
  dn_fp <- file.path(CFG$out_dir, "filtered_DEG", paste0(n, "_DOWN_filtered.csv"))
  data.frame(
    niche = n,
    n_robust_UP   = if (file.exists(up_fp)) nrow(read_csv(up_fp, show_col_types=FALSE)) else NA,
    n_robust_DOWN = if (file.exists(dn_fp)) nrow(read_csv(dn_fp, show_col_types=FALSE)) else NA
  )
}))
write_csv(filt_summary, file.path(CFG$out_dir, "filtered_DEG_summary.csv"))
cat_log("\nFiltered (robust) DEG counts per niche:\n")
print(filt_summary)

cat_log("\n=====================================================================\n")
cat_log("D4 step1 完成: ", format(Sys.time()), "\n")
cat_log("Output: ", CFG$out_dir, "\n")
cat_log("=====================================================================\n")
