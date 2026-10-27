#!/usr/bin/env Rscript
# =============================================================================
# SGS10 D3 step1 v2: Per-cell DEG (Wilcoxon + edgeR pseudobulk)
# -----------------------------------------------------------------------------
# v2 changes:
#   - Removed silent outer tryCatch on Wilcoxon (let errors surface)
#   - Improved subset ordering: subsample BEFORE creating sub-Seurat
#   - More verbose diagnostic prints at each Wilcoxon step
#   - Pseudobulk and Wilcoxon now fully independent (no shared state)
#   - Verbose tryCatch message printing (not just "ERROR")
#
# Output structure same as v1:
#   wilcoxon/<contrast>.csv
#   pseudobulk_edgeR/<contrast>.csv
#   combined/<contrast>_combined.csv
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
  make_option("--d2_dir",   type="character",
              default="/mnt/BioDisk/wangmingyuan/SGS10_analysis/D2_zone_output"),
  make_option("--out_dir",  type="character",
              default="/mnt/BioDisk/wangmingyuan/SGS10_analysis/D3_DEG_output"),
  make_option("--contrast", type="character", default="all",
              help="A1,A2,A3,B3,C4 or 'all'"),
  make_option("--min_pct", type="double", default=0.10),
  make_option("--logfc_thresh", type="double", default=0.25),
  make_option("--max_spots_per_group", type="integer", default=20000,
              help="cap spots per group for Wilcoxon (speed). Set 0 to disable"),
  make_option("--skip_wilcoxon",   action="store_true", default=FALSE),
  make_option("--skip_pseudobulk", action="store_true", default=FALSE),
  make_option("--dry_run",         action="store_true", default=FALSE)
)
opt <- parse_args(OptionParser(option_list=option_list))

CFG <- list(
  base_dir = opt$base_dir,
  enriched_dir = file.path(opt$base_dir, "spatial_results/rctd/enriched"),
  d2_dir = opt$d2_dir,
  out_dir = opt$out_dir,
  contrast = opt$contrast,
  min_pct = opt$min_pct,
  logfc_thresh = opt$logfc_thresh,
  max_spots_per_group = opt$max_spots_per_group,
  skip_wilcoxon = opt$skip_wilcoxon,
  skip_pseudobulk = opt$skip_pseudobulk,
  dry_run = opt$dry_run,
  liver_samples = c("Con_liver_1","Bla_liver_1","Liver_3","Liver_8","gan56"),
  brain_samples = c("Con_brain_1","Bla_brain_1","Brain_3","Brain_8","nao56"),
  liver_baseline = c("Con_liver_1","Bla_liver_1"),
  brain_baseline = c("Con_brain_1","Bla_brain_1")
)

dir.create(CFG$out_dir, showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(CFG$out_dir, "wilcoxon"),         showWarnings=FALSE)
dir.create(file.path(CFG$out_dir, "pseudobulk_edgeR"), showWarnings=FALSE)
dir.create(file.path(CFG$out_dir, "combined"),         showWarnings=FALSE)

cat_log <- function(...) cat(..., sep="")
sec <- function(name) {
  cat_log("\n", strrep("-", 70), "\n  ", name, "\n",
          strrep("-", 70), "\n")
}

cat_log("=====================================================================\n")
cat_log("D3 Step 1 v2: Per-cell DEG (Wilcoxon + edgeR pseudobulk)\n")
cat_log("Start: ", format(Sys.time()), "\n")
cat_log("Contrast: ", CFG$contrast, "\n")
cat_log("Seurat version: ", as.character(packageVersion("Seurat")), "\n")
cat_log("=====================================================================\n")

# ============ Define contrasts ============
contrasts_def <- list(
  A1 = list(
    organ = "liver",
    group1_filter = function(zsd) zsd$zone == "tumor_core",
    group2_filter = function(zsd) zsd$zone == "baseline_host",
    group1_label = "tumor_core",
    group2_label = "baseline_host",
    can_pseudobulk = TRUE,
    description = "Liver tumor_core (all 3 stages) vs Liver baseline_host"
  ),
  A2 = list(
    organ = "liver",
    group1_filter = function(zsd) zsd$zone == "invaded_host",
    group2_filter = function(zsd) zsd$zone == "baseline_host",
    group1_label = "invaded_host",
    group2_label = "baseline_host",
    can_pseudobulk = TRUE,
    description = "Liver invaded_host (all 3 stages) vs Liver baseline_host"
  ),
  A3 = list(
    organ = "liver",
    group1_filter = function(zsd) zsd$zone == "distant_host",
    group2_filter = function(zsd) zsd$zone == "baseline_host",
    group1_label = "distant_host",
    group2_label = "baseline_host",
    can_pseudobulk = TRUE,
    description = "Liver distant_host (all 3 stages) vs Liver baseline_host"
  ),
  B3 = list(
    organ = "liver",
    group1_filter = function(zsd) zsd$zone == "distant_host" & zsd$stage_label == "late",
    group2_filter = function(zsd) zsd$zone == "distant_host" & zsd$stage_label == "early",
    group1_label = "distant_late",
    group2_label = "distant_early",
    can_pseudobulk = FALSE,
    description = "Liver distant_host late vs early (pseudobulk skipped, N=1 each)"
  ),
  C4 = list(
    organ = "brain",
    group1_filter = function(zsd) zsd$zone == "distant_host",
    group2_filter = function(zsd) zsd$zone == "baseline_host",
    group1_label = "brain_distant",
    group2_label = "brain_baseline",
    can_pseudobulk = TRUE,
    description = "Brain distant_host (all 3 stages) vs Brain baseline_host"
  )
)

if (CFG$contrast == "all") {
  contrasts_to_run <- names(contrasts_def)
} else {
  contrasts_to_run <- strsplit(CFG$contrast, ",")[[1]]
  unknown <- setdiff(contrasts_to_run, names(contrasts_def))
  if (length(unknown) > 0) stop("Unknown contrast: ", paste(unknown, collapse=","))
}
cat_log("\nContrasts to run: ", paste(contrasts_to_run, collapse=", "), "\n")

# ============ Helper: Load mouse-only counts ============

load_organ_seurat <- function(organ) {
  samples <- if (organ == "liver") CFG$liver_samples else CFG$brain_samples
  cat_log("  Loading ", length(samples), " ", organ, " samples (mouse-only)...\n")

  so_list <- list()
  for (s in samples) {
    fp <- file.path(CFG$enriched_dir, paste0(s, "_seekspace_rctd.rds"))
    if (!file.exists(fp)) {
      warning("Missing: ", fp); next
    }
    cat_log("    ", s, " ...\n")
    so <- readRDS(fp)
    counts <- GetAssayData(so, layer = "counts")
    feats <- rownames(counts)

    if (any(grepl("^mm10---", feats))) {
      keep_idx <- grep("^mm10---", feats)
      counts <- counts[keep_idx, , drop=FALSE]
      rownames(counts) <- sub("^mm10---", "", rownames(counts))
    }

    new_so <- CreateSeuratObject(counts = counts)
    new_so$sample_id <- s
    so_list[[s]] <- new_so
    rm(so, counts); gc(verbose=FALSE)
    cat_log("      ", s, ": ", ncol(new_so), " spots, ",
            nrow(new_so), " mouse genes\n")
  }

  cat_log("  Merging ", length(so_list), " samples...\n")
  if (length(so_list) == 1) {
    merged <- so_list[[1]]
  } else {
    merged <- merge(so_list[[1]], y = so_list[-1],
                    add.cell.ids = names(so_list))
  }
  rm(so_list); gc(verbose=FALSE)
  cat_log("  merged: ", ncol(merged), " spots, ", nrow(merged), " genes\n")

  cat_log("  JoinLayers + NormalizeData ...\n")
  merged <- JoinLayers(merged)
  merged <- NormalizeData(merged, verbose = FALSE)

  zsd_path <- file.path(CFG$d2_dir, paste0("spot_zone_stage_", organ, ".rds"))
  zsd <- readRDS(zsd_path)

  matched <- intersect(colnames(merged), zsd$spot_key)
  cat_log("  matched ", length(matched), " / ", ncol(merged),
          " spots with D2 labels\n")
  if (length(matched) < ncol(merged) * 0.95) {
    warning(organ, ": only ", length(matched), " / ", ncol(merged), " matched")
  }

  merged <- merged[, matched]   # use [,cells] instead of subset()
  zsd_ord <- zsd[match(colnames(merged), zsd$spot_key), ]
  merged$zone        <- zsd_ord$zone
  merged$stage_label <- zsd_ord$stage_label
  merged$animal_id   <- zsd_ord$animal_id
  merged$human_pct   <- zsd_ord$human_pct

  cat_log("  Final: ", ncol(merged), " spots, ", nrow(merged), " genes\n")
  cat_log("  zone distribution:\n")
  print(table(merged$zone, useNA="ifany"))

  merged
}

# ============ Wilcoxon helper (v2: more robust) ============

run_wilcoxon <- function(seurat_obj, ct, group1_idx, group2_idx) {
  cat_log("    [Wilcoxon] g1=", length(group1_idx), " g2=", length(group2_idx), " spots\n")

  # Subsample FIRST (before subset)
  if (CFG$max_spots_per_group > 0) {
    if (length(group1_idx) > CFG$max_spots_per_group) {
      set.seed(42)
      group1_idx <- sample(group1_idx, CFG$max_spots_per_group)
      cat_log("    [Wilcoxon] g1 downsampled to ", length(group1_idx), "\n")
    }
    if (length(group2_idx) > CFG$max_spots_per_group) {
      set.seed(43)
      group2_idx <- sample(group2_idx, CFG$max_spots_per_group)
      cat_log("    [Wilcoxon] g2 downsampled to ", length(group2_idx), "\n")
    }
  }

  # Get cell names for selected indices
  all_cells <- colnames(seurat_obj)
  g1_cells <- all_cells[group1_idx]
  g2_cells <- all_cells[group2_idx]
  use_cells <- c(g1_cells, g2_cells)

  cat_log("    [Wilcoxon] subsetting to ", length(use_cells), " cells via [,cells]\n")
  seurat_sub <- seurat_obj[, use_cells]   # safer than subset()
  cat_log("    [Wilcoxon] subset OK: ", ncol(seurat_sub), " cells\n")

  # Build ident vector aligned with subset cells
  cat_log("    [Wilcoxon] setting Idents...\n")
  ident_vec <- ifelse(colnames(seurat_sub) %in% g1_cells,
                       ct$group1_label, ct$group2_label)
  names(ident_vec) <- colnames(seurat_sub)
  Idents(seurat_sub) <- ident_vec
  cat_log("    [Wilcoxon] Idents OK. counts: ")
  print(table(Idents(seurat_sub)))

  # Run FindMarkers
  cat_log("    [Wilcoxon] running FindMarkers...\n")
  t0 <- Sys.time()
  markers <- FindMarkers(seurat_sub,
                         ident.1 = ct$group1_label,
                         ident.2 = ct$group2_label,
                         test.use = "wilcox",
                         min.pct = CFG$min_pct,
                         logfc.threshold = CFG$logfc_thresh,
                         only.pos = FALSE,
                         verbose = FALSE)
  dt <- round(as.numeric(Sys.time()-t0, units="secs"), 1)
  cat_log("    [Wilcoxon] FindMarkers done. ", round(dt/60, 1),
          " min, ", nrow(markers), " markers returned\n")

  if (nrow(markers) == 0) {
    cat_log("    [Wilcoxon] WARNING: no markers found, returning empty\n")
    markers$gene <- character(0)
    markers$sample_consistency <- numeric(0)
    return(markers)
  }

  markers$gene <- rownames(markers)

  # Per-sample consistency (defensive)
  cat_log("    [Wilcoxon] computing per-sample consistency...\n")
  markers$sample_consistency <- NA_real_

  consistency_result <- tryCatch({
    use_genes <- rownames(markers)
    data_mat <- GetAssayData(seurat_sub, layer = "data")
    use_genes <- intersect(use_genes, rownames(data_mat))
    if (length(use_genes) == 0) stop("no genes in data layer")

    sub_data <- data_mat[use_genes, , drop = FALSE]
    sub_sample <- as.character(seurat_sub$sample_id)
    sub_group  <- as.character(Idents(seurat_sub))

    samples_all <- unique(sub_sample)
    cat_log("      [consistency] ", length(samples_all), " samples x ",
            length(use_genes), " genes\n")

    sample_means <- matrix(NA_real_,
                           nrow = length(use_genes),
                           ncol = length(samples_all),
                           dimnames = list(use_genes, samples_all))
    for (s in samples_all) {
      cells_s <- which(sub_sample == s)
      if (length(cells_s) > 0) {
        sample_means[, s] <- as.numeric(Matrix::rowMeans(sub_data[, cells_s, drop = FALSE]))
      }
    }

    samples_g1 <- unique(sub_sample[sub_group == ct$group1_label])
    samples_g2 <- unique(sub_sample[sub_group == ct$group2_label])
    cat_log("      [consistency] g1 samples: ",
            paste(samples_g1, collapse=","), "\n")
    cat_log("      [consistency] g2 samples: ",
            paste(samples_g2, collapse=","), "\n")

    if (length(samples_g1) == 0 || length(samples_g2) == 0) {
      stop("zero samples in one group")
    }

    # mean of g2 samples per gene
    if (length(samples_g2) == 1) {
      g2_avg <- sample_means[, samples_g2]
    } else {
      g2_avg <- rowMeans(sample_means[, samples_g2, drop = FALSE], na.rm = TRUE)
    }

    g1_means_mat <- sample_means[, samples_g1, drop = FALSE]
    cons_up <- rowMeans(g1_means_mat > g2_avg, na.rm = TRUE)
    cons_down <- 1 - cons_up

    cons_final <- rep(NA_real_, nrow(markers))
    matched_idx <- match(use_genes, rownames(markers))
    for (i in seq_along(use_genes)) {
      mi <- matched_idx[i]
      if (is.na(mi)) next
      if (markers$avg_log2FC[mi] > 0) {
        cons_final[mi] <- cons_up[i]
      } else {
        cons_final[mi] <- cons_down[i]
      }
    }
    cat_log("      [consistency] ", sum(!is.na(cons_final)), " / ",
            nrow(markers), " genes scored\n")
    cons_final
  }, error = function(e) {
    cat_log("      [consistency] FAILED: ", conditionMessage(e),
            " (proceeding with NA)\n")
    rep(NA_real_, nrow(markers))
  })

  markers$sample_consistency <- consistency_result
  markers <- markers[order(-markers$avg_log2FC), ]
  markers
}

# ============ Pseudobulk edgeR helper ============

run_pseudobulk_edgeR <- function(seurat_obj, ct, group1_idx, group2_idx) {
  cat_log("    [pseudobulk] aggregating by sample...\n")

  ident_vec <- rep(NA_character_, ncol(seurat_obj))
  ident_vec[group1_idx] <- ct$group1_label
  ident_vec[group2_idx] <- ct$group2_label

  counts <- GetAssayData(seurat_obj, layer = "counts")
  meta <- data.frame(
    sample_id = seurat_obj$sample_id,
    group = ident_vec,
    stringsAsFactors = FALSE
  )

  meta$grp_sample <- paste0(meta$sample_id, "__", meta$group)
  keep_spots <- !is.na(meta$group)
  meta_sub <- meta[keep_spots, ]
  counts_sub <- counts[, keep_spots, drop=FALSE]

  pb_groups <- unique(meta_sub$grp_sample)
  cat_log("    [pseudobulk] units: ", length(pb_groups), "\n")
  pb_counts <- sapply(pb_groups, function(g) {
    cells_g <- which(meta_sub$grp_sample == g)
    Matrix::rowSums(counts_sub[, cells_g, drop=FALSE])
  })

  parts <- strsplit(pb_groups, "__")
  pb_meta <- data.frame(
    pb_id = pb_groups,
    sample_id = sapply(parts, `[`, 1),
    group = sapply(parts, `[`, 2),
    n_spots = as.integer(table(meta_sub$grp_sample)[pb_groups]),
    stringsAsFactors = FALSE
  )
  pb_meta$group <- factor(pb_meta$group,
                          levels = c(ct$group2_label, ct$group1_label))
  cat_log("    [pseudobulk] table:\n")
  print(pb_meta)

  if (nrow(pb_meta) < 4 ||
      length(unique(pb_meta$sample_id[pb_meta$group == ct$group1_label])) < 2 ||
      length(unique(pb_meta$sample_id[pb_meta$group == ct$group2_label])) < 2) {
    cat_log("    [pseudobulk] SKIP: insufficient samples (need >= 2 per group)\n")
    return(NULL)
  }

  cat_log("    [pseudobulk] running edgeR glmQLFit...\n")
  y <- DGEList(counts = pb_counts, samples = pb_meta)

  keep_genes <- filterByExpr(y, group = pb_meta$group)
  cat_log("    [pseudobulk] filterByExpr keeps ", sum(keep_genes), " / ",
          nrow(y), " genes\n")
  y <- y[keep_genes, , keep.lib.sizes = FALSE]

  y <- calcNormFactors(y)
  design <- model.matrix(~ group, data = pb_meta)

  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design, robust = TRUE)
  qlf <- glmQLFTest(fit, coef = 2)

  res <- topTags(qlf, n = Inf, sort.by = "PValue")$table
  res$gene <- rownames(res)
  cat_log("    [pseudobulk] edgeR done: ", nrow(res), " genes tested, ",
          sum(res$FDR < 0.05), " FDR<0.05\n")
  res
}

# ============ Main ============

organ_needed <- unique(sapply(contrasts_to_run,
                              function(ct_id) contrasts_def[[ct_id]]$organ))
cat_log("\nWill load organs: ", paste(organ_needed, collapse=", "), "\n")

if (CFG$dry_run) {
  cat_log("\n[DRY RUN] OK\n")
  quit(save="no", status=0)
}

organ_objs <- list()
for (org in organ_needed) {
  sec(paste("Load", org))
  organ_objs[[org]] <- load_organ_seurat(org)
  gc(verbose=FALSE)
}

# ============ Per contrast ============
for (ct_id in contrasts_to_run) {
  ct <- contrasts_def[[ct_id]]
  sec(paste0("Contrast ", ct_id, ": ", ct$description))

  so <- organ_objs[[ct$organ]]
  group1_idx <- which(ct$group1_filter(so@meta.data))
  group2_idx <- which(ct$group2_filter(so@meta.data))

  cat_log("  group1 (", ct$group1_label, "): ", length(group1_idx), " spots\n")
  cat_log("  group2 (", ct$group2_label, "): ", length(group2_idx), " spots\n")

  # Print sample composition (sanity check)
  if (length(group1_idx) > 0) {
    cat_log("  group1 sample distribution: ")
    print(table(so$sample_id[group1_idx]))
  }
  if (length(group2_idx) > 0) {
    cat_log("  group2 sample distribution: ")
    print(table(so$sample_id[group2_idx]))
  }

  if (length(group1_idx) < 50 || length(group2_idx) < 50) {
    cat_log("  SKIP: too few spots\n")
    next
  }

  # ----- Wilcoxon -----
  if (!CFG$skip_wilcoxon) {
    cat_log("\n  ===== Wilcoxon =====\n")
    wilc_result <- tryCatch({
      run_wilcoxon(so, ct, group1_idx, group2_idx)
    }, error = function(e) {
      cat_log("  *** WILCOXON FAILED ***\n")
      cat_log("  message: ", conditionMessage(e), "\n")
      cat_log("  traceback:\n")
      print(sys.calls())
      NULL
    })

    if (!is.null(wilc_result) && nrow(wilc_result) > 0) {
      out_path <- file.path(CFG$out_dir, "wilcoxon",
                           paste0(ct_id, "_", ct$group1_label, "_vs_",
                                  ct$group2_label, ".csv"))
      write_csv(wilc_result, out_path)
      cat_log("  saved Wilcoxon: ", out_path,
              " (", nrow(wilc_result), " rows)\n")
    } else {
      cat_log("  Wilcoxon produced 0 rows or failed, no file saved\n")
    }
  }

  # ----- Pseudobulk -----
  if (!CFG$skip_pseudobulk && ct$can_pseudobulk) {
    cat_log("\n  ===== Pseudobulk edgeR =====\n")
    pb_result <- tryCatch({
      run_pseudobulk_edgeR(so, ct, group1_idx, group2_idx)
    }, error = function(e) {
      cat_log("  *** PSEUDOBULK FAILED ***\n")
      cat_log("  message: ", conditionMessage(e), "\n")
      NULL
    })

    if (!is.null(pb_result)) {
      out_path <- file.path(CFG$out_dir, "pseudobulk_edgeR",
                           paste0(ct_id, "_", ct$group1_label, "_vs_",
                                  ct$group2_label, ".csv"))
      write_csv(pb_result, out_path)
      cat_log("  saved pseudobulk: ", out_path, "\n")
    }
  } else if (!ct$can_pseudobulk) {
    cat_log("\n  pseudobulk SKIPPED for this contrast (N too small)\n")
  }

  gc(verbose=FALSE)
}

# ============ Combine ============
sec("Combine Wilcoxon + pseudobulk")

for (ct_id in contrasts_to_run) {
  ct <- contrasts_def[[ct_id]]
  if (!ct$can_pseudobulk) {
    cat_log("  ", ct_id, ": no pseudobulk (skip combine)\n")
    next
  }

  wilc_path <- file.path(CFG$out_dir, "wilcoxon",
                        paste0(ct_id, "_", ct$group1_label, "_vs_",
                               ct$group2_label, ".csv"))
  pb_path <- file.path(CFG$out_dir, "pseudobulk_edgeR",
                      paste0(ct_id, "_", ct$group1_label, "_vs_",
                             ct$group2_label, ".csv"))
  if (!file.exists(wilc_path) || !file.exists(pb_path)) {
    cat_log("  ", ct_id, ": missing wilc (",
            file.exists(wilc_path), ") or pb (",
            file.exists(pb_path), "), skip combine\n")
    next
  }

  wilc <- read_csv(wilc_path, show_col_types=FALSE)
  pb   <- read_csv(pb_path,   show_col_types=FALSE)

  combined <- wilc %>%
    select(gene, wilc_log2FC = avg_log2FC, wilc_p_adj = p_val_adj,
           wilc_pct1 = pct.1, wilc_pct2 = pct.2,
           sample_consistency) %>%
    inner_join(pb %>% select(gene, pb_log2FC = logFC, pb_p_adj = FDR,
                              pb_logCPM = logCPM),
               by = "gene") %>%
    mutate(both_up   = wilc_log2FC > 0 & pb_log2FC > 0 & pb_p_adj < 0.05,
           both_down = wilc_log2FC < 0 & pb_log2FC < 0 & pb_p_adj < 0.05,
           direction = case_when(both_up ~ "UP",
                                  both_down ~ "DOWN",
                                  TRUE ~ "discordant_or_ns"))
  combined <- combined %>% arrange(desc(wilc_log2FC * (pb_log2FC > 0)))

  out_path <- file.path(CFG$out_dir, "combined",
                        paste0(ct_id, "_", ct$group1_label, "_vs_",
                               ct$group2_label, "_combined.csv"))
  write_csv(combined, out_path)

  n_up <- sum(combined$both_up, na.rm=TRUE)
  n_down <- sum(combined$both_down, na.rm=TRUE)
  cat_log("  ", ct_id, ": ", nrow(combined), " genes overlap, ",
          n_up, " UP both, ", n_down, " DOWN both\n")
}

cat_log("\n=== DONE: ", format(Sys.time()), " ===\n")
cat_log("Output: ", CFG$out_dir, "\n\n")
