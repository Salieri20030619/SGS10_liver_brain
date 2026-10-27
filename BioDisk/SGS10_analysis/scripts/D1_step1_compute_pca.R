#!/usr/bin/env Rscript
# =============================================================================
# SGS10 D1 step1 v2: Compute integrated PCA on MOUSE genes only
# -----------------------------------------------------------------------------
# IMPORTANT: source RDS files have inconsistent gene naming:
#   baseline samples (Con/Bla): no prefix, mouse-only (32285 genes)
#   tumor samples: mm10---xxx + GRCh38-xxx (68886 genes)
# This script unifies to MOUSE-ONLY (strip mm10--- prefix), giving 32285
# common features across all 5 samples. PCA describes host (mouse) response.
# Human tumor genes are dropped here -- to be analyzed separately in D2.
# -----------------------------------------------------------------------------
# Splits v3.9 into two stages:
#   Step1 (this) : compute integrated PCA, save to RDS (~45 min, ONE-TIME)
#   Step2        : load cached PCA, run SNN+Leiden+smoothing (~20 min, RE-RUN)
#
# This way, parameter tuning (resolution, smoothing) only needs Step2.
#
# Output:
#   <out_dir>/cache/<organ>_integrated_pca.rds
#       contains: list(pca, coord, meta, qc_stats)
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(Matrix)
  library(dplyr)
})

option_list <- list(
  make_option("--organ",    type="character", default="liver"),
  make_option("--base_dir", type="character",
              default="/mnt/BioDisk/wangmingyuan/SGS10_analysis"),
  make_option("--out_dir",  type="character",
              default="/mnt/BioDisk/wangmingyuan/SGS10_analysis/D1_niche_output"),
  make_option("--n_hvg",    type="integer", default=2000),
  make_option("--n_pcs",    type="integer", default=30),
  make_option("--regress_vars", type="character", default="")
)
opt <- parse_args(OptionParser(option_list=option_list))

CFG <- list(
  organ = opt$organ,
  enriched_dir = file.path(opt$base_dir, "spatial_results/rctd/enriched"),
  out_dir = opt$out_dir,
  cache_dir = file.path(opt$out_dir, "cache"),
  n_hvg = opt$n_hvg,
  n_pcs = opt$n_pcs,
  regress_vars = if (nchar(opt$regress_vars) == 0) NULL
                 else strsplit(opt$regress_vars, ",")[[1]],
  liver_samples = c("Con_liver_1","Bla_liver_1","Liver_3","Liver_8","gan56"),
  brain_samples = c("Con_brain_1","Bla_brain_1","Brain_3","Brain_8","nao56"),
  meta_keep_fields = c("human_pct","human_count","mouse_count",
                       "species_class","zone","rctd_dominant",
                       "rctd_dom_score","group")
)

dir.create(CFG$cache_dir, showWarnings=FALSE, recursive=TRUE)

# Use cat() instead of message() to ensure output goes to stdout
cat_log <- function(...) cat(..., sep="")
sec <- function(name) {
  cat_log("\n", strrep("-", 70), "\n  ", name, "\n",
          strrep("-", 70), "\n")
}

cat_log("=====================================================================\n")
cat_log("D1 Step 1: Integrated PCA caching\n")
cat_log("Organ: ", CFG$organ, "\n")
cat_log("Start: ", format(Sys.time()), "\n")
cat_log("=====================================================================\n")

samples <- if (CFG$organ == "liver") CFG$liver_samples else CFG$brain_samples

# ============ Step 1: Load all samples ============
sec("Step 1: Load samples")
so_list <- list()
for (s in samples) {
  fp <- file.path(CFG$enriched_dir, paste0(s, "_seekspace_rctd.rds"))
  if (!file.exists(fp)) {
    stop("Missing: ", fp)
  }
  cat_log("  loading ", s, " ...\n")
  so <- readRDS(fp)
  counts <- GetAssayData(so, layer = "counts")
  meta   <- so@meta.data
  coord  <- if ("spatial" %in% names(so@reductions)) {
    Embeddings(so, "spatial")[, 1:2]
  } else {
    matrix(NA, nrow=ncol(so), ncol=2,
           dimnames=list(colnames(so), c("x","y")))
  }

  # ==== KEY: unify gene names to mouse-only, 32285 features ====
  # Naming inconsistency in source RDS:
  #   baseline (Con/Bla):  no prefix       (e.g. Xkr4, Gm1992)
  #   tumor (Liver_3/8/gan56): mm10---xxx + GRCh38-xxx
  # Strategy: keep only mouse genes, strip "mm10---" prefix to align with baseline.
  feats <- rownames(counts)
  has_mm10 <- any(grepl("^mm10---", feats))
  has_grch <- any(grepl("^GRCh38-", feats))

  if (has_mm10) {
    # tumor sample: keep mm10---*, strip prefix
    keep_idx <- grep("^mm10---", feats)
    counts <- counts[keep_idx, , drop=FALSE]
    rownames(counts) <- sub("^mm10---", "", rownames(counts))
    cat_log("    [", s, "] stripped mm10--- prefix, kept ",
            nrow(counts), " mouse genes (dropped ",
            sum(has_grch * (length(feats) - length(keep_idx))),
            " human genes)\n")
  } else {
    cat_log("    [", s, "] no mm10--- prefix, keep as-is (", nrow(counts), " genes)\n")
  }

  new_so <- CreateSeuratObject(counts = counts, meta.data = meta)
  new_so@misc$spatial_coords <- coord
  new_so$sample_id <- s
  so_list[[s]] <- new_so
  rm(so, counts, meta, coord); gc(verbose=FALSE)
  cat_log("    ", s, ": ", ncol(new_so), " spots, ",
          nrow(new_so), " features\n")
}

# ============ Step 2: Merge ============
sec("Step 2: Merge")
t0 <- Sys.time()
merged <- if (length(so_list) == 1) {
  so_list[[1]]
} else {
  merge(so_list[[1]], y = so_list[-1],
        add.cell.ids = names(so_list))
}
rm(so_list); gc(verbose=FALSE)
dt <- round(as.numeric(Sys.time()-t0, units="secs"), 1)
cat_log("  merged: ", ncol(merged), " spots, ",
        nrow(merged), " features. ", dt, "s\n")
cat_log("  sample distribution:\n")
print(table(merged$sample_id))

# ============ Step 3: NormalizeData ============
sec("Step 3: NormalizeData")
t0 <- Sys.time()
merged <- NormalizeData(merged, verbose = FALSE)
cat_log("  done. ",
        round(as.numeric(Sys.time()-t0, units="secs"), 1), "s\n")

# ============ Step 4: JoinLayers ============
sec("Step 4: JoinLayers (Seurat 5)")
t0 <- Sys.time()
tryCatch({
  merged <- JoinLayers(merged)
  cat_log("  done. ",
          round(as.numeric(Sys.time()-t0, units="secs"), 1), "s\n")
}, error = function(e) {
  cat_log("  not needed or failed: ", conditionMessage(e), "\n")
})

# ============ Step 5: FindVariableFeatures ============
sec("Step 5: FindVariableFeatures")
t0 <- Sys.time()
merged <- FindVariableFeatures(merged, nfeatures = CFG$n_hvg,
                                selection.method = "vst", verbose = FALSE)
hvg <- VariableFeatures(merged)
cat_log("  ", length(hvg), " HVG selected. ",
        round(as.numeric(Sys.time()-t0, units="secs"), 1), "s\n")
if (length(hvg) == 0) stop("0 HVG -- check layer state")

# ============ Step 6: ScaleData ============
sec("Step 6: ScaleData (slow, ~15-25 min)")
t0 <- Sys.time()
merged <- ScaleData(merged, features = hvg,
                    vars.to.regress = CFG$regress_vars, verbose = FALSE)
cat_log("  done. ",
        round(as.numeric(Sys.time()-t0, units="secs")/60, 1), " min\n")
gc(verbose=FALSE)

# ============ Step 7: RunPCA ============
sec("Step 7: RunPCA")
t0 <- Sys.time()
merged <- RunPCA(merged, features = hvg, npcs = CFG$n_pcs, verbose = FALSE)
cat_log("  done. ",
        round(as.numeric(Sys.time()-t0, units="secs")/60, 1), " min\n")
gc(verbose=FALSE)

# ============ Step 8: Extract & QC ============
sec("Step 8: Extract PCA + QC")
pca <- Embeddings(merged, "pca")
cat_log("  PCA dim: ", nrow(pca), " x ", ncol(pca), "\n")

cat_log("\n  PC1 mean by sample (baseline should cluster, tumor should cluster):\n")
pc1_by_sample <- tapply(pca[,1], merged$sample_id, mean)
print(round(pc1_by_sample, 2))

qc_stats <- list()
qc_stats$pc1_by_sample <- pc1_by_sample

if ("human_pct" %in% colnames(merged@meta.data)) {
  hp <- merged$human_pct
  if (sd(hp, na.rm=TRUE) > 0) {
    cors <- sapply(1:min(5, ncol(pca)), function(i) {
      cor(pca[,i], hp, use="complete.obs")
    })
    cat_log("\n  cor(PC1..PC5, human_pct) [global]:\n")
    print(round(cors, 3))
    qc_stats$cor_human_pct <- cors
  }
}

# ============ Step 9: Reload spatial coords ============
sec("Step 9: Reload spatial coords")
coord_full <- list()
for (s in samples) {
  fp <- file.path(CFG$enriched_dir, paste0(s, "_seekspace_rctd.rds"))
  if (!file.exists(fp)) next
  so <- readRDS(fp)
  if (!"spatial" %in% names(so@reductions)) {
    rm(so); next
  }
  cc <- as.data.frame(Embeddings(so, "spatial"))
  colnames(cc)[1:2] <- c("x","y")
  cc$spot_id_orig <- rownames(cc)
  cc$sample <- s
  cc$spot_key <- paste0(s, "_", cc$spot_id_orig)
  coord_full[[s]] <- cc
  rm(so); gc(verbose=FALSE)
}
coord_df_full <- do.call(rbind, coord_full)
rownames(coord_df_full) <- coord_df_full$spot_key
cat_log("  total coords: ", nrow(coord_df_full), "\n")

# ============ Step 10: Match PCA to coords ============
sec("Step 10: Match PCA to coords")
matched <- intersect(rownames(pca), coord_df_full$spot_key)
cat_log("  matched: ", length(matched), " / ", nrow(pca),
        " (", round(100*length(matched)/nrow(pca), 1), "%)\n")
if (length(matched) < nrow(pca) * 0.95) {
  warning("Many spots cannot be matched to coords (",
          length(matched), " of ", nrow(pca), ")")
}
pca <- pca[matched, , drop=FALSE]
coord_out <- coord_df_full[matched, c("spot_key","sample","spot_id_orig","x","y")]
colnames(coord_out)[3] <- "spot_id"

meta_keep <- intersect(CFG$meta_keep_fields, colnames(merged@meta.data))
meta_out <- merged@meta.data[matched, meta_keep, drop=FALSE]
meta_out$sample   <- merged$sample_id[matched]
meta_out$spot_id  <- coord_out$spot_id
meta_out$spot_key <- matched

# ============ Step 11: Save cache ============
sec("Step 11: Save PCA cache")
cache_file <- file.path(CFG$cache_dir,
                        paste0(CFG$organ, "_integrated_pca.rds"))
saveRDS(list(
  pca = pca,
  coord = coord_out,
  meta = meta_out,
  qc = qc_stats,
  config = CFG,
  date = format(Sys.time())
), cache_file)
cat_log("  saved: ", cache_file, "\n")
cat_log("  size: ", round(file.info(cache_file)$size / 1024^2, 1), " MB\n")

cat_log("\n=== DONE: ", format(Sys.time()), " ===\n")
cat_log("Next step: run D1_step2_cluster.R --organ ", CFG$organ, "\n\n")
