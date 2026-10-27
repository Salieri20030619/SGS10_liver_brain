#!/usr/bin/env Rscript
# =============================================================================
# SGS10 D1 v3.9 (ASCII): Cross-sample integrated PCA + Leiden + smoothing
# -----------------------------------------------------------------------------
# v3.9 KEY CHANGE vs v3.8:
#   v3.8 rbind'd per-sample PCAs --- but each sample's PCA had different
#   axis orientations (e.g. Liver_8 cor(PC1,human)=-0.89, gan56=+0.85),
#   making rbind mathematically invalid. Harmony cant fix this either.
#
#   v3.9 computes a single shared PCA across all samples:
#     - merge raw counts -> NormalizeData -> FindVariableFeatures (2000 HVG)
#     - ScaleData -> RunPCA on merged Seurat object
#   This produces a PCA where PC1 IS the same axis across all samples.
#
# Memory peak: ~20 GB (ScaleData on 276K x 2000 dense)
# Time: 50-65 min total (30-45 min PCA + 20 min downstream)
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(RANN)
  library(igraph)
  library(FNN)
})

option_list <- list(
  make_option("--organ",    type="character", default="both"),
  make_option("--base_dir", type="character",
              default="/mnt/BioDisk/wangmingyuan/SGS10_analysis"),
  make_option("--out_dir",  type="character",
              default="/mnt/BioDisk/wangmingyuan/SGS10_analysis/D1_niche_output"),
  make_option("--n_hvg",    type="integer", default=2000,
              help="number of highly variable genes for integrated PCA"),
  make_option("--n_pcs",    type="integer", default=30),
  make_option("--n_pcs_use",type="integer", default=20,
              help="how many PCs to use downstream"),
  make_option("--k_nn",     type="integer", default=30),
  make_option("--snn_chunk_size", type="integer", default=2000),
  make_option("--snn_prune_threshold", type="double", default=1/15),
  make_option("--target_q_liver", type="integer", default=8),
  make_option("--target_q_brain", type="integer", default=10),
  make_option("--res_scan", type="character",
              default="0.05,0.10,0.15,0.20,0.30,0.50,0.80"),
  make_option("--min_spots_per_niche", type="integer", default=100),
  make_option("--spatial_smooth",    action="store_true", default=TRUE),
  make_option("--no_spatial_smooth", action="store_true", default=FALSE),
  make_option("--smooth_k",          type="integer", default=6),
  make_option("--smooth_maxiter",    type="integer", default=3),
  make_option("--smooth_threshold",  type="double",  default=0.5),
  make_option("--regress_vars", type="character", default="",
              help="comma-sep vars to regress in ScaleData (e.g. nCount_RNA), empty=no regression"),
  make_option("--save_integrated_pca", action="store_true", default=FALSE,
              help="save merged Seurat after PCA (large, for debugging)"),
  make_option("--dry_run",  action="store_true", default=FALSE)
)
opt <- parse_args(OptionParser(option_list=option_list))

CFG_smooth <- if (opt$no_spatial_smooth) FALSE else opt$spatial_smooth

CFG <- list(
  method   = "integrated_PCA + Leiden + chunked_SNN + spatial_smooth",
  organ    = opt$organ,
  base_dir = opt$base_dir,
  rctd_enriched_dir   = file.path(opt$base_dir, "spatial_results/rctd/enriched"),
  rctd_per_sample_dir = file.path(opt$base_dir, "spatial_results/rctd/per_sample"),
  out_dir  = opt$out_dir,
  n_hvg    = opt$n_hvg,
  n_pcs    = opt$n_pcs,
  n_pcs_use= opt$n_pcs_use,
  k_nn     = opt$k_nn,
  snn_chunk_size     = opt$snn_chunk_size,
  snn_prune_threshold= opt$snn_prune_threshold,
  target_q_liver = opt$target_q_liver,
  target_q_brain = opt$target_q_brain,
  res_scan = as.numeric(strsplit(opt$res_scan, ",")[[1]]),
  min_spots_per_niche = opt$min_spots_per_niche,
  spatial_smooth  = CFG_smooth,
  smooth_k        = opt$smooth_k,
  smooth_maxiter  = opt$smooth_maxiter,
  smooth_threshold= opt$smooth_threshold,
  regress_vars = if (nchar(opt$regress_vars) == 0) NULL
                 else strsplit(opt$regress_vars, ",")[[1]],
  save_integrated_pca = opt$save_integrated_pca,
  dry_run = opt$dry_run,
  liver_samples = c("Con_liver_1","Bla_liver_1","Liver_3","Liver_8","gan56"),
  brain_samples = c("Con_brain_1","Bla_brain_1","Brain_3","Brain_8","nao56"),
  meta_keep_fields = c("human_pct","human_count","mouse_count",
                       "species_class","zone","rctd_dominant",
                       "rctd_dom_score","group")
)

dir.create(CFG$out_dir, showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(CFG$out_dir, "logs"),    showWarnings=FALSE)
dir.create(file.path(CFG$out_dir, "figures"), showWarnings=FALSE)

log_file <- file.path(CFG$out_dir, "logs", "D1_run.log")
con_log  <- file(log_file, open="wt")
sink(con_log, split=TRUE, type="output")
sink(con_log, type="message")

cat("=====================================================================\n")
cat("SGS10 D1 v3.9 (Integrated PCA + Leiden)\n")
cat("Start:", format(Sys.time()), "\n")
cat("KEY CHANGE: cross-sample PCA replaces per-sample-rbind\n")
cat("=====================================================================\n\n")

log_section <- function(msg) {
  cat("\n", strrep("-", 70), "\n  ", msg, "\n", strrep("-", 70), "\n", sep="")
}

# ============ Helpers (unchanged from v3.8) ============

build_knn_adj <- function(pca, k = 30) {
  N <- nrow(pca)
  message("    RANN::nn2 on ", N, " x ", ncol(pca), " (k=", k, ")")
  t0 <- Sys.time()
  nn <- nn2(pca, k = k + 1)
  nn_idx <- nn$nn.idx[, -1]
  dt <- round(as.numeric(Sys.time()-t0, units="secs"), 1)
  message("    kNN done. ", dt, "s")

  from <- rep(seq_len(N), each=k)
  to   <- as.vector(t(nn_idx))
  knn_mat <- sparseMatrix(i = from, j = to, x = 1,
                          dims = c(N, N), symmetric = FALSE)
  rm(nn, nn_idx, from, to); gc(verbose=FALSE)
  mem <- round(object.size(knn_mat) / 1024^3, 2)
  message("    knn_mat: ", N, " x ", N, " sparse, ", mem, " GB")
  knn_mat
}

build_snn_chunked <- function(knn_mat, k = 30,
                              chunk_size = 2000,
                              prune_threshold = 1/15) {
  N <- nrow(knn_mat)
  n_chunks <- ceiling(N / chunk_size)
  message("    Chunked SNN: ", n_chunks, " chunks of ~", chunk_size)
  min_common <- prune_threshold * k
  t0 <- Sys.time()
  i_list <- vector("list", n_chunks)
  j_list <- vector("list", n_chunks)
  x_list <- vector("list", n_chunks)
  knn_t <- t(knn_mat)

  for (ch in seq_len(n_chunks)) {
    start_row <- (ch - 1) * chunk_size + 1
    end_row   <- min(ch * chunk_size, N)
    chunk <- knn_mat[start_row:end_row, , drop=FALSE]
    snn_chunk <- chunk %*% knn_t
    snn_chunk@x[snn_chunk@x < min_common] <- 0
    snn_chunk <- drop0(snn_chunk)
    trip <- summary(snn_chunk)
    if (nrow(trip) > 0) {
      trip$i <- trip$i + start_row - 1
      trip <- trip[trip$i != trip$j, ]
      if (nrow(trip) > 0) {
        i_list[[ch]] <- trip$i
        j_list[[ch]] <- trip$j
        x_list[[ch]] <- trip$x / k
      }
    }
    if (ch %% max(1, n_chunks %/% 20) == 0 || ch == n_chunks) {
      dt <- round(as.numeric(Sys.time()-t0, units="secs"), 1)
      mem_mb <- round(sum(as.numeric(gc(verbose=FALSE)[, "(Mb)"])), 0)
      pct <- round(100 * ch / n_chunks, 1)
      message(sprintf("      chunk %d/%d (%.1f%%), %.0fs elapsed, R mem %.0f Mb",
                      ch, n_chunks, pct, dt, mem_mb))
    }
    rm(chunk, snn_chunk, trip); gc(verbose=FALSE)
  }
  rm(knn_t); gc(verbose=FALSE)

  message("    Concatenating edges...")
  all_i <- unlist(i_list); all_j <- unlist(j_list); all_x <- unlist(x_list)
  rm(i_list, j_list, x_list); gc(verbose=FALSE)
  message("      total edges before symmetrization: ", length(all_i))

  message("    Symmetrizing (mirror edges)...")
  t_sym <- Sys.time()
  sym_i <- c(all_i, all_j)
  sym_j <- c(all_j, all_i)
  sym_x <- c(all_x, all_x)
  rm(all_i, all_j, all_x); gc(verbose=FALSE)

  snn <- sparseMatrix(i = sym_i, j = sym_j, x = sym_x,
                      dims = c(N, N), repr = "T")
  diag(snn) <- 0
  snn <- drop0(snn)
  snn <- as(snn, "CsparseMatrix")
  rm(sym_i, sym_j, sym_x); gc(verbose=FALSE)
  dt_sym <- round(as.numeric(Sys.time()-t_sym, units="secs"), 1)
  message("    symmetrize done. ", dt_sym, "s")

  dt <- round(as.numeric(Sys.time()-t0, units="secs"), 1)
  mem <- round(object.size(snn) / 1024^3, 2)
  message("    SNN built: ", length(snn@x), " edges, ", mem, " GB, ",
          round(dt/60, 1), " min")
  snn
}

build_snn_graph <- function(pca, k = 30, chunk_size = 2000,
                            prune_threshold = 1/15) {
  knn_mat <- build_knn_adj(pca, k = k)
  snn <- build_snn_chunked(knn_mat, k = k,
                           chunk_size = chunk_size,
                           prune_threshold = prune_threshold)
  rm(knn_mat); gc(verbose=FALSE)
  message("    Converting to igraph...")
  t0 <- Sys.time()
  g <- graph_from_adjacency_matrix(snn, mode = "undirected",
                                    weighted = TRUE, diag = FALSE)
  dt <- round(as.numeric(Sys.time()-t0, units="secs"), 1)
  message("    igraph: ", vcount(g), " v, ", ecount(g), " e, ", dt, "s")
  g
}

run_louvain <- function(g, resolution = 1.0, seed = 42) {
  set.seed(seed)
  tryCatch({
    r <- cluster_leiden(g, resolution_parameter = resolution,
                        objective_function = "modularity",
                        n_iterations = 3, beta = 0.01)
    as.integer(membership(r))
  }, error = function(e) {
    message("    Leiden failed, fallback Louvain: ", conditionMessage(e))
    r <- cluster_louvain(g, resolution = resolution)
    as.integer(membership(r))
  })
}

spatial_majority_vote <- function(cluster_vec, coord_df,
                                  k = 6, max_iter = 3, threshold = 0.5) {
  stopifnot(length(cluster_vec) == nrow(coord_df))
  message("  Spatial majority-vote (k=", k, ", maxiter=", max_iter,
          ", thresh=", threshold, ")")
  cur <- cluster_vec
  samples <- unique(coord_df$sample)
  neighbor_list <- list()
  for (s in samples) {
    idx_s <- which(coord_df$sample == s)
    if (length(idx_s) < k + 1) next
    coords_s <- as.matrix(coord_df[idx_s, c("x","y")])
    nn <- get.knn(coords_s, k = k)
    nn_global <- matrix(idx_s[nn$nn.index], nrow = length(idx_s), ncol = k)
    neighbor_list[[s]] <- list(idx_global = idx_s, nn = nn_global)
  }
  change_log <- integer(0)
  for (it in seq_len(max_iter)) {
    new_cluster <- cur
    n_changes <- 0L
    for (s in names(neighbor_list)) {
      nl <- neighbor_list[[s]]
      nn_clusters <- matrix(cur[nl$nn], nrow = length(nl$idx_global), ncol = k)
      for (i in seq_along(nl$idx_global)) {
        tab <- table(nn_clusters[i, ])
        top_cl <- as.integer(names(tab)[which.max(tab)])
        top_frac <- max(tab) / k
        gi <- nl$idx_global[i]
        if (top_cl != cur[gi] && top_frac >= threshold) {
          new_cluster[gi] <- top_cl
          n_changes <- n_changes + 1L
        }
      }
    }
    change_log <- c(change_log, n_changes)
    message(sprintf("    iter %d: %d changed (%.2f%%)",
                    it, n_changes, 100 * n_changes / length(cur)))
    cur <- new_cluster
    if (n_changes == 0) { message("    converged"); break }
  }
  list(cluster = cur, change_log = change_log)
}

load_rctd_weights_merged <- function(samples, per_sample_dir) {
  message("Loading per_sample RCTD weights...")
  weights_list <- list()
  for (s in samples) {
    fp <- file.path(per_sample_dir, paste0(s, "_rctd_results.rds"))
    if (!file.exists(fp)) next
    r <- readRDS(fp)
    if (!("weights" %in% names(r@results))) next
    w <- r@results$weights
    rs <- Matrix::rowSums(w); rs[rs == 0] <- 1
    w <- w / rs
    rownames(w) <- paste0(s, "__", rownames(w))
    weights_list[[s]] <- w
    message("  ", s, ": ", nrow(w), " x ", ncol(w))
    rm(r); gc(verbose=FALSE)
  }
  if (length(weights_list) == 0) return(NULL)
  all_cts <- sort(unique(unlist(lapply(weights_list, colnames))))
  aligned <- lapply(weights_list, function(w) {
    missing_cts <- setdiff(all_cts, colnames(w))
    if (length(missing_cts) > 0) {
      zeros <- Matrix(0, nrow=nrow(w), ncol=length(missing_cts), sparse=TRUE,
                      dimnames=list(rownames(w), missing_cts))
      w <- cbind(w, zeros)
    }
    w[, all_cts, drop=FALSE]
  })
  do.call(rbind, aligned)
}

# ============ NEW v3.9: Integrated PCA from raw counts ============

compute_integrated_pca <- function(samples, organ_label,
                                    enriched_dir,
                                    n_hvg = 2000,
                                    n_pcs = 30,
                                    regress_vars = NULL) {
  message("[", organ_label, "] === Integrated PCA pipeline ===")
  t_total <- Sys.time()

  # Step 1: read all samples as Seurat objects, keep raw counts only
  message("  Step 1: loading raw counts from ", length(samples), " samples")
  so_list <- list()
  for (s in samples) {
    fp <- file.path(enriched_dir, paste0(s, "_seekspace_rctd.rds"))
    if (!file.exists(fp)) {
      warning("Missing: ", fp); next
    }
    message("    loading ", s, " ...")
    so <- readRDS(fp)
    # We only need: raw counts + meta + spatial coords
    counts <- GetAssayData(so, layer = "counts")
    meta   <- so@meta.data
    coord  <- if ("spatial" %in% names(so@reductions)) {
      Embeddings(so, "spatial")[, 1:2]
    } else {
      matrix(NA, nrow=ncol(so), ncol=2,
             dimnames=list(colnames(so), c("x","y")))
    }
    # Rebuild as a minimal Seurat object (drop pca etc to save memory)
    new_so <- CreateSeuratObject(counts = counts, meta.data = meta)
    # store coords as misc
    new_so@misc$spatial_coords <- coord
    new_so$sample_id <- s
    so_list[[s]] <- new_so
    rm(so, counts, meta, coord); gc(verbose=FALSE)
    message("      ", s, ": ", ncol(new_so), " spots, ",
            nrow(new_so), " features")
  }
  if (length(so_list) == 0) stop("No samples loaded for ", organ_label)

  # Step 2: merge into one Seurat object
  message("\n  Step 2: merging samples")
  t0 <- Sys.time()
  merged <- if (length(so_list) == 1) so_list[[1]]
            else merge(so_list[[1]], y = so_list[-1],
                       add.cell.ids = names(so_list))
  rm(so_list); gc(verbose=FALSE)
  dt <- round(as.numeric(Sys.time()-t0, units="secs"), 1)
  message("    merged: ", ncol(merged), " spots, ", nrow(merged),
          " features. ", dt, "s")
  message("    sample distribution: ")
  print(table(merged$sample_id))

  # Save coords from merged meta (cell names get prefixed)
  coord_df <- data.frame(
    spot_key = colnames(merged),
    sample   = merged$sample_id,
    stringsAsFactors = FALSE
  )
  # spot_id (without sample prefix)
  coord_df$spot_id <- sub("^[^_]+_[^_]+_?[0-9]*_", "", coord_df$spot_key)
  # Actually safer: use original cell name = spot_key with prefix
  # We'll use spot_key directly as ID. extract spot_id from suffix:
  # add.cell.ids=sample_name -> spot_key = "sample_origCellName"
  coord_df$spot_id_orig <- mapply(function(key, samp) {
    sub(paste0("^", samp, "_"), "", key, perl = FALSE)
  }, coord_df$spot_key, coord_df$sample)

  # Step 3: NormalizeData
  message("\n  Step 3: NormalizeData")
  t0 <- Sys.time()
  merged <- NormalizeData(merged, verbose = FALSE)
  dt <- round(as.numeric(Sys.time()-t0, units="secs"), 1)
  message("    done. ", dt, "s")

  # Step 4: JoinLayers FIRST (Seurat 5: counts/data split per-sample after merge)
  message("\n  Step 4: JoinLayers (Seurat 5 -- before HVG)")
  t0 <- Sys.time()
  tryCatch({
    merged <- JoinLayers(merged)
    message("    JoinLayers done. ",
            round(as.numeric(Sys.time()-t0, units="secs"), 1), "s")
  }, error = function(e) {
    message("    JoinLayers not needed (Seurat <5?): ",
            conditionMessage(e))
  })

  # Step 5: FindVariableFeatures on joined data
  message("\n  Step 5: FindVariableFeatures (n=", n_hvg, ")")
  t0 <- Sys.time()
  merged <- FindVariableFeatures(merged, nfeatures = n_hvg,
                                  selection.method = "vst",
                                  verbose = FALSE)
  hvg <- VariableFeatures(merged)
  dt <- round(as.numeric(Sys.time()-t0, units="secs"), 1)
  message("    done. ", length(hvg), " HVG selected. ", dt, "s")
  if (length(hvg) == 0) {
    stop("FindVariableFeatures returned 0 HVG -- check layer state")
  }

  # Step 6: ScaleData (this is the slow + memory-heavy step)
  message("\n  Step 6: ScaleData on ", length(hvg), " HVG")
  if (!is.null(regress_vars)) {
    message("    regressing out: ", paste(regress_vars, collapse=", "))
  }
  t0 <- Sys.time()
  merged <- ScaleData(merged, features = hvg,
                      vars.to.regress = regress_vars,
                      verbose = FALSE)
  dt <- round(as.numeric(Sys.time()-t0, units="secs"), 1)
  message("    done. ", round(dt/60, 1), " min")
  gc(verbose=FALSE)

  # Step 7: RunPCA
  message("\n  Step 7: RunPCA (n_pcs=", n_pcs, ")")
  t0 <- Sys.time()
  merged <- RunPCA(merged, features = hvg, npcs = n_pcs, verbose = FALSE)
  dt <- round(as.numeric(Sys.time()-t0, units="secs"), 1)
  message("    done. ", round(dt/60, 1), " min")
  gc(verbose=FALSE)

  # Get integrated PCA embedding
  pca <- Embeddings(merged, "pca")
  message("\n  Integrated PCA: ", nrow(pca), " spots x ", ncol(pca), " PCs")

  # Quick QC: check PC1 vs sample_id
  message("  PC1 mean by sample (should differ if tumor signal is in PC1):")
  pc1_by_sample <- tapply(pca[,1], merged$sample_id, mean)
  print(round(pc1_by_sample, 2))

  if ("human_pct" %in% colnames(merged@meta.data)) {
    hp <- merged$human_pct
    if (sd(hp, na.rm=TRUE) > 0) {
      cors <- sapply(1:min(5, ncol(pca)), function(i) {
        cor(pca[,i], hp, use="complete.obs")
      })
      message("  cor(PC1..PC5, human_pct):")
      print(round(cors, 3))
    }
  }

  # Get spatial coords from misc (need to merge per-sample coords into df)
  # Actually we stored them in misc which got dropped after merge.
  # Reload coords from enriched RDS files (faster than re-running PCA on fail)
  message("\n  Step 8: reloading spatial coords")
  coord_full <- list()
  for (s in samples) {
    fp <- file.path(enriched_dir, paste0(s, "_seekspace_rctd.rds"))
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

  # Match PCA spots to coord
  matched <- intersect(rownames(pca), coord_df_full$spot_key)
  message("    matched ", length(matched), " / ", nrow(pca), " spots")
  if (length(matched) < nrow(pca) * 0.95) {
    warning("Many spots cannot be matched to coords -- check naming")
  }
  pca <- pca[matched, , drop=FALSE]
  coord_out <- coord_df_full[matched, c("spot_key","sample","spot_id_orig","x","y")]
  colnames(coord_out)[3] <- "spot_id"

  # Get meta
  meta_keep <- intersect(CFG$meta_keep_fields, colnames(merged@meta.data))
  meta_out <- merged@meta.data[matched, meta_keep, drop=FALSE]
  meta_out$sample   <- merged$sample_id[matched]
  meta_out$spot_id  <- coord_out$spot_id
  meta_out$spot_key <- matched

  dt_total <- round(as.numeric(Sys.time()-t_total, units="secs"), 1)
  message("\n[", organ_label, "] integrated PCA total: ",
          round(dt_total/60, 1), " min")

  rm(merged); gc(verbose=FALSE)

  list(pca = pca, coord = coord_out, meta = meta_out)
}

# ============ Environment check ============

log_section("Environment check")
check_inputs <- function() {
  samples <- switch(CFG$organ,
                    liver = CFG$liver_samples,
                    brain = CFG$brain_samples,
                    both  = c(CFG$liver_samples, CFG$brain_samples))
  miss <- character()
  for (s in samples) {
    fp <- file.path(CFG$rctd_enriched_dir, paste0(s, "_seekspace_rctd.rds"))
    if (!file.exists(fp)) miss <- c(miss, fp)
  }
  if (length(miss) > 0) {
    cat("  [MISSING]:\n"); for (m in miss) cat("    -", m, "\n")
    stop("Missing")
  }
  cat("  [OK] inputs + Seurat", as.character(packageVersion("Seurat")),
      "+ RANN + igraph", as.character(packageVersion("igraph")), "\n")
}
check_inputs()

if (CFG$dry_run) {
  cat("\n[DRY RUN] OK\n")
  sink(type="message"); sink(type="output"); close(con_log)
  quit(save="no", status=0)
}

# ============ Section: Integrated PCA per organ ============

log_section("Integrated PCA computation")
organs_to_run <- switch(CFG$organ, liver="liver", brain="brain",
                        both=c("liver","brain"))

organ_data <- list()
for (org in organs_to_run) {
  samples <- if (org == "liver") CFG$liver_samples else CFG$brain_samples
  combined <- compute_integrated_pca(samples, org,
                                      CFG$rctd_enriched_dir,
                                      n_hvg = CFG$n_hvg,
                                      n_pcs = CFG$n_pcs,
                                      regress_vars = CFG$regress_vars)
  # Trim to n_pcs_use
  combined$pca <- combined$pca[, seq_len(min(CFG$n_pcs_use, ncol(combined$pca))), drop=FALSE]
  organ_data[[org]] <- list(organ=org, samples=samples, combined=combined)
  gc(verbose=FALSE)
}

# ============ Graph clustering (same as v3.8) ============

log_section("Graph clustering (chunked SNN + Leiden + smoothing)")

run_graph_clustering <- function(org_dat) {
  org    <- org_dat$organ
  pca    <- org_dat$combined$pca
  coord  <- org_dat$combined$coord
  target_q <- if (org == "liver") CFG$target_q_liver else CFG$target_q_brain

  message("[", org, "] ", nrow(pca), " spots x ", ncol(pca),
          " PCs, target q=", target_q)
  t_all <- Sys.time()

  g <- build_snn_graph(pca, k = CFG$k_nn,
                       chunk_size = CFG$snn_chunk_size,
                       prune_threshold = CFG$snn_prune_threshold)
  gc(verbose=FALSE)

  message("  Leiden resolution scan...")
  scan_res <- data.frame(resolution=numeric(), n_clusters=integer(),
                          n_small=integer())
  best_diff <- Inf; best_res <- NA; best_clusters <- NULL
  for (res in CFG$res_scan) {
    t0 <- Sys.time()
    cl <- run_louvain(g, resolution = res)
    dt <- round(as.numeric(Sys.time()-t0, units="secs"), 1)
    n_cl <- length(unique(cl))
    n_small <- sum(table(cl) < CFG$min_spots_per_niche)
    scan_res <- rbind(scan_res,
                      data.frame(resolution=res, n_clusters=n_cl, n_small=n_small))
    message(sprintf("    res=%.3f: n_cl=%d, small=%d, %.1fs",
                    res, n_cl, n_small, dt))
    diff_val <- abs(n_cl - target_q) + 0.3 * n_small
    if (diff_val < best_diff) {
      best_diff <- diff_val; best_res <- res; best_clusters <- cl
    }
  }
  message("  -> best res = ", best_res, " (q = ", length(unique(best_clusters)), ")")
  scan_res$is_best <- scan_res$resolution == best_res
  p <- ggplot(scan_res, aes(x=resolution, y=n_clusters)) +
    geom_line(color="steelblue") +
    geom_point(aes(color=is_best), size=3) +
    scale_color_manual(values=c("FALSE"="gray50","TRUE"="red"), guide="none") +
    geom_hline(yintercept=target_q, linetype="dashed", color="red", alpha=0.5) +
    labs(title=paste0(org, " resolution scan"),
         x="resolution", y="n_clusters") +
    theme_minimal(base_size=11)
  ggsave(file.path(CFG$out_dir, "figures", paste0("res_scan_", org, ".pdf")),
         p, width=6, height=4)
  rm(g); gc(verbose=FALSE)

  smoothed_clusters <- best_clusters
  change_log <- integer(0)
  if (CFG$spatial_smooth) {
    message("\n  Spatial smoothing")
    sr <- spatial_majority_vote(best_clusters, coord,
                                k = CFG$smooth_k,
                                max_iter = CFG$smooth_maxiter,
                                threshold = CFG$smooth_threshold)
    smoothed_clusters <- sr$cluster
    change_log        <- sr$change_log
    message(sprintf("  -> total %d changed (%.2f%%)",
                    sum(change_log),
                    100 * sum(change_log) / length(best_clusters)))
  }

  prefix <- toupper(substr(org, 1, 1))
  niches <- data.frame(
    spot_key = rownames(coord),
    sample   = coord$sample,
    spot_id  = coord$spot_id,
    x = coord$x, y = coord$y,
    niche_raw = sprintf("%s-N%02d", prefix, best_clusters),
    niche     = sprintf("%s-N%02d", prefix, smoothed_clusters)
  )

  dt_total <- round(as.numeric(Sys.time()-t_all, units="secs"), 1)
  message("  Clustering total: ", round(dt_total/60, 1), " min")

  list(niches=niches,
       best_q=length(unique(smoothed_clusters)),
       best_res=best_res,
       scan_df=scan_res,
       smooth_change_log=change_log)
}

niche_results <- lapply(organ_data, run_graph_clustering)

# ============ QC + outputs (same as v3.8) ============

log_section("Per-niche QC")
spot_count_qc <- list()
for (org in names(niche_results)) {
  ni <- niche_results[[org]]$niches
  tab <- ni %>% count(niche, sample) %>%
    pivot_wider(names_from=sample, values_from=n, values_fill=0) %>%
    mutate(total = rowSums(across(where(is.numeric)))) %>%
    arrange(desc(total))
  tab$organ <- org
  spot_count_qc[[org]] <- tab
  cat("\n[", org, "]\n")
  for (i in seq_len(min(nrow(tab), 15))) {
    mk <- if (tab$total[i] < CFG$min_spots_per_niche) " [SMALL]" else ""
    cat(sprintf("    %s : %6d%s\n", tab$niche[i], tab$total[i], mk))
  }
}
organ_suffix <- if (CFG$organ == "both") "" else paste0("_", CFG$organ)
write_csv(bind_rows(spot_count_qc),
          file.path(CFG$out_dir, paste0("niche_spot_count_qc", organ_suffix, ".csv")))

log_section("Celltype composition (RCTD weights)")
for (org in names(niche_results)) {
  samples <- if (org == "liver") CFG$liver_samples else CFG$brain_samples
  weights <- load_rctd_weights_merged(samples, CFG$rctd_per_sample_dir)
  if (is.null(weights)) next
  ni <- niche_results[[org]]$niches
  common <- intersect(ni$spot_key, rownames(weights))
  message("[", org, "] intersect: ", length(common),
          " (", round(100*length(common)/nrow(ni),1), "%)")
  if (length(common) < 100) next
  w_sub <- as.matrix(weights[common, , drop=FALSE])
  ni_sub <- ni[match(common, ni$spot_key), ]
  df <- as.data.frame(w_sub); df$niche <- ni_sub$niche
  comp <- df %>% group_by(niche) %>%
    summarise(n_rctd_spots = n(),
              across(everything(), mean), .groups="drop")
  comp$organ <- org
  write_csv(comp, file.path(CFG$out_dir,
                            paste0("niche_celltype_composition_", org, ".csv")))
  top3 <- comp %>% select(-organ, -n_rctd_spots) %>%
    pivot_longer(-niche, names_to="celltype", values_to="proportion") %>%
    group_by(niche) %>% arrange(desc(proportion)) %>%
    slice_head(n=3) %>% mutate(rank=row_number()) %>%
    select(niche, rank, celltype, proportion) %>%
    pivot_wider(names_from=rank, values_from=c(celltype, proportion))
  write_csv(top3, file.path(CFG$out_dir,
                            paste0("niche_top3_celltypes_", org, ".csv")))
  rm(weights, w_sub, df); gc(verbose=FALSE)
}

log_section("Tumor-feature preview")
tumor_rows <- list()
for (org in names(niche_results)) {
  ni   <- niche_results[[org]]$niches
  meta <- organ_data[[org]]$combined$meta
  merged_df <- merge(ni, meta[, c("spot_key", "human_pct")],
                     by="spot_key", all.x=TRUE)
  if ("human_pct" %in% colnames(merged_df)) {
    stats <- merged_df %>% group_by(niche) %>%
      summarise(
        n_spots = n(),
        mean_human_pct = mean(human_pct, na.rm=TRUE),
        med_human_pct  = median(human_pct, na.rm=TRUE),
        pct_tumor_spot = sum(human_pct > 50, na.rm=TRUE) / n(),
        pct_mixed_spot = sum(human_pct > 5 & human_pct <= 50, na.rm=TRUE) / n(),
        pct_host_spot  = sum(human_pct <= 5, na.rm=TRUE) / n(),
        .groups="drop")
    stats$organ <- org
    tumor_rows[[org]] <- stats
    cat("\n[", org, "] top by pct_tumor_spot:\n")
    stats %>% arrange(desc(pct_tumor_spot)) %>% head(15) %>%
      as.data.frame() %>% print()
  }
}
if (length(tumor_rows) > 0) {
  write_csv(bind_rows(tumor_rows),
            file.path(CFG$out_dir,
                      paste0("niche_tumor_feature_preview", organ_suffix, ".csv")))
}

log_section("Main outputs")
all_ni <- bind_rows(lapply(niche_results, function(r) r$niches))
niches_by_sample <- split(all_ni, all_ni$sample)
attr(niches_by_sample, "method") <- CFG$method
attr(niches_by_sample, "q_used") <- sapply(niche_results, function(r) r$best_q)
attr(niches_by_sample, "best_res") <- sapply(niche_results, function(r) r$best_res)
saveRDS(niches_by_sample,
        file.path(CFG$out_dir, paste0("niches_by_sample", organ_suffix, ".rds")))

prop <- all_ni %>% count(sample, niche) %>%
  group_by(sample) %>% mutate(prop = n / sum(n)) %>% ungroup() %>%
  select(niche, sample, prop) %>%
  pivot_wider(names_from=sample, values_from=prop, values_fill=0)
write_csv(prop, file.path(CFG$out_dir,
                          paste0("niche_proportion_per_sample", organ_suffix, ".csv")))

scan_all <- bind_rows(lapply(names(niche_results), function(org) {
  df <- niche_results[[org]]$scan_df; df$organ <- org; df
}))
write_csv(scan_all, file.path(CFG$out_dir,
                              paste0("niche_resolution_scan", organ_suffix, ".csv")))

if (CFG$spatial_smooth) {
  smooth_log <- bind_rows(lapply(names(niche_results), function(org) {
    cl <- niche_results[[org]]$smooth_change_log
    data.frame(organ=org, iter=seq_along(cl), n_changed=cl)
  }))
  write_csv(smooth_log, file.path(CFG$out_dir,
                                  paste0("spatial_smoothing_log", organ_suffix, ".csv")))
}

log_section("Spatial viz")
for (org in names(niche_results)) {
  ni <- niche_results[[org]]$niches
  samples <- unique(ni$sample)
  p_list <- lapply(samples, function(s) {
    sub <- ni[ni$sample==s, ]
    ggplot(sub, aes(x=x, y=y, color=niche)) +
      geom_point(size=0.3, alpha=0.8) +
      coord_equal() + theme_void() +
      labs(title=paste0(s, " (n=", nrow(sub), ")")) +
      theme(legend.position="none",
            plot.title=element_text(size=10, hjust=0.5))
  })
  p <- wrap_plots(p_list, ncol=min(5, length(samples))) +
    plot_annotation(
      title=paste(org, "niche map (integrated PCA, smoothed)"),
      subtitle=paste0("q=", niche_results[[org]]$best_q,
                      " | res=", niche_results[[org]]$best_res))
  ggsave(file.path(CFG$out_dir, "figures", paste0("niche_map_", org, ".pdf")),
         p, width=16, height=5, limitsize=FALSE)

  if (CFG$spatial_smooth && "niche_raw" %in% colnames(ni)) {
    p_list_raw <- lapply(samples, function(s) {
      sub <- ni[ni$sample==s, ]
      ggplot(sub, aes(x=x, y=y, color=niche_raw)) +
        geom_point(size=0.3, alpha=0.8) + coord_equal() + theme_void() +
        labs(title=paste0(s, " (raw)")) +
        theme(legend.position="none",
              plot.title=element_text(size=9, hjust=0.5))
    })
    p_list_sm <- lapply(samples, function(s) {
      sub <- ni[ni$sample==s, ]
      ggplot(sub, aes(x=x, y=y, color=niche)) +
        geom_point(size=0.3, alpha=0.8) + coord_equal() + theme_void() +
        labs(title=paste0(s, " (smoothed)")) +
        theme(legend.position="none",
              plot.title=element_text(size=9, hjust=0.5))
    })
    p_cmp <- wrap_plots(c(p_list_raw, p_list_sm),
                        ncol=length(samples), nrow=2) +
      plot_annotation(title=paste(org, "niche: raw (top) vs smoothed (bottom)"))
    ggsave(file.path(CFG$out_dir, "figures",
                     paste0("niche_map_raw_vs_smooth_", org, ".pdf")),
           p_cmp, width=16, height=10, limitsize=FALSE)
  }
}

log_section("Summary")
for (org in names(niche_results)) {
  cat("\n[", org, "] res=", niche_results[[org]]$best_res,
      " q=", niche_results[[org]]$best_q, "\n")
  tab <- sort(table(niche_results[[org]]$niches$niche), decreasing=TRUE)
  cat("  total:", sum(tab), " | range:", min(tab), "-", max(tab), "\n")
  if (CFG$spatial_smooth) {
    cl <- niche_results[[org]]$smooth_change_log
    if (length(cl) > 0) {
      cat("  smooth: ", length(cl), "iter, changed",
          sprintf("%.2f%%", 100*sum(cl)/sum(tab)), "\n")
    }
  }
}

cat("\nDone:", format(Sys.time()), "\n")
sink(type="message"); sink(type="output"); close(con_log)
