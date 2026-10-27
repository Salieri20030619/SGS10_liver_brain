#!/usr/bin/env Rscript
# =============================================================================
# SGS10 D1 step2: From cached PCA, run SNN+Leiden+smoothing+output
# -----------------------------------------------------------------------------
# Input: <out_dir>/cache/<organ>_integrated_pca.rds (from step1)
# Output: niches_by_sample_<organ>.rds, all CSVs, figures
#
# Time: ~20 min (no PCA recomputation)
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(RANN)
  library(igraph)
  library(FNN)
  library(Seurat)  # for RCTD weights handling
})

option_list <- list(
  make_option("--organ",    type="character", default="liver"),
  make_option("--base_dir", type="character",
              default="/mnt/BioDisk/wangmingyuan/SGS10_analysis"),
  make_option("--out_dir",  type="character",
              default="/mnt/BioDisk/wangmingyuan/SGS10_analysis/D1_niche_output"),
  make_option("--n_pcs_use",type="integer", default=20),
  make_option("--k_nn",     type="integer", default=30),
  make_option("--snn_chunk_size", type="integer", default=2000),
  make_option("--snn_prune_threshold", type="double", default=1/15),
  make_option("--target_q_liver", type="integer", default=8),
  make_option("--target_q_brain", type="integer", default=10),
  make_option("--res_scan", type="character",
              default="0.05,0.10,0.15,0.20,0.30,0.50,0.80"),
  make_option("--min_spots_per_niche", type="integer", default=100),
  make_option("--no_spatial_smooth", action="store_true", default=FALSE),
  make_option("--smooth_k",          type="integer", default=6),
  make_option("--smooth_maxiter",    type="integer", default=3),
  make_option("--smooth_threshold",  type="double",  default=0.5)
)
opt <- parse_args(OptionParser(option_list=option_list))

CFG <- list(
  organ = opt$organ,
  base_dir = opt$base_dir,
  out_dir = opt$out_dir,
  cache_dir = file.path(opt$out_dir, "cache"),
  rctd_per_sample_dir = file.path(opt$base_dir, "spatial_results/rctd/per_sample"),
  n_pcs_use = opt$n_pcs_use,
  k_nn = opt$k_nn,
  snn_chunk_size = opt$snn_chunk_size,
  snn_prune_threshold = opt$snn_prune_threshold,
  target_q_liver = opt$target_q_liver,
  target_q_brain = opt$target_q_brain,
  res_scan = as.numeric(strsplit(opt$res_scan, ",")[[1]]),
  min_spots_per_niche = opt$min_spots_per_niche,
  spatial_smooth = !opt$no_spatial_smooth,
  smooth_k = opt$smooth_k,
  smooth_maxiter = opt$smooth_maxiter,
  smooth_threshold = opt$smooth_threshold,
  liver_samples = c("Con_liver_1","Bla_liver_1","Liver_3","Liver_8","gan56"),
  brain_samples = c("Con_brain_1","Bla_brain_1","Brain_3","Brain_8","nao56")
)

dir.create(file.path(CFG$out_dir, "figures"), showWarnings=FALSE, recursive=TRUE)

cat_log <- function(...) cat(..., sep="")
sec <- function(name) {
  cat_log("\n", strrep("-", 70), "\n  ", name, "\n",
          strrep("-", 70), "\n")
}

cat_log("=====================================================================\n")
cat_log("D1 Step 2: Cluster on cached PCA\n")
cat_log("Organ: ", CFG$organ, "\n")
cat_log("Start: ", format(Sys.time()), "\n")
cat_log("=====================================================================\n")

# ============ Load cache ============
sec("Load PCA cache")
cache_file <- file.path(CFG$cache_dir,
                        paste0(CFG$organ, "_integrated_pca.rds"))
if (!file.exists(cache_file)) {
  stop("Cache not found: ", cache_file,
       "\n  Run D1_step1_compute_pca.R --organ ", CFG$organ, " first.")
}
cache <- readRDS(cache_file)
pca   <- cache$pca[, seq_len(min(CFG$n_pcs_use, ncol(cache$pca))), drop=FALSE]
coord <- cache$coord
meta  <- cache$meta
cat_log("  PCA: ", nrow(pca), " x ", ncol(pca), "\n")
cat_log("  cache date: ", cache$date, "\n")

if (!is.null(cache$qc$cor_human_pct)) {
  cat_log("  cor(PC1..PC5, human_pct):\n")
  print(round(cache$qc$cor_human_pct, 3))
}

# ============ Helper functions (same as v3.9) ============

build_knn_adj <- function(pca, k = 30) {
  N <- nrow(pca)
  cat_log("    RANN::nn2 on ", N, " x ", ncol(pca), " (k=", k, ")\n")
  t0 <- Sys.time()
  nn <- nn2(pca, k = k + 1)
  nn_idx <- nn$nn.idx[, -1]
  cat_log("    kNN done. ",
          round(as.numeric(Sys.time()-t0, units="secs"), 1), "s\n")
  from <- rep(seq_len(N), each=k)
  to   <- as.vector(t(nn_idx))
  knn_mat <- sparseMatrix(i = from, j = to, x = 1,
                          dims = c(N, N), symmetric = FALSE)
  rm(nn, nn_idx, from, to); gc(verbose=FALSE)
  cat_log("    knn_mat: ", N, " x ", N, " sparse, ",
          round(object.size(knn_mat) / 1024^3, 2), " GB\n")
  knn_mat
}

build_snn_chunked <- function(knn_mat, k = 30,
                              chunk_size = 2000,
                              prune_threshold = 1/15) {
  N <- nrow(knn_mat)
  n_chunks <- ceiling(N / chunk_size)
  cat_log("    Chunked SNN: ", n_chunks, " chunks of ~", chunk_size, "\n")
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
      cat_log(sprintf("      chunk %d/%d (%.1f%%), %.0fs, mem %.0f Mb\n",
                      ch, n_chunks, 100*ch/n_chunks, dt, mem_mb))
    }
    rm(chunk, snn_chunk, trip); gc(verbose=FALSE)
  }
  rm(knn_t); gc(verbose=FALSE)

  cat_log("    Concat edges...\n")
  all_i <- unlist(i_list); all_j <- unlist(j_list); all_x <- unlist(x_list)
  rm(i_list, j_list, x_list); gc(verbose=FALSE)
  cat_log("      total edges: ", length(all_i), "\n")

  cat_log("    Symmetrizing...\n")
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
  cat_log("    sym done. ",
          round(as.numeric(Sys.time()-t_sym, units="secs"), 1), "s\n")
  cat_log("    SNN built: ", length(snn@x), " edges, ",
          round(object.size(snn) / 1024^3, 2), " GB\n")
  snn
}

build_snn_graph <- function(pca, k = 30, chunk_size = 2000,
                            prune_threshold = 1/15) {
  knn_mat <- build_knn_adj(pca, k = k)
  snn <- build_snn_chunked(knn_mat, k = k,
                           chunk_size = chunk_size,
                           prune_threshold = prune_threshold)
  rm(knn_mat); gc(verbose=FALSE)
  cat_log("    Convert to igraph...\n")
  t0 <- Sys.time()
  g <- graph_from_adjacency_matrix(snn, mode = "undirected",
                                    weighted = TRUE, diag = FALSE)
  cat_log("    igraph: ", vcount(g), " v, ", ecount(g), " e, ",
          round(as.numeric(Sys.time()-t0, units="secs"), 1), "s\n")
  g
}

run_leiden <- function(g, resolution = 1.0, seed = 42) {
  set.seed(seed)
  tryCatch({
    r <- cluster_leiden(g, resolution_parameter = resolution,
                        objective_function = "modularity",
                        n_iterations = 3, beta = 0.01)
    as.integer(membership(r))
  }, error = function(e) {
    cat_log("    Leiden failed, fallback Louvain: ", conditionMessage(e), "\n")
    r <- cluster_louvain(g, resolution = resolution)
    as.integer(membership(r))
  })
}

spatial_majority_vote <- function(cluster_vec, coord_df,
                                  k = 6, max_iter = 3, threshold = 0.5) {
  cat_log("  Smoothing (k=", k, ", maxiter=", max_iter,
          ", thresh=", threshold, ")\n")
  cur <- cluster_vec
  samples_local <- unique(coord_df$sample)
  neighbor_list <- list()
  for (s in samples_local) {
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
    cat_log(sprintf("    iter %d: %d changed (%.2f%%)\n",
                    it, n_changes, 100 * n_changes / length(cur)))
    cur <- new_cluster
    if (n_changes == 0) { cat_log("    converged\n"); break }
  }
  list(cluster = cur, change_log = change_log)
}

# ============ Build SNN ============
sec("Build SNN graph")
g <- build_snn_graph(pca, k = CFG$k_nn,
                     chunk_size = CFG$snn_chunk_size,
                     prune_threshold = CFG$snn_prune_threshold)

# ============ Resolution scan ============
sec("Leiden resolution scan")
target_q <- if (CFG$organ == "liver") CFG$target_q_liver else CFG$target_q_brain
cat_log("  target q = ", target_q, "\n")
scan_res <- data.frame(resolution=numeric(), n_clusters=integer(),
                        n_small=integer())
best_diff <- Inf; best_res <- NA; best_clusters <- NULL
for (res in CFG$res_scan) {
  t0 <- Sys.time()
  cl <- run_leiden(g, resolution = res)
  dt <- round(as.numeric(Sys.time()-t0, units="secs"), 1)
  n_cl <- length(unique(cl))
  n_small <- sum(table(cl) < CFG$min_spots_per_niche)
  scan_res <- rbind(scan_res,
                    data.frame(resolution=res, n_clusters=n_cl, n_small=n_small))
  cat_log(sprintf("    res=%.3f: n_cl=%d, small=%d, %.1fs\n",
                  res, n_cl, n_small, dt))
  diff_val <- abs(n_cl - target_q) + 0.3 * n_small
  if (diff_val < best_diff) {
    best_diff <- diff_val; best_res <- res; best_clusters <- cl
  }
}
cat_log("  best res = ", best_res, " (q = ",
        length(unique(best_clusters)), ")\n")

scan_res$is_best <- scan_res$resolution == best_res
p <- ggplot(scan_res, aes(x=resolution, y=n_clusters)) +
  geom_line(color="steelblue") +
  geom_point(aes(color=is_best), size=3) +
  scale_color_manual(values=c("FALSE"="gray50","TRUE"="red"), guide="none") +
  geom_hline(yintercept=target_q, linetype="dashed", color="red", alpha=0.5) +
  labs(title=paste0(CFG$organ, " resolution scan"),
       x="resolution", y="n_clusters") +
  theme_minimal(base_size=11)
ggsave(file.path(CFG$out_dir, "figures", paste0("res_scan_", CFG$organ, ".pdf")),
       p, width=6, height=4)
rm(g); gc(verbose=FALSE)

# ============ Spatial smoothing ============
sec("Spatial smoothing")
smoothed_clusters <- best_clusters
change_log <- integer(0)
if (CFG$spatial_smooth) {
  sr <- spatial_majority_vote(best_clusters, coord,
                              k = CFG$smooth_k,
                              max_iter = CFG$smooth_maxiter,
                              threshold = CFG$smooth_threshold)
  smoothed_clusters <- sr$cluster
  change_log <- sr$change_log
  cat_log(sprintf("  total %d changed (%.2f%%)\n",
                  sum(change_log),
                  100 * sum(change_log) / length(best_clusters)))
}

# ============ Build niches data ============
prefix <- toupper(substr(CFG$organ, 1, 1))
niches <- data.frame(
  spot_key = rownames(coord),
  sample   = coord$sample,
  spot_id  = coord$spot_id,
  x = coord$x, y = coord$y,
  niche_raw = sprintf("%s-N%02d", prefix, best_clusters),
  niche     = sprintf("%s-N%02d", prefix, smoothed_clusters)
)

# ============ Per-niche QC ============
sec("Per-niche QC")
tab <- niches %>%
  count(niche, sample) %>%
  pivot_wider(names_from=sample, values_from=n, values_fill=0) %>%
  mutate(total = rowSums(across(where(is.numeric)))) %>%
  arrange(desc(total))
tab$organ <- CFG$organ
cat_log("\n  Top 15 niches by spot count:\n")
for (i in seq_len(min(nrow(tab), 15))) {
  mk <- if (tab$total[i] < CFG$min_spots_per_niche) " [SMALL]" else ""
  cat_log(sprintf("    %s : %6d%s\n", tab$niche[i], tab$total[i], mk))
}
write_csv(tab, file.path(CFG$out_dir,
                          paste0("niche_spot_count_qc_", CFG$organ, ".csv")))

# ============ Tumor-feature preview ============
sec("Tumor-feature preview")
merged_meta <- merge(niches, meta[, c("spot_key", "human_pct")],
                     by="spot_key", all.x=TRUE)
if ("human_pct" %in% colnames(merged_meta)) {
  stats <- merged_meta %>% group_by(niche) %>%
    summarise(
      n_spots = n(),
      mean_human_pct = mean(human_pct, na.rm=TRUE),
      med_human_pct  = median(human_pct, na.rm=TRUE),
      pct_tumor_spot = sum(human_pct > 50, na.rm=TRUE) / n(),
      pct_mixed_spot = sum(human_pct > 5 & human_pct <= 50, na.rm=TRUE) / n(),
      pct_host_spot  = sum(human_pct <= 5, na.rm=TRUE) / n(),
      .groups="drop")
  stats$organ <- CFG$organ
  cat_log("\n  Top by pct_tumor_spot:\n")
  print(as.data.frame(stats %>% arrange(desc(pct_tumor_spot))))
  write_csv(stats, file.path(CFG$out_dir,
                             paste0("niche_tumor_feature_preview_",
                                    CFG$organ, ".csv")))
}

# ============ Celltype composition ============
sec("Celltype composition (RCTD weights)")
samples <- unique(niches$sample)
weights_list <- list()
for (s in samples) {
  fp <- file.path(CFG$rctd_per_sample_dir, paste0(s, "_rctd_results.rds"))
  if (!file.exists(fp)) next
  r <- readRDS(fp)
  if (!("weights" %in% names(r@results))) next
  w <- r@results$weights
  rs <- Matrix::rowSums(w); rs[rs == 0] <- 1
  w <- w / rs
  # IMPORTANT: niches$spot_key is "sample_origname" (one underscore) from merge add.cell.ids
  # but RCTD weights rownames are origname. Need to add prefix matching the merge format:
  rownames(w) <- paste0(s, "_", rownames(w))
  weights_list[[s]] <- w
  cat_log("  ", s, ": ", nrow(w), " x ", ncol(w), "\n")
  rm(r); gc(verbose=FALSE)
}

if (length(weights_list) > 0) {
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
  weights <- do.call(rbind, aligned)

  common <- intersect(niches$spot_key, rownames(weights))
  cat_log("  intersect: ", length(common), " (",
          round(100*length(common)/nrow(niches),1), "%)\n")

  if (length(common) >= 100) {
    w_sub <- as.matrix(weights[common, , drop=FALSE])
    ni_sub <- niches[match(common, niches$spot_key), ]
    df <- as.data.frame(w_sub); df$niche <- ni_sub$niche
    comp <- df %>% group_by(niche) %>%
      summarise(n_rctd_spots = n(),
                across(everything(), mean), .groups="drop")
    comp$organ <- CFG$organ
    write_csv(comp, file.path(CFG$out_dir,
                              paste0("niche_celltype_composition_",
                                     CFG$organ, ".csv")))
    top3 <- comp %>% select(-organ, -n_rctd_spots) %>%
      pivot_longer(-niche, names_to="celltype", values_to="proportion") %>%
      group_by(niche) %>% arrange(desc(proportion)) %>%
      slice_head(n=3) %>% mutate(rank=row_number()) %>%
      select(niche, rank, celltype, proportion) %>%
      pivot_wider(names_from=rank, values_from=c(celltype, proportion))
    write_csv(top3, file.path(CFG$out_dir,
                              paste0("niche_top3_celltypes_",
                                     CFG$organ, ".csv")))
    cat_log("  written\n")
  }
  rm(weights); gc(verbose=FALSE)
}

# ============ Main outputs ============
sec("Main outputs")
niches_by_sample <- split(niches, niches$sample)
attr(niches_by_sample, "method") <- "integrated_PCA + Leiden + smoothing"
attr(niches_by_sample, "q_used") <- length(unique(smoothed_clusters))
attr(niches_by_sample, "best_res") <- best_res
attr(niches_by_sample, "n_pcs") <- ncol(pca)
saveRDS(niches_by_sample,
        file.path(CFG$out_dir, paste0("niches_by_sample_", CFG$organ, ".rds")))

prop <- niches %>% count(sample, niche) %>%
  group_by(sample) %>% mutate(prop = n / sum(n)) %>% ungroup() %>%
  select(niche, sample, prop) %>%
  pivot_wider(names_from=sample, values_from=prop, values_fill=0)
write_csv(prop, file.path(CFG$out_dir,
                          paste0("niche_proportion_per_sample_", CFG$organ, ".csv")))

scan_res$organ <- CFG$organ
write_csv(scan_res, file.path(CFG$out_dir,
                              paste0("niche_resolution_scan_", CFG$organ, ".csv")))

if (CFG$spatial_smooth) {
  smooth_log <- data.frame(organ=CFG$organ,
                           iter=seq_along(change_log),
                           n_changed=change_log)
  write_csv(smooth_log, file.path(CFG$out_dir,
                                  paste0("spatial_smoothing_log_", CFG$organ, ".csv")))
}

# ============ Spatial viz ============
sec("Spatial viz")
samples <- unique(niches$sample)
p_list <- lapply(samples, function(s) {
  sub <- niches[niches$sample==s, ]
  ggplot(sub, aes(x=x, y=y, color=niche)) +
    geom_point(size=0.3, alpha=0.8) +
    coord_equal() + theme_void() +
    labs(title=paste0(s, " (n=", nrow(sub), ")")) +
    theme(legend.position="none",
          plot.title=element_text(size=10, hjust=0.5))
})
p <- wrap_plots(p_list, ncol=min(5, length(samples))) +
  plot_annotation(
    title=paste(CFG$organ, "niche map (integrated PCA, smoothed)"),
    subtitle=paste0("q=", length(unique(smoothed_clusters)),
                    " | res=", best_res))
ggsave(file.path(CFG$out_dir, "figures",
                 paste0("niche_map_", CFG$organ, ".pdf")),
       p, width=16, height=5, limitsize=FALSE)

if (CFG$spatial_smooth) {
  p_list_raw <- lapply(samples, function(s) {
    sub <- niches[niches$sample==s, ]
    ggplot(sub, aes(x=x, y=y, color=niche_raw)) +
      geom_point(size=0.3, alpha=0.8) + coord_equal() + theme_void() +
      labs(title=paste0(s, " (raw)")) +
      theme(legend.position="none",
            plot.title=element_text(size=9, hjust=0.5))
  })
  p_list_sm <- lapply(samples, function(s) {
    sub <- niches[niches$sample==s, ]
    ggplot(sub, aes(x=x, y=y, color=niche)) +
      geom_point(size=0.3, alpha=0.8) + coord_equal() + theme_void() +
      labs(title=paste0(s, " (smoothed)")) +
      theme(legend.position="none",
            plot.title=element_text(size=9, hjust=0.5))
  })
  p_cmp <- wrap_plots(c(p_list_raw, p_list_sm),
                      ncol=length(samples), nrow=2) +
    plot_annotation(title=paste(CFG$organ, "raw (top) vs smoothed (bottom)"))
  ggsave(file.path(CFG$out_dir, "figures",
                   paste0("niche_map_raw_vs_smooth_", CFG$organ, ".pdf")),
         p_cmp, width=16, height=10, limitsize=FALSE)
}

cat_log("\n=== DONE: ", format(Sys.time()), " ===\n")
cat_log("[", CFG$organ, "] q=", length(unique(smoothed_clusters)),
        " res=", best_res, " smoothing=",
        sprintf("%.1f%%", 100*sum(change_log)/length(best_clusters)), "\n\n")
