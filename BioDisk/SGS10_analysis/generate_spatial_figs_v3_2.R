#!/usr/bin/env Rscript
# =============================================================================
# SGS10 Task B Spatial Figures v3.1 (cross-sample merged UMAP + group overlay)
# Date: 2026-04-20
# Author: Claude
#
# 主要改动 vs v2:
#   - UMAP: 改为跨样本合并 (5 liver 合一, 5 brain 合一)
#           产出 3 张 merged UMAP: by niche / by sample / by group (tumor vs control)
#   - k-means: 在合并的 weights 矩阵上全局跑, 所有样本共享 niche 标签
#   - Spatial Domain: 单样本图, 但 niche 颜色跨样本统一
#   - Overlay: 保持 v2 的 4-color argmax 逻辑
# 产出: 3 merged UMAP × 2 organ + 10 × 2 单样本 = 26 张 PNG
#
# Input:
#   ~/SGS10_analysis/spatial_results/rctd/per_sample/<sample>_rctd_results.rds
#
# Output:
#   ~/SGS10_analysis/taskB_spatial_figs_2026-04-20_v3/{liver,brain}/*
#
# Usage:
#   mkdir -p ~/SGS10_analysis/taskB_spatial_figs_v3
#   # scp generate_spatial_figs_v3.R to that dir
#   cd ~/SGS10_analysis/
#   screen -S taskB_figs_v3
#   Rscript taskB_spatial_figs_v3/generate_spatial_figs_v3.R 2>&1 | tee run_v3.log
# =============================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(uwot)
  library(FNN)
  library(readr)
  library(stringr)
  library(RColorBrewer)
  library(scales)
})

# ---------- 全局配置 ----------
RCTD_DIR <- "~/SGS10_analysis/spatial_results/rctd/per_sample"
OUT_ROOT <- "~/SGS10_analysis/taskB_spatial_figs_2026-04-20_v3"
DPI      <- 300
K_MEANS  <- 10
KNN_SMOOTH <- 15
SET_SEED <- 42

dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)
LOG_FILE <- file.path(OUT_ROOT, "run.log")
writeLines(c(paste("# SGS10 Task B spatial figs v3 run log -", Sys.time()), ""),
           LOG_FILE)

log_msg <- function(msg) {
  line <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", msg)
  cat(line, "\n")
  cat(line, "\n", file = LOG_FILE, append = TRUE, sep = "")
}

# ---------- 样本定义 ----------
SAMPLES_LIVER <- c(
  "Bla_liver_1" = "baseline",
  "Con_liver_1" = "baseline",
  "Liver_3"     = "tumor",
  "Liver_8"     = "tumor",
  "gan56"       = "tumor"
)
SAMPLES_BRAIN <- c(
  "Bla_brain_1" = "baseline",
  "Con_brain_1" = "baseline",
  "Brain_3"     = "tumor",
  "Brain_8"     = "tumor",
  "nao56"       = "tumor"
)

CORE_CT_LIVER <- c("HSC_Activated_ECM", "Macrophage_M1",
                   "Kupffer",           "Hepatocyte_PV")
CORE_CT_BRAIN <- c("Microglia_Homeostatic", "Pericytes",
                   "Astrocyte_Resting",     "Endothelial")

CORE_COLORS <- c("#E31A1C", "#1F78B4", "#33A02C", "#FF7F00")

# 10 色调色板 (niche)
NICHE_PAL <- c(brewer.pal(9, "Set1"), "#666666")

# 样本配色 (5 色)
SAMPLE_PAL <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")

# ---------- paper theme ----------
paper_theme <- theme_bw(base_size = 11) +
  theme(panel.grid = element_blank(),
        axis.text  = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, color = "grey30"),
        legend.position = "right",
        legend.title = element_text(face = "bold"),
        aspect.ratio = 1)

# ---------- 工具: 从 RCTD 对象提取 weights + coords ----------
extract_rctd <- function(rctd_file) {
  rctd <- readRDS(rctd_file)

  w <- NULL
  if (!is.null(rctd@results$weights)) {
    w <- as.matrix(rctd@results$weights)
  } else if (!is.null(rctd@results$weights_doublet)) {
    if (!is.null(rctd@results$norm_weights)) {
      w <- as.matrix(rctd@results$norm_weights)
    } else {
      w <- as.matrix(rctd@results$weights_doublet)
    }
  }
  if (is.null(w)) stop("Cannot find weights in RCTD object")

  coords <- rctd@spatialRNA@coords
  common <- intersect(rownames(w), rownames(coords))
  if (length(common) == 0) {
    common <- intersect(colnames(rctd@spatialRNA@counts), rownames(w))
  }
  w <- w[common, , drop = FALSE]
  coords <- coords[common, , drop = FALSE]

  row_sums <- rowSums(w)
  row_sums[row_sums == 0] <- 1
  w <- w / row_sums

  list(weights = w, coords = coords, n_spots = nrow(w),
       n_ct = ncol(w), celltypes = colnames(w))
}

# ---------- KNN 空间平滑 ----------
smooth_domain_knn <- function(coords, cluster, k = KNN_SMOOTH, iters = 2) {
  xy <- as.matrix(coords[, 1:2])
  if (!is.numeric(xy[1, 1])) xy <- apply(xy, 2, as.numeric)

  nn <- FNN::get.knn(xy, k = k)$nn.index

  current <- cluster
  for (it in seq_len(iters)) {
    new_cluster <- integer(length(current))
    for (i in seq_along(current)) {
      neigh_labels <- c(current[i], current[nn[i, ]])
      tab <- table(neigh_labels)
      new_cluster[i] <- as.integer(names(tab)[which.max(tab)])
    }
    current <- new_cluster
  }
  current
}

# ---------- 自动 domain 命名 (based on cluster centers) ----------
name_domains <- function(centers) {
  k <- nrow(centers)
  out <- character(k)
  for (i in seq_len(k)) {
    top2 <- sort(centers[i, ], decreasing = TRUE)[1:2]
    out[i] <- sprintf("D%02d [%s+%s]", i, names(top2)[1], names(top2)[2])
  }
  out
}

# =============================================================================
# 主处理流程 - 按器官分别跑
# =============================================================================
process_organ <- function(organ, samples_dict, core_cts) {
  t0 <- Sys.time()
  log_msg(sprintf("\n================ ORGAN: %s ================", organ))

  # 确保 organ 子目录存在 (fix for v3.1)
  dir.create(file.path(OUT_ROOT, organ), recursive = TRUE, showWarnings = FALSE)

  # ---------- 阶段 1: 读入所有样本 ----------
  log_msg("Phase 1: loading all RCTD objects...")
  all_data <- list()
  for (s in names(samples_dict)) {
    rctd_file <- file.path(RCTD_DIR,
                           sprintf("%s_rctd_results.rds", s))
    rctd_file <- path.expand(rctd_file)
    if (!file.exists(rctd_file)) {
      log_msg(sprintf("  [SKIP] %s not found", s))
      next
    }
    log_msg(sprintf("  reading %s...", s))
    d <- tryCatch(extract_rctd(rctd_file),
                  error = function(e) {
                    log_msg(sprintf("  [ERR] %s: %s", s, e$message))
                    NULL
                  })
    if (is.null(d)) next
    all_data[[s]] <- d
  }
  if (length(all_data) == 0) {
    log_msg(sprintf("[FATAL] no %s data loaded", organ))
    return(FALSE)
  }

  # ---------- 阶段 2: 合并 weights 矩阵 (按共同 celltype) ----------
  log_msg("Phase 2: merging weights across samples (common celltypes)...")
  common_cts <- Reduce(intersect, lapply(all_data, function(d) d$celltypes))
  log_msg(sprintf("  common celltypes: %d", length(common_cts)))

  merged_w <- do.call(rbind, lapply(all_data, function(d) d$weights[, common_cts]))
  # 记录每行来自哪个样本
  sample_labels <- unlist(lapply(names(all_data), function(s)
                         rep(s, all_data[[s]]$n_spots)))
  log_msg(sprintf("  merged matrix: %d spots × %d celltypes",
                  nrow(merged_w), ncol(merged_w)))

  # ---------- 阶段 3: 全局 k-means ----------
  log_msg("Phase 3: global k-means on merged matrix...")
  set.seed(SET_SEED)
  km_global <- kmeans(merged_w, centers = K_MEANS,
                      iter.max = 50, nstart = 10)
  domain_names <- name_domains(km_global$centers)
  log_msg(sprintf("  k-means done (%d clusters)", K_MEANS))
  for (i in seq_len(K_MEANS)) {
    log_msg(sprintf("    %s (n=%d)", domain_names[i],
                    sum(km_global$cluster == i)))
  }

  # ---------- 阶段 4: 合并 UMAP ----------
  log_msg("Phase 4: running merged UMAP (30-120s)...")
  set.seed(SET_SEED)
  umap_res <- uwot::umap(merged_w,
                         n_neighbors = 50,
                         min_dist    = 0.1,
                         spread      = 1.5,
                         metric      = "cosine",
                         verbose     = FALSE)
  colnames(umap_res) <- c("UMAP1", "UMAP2")

  umap_df <- data.frame(umap_res,
                        sample = sample_labels,
                        domain = domain_names[km_global$cluster],
                        group = samples_dict[sample_labels])
  umap_df$domain <- factor(umap_df$domain, levels = domain_names)
  umap_df$sample <- factor(umap_df$sample, levels = names(samples_dict))

  # ----- 画图 0a: merged UMAP by niche -----
  log_msg("  plot 0a: merged UMAP by niche")
  p1 <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = domain)) +
    geom_point(size = 0.8, alpha = 0.5) +
    scale_color_manual(values = NICHE_PAL) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1),
                                ncol = 1)) +
    labs(title = sprintf("Merged Niche UMAP - %s (5 samples)", toupper(organ)),
         subtitle = sprintf("n = %d spots, k-means k=%d on merged RCTD weights",
                           nrow(umap_df), K_MEANS),
         x = "UMAP1", y = "UMAP2", color = "Niche") +
    paper_theme +
    theme(legend.text = element_text(size = 7),
          legend.key.size = unit(0.4, "cm"))

  ggsave(file.path(OUT_ROOT, organ, "00_merged_niche_umap.png"),
         p1, width = 9, height = 6.5, dpi = DPI, bg = "white")

  # ----- 画图 0b: merged UMAP by sample (诊断 batch) -----
  log_msg("  plot 0b: merged UMAP by sample")
  p2 <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = sample)) +
    geom_point(size = 0.8, alpha = 0.5) +
    scale_color_manual(values = SAMPLE_PAL) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1),
                                ncol = 1)) +
    labs(title = sprintf("Merged UMAP by sample - %s", toupper(organ)),
         subtitle = "Color = sample (diagnostic for batch effect)",
         x = "UMAP1", y = "UMAP2", color = "Sample") +
    paper_theme +
    theme(legend.text = element_text(size = 8))

  ggsave(file.path(OUT_ROOT, organ, "00_merged_umap_by_sample.png"),
         p2, width = 9, height = 6.5, dpi = DPI, bg = "white")

  # ----- 画图 0c: merged UMAP by group (tumor vs control) - 核心对比图 -----
  log_msg("  plot 0c: merged UMAP by group (tumor vs control)")
  umap_df$group <- factor(umap_df$group, levels = c("baseline", "tumor"))

  # 计算 tumor-enriched / control-enriched / shared regions
  # 用 bin2d 粗格 + tumor 比例计算, 只在图说明里给总数
  n_tumor    <- sum(umap_df$group == "tumor")
  n_baseline <- sum(umap_df$group == "baseline")

  # 优化叠绘: tumor 在上层 (因为关心 tumor 特异 spot 的分布)
  p3 <- ggplot() +
    geom_point(data = umap_df %>% filter(group == "baseline"),
               aes(x = UMAP1, y = UMAP2),
               color = "#3690C0", size = 0.8, alpha = 0.35) +
    geom_point(data = umap_df %>% filter(group == "tumor"),
               aes(x = UMAP1, y = UMAP2),
               color = "#E31A1C", size = 0.8, alpha = 0.45) +
    # 用一个 geom_point 填 legend
    geom_point(data = data.frame(UMAP1 = c(Inf, Inf),
                                 UMAP2 = c(Inf, Inf),
                                 group = factor(c("baseline", "tumor"),
                                                levels = c("baseline", "tumor"))),
               aes(x = UMAP1, y = UMAP2, color = group), size = 3) +
    scale_color_manual(values = c("baseline" = "#3690C0",
                                  "tumor"    = "#E31A1C")) +
    guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +
    labs(title = sprintf("Merged UMAP by group - %s", toupper(organ)),
         subtitle = sprintf("Red = tumor (n=%d spots, 3 samples), Blue = control (n=%d spots, 2 samples)",
                            n_tumor, n_baseline),
         x = "UMAP1", y = "UMAP2", color = "Group") +
    paper_theme +
    theme(legend.text = element_text(size = 9))

  ggsave(file.path(OUT_ROOT, organ, "00_merged_umap_by_group.png"),
         p3, width = 9, height = 6.5, dpi = DPI, bg = "white")

  # ---------- 阶段 5: 按样本切分 cluster labels, 画单样本 Spatial Domain ----------
  log_msg("Phase 5: per-sample spatial domain plots (unified niche colors)...")
  # 按样本把 cluster 切回
  row_idx <- 0
  for (s in names(all_data)) {
    t_sample <- Sys.time()
    d <- all_data[[s]]
    n_s <- d$n_spots

    # 对应的 global cluster slice
    cluster_s <- km_global$cluster[(row_idx + 1) : (row_idx + n_s)]
    row_idx <- row_idx + n_s

    # KNN 空间平滑 (用该样本自己的坐标)
    log_msg(sprintf("  %s: KNN smoothing...", s))
    smoothed <- smooth_domain_knn(d$coords, cluster_s,
                                  k = KNN_SMOOTH, iters = 2)

    out_dir <- file.path(OUT_ROOT, organ, s)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    # ----- Spatial Domain 图 (统一 niche 颜色) -----
    coord_df <- as.data.frame(d$coords)
    colnames(coord_df)[1:2] <- c("x", "y")
    coord_df$domain <- factor(domain_names[smoothed], levels = domain_names)

    grp <- samples_dict[s]
    p_sd <- ggplot(coord_df, aes(x = x, y = y, color = domain)) +
      geom_point(size = 1.2, alpha = 0.8) +
      scale_color_manual(values = NICHE_PAL, drop = FALSE) +
      scale_y_reverse() + coord_fixed() +
      guides(color = guide_legend(override.aes = list(size = 3, alpha = 1),
                                  ncol = 1)) +
      labs(title = sprintf("Spatial Domain - %s (%s)", s, grp),
           subtitle = sprintf("%d spots, niche colors unified across %s samples (merged k-means)",
                              n_s, organ),
           x = "Spatial x", y = "Spatial y", color = "Niche") +
      paper_theme +
      theme(legend.text = element_text(size = 7),
            legend.key.size = unit(0.4, "cm"),
            axis.text = element_blank(),
            axis.ticks = element_blank())

    ggsave(file.path(out_dir, "02_spatial_domain.png"),
           p_sd, width = 9, height = 7, dpi = DPI, bg = "white")

    # ----- 4-color Overlay (保持 v2 逻辑) -----
    avail_cts <- intersect(core_cts, colnames(d$weights))
    if (length(avail_cts) > 0) {
      w_core <- d$weights[, avail_cts, drop = FALSE]
      argmax_idx <- max.col(w_core, ties.method = "first")
      argmax_weight <- w_core[cbind(seq_len(nrow(w_core)), argmax_idx)]
      argmax_ct <- ifelse(argmax_weight >= 0.1,
                          avail_cts[argmax_idx], "Other")

      df_ov <- as.data.frame(d$coords)
      colnames(df_ov)[1:2] <- c("x", "y")
      df_ov$TopCT <- factor(argmax_ct, levels = c(core_cts, "Other"))

      df_other <- df_ov %>% filter(TopCT == "Other")
      df_core  <- df_ov %>% filter(TopCT != "Other")

      pal_ov <- setNames(c(CORE_COLORS[1:length(core_cts)], "grey85"),
                         c(core_cts, "Other"))

      p_ov <- ggplot() +
        geom_point(data = df_other, aes(x = x, y = y),
                   color = "grey85", size = 0.8, alpha = 0.6) +
        geom_point(data = df_core,
                   aes(x = x, y = y, color = TopCT),
                   size = 1.3, alpha = 0.9) +
        scale_color_manual(values = pal_ov, drop = FALSE) +
        scale_y_reverse() + coord_fixed() +
        guides(color = guide_legend(override.aes = list(size = 3, alpha = 1),
                                    ncol = 1)) +
        labs(title = sprintf("Cell-type Overlay - %s (%s)", s, grp),
             subtitle = "argmax of 4 core CT (weight >= 0.10), grey = other",
             x = "Spatial x", y = "Spatial y", color = "Top cell type") +
        paper_theme +
        theme(legend.text = element_text(size = 8),
              legend.key.size = unit(0.4, "cm"),
              axis.text = element_blank(),
              axis.ticks = element_blank())

      ggsave(file.path(out_dir, "03_celltype_overlay_4color.png"),
             p_ov, width = 9, height = 7, dpi = DPI, bg = "white")
    }

    dt_s <- as.numeric(difftime(Sys.time(), t_sample, units = "secs"))
    log_msg(sprintf("    %s done in %.1f s", s, dt_s))
  }

  # ---------- 保存 cluster 标签到 CSV (供后续分析) ----------
  out_csv <- file.path(OUT_ROOT, organ,
                       sprintf("%s_merged_cluster_labels.csv", organ))
  write.csv(data.frame(sample = sample_labels,
                       group = samples_dict[sample_labels],
                       global_cluster = km_global$cluster,
                       domain = domain_names[km_global$cluster]),
            out_csv, row.names = FALSE)
  log_msg(sprintf("Saved cluster labels to: %s", basename(out_csv)))

  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  log_msg(sprintf("================ %s done in %.1f s ================", organ, dt))
  TRUE
}

# ---------- 主调度 ----------
log_msg(paste0("RCTD_DIR: ", path.expand(RCTD_DIR)))
log_msg(paste0("OUT_ROOT: ", path.expand(OUT_ROOT)))
log_msg("Strategy: merged cross-sample UMAP + k-means, unified niche colors")
log_msg("")

ok_liver <- process_organ("liver", SAMPLES_LIVER, CORE_CT_LIVER)
ok_brain <- process_organ("brain", SAMPLES_BRAIN, CORE_CT_BRAIN)

# ---------- 打包 ----------
log_msg("")
log_msg("Packaging tar.gz...")
system(sprintf("cd %s && cd .. && tar -czf taskB_spatial_figs_2026-04-20_v3.tar.gz taskB_spatial_figs_2026-04-20_v3/",
               path.expand(OUT_ROOT)))
tar_path <- "~/SGS10_analysis/taskB_spatial_figs_2026-04-20_v3.tar.gz"
log_msg(sprintf("Tarball: %s", path.expand(tar_path)))
system(sprintf("ls -lh %s", path.expand(tar_path)))
log_msg(paste0("End: ", Sys.time()))
