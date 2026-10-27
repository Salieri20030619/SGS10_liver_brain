#!/usr/bin/env Rscript
# ============================================================
# pseudotime_v4_core.R --- v4 注释拟时序分析 (合并分析 + 核心细类)
# ============================================================
# 用法:
#   cd ~/SGS10_analysis
#   Rscript scripts/pseudotime_v4_core.R 2>&1 | tee scripts/logs/pseudotime_v4.log
#
# 和原 pseudotime_multicurve_optimized 脚本的关键区别:
#   1. 使用 v4_refined 注释 (而不是 final_celltype)
#   2. 使用 v4 RDS (通过 diff_common.R 的 resolve_rds_path 解析路径)
#   3. 合并所有样本 + Harmony 整合, 得到统一伪时间坐标系
#      (不再分组独立跑, 避免组间 pseudotime 不可比的方法学问题)
#   4. 只对 10 个核心细类做拟时序
#      (跳过 All_Cells×2 和 Macrophage_NAM, 因 NAM 细胞数不足)
#   5. tradeSeq 正式版: fitGAM + associationTest + patternTest
#      (patternTest 用于比较 tumor 和 control 细胞在同一轨迹上的动态差异)
#   6. 复用 Stage 1 的 diff_common.R 模块
#
# 前置:
#   - Stage 1/2/3 已完成, diff_common.R 已在 scripts/ 下
#   - v4 RDS 存在于 <sample>/mouse/ 子目录
#
# 输出:
#   pseudotime_v4_results/
#     integrated/          两大整合对象的 UMAP 概览
#     {liver,brain}/       每个细类一个子目录
#       <virtual_ct>/
#         trajectory.{pdf,png}
#         curve_assignment.{pdf,png}
#         pseudotime_by_group.{pdf,png}    ← tumor vs ctrl 在同一轨迹上的分布
#         spearman_correlation.csv
#         heatmap_top_genes.pdf
#         dynamic_core_genes.{pdf,png}     ← 三轴核心基因动态
#         tradeseq_results/
#           association_test.csv
#           pattern_test.csv               ← tumor 曲线 vs control 曲线差异
#           core_panel_fit.rds
#     summary/
#       task_status.csv
#       pseudotime_summary.md
#       error_log.txt
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(slingshot)
  library(tradeSeq)
  library(harmony)
  library(BiocParallel)
  library(FNN)
  library(parallel)
  library(future)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(RColorBrewer)
  library(viridis)
  library(pheatmap)
})

# dplyr namespace lock
select    <- dplyr::select
rename    <- dplyr::rename
filter    <- dplyr::filter
mutate    <- dplyr::mutate
arrange   <- dplyr::arrange
summarise <- dplyr::summarise
group_by  <- dplyr::group_by
bind_rows <- dplyr::bind_rows
pull      <- dplyr::pull

# ============================================================
# 0. 路径 & 共用模块
# ============================================================
PROJECT_ROOT <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
SCRIPT_DIR   <- file.path(PROJECT_ROOT, "scripts")
setwd(PROJECT_ROOT)

source(file.path(SCRIPT_DIR, "diff_common.R"))

OUT_DIR <- file.path(PROJECT_ROOT, "pseudotime_v4_results")
for (d in c(
  "integrated",
  "liver", "brain",
  "summary"
)) dir.create(file.path(OUT_DIR, d), recursive = TRUE, showWarnings = FALSE)

ERROR_LOG <- file.path(OUT_DIR, "summary", "error_log.txt")
writeLines(sprintf("== pseudotime v4 error log - %s ==",
                   format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ERROR_LOG)

write_error <- function(context, msg) {
  line <- sprintf("[%s] %s => %s",
                  format(Sys.time(), "%H:%M:%S"), context, msg)
  cat("  ⚠", line, "\n")
  cat(line, "\n", file = ERROR_LOG, append = TRUE)
}

# ============================================================
# 1. 参数
# ============================================================

# 要跑的核心细类 (10 个, 去掉 All_Cells×2 + NAM)
CELL_TYPES_TO_RUN <- list(
  liver = c("Hepatocyte_PP", "Hepatocyte_PV", "HSC",
            "HSC_Activated_ECM", "Kupffer", "DC"),   # 6 个
  brain = c("Endothelial", "Microglia_Homeostatic",
            "Microglia_DAM", "Astrocyte_Gfap")        # 4 个
)
# 总计: 6 + 4 = 10 个细类

# 每个细类独立处理时的参数
HARMONY_DIMS <- 30
UMAP_DIMS    <- 30
KMEANS_K     <- 15       # slingshot 用的 kmeans 聚类数
KNN_K        <- 10       # pseudotime NA 填充的邻居数

# tradeSeq 参数
N_KNOTS_GAM        <- 6       # GAM 平滑节点数
TOP_GENES_FOR_GAM  <- 500     # Spearman top 500 跑 GAM
MIN_CELLS_FOR_RUN  <- 500     # 细类最少细胞数, 低于此跳过

# 核心基因面板 (用于 dynamic_core_genes 图 + tradeSeq 必跑)
CORE_PANEL <- unique(c(
  # 轴 1: BBB + 脑内皮 MHC-I
  "H2-D1", "H2-K1", "B2m", "Cdh5", "Cldn5", "Klf4", "Kdr", "Flt1", "Bsg",
  # 轴 2: Efferocytosis
  "Elmo1", "Dock4", "Rac1", "Itgb5", "Ifitm10", "Mertk", "Gas6", "Axl",
  "Lrp1", "Pld4", "Marcks", "Arhgdia", "Arhgap17", "Cx3cr1",
  # 轴 3: TGFb/IL-6
  "Smad3", "Tgfb1", "Il6st", "Il6", "Stat3", "Jak1", "Jak2", "Socs3",
  "Nfkb1", "Rela",
  # Microglia / Astrocyte 标志
  "P2ry12", "Tmem119", "Trem2", "Cst7", "Apoe", "Gfap", "Slc1a2", "S100b",
  # Hepatocyte zonation
  "Cyp2e1", "Cyp1a2", "Glul", "Pck1", "Alb",
  # HSC activation
  "Lrat", "Acta2", "Col1a1", "Col1a2", "Des",
  # Kupffer / DC
  "Clec4f", "Flt3", "Itgax", "Irf8"
))

# 每个细类的 root marker (用于自动判定 pseudotime 方向)
ROOT_MARKERS <- list(
  Hepatocyte_PP         = "Cyp2e1",  # PP zone 高表达 Cyp2e1 (实际 PP 低，作为方向参考)
  Hepatocyte_PV         = "Cyp2e1",  # 同上
  HSC                   = "Lrat",    # quiescent HSC 高表达 Lrat
  HSC_Activated_ECM     = "Lrat",    # activated HSC 低 Lrat
  Kupffer               = "Clec4f",  # resident Kupffer 高表达 Clec4f
  DC                    = "Flt3",    # precursor DC 高表达 Flt3
  Endothelial           = "Cldn5",   # 静息内皮高 Cldn5, 激活后改变
  Microglia_Homeostatic = "P2ry12",  # 稳态标志
  Microglia_DAM         = "P2ry12",  # 同上 (反向作为 DAM 起点)
  Astrocyte_Gfap        = "Slc1a2"   # 静息 astrocyte 高表达 Slc1a2
)

# 组颜色 (和其他 stages 一致)
GROUP_COLORS <- c("Tumor"="#E41A1C", "Blank"="#FF7F00", "Control"="#377EB8")

# 并发控制 (硬限制, 不按 CPU 核数算)
TASKS_CONCURRENT   <- 4       # 同时跑的任务数上限 (脚本是顺序执行, 此参数作为文档)
CORES_PER_TASK     <- 8       # tradeSeq fitGAM 的 BiocParallel worker 数
                              # 每个 worker fork 会复制父进程 ~50 GB, 8 × 50 = 400 GB 上限
                              # 从原 40 降到 8 是为了避免爆内存 (32 × 50 GB = 1.6 TB)
MAX_PARALLEL_CORES <- 128     # 全脚本 CPU 使用上限 (仅文档用)

# tradeSeq 开关
# 2026-04-12 经验教训: multi-condition fitGAM 在 50K+ 细胞 + 3 conditions × 3 lineages
#   场景下收敛极慢 (超过 2 小时仍在跑, TIME+ 一直增长但看不到结束点)
# 降级方案: 跳过 tradeSeq, 只保留 Spearman 相关 + lowess 分组动态图作为时间维度证据
# 如果未来数据更小或 tradeSeq 有性能改进, 可以重新开启
ENABLE_TRADESEQ <- FALSE

# 内存硬限制
MIN_FREE_GB_TO_START <- 200   # 启动单个任务前要求的最小可用内存
MIN_FREE_GB_ABORT    <- 100   # 运行中低于此立即 abort 整个脚本
MAX_R_PROCS_ALLOWED  <- 2     # 脚本启动时允许的最大残留 R 进程数 (自身 + 1)

cat("========================================\n")
cat("Pseudotime v4: core cell types + merged\n")
cat("Start:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================\n\n")

options(future.globals.maxSize = 256 * 1024^3)

# ============================================================
# 1.5 硬限制安全检查函数
# ============================================================

# 检查可用内存 (GB)
get_free_mem_gb <- function() {
  tryCatch({
    out <- system("free -g | awk 'NR==2 {print $7}'", intern = TRUE)
    as.numeric(out)
  }, error = function(e) NA)
}

# 检查 swap 使用 (GB)
get_swap_used_gb <- function() {
  tryCatch({
    out <- system("free -g | awk 'NR==3 {print $3}'", intern = TRUE)
    as.numeric(out)
  }, error = function(e) NA)
}

# 检查残留 R 进程数 (当前用户)
get_r_proc_count <- function() {
  tryCatch({
    out <- system("pgrep -u $USER -c '^(R|Rscript)$' 2>/dev/null", intern = TRUE)
    n <- suppressWarnings(as.integer(out))
    if (is.na(n)) 0 else n
  }, error = function(e) 0)
}

# 启动前检查 (失败就 stop)
preflight_check <- function() {
  cat("========== Pre-run safety check ==========\n")

  n_r <- get_r_proc_count()
  cat(sprintf("  Residual R processes: %d (limit: %d)\n",
              n_r, MAX_R_PROCS_ALLOWED))
  if (n_r > MAX_R_PROCS_ALLOWED) {
    stop(sprintf("ABORT: %d R 进程残留 (限制 %d)。请先: pkill -9 -u $USER R",
                 n_r, MAX_R_PROCS_ALLOWED))
  }

  free_gb <- get_free_mem_gb()
  cat(sprintf("  Available memory: %.0f GB (min required: %d GB)\n",
              if (is.na(free_gb)) 0 else free_gb, MIN_FREE_GB_TO_START))
  if (is.na(free_gb) || free_gb < MIN_FREE_GB_TO_START) {
    stop(sprintf("ABORT: 可用内存 %.0f GB < %d GB。请先清理进程或等待 swap 释放",
                 if (is.na(free_gb)) 0 else free_gb, MIN_FREE_GB_TO_START))
  }

  swap_used <- get_swap_used_gb()
  if (!is.na(swap_used) && swap_used > 1) {
    cat(sprintf("  ⚠ Swap in use: %.1f GB (服务器可能还在恢复, 不阻止运行)\n",
                swap_used))
  } else {
    cat(sprintf("  Swap used: %.1f GB (OK)\n",
                if (is.na(swap_used)) 0 else swap_used))
  }

  cat("  ✓ Safety check passed\n\n")
}

# 运行中检查 (启动新任务前调用, 低于阈值就 stop)
runtime_check <- function(task_idx, task_key) {
  free_gb <- get_free_mem_gb()
  if (is.na(free_gb)) {
    cat("  ⚠ runtime_check: 无法获取内存状态, skip\n")
    return(invisible(NULL))
  }
  cat(sprintf("  [runtime_check] task %d (%s) starting, free mem: %.0f GB\n",
              task_idx, task_key, free_gb))
  if (free_gb < MIN_FREE_GB_ABORT) {
    stop(sprintf("ABORT at task %d (%s): free mem %.0f GB < %d GB threshold",
                 task_idx, task_key, free_gb, MIN_FREE_GB_ABORT))
  }
  invisible(NULL)
}

# 立即做启动前检查
preflight_check()

# ============================================================
# 2. sample → group 的映射
# ============================================================
# 这里用 diff_common.R 定义的样本列表, 拆成三组
LIVER_SAMPLES_BY_GROUP <- list(
  Tumor   = c("liver3", "liver8", "gan56"),
  Blank   = c("bla_liver"),
  Control = c("con_liver")
)
BRAIN_SAMPLES_BY_GROUP <- list(
  Tumor   = c("brain3", "brain8", "nao56"),
  Blank   = c("bla_brain"),
  Control = c("con_brain")
)

get_sample_group <- function(sample, organ) {
  grp_map <- if (organ == "liver") LIVER_SAMPLES_BY_GROUP else BRAIN_SAMPLES_BY_GROUP
  for (g in names(grp_map)) {
    if (sample %in% grp_map[[g]]) return(g)
  }
  return(NA)
}

# ============================================================
# 3. 加载并合并 v4 RDS (肝 + 脑分开)
# ============================================================
cat("========== Loading v4 RDS and merging by organ ==========\n")

load_and_merge_organ <- function(organ) {
  samples <- unlist(
    if (organ == "liver") LIVER_SAMPLES_BY_GROUP else BRAIN_SAMPLES_BY_GROUP
  )
  cat(sprintf("[%s] loading %d samples\n", organ, length(samples)))

  objs <- list()
  for (s in samples) {
    fp <- resolve_rds_path(s, "v4")
    if (!file.exists(fp)) {
      write_error(s, paste0("missing v4 RDS: ", fp))
      next
    }
    o <- tryCatch({
      x <- readRDS(fp)
      if (inherits(x, "Seurat")) x <- tryCatch(JoinLayers(x), error = function(e) x)
      x <- RenameCells(x, add.cell.id = s)
      x$sample <- s
      x$group  <- get_sample_group(s, organ)
      x$organ  <- organ
      x
    }, error = function(e) {
      write_error(s, e$message); NULL
    })
    if (!is.null(o)) {
      objs[[s]] <- o
      cat(sprintf("  ✓ %-12s  cells=%6d  group=%s\n", s, ncol(o), o$group[1]))
    }
  }

  if (length(objs) < 2) {
    write_error(organ, "fewer than 2 samples loaded, cannot merge")
    return(NULL)
  }

  cat(sprintf("  merging %d %s samples...\n", length(objs), organ))
  merged <- merge(objs[[1]], objs[-1])
  merged <- tryCatch(JoinLayers(merged), error = function(e) merged)
  cat(sprintf("  [%s] merged cells: %d\n", organ, ncol(merged)))

  rm(objs); gc(verbose = FALSE)
  return(merged)
}

seurat_liver <- NULL  # 按需加载, 避免同时持有两个 organ 的合并对象
seurat_brain <- NULL

# ============================================================
# 4. Integrated UMAP (Harmony on organ-level merged object)
#    用于全景概览, 单个细类后续会独立重跑 PCA/UMAP
#    这一步在主循环里按需调用, 不在脚本加载阶段自动跑
# ============================================================
cat("\n========== Integrated UMAP (will be built on demand) ==========\n")

build_integrated_umap <- function(obj, organ) {
  if (is.null(obj)) return(NULL)
  cat(sprintf("[%s overview] Normalize → VariableFeatures → Scale → PCA → Harmony → UMAP\n",
              organ))
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = 50, verbose = FALSE)
  obj <- tryCatch(
    harmony::RunHarmony(obj, group.by.vars = "sample", verbose = FALSE),
    error = function(e) {
      write_error(paste0(organ, "_harmony"), e$message)
      obj  # 回退到 PCA
    }
  )
  red_for_umap <- if ("harmony" %in% names(obj@reductions)) "harmony" else "pca"
  obj <- RunUMAP(obj, reduction = red_for_umap, dims = 1:UMAP_DIMS, verbose = FALSE)

  # 保存整合后 UMAP 概览图
  # 兼容不同命名: celltype_refined (v4 常用) 或 v4_refined 或 final_celltype
  ct_col <- NULL
  for (candidate in c("celltype_refined", "v4_refined", "final_celltype")) {
    if (candidate %in% colnames(obj@meta.data)) { ct_col <- candidate; break }
  }
  if (is.null(ct_col)) stop("没找到细胞类型列 (celltype_refined/v4_refined/final_celltype)")
  cat(sprintf("  using cell type column: %s\n", ct_col))

  df <- data.frame(
    UMAP1 = Embeddings(obj, "umap")[, 1],
    UMAP2 = Embeddings(obj, "umap")[, 2],
    sample = obj$sample,
    group  = obj$group,
    celltype = obj@meta.data[[ct_col]]
  )
  p1 <- ggplot(df, aes(UMAP1, UMAP2, color = group)) +
    geom_point(size = 0.15, alpha = 0.6) +
    scale_color_manual(values = GROUP_COLORS) +
    labs(title = paste0(organ, " integrated UMAP (by group)")) +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5))
  p2 <- ggplot(df, aes(UMAP1, UMAP2, color = celltype)) +
    geom_point(size = 0.15, alpha = 0.6) +
    labs(title = paste0(organ, " integrated UMAP (by v4 cell type)")) +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                       legend.text = element_text(size = 7))
  base <- file.path(OUT_DIR, "integrated", paste0(organ, "_integrated_umap"))
  ggsave(paste0(base, "_bygroup.pdf"),    p1, width = 9, height = 7, dpi = 300)
  ggsave(paste0(base, "_bygroup.png"),    p1, width = 9, height = 7, dpi = 200)
  ggsave(paste0(base, "_bycelltype.pdf"), p2, width = 11, height = 7, dpi = 300)
  ggsave(paste0(base, "_bycelltype.png"), p2, width = 11, height = 7, dpi = 200)

  return(obj)
}

# 注意: 不在这里立即调用 build_integrated_umap
# 而是在主循环里按 organ 分段运行

# ============================================================
# 5. 辅助: KNN 填充 pseudotime 的 NA
# ============================================================
fill_pt_na_knn <- function(umap_mat, pseudotime_vec, k = KNN_K) {
  na_idx <- which(is.na(pseudotime_vec))
  nona_idx <- which(!is.na(pseudotime_vec))
  if (length(na_idx) == 0) return(pseudotime_vec)
  if (length(nona_idx) < k) k <- max(1, floor(length(nona_idx) / 2))

  knn_res <- FNN::get.knnx(
    data  = umap_mat[nona_idx, , drop = FALSE],
    query = umap_mat[na_idx,   , drop = FALSE],
    k     = min(k, length(nona_idx))
  )
  filled <- numeric(length(na_idx))
  for (i in seq_along(na_idx)) {
    neighbors <- knn_res$nn.index[i, ]
    dists     <- knn_res$nn.dist[i, ]
    if (all(dists == 0)) dists <- rep(1, length(dists))
    else dists[dists == 0] <- min(dists[dists > 0]) / 10
    w <- 1 / dists; w <- w / sum(w)
    filled[i] <- sum(pseudotime_vec[nona_idx[neighbors]] * w)
  }
  pseudotime_vec[na_idx] <- filled
  pseudotime_vec
}

# ============================================================
# 6. 核心: 单细类拟时序分析流程
# ============================================================
run_one_cell_type <- function(parent_obj, organ, virtual_ct) {
  key <- paste0(organ, "_", virtual_ct)
  ct_dir <- file.path(OUT_DIR, organ, virtual_ct)
  dir.create(ct_dir, recursive = TRUE, showWarnings = FALSE)
  ts_dir <- file.path(ct_dir, "tradeseq_results")
  dir.create(ts_dir, recursive = TRUE, showWarnings = FALSE)

  result <- list(
    organ = organ, virtual_ct = virtual_ct, key = key,
    status = "error", n_cells = NA, n_curves = NA,
    n_spearman_sig = NA, n_associated_gam = NA, n_pattern_sig = NA,
    elapsed_min = NA
  )
  t_begin <- proc.time()

  tryCatch({
    # 1. 抽出该细类的细胞 (extract_cells_v2 内部处理列名)
    cell_type_map <- if (organ == "liver") LIVER_CELL_TYPES else BRAIN_CELL_TYPES
    cells <- extract_cells_v2(parent_obj, virtual_ct, "v4_refined", cell_type_map)
    if (length(cells) < MIN_CELLS_FOR_RUN) {
      result$status <- "insufficient_cells"
      result$n_cells <- length(cells)
      return(result)
    }

    sub <- subset(parent_obj, cells = cells)
    sub <- tryCatch(JoinLayers(sub), error = function(e) sub)
    result$n_cells <- ncol(sub)
    cat(sprintf("\n[%s] n_cells=%d, groups: %s\n",
                key, ncol(sub),
                paste0(names(table(sub$group)), "=",
                       as.vector(table(sub$group)), collapse = ", ")))

    # 2. 子集内部重新做 PCA / Harmony / UMAP
    #    (因为亚群细胞组成和整合 umap 不完全一样)
    sub <- NormalizeData(sub, verbose = FALSE)
    sub <- FindVariableFeatures(sub, nfeatures = 2000, verbose = FALSE)
    sub <- ScaleData(sub, verbose = FALSE)
    sub <- RunPCA(sub, npcs = 30, verbose = FALSE)
    sub <- tryCatch(
      harmony::RunHarmony(sub, group.by.vars = "sample", verbose = FALSE),
      error = function(e) { write_error(paste0(key, "_harmony"), e$message); sub }
    )
    red_use <- if ("harmony" %in% names(sub@reductions)) "harmony" else "pca"
    sub <- RunUMAP(sub, reduction = red_use, dims = 1:min(20, UMAP_DIMS),
                   verbose = FALSE)

    umap_mat <- Embeddings(sub, "umap")
    colnames(umap_mat) <- c("UMAP1", "UMAP2")

    # 3. kmeans clustering for slingshot
    set.seed(42)
    n_clust <- min(KMEANS_K, max(3, floor(ncol(sub) / 200)))
    kc <- kmeans(umap_mat, centers = n_clust, iter.max = 50)
    sub$sling_cluster <- as.factor(kc$cluster)

    # 4. 转 SCE
    counts <- GetAssayData(sub, assay = "RNA", layer = "counts")
    if (!inherits(counts, "dgCMatrix")) counts <- as(counts, "dgCMatrix")
    sce <- SingleCellExperiment(
      assays  = list(counts = counts),
      colData = sub@meta.data
    )
    reducedDim(sce, "UMAP") <- umap_mat

    # 5. 自动找 root cluster (基于 root marker 基因表达最高的 cluster)
    root_marker <- ROOT_MARKERS[[virtual_ct]]
    start_cluster <- NULL
    if (!is.null(root_marker) && root_marker %in% rownames(counts)) {
      marker_expr <- as.numeric(counts[root_marker, ])
      cluster_mean <- tapply(marker_expr, sub$sling_cluster, mean)
      start_cluster <- as.character(which.max(cluster_mean))
      cat(sprintf("  [%s] root marker=%s, start_cluster=%s (mean_expr=%.3f)\n",
                  key, root_marker, start_cluster, max(cluster_mean)))
    } else {
      cat(sprintf("  [%s] no root marker or not expressed, auto-detect\n", key))
    }

    # 6. 运行 slingshot
    sce <- slingshot::slingshot(
      sce,
      clusterLabels = sub$sling_cluster,
      reducedDim    = "UMAP",
      start.clus    = start_cluster,
      allow.breaks  = TRUE
    )

    pt_mat <- slingPseudotime(sce)
    cw_mat <- slingCurveWeights(sce)
    if (is.null(pt_mat)) stop("slingshot 未生成 pseudotime")

    n_curves <- ncol(pt_mat)
    result$n_curves <- n_curves
    cat(sprintf("  [%s] detected %d curves\n", key, n_curves))

    # 7. NA 填充 (使用 KNN)
    pt_main <- pt_mat[, 1]
    na_pct <- 100 * sum(is.na(pt_main)) / length(pt_main)
    cat(sprintf("  [%s] main curve NA: %.1f%%\n", key, na_pct))
    if (na_pct > 0) {
      pt_main <- fill_pt_na_knn(umap_mat, pt_main, k = KNN_K)
    }
    sub$pseudotime <- pt_main
    sce$pseudotime <- pt_main

    # 8. 可视化: trajectory
    # 探测 sub 里的细胞类型列 (用于 df_traj)
    ct_col_sub <- NULL
    for (candidate in c("celltype_refined", "v4_refined", "final_celltype")) {
      if (candidate %in% colnames(sub@meta.data)) { ct_col_sub <- candidate; break }
    }

    df_traj <- data.frame(
      UMAP1      = umap_mat[, 1],
      UMAP2      = umap_mat[, 2],
      pseudotime = pt_main,
      group      = sub$group,
      sample     = sub$sample,
      celltype   = if (!is.null(ct_col_sub)) as.character(sub@meta.data[[ct_col_sub]])
                   else rep(virtual_ct, ncol(sub)),
      sling_cluster = as.character(sub$sling_cluster)
    )

    p_traj <- ggplot(df_traj, aes(UMAP1, UMAP2, color = pseudotime)) +
      geom_point(size = 0.5, alpha = 0.8) +
      scale_color_viridis_c(option = "plasma") +
      labs(title = sprintf("%s · %s · pseudotime", organ, virtual_ct),
           subtitle = sprintf("n_cells=%d, n_curves=%d", ncol(sub), n_curves)) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5))
    # 叠加曲线
    curves <- slingshot::slingCurves(sce)
    for (cv in curves) {
      cdf <- as.data.frame(cv$s[cv$tag, ])
      colnames(cdf) <- c("UMAP1", "UMAP2")
      p_traj <- p_traj + geom_path(data = cdf,
                                    aes(UMAP1, UMAP2),
                                    color = "black", linewidth = 0.8,
                                    inherit.aes = FALSE)
    }
    ggsave(file.path(ct_dir, "trajectory.pdf"), p_traj, width = 9, height = 7, dpi = 300)
    ggsave(file.path(ct_dir, "trajectory.png"), p_traj, width = 9, height = 7, dpi = 200)

    # 可视化: 按 group 的 pseudotime 分布 (关键图: 合并分析的独特价值)
    p_by_group <- ggplot(df_traj, aes(x = pseudotime, fill = group)) +
      geom_density(alpha = 0.55) +
      scale_fill_manual(values = GROUP_COLORS) +
      labs(title = sprintf("%s · %s · pseudotime by group", organ, virtual_ct),
           subtitle = "同一轨迹上 tumor vs blank vs control 细胞的分布") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5))
    ggsave(file.path(ct_dir, "pseudotime_by_group.pdf"),
           p_by_group, width = 9, height = 6, dpi = 300)
    ggsave(file.path(ct_dir, "pseudotime_by_group.png"),
           p_by_group, width = 9, height = 6, dpi = 200)

    # 可视化: 曲线分配
    df_traj$curve_assign <- as.character(apply(cw_mat, 1, which.max))
    p_curve <- ggplot(df_traj, aes(UMAP1, UMAP2, color = curve_assign)) +
      geom_point(size = 0.5, alpha = 0.8) +
      scale_color_brewer(palette = "Set1") +
      labs(title = sprintf("%s · %s · curve assignment", organ, virtual_ct)) +
      theme_bw() + theme(plot.title = element_text(hjust = 0.5))
    ggsave(file.path(ct_dir, "curve_assignment.pdf"),
           p_curve, width = 9, height = 7, dpi = 300)
    ggsave(file.path(ct_dir, "curve_assignment.png"),
           p_curve, width = 9, height = 7, dpi = 200)

    # 导出 pseudotime 值
    df_pt_out <- df_traj
    for (i in seq_len(n_curves)) {
      df_pt_out[[paste0("curve_", i, "_pt")]] <- pt_mat[, i]
    }
    write.csv(df_pt_out, file.path(ct_dir, "pseudotime_values.csv"),
              row.names = FALSE)

    # 9. Spearman 全基因组相关
    cat(sprintf("  [%s] Spearman correlation on main curve...\n", key))
    log_counts <- log1p(as.matrix(counts))
    min_cells_expr <- max(50, ncol(sub) * 0.05)
    genes_keep <- rowSums(counts > 0) >= min_cells_expr
    log_counts_f <- log_counts[genes_keep, ]
    cat(sprintf("  [%s] Spearman on %d genes × %d cells\n",
                key, nrow(log_counts_f), ncol(log_counts_f)))

    cor_df <- data.frame(
      gene = rownames(log_counts_f),
      rho  = NA_real_,
      pval = NA_real_,
      stringsAsFactors = FALSE
    )
    for (i in seq_len(nrow(log_counts_f))) {
      tryCatch({
        ct <- cor.test(log_counts_f[i, ], pt_main, method = "spearman",
                       exact = FALSE)
        cor_df$rho[i]  <- ct$estimate
        cor_df$pval[i] <- ct$p.value
      }, error = function(e) NULL)
    }
    cor_df$padj <- p.adjust(cor_df$pval, method = "BH")
    cor_df <- cor_df[order(abs(cor_df$rho), decreasing = TRUE), ]
    write.csv(cor_df, file.path(ct_dir, "spearman_correlation.csv"),
              row.names = FALSE)

    n_sig <- sum(cor_df$padj < 0.05 & abs(cor_df$rho) > 0.3, na.rm = TRUE)
    result$n_spearman_sig <- n_sig
    cat(sprintf("  [%s] Spearman sig (padj<0.05, |rho|>0.3): %d\n", key, n_sig))

    # 10. Top gene 热图
    top_genes_heatmap <- head(cor_df$gene[cor_df$padj < 0.05 &
                                          abs(cor_df$rho) > 0.3], 80)
    if (length(top_genes_heatmap) >= 10) {
      n_cells_hm <- min(2000, ncol(sub))
      set.seed(42)
      cell_idx_hm <- sample(ncol(sub), n_cells_hm)
      pt_hm <- pt_main[cell_idx_hm]
      ord_hm <- order(pt_hm)
      mat_hm <- log_counts_f[top_genes_heatmap, cell_idx_hm][, ord_hm]
      mat_hm <- t(scale(t(mat_hm)))
      mat_hm[mat_hm > 2] <- 2
      mat_hm[mat_hm < -2] <- -2
      mat_hm <- mat_hm[complete.cases(mat_hm), ]
      if (nrow(mat_hm) >= 5) {
        col_anno <- data.frame(
          pseudotime = pt_hm[ord_hm],
          group      = df_traj$group[cell_idx_hm][ord_hm]
        )
        rownames(col_anno) <- colnames(mat_hm)
        tryCatch({
          pheatmap(mat_hm,
                   cluster_rows = TRUE, cluster_cols = FALSE,
                   annotation_col = col_anno,
                   show_colnames = FALSE,
                   main = sprintf("%s · %s · top %d Spearman genes",
                                  organ, virtual_ct, nrow(mat_hm)),
                   filename = file.path(ct_dir, "heatmap_top_genes.pdf"),
                   width = 11, height = 12)
        }, error = function(e) write_error(paste0(key, "_heatmap"), e$message))
      }
    }

    # 11. Core panel dynamic genes plot (lowess smoothing by group)
    core_in_data <- intersect(CORE_PANEL, rownames(log_counts_f))
    cat(sprintf("  [%s] core panel genes in data: %d / %d\n",
                key, length(core_in_data), length(CORE_PANEL)))
    if (length(core_in_data) >= 3) {
      panel_dfs <- list()
      for (g in core_in_data) {
        expr <- log_counts_f[g, ]
        for (grp in unique(df_traj$group)) {
          idx <- which(df_traj$group == grp)
          if (length(idx) < 20) next
          ord <- idx[order(pt_main[idx])]
          sm <- tryCatch(
            stats::lowess(pt_main[ord], expr[ord], f = 0.2),
            error = function(e) NULL
          )
          if (!is.null(sm)) {
            panel_dfs[[length(panel_dfs) + 1]] <- data.frame(
              gene = g, group = grp,
              pseudotime = sm$x, expression = sm$y
            )
          }
        }
      }
      if (length(panel_dfs) > 0) {
        pdf_panel <- dplyr::bind_rows(panel_dfs)
        p_dyn <- ggplot(pdf_panel,
                        aes(pseudotime, expression, color = group)) +
          geom_line(linewidth = 0.8) +
          facet_wrap(~ gene, scales = "free_y", ncol = 4) +
          scale_color_manual(values = GROUP_COLORS) +
          labs(title = sprintf("%s · %s · core panel dynamics",
                               organ, virtual_ct),
               x = "Pseudotime", y = "Log1p(Expression) lowess") +
          theme_bw(base_size = 9) +
          theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                strip.text = element_text(size = 7),
                legend.position = "bottom")
        n_genes_pdyn <- length(unique(pdf_panel$gene))
        ncol_dyn <- 4
        nrow_dyn <- ceiling(n_genes_pdyn / ncol_dyn)
        ggsave(file.path(ct_dir, "dynamic_core_genes.pdf"),
               p_dyn,
               width = ncol_dyn * 3.2,
               height = nrow_dyn * 2.0, dpi = 300, limitsize = FALSE)
        ggsave(file.path(ct_dir, "dynamic_core_genes.png"),
               p_dyn,
               width = ncol_dyn * 3.2,
               height = nrow_dyn * 2.0, dpi = 200, limitsize = FALSE)
      }
    }

    # 12. tradeSeq: 选基因 (可通过 ENABLE_TRADESEQ 跳过)
    if (!ENABLE_TRADESEQ) {
      cat(sprintf("  [%s] SKIP tradeSeq (ENABLE_TRADESEQ=FALSE)\n", key))
      result$status <- "ok_no_gam"
      result$elapsed_min <- round((proc.time() - t_begin)["elapsed"] / 60, 1)
      cat(sprintf("  ✓ [%s] DONE  n=%d  curves=%d  sp_sig=%s  (%.1f min, no fitGAM)\n",
                  key, result$n_cells, result$n_curves,
                  result$n_spearman_sig, result$elapsed_min))
      return(result)
    }

    top_spearman_genes <- head(
      cor_df$gene[cor_df$padj < 0.05 & abs(cor_df$rho) > 0.25],
      TOP_GENES_FOR_GAM
    )
    genes_for_gam <- unique(c(top_spearman_genes, core_in_data))
    genes_for_gam <- intersect(genes_for_gam, rownames(counts))
    cat(sprintf("  [%s] genes for tradeSeq fitGAM: %d\n",
                key, length(genes_for_gam)))

    if (length(genes_for_gam) < 20) {
      cat(sprintf("  [%s] too few genes for tradeSeq, skip\n", key))
      result$status <- "ok_no_gam"
      result$elapsed_min <- round((proc.time() - t_begin)["elapsed"] / 60, 1)
      return(result)
    }

    # 13. 准备 tradeSeq 输入
    sce_for_gam <- sce[genes_for_gam, ]
    cw <- slingshot::slingCurveWeights(sce)
    pt_for_gam <- slingshot::slingPseudotime(sce)
    # 填充 pt NA
    for (ci in seq_len(ncol(pt_for_gam))) {
      pt_for_gam[, ci] <- fill_pt_na_knn(umap_mat, pt_for_gam[, ci], k = KNN_K)
    }

    # 构造 group conditions (tumor/blank/control)
    cond_vec <- factor(sub$group, levels = c("Control", "Blank", "Tumor"))

    # 14. fitGAM 前的内存释放 (关键!)
    # 在 fork 8 worker 之前, 尽可能 rm 掉父进程的大对象
    # 只保留 fitGAM 必需的: counts_for_gam, pt_for_gam, cw, cond_vec, genes_for_gam
    counts_for_gam <- counts[genes_for_gam, , drop = FALSE]
    # sce_for_gam 是之前没用的子集, counts / log_counts / log_counts_f / sub / sce 都可以释放
    rm(log_counts, log_counts_f, sub, sce, sce_for_gam, counts)
    if (exists("expr_mat")) rm(expr_mat)
    if (exists("expr_mat_clean")) rm(expr_mat_clean)
    if (exists("mat_hm")) rm(mat_hm)
    if (exists("cor_df")) rm(cor_df)
    for (j in 1:3) gc(verbose = FALSE)

    free_before_gam <- get_free_mem_gb()
    cat(sprintf("  [%s] free mem before fitGAM: %.0f GB\n",
                key, if (is.na(free_before_gam)) -1 else free_before_gam))

    # 15. fitGAM — 使用 SerialParam (完全不 fork worker)
    # 原因: 并行 fork 会让每个 worker 复制父进程大对象 (seurat_liver/brain 等)
    # 导致 N × 大对象 内存爆炸. 串行反而更快, 因为 mgcv GAM 拟合本身很快
    # fork 开销占比高, SerialParam 避免所有 fork 成本
    cat(sprintf("  [%s] fitGAM (%d knots, SerialParam - no fork)...\n",
                key, N_KNOTS_GAM))
    bpp <- BiocParallel::SerialParam()
    gam_fit <- tryCatch(
      tradeSeq::fitGAM(
        counts    = counts_for_gam,
        pseudotime = pt_for_gam,
        cellWeights = cw,
        conditions = cond_vec,
        nknots    = N_KNOTS_GAM,
        parallel  = FALSE,
        BPPARAM   = bpp,
        verbose   = FALSE
      ),
      error = function(e) {
        write_error(paste0(key, "_fitGAM"), e$message); NULL
      }
    )

    if (is.null(gam_fit)) {
      result$status <- "ok_gam_failed"
      result$elapsed_min <- round((proc.time() - t_begin)["elapsed"] / 60, 1)
      return(result)
    }

    saveRDS(gam_fit, file.path(ts_dir, "core_panel_fit.rds"))
    rm(counts_for_gam); gc(verbose = FALSE)

    # 16. associationTest: 基因是否沿伪时间变化
    assoc <- tryCatch(
      tradeSeq::associationTest(gam_fit, global = TRUE, lineages = FALSE),
      error = function(e) { write_error(paste0(key, "_assoc"), e$message); NULL }
    )
    if (!is.null(assoc)) {
      assoc$gene <- rownames(assoc)
      assoc$padj <- p.adjust(assoc$pvalue, method = "BH")
      assoc <- assoc[order(assoc$padj), ]
      write.csv(assoc,
                file.path(ts_dir, "association_test.csv"),
                row.names = FALSE)
      result$n_associated_gam <- sum(assoc$padj < 0.05, na.rm = TRUE)
    }

    # 16. patternTest: 条件间的轨迹形状差异 (关键!)
    ptn <- tryCatch(
      tradeSeq::conditionTest(gam_fit, global = TRUE, pairwise = TRUE),
      error = function(e) { write_error(paste0(key, "_pattern"), e$message); NULL }
    )
    if (!is.null(ptn)) {
      ptn$gene <- rownames(ptn)
      ptn$padj <- p.adjust(ptn$pvalue, method = "BH")
      ptn <- ptn[order(ptn$padj), ]
      write.csv(ptn,
                file.path(ts_dir, "pattern_test.csv"),
                row.names = FALSE)
      result$n_pattern_sig <- sum(ptn$padj < 0.05, na.rm = TRUE)
    }

    result$status <- "ok"
    result$elapsed_min <- round((proc.time() - t_begin)["elapsed"] / 60, 1)
    cat(sprintf("  ✓ [%s] DONE  n=%d  curves=%d  sp_sig=%s  assoc=%s  pattern=%s  (%.1f min)\n",
                key, result$n_cells, result$n_curves,
                result$n_spearman_sig, result$n_associated_gam,
                result$n_pattern_sig, result$elapsed_min))

  }, error = function(e) {
    result$status <- "error"
    result$elapsed_min <- round((proc.time() - t_begin)["elapsed"] / 60, 1)
    write_error(key, e$message)
  })

  gc(verbose = FALSE)
  return(result)
}

# ============================================================
# 7. 跑所有 10 个核心细类
# ============================================================
cat("\n========== Running pseudotime analysis (10 core cell types) ==========\n")

all_tasks <- list()
for (ct in CELL_TYPES_TO_RUN$liver) {
  all_tasks[[length(all_tasks) + 1]] <- list(organ = "liver", virtual_ct = ct)
}
for (ct in CELL_TYPES_TO_RUN$brain) {
  all_tasks[[length(all_tasks) + 1]] <- list(organ = "brain", virtual_ct = ct)
}
cat(sprintf("Total tasks: %d (%d liver + %d brain)\n\n",
            length(all_tasks),
            length(CELL_TYPES_TO_RUN$liver),
            length(CELL_TYPES_TO_RUN$brain)))

# 按器官分段执行: 先跑 liver 全部任务, 释放 seurat_liver, 再跑 brain
# 同时在每个 organ 段内, 先加载对象 + 建整合 UMAP, 然后跑该 organ 的全部任务
t_all_start <- proc.time()
task_results <- list()

process_organ_tasks <- function(organ_name, task_indices_in_all) {
  cat(sprintf("\n\n========== Processing all %s tasks ==========\n",
              toupper(organ_name)))

  # 1. 加载该器官的合并对象
  cat(sprintf("\n-- Loading and merging %s samples --\n", organ_name))
  organ_obj <- load_and_merge_organ(organ_name)
  if (is.null(organ_obj)) {
    write_error(organ_name, "load_and_merge_organ returned NULL")
    return(invisible(NULL))
  }

  free_after_load <- get_free_mem_gb()
  cat(sprintf("  [%s loaded] free mem: %.0f GB\n",
              organ_name,
              if (is.na(free_after_load)) -1 else free_after_load))

  # 2. 建整合 UMAP (需要父对象)
  cat(sprintf("\n-- Building %s integrated UMAP --\n", organ_name))
  organ_obj <- build_integrated_umap(organ_obj, organ_name)

  free_after_umap <- get_free_mem_gb()
  cat(sprintf("  [%s umap done] free mem: %.0f GB\n",
              organ_name,
              if (is.na(free_after_umap)) -1 else free_after_umap))

  # 3. 跑该器官的全部 task
  for (i in task_indices_in_all) {
    task <- all_tasks[[i]]
    task_key <- paste0(task$organ, "_", task$virtual_ct)

    cat(sprintf("\n### Task %d/%d: %s · %s ###\n",
                i, length(all_tasks), task$organ, task$virtual_ct))

    runtime_check(i, task_key)

    task_results[[i]] <<- run_one_cell_type(organ_obj, task$organ, task$virtual_ct)

    for (j in 1:3) gc(verbose = FALSE)

    free_gb_after <- get_free_mem_gb()
    if (!is.na(free_gb_after)) {
      cat(sprintf("  [after task %d] free mem: %.0f GB\n", i, free_gb_after))
    }
  }

  # 4. 释放该器官的对象 (关键! 避免下一个 organ 时同时持有两个)
  cat(sprintf("\n-- Releasing %s merged object --\n", organ_name))
  rm(organ_obj)
  for (j in 1:5) gc(verbose = FALSE)  # 5 次 gc 充分释放

  free_after_release <- get_free_mem_gb()
  cat(sprintf("  [%s released] free mem: %.0f GB\n",
              organ_name,
              if (is.na(free_after_release)) -1 else free_after_release))
}

# 分离 liver 任务和 brain 任务的 index
liver_task_idx <- which(sapply(all_tasks, function(t) t$organ == "liver"))
brain_task_idx <- which(sapply(all_tasks, function(t) t$organ == "brain"))

# 跑 liver 段
process_organ_tasks("liver", liver_task_idx)

# 跑 brain 段
process_organ_tasks("brain", brain_task_idx)

t_all_elapsed <- (proc.time() - t_all_start)["elapsed"]
cat(sprintf("\n\nAll tasks elapsed: %.1f min\n", t_all_elapsed / 60))

# ============================================================
# 8. 汇总 task_status.csv
# ============================================================
cat("\n========== Writing task_status.csv ==========\n")
ts_df <- dplyr::bind_rows(lapply(task_results, function(r) {
  data.frame(
    key              = r$key,
    organ            = r$organ,
    virtual_ct       = r$virtual_ct,
    status           = r$status,
    n_cells          = if (is.null(r$n_cells)) NA else r$n_cells,
    n_curves         = if (is.null(r$n_curves)) NA else r$n_curves,
    n_spearman_sig   = if (is.null(r$n_spearman_sig)) NA else r$n_spearman_sig,
    n_associated_gam = if (is.null(r$n_associated_gam)) NA else r$n_associated_gam,
    n_pattern_sig    = if (is.null(r$n_pattern_sig)) NA else r$n_pattern_sig,
    elapsed_min      = if (is.null(r$elapsed_min)) NA else r$elapsed_min,
    stringsAsFactors = FALSE
  )
}))
write.csv(ts_df, file.path(OUT_DIR, "summary", "task_status.csv"),
          row.names = FALSE)
print(ts_df)

# ============================================================
# 9. 跨细类对比: 核心基因 pattern_test 显著性汇总
# ============================================================
cat("\n========== Aggregating pattern_test results across cell types ==========\n")
pattern_agg <- list()
for (task in all_tasks) {
  key <- paste0(task$organ, "_", task$virtual_ct)
  fp <- file.path(OUT_DIR, task$organ, task$virtual_ct,
                  "tradeseq_results", "pattern_test.csv")
  if (!file.exists(fp)) next
  df <- tryCatch(read.csv(fp, stringsAsFactors = FALSE),
                 error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) next
  df <- df[df$gene %in% CORE_PANEL, , drop = FALSE]
  if (nrow(df) == 0) next
  df$key <- key
  df$organ <- task$organ
  df$virtual_ct <- task$virtual_ct
  pattern_agg[[length(pattern_agg) + 1]] <- df
}

if (length(pattern_agg) > 0) {
  pattern_all <- dplyr::bind_rows(pattern_agg)
  pattern_all <- pattern_all %>%
    dplyr::select(dplyr::any_of(c("gene", "organ", "virtual_ct", "key",
                                  "waldStat", "df", "pvalue", "padj")))
  write.csv(pattern_all,
            file.path(OUT_DIR, "summary", "core_panel_pattern_test_all.csv"),
            row.names = FALSE)
  cat(sprintf("  core panel pattern test rows: %d\n", nrow(pattern_all)))
} else {
  cat("  no pattern_test results aggregated\n")
}

# ============================================================
# 10. 写 summary md
# ============================================================
cat("\n========== Writing pseudotime_summary.md ==========\n")

summary_lines <- c(
  "# Pseudotime v4: Core cell type analysis summary",
  "",
  sprintf("**Generated:** %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  sprintf("**Total elapsed:** %.1f min", t_all_elapsed / 60),
  "",
  "## Method overview",
  "",
  "- **v4 annotation** (`v4_refined`), 复用 Stage 1-3 的细类定义",
  "- **合并分析策略**：每个器官 (liver, brain) 合并所有样本 (4 liver, 4 brain)",
  "- **Harmony batch integration** (by `sample`)",
  "- **Per cell type**: subset → re-PCA → re-Harmony → re-UMAP → kmeans clustering → slingshot",
  "- **Root selection**: 基于 marker gene expression (各细类不同)",
  "- **NA handling**: KNN (k=10) 邻居加权填充",
  "- **Differential analysis**:",
  "  - Spearman correlation (全基因组, fast pass)",
  "  - tradeSeq fitGAM (top 500 Spearman + core panel ~80)",
  "  - tradeSeq associationTest (global, 检测沿伪时间变化)",
  "  - tradeSeq conditionTest (tumor vs control vs blank 轨迹形状差异)",
  "",
  "## Task status",
  "",
  "| Key | n_cells | n_curves | Spearman sig | GAM assoc | Pattern sig | Elapsed (min) |",
  "|---|---:|---:|---:|---:|---:|---:|"
)
for (i in seq_len(nrow(ts_df))) {
  r <- ts_df[i, ]
  summary_lines <- c(summary_lines,
    sprintf("| %s | %s | %s | %s | %s | %s | %s |",
            r$key,
            if (is.na(r$n_cells))          "-" else as.character(r$n_cells),
            if (is.na(r$n_curves))         "-" else as.character(r$n_curves),
            if (is.na(r$n_spearman_sig))   "-" else as.character(r$n_spearman_sig),
            if (is.na(r$n_associated_gam)) "-" else as.character(r$n_associated_gam),
            if (is.na(r$n_pattern_sig))    "-" else as.character(r$n_pattern_sig),
            if (is.na(r$elapsed_min))      "-" else sprintf("%.1f", r$elapsed_min)))
}

summary_lines <- c(summary_lines, "",
  "## Files",
  "",
  "- `integrated/{liver,brain}_integrated_umap_bygroup.{pdf,png}` — 两大整合对象概览",
  "- `{liver,brain}/<ct>/trajectory.{pdf,png}` — 每个细类的轨迹图",
  "- `{liver,brain}/<ct>/pseudotime_by_group.{pdf,png}` — tumor/blank/control 在同一轨迹上的分布 **(关键)**",
  "- `{liver,brain}/<ct>/spearman_correlation.csv` — 全基因 Spearman 结果",
  "- `{liver,brain}/<ct>/heatmap_top_genes.pdf` — top 基因热图",
  "- `{liver,brain}/<ct>/dynamic_core_genes.{pdf,png}` — 核心基因沿伪时间分组曲线 **(关键)**",
  "- `{liver,brain}/<ct>/tradeseq_results/association_test.csv` — 沿伪时间变化的基因",
  "- `{liver,brain}/<ct>/tradeseq_results/pattern_test.csv` — tumor vs control 轨迹差异基因 **(关键)**",
  "- `summary/core_panel_pattern_test_all.csv` — 核心基因在所有细类的 pattern test 汇总",
  "- `summary/task_status.csv` — 每个任务的状态和统计",
  "- `summary/error_log.txt` — 错误日志"
)

writeLines(summary_lines, file.path(OUT_DIR, "summary", "pseudotime_summary.md"))

cat("\n========================================\n")
cat("Pseudotime v4 analysis complete\n")
cat(sprintf("  Output dir : %s\n", OUT_DIR))
cat(sprintf("  Summary    : %s\n", file.path(OUT_DIR, "summary", "pseudotime_summary.md")))
cat(sprintf("  Elapsed    : %.1f min\n", t_all_elapsed / 60))
cat("========================================\n")
