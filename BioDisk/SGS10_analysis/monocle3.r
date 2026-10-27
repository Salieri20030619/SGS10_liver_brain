# ==================== 拟时序分析 Monocle3（项目定制版）====================
# 项目：26年第82个项目 / 26.2.25第2个项目
# 用法：conda activate R_4.5.2 && Rscript pseudotime_monocle3.R
#
# 【样本结构】
#   肝脏：gan56 / liver3 / liver8 / bla_liver / con_liver
#   脑部：nao56 / brain3 / brain8 / bla_brain / con_brain
#   每个样本路径：<样本名>/<样本名>_manual_annotated.rds
#
# 【分析策略】
#   每个器官把所有样本合并后统一跑一次 Monocle3
#   保留 sample / group / final_celltype 信息，用于着色和对比
# =============================================================================

cat("========================================\n")
cat("拟时序分析 Monocle3（项目定制版）\n")
cat("========================================\n\n")

# ==================== 0. 环境 ====================
.libPaths("/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/lib/R/library")

suppressPackageStartupMessages({
  library(igraph)
  library(Seurat)
  library(monocle3)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(RColorBrewer)
  library(viridis)
})

select    <- dplyr::select
filter    <- dplyr::filter
mutate    <- dplyr::mutate
arrange   <- dplyr::arrange
group_by  <- dplyr::group_by
summarise <- dplyr::summarise
bind_rows <- dplyr::bind_rows

cat("✓ 包加载完成\n\n")

# ==================== 1. 参数（按需修改）====================
setwd("/mnt/BioDisk/wangmingyuan/SGS10_analysis")

# ---- 样本定义 ----
LIVER_SAMPLES <- list(
  list(name="gan56",    group="Tumor",   organ="Liver"),
  list(name="liver3",   group="Tumor",   organ="Liver"),
  list(name="liver8",   group="Tumor",   organ="Liver"),
  list(name="bla_liver",group="Blank",   organ="Liver"),
  list(name="con_liver",group="Control", organ="Liver")
)

BRAIN_SAMPLES <- list(
  list(name="nao56",    group="Tumor",   organ="Brain"),
  list(name="brain3",   group="Tumor",   organ="Brain"),
  list(name="brain8",   group="Tumor",   organ="Brain"),
  list(name="bla_brain",group="Blank",   organ="Brain"),
  list(name="con_brain",group="Control", organ="Brain")
)

# ---- 列名（与注释脚本保持一致）----
CELLTYPE_COL <- "final_celltype"

# ---- 拟时序根细胞类型（生物学先验）----
# 肝脏：从正常肝细胞出发，沿肿瘤化方向
LIVER_ROOT_CELLTYPE <- "Hepatocytes"
# 脑部：从稳态小胶质细胞出发
BRAIN_ROOT_CELLTYPE <- "Microglia_Homeostatic"

# ---- 其他参数 ----
CORES      <- 8
TOP_N_GENE <- 200

# 分组配色
GROUP_COLORS  <- c("Tumor"="#E41A1C", "Blank"="#FF7F00", "Control"="#377EB8")
SAMPLE_COLORS <- c(
  "gan56"="#E41A1C",  "liver3"="#FB9A99",  "liver8"="#FDBF6F",
  "bla_liver"="#FF7F00", "con_liver"="#377EB8",
  "nao56"="#A6CEE3",  "brain3"="#1F78B4",  "brain8"="#B2DF8A",
  "bla_brain"="#33A02C", "con_brain"="#6A3D9A"
)

OUTPUT_ROOT <- "pseudotime_results"
if (!dir.exists(OUTPUT_ROOT)) dir.create(OUTPUT_ROOT, recursive=TRUE)

# ==================== 2. 工具函数 ====================

load_one_sample <- function(s_info) {
  fp <- file.path(s_info$name, paste0(s_info$name, "_manual_annotated.rds"))
  if (!file.exists(fp)) {
    cat(sprintf("  ✗ %-15s 文件不存在: %s\n", s_info$name, fp))
    return(NULL)
  }
  obj <- tryCatch(readRDS(fp), error=function(e) {
    cat(sprintf("  ✗ %-15s 读取失败: %s\n", s_info$name, e$message)); NULL
  })
  if (is.null(obj)) return(NULL)
  obj <- JoinLayers(obj)
  obj <- RenameCells(obj, add.cell.id=s_info$name)
  obj$sample <- s_info$name
  obj$group  <- s_info$group
  obj$organ  <- s_info$organ
  cat(sprintf("  ✓ %-15s 细胞数: %5d  分组: %s\n",
              s_info$name, ncol(obj), s_info$group))
  return(obj)
}

merge_samples <- function(sample_list) {
  objs <- Filter(Negate(is.null), lapply(sample_list, load_one_sample))
  if (length(objs) == 0) stop("没有成功加载的样本")
  cat(sprintf("  合并 %d 个样本，共 %d 个细胞\n",
              length(objs), sum(sapply(objs, ncol))))
  JoinLayers(Reduce(merge, objs))
}

seurat_to_cds <- function(seurat_obj) {
  counts <- tryCatch(
    GetAssayData(seurat_obj, assay="RNA", layer="counts"),
    error=function(e) GetAssayData(seurat_obj, assay="RNA", layer="data")
  )
  if (!inherits(counts, "dgCMatrix")) counts <- as(counts, "dgCMatrix")

  gene_meta <- data.frame(
    gene_short_name = rownames(counts),
    row.names       = rownames(counts)
  )
  cds <- new_cell_data_set(
    counts,
    cell_metadata = seurat_obj@meta.data,
    gene_metadata = gene_meta
  )

  # 注入已有 UMAP 坐标
  umap_key <- intersect(c("umap","UMAP"), names(seurat_obj@reductions))
  if (length(umap_key) > 0) {
    embed <- Embeddings(seurat_obj, umap_key[1])[, 1:2]
    reducedDim(cds, "UMAP") <- embed
    cat("  ✓ 注入已有 UMAP 坐标\n")
  } else {
    cat("  ⚠ 未找到 UMAP，重新计算...\n")
    cds <- preprocess_cds(cds, num_dim=50)
    cds <- reduce_dimension(cds, reduction_method="UMAP")
  }
  return(cds)
}

set_root_node <- function(cds, root_celltype) {
  cell_ids <- colnames(cds)[colData(cds)[[CELLTYPE_COL]] == root_celltype]
  if (length(cell_ids) < 3) {
    cat(sprintf("  ⚠ '%s' 细胞数不足，使用默认起点\n", root_celltype))
    return(order_cells(cds))
  }
  center <- colMeans(reducedDim(cds, "UMAP")[cell_ids, , drop=FALSE])
  nodes  <- t(cds@principal_graph_aux[["UMAP"]]$dp_mst)
  dists  <- sqrt(rowSums(sweep(nodes, 2, center)^2))
  root   <- names(which.min(dists))
  cds    <- order_cells(cds, root_pr_nodes=root)
  cat(sprintf("  ✓ 根节点: '%s' → 轨迹节点 %s\n", root_celltype, root))
  return(cds)
}

save_fig <- function(p, path, w=10, h=8) {
  ggsave(paste0(path,".pdf"), p, width=w, height=h, dpi=300)
  ggsave(paste0(path,".png"), p, width=w, height=h, dpi=200)
}

# ==================== 3. 可视化套件 ====================

plot_suite <- function(cds, organ, out_dir) {
  ct_vals <- as.character(unique(colData(cds)[[CELLTYPE_COL]]))
  n_ct    <- length(ct_vals)
  ct_pal  <- setNames(
    colorRampPalette(brewer.pal(min(12,n_ct),"Paired"))(n_ct),
    ct_vals
  )
  grp_vals <- as.character(unique(colData(cds)[["group"]]))
  smp_vals <- as.character(unique(colData(cds)[["sample"]]))

  base_theme <- theme_bw(base_size=11) +
    theme(plot.title=element_text(hjust=0.5, face="bold", size=12),
          legend.key.size=unit(0.4,"cm"),
          legend.text=element_text(size=8))

  # 1. 伪时间轨迹
  p1 <- plot_cells(cds,
    color_cells_by="pseudotime", label_cell_groups=FALSE,
    label_leaves=TRUE, label_branch_points=TRUE,
    graph_label_size=3, cell_size=0.5,
    trajectory_graph_color="grey35",
    trajectory_graph_segment_size=0.7
  ) + scale_color_viridis_c(option="plasma", name="Pseudotime") +
    labs(title=paste0(organ, " · 伪时间轨迹")) + base_theme
  save_fig(p1, file.path(out_dir,"01_pseudotime_trajectory"))

  # 2. 细胞类型
  p2 <- plot_cells(cds,
    color_cells_by=CELLTYPE_COL, label_cell_groups=FALSE,
    label_leaves=FALSE, label_branch_points=FALSE,
    cell_size=0.5, trajectory_graph_color="grey50",
    trajectory_graph_segment_size=0.5
  ) + scale_color_manual(values=ct_pal, name="细胞类型") +
    guides(color=guide_legend(override.aes=list(size=3), ncol=1)) +
    labs(title=paste0(organ, " · 细胞类型")) + base_theme
  save_fig(p2, file.path(out_dir,"02_celltype"))

  # 3. 分组
  p3 <- plot_cells(cds,
    color_cells_by="group", label_cell_groups=FALSE,
    label_leaves=FALSE, label_branch_points=FALSE,
    cell_size=0.5, trajectory_graph_color="grey50",
    trajectory_graph_segment_size=0.5
  ) + scale_color_manual(values=GROUP_COLORS[grp_vals], name="分组") +
    guides(color=guide_legend(override.aes=list(size=3))) +
    labs(title=paste0(organ, " · 样本分组")) + base_theme
  save_fig(p3, file.path(out_dir,"03_group"))

  # 4. 样本来源
  p4 <- plot_cells(cds,
    color_cells_by="sample", label_cell_groups=FALSE,
    label_leaves=FALSE, label_branch_points=FALSE,
    cell_size=0.5, trajectory_graph_color="grey50",
    trajectory_graph_segment_size=0.5
  ) + scale_color_manual(values=SAMPLE_COLORS[smp_vals], name="样本") +
    guides(color=guide_legend(override.aes=list(size=3))) +
    labs(title=paste0(organ, " · 样本来源")) + base_theme
  save_fig(p4, file.path(out_dir,"04_sample"))

  # 5. 四图拼合
  save_fig((p1|p2)/(p3|p4), file.path(out_dir,"05_combined_4panel"), w=20, h=16)
  cat("  ✓ 01-05 轨迹图\n")

  # 伪时间值数据框
  pt_df <- data.frame(
    pseudotime = pseudotime(cds),
    group      = as.character(colData(cds)[["group"]]),
    sample     = as.character(colData(cds)[["sample"]]),
    celltype   = as.character(colData(cds)[[CELLTYPE_COL]])
  ) %>% dplyr::filter(!is.infinite(pseudotime), !is.na(pseudotime))

  # 6. 各组分布
  p6 <- ggplot(pt_df, aes(x=group, y=pseudotime, fill=group)) +
    geom_violin(alpha=0.75, trim=FALSE) +
    geom_boxplot(width=0.14, fill="white",
                 outlier.size=0.3, outlier.alpha=0.3) +
    scale_fill_manual(values=GROUP_COLORS) +
    labs(title=paste0(organ," · 各组伪时间分布"), x=NULL, y="Pseudotime") +
    base_theme + theme(legend.position="none")
  save_fig(p6, file.path(out_dir,"06_pseudotime_by_group"), w=7, h=6)

  # 7. 各细胞类型分布
  med_ord <- pt_df %>%
    dplyr::group_by(celltype) %>%
    dplyr::summarise(med=median(pseudotime),.groups="drop") %>%
    dplyr::arrange(med) %>% dplyr::pull(celltype)
  pt_df$celltype <- factor(pt_df$celltype, levels=med_ord)

  p7 <- ggplot(pt_df, aes(x=celltype, y=pseudotime, fill=celltype)) +
    geom_violin(alpha=0.75, trim=FALSE) +
    geom_boxplot(width=0.14, fill="white",
                 outlier.size=0.2, outlier.alpha=0.3) +
    scale_fill_manual(values=ct_pal) +
    coord_flip() +
    labs(title=paste0(organ," · 各细胞类型伪时间"),
         x=NULL, y="Pseudotime") +
    base_theme + theme(legend.position="none")
  save_fig(p7, file.path(out_dir,"07_pseudotime_by_celltype"),
           w=9, h=max(5, n_ct*0.55+2))

  # 8. 细胞类型 × 分组分面
  p8 <- ggplot(pt_df, aes(x=group, y=pseudotime, fill=group)) +
    geom_violin(alpha=0.7, trim=FALSE) +
    geom_boxplot(width=0.18, fill="white", outlier.size=0.2) +
    scale_fill_manual(values=GROUP_COLORS) +
    facet_wrap(~celltype, scales="free_y", ncol=4) +
    labs(title=paste0(organ," · 细胞类型 × 分组"),
         x=NULL, y="Pseudotime") +
    theme_bw(base_size=10) +
    theme(plot.title=element_text(hjust=0.5, face="bold"),
          axis.text.x=element_text(angle=30, hjust=1, size=8),
          strip.text=element_text(size=8),
          legend.position="none")
  nc <- min(4, n_ct); nr <- ceiling(n_ct/nc)
  save_fig(p8, file.path(out_dir,"08_celltype_x_group"),
           w=nc*3.5, h=nr*3+1)
  cat("  ✓ 06-08 分布图\n")

  return(pt_df)
}

# ==================== 4. 伪时间差异基因 ====================

run_graph_test <- function(cds, out_dir) {
  cat("  graph_test 检验中（耗时较长）...\n")
  cds_f <- cds[rowSums(exprs(cds) > 0) >= 10, ]
  cat(sprintf("  过滤后基因数: %d\n", nrow(cds_f)))

  pr_res <- tryCatch(
    graph_test(cds_f, neighbor_graph="principal_graph", cores=CORES),
    error=function(e) {
      cat(sprintf("  ✗ graph_test 失败: %s\n", e$message)); NULL
    }
  )
  if (is.null(pr_res)) return(NULL)

  pr_res <- pr_res %>% dplyr::arrange(q_value, desc(morans_I))
  sig    <- pr_res %>% dplyr::filter(q_value < 0.05)
  top_n  <- head(sig, TOP_N_GENE)

  write.csv(pr_res, file.path(out_dir,"gene_test_all.csv"),    row.names=FALSE)
  write.csv(top_n,  file.path(out_dir,"gene_test_top_sig.csv"),row.names=FALSE)
  cat(sprintf("  ✓ 显著基因 %d 个（q<0.05），Top%d 已保存\n",
              nrow(sig), nrow(top_n)))

  if (nrow(top_n) >= 5) {
    genes_show <- head(top_n$gene_short_name, 50)
    tryCatch({
      p_dyn <- plot_genes_in_pseudotime(
        cds_f[genes_show,], color_cells_by="pseudotime",
        min_expr=0.5, ncol=5
      )
      nh <- ceiling(length(genes_show)/5)*3+2
      ggsave(file.path(out_dir,"09_dynamic_genes.pdf"),
             p_dyn, width=22, height=nh, dpi=300)
      ggsave(file.path(out_dir,"09_dynamic_genes.png"),
             p_dyn, width=22, height=nh, dpi=150)
      cat("  ✓ 09 动态基因图（Top50）\n")
    }, error=function(e)
      cat(sprintf("  ⚠ 动态基因图失败: %s\n", e$message)))
  }
  return(top_n)
}

# ==================== 5. 单器官主流程 ====================

run_organ <- function(organ_name, sample_list, root_celltype) {
  cat(sprintf("\n%s\n器官: %s\n%s\n\n",
              strrep("=",60), organ_name, strrep("=",60)))
  out_dir <- file.path(OUTPUT_ROOT, organ_name)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE)

  # 5.1 加载合并
  cat("加载样本...\n")
  seurat_m <- tryCatch(merge_samples(sample_list), error=function(e) {
    cat(sprintf("✗ 合并失败: %s\n", e$message)); NULL
  })
  if (is.null(seurat_m)) return(NULL)

  has_umap <- length(intersect(c("umap","UMAP"),
                               names(seurat_m@reductions))) > 0
  cat(sprintf("  UMAP: %s\n", ifelse(has_umap,"已存在✓","将重新计算")))

  # 5.2 转换
  cat("\nSeurat → CellDataSet...\n")
  cds <- tryCatch(seurat_to_cds(seurat_m), error=function(e) {
    cat(sprintf("✗ 转换失败: %s\n", e$message)); NULL
  })
  if (is.null(cds)) return(NULL)
  cat(sprintf("  细胞数: %d  基因数: %d\n", ncol(cds), nrow(cds)))

  # 5.3 轨迹学习
  cat("\n学习轨迹...\n")
  cds <- tryCatch({
    cds <- cluster_cells(cds, reduction_method="UMAP",
                         resolution=1e-3, random_seed=42)
    cat("  ✓ cluster_cells\n")
    cds <- learn_graph(cds, use_partition=TRUE,
                       learn_graph_control=list(
                         minimal_branch_len=10,
                         geodesic_distance_ratio=0.3))
    cat("  ✓ learn_graph\n")
    cds
  }, error=function(e) {
    cat(sprintf("✗ 轨迹学习失败: %s\n", e$message)); NULL
  })
  if (is.null(cds)) return(NULL)

  # 5.4 设根节点
  cat("\n设定根节点...\n")
  cds <- tryCatch(
    set_root_node(cds, root_celltype),
    error=function(e) {
      cat(sprintf("⚠ 根节点失败(%s)，使用默认\n", e$message))
      order_cells(cds)
    }
  )

  # 5.5 保存伪时间值
  pt_csv <- data.frame(
    cell       = colnames(cds),
    pseudotime = pseudotime(cds),
    sample     = as.character(colData(cds)[["sample"]]),
    group      = as.character(colData(cds)[["group"]]),
    celltype   = as.character(colData(cds)[[CELLTYPE_COL]])
  )
  write.csv(pt_csv,
            file.path(out_dir,"pseudotime_values.csv"), row.names=FALSE)
  n_fin <- sum(!is.infinite(pt_csv$pseudotime))
  cat(sprintf("  ✓ 伪时间值保存（有效: %d / %d 细胞）\n",
              n_fin, nrow(pt_csv)))

  # 5.6 可视化
  cat("\n生成可视化图...\n")
  pt_df <- tryCatch(plot_suite(cds, organ_name, out_dir),
    error=function(e) {
      cat(sprintf("✗ 可视化失败: %s\n", e$message)); NULL
    })

  # 5.7 差异基因
  cat("\n伪时间差异基因...\n")
  sig_genes <- tryCatch(run_graph_test(cds, out_dir),
    error=function(e) {
      cat(sprintf("✗ 差异基因失败: %s\n", e$message)); NULL
    })

  # 5.8 保存 CDS
  rds_out <- file.path(out_dir, paste0(organ_name,"_pseudotime.rds"))
  saveRDS(cds, rds_out)
  cat(sprintf("  ✓ CDS 保存: %s\n", rds_out))

  cat(sprintf("\n【%s 完成】输出文件:\n", organ_name))
  for (f in sort(list.files(out_dir, pattern="\\.(pdf|png|csv|rds)$")))
    cat(sprintf("  %s\n", f))

  return(list(cds=cds, pt_df=pt_df, sig_genes=sig_genes))
}

# ==================== 6. 执行两个器官 ====================

t_start <- proc.time()

res_liver <- run_organ("Liver", LIVER_SAMPLES, LIVER_ROOT_CELLTYPE)
res_brain <- run_organ("Brain", BRAIN_SAMPLES, BRAIN_ROOT_CELLTYPE)

# ==================== 7. 跨器官对比图 ====================

if (!is.null(res_liver$pt_df) && !is.null(res_brain$pt_df)) {
  cat("\n========== 跨器官伪时间对比图 ==========\n")
  tryCatch({
    cmp_dir <- file.path(OUTPUT_ROOT, "comparison")
    if (!dir.exists(cmp_dir)) dir.create(cmp_dir)

    pt_all <- dplyr::bind_rows(
      res_liver$pt_df %>% dplyr::mutate(organ="Liver"),
      res_brain$pt_df %>% dplyr::mutate(organ="Brain")
    ) %>% dplyr::filter(!is.infinite(pseudotime), !is.na(pseudotime))

    # 分组对比
    p_cmp <- ggplot(pt_all, aes(x=group, y=pseudotime, fill=group)) +
      geom_violin(alpha=0.75, trim=FALSE) +
      geom_boxplot(width=0.14, fill="white", outlier.size=0.3) +
      scale_fill_manual(values=GROUP_COLORS) +
      facet_wrap(~organ, scales="free") +
      labs(title="两器官各组伪时间分布对比",
           x=NULL, y="Pseudotime") +
      theme_bw(base_size=12) +
      theme(plot.title=element_text(hjust=0.5, face="bold"),
            strip.text=element_text(size=12, face="bold"),
            axis.text.x=element_text(angle=20, hjust=1))
    save_fig(p_cmp,
             file.path(cmp_dir,"cross_organ_pseudotime"), w=12, h=6)

    # 样本级对比
    p_smp <- ggplot(pt_all, aes(x=sample, y=pseudotime, fill=group)) +
      geom_violin(alpha=0.7, trim=FALSE) +
      geom_boxplot(width=0.14, fill="white", outlier.size=0.2) +
      scale_fill_manual(values=GROUP_COLORS) +
      facet_wrap(~organ, scales="free_x") +
      labs(title="两器官各样本伪时间分布",
           x=NULL, y="Pseudotime") +
      theme_bw(base_size=11) +
      theme(plot.title=element_text(hjust=0.5, face="bold"),
            axis.text.x=element_text(angle=30, hjust=1, size=9),
            strip.text=element_text(size=11, face="bold"))
    save_fig(p_smp,
             file.path(cmp_dir,"cross_organ_by_sample"), w=14, h=6)
    cat("  ✓ 跨器官对比图已保存\n")
  }, error=function(e)
    cat(sprintf("✗ 对比图失败: %s\n", e$message)))
}

# ==================== 8. 总结 ====================

cat(sprintf("\n%s\n全部完成！耗时: %.1f 秒\n%s\n",
            strrep("=",40),
            (proc.time()-t_start)["elapsed"],
            strrep("=",40)))
cat(sprintf("\n输出目录: %s/\n", OUTPUT_ROOT))
cat("├── Liver/\n")
cat("│   ├── 01_pseudotime_trajectory.pdf/png\n")
cat("│   ├── 02_celltype.pdf/png\n")
cat("│   ├── 03_group.pdf/png\n")
cat("│   ├── 04_sample.pdf/png\n")
cat("│   ├── 05_combined_4panel.pdf/png\n")
cat("│   ├── 06_pseudotime_by_group.pdf/png\n")
cat("│   ├── 07_pseudotime_by_celltype.pdf/png\n")
cat("│   ├── 08_celltype_x_group.pdf/png\n")
cat("│   ├── 09_dynamic_genes.pdf/png\n")
cat("│   ├── pseudotime_values.csv\n")
cat("│   ├── gene_test_all.csv\n")
cat("│   ├── gene_test_top_sig.csv\n")
cat("│   └── Liver_pseudotime.rds\n")
cat("├── Brain/（同上结构）\n")
cat("└── comparison/\n")
cat("    ├── cross_organ_pseudotime.pdf/png\n")
cat("    └── cross_organ_by_sample.pdf/png\n")
