# ==================== 肝脏肿瘤样本人鼠细胞重建表 ====================
# 文件名: rebuild_species_table.R
# 用法:   Rscript rebuild_species_table.R
#
# 【目标】
#   从原始 10X features.tsv.gz 重新识别人鼠基因，
#   计算每个 barcode 的人鼠表达比例，输出细胞物种分配表，
#   同时注入现有 RDS 对象（如存在）供后续 TF 分析使用
#
# 【仅处理三个肿瘤样本（人鼠混合）】
#   gan56 / liver3 / liver8
# =============================================================================

cat("========================================\n")
cat("肝脏肿瘤样本人鼠细胞重建\n")
cat("========================================\n\n")

.libPaths("/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/lib/R/library")
suppressPackageStartupMessages({
  library(Seurat); library(Matrix)
  library(ggplot2); library(dplyr); library(patchwork)
})

setwd("/mnt/BioDisk/wangmingyuan/SGS10_analysis")

# ==================== 1. 配置 ====================

SAMPLES <- list(
  gan56  = list(
    data_dir = "data/gan56/gan_56_filtered_feature_bc_matrix",
    rds_path = "gan56/gan56_manual_annotated.rds"
  ),
  liver3 = list(
    data_dir = "data/liver_samples/25090559_Liver_3_expression/Liver_3_filtered_feature_bc_matrix",
    rds_path = "liver3/liver3_manual_annotated.rds"
  ),
  liver8 = list(
    data_dir = "data/liver_samples/25090559_Liver_8_expression/Liver_8_filtered_feature_bc_matrix",
    rds_path = "liver8/liver8_manual_annotated.rds"
  )
)

# 原脚本的同源映射（仅用于基因名还原，不影响物种判断）
HOMOLOG_MAP <- c(
  "ALB"="Alb","AFP"="Afp","TTR"="Ttr","HNF4A"="Hnf4a","HNF1A"="Hnf1a",
  "CYP2E1"="Cyp2e1","CYP1A2"="Cyp1a2","CYP3A4"="Cyp3a11",
  "PECAM1"="Pecam1","CDH5"="Cdh5","VWF"="Vwf","KDR"="Kdr",
  "CLEC4G"="Clec4g","LYVE1"="Lyve1",
  "PDGFRA"="Pdgfra","PDGFRB"="Pdgfrb","DES"="Des","ACTA2"="Acta2",
  "COL1A1"="Col1a1","LRAT"="Lrat",
  "CD68"="Cd68","CD163"="Cd163","VSIG4"="Vsig4","C1QA"="C1qa",
  "MRC1"="Mrc1","TREM2"="Trem2","GPNMB"="Gpnmb",
  "CD3D"="Cd3d","CD3E"="Cd3e","CD4"="Cd4","CD8A"="Cd8a",
  "FOXP3"="Foxp3","CD19"="Cd19","MS4A1"="Ms4a1","CD79A"="Cd79a",
  "NKG7"="Nkg7","KLRB1"="Klrb1",
  "CD14"="Cd14","LYZ"="Lyz","S100A8"="S100a8","S100A9"="S100a9",
  "KRT7"="Krt7","KRT8"="Krt8","KRT19"="Krt19","EPCAM"="Epcam",
  "GAPDH"="Gapdh","ACTB"="Actb"
)

OUTPUT_DIR <- "species_assignment"
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR)

# ==================== 2. 核心函数 ====================

# 清理基因名（与原脚本完全一致）
clean_gene_name <- function(x) {
  x <- gsub("^GRCh38[-_]","",x)
  x <- gsub("^mm10[-_]","",x)
  x <- gsub("\\.[0-9]+$","",x)
  return(x)
}

# 从原始 features.tsv.gz 读取基因物种信息
read_features_with_species <- function(data_dir) {
  ft_file <- file.path(data_dir, "features.tsv.gz")
  if (!file.exists(ft_file)) stop("features.tsv.gz 不存在: ", ft_file)

  ft <- read.table(ft_file, header=FALSE, sep="\t", stringsAsFactors=FALSE)
  colnames(ft) <- c("gene_id","gene_name","feature_type")

  # ---- 物种判断（在清理前，用原始前缀）----
  ft$species <- dplyr::case_when(
    grepl("^GRCh38[-_]", ft$gene_name)                               ~ "human",
    grepl("^mm10[-_]",   ft$gene_name)                               ~ "mouse",
    # 没有前缀：全大写（含数字/连字符）→ 人；首大写后跟小写 → 鼠
    grepl("^[A-Z][A-Z0-9\\-]+$", ft$gene_name) &
      !grepl("[a-z]", ft$gene_name)                                  ~ "human",
    grepl("^[A-Z][a-z]", ft$gene_name)                               ~ "mouse",
    TRUE                                                              ~ "ambiguous"
  )

  # ---- 清理基因名（去前缀）----
  ft$clean <- clean_gene_name(ft$gene_name)

  # ---- 应用同源映射（还原成小鼠格式，与原脚本一致）----
  ft$final <- ft$clean
  for (hg in names(HOMOLOG_MAP)) {
    idx <- which(ft$clean == hg)
    if (length(idx)>0) ft$final[idx] <- HOMOLOG_MAP[hg]
  }

  # ---- 处理重复（去重后 final_name → species）----
  # 如果同一 final_name 对应人鼠两个基因（同源映射导致），保留最先出现的物种
  ft_dedup <- ft[!duplicated(ft$final), ]
  species_map <- setNames(ft_dedup$species, ft_dedup$final)

  stats <- table(ft$species)
  cat(sprintf("    features: 人类=%d 小鼠=%d 模糊=%d\n",
              sum(ft$species=="human"), sum(ft$species=="mouse"),
              sum(ft$species=="ambiguous")))

  return(list(ft=ft, species_map=species_map))
}

# 读取 10X count 矩阵并重命名
read_counts_renamed <- function(data_dir, ft) {
  counts <- Read10X(data_dir)
  if (nrow(counts) == nrow(ft)) {
    # 处理重复 final_name：保留表达量最高的行
    rnames <- ft$final
    dup    <- unique(rnames[duplicated(rnames)])
    if (length(dup)>0) {
      row_sum <- Matrix::rowSums(counts)
      keep    <- rep(TRUE, nrow(counts))
      for (g in dup) {
        idx  <- which(rnames==g)
        best <- idx[which.max(row_sum[idx])]
        keep[setdiff(idx,best)] <- FALSE
      }
      counts <- counts[keep,]
      ft     <- ft[keep,]
    }
    rownames(counts) <- ft$final[1:nrow(counts)]
  } else {
    cat(sprintf("    ⚠ 行数不匹配（矩阵%d vs features%d），使用原始基因名\n",
                nrow(counts), nrow(ft)))
  }
  return(counts)
}

# 计算每个 barcode 的人鼠比例
calc_species_ratio <- function(counts, species_map) {
  all_genes   <- rownames(counts)
  human_genes <- intersect(all_genes, names(species_map)[species_map=="human"])
  mouse_genes <- intersect(all_genes, names(species_map)[species_map=="mouse"])

  cat(sprintf("    人类基因数: %d  小鼠基因数: %d\n",
              length(human_genes), length(mouse_genes)))

  h_sum <- Matrix::colSums(counts[human_genes,,drop=FALSE])
  m_sum <- Matrix::colSums(counts[mouse_genes,,drop=FALSE])
  total <- h_sum + m_sum

  human_ratio <- ifelse(total>0, h_sum/total, 0)
  return(list(human_ratio=human_ratio, human_sum=h_sum,
              mouse_sum=m_sum, total=total))
}

# 自动检测阈值（双峰谷值）
auto_threshold <- function(human_ratio, default=0.30) {
  h <- hist(human_ratio, breaks=100, plot=FALSE)
  # 找局部最小值（谷值），限定在0.05~0.95之间
  cnts    <- h$counts
  valleys <- which(diff(sign(diff(cnts)))>0) + 1
  valid   <- valleys[h$mids[valleys]>0.05 & h$mids[valleys]<0.95]

  if (length(valid)>0) {
    # 取最深的谷（绝对值最小的 count）
    best_v <- valid[which.min(cnts[valid])]
    thr    <- h$mids[best_v]
    cat(sprintf("    双峰谷值检测: %.3f\n", thr))
  } else {
    thr <- default
    cat(sprintf("    未检测到双峰，使用默认阈值: %.3f\n", thr))
  }
  return(thr)
}

# ==================== 3. 逐样本处理 ====================

all_species_dfs <- list()

for (sname in names(SAMPLES)) {
  cfg <- SAMPLES[[sname]]
  cat(sprintf("\n%s\n样本: %s\n%s\n",
              strrep("=",50), sname, strrep("=",50)))

  # 3.1 检查原始数据
  if (!dir.exists(cfg$data_dir)) {
    cat(sprintf("  ✗ 数据目录不存在: %s\n", cfg$data_dir)); next
  }

  # 3.2 读取 features 物种信息
  cat("  读取 features.tsv.gz...\n")
  ft_result <- tryCatch(
    read_features_with_species(cfg$data_dir),
    error=function(e){ cat(sprintf("  ✗ %s\n",e$message)); NULL })
  if (is.null(ft_result)) next
  ft         <- ft_result$ft
  species_map<- ft_result$species_map

  # 3.3 读取 count 矩阵
  cat("  读取 count 矩阵...\n")
  counts <- tryCatch(
    read_counts_renamed(cfg$data_dir, ft),
    error=function(e){ cat(sprintf("  ✗ %s\n",e$message)); NULL })
  if (is.null(counts)) next
  cat(sprintf("  矩阵: %d genes × %d barcodes\n",
              nrow(counts), ncol(counts)))

  # 3.4 计算人鼠比例
  cat("  计算人鼠比例...\n")
  ratio_res  <- calc_species_ratio(counts, species_map)
  human_ratio <- ratio_res$human_ratio

  cat(sprintf("  比例统计: min=%.3f  median=%.3f  max=%.3f\n",
              min(human_ratio), median(human_ratio), max(human_ratio)))

  # 3.5 确定阈值
  thr <- auto_threshold(human_ratio)

  # 3.6 分配物种
  cell_species <- ifelse(human_ratio > thr, "human","mouse")
  n_human <- sum(cell_species=="human")
  n_mouse <- sum(cell_species=="mouse")
  cat(sprintf("  分配: 人类=%d (%.1f%%)  小鼠=%d (%.1f%%)\n",
              n_human, n_human/length(cell_species)*100,
              n_mouse, n_mouse/length(cell_species)*100))

  # 3.7 可视化分布
  pdf(file.path(OUTPUT_DIR, paste0(sname,"_ratio_dist.pdf")),
      width=8, height=5)
  hist(human_ratio, breaks=100, col="#6BAED6", border="white",
       main=paste0(sname," · 人类基因比例分布（",
                   n_human,"人 / ",n_mouse,"鼠）"),
       xlab="Human gene ratio", ylab="Cell count")
  abline(v=thr, col="#E41A1C", lwd=2, lty=2)
  legend("topright", legend=sprintf("cutoff=%.3f",thr),
         col="#E41A1C", lty=2, lwd=2)
  dev.off()

  # 3.8 构建输出表（以原始 barcode 为 key）
  sp_df <- data.frame(
    barcode      = names(human_ratio),
    sample       = sname,
    group        = "Tumor",
    human_count  = ratio_res$human_sum,
    mouse_count  = ratio_res$mouse_sum,
    total_count  = ratio_res$total,
    human_ratio  = round(human_ratio, 4),
    cell_species = cell_species,
    threshold    = thr,
    stringsAsFactors=FALSE
  )
  all_species_dfs[[sname]] <- sp_df

  # 3.9 注入现有 RDS（如果存在）
  if (file.exists(cfg$rds_path)) {
    cat("  注入 RDS 对象...\n")
    seurat_obj <- tryCatch(readRDS(cfg$rds_path),
                           error=function(e){ cat("  ✗ RDS读取失败\n"); NULL })
    if (!is.null(seurat_obj)) {
      seurat_barcodes <- colnames(seurat_obj)
      # 尝试匹配（去样本前缀）
      bc_clean <- sub(paste0("^",sname,"_"),"", seurat_barcodes)
      m <- match(bc_clean, sp_df$barcode)
      valid <- !is.na(m)
      cat(sprintf("  匹配: %d / %d\n", sum(valid), length(seurat_barcodes)))

      seurat_obj$human_ratio  <- NA_real_
      seurat_obj$cell_species <- NA_character_
      seurat_obj$human_ratio[valid]  <- sp_df$human_ratio[m[valid]]
      seurat_obj$cell_species[valid] <- sp_df$cell_species[m[valid]]

      # 未匹配到的用基因名规则补充
      unmatched <- which(is.na(seurat_obj$human_ratio))
      if (length(unmatched)>0) {
        cat(sprintf("  %d 个细胞未匹配，用基因名规则补充...\n",
                    length(unmatched)))
        expr_sub <- GetAssayData(seurat_obj[,unmatched], layer="counts")
        g_all <- rownames(expr_sub)
        g_h   <- g_all[grepl("^[A-Z][A-Z0-9\\-]+$",g_all)&!grepl("[a-z]",g_all)]
        g_m   <- g_all[grepl("^[A-Z][a-z]",g_all)]
        hs <- Matrix::colSums(expr_sub[intersect(g_h,rownames(expr_sub)),,drop=FALSE])
        ms <- Matrix::colSums(expr_sub[intersect(g_m,rownames(expr_sub)),,drop=FALSE])
        tot<- hs+ms; hr <- ifelse(tot>0,hs/tot,0)
        seurat_obj$human_ratio[unmatched]  <- round(hr,4)
        seurat_obj$cell_species[unmatched] <- ifelse(hr>thr,"human","mouse")
      }

      # 在 RDS 里也存一份含样本前缀的 barcode（供 TF 脚本直接用）
      seurat_obj$cell_species[is.na(seurat_obj$cell_species)] <- "mouse"

      # 保存更新后的 RDS
      saveRDS(seurat_obj, cfg$rds_path)
      cat(sprintf("  ✓ RDS 已更新: %s\n", cfg$rds_path))

      # 生成 UMAP 验证图（如果有 UMAP）
      if ("umap" %in% names(seurat_obj@reductions) ||
          "UMAP" %in% names(seurat_obj@reductions)) {
        red <- ifelse("umap" %in% names(seurat_obj@reductions),"umap","UMAP")
        pdf(file.path(OUTPUT_DIR, paste0(sname,"_species_umap.pdf")),
            width=15, height=5)
        p1 <- DimPlot(seurat_obj, reduction=red, group.by="cell_species",
                      cols=c("human"="#E41A1C","mouse"="#377EB8"),
                      pt.size=0.3) +
          labs(title=paste0(sname," · 人鼠细胞")) +
          theme(plot.title=element_text(hjust=0.5,face="bold"))
        p2 <- FeaturePlot(seurat_obj, features="human_ratio",
                          reduction=red, order=TRUE, pt.size=0.3) +
          scale_color_gradient(low="grey90",high="#E41A1C",name="ratio") +
          labs(title=paste0(sname," · 人类基因比例")) +
          theme(plot.title=element_text(hjust=0.5,face="bold"))
        p3 <- DimPlot(seurat_obj, reduction=red, group.by="final_celltype",
                      label=TRUE, label.size=2.8, repel=TRUE, pt.size=0.3) +
          labs(title=paste0(sname," · 细胞类型")) +
          theme(plot.title=element_text(hjust=0.5,face="bold"))
        print(p1|p2|p3)
        dev.off()
        cat(sprintf("  ✓ UMAP验证图: %s\n",
                    file.path(OUTPUT_DIR,paste0(sname,"_species_umap.pdf"))))
      }
    }
  } else {
    cat(sprintf("  ⚠ RDS 不存在，跳过注入: %s\n", cfg$rds_path))
  }

  # 单样本 CSV
  csv_out <- file.path(OUTPUT_DIR, paste0(sname,"_species_assignment.csv"))
  write.csv(sp_df, csv_out, row.names=FALSE)
  cat(sprintf("  ✓ CSV: %s\n", csv_out))
}

# ==================== 4. 合并三个样本输出汇总表 ====================

if (length(all_species_dfs)>0) {
  cat("\n合并三样本汇总表...\n")
  combined_df <- dplyr::bind_rows(all_species_dfs)

  # 汇总表
  summary_df <- combined_df %>%
    dplyr::group_by(sample, cell_species) %>%
    dplyr::summarise(n=n(), .groups="drop") %>%
    tidyr::pivot_wider(names_from=cell_species, values_from=n, values_fill=0) %>%
    dplyr::mutate(
      total       = human + mouse,
      human_pct   = round(human/total*100,1),
      mouse_pct   = round(mouse/total*100,1)
    )

  cat("\n=== 三样本人鼠细胞汇总 ===\n")
  print(as.data.frame(summary_df))

  # 保存
  write.csv(combined_df,
            file.path(OUTPUT_DIR,"all_tumor_samples_species.csv"),
            row.names=FALSE)
  write.csv(summary_df,
            file.path(OUTPUT_DIR,"species_summary.csv"),
            row.names=FALSE)
  cat(sprintf("\n✓ 汇总 CSV: %s\n",
              file.path(OUTPUT_DIR,"all_tumor_samples_species.csv")))

  # 汇总柱状图
  pdf(file.path(OUTPUT_DIR,"species_summary_barplot.pdf"), width=8, height=5)
  plot_df <- combined_df %>%
    dplyr::count(sample, cell_species) %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(pct=n/sum(n)*100)

  p <- ggplot(plot_df,
              aes(x=sample, y=pct, fill=cell_species)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values=c("human"="#E41A1C","mouse"="#377EB8"),
                      name="物种") +
    geom_text(aes(label=sprintf("%.1f%%\n(n=%d)",pct,n)),
              position=position_stack(vjust=0.5),
              size=3.5, color="white", fontface="bold") +
    labs(title="肝脏肿瘤样本人鼠细胞比例",
         x=NULL, y="百分比 (%)") +
    theme_bw(base_size=12) +
    theme(plot.title=element_text(hjust=0.5,face="bold"))
  print(p)
  dev.off()
  cat(sprintf("✓ 汇总柱状图: %s\n",
              file.path(OUTPUT_DIR,"species_summary_barplot.pdf")))
}

# ==================== 5. 输出使用说明 ====================

cat("\n========================================\n")
cat("完成！输出文件:\n")
cat("========================================\n")
cat(sprintf("目录: %s/\n", OUTPUT_DIR))
cat("├── gan56_ratio_dist.pdf          # 人鼠比例分布直方图\n")
cat("├── liver3_ratio_dist.pdf\n")
cat("├── liver8_ratio_dist.pdf\n")
cat("├── gan56_species_umap.pdf        # UMAP验证图（人鼠+细胞类型）\n")
cat("├── liver3_species_umap.pdf\n")
cat("├── liver8_species_umap.pdf\n")
cat("├── gan56_species_assignment.csv  # 单样本细胞物种表\n")
cat("├── liver3_species_assignment.csv\n")
cat("├── liver8_species_assignment.csv\n")
cat("├── all_tumor_samples_species.csv # 三样本合并表\n")
cat("├── species_summary.csv           # 汇总统计\n")
cat("└── species_summary_barplot.pdf   # 汇总柱状图\n")
cat("\n列说明（*_species_assignment.csv）:\n")
cat("  barcode      : 原始细胞barcode\n")
cat("  sample       : 样本名\n")
cat("  human_count  : 人类基因总 UMI 数\n")
cat("  mouse_count  : 小鼠基因总 UMI 数\n")
cat("  human_ratio  : 人类基因比例（0~1）\n")
cat("  cell_species : human / mouse\n")
cat("  threshold    : 使用的切分阈值\n")
cat("\nRDS 对象新增列:\n")
cat("  $human_ratio    可用于 FeaturePlot 可视化\n")
cat("  $cell_species   可用于 DimPlot 或 subset 过滤\n")
cat("\n后续 TF 分析直接读取 RDS，cell_species 已就绪。\n")
cat(sprintf("完成时间: %s\n",format(Sys.time(),"%Y-%m-%d %H:%M:%S")))
