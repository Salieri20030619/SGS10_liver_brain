# ==================== PDX 通用工具函数 ====================
# 文件名: ~/SGS10_analysis/scripts/pdx_utils.R
# 功能:
#   1. 读取 10X 数据并保留物种前缀信息
#   2. 按 GRCh38- / mm10- 前缀拆分人鼠矩阵
#   3. 按 UMI 占比为每个 barcode 分类 Human / Mouse / Doublet
#   4. 基因名清理(去前缀、去版本号、去重)
#   5. 人鼠同源基因映射(仅在需要可视化对比时使用,默认不混合)

# ----- 物种前缀检测 -----
# CellRanger 混合参考基因名通常形如:
#   GRCh38-GENE / GRCh38_GENE / hg38-GENE
#   mm10-Gene  / mm10_Gene  / GRCm38-Gene
HUMAN_PREFIX_RE <- "^(GRCh38|hg38|HUMAN)[-_]"
MOUSE_PREFIX_RE <- "^(mm10|GRCm38|GRCm39|MOUSE)[-_]"

detect_species_in_features <- function(gene_names) {
  n_human <- sum(grepl(HUMAN_PREFIX_RE, gene_names, ignore.case = TRUE))
  n_mouse <- sum(grepl(MOUSE_PREFIX_RE, gene_names, ignore.case = TRUE))
  n_total <- length(gene_names)
  cat("  基因名物种检测:\n")
  cat("    总基因数:", n_total, "\n")
  cat("    人源前缀:", n_human, sprintf("(%.1f%%)\n", 100 * n_human / n_total))
  cat("    鼠源前缀:", n_mouse, sprintf("(%.1f%%)\n", 100 * n_mouse / n_total))
  list(
    n_human = n_human,
    n_mouse = n_mouse,
    is_mixed = (n_human > 100 & n_mouse > 100)
  )
}

# ----- 清理基因名 -----
clean_gene_name_human <- function(x) {
  x <- gsub(HUMAN_PREFIX_RE, "", x, ignore.case = TRUE)
  x <- gsub("\\.[0-9]+$", "", x)
  x
}
clean_gene_name_mouse <- function(x) {
  x <- gsub(MOUSE_PREFIX_RE, "", x, ignore.case = TRUE)
  x <- gsub("\\.[0-9]+$", "", x)
  x
}

# ----- 矩阵去重(同名基因取表达和最大那一行) -----
deduplicate_genes <- function(mat) {
  dups <- unique(rownames(mat)[duplicated(rownames(mat))])
  if (length(dups) == 0) return(mat)
  cat("    发现", length(dups), "个重复基因名,合并取最大表达行\n")
  row_max <- sparseMatrixStats::rowMaxs(mat)
  if (is.null(row_max)) row_max <- apply(mat, 1, max)
  keep <- rep(TRUE, nrow(mat))
  for (g in dups) {
    idx <- which(rownames(mat) == g)
    best <- idx[which.max(row_max[idx])]
    keep[setdiff(idx, best)] <- FALSE
  }
  mat[keep, , drop = FALSE]
}

# ----- 加载并拆分 PDX 数据 -----
# 返回 list(human = matrix or NULL, mouse = matrix or NULL, mixed = logical)
load_and_split_pdx <- function(data_dir, sample_name) {
  cat("\n[load_and_split_pdx] 样本:", sample_name, "\n")
  cat("  数据目录:", data_dir, "\n")
  if (!dir.exists(data_dir)) stop("数据目录不存在: ", data_dir)

  # 用 Read10X 直接读
  raw <- Seurat::Read10X(data_dir)
  if (is.list(raw)) {
    # 如果是多模态返回的 list, 取 Gene Expression
    if ("Gene Expression" %in% names(raw)) raw <- raw[["Gene Expression"]] else raw <- raw[[1]]
  }
  cat("  原始矩阵维度:", nrow(raw), "基因 ×", ncol(raw), "细胞\n")

  gene_names <- rownames(raw)
  info <- detect_species_in_features(gene_names)

  human_idx <- grepl(HUMAN_PREFIX_RE, gene_names, ignore.case = TRUE)
  mouse_idx <- grepl(MOUSE_PREFIX_RE, gene_names, ignore.case = TRUE)

  result <- list(human = NULL, mouse = NULL, is_mixed_ref = info$is_mixed)

  if (info$is_mixed) {
    cat("  ✓ 检测到人鼠混合参考,执行拆分\n")
    # 人源矩阵
    h_mat <- raw[human_idx, , drop = FALSE]
    rownames(h_mat) <- clean_gene_name_human(rownames(h_mat))
    h_mat <- deduplicate_genes(h_mat)
    cat("    人源矩阵:", nrow(h_mat), "基因 ×", ncol(h_mat), "细胞\n")

    # 鼠源矩阵
    m_mat <- raw[mouse_idx, , drop = FALSE]
    rownames(m_mat) <- clean_gene_name_mouse(rownames(m_mat))
    m_mat <- deduplicate_genes(m_mat)
    cat("    鼠源矩阵:", nrow(m_mat), "基因 ×", ncol(m_mat), "细胞\n")

    result$human <- h_mat
    result$mouse <- m_mat
  } else {
    cat("  ✓ 单物种参考,不拆分\n")
    # 不是混合参考,猜物种: 鼠基因名首字母大写,人全大写
    is_likely_mouse <- mean(grepl("^[A-Z][a-z]", gene_names)) > 0.5
    rn <- gene_names
    rn <- gsub("\\.[0-9]+$", "", rn)
    rownames(raw) <- rn
    raw <- deduplicate_genes(raw)
    if (is_likely_mouse) {
      result$mouse <- raw
      cat("    判定为鼠源单物种数据\n")
    } else {
      result$human <- raw
      cat("    判定为人源单物种数据\n")
    }
  }
  result
}

# ----- 按 UMI 占比为 barcode 分类 -----
# 返回 data.frame(barcode, human_umi, mouse_umi, total, human_frac, species)
classify_barcodes_by_species <- function(human_mat, mouse_mat,
                                         pure_threshold = 0.90,
                                         min_total_umi = 200) {
  if (is.null(human_mat) || is.null(mouse_mat)) {
    stop("classify_barcodes_by_species 需要同时提供人鼠矩阵")
  }
  common_bc <- intersect(colnames(human_mat), colnames(mouse_mat))
  cat("\n[classify_barcodes_by_species]\n")
  cat("  共有 barcode 数:", length(common_bc), "\n")
  cat("  纯物种阈值:", pure_threshold, "\n")
  cat("  最小总 UMI:", min_total_umi, "\n")

  human_mat <- human_mat[, common_bc, drop = FALSE]
  mouse_mat <- mouse_mat[, common_bc, drop = FALSE]

  human_umi <- Matrix::colSums(human_mat)
  mouse_umi <- Matrix::colSums(mouse_mat)
  total_umi <- human_umi + mouse_umi
  human_frac <- ifelse(total_umi > 0, human_umi / total_umi, 0)

  species <- rep("Doublet", length(common_bc))
  species[human_frac >= pure_threshold] <- "Human"
  species[(1 - human_frac) >= pure_threshold] <- "Mouse"
  species[total_umi < min_total_umi] <- "LowUMI"

  df <- data.frame(
    barcode    = common_bc,
    human_umi  = human_umi,
    mouse_umi  = mouse_umi,
    total      = total_umi,
    human_frac = human_frac,
    species    = species,
    stringsAsFactors = FALSE
  )

  cat("  分类结果:\n")
  print(table(df$species))
  df
}

# ----- 绘制 barnyard 图 (人 UMI vs 鼠 UMI 散点) -----
plot_barnyard <- function(class_df, sample_name, out_pdf) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  library(ggplot2)
  p <- ggplot(class_df, aes(x = mouse_umi + 1, y = human_umi + 1, color = species)) +
    geom_point(size = 0.4, alpha = 0.6) +
    scale_x_log10() + scale_y_log10() +
    scale_color_manual(values = c(
      Human = "#E64B35", Mouse = "#4DBBD5",
      Doublet = "#7E6148", LowUMI = "grey70"
    )) +
    labs(title = paste0(sample_name, " - Barnyard plot"),
         x = "Mouse UMI (log10)", y = "Human UMI (log10)") +
    theme_bw() +
    theme(legend.position = "right")
  pdf(out_pdf, width = 7, height = 6)
  print(p)
  dev.off()
  cat("  ✓ Barnyard 图已保存:", out_pdf, "\n")
}

cat("[pdx_utils.R] 工具函数已加载\n")
