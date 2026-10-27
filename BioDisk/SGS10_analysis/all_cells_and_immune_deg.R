#!/usr/bin/env Rscript
# =============================================================================
# 文件名: all_cells_and_immune_deg.R
# 用途:   分析所有细胞和免疫细胞中 Il6st 的表达差异
# 用法:   conda activate R_4.5.2 && Rscript all_cells_and_immune_deg.R
# =============================================================================

cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║     所有细胞 & 免疫细胞 DEG 分析（含 Il6st 专项诊断）                ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n\n")

# =============================================================================
# 1. 环境配置
# =============================================================================
.libPaths("/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/lib/R/library")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

# 解决函数冲突
select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
rename <- dplyr::rename

cat("✓ 包加载完成\n\n")

# =============================================================================
# 2. 参数设置
# =============================================================================
setwd("/mnt/BioDisk/wangmingyuan/SGS10_analysis")

FC_THRESHOLD <- 0.20      # log2FC 阈值
PADJ_CUT     <- 0.10      # 校正后p值阈值
MIN_PCT      <- 0.10      # 最小检测比例

output_dir <- "all_cells_immune_deg_results"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive=TRUE)

cat(sprintf("输出目录: %s\n", output_dir))
cat(sprintf("DEG参数: |log2FC| > %.2f, p_adj < %.2f, min.pct = %.2f\n\n",
            FC_THRESHOLD, PADJ_CUT, MIN_PCT))

# =============================================================================
# 3. 加载样本数据
# =============================================================================
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("步骤1: 加载样本数据\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

load_sample <- function(sample_name) {
  fp <- file.path(sample_name, paste0(sample_name, "_manual_annotated.rds"))
  
  if (!file.exists(fp)) {
    cat(sprintf("  ✗ 文件不存在: %s\n", fp))
    return(NULL)
  }
  
  obj <- tryCatch(readRDS(fp), error = function(e) {
    cat(sprintf("  ✗ 读取失败: %s - %s\n", sample_name, e$message))
    return(NULL)
  })
  
  if (is.null(obj)) return(NULL)
  
  obj <- JoinLayers(obj)
  obj <- RenameCells(obj, add.cell.id = sample_name)
  cat(sprintf("  ✓ %-12s 总细胞数: %d\n", sample_name, ncol(obj)))
  
  return(obj)
}

# 样本列表
tumor_samples <- c("nao56", "brain3", "brain8")
control_sample <- "con_brain"
blank_sample   <- "bla_brain"

# 加载所有样本
all_objs <- list()
for (s in c(tumor_samples, control_sample, blank_sample)) {
  all_objs[[s]] <- load_sample(s)
}
all_objs <- Filter(Negate(is.null), all_objs)
cat(sprintf("\n✓ 成功加载 %d 个样本\n\n", length(all_objs)))

# =============================================================================
# 4. 定义免疫细胞类型
# =============================================================================
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("步骤2: 定义免疫细胞类型\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

# 根据您的注释，免疫细胞类型包括：
immune_cell_types <- c(
  "Microglia",      # 小胶质细胞
  "Macrophage",     # 巨噬细胞
  "Monocyte",       # 单核细胞
  "DC",             # 树突状细胞
  "T cell",         # T细胞
  "B cell",         # B细胞
  "NK cell",        # NK细胞
  "Neutrophil",     # 中性粒细胞
  "Mast cell"       # 肥大细胞
)

# 也可以使用更广泛的定义：所有非上皮/非内皮/非基质细胞
# 或者检查实际存在的细胞类型
cat("免疫细胞类型列表:\n")
for (ct in immune_cell_types) {
  cat(sprintf("  - %s\n", ct))
}
cat("\n")

# =============================================================================
# 5. 提取各组细胞
# =============================================================================
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("步骤3: 提取各组细胞\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

extract_cells <- function(obj, group_label, cell_type_filter = NULL) {
  # 如果指定了细胞类型，则过滤
  if (!is.null(cell_type_filter)) {
    if ("final_celltype" %in% colnames(obj@meta.data)) {
      cells <- WhichCells(obj, expression = final_celltype %in% cell_type_filter)
    } else if ("celltype" %in% colnames(obj@meta.data)) {
      cells <- WhichCells(obj, expression = celltype %in% cell_type_filter)
    } else {
      cat(sprintf("  ⚠️ %s: 未找到细胞类型注释，使用全部细胞\n", group_label))
      cells <- colnames(obj)
    }
  } else {
    # 所有细胞
    cells <- colnames(obj)
  }
  
  if (length(cells) == 0) {
    cat(sprintf("  ✗ %s: 未找到符合条件的细胞\n", group_label))
    return(NULL)
  }
  
  sub_obj <- subset(obj, cells = cells)
  sub_obj$group <- group_label
  
  cat(sprintf("  ✓ %-12s %s: %d 个细胞\n", 
              group_label, 
              ifelse(is.null(cell_type_filter), "所有细胞", "免疫细胞"),
              ncol(sub_obj)))
  
  return(sub_obj)
}

# 创建分组对象
create_group_objects <- function(sample_list, group_name, cell_types = NULL) {
  objs <- list()
  for (s in names(sample_list)) {
    obj <- extract_cells(sample_list[[s]], group_name, cell_types)
    if (!is.null(obj)) objs[[s]] <- obj
  }
  
  if (length(objs) == 0) return(NULL)
  
  # 合并同一组的所有样本
  merged <- Reduce(merge, objs)
  merged <- JoinLayers(merged)
  return(merged)
}

# 准备肿瘤组、对照组、空白组的对象
prepare_objects <- function(cell_types = NULL) {
  cell_desc <- ifelse(is.null(cell_types), "all_cells", "immune_cells")
  
  # 肿瘤组（合并三个样本）
  tumor_objs <- all_objs[tumor_samples]
  tumor_merged <- create_group_objects(tumor_objs, "Tumor", cell_types)
  
  # 对照组
  control_obj <- extract_cells(all_objs[[control_sample]], "Control", cell_types)
  
  # 空白组
  blank_obj <- extract_cells(all_objs[[blank_sample]], "Blank", cell_types)
  
  return(list(
    tumor = tumor_merged,
    control = control_obj,
    blank = blank_obj,
    desc = cell_desc
  ))
}

# =============================================================================
# 6. Il6st 表达诊断函数
# =============================================================================
diagnose_il6st <- function(obj, group_name, cell_type_name) {
  if (is.null(obj)) return(NULL)
  if (!"Il6st" %in% rownames(obj)) {
    cat(sprintf("  ⚠️ %s - %s: Il6st 不在基因列表中\n", cell_type_name, group_name))
    return(NULL)
  }
  
  exp <- FetchData(obj, vars = "Il6st")$Il6st
  pct_expressed <- mean(exp > 0) * 100
  mean_exp <- mean(exp)
  median_exp <- median(exp)
  
  return(list(
    group = group_name,
    pct = pct_expressed,
    mean = mean_exp,
    median = median_exp,
    n_cells = length(exp),
    exp_vector = exp
  ))
}

# =============================================================================
# 7. DEG 分析函数
# =============================================================================
run_deg <- function(obj1, obj2, name1, name2, comp_name, cell_type_name) {
  cat(sprintf("\n【%s - %s】\n", cell_type_name, comp_name))
  
  if (is.null(obj1) || is.null(obj2)) {
    cat("  ✗ 对象为空，跳过\n")
    return(NULL)
  }
  
  cat(sprintf("  %s: %d 个细胞\n", name1, ncol(obj1)))
  cat(sprintf("  %s: %d 个细胞\n", name2, ncol(obj2)))
  
  # 合并对象
  merged <- merge(obj1, obj2)
  merged <- JoinLayers(merged)
  Idents(merged) <- merged$group
  
  # 运行 FindMarkers
  cat("  运行 FindMarkers (Wilcoxon)...\n")
  
  deg_results <- tryCatch({
    FindMarkers(merged,
                ident.1 = name1,
                ident.2 = name2,
                test.use = "wilcox",
                min.pct = MIN_PCT,
                logfc.threshold = 0,
                only.pos = FALSE,
                verbose = FALSE)
  }, error = function(e) {
    cat(sprintf("  ✗ 错误: %s\n", e$message))
    return(NULL)
  })
  
  if (is.null(deg_results) || nrow(deg_results) == 0) {
    cat("  ✗ 无结果返回\n")
    return(NULL)
  }
  
  # 整理结果
  deg_results$gene <- rownames(deg_results)
  deg_results <- deg_results %>%
    rename(pct1 = pct.1, pct2 = pct.2) %>%
    select(gene, avg_log2FC, pct1, pct2, p_val, p_val_adj) %>%
    arrange(p_val_adj, desc(abs(avg_log2FC)))
  
  # 标记显著性
  deg_results$sig <- case_when(
    deg_results$avg_log2FC > FC_THRESHOLD & deg_results$p_val_adj < PADJ_CUT ~ "Up",
    deg_results$avg_log2FC < -FC_THRESHOLD & deg_results$p_val_adj < PADJ_CUT ~ "Down",
    TRUE ~ "NS"
  )
  
  n_up <- sum(deg_results$sig == "Up")
  n_down <- sum(deg_results$sig == "Down")
  
  cat(sprintf("  ✓ 完成: 总基因=%d, Up=%d, Down=%d\n",
              nrow(deg_results), n_up, n_down))
  
  # Il6st 结果
  if ("Il6st" %in% deg_results$gene) {
    il6st_res <- deg_results[deg_results$gene == "Il6st", ]
    cat(sprintf("\n  ★ Il6st:\n"))
    cat(sprintf("    avg_log2FC = %.4f\n", il6st_res$avg_log2FC))
    cat(sprintf("    p_val      = %.2e\n", il6st_res$p_val))
    cat(sprintf("    p_val_adj  = %.2e\n", il6st_res$p_val_adj))
    cat(sprintf("    pct1(Tumor)= %.1f%%\n", il6st_res$pct1 * 100))
    cat(sprintf("    pct2(Cont) = %.1f%%\n", il6st_res$pct2 * 100))
    cat(sprintf("    显著性     = %s\n", il6st_res$sig))
  } else {
    cat(sprintf("\n  ⚠️ Il6st 不在DEG结果中！\n"))
  }
  
  return(deg_results)
}

# =============================================================================
# 8. 火山图函数
# =============================================================================
save_volcano <- function(res, comp_name, cell_type_name) {
  if (is.null(res) || nrow(res) == 0) return()
  
  res$log10p <- -log10(res$p_val_adj + 1e-300)
  
  # 标记显著基因
  top_genes <- bind_rows(
    res %>% filter(sig == "Up") %>% arrange(p_val_adj) %>% head(10),
    res %>% filter(sig == "Down") %>% arrange(p_val_adj) %>% head(10)
  )
  
  # 标记 Il6st
  res$is_il6st <- res$gene == "Il6st"
  
  p <- ggplot(res, aes(x = avg_log2FC, y = log10p, color = sig)) +
    geom_point(alpha = 0.5, size = 0.8) +
    scale_color_manual(values = c("Up" = "#E41A1C", "Down" = "#377EB8", 
                                   "NS" = "grey80")) +
    geom_vline(xintercept = c(-FC_THRESHOLD, FC_THRESHOLD), 
               linetype = "dashed", linewidth = 0.3) +
    geom_hline(yintercept = -log10(PADJ_CUT), 
               linetype = "dashed", linewidth = 0.3) +
    geom_text_repel(data = top_genes, aes(label = gene), 
                    size = 2.5, max.overlaps = 30, segment.size = 0.2) +
    geom_point(data = res %>% filter(is_il6st),
               aes(x = avg_log2FC, y = log10p),
               color = "#FF7F00", size = 3, shape = 17) +
    geom_text_repel(data = res %>% filter(is_il6st),
                    aes(label = gene), size = 4, color = "#FF7F00", 
                    fontface = "bold", nudge_x = 0.3) +
    labs(title = paste0(cell_type_name, ": ", comp_name),
         subtitle = sprintf("|log2FC| > %.2f, p_adj < %.2f | Up=%d, Down=%d",
                           FC_THRESHOLD, PADJ_CUT,
                           sum(res$sig == "Up"), sum(res$sig == "Down")),
         x = "avg_log2FC", y = "-log10(p_val_adj)") +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 8))
  
  file_prefix <- paste0(cell_type_name, "_", comp_name)
  ggsave(file.path(output_dir, paste0(file_prefix, "_volcano.pdf")), 
         p, width = 8, height = 6)
  ggsave(file.path(output_dir, paste0(file_prefix, "_volcano.png")), 
         p, width = 8, height = 6, dpi = 200)
  cat(sprintf("  ✓ 火山图: %s_volcano.pdf/png\n", file_prefix))
}

# =============================================================================
# 9. 执行分析 - 所有细胞
# =============================================================================
cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║                    分析1: 所有细胞 (All Cells)                       ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n")

all_cells <- prepare_objects(cell_types = NULL)

# 诊断 Il6st 表达
cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("Il6st 表达诊断 - 所有细胞\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

all_diag <- list()
all_diag[["Tumor"]] <- diagnose_il6st(all_cells$tumor, "Tumor", "所有细胞")
all_diag[["Control"]] <- diagnose_il6st(all_cells$control, "Control", "所有细胞")
all_diag[["Blank"]] <- diagnose_il6st(all_cells$blank, "Blank", "所有细胞")

cat("\n所有细胞 Il6st 检测率:\n")
for (g in names(all_diag)) {
  if (!is.null(all_diag[[g]])) {
    cat(sprintf("  %s: %.1f%% (%d 细胞)\n", 
                g, all_diag[[g]]$pct, all_diag[[g]]$n_cells))
  }
}

# DEG 分析 - 所有细胞
cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("DEG 分析 - 所有细胞\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")

all_results <- list()

# Tumor vs Control
all_results[["tumor_vs_control"]] <- run_deg(
  all_cells$tumor, all_cells$control, 
  "Tumor", "Control", "tumor_vs_control", "AllCells"
)

# Tumor vs Blank
all_results[["tumor_vs_blank"]] <- run_deg(
  all_cells$tumor, all_cells$blank,
  "Tumor", "Blank", "tumor_vs_blank", "AllCells"
)

# Blank vs Control
all_results[["blank_vs_control"]] <- run_deg(
  all_cells$blank, all_cells$control,
  "Blank", "Control", "blank_vs_control", "AllCells"
)

# 保存结果和火山图
for (comp in names(all_results)) {
  if (!is.null(all_results[[comp]])) {
    write.csv(all_results[[comp]], 
              file.path(output_dir, paste0("AllCells_", comp, ".csv")), 
              row.names = FALSE)
    save_volcano(all_results[[comp]], comp, "AllCells")
  }
}

# =============================================================================
# 10. 执行分析 - 免疫细胞
# =============================================================================
cat("\n\n")
cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║                  分析2: 免疫细胞 (Immune Cells)                      ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n")

immune_cells <- prepare_objects(cell_types = immune_cell_types)

# 诊断 Il6st 表达
cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("Il6st 表达诊断 - 免疫细胞\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

immune_diag <- list()
immune_diag[["Tumor"]] <- diagnose_il6st(immune_cells$tumor, "Tumor", "免疫细胞")
immune_diag[["Control"]] <- diagnose_il6st(immune_cells$control, "Control", "免疫细胞")
immune_diag[["Blank"]] <- diagnose_il6st(immune_cells$blank, "Blank", "免疫细胞")

cat("\n免疫细胞 Il6st 检测率:\n")
for (g in names(immune_diag)) {
  if (!is.null(immune_diag[[g]])) {
    cat(sprintf("  %s: %.1f%% (%d 细胞)\n", 
                g, immune_diag[[g]]$pct, immune_diag[[g]]$n_cells))
  }
}

# DEG 分析 - 免疫细胞
cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("DEG 分析 - 免疫细胞\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")

immune_results <- list()

# Tumor vs Control
immune_results[["tumor_vs_control"]] <- run_deg(
  immune_cells$tumor, immune_cells$control,
  "Tumor", "Control", "tumor_vs_control", "ImmuneCells"
)

# Tumor vs Blank
immune_results[["tumor_vs_blank"]] <- run_deg(
  immune_cells$tumor, immune_cells$blank,
  "Tumor", "Blank", "tumor_vs_blank", "ImmuneCells"
)

# Blank vs Control
immune_results[["blank_vs_control"]] <- run_deg(
  immune_cells$blank, immune_cells$control,
  "Blank", "Control", "blank_vs_control", "ImmuneCells"
)

# 保存结果和火山图
for (comp in names(immune_results)) {
  if (!is.null(immune_results[[comp]])) {
    write.csv(immune_results[[comp]], 
              file.path(output_dir, paste0("ImmuneCells_", comp, ".csv")), 
              row.names = FALSE)
    save_volcano(immune_results[[comp]], comp, "ImmuneCells")
  }
}

# =============================================================================
# 11. 汇总报告
# =============================================================================
cat("\n\n")
cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║                        汇总报告                                      ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n")

# 创建汇总表
summary_list <- list()

# 所有细胞
for (comp in names(all_results)) {
  if (!is.null(all_results[[comp]])) {
    res <- all_results[[comp]]
    il6st_entry <- res[res$gene == "Il6st", ]
    summary_list[[paste0("AllCells_", comp)]] <- data.frame(
      CellType = "AllCells",
      Comparison = comp,
      Total_genes = nrow(res),
      Up = sum(res$sig == "Up"),
      Down = sum(res$sig == "Down"),
      Il6st_FC = ifelse(nrow(il6st_entry) > 0, round(il6st_entry$avg_log2FC, 4), NA),
      Il6st_p_val = ifelse(nrow(il6st_entry) > 0, sprintf("%.2e", il6st_entry$p_val), NA),
      Il6st_p_adj = ifelse(nrow(il6st_entry) > 0, sprintf("%.2e", il6st_entry$p_val_adj), NA),
      Il6st_sig = ifelse(nrow(il6st_entry) > 0, il6st_entry$sig, "Not Found"),
      Il6st_pct_tumor = ifelse(nrow(il6st_entry) > 0, sprintf("%.1f%%", il6st_entry$pct1 * 100), NA),
      Il6st_pct_control = ifelse(nrow(il6st_entry) > 0, sprintf("%.1f%%", il6st_entry$pct2 * 100), NA),
      stringsAsFactors = FALSE
    )
  }
}

# 免疫细胞
for (comp in names(immune_results)) {
  if (!is.null(immune_results[[comp]])) {
    res <- immune_results[[comp]]
    il6st_entry <- res[res$gene == "Il6st", ]
    summary_list[[paste0("ImmuneCells_", comp)]] <- data.frame(
      CellType = "ImmuneCells",
      Comparison = comp,
      Total_genes = nrow(res),
      Up = sum(res$sig == "Up"),
      Down = sum(res$sig == "Down"),
      Il6st_FC = ifelse(nrow(il6st_entry) > 0, round(il6st_entry$avg_log2FC, 4), NA),
      Il6st_p_val = ifelse(nrow(il6st_entry) > 0, sprintf("%.2e", il6st_entry$p_val), NA),
      Il6st_p_adj = ifelse(nrow(il6st_entry) > 0, sprintf("%.2e", il6st_entry$p_val_adj), NA),
      Il6st_sig = ifelse(nrow(il6st_entry) > 0, il6st_entry$sig, "Not Found"),
      Il6st_pct_tumor = ifelse(nrow(il6st_entry) > 0, sprintf("%.1f%%", il6st_entry$pct1 * 100), NA),
      Il6st_pct_control = ifelse(nrow(il6st_entry) > 0, sprintf("%.1f%%", il6st_entry$pct2 * 100), NA),
      stringsAsFactors = FALSE
    )
  }
}

summary_df <- do.call(rbind, summary_list)
write.csv(summary_df, file.path(output_dir, "summary.csv"), row.names = FALSE)

cat("\n【Il6st 结果汇总】\n")
cat("──────────────────────────────────────────────────────────────────────\n")
print(summary_df[, c("CellType", "Comparison", "Il6st_FC", "Il6st_p_adj", "Il6st_sig", 
                     "Il6st_pct_tumor", "Il6st_pct_control")])
cat("──────────────────────────────────────────────────────────────────────\n")

# =============================================================================
# 12. Il6st 表达小提琴图
# =============================================================================
cat("\n生成 Il6st 表达分布图...\n")

# 收集所有细胞的 Il6st 表达数据
plot_data_all <- data.frame()
for (g in c("Tumor", "Control", "Blank")) {
  if (!is.null(all_diag[[g]])) {
    tmp <- data.frame(
      Expression = all_diag[[g]]$exp_vector,
      Group = g,
      CellType = "All Cells"
    )
    plot_data_all <- rbind(plot_data_all, tmp)
  }
}

# 收集免疫细胞的 Il6st 表达数据
plot_data_immune <- data.frame()
for (g in c("Tumor", "Control", "Blank")) {
  if (!is.null(immune_diag[[g]])) {
    tmp <- data.frame(
      Expression = immune_diag[[g]]$exp_vector,
      Group = g,
      CellType = "Immune Cells"
    )
    plot_data_immune <- rbind(plot_data_immune, tmp)
  }
}

plot_data <- rbind(plot_data_all, plot_data_immune)
plot_data$Group <- factor(plot_data$Group, levels = c("Tumor", "Control", "Blank"))

p_combined <- ggplot(plot_data, aes(x = Group, y = Expression, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.2, outlier.size = 0.3) +
  facet_wrap(~CellType, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = c("Tumor" = "#E41A1C", "Control" = "#377EB8", 
                                "Blank" = "#4DAF4A")) +
  labs(title = "Il6st Expression Across Groups",
       subtitle = paste0("Tumor (n=", ncol(all_cells$tumor), " cells), ",
                         "Control (n=", ncol(all_cells$control), "), ",
                         "Blank (n=", ncol(all_cells$blank), ")"),
       y = "Il6st Expression Level", x = NULL) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none",
        strip.background = element_rect(fill = "grey90"),
        strip.text = element_text(face = "bold"))

ggsave(file.path(output_dir, "Il6st_expression_across_groups.pdf"), 
       p_combined, width = 8, height = 5)
ggsave(file.path(output_dir, "Il6st_expression_across_groups.png"), 
       p_combined, width = 8, height = 5, dpi = 200)

cat("  ✓ 已保存: Il6st_expression_across_groups.pdf/png\n")

# =============================================================================
# 13. 完成
# =============================================================================
cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║                        分析完成！                                    ║\n")
cat("╠══════════════════════════════════════════════════════════════════════╣\n")
cat(sprintf("║  输出目录: %-50s║\n", output_dir))
cat("╠══════════════════════════════════════════════════════════════════════╣\n")
cat("║  生成文件:                                                          ║\n")
cat("║    - AllCells_*.csv              所有细胞DEG结果                     ║\n")
cat("║    - ImmuneCells_*.csv           免疫细胞DEG结果                    ║\n")
cat("║    - *_volcano.pdf/png           火山图                              ║\n")
cat("║    - summary.csv                 汇总表                              ║\n")
cat("║    - Il6st_expression_across_groups.pdf/png  Il6st表达分布          ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n")
