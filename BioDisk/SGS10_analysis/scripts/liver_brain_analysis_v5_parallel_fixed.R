#!/usr/bin/env Rscript
# =============================================================================
# 文件名: liver_brain_analysis_v5_parallel_fixed.R
# 用途:   肝脑肿瘤分析 - 稳态小胶质细胞 + 完整富集分析（并行修复版）
# =============================================================================

cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║   肝脑肿瘤分析 v5 - 稳态小胶质细胞 + 并行加速（修复版）              ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n\n")

# ==================== 0. 环境检查和包加载 ====================
cat("========================================\n")
cat("检查R环境和依赖包...\n")
cat("========================================\n")

.libPaths("/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/lib/R/library")
cat("R库路径:", .libPaths()[1], "\n")

required_packages <- c(
  "Seurat", "parallel", "future", "ggplot2", "ggrepel", "dplyr", "tidyr",
  "clusterProfiler", "org.Mm.eg.db", "enrichplot"
)

for (pkg in required_packages) {
  if (!suppressPackageStartupMessages(require(pkg, character.only = TRUE, quietly = TRUE))) {
    cat("  ✗", pkg, "未安装\n")
    quit(status = 1)
  }
  cat("  ✓", pkg, "已加载\n")
}

suppressPackageStartupMessages({
  library(Seurat)
  library(parallel)
  library(future)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tidyr)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
})

cat("✓ 所有必需包已加载\n\n")

# ==================== 并行核心数配置 ====================
available_cores <- detectCores()
cat("检测到可用核心数:", available_cores, "\n")

CORES_FINDMARKERS <- min(16, available_cores)  # FindMarkers内部并行
CORES_TASK        <- min(8, available_cores)   # 任务级并行
CORES_ENRICH      <- min(4, available_cores)   # 富集分析并行

cat("使用核心数配置:\n")
cat("  FindMarkers并行:", CORES_FINDMARKERS, "核\n")
cat("  任务级并行:     ", CORES_TASK, "核\n")
cat("  富集分析并行:   ", CORES_ENRICH, "核\n\n")

options(future.globals.maxSize = 32 * 1024^3)
options(future.rng.onMisuse = "ignore")

# ==================== 工作目录 ====================
setwd("/mnt/BioDisk/wangmingyuan/SGS10_analysis")

output_dir <- "liver_brain_v5_homeostatic_results"
dirs <- c(
  output_dir,
  file.path(output_dir, "01_liver_DEG"),
  file.path(output_dir, "02_brain_DEG"),
  file.path(output_dir, "03_liver_enrichment"),
  file.path(output_dir, "04_brain_enrichment"),
  file.path(output_dir, "05_crosstalk"),
  file.path(output_dir, "06_plots/volcano"),
  file.path(output_dir, "07_summary"),
  file.path(output_dir, "08_enrichment_detailed")
)
for (d in dirs) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

# ==================== DEG参数 ====================
FC_THRESHOLD <- 0.20
PADJ_CUT     <- 0.10
MIN_PCT      <- 0.10
MAX_CELLS    <- 500
MIN_CELLS    <- 5

# 细胞类型 - 关键修改：将 Microglia 改为 Microglia_Homeostatic
liver_cell_types <- c("All_Cells", "Hepatocytes", "Endothelial", "Kupffer")
brain_cell_types <- c("All_Cells", "Endothelial", "Microglia_Homeostatic")

cat("========================================\n")
cat("肝脑肿瘤分析 v5（稳态小胶质细胞 + 并行）\n")
cat("========================================\n")
cat(sprintf("  脑部细胞类型: %s\n", paste(brain_cell_types, collapse=", ")))
cat(sprintf("  avg_log2FC:   > %.2f\n", FC_THRESHOLD))
cat(sprintf("  p_val_adj:    < %.2f\n\n", PADJ_CUT))

# ==================== 1. 加载样本 ====================
cat("========== 1. 加载样本 ==========\n")

samples_to_load <- c("gan56","liver3","liver8",
                     "nao56","brain3","brain8",
                     "con_brain","bla_brain","bla_liver","con_liver")

seurat_objs <- list()
for (s in samples_to_load) {
  fp <- file.path(s, paste0(s, "_manual_annotated.rds"))
  if (!file.exists(fp)) {
    cat(sprintf("  ✗ %-15s  文件不存在\n", s))
    next
  }
  obj <- tryCatch({
    readRDS(fp)
  }, error = function(e) {
    cat(sprintf("  ✗ %-15s  读取失败: %s\n", s, e$message))
    return(NULL)
  })
  
  if (is.null(obj)) next
  
  obj <- JoinLayers(obj)
  obj <- RenameCells(obj, add.cell.id = s)
  cat(sprintf("  ✓ %-15s  细胞数: %d\n", s, ncol(obj)))
  seurat_objs[[s]] <- obj
}
cat(sprintf("  成功加载: %d / %d 个样本\n\n", length(seurat_objs), length(samples_to_load)))

# ==================== 2. 辅助函数 ====================

extract_cells <- function(obj, cell_type = "All_Cells") {
  if (cell_type != "All_Cells" && "final_celltype" %in% colnames(obj@meta.data)) {
    return(WhichCells(obj, expression = final_celltype == cell_type))
  }
  return(colnames(obj))
}

merge_two_groups <- function(obj1, obj2, cell_type, label1, label2) {
  cells1 <- extract_cells(obj1, cell_type)
  cells2 <- extract_cells(obj2, cell_type)
  
  if (length(cells1) < MIN_CELLS || length(cells2) < MIN_CELLS) return(NULL)
  
  s1 <- subset(obj1, cells = cells1)
  s1$group <- label1
  
  s2 <- subset(obj2, cells = cells2)
  s2$group <- label2
  
  merged <- merge(s1, s2)
  merged <- JoinLayers(merged)
  return(merged)
}

merge_to_group <- function(obj_list, group_label, cell_type) {
  sub_list <- list()
  for (obj in obj_list) {
    cells <- extract_cells(obj, cell_type)
    if (length(cells) >= MIN_CELLS) {
      s <- subset(obj, cells = cells)
      s$group <- group_label
      sub_list[[length(sub_list) + 1]] <- s
    }
  }
  
  if (length(sub_list) == 0) return(NULL)
  
  merged <- Reduce(merge, sub_list)
  merged <- JoinLayers(merged)
  return(merged)
}

# ==================== 3. 完整富集分析函数 ====================
run_complete_enrichment <- function(deg_res, key, enrich_dir) {
  
  detailed_dir <- file.path(enrich_dir, "detailed")
  if (!dir.exists(detailed_dir)) dir.create(detailed_dir, recursive = TRUE)
  
  sig_up <- deg_res %>% filter(sig == "Up") %>% pull(gene)
  sig_down <- deg_res %>% filter(sig == "Down") %>% pull(gene)
  sig_all <- deg_res %>% filter(sig != "NS") %>% pull(gene)
  
  if (length(sig_all) < 5) {
    return(NULL)
  }
  
  # 转换基因ID
  gene_df_all <- tryCatch(
    bitr(sig_all, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Mm.eg.db),
    error = function(e) NULL
  )
  
  if (is.null(gene_df_all) || nrow(gene_df_all) < 3) return(NULL)
  
  entrez_all <- unique(gene_df_all$ENTREZID)
  
  result <- list(key = key, n_genes = length(sig_all), 
                 n_up = length(sig_up), n_down = length(sig_down),
                 kegg_n = 0, go_n = 0)
  
  # KEGG 富集
  kegg_all <- tryCatch(
    enrichKEGG(gene = entrez_all, organism = "mmu", pvalueCutoff = 0.05),
    error = function(e) NULL
  )
  
  if (!is.null(kegg_all) && nrow(kegg_all@result) > 0) {
    write.csv(as.data.frame(kegg_all),
              file.path(detailed_dir, paste0(key, "_KEGG.csv")), row.names = FALSE)
    result$kegg_n <- sum(kegg_all@result$p.adjust < 0.05, na.rm = TRUE)
    
    if (result$kegg_n > 0) {
      p <- dotplot(kegg_all, showCategory = 20, title = paste0("KEGG: ", key))
      ggsave(file.path(detailed_dir, paste0(key, "_KEGG_bubble.pdf")), p, width = 10, height = 8)
    }
  }
  
  # GO-BP 富集
  go_all <- tryCatch(
    enrichGO(gene = entrez_all, OrgDb = org.Mm.eg.db, ont = "BP",
             pvalueCutoff = 0.05, readable = TRUE),
    error = function(e) NULL
  )
  
  if (!is.null(go_all) && nrow(go_all@result) > 0) {
    write.csv(as.data.frame(go_all),
              file.path(detailed_dir, paste0(key, "_GO_BP.csv")), row.names = FALSE)
    result$go_n <- sum(go_all@result$p.adjust < 0.05, na.rm = TRUE)
    
    if (result$go_n > 0) {
      p <- dotplot(go_all, showCategory = 20, title = paste0("GO-BP: ", key))
      ggsave(file.path(detailed_dir, paste0(key, "_GO_BP_bubble.pdf")), p, width = 10, height = 8)
    }
  }
  
  return(result)
}

# ==================== 4. 单任务DEG函数 ====================
run_one_deg <- function(merged_obj, comp_name, group1, group2, ct_name, key) {

  if (is.null(merged_obj)) return(NULL)
  
  n1 <- sum(merged_obj$group == group1)
  n2 <- sum(merged_obj$group == group2)
  
  if (n1 < MIN_CELLS || n2 < MIN_CELLS) return(NULL)

  # DownSample
  set.seed(42)
  cg1 <- colnames(merged_obj)[merged_obj$group == group1]
  cg2 <- colnames(merged_obj)[merged_obj$group == group2]
  
  if (length(cg1) > MAX_CELLS) cg1 <- sample(cg1, MAX_CELLS)
  if (length(cg2) > MAX_CELLS) cg2 <- sample(cg2, MAX_CELLS)
  
  obj_sub <- subset(merged_obj, cells = c(cg1, cg2))
  Idents(obj_sub) <- obj_sub$group

  # 设置并行
  plan("multicore", workers = CORES_FINDMARKERS)

  res <- tryCatch({
    FindMarkers(obj_sub,
                ident.1 = group1,
                ident.2 = group2,
                test.use = "wilcox",
                min.pct = MIN_PCT,
                logfc.threshold = 0,
                only.pos = FALSE,
                verbose = FALSE)
  }, error = function(e) {
    return(NULL)
  })

  plan("sequential")

  if (is.null(res) || nrow(res) == 0) return(NULL)

  res$gene <- rownames(res)
  res <- res %>%
    select(gene, avg_log2FC, pct.1, pct.2, p_val, p_val_adj) %>%
    rename(pct1 = pct.1, pct2 = pct.2) %>%
    arrange(p_val_adj, desc(abs(avg_log2FC)))

  res$sig <- ifelse(
    abs(res$avg_log2FC) > FC_THRESHOLD & res$p_val_adj < PADJ_CUT,
    ifelse(res$avg_log2FC > 0, "Up", "Down"), "NS"
  )

  cat(sprintf("  ✓ %-50s Up=%-4d Down=%d\n", key, sum(res$sig=="Up"), sum(res$sig=="Down")))
  return(res)
}

# ==================== 5. 构建所有任务 ====================
cat("\n========== 2. 构建任务列表 ==========\n")

build_all_tasks <- function() {
  tasks <- list()
  
  # 肝脏任务
  liver_comparisons <- list(
    list(comp="liver_tumor_vs_blank", g1_samples=c("gan56","liver3","liver8"), 
         g1_label="Tumor", g2_obj="bla_liver", g2_label="Blank"),
    list(comp="liver_tumor_vs_control", g1_samples=c("gan56","liver3","liver8"), 
         g1_label="Tumor", g2_obj="con_liver", g2_label="Control"),
    list(comp="liver_blank_vs_control", g1_obj="bla_liver", g1_label="Blank", 
         g2_obj="con_liver", g2_label="Control")
  )
  
  for (comp in liver_comparisons) {
    for (ct_name in liver_cell_types) {
      tasks[[length(tasks) + 1]] <- list(
        type = "liver",
        comp = comp,
        ct_name = ct_name,
        key = paste0(comp$comp, "_", ct_name)
      )
    }
  }
  
  # 脑部任务
  brain_comparisons <- list(
    list(comp="brain_tumor_vs_blank", g1_samples=c("nao56","brain3","brain8"), 
         g1_label="TumorBrain", g2_obj="bla_brain", g2_label="Blank"),
    list(comp="brain_tumor_vs_control", g1_samples=c("nao56","brain3","brain8"), 
         g1_label="TumorBrain", g2_obj="con_brain", g2_label="Control"),
    list(comp="brain_blank_vs_control", g1_obj="bla_brain", g1_label="Blank", 
         g2_obj="con_brain", g2_label="Control")
  )
  
  for (comp in brain_comparisons) {
    for (ct_name in brain_cell_types) {
      tasks[[length(tasks) + 1]] <- list(
        type = "brain",
        comp = comp,
        ct_name = ct_name,
        key = paste0(comp$comp, "_", ct_name)
      )
    }
  }
  
  return(tasks)
}

# ==================== 6. 准备合并对象函数 ====================
prepare_merged_object <- function(task, seurat_objs) {
  
  comp <- task$comp
  ct_name <- task$ct_name
  
  if (!is.null(comp$g1_samples)) {
    # 多样本组 vs 单样本组
    g1_objs <- list()
    for (s in comp$g1_samples) {
      if (!is.null(seurat_objs[[s]])) {
        g1_objs[[length(g1_objs) + 1]] <- seurat_objs[[s]]
      }
    }
    if (length(g1_objs) == 0) return(NULL)
    
    g1_merged <- merge_to_group(g1_objs, comp$g1_label, ct_name)
    if (is.null(g1_merged)) return(NULL)
    
    g2_obj <- seurat_objs[[comp$g2_obj]]
    if (is.null(g2_obj)) return(NULL)
    
    g2_cells <- extract_cells(g2_obj, ct_name)
    if (length(g2_cells) < MIN_CELLS) return(NULL)
    
    g2_sub <- subset(g2_obj, cells = g2_cells)
    g2_sub$group <- comp$g2_label
    
    merged <- merge(g1_merged, g2_sub)
    merged <- JoinLayers(merged)
    
  } else {
    # 单样本 vs 单样本
    merged <- merge_two_groups(
      seurat_objs[[comp$g1_obj]], seurat_objs[[comp$g2_obj]],
      ct_name, comp$g1_label, comp$g2_label
    )
  }
  
  return(merged)
}

# ==================== 7. 并行执行DEG ====================
cat("\n========== 3. 并行执行DEG ==========\n")

all_tasks <- build_all_tasks()
cat(sprintf("  总任务数: %d\n", length(all_tasks)))

# 设置并行任务数
n_cores <- min(CORES_TASK, length(all_tasks))
cat(sprintf("  使用 %d 个核心并行执行\n", n_cores))

t_start <- proc.time()

# 使用 mclapply 并行执行，正确捕获结果
results_list <- mclapply(seq_along(all_tasks), function(i) {
  task <- all_tasks[[i]]
  
  # 准备合并对象
  merged <- prepare_merged_object(task, seurat_objs)
  
  if (is.null(merged)) {
    return(list(key = task$key, success = FALSE, error = "merged_null"))
  }
  
  # 运行DEG
  res <- run_one_deg(merged, task$comp$comp, task$comp$g1_label, 
                     task$comp$g2_label, task$ct_name, task$key)
  
  if (is.null(res)) {
    return(list(key = task$key, success = FALSE, error = "deg_failed"))
  }
  
  # 保存结果
  deg_dir <- ifelse(task$type == "liver", 
                    file.path(output_dir, "01_liver_DEG"),
                    file.path(output_dir, "02_brain_DEG"))
  write.csv(res, file.path(deg_dir, paste0(task$key, ".csv")), row.names = FALSE)
  
  # 保存火山图
  tryCatch({
    res_v <- res
    res_v$log10p <- -log10(res_v$p_val_adj + 1e-300)
    top_genes <- bind_rows(
      res_v %>% filter(sig=="Up")   %>% arrange(p_val_adj) %>% head(10),
      res_v %>% filter(sig=="Down") %>% arrange(p_val_adj) %>% head(10)
    )
    n_up <- sum(res_v$sig == "Up")
    n_down <- sum(res_v$sig == "Down")
    p_v <- ggplot(res_v, aes(x=avg_log2FC, y=log10p, color=sig)) +
      geom_point(alpha=0.5, size=0.8) +
      scale_color_manual(values=c("Up"="#E41A1C","Down"="#377EB8","NS"="grey75")) +
      geom_vline(xintercept=c(-FC_THRESHOLD, FC_THRESHOLD), linetype="dashed") +
      geom_hline(yintercept=-log10(PADJ_CUT), linetype="dashed") +
      geom_text_repel(data=top_genes, aes(label=gene), size=2.5, max.overlaps=30) +
      labs(title = paste0(task$comp$comp, " [", task$ct_name, "]"),
           subtitle = sprintf("Up=%d, Down=%d", n_up, n_down),
           x="avg_log2FC", y="-log10(p_val_adj)") +
      theme_bw()
    
    volcano_dir <- file.path(output_dir, "06_plots/volcano")
    ggsave(file.path(volcano_dir, paste0(task$key, "_volcano.pdf")), p_v, width=8, height=6)
  }, error = function(e) NULL)
  
  return(list(key = task$key, success = TRUE, result = res))
  
}, mc.cores = n_cores, mc.preschedule = FALSE, mc.cleanup = TRUE)

t_elapsed <- proc.time() - t_start
cat(sprintf("\n  ✓ DEG完成，耗时 %.1f 秒\n", t_elapsed["elapsed"]))

# 收集成功的结果
all_results <- list()
for (r in results_list) {
  if (is.list(r) && r$success && !is.null(r$result)) {
    all_results[[r$key]] <- r$result
  }
}

cat(sprintf("  成功完成: %d / %d 个任务\n", length(all_results), length(all_tasks)))

# ==================== 8. 并行富集分析 ====================
cat("\n\n========== 4. 并行完整富集分析 ==========\n")

if (length(all_results) > 0) {
  
  # 准备富集任务
  enrich_tasks <- list()
  for (key in names(all_results)) {
    enrich_tasks[[length(enrich_tasks) + 1]] <- list(
      key = key,
      deg_res = all_results[[key]],
      enrich_dir = ifelse(grepl("liver", key),
                          file.path(output_dir, "03_liver_enrichment"),
                          file.path(output_dir, "04_brain_enrichment"))
    )
  }
  
  n_enrich_cores <- min(CORES_ENRICH, length(enrich_tasks))
  cat(sprintf("  使用 %d 个核心并行富集分析\n", n_enrich_cores))
  
  enrich_results_list <- mclapply(enrich_tasks, function(task) {
    res <- run_complete_enrichment(task$deg_res, task$key, task$enrich_dir)
    return(res)
  }, mc.cores = n_enrich_cores, mc.preschedule = FALSE)
  
  # 收集富集结果
  enrich_results <- list()
  for (r in enrich_results_list) {
    if (!is.null(r)) {
      enrich_results[[r$key]] <- r
    }
  }
  
  cat(sprintf("  富集完成: %d / %d 组有显著结果\n", length(enrich_results), length(enrich_tasks)))
}

# ==================== 9. 生成汇总报告 ====================
cat("\n\n========== 5. 生成汇总报告 ==========\n")

if (length(all_results) > 0) {
  deg_summary <- do.call(rbind, lapply(names(all_results), function(nm) {
    res <- all_results[[nm]]
    
    if ("Il6st" %in% res$gene) {
      il6st <- res[res$gene == "Il6st", ]
      il6st_fc <- il6st$avg_log2FC
      il6st_padj <- il6st$p_val_adj
      il6st_sig <- il6st$sig
    } else {
      il6st_fc <- NA
      il6st_padj <- NA
      il6st_sig <- "Not Found"
    }
    
    data.frame(
      Comparison = nm,
      CellType = ifelse(grepl("Microglia_Homeostatic", nm), "Microglia_Homeostatic",
                       ifelse(grepl("Endothelial", nm), "Endothelial", "Other")),
      Up = sum(res$sig=="Up", na.rm=TRUE),
      Down = sum(res$sig=="Down", na.rm=TRUE),
      Total = sum(res$sig!="NS", na.rm=TRUE),
      Il6st_FC = round(il6st_fc, 4),
      Il6st_padj = sprintf("%.2e", il6st_padj),
      Il6st_sig = il6st_sig,
      stringsAsFactors = FALSE
    )
  }))
  
  write.csv(deg_summary, file.path(output_dir, "07_summary", "DEG_count_summary.csv"), row.names = FALSE)
  cat("✓ DEG汇总已保存\n")
  
  # 打印 Il6st 结果
  cat("\n【Il6st 结果汇总】\n")
  cat("────────────────────────────────────────────────────────────────────────────\n")
  il6st_summary <- deg_summary[!is.na(deg_summary$Il6st_FC) & deg_summary$Il6st_sig != "Not Found", 
                                c("Comparison", "CellType", "Il6st_FC", "Il6st_padj", "Il6st_sig")]
  if (nrow(il6st_summary) > 0) {
    print(il6st_summary)
  } else {
    cat("  未检测到 Il6st 显著差异\n")
  }
  cat("────────────────────────────────────────────────────────────────────────────\n")
}

# 富集结果汇总
if (exists("enrich_results") && length(enrich_results) > 0) {
  enrich_summary <- do.call(rbind, lapply(names(enrich_results), function(nm) {
    er <- enrich_results[[nm]]
    data.frame(
      Comparison = nm,
      Total_Genes = er$n_genes,
      Up_Genes = er$n_up,
      Down_Genes = er$n_down,
      KEGG_Pathways = er$kegg_n,
      GO_BP_Terms = er$go_n,
      stringsAsFactors = FALSE
    )
  }))
  write.csv(enrich_summary, file.path(output_dir, "08_enrichment_detailed", "enrichment_summary.csv"), row.names = FALSE)
  cat("✓ 富集汇总已保存\n")
}

# 生成分析报告
report_path <- file.path(output_dir, "07_summary", "analysis_report.txt")
sink(report_path)
cat("========================================\n")
cat("肝脑肿瘤分析报告 v5（稳态小胶质细胞 + 并行）\n")
cat("========================================\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("【并行配置】\n")
cat(sprintf("  可用核心数:       %d\n", available_cores))
cat(sprintf("  FindMarkers并行:  %d 核\n", CORES_FINDMARKERS))
cat(sprintf("  任务级并行:       %d 核\n", CORES_TASK))
cat(sprintf("  富集分析并行:     %d 核\n", CORES_ENRICH))
cat(sprintf("  实际使用核心:     %d 核\n", n_cores))
cat(sprintf("  总耗时:           %.1f 秒\n\n", t_elapsed["elapsed"]))
cat("【DEG参数】\n")
cat(sprintf("  avg_log2FC: > %.2f\n", FC_THRESHOLD))
cat(sprintf("  p_val_adj:  < %.2f\n", PADJ_CUT))
cat(sprintf("  min.pct:    %.2f\n", MIN_PCT))
cat(sprintf("  DownSample: 每组最多 %d 个细胞\n\n", MAX_CELLS))
cat("【细胞类型配置】\n")
cat(sprintf("  肝脏: %s\n", paste(liver_cell_types, collapse=", ")))
cat(sprintf("  脑部: %s\n\n", paste(brain_cell_types, collapse=", ")))
cat("【任务统计】\n")
cat(sprintf("  总任务数:       %d\n", length(all_tasks)))
cat(sprintf("  成功完成:       %d\n", length(all_results)))
cat(sprintf("  富集成功:       %d\n\n", ifelse(exists("enrich_results"), length(enrich_results), 0)))
cat("【DEG数量统计】\n")
if (exists("deg_summary") && nrow(deg_summary) > 0) print(deg_summary) else cat("  无显著结果\n")
sink()

cat("\n✓ 报告已保存:", report_path, "\n")
cat("\n========================================\n")
cat("分析完成！\n")
cat(sprintf("  DEG组数: %d\n", length(all_results)))
cat(sprintf("  富集组数: %d\n", ifelse(exists("enrich_results"), length(enrich_results), 0)))
cat(sprintf("  总耗时: %.1f 秒\n", t_elapsed["elapsed"]))
cat(sprintf("  输出目录: %s\n", output_dir))
cat("========================================\n")
