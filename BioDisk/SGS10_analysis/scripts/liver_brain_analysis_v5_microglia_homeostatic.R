#!/usr/bin/env Rscript
# =============================================================================
# 文件名: liver_brain_analysis_v5_microglia_homeostatic.R
# 用途:   肝脑肿瘤分析 - 使用稳态小胶质细胞 + 完整富集分析
# 用法:   conda activate R_4.5.2 && Rscript liver_brain_analysis_v5_microglia_homeostatic.R
# =============================================================================

cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║     肝脑肿瘤分析 v5 - 稳态小胶质细胞 + 完整富集分析                  ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n\n")

# ==================== 0. 环境检查和包加载 ====================
cat("========================================\n")
cat("检查R环境和依赖包...\n")
cat("========================================\n")

# 设置库路径
.libPaths("/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/lib/R/library")
cat("R库路径:", .libPaths()[1], "\n")

# 定义需要的包
required_packages <- c(
  "Seurat", "parallel", "future", "ggplot2", "ggrepel", "dplyr", "tidyr",
  "ComplexHeatmap", "circlize", "clusterProfiler", "org.Mm.eg.db", 
  "enrichplot", "VennDiagram", "grid", "pheatmap", "RColorBrewer"
)

# 检查包是否可用
missing_packages <- c()
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    missing_packages <- c(missing_packages, pkg)
    cat("  ✗", pkg, "未安装\n")
  } else {
    cat("  ✓", pkg, "已加载\n")
  }
}

if (length(missing_packages) > 0) {
  cat("\n⚠️ 缺少以下包，请先安装:\n")
  cat(paste("  ", missing_packages, collapse = "\n"), "\n")
  quit(status = 1)
}

cat("✓ 所有必需包已加载\n\n")

suppressPackageStartupMessages({
  library(Seurat)
  library(parallel)
  library(future)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tidyr)
  library(ComplexHeatmap)
  library(circlize)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
  library(VennDiagram)
  library(grid)
  library(pheatmap)
  library(RColorBrewer)
})

# ==================== 并行核心数配置 ====================
available_cores <- detectCores()
cat("检测到可用核心数:", available_cores, "\n")

CORES_FINDMARKERS <- min(32, available_cores)
CORES_TASK        <- min(8, available_cores)
CORES_ENRICH      <- min(4, available_cores)

cat("使用核心数配置:\n")
cat("  FindMarkers并行:", CORES_FINDMARKERS, "核\n")
cat("  任务级并行:     ", CORES_TASK, "核\n")
cat("  富集分析并行:   ", CORES_ENRICH, "核\n\n")

options(future.globals.maxSize = 32 * 1024^3)

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
  file.path(output_dir, "06_plots/heatmap"),
  file.path(output_dir, "06_plots/venn"),
  file.path(output_dir, "06_plots/bubble"),
  file.path(output_dir, "06_plots/IL6_mechanism"),
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
brain_cell_types <- c("All_Cells", "Endothelial", "Microglia_Homeostatic")  # ← 修改此处

# 配对样本
paired_samples <- list(
  list(liver = "gan56",  brain = "nao56",  id = "mouse56"),
  list(liver = "liver3", brain = "brain3", id = "mouse3"),
  list(liver = "liver8", brain = "brain8", id = "mouse8")
)

cat("========================================\n")
cat("肝脑肿瘤分析 v5（稳态小胶质细胞版）\n")
cat("========================================\n")
cat(sprintf("  脑部细胞类型: %s\n", paste(brain_cell_types, collapse=", ")))
cat(sprintf("  FindMarkers并行: %d 核\n", CORES_FINDMARKERS))
cat(sprintf("  任务级并行:      %d 核\n", CORES_TASK))
cat(sprintf("  avg_log2FC:      > %.2f\n", FC_THRESHOLD))
cat(sprintf("  p_val_adj:       < %.2f\n", PADJ_CUT))
cat(sprintf("  DownSample:      每组最多 %d 个细胞\n\n", MAX_CELLS))

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

merge_to_group <- function(obj_list, group_label, cell_type) {
  sub_list <- lapply(obj_list, function(obj) {
    cells <- extract_cells(obj, cell_type)
    if (length(cells) < MIN_CELLS) return(NULL)
    s <- subset(obj, cells = cells)
    s$group <- group_label
    return(s)
  })
  sub_list <- Filter(Negate(is.null), sub_list)
  if (length(sub_list) == 0) return(NULL)
  merged <- Reduce(merge, sub_list)
  merged <- JoinLayers(merged)
  return(merged)
}

merge_two_groups <- function(obj1, obj2, cell_type, label1, label2) {
  cells1 <- extract_cells(obj1, cell_type)
  cells2 <- extract_cells(obj2, cell_type)
  if (length(cells1) < MIN_CELLS || length(cells2) < MIN_CELLS) return(NULL)
  s1 <- subset(obj1, cells = cells1); s1$group <- label1
  s2 <- subset(obj2, cells = cells2); s2$group <- label2
  merged <- merge(s1, s2)
  merged <- JoinLayers(merged)
  return(merged)
}

# ==================== 3. 完整富集分析函数 ====================
run_complete_enrichment <- function(deg_res, key, enrich_dir, organ_tag) {
  
  # 创建详细富集目录
  detailed_dir <- file.path(enrich_dir, "detailed")
  if (!dir.exists(detailed_dir)) dir.create(detailed_dir, recursive = TRUE)
  
  # 提取显著基因
  sig_up <- deg_res %>% filter(sig == "Up") %>% pull(gene)
  sig_down <- deg_res %>% filter(sig == "Down") %>% pull(gene)
  sig_all <- deg_res %>% filter(sig != "NS") %>% pull(gene)
  
  if (length(sig_all) < 5) {
    cat(sprintf("  ⚠️ [%s] 显著基因不足5个，跳过富集\n", key))
    return(NULL)
  }
  
  # 转换基因ID
  gene_df_all <- tryCatch(
    bitr(sig_all, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Mm.eg.db),
    error = function(e) NULL
  )
  
  if (is.null(gene_df_all) || nrow(gene_df_all) < 3) {
    cat(sprintf("  ⚠️ [%s] 基因ID转换失败\n", key))
    return(NULL)
  }
  
  entrez_all <- unique(gene_df_all$ENTREZID)
  entrez_up <- NULL
  entrez_down <- NULL
  
  if (length(sig_up) >= 3) {
    gene_df_up <- bitr(sig_up, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Mm.eg.db)
    if (!is.null(gene_df_up)) entrez_up <- unique(gene_df_up$ENTREZID)
  }
  if (length(sig_down) >= 3) {
    gene_df_down <- bitr(sig_down, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Mm.eg.db)
    if (!is.null(gene_df_down)) entrez_down <- unique(gene_df_down$ENTREZID)
  }
  
  # ===== KEGG 富集 =====
  cat(sprintf("  [%s] 运行KEGG富集...\n", key))
  kegg_all <- tryCatch(
    enrichKEGG(gene = entrez_all, organism = "mmu", 
               pvalueCutoff = 0.05, qvalueCutoff = 0.10),
    error = function(e) NULL
  )
  
  if (!is.null(kegg_all) && nrow(kegg_all@result) > 0) {
    # 保存完整结果
    write.csv(as.data.frame(kegg_all),
              file.path(detailed_dir, paste0(key, "_KEGG_all.csv")),
              row.names = FALSE)
    
    # 筛选显著通路
    kegg_sig <- kegg_all@result[kegg_all@result$p.adjust < 0.05, ]
    if (nrow(kegg_sig) > 0) {
      # 气泡图
      if (nrow(kegg_sig) > 20) kegg_sig <- kegg_sig[1:20, ]
      p_bubble <- dotplot(kegg_all, showCategory = 20, title = paste0("KEGG: ", key))
      ggsave(file.path(detailed_dir, paste0(key, "_KEGG_bubble.pdf")),
             p_bubble, width = 10, height = 8)
      ggsave(file.path(detailed_dir, paste0(key, "_KEGG_bubble.png")),
             p_bubble, width = 10, height = 8, dpi = 200)
      
      # 通路网络图
      tryCatch({
        p_cnet <- cnetplot(kegg_all, showCategory = 10, circular = FALSE)
        ggsave(file.path(detailed_dir, paste0(key, "_KEGG_cnet.pdf")),
               p_cnet, width = 12, height = 10)
      }, error = function(e) NULL)
    }
  }
  
  # ===== GO-BP 富集 =====
  cat(sprintf("  [%s] 运行GO-BP富集...\n", key))
  go_all <- tryCatch(
    enrichGO(gene = entrez_all, OrgDb = org.Mm.eg.db, ont = "BP",
             pvalueCutoff = 0.05, qvalueCutoff = 0.10, readable = TRUE),
    error = function(e) NULL
  )
  
  if (!is.null(go_all) && nrow(go_all@result) > 0) {
    write.csv(as.data.frame(go_all),
              file.path(detailed_dir, paste0(key, "_GO_BP_all.csv")),
              row.names = FALSE)
    
    go_sig <- go_all@result[go_all@result$p.adjust < 0.05, ]
    if (nrow(go_sig) > 0) {
      # 气泡图
      if (nrow(go_sig) > 20) go_sig <- go_sig[1:20, ]
      p_go <- dotplot(go_all, showCategory = 20, title = paste0("GO-BP: ", key))
      ggsave(file.path(detailed_dir, paste0(key, "_GO_BP_bubble.pdf")),
             p_go, width = 10, height = 8)
      
      # 基因概念网络图
      tryCatch({
        p_gocnet <- cnetplot(go_all, showCategory = 10, circular = FALSE)
        ggsave(file.path(detailed_dir, paste0(key, "_GO_BP_cnet.pdf")),
               p_gocnet, width = 12, height = 10)
      }, error = function(e) NULL)
    }
  }
  
  # ===== 上调基因特异性富集 =====
  if (!is.null(entrez_up) && length(entrez_up) >= 5) {
    cat(sprintf("  [%s] 上调基因KEGG...\n", key))
    kegg_up <- tryCatch(
      enrichKEGG(gene = entrez_up, organism = "mmu", 
                 pvalueCutoff = 0.05, qvalueCutoff = 0.10),
      error = function(e) NULL
    )
    if (!is.null(kegg_up) && nrow(kegg_up@result) > 0) {
      write.csv(as.data.frame(kegg_up),
                file.path(detailed_dir, paste0(key, "_KEGG_up.csv")),
                row.names = FALSE)
    }
    
    go_up <- tryCatch(
      enrichGO(gene = entrez_up, OrgDb = org.Mm.eg.db, ont = "BP",
               pvalueCutoff = 0.05, qvalueCutoff = 0.10, readable = TRUE),
      error = function(e) NULL
    )
    if (!is.null(go_up) && nrow(go_up@result) > 0) {
      write.csv(as.data.frame(go_up),
                file.path(detailed_dir, paste0(key, "_GO_BP_up.csv")),
                row.names = FALSE)
    }
  }
  
  # ===== 下调基因特异性富集 =====
  if (!is.null(entrez_down) && length(entrez_down) >= 5) {
    cat(sprintf("  [%s] 下调基因KEGG...\n", key))
    kegg_down <- tryCatch(
      enrichKEGG(gene = entrez_down, organism = "mmu", 
                 pvalueCutoff = 0.05, qvalueCutoff = 0.10),
      error = function(e) NULL
    )
    if (!is.null(kegg_down) && nrow(kegg_down@result) > 0) {
      write.csv(as.data.frame(kegg_down),
                file.path(detailed_dir, paste0(key, "_KEGG_down.csv")),
                row.names = FALSE)
    }
    
    go_down <- tryCatch(
      enrichGO(gene = entrez_down, OrgDb = org.Mm.eg.db, ont = "BP",
               pvalueCutoff = 0.05, qvalueCutoff = 0.10, readable = TRUE),
      error = function(e) NULL
    )
    if (!is.null(go_down) && nrow(go_down@result) > 0) {
      write.csv(as.data.frame(go_down),
                file.path(detailed_dir, paste0(key, "_GO_BP_down.csv")),
                row.names = FALSE)
    }
  }
  
  # 返回富集结果摘要
  return(list(
    key = key,
    n_genes = length(sig_all),
    n_up = length(sig_up),
    n_down = length(sig_down),
    kegg_n = ifelse(!is.null(kegg_all) && !is.null(kegg_all@result), 
                    sum(kegg_all@result$p.adjust < 0.05, na.rm = TRUE), 0),
    go_n = ifelse(!is.null(go_all) && !is.null(go_all@result), 
                  sum(go_all@result$p.adjust < 0.05, na.rm = TRUE), 0)
  ))
}

# ==================== 4. 单任务DEG函数 ====================
run_one_deg <- function(merged_obj, comp_name, group1, group2, ct_name,
                        fm_cores = CORES_FINDMARKERS) {

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

  plan("multicore", workers = fm_cores)

  res <- tryCatch({
    FindMarkers(obj_sub,
                ident.1         = group1,
                ident.2         = group2,
                test.use        = "wilcox",
                min.pct         = MIN_PCT,
                logfc.threshold = 0,
                only.pos        = FALSE,
                verbose         = FALSE)
  }, error = function(e) {
    message(sprintf("    ✗ FindMarkers失败 [%s|%s]: %s", comp_name, ct_name, e$message))
    NULL
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

  message(sprintf("  ✓ %-48s [%s] Up=%-4d Down=%d",
                  comp_name, ct_name,
                  sum(res$sig=="Up"), sum(res$sig=="Down")))
  return(res)
}

# ==================== 5. 任务构建 ====================

build_task_list <- function(comparisons, cell_types, seurat_objs,
                             deg_dir, enrich_dir, plot_dir, organ_tag) {
  tasks <- list()
  for (comp in comparisons) {
    for (ct_name in cell_types) {

      merged <- NULL
      if (!is.null(comp$g1_samples)) {
        g1_objs   <- Filter(Negate(is.null), seurat_objs[comp$g1_samples])
        g1_merged <- merge_to_group(g1_objs, comp$g1_label, ct_name)
        g2_cells  <- extract_cells(seurat_objs[[comp$g2_obj]], ct_name)
        if (!is.null(g1_merged) && length(g2_cells) >= MIN_CELLS) {
          g2_sub <- subset(seurat_objs[[comp$g2_obj]], cells = g2_cells)
          g2_sub$group <- comp$g2_label
          merged <- merge(g1_merged, g2_sub)
          merged <- JoinLayers(merged)
        }
      } else {
        merged <- merge_two_groups(
          seurat_objs[[comp$g1_obj]], seurat_objs[[comp$g2_obj]],
          ct_name, comp$g1_label, comp$g2_label
        )
      }

      if (is.null(merged)) next

      tasks[[length(tasks) + 1]] <- list(
        merged     = merged,
        comp_name  = comp$comp,
        group1     = comp$g1_label,
        group2     = comp$g2_label,
        ct_name    = ct_name,
        key        = paste0(comp$comp, "_", ct_name),
        deg_dir    = deg_dir,
        enrich_dir = enrich_dir,
        plot_dir   = plot_dir,
        organ_tag  = organ_tag
      )
    }
  }
  return(tasks)
}

# ==================== 6. 定义比较组 ====================

liver_comparisons <- list(
  list(comp="liver_tumor_vs_blank",
       g1_samples=c("gan56","liver3","liver8"), g1_label="Tumor",
       g2_obj="bla_liver", g2_label="Blank"),
  list(comp="liver_tumor_vs_control",
       g1_samples=c("gan56","liver3","liver8"), g1_label="Tumor",
       g2_obj="con_liver", g2_label="Control"),
  list(comp="liver_blank_vs_control",
       g1_obj="bla_liver", g1_label="Blank",
       g2_obj="con_liver", g2_label="Control")
)

brain_comparisons <- list(
  list(comp="brain_tumor_vs_blank",
       g1_samples=c("nao56","brain3","brain8"), g1_label="TumorBrain",
       g2_obj="bla_brain", g2_label="Blank"),
  list(comp="brain_tumor_vs_control",
       g1_samples=c("nao56","brain3","brain8"), g1_label="TumorBrain",
       g2_obj="con_brain", g2_label="Control"),
  list(comp="brain_blank_vs_control",
       g1_obj="bla_brain", g1_label="Blank",
       g2_obj="con_brain", g2_label="Control")
)

# ==================== 7. 构建所有任务 ====================
cat("\n========== 2. 构建任务列表 ==========\n")

liver_tasks <- build_task_list(
  liver_comparisons, liver_cell_types, seurat_objs,
  file.path(output_dir,"01_liver_DEG"),
  file.path(output_dir,"03_liver_enrichment"),
  file.path(output_dir,"06_plots"),
  organ_tag = "liver"
)

brain_tasks <- build_task_list(
  brain_comparisons, brain_cell_types, seurat_objs,
  file.path(output_dir,"02_brain_DEG"),
  file.path(output_dir,"04_brain_enrichment"),
  file.path(output_dir,"06_plots"),
  organ_tag = "brain"
)

all_tasks <- c(liver_tasks, brain_tasks)
cat(sprintf("  总任务数: %d（肝脏 %d + 脑部 %d）\n\n",
            length(all_tasks), length(liver_tasks), length(brain_tasks)))

# ==================== 8. 并行执行所有任务 ====================
cat("========== 3. 并行执行DEG ==========\n")

CONCURRENT_TASKS <- min(CORES_TASK, length(all_tasks), 4)
cat(sprintf("  同时运行 %d 个任务\n", CONCURRENT_TASKS))

t_start <- proc.time()

results_raw <- mclapply(all_tasks, function(task) {
  res <- run_one_deg(task$merged, task$comp_name,
                     task$group1, task$group2, task$ct_name, 
                     CORES_FINDMARKERS)
  
  if (!is.null(res)) {
    # 保存 DEG CSV
    write.csv(res, file.path(task$deg_dir, paste0(task$key, ".csv")),
              row.names = FALSE)
    
    # 火山图
    tryCatch({
      res_v <- res
      res_v$log10p <- -log10(res_v$p_val_adj + 1e-300)
      top_genes <- bind_rows(
        res_v %>% filter(sig=="Up")   %>% arrange(p_val_adj) %>% head(10),
        res_v %>% filter(sig=="Down") %>% arrange(p_val_adj) %>% head(10)
      )
      n_up   <- sum(res_v$sig == "Up")
      n_down <- sum(res_v$sig == "Down")
      p_v <- ggplot(res_v, aes(x=avg_log2FC, y=log10p, color=sig)) +
        geom_point(alpha=0.5, size=0.8) +
        scale_color_manual(values=c("Up"="#E41A1C","Down"="#377EB8","NS"="grey75")) +
        geom_vline(xintercept=c(-FC_THRESHOLD, FC_THRESHOLD),
                   linetype="dashed", linewidth=0.4) +
        geom_hline(yintercept=-log10(PADJ_CUT),
                   linetype="dashed", linewidth=0.4) +
        geom_text_repel(data=top_genes, aes(label=gene),
                        size=2.5, max.overlaps=30, segment.size=0.2) +
        labs(title    = paste0(task$comp_name, "  [", task$ct_name, "]"),
             subtitle = sprintf("Wilcoxon | FC>%.2f | p_adj<%.2f | Up=%d Down=%d",
                                FC_THRESHOLD, PADJ_CUT, n_up, n_down),
             x="avg_log2FC", y="-log10(p_val_adj)", color=NULL) +
        theme_bw(base_size=11) +
        theme(plot.title    = element_text(hjust=0.5, size=10, face="bold"),
              plot.subtitle = element_text(hjust=0.5, size=8, color="grey40"))
      base_v <- file.path(task$plot_dir, "volcano",
                          paste0(task$comp_name, "_", task$ct_name))
      ggsave(paste0(base_v,"_volcano.pdf"), p_v, width=8, height=6, dpi=300)
      ggsave(paste0(base_v,"_volcano.png"), p_v, width=8, height=6, dpi=200)
    }, error = function(e) NULL)
  }
  
  return(list(key=task$key, res=res))
}, mc.cores = CONCURRENT_TASKS,
   mc.preschedule = FALSE,
   mc.cleanup = TRUE)

t_elapsed <- proc.time() - t_start
cat(sprintf("\n  ✓ DEG完成，耗时 %.1f 秒\n", t_elapsed["elapsed"]))

# 整理结果
all_deg <- list()
for (r in results_raw) {
  if (!is.null(r) && !is.null(r$res)) {
    all_deg[[r$key]] <- r$res
  }
}

liver_deg <- all_deg[grepl("^liver", names(all_deg))]
brain_deg  <- all_deg[grepl("^brain", names(all_deg))]
cat(sprintf("  成功完成: 肝脏 %d 组，脑部 %d 组\n",
            length(liver_deg), length(brain_deg)))

# ==================== 9. 完整富集分析 ====================
cat("\n\n========== 4. 完整富集分析 ==========\n")

enrich_results <- list()

# 收集所有需要富集的结果
all_for_enrich <- c(liver_deg, brain_deg)

if (length(all_for_enrich) > 0) {
  enrich_tasks <- lapply(names(all_for_enrich), function(key) {
    list(
      deg_res = all_for_enrich[[key]],
      key = key,
      enrich_dir = ifelse(grepl("^liver", key),
                          file.path(output_dir, "03_liver_enrichment"),
                          file.path(output_dir, "04_brain_enrichment")),
      organ_tag = ifelse(grepl("^liver", key), "liver", "brain")
    )
  })
  
  # 并行执行富集分析
  enrich_raw <- mclapply(enrich_tasks, function(task) {
    run_complete_enrichment(task$deg_res, task$key, 
                            task$enrich_dir, task$organ_tag)
  }, mc.cores = min(CORES_ENRICH, length(enrich_tasks), 4),
     mc.preschedule = FALSE)
  
  for (i in seq_along(enrich_raw)) {
    if (!is.null(enrich_raw[[i]])) {
      enrich_results[[names(all_for_enrich)[i]]] <- enrich_raw[[i]]
    }
  }
  
  cat(sprintf("\n  ✓ 富集分析完成: %d 组有显著结果\n", length(enrich_results)))
} else {
  cat("  ⚠️ 无DEG结果，跳过富集分析\n")
}

# ==================== 10. 配对样本分析 ====================
cat("\n\n========== 5. 配对样本分析 ==========\n")

paired_liver_cts <- c("All_Cells","Endothelial","Kupffer","Hepatocytes")
paired_brain_cts <- c("All_Cells","Endothelial","Microglia_Homeostatic")  # ← 修改

paired_tasks <- list()
for (pair in paired_samples) {
  lo <- seurat_objs[[pair$liver]]
  bo <- seurat_objs[[pair$brain]]
  if (is.null(lo) || is.null(bo)) next

  for (ct in paired_liver_cts) {
    m <- merge_two_groups(lo, seurat_objs[["con_liver"]],
                           ct, "TumorLiver","Control")
    if (!is.null(m)) {
      paired_tasks[[length(paired_tasks)+1]] <- list(
        merged     = m,
        comp_name  = paste0("paired_",pair$id,"_liver"),
        group1     = "TumorLiver", group2 = "Control",
        ct_name    = ct,
        key        = paste0("paired_",pair$id,"_liver_",ct),
        deg_dir    = file.path(output_dir,"05_crosstalk"),
        enrich_dir = file.path(output_dir,"05_crosstalk"),
        plot_dir   = file.path(output_dir,"06_plots")
      )
    }
  }
  for (ct in paired_brain_cts) {
    m <- merge_two_groups(bo, seurat_objs[["con_brain"]],
                           ct, "TumorBrain","Control")
    if (!is.null(m)) {
      paired_tasks[[length(paired_tasks)+1]] <- list(
        merged     = m,
        comp_name  = paste0("paired_",pair$id,"_brain"),
        group1     = "TumorBrain", group2 = "Control",
        ct_name    = ct,
        key        = paste0("paired_",pair$id,"_brain_",ct),
        deg_dir    = file.path(output_dir,"05_crosstalk"),
        enrich_dir = file.path(output_dir,"05_crosstalk"),
        plot_dir   = file.path(output_dir,"06_plots")
      )
    }
  }
}

cat(sprintf("  配对任务数: %d\n", length(paired_tasks)))

if (length(paired_tasks) > 0) {
  paired_results <- mclapply(paired_tasks, function(task) {
    res <- run_one_deg(task$merged, task$comp_name,
                       task$group1, task$group2, task$ct_name,
                       CORES_FINDMARKERS)
    if (!is.null(res)) {
      write.csv(res, file.path(task$deg_dir, paste0(task$key, ".csv")),
                row.names = FALSE)
    }
    return(list(key=task$key, res=res))
  }, mc.cores = min(CORES_TASK, length(paired_tasks), 2),
     mc.preschedule = FALSE)

  paired_deg <- list()
  for (r in paired_results) {
    if (!is.null(r) && !is.null(r$res)) paired_deg[[r$key]] <- r$res
  }
  cat(sprintf("  配对分析完成: %d 组\n", length(paired_deg)))
} else {
  paired_deg <- list()
}

# ==================== 11. 生成汇总报告 ====================
cat("\n\n========== 6. 生成汇总报告 ==========\n")

all_deg_final <- c(liver_deg, brain_deg, paired_deg)

if (length(all_deg_final) > 0) {
  deg_summary <- do.call(rbind, lapply(names(all_deg_final), function(nm) {
    res <- all_deg_final[[nm]]
    # 提取 Il6st 信息
    il6st_info <- if ("Il6st" %in% res$gene) {
      il6st <- res[res$gene == "Il6st", ]
      c(il6st$avg_log2FC, il6st$p_val_adj, il6st$sig)
    } else {
      c(NA, NA, "Not Found")
    }
    
    data.frame(
      Comparison = nm,
      CellType = ifelse(grepl("Microglia_Homeostatic", nm), "Microglia_Homeostatic",
                       ifelse(grepl("Microglia", nm), "Microglia", "Other")),
      Up         = sum(res$sig=="Up",   na.rm=TRUE),
      Down       = sum(res$sig=="Down", na.rm=TRUE),
      Total      = sum(res$sig!="NS",   na.rm=TRUE),
      Il6st_FC   = as.numeric(il6st_info[1]),
      Il6st_padj = as.numeric(il6st_info[2]),
      Il6st_sig  = il6st_info[3],
      Method     = "Wilcoxon_FindMarkers",
      FC_cut     = FC_THRESHOLD,
      Padj_cut   = PADJ_CUT,
      stringsAsFactors = FALSE
    )
  }))
  write.csv(deg_summary,
            file.path(output_dir,"07_summary","DEG_count_summary.csv"),
            row.names=FALSE)
  cat("✓ DEG汇总已保存\n")
  
  # 打印 Il6st 结果
  cat("\n【Il6st 结果汇总】\n")
  cat("──────────────────────────────────────────────────────────────────────\n")
  il6st_summary <- deg_summary[!is.na(deg_summary$Il6st_FC), 
                                c("Comparison", "CellType", "Il6st_FC", "Il6st_padj", "Il6st_sig")]
  print(il6st_summary)
  cat("──────────────────────────────────────────────────────────────────────\n")
}

# 富集结果汇总
if (length(enrich_results) > 0) {
  enrich_summary <- do.call(rbind, lapply(names(enrich_results), function(nm) {
    er <- enrich_results[[nm]]
    data.frame(
      Comparison = nm,
      Total_Genes = er$n_genes,
      Up_Genes = er$n_up,
      Down_Genes = er$n_down,
      KEGG_Sig_Pathways = er$kegg_n,
      GO_BP_Sig_Pathways = er$go_n,
      stringsAsFactors = FALSE
    )
  }))
  write.csv(enrich_summary,
            file.path(output_dir,"08_enrichment_detailed","enrichment_summary.csv"),
            row.names=FALSE)
  cat("✓ 富集汇总已保存\n")
}

# 生成分析报告
report_path <- file.path(output_dir,"07_summary","analysis_report.txt")
sink(report_path)
cat("========================================\n")
cat("肝脑肿瘤分析报告 v5（稳态小胶质细胞版）\n")
cat("========================================\n")
cat("时间:", format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"\n\n")
cat("【系统配置】\n")
cat(sprintf("  可用核心数:      %d\n", available_cores))
cat(sprintf("  FindMarkers:     %d 核\n", CORES_FINDMARKERS))
cat(sprintf("  任务级并行:      %d 核\n", CORES_TASK))
cat(sprintf("  富集分析并行:    %d 核\n", CORES_ENRICH))
cat(sprintf("  同时运行任务:    %d 个\n", CONCURRENT_TASKS))
cat(sprintf("  总耗时:          %.1f 秒\n\n", t_elapsed["elapsed"]))
cat("【细胞类型配置】\n")
cat(sprintf("  肝脏: %s\n", paste(liver_cell_types, collapse=", ")))
cat(sprintf("  脑部: %s\n", paste(brain_cell_types, collapse=", ")))
cat("\n【DEG参数】\n")
cat(sprintf("  方法:           Wilcoxon (Seurat::FindMarkers)\n"))
cat(sprintf("  avg_log2FC:     > %.2f\n", FC_THRESHOLD))
cat(sprintf("  p_val_adj:      < %.2f\n", PADJ_CUT))
cat(sprintf("  DownSample:     每组最多 %d 个细胞\n", MAX_CELLS))
cat(sprintf("  min.pct:        %.2f\n\n", MIN_PCT))
cat("【DEG数量统计】\n")
if (nrow(deg_summary) > 0) print(deg_summary) else cat("  无显著结果\n")
cat("\n【Il6st 结果】\n")
if (exists("il6st_summary") && nrow(il6st_summary) > 0) print(il6st_summary) else cat("  未检测到 Il6st\n")
sink()

cat("✓ 报告已保存:", report_path,"\n")
cat("\n========================================\n")
cat("分析完成！\n")
cat(sprintf("  肝脏: %d 组  脑部: %d 组  配对: %d 组  富集: %d 组\n",
            length(liver_deg), length(brain_deg), length(paired_deg), length(enrich_results)))
cat(sprintf("  总耗时: %.1f 秒\n", t_elapsed["elapsed"]))
cat(sprintf("  输出目录: %s\n", output_dir))
cat("========================================\n")
