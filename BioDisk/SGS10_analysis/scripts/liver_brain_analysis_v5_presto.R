# ==================== 肝脑肿瘤分析脚本 v5（冲突修复最终版） ====================
# 文件名: liver_brain_analysis_v5_presto.R
# 用法:   conda activate R_4.5.2 && Rscript liver_brain_analysis_v5_presto.R
#
# 【并行策略】
#   每个任务 FindMarkers 内部用 4 核（future multicore）
#   同时并发 = 128 ÷ 4 = 32 个任务
#   总占用   = 32 × 4 = 128 核
#
# 【修复】
#   AnnotationDbi::select 与 dplyr::select 冲突
#   → 全局绑定 dplyr 动词 + 子进程内所有调用显式用 dplyr:: 前缀
#
# 【修改】
#   脑部小胶质细胞: Microglia → Microglia_Homeostatic
# =============================================================================

# ==================== 0. 包加载 ====================
cat("========================================\n")
cat("肝脑肿瘤分析 v5（稳态小胶质细胞版）\n")
cat("========================================\n")
cat("加载依赖包...\n")

.libPaths("/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/lib/R/library")

suppressPackageStartupMessages({
  library(igraph)           # 必须在 Seurat 前
  library(Seurat)
  library(presto)           # Wilcoxon 加速
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
})

# ---- 关键：全局绑定 dplyr 动词，防止 AnnotationDbi 覆盖 ----
select    <- dplyr::select
rename    <- dplyr::rename
filter    <- dplyr::filter
mutate    <- dplyr::mutate
arrange   <- dplyr::arrange
summarise <- dplyr::summarise
group_by  <- dplyr::group_by
left_join <- dplyr::left_join
inner_join<- dplyr::inner_join
bind_rows <- dplyr::bind_rows
distinct  <- dplyr::distinct

cat("✓ 所有包加载完成，dplyr 命名空间已锁定\n\n")

# ==================== 1. 全局参数 ====================

available_cores  <- detectCores()
CORES_PER_TASK   <- 4
CORES_TOTAL      <- 128
CONCURRENT_TASKS <- floor(CORES_TOTAL / CORES_PER_TASK)   # = 32

cat(sprintf("服务器: %d核  策略: 每任务%d核 × %d并发 = %d核\n\n",
            available_cores, CORES_PER_TASK,
            CONCURRENT_TASKS, CORES_PER_TASK * CONCURRENT_TASKS))

options(future.globals.maxSize = 128 * 1024^3)

FC_THRESHOLD <- 0.20
PADJ_CUT     <- 0.10
MIN_PCT      <- 0.10
MAX_CELLS    <- 500
MIN_CELLS    <- 5

liver_cell_types <- c("All_Cells", "Hepatocytes", "Endothelial", "Kupffer")

# ==================== 关键修改：Microglia → Microglia_Homeostatic ====================
brain_cell_types <- c("All_Cells", "Endothelial", "Microglia_Homeostatic")

paired_samples <- list(
  list(liver="gan56",  brain="nao56",  id="mouse56"),
  list(liver="liver3", brain="brain3", id="mouse3"),
  list(liver="liver8", brain="brain8", id="mouse8")
)

mechanism_genes <- list(
  IL6_JAK_STAT3  = c("Il6","Il6ra","Il6st","Jak1","Jak2","Stat3","Socs3","Gp130"),
  TGFb_SMAD      = c("Tgfb1","Tgfb2","Tgfbr1","Tgfbr2","Smad2","Smad3","Smad7"),
  NFkB           = c("Nfkb1","Nfkb2","Rela","Relb","Tnf","Il1b"),
  PI3K_AKT       = c("Pik3ca","Pik3cb","Akt1","Akt2","Mtor","Pten"),
  MAPK           = c("Mapk1","Mapk3","Mapk8","Mapk14","Map2k1","Map2k2"),
  Microglia_activ= c("Cd68","Cx3cr1","Tmem119","P2ry12","Hexb",
                     "Tspan18","Basp1","Il6st","Cd74","Smad3",
                     "Cd180","Ccr5","Nav2","Creb5","Ifitm10")
)
all_mech_genes <- unique(unlist(mechanism_genes))

# ==================== 2. 输出目录 ====================
setwd("/mnt/BioDisk/wangmingyuan/SGS10_analysis")

output_dir <- "liver_brain_v5_homeostatic_results"
dirs <- c(output_dir, file.path(output_dir, c(
  "01_liver_DEG","02_brain_DEG",
  "03_liver_enrichment","04_brain_enrichment",
  "05_crosstalk",
  "06_plots/volcano","06_plots/heatmap",
  "06_plots/venn","06_plots/bubble",
  "06_plots/IL6_mechanism",
  "07_summary"
)))
for (d in dirs) if (!dir.exists(d)) dir.create(d, recursive=TRUE)

ERROR_LOG <- file.path(output_dir, "07_summary", "error_log.txt")
writeLines(sprintf(
  "========================================\n错误日志 - %s\n========================================",
  format(Sys.time(),"%Y-%m-%d %H:%M:%S")), ERROR_LOG)

write_error <- function(context, msg) {
  line <- sprintf("[%s] %s => %s", format(Sys.time(),"%H:%M:%S"), context, msg)
  cat("  ⚠", line, "\n")
  cat(line, "\n", file=ERROR_LOG, append=TRUE)
}

cat(sprintf("输出目录: %s\n错误日志: %s\n\n", output_dir, ERROR_LOG))

# ==================== 3. 加载样本 ====================
cat("========== 加载样本 ==========\n")

samples_to_load <- c("gan56","liver3","liver8",
                     "nao56","brain3","brain8",
                     "con_brain","bla_brain","bla_liver","con_liver")

seurat_objs <- list()
for (s in samples_to_load) {
  fp <- file.path(s, paste0(s,"_manual_annotated.rds"))
  if (!file.exists(fp)) { write_error(s,"文件不存在"); next }
  obj <- tryCatch(readRDS(fp), error=function(e){ write_error(s,e$message); NULL })
  if (is.null(obj)) next
  obj <- JoinLayers(obj)
  obj <- RenameCells(obj, add.cell.id=s)
  seurat_objs[[s]] <- obj
  cat(sprintf("  ✓ %-15s 细胞数: %d\n", s, ncol(obj)))
}
cat(sprintf("  成功加载: %d / %d\n\n", length(seurat_objs), length(samples_to_load)))

# ==================== 4. 辅助函数 ====================

extract_cells <- function(obj, cell_type="All_Cells") {
  if (cell_type != "All_Cells" && "final_celltype" %in% colnames(obj@meta.data))
    return(WhichCells(obj, expression=final_celltype == cell_type))
  colnames(obj)
}

merge_to_group <- function(obj_list, group_label, cell_type) {
  sub_list <- lapply(obj_list, function(obj) {
    cells <- extract_cells(obj, cell_type)
    if (length(cells) < MIN_CELLS) return(NULL)
    s <- subset(obj, cells=cells); s$group <- group_label; s
  })
  sub_list <- Filter(Negate(is.null), sub_list)
  if (length(sub_list) == 0) return(NULL)
  JoinLayers(Reduce(merge, sub_list))
}

merge_two_groups <- function(obj1, obj2, cell_type, label1, label2) {
  c1 <- extract_cells(obj1, cell_type)
  c2 <- extract_cells(obj2, cell_type)
  if (length(c1) < MIN_CELLS || length(c2) < MIN_CELLS) return(NULL)
  s1 <- subset(obj1, cells=c1); s1$group <- label1
  s2 <- subset(obj2, cells=c2); s2$group <- label2
  JoinLayers(merge(s1, s2))
}

collect_results <- function(raw_list, stage_name) {
  ok <- list(); n_err <- 0
  for (i in seq_along(raw_list)) {
    r <- raw_list[[i]]
    if (is.null(r)) { n_err <- n_err+1; next }
    if (inherits(r,"try-error") || inherits(r,"error") || is.character(r)) {
      msg <- if (is.character(r)) r else conditionMessage(r)
      cat(sprintf("[%s] 任务%d失败: %s\n", stage_name, i, msg),
          file=ERROR_LOG, append=TRUE)
      n_err <- n_err+1; next
    }
    if (!is.list(r) || is.null(r$key) || is.null(r$res)) { n_err <- n_err+1; next }
    ok[[r$key]] <- r$res
  }
  if (n_err > 0)
    cat(sprintf("  ⚠ %s: %d个任务失败（详见error_log.txt）\n", stage_name, n_err))
  cat(sprintf("  %s: 成功 %d 个\n", stage_name, length(ok)))
  ok
}

# ==================== 5. DEG 核心函数 ====================
run_one_deg <- function(merged_obj, comp_name, group1, group2, ct_name,
                        fm_cores=CORES_PER_TASK) {
  n1 <- sum(merged_obj$group == group1)
  n2 <- sum(merged_obj$group == group2)
  if (n1 < MIN_CELLS || n2 < MIN_CELLS) return(NULL)

  set.seed(42)
  cg1 <- colnames(merged_obj)[merged_obj$group == group1]
  cg2 <- colnames(merged_obj)[merged_obj$group == group2]
  if (length(cg1) > MAX_CELLS) cg1 <- sample(cg1, MAX_CELLS)
  if (length(cg2) > MAX_CELLS) cg2 <- sample(cg2, MAX_CELLS)
  obj_sub <- subset(merged_obj, cells=c(cg1,cg2))
  Idents(obj_sub) <- obj_sub$group

  plan("multicore", workers=fm_cores)
  res <- tryCatch(
    FindMarkers(obj_sub,
                ident.1=group1, ident.2=group2,
                test.use="wilcox", min.pct=MIN_PCT,
                logfc.threshold=0, only.pos=FALSE, verbose=FALSE),
    error=function(e) stop(sprintf("[%s|%s] %s", comp_name, ct_name, e$message))
  )
  plan("sequential")

  if (is.null(res) || nrow(res)==0) return(NULL)

  res$gene <- rownames(res)
  res <- res %>%
    dplyr::select(gene, avg_log2FC, pct.1, pct.2, p_val, p_val_adj) %>%
    dplyr::rename(pct1=pct.1, pct2=pct.2) %>%
    dplyr::arrange(p_val_adj, dplyr::desc(abs(avg_log2FC)))

  res$sig <- ifelse(
    abs(res$avg_log2FC) > FC_THRESHOLD & res$p_val_adj < PADJ_CUT,
    ifelse(res$avg_log2FC > 0, "Up", "Down"), "NS"
  )

  message(sprintf("  ✓ %-50s [%s] Up=%-4d Down=%d",
                  comp_name, ct_name,
                  sum(res$sig=="Up"), sum(res$sig=="Down")))
  return(res)
}

# ==================== 6. 富集分析 ====================
run_enrichment <- function(res, key, enrich_dir) {
  tryCatch({
    sig_genes <- res %>% dplyr::filter(sig != "NS") %>% dplyr::pull(gene)
    if (length(sig_genes) < 10) return(invisible(NULL))

    gene_df <- tryCatch(
      bitr(sig_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Mm.eg.db),
      error=function(e) NULL
    )
    if (is.null(gene_df) || nrow(gene_df) < 5) return(invisible(NULL))
    entrez <- unique(gene_df$ENTREZID)

    kegg <- tryCatch(
      enrichKEGG(gene=entrez, organism="mmu", pvalueCutoff=0.05),
      error=function(e) NULL
    )
    if (!is.null(kegg) && nrow(kegg@result) > 0) {
      write.csv(as.data.frame(kegg),
                file.path(enrich_dir, paste0(key,"_KEGG.csv")), row.names=FALSE)
      n_sig <- sum(kegg@result$p.adjust < 0.05)
      if (n_sig >= 5) {
        p_dot <- dotplot(kegg, showCategory=20, title=paste0("KEGG: ",key))
        bk <- file.path(enrich_dir, paste0(key,"_KEGG_dotplot"))
        ggsave(paste0(bk,".pdf"), p_dot, width=10, height=8, dpi=300)
        ggsave(paste0(bk,".png"), p_dot, width=10, height=8, dpi=200)
      }
    }

    go <- tryCatch(
      enrichGO(gene=entrez, OrgDb=org.Mm.eg.db, ont="BP",
               pvalueCutoff=0.05, readable=TRUE),
      error=function(e) NULL
    )
    if (!is.null(go) && nrow(go@result) > 0)
      write.csv(as.data.frame(go),
                file.path(enrich_dir, paste0(key,"_GO_BP.csv")), row.names=FALSE)

  }, error=function(e) write_error(paste("富集",key), e$message))
}

# ==================== 7. 可视化函数 ====================

save_volcano <- function(res, comp_name, ct_name, plot_dir) {
  tryCatch({
    rv <- res
    rv$log10p <- -log10(rv$p_val_adj + 1e-300)
    tg <- dplyr::bind_rows(
      rv %>% dplyr::filter(sig=="Up")   %>% dplyr::arrange(p_val_adj) %>% head(10),
      rv %>% dplyr::filter(sig=="Down") %>% dplyr::arrange(p_val_adj) %>% head(10)
    )
    
    # 标记 IL-6 通路基因
    il6_genes <- c("Il6st","Jak1","Jak2","Stat3","Socs3","Il6ra")
    rv$is_il6 <- rv$gene %in% il6_genes
    
    pv <- ggplot(rv, aes(x=avg_log2FC, y=log10p, color=sig)) +
      geom_point(alpha=0.5, size=0.8) +
      scale_color_manual(values=c("Up"="#E41A1C","Down"="#377EB8","NS"="grey75")) +
      geom_vline(xintercept=c(-FC_THRESHOLD,FC_THRESHOLD),
                 linetype="dashed", linewidth=0.4) +
      geom_hline(yintercept=-log10(PADJ_CUT), linetype="dashed", linewidth=0.4) +
      geom_text_repel(data=tg, aes(label=gene),
                      size=2.5, max.overlaps=30, segment.size=0.2) +
      geom_point(data=rv %>% dplyr::filter(is_il6 & sig != "NS"),
                 aes(x=avg_log2FC, y=log10p), color="#FF7F00", size=2, shape=17) +
      labs(title    = paste0(comp_name,"  [",ct_name,"]"),
           subtitle = sprintf("Wilcoxon(presto)|FC>%.2f|p_adj<%.2f|Up=%d Down=%d\n橙色三角=IL-6通路基因",
                              FC_THRESHOLD, PADJ_CUT,
                              sum(rv$sig=="Up"), sum(rv$sig=="Down")),
           x="avg_log2FC", y="-log10(p_val_adj)", color=NULL) +
      theme_bw(base_size=11) +
      theme(plot.title    = element_text(hjust=0.5, size=10, face="bold"),
            plot.subtitle = element_text(hjust=0.5, size=8,  color="grey40"))
    bv <- file.path(plot_dir,"volcano",paste0(comp_name,"_",ct_name))
    ggsave(paste0(bv,"_volcano.pdf"), pv, width=9, height=7, dpi=300)
    ggsave(paste0(bv,"_volcano.png"), pv, width=9, height=7, dpi=200)
  }, error=function(e) write_error(paste("火山图",comp_name,ct_name), e$message))
}

save_heatmap <- function(res, merged_obj, comp_name, ct_name, plot_dir) {
  tryCatch({
    sh <- res %>%
      dplyr::filter(sig!="NS") %>%
      dplyr::arrange(dplyr::desc(abs(avg_log2FC))) %>%
      head(40) %>% dplyr::pull(gene)
    if (length(sh) < 3) return()
    set.seed(42)
    cu <- colnames(merged_obj)
    if (length(cu) > 300) cu <- sample(cu, 300)
    os <- subset(merged_obj, cells=cu)
    ex <- GetAssayData(os, layer="data")
    gu <- intersect(sh, rownames(ex))
    if (length(gu) < 3) return()
    ez <- t(scale(t(as.matrix(ex[gu,])))); ez[is.na(ez)] <- 0
    ad <- data.frame(Group=os$group, row.names=colnames(ez))
    cf <- colorRamp2(c(-2,0,2), c("#4575B4","white","#D73027"))
    bh <- file.path(plot_dir,"heatmap",paste0(comp_name,"_",ct_name,"_heatmap"))
    for (fmt in c("pdf","png")) {
      if (fmt=="pdf") pdf(paste0(bh,".pdf"), width=12,
                          height=max(6, length(gu)*0.3+2))
      else            png(paste0(bh,".png"), width=2400,
                          height=max(1200, length(gu)*60+400), res=200)
      draw(Heatmap(ez, name="Z-score",
                   top_annotation=HeatmapAnnotation(df=ad),
                   col=cf, show_column_names=FALSE,
                   row_names_gp=gpar(fontsize=8),
                   column_title=paste0(comp_name," (",ct_name,")"),
                   cluster_rows=TRUE, cluster_columns=TRUE))
      dev.off()
    }
  }, error=function(e) write_error(paste("热图",comp_name,ct_name), e$message))
}

# ==================== 8. 单任务完整执行 ====================
run_one_task <- function(task) {
  tryCatch({
    # 子进程重新加载包
    suppressPackageStartupMessages({
      library(igraph); library(Seurat); library(presto)
      library(future); library(dplyr); library(ggplot2); library(ggrepel)
      library(ComplexHeatmap); library(circlize)
      library(clusterProfiler); library(org.Mm.eg.db); library(enrichplot)
    })
    # 子进程内同样锁定 dplyr 命名空间
    select    <- dplyr::select
    rename    <- dplyr::rename
    filter    <- dplyr::filter
    mutate    <- dplyr::mutate
    arrange   <- dplyr::arrange
    summarise <- dplyr::summarise
    group_by  <- dplyr::group_by
    left_join <- dplyr::left_join
    inner_join<- dplyr::inner_join
    bind_rows <- dplyr::bind_rows
    distinct  <- dplyr::distinct

    res <- run_one_deg(task$merged, task$comp_name,
                       task$group1, task$group2, task$ct_name,
                       fm_cores=task$fm_cores)
    if (is.null(res)) return(NULL)

    write.csv(res, file.path(task$deg_dir, paste0(task$key,".csv")),
              row.names=FALSE)
    run_enrichment(res, task$key, task$enrich_dir)
    save_volcano(res, task$comp_name, task$ct_name, task$plot_dir)
    save_heatmap(res, task$merged, task$comp_name, task$ct_name, task$plot_dir)

    list(key=task$key, res=res)

  }, error=function(e) {
    cat(sprintf("[%s] 任务[%s|%s]失败: %s\n",
                format(Sys.time(),"%H:%M:%S"),
                task$comp_name, task$ct_name, e$message),
        file=ERROR_LOG, append=TRUE)
    NULL
  })
}

# ==================== 9. 任务构建 ====================
build_task_list <- function(comparisons, cell_types, seurat_objs,
                             deg_dir, enrich_dir, plot_dir) {
  tasks <- list()
  for (comp in comparisons) {
    for (ct_name in cell_types) {
      merged <- NULL
      if (!is.null(comp$g1_samples)) {
        g1_objs   <- Filter(Negate(is.null), seurat_objs[comp$g1_samples])
        g1_merged <- merge_to_group(g1_objs, comp$g1_label, ct_name)
        g2_cells  <- extract_cells(seurat_objs[[comp$g2_obj]], ct_name)
        if (!is.null(g1_merged) && length(g2_cells) >= MIN_CELLS) {
          g2_sub <- subset(seurat_objs[[comp$g2_obj]], cells=g2_cells)
          g2_sub$group <- comp$g2_label
          merged <- JoinLayers(merge(g1_merged, g2_sub))
        }
      } else {
        merged <- merge_two_groups(
          seurat_objs[[comp$g1_obj]], seurat_objs[[comp$g2_obj]],
          ct_name, comp$g1_label, comp$g2_label
        )
      }
      if (is.null(merged)) next
      tasks[[length(tasks)+1]] <- list(
        merged    = merged,
        comp_name = comp$comp,
        group1    = comp$g1_label,
        group2    = comp$g2_label,
        ct_name   = ct_name,
        key       = paste0(comp$comp,"_",ct_name),
        fm_cores  = CORES_PER_TASK,
        deg_dir   = deg_dir,
        enrich_dir= enrich_dir,
        plot_dir  = plot_dir
      )
    }
  }
  tasks
}

# ==================== 10. 定义比较组 ====================
liver_comparisons <- list(
  list(comp="liver_tumor_vs_blank",
       g1_samples=c("gan56","liver3","liver8"), g1_label="Tumor",
       g2_obj="bla_liver",  g2_label="Blank"),
  list(comp="liver_tumor_vs_control",
       g1_samples=c("gan56","liver3","liver8"), g1_label="Tumor",
       g2_obj="con_liver",  g2_label="Control"),
  list(comp="liver_blank_vs_control",
       g1_obj="bla_liver",  g1_label="Blank",
       g2_obj="con_liver",  g2_label="Control")
)
brain_comparisons <- list(
  list(comp="brain_tumor_vs_blank",
       g1_samples=c("nao56","brain3","brain8"), g1_label="TumorBrain",
       g2_obj="bla_brain",  g2_label="Blank"),
  list(comp="brain_tumor_vs_control",
       g1_samples=c("nao56","brain3","brain8"), g1_label="TumorBrain",
       g2_obj="con_brain",  g2_label="Control"),
  list(comp="brain_blank_vs_control",
       g1_obj="bla_brain",  g1_label="Blank",
       g2_obj="con_brain",  g2_label="Control")
)

# ==================== 11. 构建并执行主分析 ====================
cat("========== 构建任务列表 ==========\n")

liver_tasks <- build_task_list(
  liver_comparisons, liver_cell_types, seurat_objs,
  file.path(output_dir,"01_liver_DEG"),
  file.path(output_dir,"03_liver_enrichment"),
  file.path(output_dir,"06_plots")
)
brain_tasks <- build_task_list(
  brain_comparisons, brain_cell_types, seurat_objs,
  file.path(output_dir,"02_brain_DEG"),
  file.path(output_dir,"04_brain_enrichment"),
  file.path(output_dir,"06_plots")
)
all_tasks <- c(liver_tasks, brain_tasks)
n_main    <- min(CONCURRENT_TASKS, length(all_tasks))

cat(sprintf("  主分析任务: %d（每任务%d核，%d个并发，共%d核）\n\n",
            length(all_tasks), CORES_PER_TASK, n_main,
            n_main * CORES_PER_TASK))

cat("========== 并行执行主分析 ==========\n")
plan("sequential")
t_start <- proc.time()

results_raw <- mclapply(
  all_tasks, run_one_task,
  mc.cores=n_main, mc.preschedule=FALSE, mc.cleanup=TRUE
)

t_main <- proc.time() - t_start
cat(sprintf("  主分析耗时: %.1f 秒\n", t_main["elapsed"]))

all_deg   <- collect_results(results_raw, "主分析")
liver_deg <- all_deg[grepl("^liver", names(all_deg))]
brain_deg  <- all_deg[grepl("^brain", names(all_deg))]
cat(sprintf("  肝脏: %d 组  脑部: %d 组\n",
            length(liver_deg), length(brain_deg)))

# ==================== 12. 配对样本分析 ====================
cat("\n========== 配对样本分析 ==========\n")

paired_liver_cts <- c("All_Cells","Endothelial","Kupffer","Hepatocytes")

# ==================== 关键修改：配对分析也使用 Microglia_Homeostatic ====================
paired_brain_cts <- c("All_Cells","Endothelial","Microglia_Homeostatic")

paired_tasks <- list()
for (pair in paired_samples) {
  lo <- seurat_objs[[pair$liver]]; bo <- seurat_objs[[pair$brain]]
  if (is.null(lo)||is.null(bo)) next
  for (ct in paired_liver_cts) {
    m <- merge_two_groups(lo, seurat_objs[["con_liver"]], ct, "TumorLiver","Control")
    if (!is.null(m)) paired_tasks[[length(paired_tasks)+1]] <- list(
      merged=m, comp_name=paste0("paired_",pair$id,"_liver"),
      group1="TumorLiver", group2="Control", ct_name=ct,
      key=paste0("paired_",pair$id,"_liver_",ct),
      fm_cores=CORES_PER_TASK,
      deg_dir    = file.path(output_dir,"05_crosstalk"),
      enrich_dir = file.path(output_dir,"05_crosstalk"),
      plot_dir   = file.path(output_dir,"06_plots")
    )
  }
  for (ct in paired_brain_cts) {
    m <- merge_two_groups(bo, seurat_objs[["con_brain"]], ct, "TumorBrain","Control")
    if (!is.null(m)) paired_tasks[[length(paired_tasks)+1]] <- list(
      merged=m, comp_name=paste0("paired_",pair$id,"_brain"),
      group1="TumorBrain", group2="Control", ct_name=ct,
      key=paste0("paired_",pair$id,"_brain_",ct),
      fm_cores=CORES_PER_TASK,
      deg_dir    = file.path(output_dir,"05_crosstalk"),
      enrich_dir = file.path(output_dir,"05_crosstalk"),
      plot_dir   = file.path(output_dir,"06_plots")
    )
  }
}

paired_deg <- list()
if (length(paired_tasks) > 0) {
  np <- min(CONCURRENT_TASKS, length(paired_tasks))
  cat(sprintf("  配对任务: %d（%d个并发）\n", length(paired_tasks), np))
  paired_raw <- mclapply(paired_tasks, run_one_task,
                         mc.cores=np, mc.preschedule=FALSE)
  paired_deg <- collect_results(paired_raw, "配对分析")
}

# ==================== 13. Il6st 专项输出 ====================
cat("\n========== Il6st 专项结果 ==========\n")

# 提取所有 brain_tumor_vs_control 结果中的 Il6st
brain_tumor_ctrl <- brain_deg[grepl("brain_tumor_vs_control", names(brain_deg))]

if (length(brain_tumor_ctrl) > 0) {
  cat("\n【脑部 tumor_vs_control 中 Il6st 表达】\n")
  cat("──────────────────────────────────────────────────────────────────────\n")
  for (nm in names(brain_tumor_ctrl)) {
    res <- brain_tumor_ctrl[[nm]]
    if ("Il6st" %in% res$gene) {
      il6st <- res[res$gene == "Il6st", ]
      cat(sprintf("  %s:\n", nm))
      cat(sprintf("    avg_log2FC = %.4f\n", il6st$avg_log2FC))
      cat(sprintf("    p_val_adj  = %.2e\n", il6st$p_val_adj))
      cat(sprintf("    pct1(Tumor)= %.1f%%\n", il6st$pct1 * 100))
      cat(sprintf("    pct2(Control)= %.1f%%\n", il6st$pct2 * 100))
      cat(sprintf("    显著性     = %s\n\n", il6st$sig))
    } else {
      cat(sprintf("  %s: Il6st 不在结果中\n\n", nm))
    }
  }
  cat("──────────────────────────────────────────────────────────────────────\n")
}

# ==================== 14. 机制验证热图 ====================
cat("\n========== 机制验证热图 ==========\n")

extract_fc_matrix <- function(deg_list, gene_list, pattern="tumor_vs_control") {
  mat <- data.frame(gene=gene_list, stringsAsFactors=FALSE)
  for (nm in names(deg_list)) {
    if (!grepl(pattern, nm)) next
    res <- deg_list[[nm]]
    tmp <- res %>%
      dplyr::filter(gene %in% gene_list) %>%
      dplyr::select(gene, avg_log2FC) %>%
      dplyr::rename(!!nm := avg_log2FC)
    mat <- dplyr::left_join(mat, tmp, by="gene")
  }
  rownames(mat) <- mat$gene
  m <- as.matrix(mat[,-1,drop=FALSE]); m[is.na(m)] <- 0; m
}

fc_liver <- extract_fc_matrix(liver_deg, all_mech_genes)
fc_brain  <- extract_fc_matrix(brain_deg,  all_mech_genes)

if (ncol(fc_liver)>0) {
  df_l <- as.data.frame(fc_liver); df_l$gene <- rownames(df_l)
  write.csv(df_l, file.path(output_dir,"06_plots/IL6_mechanism",
                             "liver_mechanism_FC.csv"), row.names=FALSE)
}
if (ncol(fc_brain)>0) {
  df_b <- as.data.frame(fc_brain); df_b$gene <- rownames(df_b)
  write.csv(df_b, file.path(output_dir,"06_plots/IL6_mechanism",
                             "brain_mechanism_FC.csv"), row.names=FALSE)
}

tryCatch({
  if (ncol(fc_liver)>0 && ncol(fc_brain)>0) {
    gene_show <- union(
      rownames(fc_liver)[apply(abs(fc_liver),1,max)>0.1],
      rownames(fc_brain) [apply(abs(fc_brain), 1,max)>0.1]
    )
    gene_show <- gene_show[gene_show %in% all_mech_genes]

    if (length(gene_show) >= 3) {
      fc_comb  <- cbind(fc_liver[gene_show,,drop=FALSE],
                        fc_brain [gene_show,,drop=FALSE])
      col_anno <- data.frame(
        Organ=c(rep("Liver",ncol(fc_liver)),rep("Brain",ncol(fc_brain))),
        row.names=colnames(fc_comb)
      )
      pw_map <- data.frame(
        gene=all_mech_genes,
        Pathway=rep(names(mechanism_genes),sapply(mechanism_genes,length)),
        stringsAsFactors=FALSE
      )
      row_anno <- pw_map %>%
        dplyr::filter(gene %in% gene_show) %>%
        dplyr::distinct(gene,.keep_all=TRUE)
      rownames(row_anno) <- row_anno$gene
      row_anno <- row_anno[gene_show,,drop=FALSE]

      pw_col <- c(IL6_JAK_STAT3="#E41A1C", TGFb_SMAD="#FF7F00",
                  NFkB="#984EA3",           PI3K_AKT="#4DAF4A",
                  MAPK="#377EB8",           Microglia_activ="#A65628")
      cf_fc  <- colorRamp2(c(-1.5,0,1.5), c("#377EB8","white","#E41A1C"))

      for (fmt in c("pdf","png")) {
        bm <- file.path(output_dir,"06_plots/IL6_mechanism","mechanism_heatmap")
        if (fmt=="pdf") pdf(paste0(bm,".pdf"), width=16,
                            height=max(7,length(gene_show)*0.45+3))
        else            png(paste0(bm,".png"), width=3200,
                            height=max(1400,length(gene_show)*90+600), res=200)
        draw(Heatmap(
          fc_comb, name="avg_log2FC",
          top_annotation = HeatmapAnnotation(
            df=col_anno,
            col=list(Organ=c("Liver"="#E41A1C","Brain"="#377EB8")),
            annotation_name_side="left"),
          left_annotation = rowAnnotation(
            df=data.frame(Pathway=row_anno$Pathway,
                          row.names=rownames(row_anno)),
            col=list(Pathway=pw_col), width=unit(3,"mm")),
          col=cf_fc,
          cluster_rows=FALSE, cluster_columns=FALSE,
          row_names_side="left", row_names_gp=gpar(fontsize=9),
          column_names_gp=gpar(fontsize=7), column_names_rot=45,
          column_title="机制基因FC | 肝脏（左）vs 脑（右）| tumor vs control",
          column_title_gp=gpar(fontsize=11,fontface="bold"),
          cell_fun=function(j,i,x,y,width,height,fill) {
            v <- fc_comb[i,j]
            if (abs(v)>FC_THRESHOLD)
              grid.text(sprintf("%.2f",v),x,y,
                        gp=gpar(fontsize=6.5,col="black"))
          }
        ))
        dev.off()
      }
      cat("  ✓ 机制验证热图已保存\n")
    }
  }
}, error=function(e) write_error("机制热图", e$message))

# ==================== 15. 肝脑 KEGG 共同通路 ====================
cat("\n========== 肝脑通路联系 ==========\n")

read_kegg_dir <- function(dir) {
  files <- list.files(dir, pattern="_KEGG\\.csv$", full.names=TRUE)
  res <- list()
  for (f in files) {
    key <- gsub("_KEGG\\.csv$","",basename(f))
    df  <- tryCatch(read.csv(f), error=function(e) NULL)
    if (!is.null(df) && nrow(df)>0) res[[key]] <- df
  }
  res
}

tryCatch({
  lk <- read_kegg_dir(file.path(output_dir,"03_liver_enrichment"))
  bk <- read_kegg_dir(file.path(output_dir,"04_brain_enrichment"))
  cat(sprintf("  肝脏KEGG: %d 组  脑部KEGG: %d 组\n", length(lk), length(bk)))

  if (length(lk)>0 && length(bk)>0) {
    pool <- function(kl,tag) dplyr::bind_rows(lapply(names(kl), function(nm) {
      df <- kl[[nm]]; df$source<-nm; df$organ<-tag; df
    }))
    lp <- pool(lk,"Liver"); bp <- pool(bk,"Brain")
    lt <- lp %>% dplyr::filter(grepl("tumor_vs_control",source))
    bt <- bp %>% dplyr::filter(grepl("tumor_vs_control",source))
    li <- unique(lt$ID); bi <- unique(bt$ID); si <- intersect(li,bi)
    cat(sprintf("  肝脏通路: %d  脑部通路: %d  共同: %d\n",
                length(li),length(bi),length(si)))

    # Venn 图
    if (length(li)>0 || length(bi)>0) {
      vd <- venn.diagram(list(Liver=li,Brain=bi), filename=NULL,
                         fill=c("#E41A1C","#377EB8"), alpha=0.5,
                         cex=1.5, cat.cex=1.5, cat.pos=c(-20,20),
                         main="肿瘤 vs 对照  肝脑共同KEGG通路", main.cex=1.2)
      vb <- file.path(output_dir,"06_plots/venn","liver_brain_pathway_venn")
      pdf(paste0(vb,".pdf"),width=7,height=6); grid.draw(vd); dev.off()
      png(paste0(vb,".png"),width=1400,height=1200,res=200); grid.draw(vd); dev.off()
      cat("  ✓ Venn图已保存\n")
    }

    if (length(si)>0) {
      ls_s <- lt %>%
        dplyr::filter(ID%in%si) %>%
        dplyr::group_by(ID,Description) %>%
        dplyr::summarise(Liver_minP=min(p.adjust),.groups="drop")
      bs_s <- bt %>%
        dplyr::filter(ID%in%si) %>%
        dplyr::group_by(ID,Description) %>%
        dplyr::summarise(Brain_minP=min(p.adjust),.groups="drop")
      shared_df <- dplyr::inner_join(ls_s,bs_s,by=c("ID","Description")) %>%
        dplyr::mutate(combined_score=-log10(Liver_minP)-log10(Brain_minP)) %>%
        dplyr::arrange(dplyr::desc(combined_score))

      write.csv(shared_df,
                file.path(output_dir,"05_crosstalk","shared_KEGG_pathways.csv"),
                row.names=FALSE)
      cat("\n  【共同通路 Top20】\n")
      print(head(shared_df %>%
                   dplyr::select(ID,Description,Liver_minP,Brain_minP,combined_score),20))

      if (nrow(shared_df)>=3) {
        pd <- shared_df %>% head(25) %>%
          dplyr::mutate(
            Liver_logP=-log10(Liver_minP), Brain_logP=-log10(Brain_minP),
            Desc_short=ifelse(nchar(Description)>45,
                              paste0(substr(Description,1,42),"..."),Description)
          )

        p_sc <- ggplot(pd, aes(x=Liver_logP,y=Brain_logP,label=Desc_short)) +
          geom_point(aes(size=combined_score,color=combined_score),alpha=0.8) +
          scale_color_gradientn(
            colors=c("#4575B4","#91BFDB","#FFFFBF","#FC8D59","#D73027")) +
          geom_text_repel(size=2.8,max.overlaps=20,segment.size=0.2) +
          labs(title="肝脑共同KEGG通路",
               subtitle="X轴：肝脏显著性  Y轴：脑部显著性",
               x="-log10 P (Liver)",y="-log10 P (Brain)") +
          theme_bw(base_size=12) +
          theme(plot.title=element_text(hjust=0.5,face="bold"),
                plot.subtitle=element_text(hjust=0.5,color="grey40"))
        ggsave(file.path(output_dir,"06_plots/bubble","shared_scatter.pdf"),
               p_sc,width=10,height=8,dpi=300)
        ggsave(file.path(output_dir,"06_plots/bubble","shared_scatter.png"),
               p_sc,width=10,height=8,dpi=200)

        pl <- pd %>%
          dplyr::select(Desc_short,Liver_logP,Brain_logP) %>%
          tidyr::pivot_longer(c(Liver_logP,Brain_logP),
                              names_to="Organ",values_to="logP") %>%
          dplyr::mutate(Organ=dplyr::recode(
            Organ,"Liver_logP"="Liver","Brain_logP"="Brain"))
        p_bar <- ggplot(pl,aes(x=reorder(Desc_short,logP),y=logP,fill=Organ)) +
          geom_bar(stat="identity",position="dodge") +
          scale_fill_manual(values=c("Liver"="#E41A1C","Brain"="#377EB8")) +
          coord_flip() +
          labs(title="肝脑共同富集通路显著性比较",
               x="",y="-log10 p.adjust",fill="器官") +
          theme_bw(base_size=11) +
          theme(plot.title=element_text(hjust=0.5,face="bold"),
                axis.text.y=element_text(size=8))
        h <- max(6,nrow(pd)*0.4+2)
        ggsave(file.path(output_dir,"06_plots/bubble","shared_barplot.pdf"),
               p_bar,width=12,height=h,dpi=300)
        ggsave(file.path(output_dir,"06_plots/bubble","shared_barplot.png"),
               p_bar,width=12,height=h,dpi=200)
        cat("  ✓ 共同通路图已保存\n")
      }
    }
  }
}, error=function(e) write_error("通路联系分析", e$message))

# ==================== 16. 汇总报告 ====================
cat("\n========== 生成汇总报告 ==========\n")

all_deg_final <- c(liver_deg, brain_deg, paired_deg)
t_total       <- proc.time() - t_start

deg_summary <- if (length(all_deg_final)>0) {
  do.call(rbind, lapply(names(all_deg_final), function(nm) {
    res <- all_deg_final[[nm]]
    
    # 提取 Il6st 信息
    if ("Il6st" %in% res$gene) {
      il6st <- res[res$gene == "Il6st", ]
      il6st_fc <- il6st$avg_log2FC
      il6st_padj <- il6st$p_val_adj
      il6st_sig <- il6st$sig
      il6st_pct1 <- il6st$pct1 * 100
      il6st_pct2 <- il6st$pct2 * 100
    } else {
      il6st_fc <- NA
      il6st_padj <- NA
      il6st_sig <- "Not Found"
      il6st_pct1 <- NA
      il6st_pct2 <- NA
    }
    
    data.frame(
      Comparison = nm,
      CellType = ifelse(grepl("Microglia_Homeostatic", nm), "Microglia_Homeostatic",
                       ifelse(grepl("Microglia", nm), "Microglia", 
                             ifelse(grepl("Endothelial", nm), "Endothelial", "Other"))),
      Up         = sum(res$sig=="Up",  na.rm=TRUE),
      Down       = sum(res$sig=="Down",na.rm=TRUE),
      Total      = sum(res$sig!="NS",  na.rm=TRUE),
      Il6st_FC   = round(il6st_fc, 4),
      Il6st_padj = sprintf("%.2e", il6st_padj),
      Il6st_sig  = il6st_sig,
      Il6st_pct1 = sprintf("%.1f%%", il6st_pct1),
      Il6st_pct2 = sprintf("%.1f%%", il6st_pct2),
      Method     = "Wilcoxon_presto",
      FC_cut     = FC_THRESHOLD,
      Padj_cut   = PADJ_CUT,
      stringsAsFactors=FALSE
    )
  }))
} else {
  write_error("汇总","没有任何成功完成的DEG分析")
  data.frame()
}

if (nrow(deg_summary)>0) {
  write.csv(deg_summary,
            file.path(output_dir,"07_summary","DEG_count_summary.csv"),
            row.names=FALSE)
  cat("✓ DEG汇总已保存\n")
  
  # 打印 Il6st 结果
  cat("\n【Il6st 结果汇总】\n")
  cat("────────────────────────────────────────────────────────────────────────────\n")
  il6st_summary <- deg_summary[!is.na(deg_summary$Il6st_FC) & deg_summary$Il6st_sig != "Not Found", 
                                c("Comparison", "CellType", "Il6st_FC", "Il6st_padj", "Il6st_sig", "Il6st_pct1", "Il6st_pct2")]
  if (nrow(il6st_summary) > 0) {
    print(il6st_summary)
  } else {
    cat("  未检测到 Il6st 显著差异\n")
  }
  cat("────────────────────────────────────────────────────────────────────────────\n")
}

report_path <- file.path(output_dir,"07_summary","analysis_report.txt")
sink(report_path)
cat("========================================\n")
cat("肝脑肿瘤分析报告 v5（稳态小胶质细胞版）\n")
cat("========================================\n")
cat("时间:", format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"\n\n")
cat("【运行配置】\n")
cat(sprintf("  服务器核心:   %d\n",  available_cores))
cat(sprintf("  每任务核数:   %d 核\n", CORES_PER_TASK))
cat(sprintf("  并发任务数:   %d 个\n", n_main))
cat(sprintf("  总占用核数:   %d 核\n", n_main*CORES_PER_TASK))
cat(sprintf("  总耗时:       %.1f 秒\n\n", t_total["elapsed"]))
cat("【DEG参数】\n")
cat(sprintf("  方法:         Wilcoxon (presto加速)\n"))
cat(sprintf("  avg_log2FC:   > %.2f\n", FC_THRESHOLD))
cat(sprintf("  p_val_adj:    < %.2f\n", PADJ_CUT))
cat(sprintf("  DownSample:   每组最多 %d 个细胞\n", MAX_CELLS))
cat(sprintf("  min.pct:      %.2f\n\n", MIN_PCT))
cat("【分析细胞类型】\n")
cat("  肝脏:", paste(liver_cell_types,collapse=", "),"\n")
cat("  脑部:", paste(brain_cell_types,collapse=", "),"\n\n")
cat("【DEG数量统计】\n")
if (nrow(deg_summary)>0) print(deg_summary) else cat("  无显著结果\n")
sink()

cat("\n========================================\n")
cat("全部分析完成！\n")
cat(sprintf("  肝脏: %d 组  脑部: %d 组  配对: %d 组\n",
            length(liver_deg), length(brain_deg), length(paired_deg)))
cat(sprintf("  总耗时:   %.1f 秒\n",  t_total["elapsed"]))
cat(sprintf("  输出目录: %s\n",       output_dir))
cat(sprintf("  错误日志: %s\n",       ERROR_LOG))
cat("========================================\n")
