# ==================== 转录因子分析（读取已有cell_species版）====================
# 项目：26年第82个项目 / 26.2.25第2个项目
# 用法：conda activate R_4.5.2 && Rscript tf_analysis.R
#
# 【前提】
#   rebuild_species_table.R 已运行完成
#   gan56/liver3/liver8 的 RDS 中已有 $cell_species 列
#
# 【分析框架】
#  肝脏
#  ├── 策略2：小鼠间质细胞 → Tumor vs Control/Blank → 宿主微环境
#  └── 策略1（描述性）：人源肿瘤细胞 → 肿瘤内TF网络+拟时序
#  脑部（纯小鼠）→ Tumor vs Control/Blank → 远端效应
# =============================================================================

cat("========================================\n")
cat("转录因子分析（cell_species已注入版）\n")
cat("========================================\n\n")

# ==================== 0. 环境 ====================
.libPaths("/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/lib/R/library")

suppressPackageStartupMessages({
  library(igraph); library(Seurat); library(dorothea); library(viper)
  library(dplyr); library(tidyr); library(ggplot2); library(ggrepel)
  library(patchwork); library(RColorBrewer)
  library(ComplexHeatmap); library(circlize); library(ggraph)
})

select    <- dplyr::select;  filter    <- dplyr::filter
mutate    <- dplyr::mutate;  arrange   <- dplyr::arrange
group_by  <- dplyr::group_by; summarise <- dplyr::summarise
bind_rows <- dplyr::bind_rows; rename   <- dplyr::rename

cat("✓ 包加载完成\n\n")

# ==================== 1. 参数 ====================
setwd("/mnt/BioDisk/wangmingyuan/SGS10_analysis")

LIVER_TUMOR_SAMPLES <- list(
  list(name="gan56",  group="Tumor", organ="Liver"),
  list(name="liver3", group="Tumor", organ="Liver"),
  list(name="liver8", group="Tumor", organ="Liver")
)
LIVER_CTRL_SAMPLES <- list(
  list(name="bla_liver", group="Blank",   organ="Liver"),
  list(name="con_liver", group="Control", organ="Liver")
)
BRAIN_SAMPLES <- list(
  list(name="nao56",    group="Tumor",   organ="Brain"),
  list(name="brain3",   group="Tumor",   organ="Brain"),
  list(name="brain8",   group="Tumor",   organ="Brain"),
  list(name="bla_brain",group="Blank",   organ="Brain"),
  list(name="con_brain",group="Control", organ="Brain")
)

CELLTYPE_COL    <- "final_celltype"
DOROTHEA_LEVELS <- c("A","B","C")
MIN_CELLS       <- 10
TOP_TF          <- 30
TOP_TF_LABEL    <- 20
PSEUDOTIME_DIR  <- "pseudotime_multicurve_optimized"
GROUP_COLORS    <- c("Tumor"="#E41A1C","Blank"="#FF7F00","Control"="#377EB8")
OUTPUT_ROOT     <- "tf_analysis_results"

for (d in c(OUTPUT_ROOT,
            file.path(OUTPUT_ROOT,
                      c("Liver_S2_mouse","Liver_S1_human","Brain","comparison"))))
  if (!dir.exists(d)) dir.create(d, recursive=TRUE)

# ==================== 2. 基础工具函数 ====================

load_one_sample <- function(s) {
  fp <- file.path(s$name, paste0(s$name,"_manual_annotated.rds"))
  if (!file.exists(fp)) { cat(sprintf("  ✗ %-12s 不存在\n",s$name)); return(NULL) }
  obj <- tryCatch(readRDS(fp), error=function(e){
    cat(sprintf("  ✗ %-12s %s\n",s$name,e$message)); NULL })
  if (is.null(obj)) return(NULL)
  obj <- JoinLayers(obj)
  obj <- RenameCells(obj, add.cell.id=s$name)
  obj$sample <- s$name; obj$group <- s$group; obj$organ <- s$organ
  cat(sprintf("  ✓ %-12s n=%5d  %s\n", s$name, ncol(obj), s$group))
  return(obj)
}

merge_samples <- function(sl) {
  objs <- Filter(Negate(is.null), lapply(sl, load_one_sample))
  if (length(objs)==0) stop("无可用样本")
  cat(sprintf("  合并%d个样本，共%d个细胞\n",length(objs),sum(sapply(objs,ncol))))
  JoinLayers(Reduce(merge, objs))
}

save_fig <- function(p, path, w=10, h=8) {
  ggsave(paste0(path,".pdf"), p, width=w, height=h, dpi=300)
  ggsave(paste0(path,".png"), p, width=w, height=h, dpi=200)
}

get_expr <- function(obj)
  as.matrix(tryCatch(
    GetAssayData(obj, assay="RNA", layer="data"),
    error=function(e) GetAssayData(obj, assay="RNA", layer="counts")))

# 基因名物种分类（备用，主要用 cell_species 列）
classify_genes <- function(genes) {
  list(
    human = genes[grepl("^[A-Z][A-Z0-9\\-]+$",genes) & !grepl("^[A-Z][a-z]",genes)],
    mouse = genes[grepl("^[A-Z][a-z]",genes)]
  )
}

# ==================== 3. 读取 cell_species（核心改动）====================

# 从已有 meta.data 直接读取 cell_species
# 注意：RenameCells 后 colnames 加了样本前缀，meta.data 里的列还在
get_species_cells <- function(obj, species="mouse") {
  if (!"cell_species" %in% colnames(obj@meta.data)) {
    cat(sprintf("  ⚠ 未找到 cell_species 列，用基因名规则推断...\n"))
    # 降级方案：用基因名比例
    gc   <- classify_genes(rownames(obj))
    expr <- GetAssayData(obj, assay="RNA", layer="counts")
    h_s  <- Matrix::colSums(expr[intersect(gc$human,rownames(expr)),,drop=FALSE])
    m_s  <- Matrix::colSums(expr[intersect(gc$mouse,rownames(expr)),,drop=FALSE])
    tot  <- h_s + m_s
    hr   <- ifelse(tot>0, h_s/tot, 0)
    obj$cell_species <- ifelse(hr>0.30,"human","mouse")
    obj$human_ratio  <- hr
    cat(sprintf("  推断: 人类=%d 小鼠=%d\n",
                sum(obj$cell_species=="human"), sum(obj$cell_species=="mouse")))
  }
  cells <- colnames(obj)[obj$cell_species == species]
  return(cells)
}

# 诊断：打印 cell_species 分布
print_species_stats <- function(obj, sample_name) {
  if ("cell_species" %in% colnames(obj@meta.data)) {
    tbl <- table(obj$cell_species)
    cat(sprintf("  [%s] cell_species: %s\n", sample_name,
                paste(paste0(names(tbl),"=",as.integer(tbl)), collapse=" | ")))
  } else {
    cat(sprintf("  [%s] ⚠ 无 cell_species 列\n", sample_name))
  }
}

# ==================== 4. Regulon 加载 ====================

cat("加载 DoRothEA regulon...\n")
data("dorothea_hs", package="dorothea")
data("dorothea_mm", package="dorothea")
reg_hs <- dorothea_hs %>% dplyr::filter(confidence %in% DOROTHEA_LEVELS)
reg_mm <- dorothea_mm %>% dplyr::filter(confidence %in% DOROTHEA_LEVELS)

to_viper_reg <- function(reg) {
  tfs <- unique(reg$tf)
  lapply(setNames(tfs,tfs), function(tf) {
    s <- reg[reg$tf==tf,]
    list(tfmode=setNames(s$mor,s$target), likelihood=rep(1,nrow(s)))
  })
}
viper_hs <- to_viper_reg(reg_hs)
viper_mm <- to_viper_reg(reg_mm)
cat(sprintf("  hs: %d TF  |  mm: %d TF\n",
            length(unique(reg_hs$tf)), length(unique(reg_mm$tf))))
cat("✓ Regulon 准备完成\n\n")

# ==================== 5. VIPER 推断 ====================

run_viper <- function(expr, vreg, tag="") {
  targets  <- unique(unlist(lapply(vreg, function(r) names(r$tfmode))))
  n_match  <- length(intersect(rownames(expr), targets))
  coverage <- n_match / length(targets) * 100
  cat(sprintf("  [%s] 基因覆盖度: %d/%d (%.1f%%)\n",
              tag, n_match, length(targets), coverage))
  if (coverage < 5) { cat("  ⚠ 覆盖度过低，跳过\n"); return(NULL) }
  res <- tryCatch(
    viper(eset=expr, regulon=vreg, nes=TRUE, method="scale",
          minsize=4, adaptive.size=FALSE, verbose=FALSE),
    error=function(e){ cat(sprintf("  ✗ VIPER: %s\n",e$message)); NULL })
  if (!is.null(res))
    cat(sprintf("  ✓ %d TF × %d 细胞\n", nrow(res), ncol(res)))
  return(res)
}

# ==================== 6. 可视化函数 ====================

plot_tf_hm <- function(tf_mat, meta, col, title="", out=NULL) {
  groups <- unique(as.character(meta[[col]]))
  groups <- groups[!is.na(groups)]
  gm <- sapply(groups, function(g) {
    cells <- intersect(rownames(meta)[as.character(meta[[col]])==g], colnames(tf_mat))
    if (length(cells)==0) return(rep(NA,nrow(tf_mat)))
    rowMeans(tf_mat[,cells,drop=FALSE])
  })
  colnames(gm) <- groups
  tv  <- apply(gm,1,function(x) var(x,na.rm=TRUE))
  tfs <- names(sort(tv,decreasing=TRUE))[1:min(TOP_TF,sum(!is.na(tv)))]
  mat <- t(gm[tfs,,drop=FALSE])

  cf  <- colorRamp2(c(-2,0,2),c("#4575B4","white","#D73027"))
  cols <- if (col=="group") GROUP_COLORS[groups]
          else setNames(
            colorRampPalette(brewer.pal(min(12,length(groups)),"Paired"))(length(groups)),
            groups)
  ha  <- HeatmapAnnotation(
    df=setNames(data.frame(groups,row.names=groups),col),
    col=setNames(list(cols),col),
    annotation_name_side="left")
  ht  <- Heatmap(t(mat), name="NES", top_annotation=ha, col=cf,
                 cluster_rows=TRUE, cluster_columns=FALSE,
                 show_column_names=TRUE,
                 column_names_gp=gpar(fontsize=9),
                 row_names_gp=gpar(fontsize=8),
                 column_title=title,
                 column_title_gp=gpar(fontsize=11,fontface="bold"))
  if (!is.null(out)) {
    h <- max(8,length(tfs)*0.32+3); w <- max(6,length(groups)*1.5+3)
    pdf(paste0(out,".pdf"),width=w,height=h); draw(ht); dev.off()
    png(paste0(out,".png"),width=round(w*200),height=round(h*200),res=200)
    draw(ht); dev.off()
  }
  return(ht)
}

diff_tf <- function(tf_mat, meta, g1="Tumor", g2="Control", ct_col=NULL) {
  groups <- if (!is.null(ct_col)) unique(as.character(meta[[ct_col]])) else "ALL"
  res <- lapply(groups, function(ct) {
    if (ct=="ALL") {
      c1 <- intersect(rownames(meta)[meta$group==g1],colnames(tf_mat))
      c2 <- intersect(rownames(meta)[meta$group==g2],colnames(tf_mat))
    } else {
      c1 <- intersect(rownames(meta)[as.character(meta[[ct_col]])==ct&meta$group==g1],colnames(tf_mat))
      c2 <- intersect(rownames(meta)[as.character(meta[[ct_col]])==ct&meta$group==g2],colnames(tf_mat))
    }
    if (length(c1)<MIN_CELLS||length(c2)<MIN_CELLS) return(NULL)
    delta <- rowMeans(tf_mat[,c1,drop=FALSE]) - rowMeans(tf_mat[,c2,drop=FALSE])
    pv    <- sapply(rownames(tf_mat), function(tf)
      tryCatch(wilcox.test(tf_mat[tf,c1],tf_mat[tf,c2])$p.value,error=function(e)1))
    data.frame(TF=rownames(tf_mat),delta_NES=delta,
               pval=pv,padj=p.adjust(pv,"BH"),
               celltype=ct,comp=paste0(g1,"_vs_",g2),
               stringsAsFactors=FALSE) %>%
      dplyr::arrange(padj,desc(abs(delta_NES)))
  })
  d <- dplyr::bind_rows(Filter(Negate(is.null),res))
  if (nrow(d)==0) return(NULL)
  d
}

plot_vol <- function(d, g1="Tumor", g2="Control", title="", out=NULL) {
  s <- d %>% dplyr::group_by(TF) %>%
    dplyr::summarise(md=mean(delta_NES),mp=min(padj),.groups="drop") %>%
    dplyr::mutate(lp=-log10(mp+1e-300),
                  sig=dplyr::case_when(
                    mp<0.05&md> 0.5~paste0("Up in ",g1),
                    mp<0.05&md<-0.5~paste0("Up in ",g2),
                    TRUE~"NS"))
  lb <- dplyr::bind_rows(
    s%>%dplyr::filter(sig==paste0("Up in ",g1))%>%dplyr::arrange(desc(lp))%>%head(TOP_TF_LABEL),
    s%>%dplyr::filter(sig==paste0("Up in ",g2))%>%dplyr::arrange(desc(lp))%>%head(TOP_TF_LABEL))
  vc <- setNames(c("#E41A1C","#377EB8","grey70"),
                 c(paste0("Up in ",g1),paste0("Up in ",g2),"NS"))
  p <- ggplot(s,aes(md,lp,color=sig))+
    geom_point(alpha=0.6,size=1.5)+
    scale_color_manual(values=vc,name=NULL)+
    geom_text_repel(data=lb,aes(label=TF),size=2.8,max.overlaps=30,segment.size=0.2)+
    geom_vline(xintercept=c(-0.5,0.5),linetype="dashed",linewidth=0.4,color="grey40")+
    geom_hline(yintercept=-log10(0.05),linetype="dashed",linewidth=0.4,color="grey40")+
    labs(title=title,x=paste0("ΔNES（",g1,"−",g2,"）"),y="-log10(padj)")+
    theme_bw(base_size=11)+
    theme(plot.title=element_text(hjust=0.5,face="bold"),legend.position="top")
  if (!is.null(out)) save_fig(p,out,w=9,h=7)
  return(p)
}

plot_grn <- function(tf_mat, reg, n_tf=20, title="", out=NULL) {
  tm   <- rowMeans(tf_mat)
  ctfs <- names(sort(tm,decreasing=TRUE))[1:min(n_tf,length(tm))]
  e    <- reg %>% dplyr::filter(tf%in%ctfs) %>%
    dplyr::group_by(tf) %>% dplyr::slice_max(abs(mor),n=8) %>% dplyr::ungroup() %>%
    dplyr::mutate(type=ifelse(mor>0,"Activation","Repression"))
  if (nrow(e)<5) return(NULL)
  g <- igraph::graph_from_data_frame(e[,c("tf","target","type","mor")],directed=TRUE)
  igraph::V(g)$is_tf <- igraph::V(g)$name %in% ctfs
  p <- ggraph(g,layout="fr")+
    geom_edge_link(aes(color=type,width=abs(mor)),alpha=0.45,
                   arrow=arrow(length=unit(1.5,"mm"),type="closed"),
                   end_cap=circle(3,"mm"))+
    geom_node_point(aes(color=is_tf,size=is_tf))+
    geom_node_label(aes(label=ifelse(is_tf,name,"")),size=2.8,repel=TRUE,max.overlaps=30)+
    scale_color_manual(values=c("TRUE"="#E41A1C","FALSE"="#6BAED6"),
                       labels=c("靶基因","TF"),name="节点")+
    scale_size_manual(values=c("TRUE"=5,"FALSE"=2),labels=c("靶基因","TF"),name="节点")+
    scale_edge_color_manual(values=c("Activation"="#E41A1C","Repression"="#377EB8"),name="调控")+
    scale_edge_width_continuous(range=c(0.3,1.5),name="|MOR|")+
    labs(title=title)+theme_graph(base_size=11)+
    theme(plot.title=element_text(hjust=0.5,face="bold"))
  if (!is.null(out)) save_fig(p,out,w=14,h=11)
  return(p)
}

run_tf_pt <- function(tf_mat, obj, organ, out_dir) {
  cat("  读取拟时序...\n")
  pt_paths <- c(
    file.path(PSEUDOTIME_DIR,organ,"Tumor","pseudotime_values.csv"),
    file.path(PSEUDOTIME_DIR,organ,"pseudotime_values.csv"),
    file.path("pseudotime_results",organ,"pseudotime_values.csv")
  )
  pf <- Filter(file.exists,pt_paths)[1]
  if (is.na(pf)) { cat("  ⚠ 未找到拟时序文件\n"); return(NULL) }
  pt  <- read.csv(pf)
  pt  <- pt[!is.infinite(pt$pseudotime)&!is.na(pt$pseudotime),]
  com <- intersect(pt$cell,colnames(tf_mat))
  if (length(com)<50) {
    bc_tf <- sub("^.*_","",colnames(tf_mat))
    bc_pt <- sub("^.*_","",pt$cell)
    m <- match(bc_pt,bc_tf); v <- which(!is.na(m))
    com <- colnames(tf_mat)[m[v]]; pt <- pt[v,]; pt$cell <- com
  }
  if (length(com)<50) { cat("  ⚠ 匹配不足，跳过\n"); return(NULL) }
  cat(sprintf("  ✓ 匹配 %d 个细胞\n",length(com)))
  pv <- setNames(pt$pseudotime,com)
  ts <- tf_mat[,com,drop=FALSE]
  cr <- apply(ts,1,function(x) cor(x,pv,method="spearman",use="complete.obs"))
  cd <- data.frame(TF=names(cr),cor=cr) %>%
    dplyr::arrange(desc(abs(cor))) %>% dplyr::filter(!is.na(cor))
  write.csv(cd,file.path(out_dir,"TF_pseudotime_correlation.csv"),row.names=FALSE)

  top20 <- cd %>%
    dplyr::slice(c(1:min(10,n()),(n()-min(10,n())+1):n())) %>%
    dplyr::mutate(dir=ifelse(cor>0,"随伪时间升高","随伪时间降低"))
  p_bar <- ggplot(top20,aes(x=reorder(TF,cor),y=cor,fill=dir))+
    geom_bar(stat="identity")+
    scale_fill_manual(values=c("随伪时间升高"="#E41A1C","随伪时间降低"="#377EB8"))+
    coord_flip()+
    labs(title=paste0(organ," · TF-伪时间相关性（Top20）"),x=NULL,y="Spearman r",fill=NULL)+
    theme_bw(base_size=11)+theme(plot.title=element_text(hjust=0.5,face="bold"),
                                  legend.position="top")
  save_fig(p_bar,file.path(out_dir,"TF_pseudotime_bar"),w=8,h=7)

  tp <- cd%>%dplyr::filter(cor>0.1)%>%head(12)%>%dplyr::pull(TF)
  tn <- cd%>%dplyr::filter(cor< -0.1)%>%head(12)%>%dplyr::pull(TF)
  dy <- c(tp,tn)
  if (length(dy)>=4) {
    ord <- order(pv); co <- com[ord]; to <- ts[dy,co,drop=FALSE]
    if (ncol(to)>2000) { idx <- round(seq(1,ncol(to),length.out=2000)); to <- to[,idx]; co <- co[idx] }
    ha <- HeatmapAnnotation(Pseudotime=pv[co],
                             col=list(Pseudotime=colorRamp2(range(pv[co]),c("#FFFFCC","#800026"))))
    cf <- colorRamp2(c(-2,0,2),c("#4575B4","white","#D73027"))
    for (fmt in c("pdf","png")) {
      bm <- file.path(out_dir,"TF_pseudotime_heatmap")
      if(fmt=="pdf") pdf(paste0(bm,".pdf"),width=16,height=max(6,length(dy)*0.5+3))
      else png(paste0(bm,".png"),width=3200,height=max(1200,length(dy)*100+600),res=200)
      draw(Heatmap(to,name="NES",top_annotation=ha,col=cf,
                   cluster_rows=TRUE,cluster_columns=FALSE,show_column_names=FALSE,
                   row_names_gp=gpar(fontsize=9),
                   column_title=paste0(organ," · 动态TF × 伪时间"),
                   column_title_gp=gpar(fontsize=11,fontface="bold")))
      dev.off()
    }
    cat("  ✓ 动态TF热图\n")
  }
  return(cd)
}

# ==================== 7. 策略2：小鼠间质细胞 ====================

run_s2_liver <- function() {
  cat("\n",strrep("=",60),"\n策略2：肝脏小鼠间质细胞（宿主微环境）\n",strrep("=",60),"\n\n")
  od <- file.path(OUTPUT_ROOT,"Liver_S2_mouse")

  # 加载全部肝脏样本
  cat("加载肝脏样本...\n")
  obj_all <- tryCatch(
    merge_samples(c(LIVER_TUMOR_SAMPLES,LIVER_CTRL_SAMPLES)),
    error=function(e){cat("✗",e$message,"\n");NULL})
  if (is.null(obj_all)) return(NULL)

  # 诊断 cell_species 分布
  cat("cell_species 分布:\n")
  for (s in c("gan56","liver3","liver8")) {
    cells_s <- colnames(obj_all)[obj_all$sample==s]
    if (length(cells_s)>0) {
      sub_obj <- obj_all[,cells_s]
      print_species_stats(sub_obj, s)
    }
  }
  cat("  （bla_liver/con_liver 为纯小鼠，不含 human 细胞）\n")

  # 提取小鼠细胞
  # con_liver/bla_liver 全是小鼠；gan56/liver3/liver8 用 cell_species=="mouse"
  ctrl_cells  <- colnames(obj_all)[obj_all$sample %in% c("bla_liver","con_liver")]
  tumor_mouse <- get_species_cells(obj_all[,colnames(obj_all)[
    obj_all$sample %in% c("gan56","liver3","liver8")]], "mouse")
  mouse_cells <- c(ctrl_cells, tumor_mouse)

  cat(sprintf("\n小鼠细胞汇总:\n"))
  cat(sprintf("  肿瘤样本小鼠细胞: %d\n", length(tumor_mouse)))
  cat(sprintf("  对照/空白细胞:    %d\n", length(ctrl_cells)))
  cat(sprintf("  合计:             %d\n", length(mouse_cells)))

  if (length(mouse_cells)<50) { cat("⚠ 细胞数不足\n"); return(NULL) }
  obj_m <- subset(obj_all, cells=mouse_cells)

  cat("\n分组 × 细胞类型分布:\n")
  print(table(obj_m@meta.data[[CELLTYPE_COL]], obj_m$group))

  # 只用小鼠基因做 VIPER
  gc     <- classify_genes(rownames(obj_m))
  expr   <- get_expr(obj_m)
  expr_m <- expr[intersect(gc$mouse,rownames(expr)),,drop=FALSE]
  cat(sprintf("小鼠基因数: %d\n",nrow(expr_m)))

  cat("\nVIPER推断（mouse regulon）...\n")
  tf <- run_viper(expr_m, viper_mm, "S2_mouse")
  if (is.null(tf)) return(NULL)

  meta <- obj_m@meta.data

  # 热图
  cat("热图...\n")
  plot_tf_hm(tf,meta,"group",
    "肝脏小鼠间质细胞 · 各组TF活性",
    file.path(od,"01_TF_heatmap_group"))
  if (length(unique(meta[[CELLTYPE_COL]]))>1)
    plot_tf_hm(tf,meta,CELLTYPE_COL,
      "肝脏小鼠间质细胞 · 各细胞类型TF活性",
      file.path(od,"02_TF_heatmap_celltype"))
  cat("  ✓\n")

  # 差异TF
  cat("差异TF...\n")
  dc <- diff_tf(tf,meta,"Tumor","Control",CELLTYPE_COL)
  if (!is.null(dc)) {
    write.csv(dc,file.path(od,"03_diff_TF_TvC.csv"),row.names=FALSE)
    plot_vol(dc,"Tumor","Control","肝脏间质细胞 · Tumor vs Control",
             file.path(od,"04_volcano_TvC"))
    cat("  ✓ TvC\n")
  }
  db <- diff_tf(tf,meta,"Tumor","Blank",CELLTYPE_COL)
  if (!is.null(db)) {
    write.csv(db,file.path(od,"03b_diff_TF_TvB.csv"),row.names=FALSE)
    plot_vol(db,"Tumor","Blank","肝脏间质细胞 · Tumor vs Blank",
             file.path(od,"04b_volcano_TvB"))
    cat("  ✓ TvB\n")
  }

  # GRN
  cat("GRN...\n")
  plot_grn(tf,reg_mm,20,"肝脏小鼠间质细胞 · TF调控网络",
           file.path(od,"05_GRN"))
  cells_t <- intersect(colnames(obj_m)[obj_m$group=="Tumor"],colnames(tf))
  cells_c <- intersect(colnames(obj_m)[obj_m$group=="Control"],colnames(tf))
  if (length(cells_t)>MIN_CELLS && length(cells_c)>MIN_CELLS) {
    dlt <- rowMeans(tf[,cells_t,drop=FALSE]) - rowMeans(tf[,cells_c,drop=FALSE])
    tt  <- names(sort(dlt,decreasing=TRUE))[1:min(15,length(dlt))]
    plot_grn(tf[tt,,drop=FALSE],reg_mm%>%dplyr::filter(tf%in%tt),
             length(tt),"肝脏间质细胞 · 肿瘤特异TF网络",
             file.path(od,"05b_GRN_tumor"))
  }
  cat("  ✓\n")

  run_tf_pt(tf,obj_m,"Liver",od)
  write.csv(as.data.frame(tf),file.path(od,"TF_activity_matrix.csv"))
  cat(sprintf("✓ 策略2完成 → %s\n",od))
  return(list(tf=tf,meta=meta,diff_ctrl=dc,diff_blank=db))
}

# ==================== 8. 策略1：人源肿瘤细胞（描述性）====================

run_s1_liver <- function() {
  cat("\n",strrep("=",60),"\n策略1：肝脏人源肿瘤细胞（描述性）\n",strrep("=",60),"\n\n")
  od <- file.path(OUTPUT_ROOT,"Liver_S1_human")

  cat("加载肿瘤样本...\n")
  obj <- tryCatch(merge_samples(LIVER_TUMOR_SAMPLES),
                  error=function(e){cat("✗",e$message,"\n");NULL})
  if (is.null(obj)) return(NULL)

  # 直接读取 cell_species
  human_cells <- get_species_cells(obj, "human")
  cat(sprintf("人类细胞: %d / %d\n", length(human_cells), ncol(obj)))
  if (length(human_cells)<50) { cat("⚠ 人类细胞不足50，跳过\n"); return(NULL) }

  obj_h  <- subset(obj, cells=human_cells)
  cat("人类细胞分布（各样本 × 细胞类型）:\n")
  print(table(obj_h@meta.data[[CELLTYPE_COL]], obj_h$sample))

  # 用人类基因做 VIPER（human regulon）
  gc     <- classify_genes(rownames(obj_h))
  expr   <- get_expr(obj_h)
  expr_h <- expr[intersect(gc$human,rownames(expr)),,drop=FALSE]
  cat(sprintf("人类基因数: %d\n",nrow(expr_h)))

  cat("\nVIPER推断（human regulon）...\n")
  tf <- run_viper(expr_h, viper_hs, "S1_human")
  if (is.null(tf)) return(NULL)

  meta    <- obj_h@meta.data
  ct_vals <- unique(as.character(meta[[CELLTYPE_COL]]))

  # 热图
  cat("热图...\n")
  if (length(ct_vals)>1)
    plot_tf_hm(tf,meta,CELLTYPE_COL,
      "人源肿瘤细胞 · 各亚群TF活性",
      file.path(od,"01_TF_heatmap_celltype"))
  plot_tf_hm(tf,meta,"sample",
    "人源肿瘤细胞 · 各样本TF活性",
    file.path(od,"02_TF_heatmap_sample"))
  cat("  ✓\n")

  # 亚群均值 + 气泡图
  if (length(ct_vals)>1) {
    ct_mean <- sapply(ct_vals, function(ct) {
      cells <- intersect(rownames(meta)[as.character(meta[[CELLTYPE_COL]])==ct],colnames(tf))
      if (length(cells)==0) return(rep(NA,nrow(tf)))
      rowMeans(tf[,cells,drop=FALSE])
    })
    colnames(ct_mean) <- ct_vals
    write.csv(as.data.frame(ct_mean),
              file.path(od,"03_TF_mean_by_celltype.csv"))

    top_ct <- as.data.frame(ct_mean) %>%
      dplyr::mutate(TF=rownames(ct_mean)) %>%
      tidyr::pivot_longer(-TF,names_to="celltype",values_to="NES") %>%
      dplyr::group_by(celltype) %>% dplyr::slice_max(NES,n=5) %>% dplyr::ungroup()

    if (nrow(top_ct)>=3) {
      p_b <- ggplot(top_ct,aes(x=celltype,y=TF,size=NES,color=NES))+
        geom_point(alpha=0.85)+
        scale_color_gradient2(low="#377EB8",mid="white",high="#E41A1C",
                               midpoint=0,name="NES")+
        scale_size_continuous(range=c(2,8))+
        labs(title="人源肿瘤细胞 · 各亚群Top TF",x=NULL,y=NULL)+
        theme_bw(base_size=11)+
        theme(plot.title=element_text(hjust=0.5,face="bold"),
              axis.text.x=element_text(angle=35,hjust=1))
      save_fig(p_b,file.path(od,"04_TF_bubble"),
               w=max(7,length(ct_vals)*1.2+3),h=8)
      cat("  ✓ 亚群TF气泡图\n")
    }
  }

  # GRN
  cat("GRN（肿瘤细胞调控程序）...\n")
  plot_grn(tf,reg_hs,20,"人源肿瘤细胞 · TF调控网络",
           file.path(od,"05_GRN"))
  cat("  ✓\n")

  run_tf_pt(tf,obj_h,"Liver",od)
  write.csv(as.data.frame(tf),file.path(od,"TF_activity_matrix.csv"))
  cat(sprintf("✓ 策略1完成 → %s\n",od))
  return(list(tf=tf,meta=meta))
}

# ==================== 9. 脑部（纯小鼠）====================

run_brain_tf <- function() {
  cat("\n",strrep("=",60),"\n脑部转录因子分析（纯小鼠）\n",strrep("=",60),"\n\n")
  od <- file.path(OUTPUT_ROOT,"Brain")

  cat("加载脑部样本...\n")
  obj <- tryCatch(merge_samples(BRAIN_SAMPLES),
                  error=function(e){cat("✗",e$message,"\n");NULL})
  if (is.null(obj)) return(NULL)

  gc   <- classify_genes(rownames(obj))
  expr <- get_expr(obj)
  expr_m <- expr[intersect(gc$mouse,rownames(expr)),,drop=FALSE]
  cat(sprintf("小鼠基因数: %d\n",nrow(expr_m)))

  cat("\nVIPER推断（mouse regulon）...\n")
  tf <- run_viper(expr_m, viper_mm, "brain")
  if (is.null(tf)) return(NULL)

  meta <- obj@meta.data

  # 热图
  cat("热图...\n")
  plot_tf_hm(tf,meta,"group","脑部 · 各组TF活性",
             file.path(od,"01_TF_heatmap_group"))
  if (length(unique(as.character(meta[[CELLTYPE_COL]])))>1)
    plot_tf_hm(tf,meta,CELLTYPE_COL,"脑部 · 各细胞类型TF活性",
               file.path(od,"02_TF_heatmap_celltype"))
  cat("  ✓\n")

  # 差异TF
  cat("差异TF...\n")
  dc <- diff_tf(tf,meta,"Tumor","Control",CELLTYPE_COL)
  if (!is.null(dc)) {
    write.csv(dc,file.path(od,"03_diff_TF_TvC.csv"),row.names=FALSE)
    plot_vol(dc,"Tumor","Control","脑部 · Tumor vs Control",
             file.path(od,"04_volcano_TvC"))
    cat("  ✓ TvC\n")
  }
  db <- diff_tf(tf,meta,"Tumor","Blank",CELLTYPE_COL)
  if (!is.null(db)) {
    write.csv(db,file.path(od,"03b_diff_TF_TvB.csv"),row.names=FALSE)
    plot_vol(db,"Tumor","Blank","脑部 · Tumor vs Blank",
             file.path(od,"04b_volcano_TvB"))
    cat("  ✓ TvB\n")
  }

  # GRN
  cat("GRN...\n")
  plot_grn(tf,reg_mm,20,"脑部 · TF调控网络",file.path(od,"05_GRN"))
  cells_t <- intersect(colnames(obj)[obj$group=="Tumor"],colnames(tf))
  cells_c <- intersect(colnames(obj)[obj$group=="Control"],colnames(tf))
  if (length(cells_t)>MIN_CELLS && length(cells_c)>MIN_CELLS) {
    dlt <- rowMeans(tf[,cells_t,drop=FALSE]) - rowMeans(tf[,cells_c,drop=FALSE])
    tt  <- names(sort(dlt,decreasing=TRUE))[1:min(15,length(dlt))]
    plot_grn(tf[tt,,drop=FALSE],reg_mm%>%dplyr::filter(tf%in%tt),
             length(tt),"脑部 · 肿瘤特异TF网络",
             file.path(od,"05b_GRN_tumor"))
  }
  cat("  ✓\n")

  run_tf_pt(tf,obj,"Brain",od)
  write.csv(as.data.frame(tf),file.path(od,"TF_activity_matrix.csv"))
  cat(sprintf("✓ 脑部完成 → %s\n",od))
  return(list(tf=tf,meta=meta,diff_ctrl=dc,diff_blank=db))
}

# ==================== 10. 跨器官对比 ====================

run_comparison <- function(rs2, rb) {
  cat("\n",strrep("=",60),"\n跨器官TF对比\n",strrep("=",60),"\n\n")
  od <- file.path(OUTPUT_ROOT,"comparison")

  get_sum <- function(d) {
    if (is.null(d)) return(NULL)
    d %>% dplyr::group_by(TF) %>%
      dplyr::summarise(md=mean(delta_NES),mp=min(padj),.groups="drop")
  }
  sl <- get_sum(rs2$diff_ctrl); sb <- get_sum(rb$diff_ctrl)
  if (is.null(sl)||is.null(sb)) { cat("⚠ 数据不足\n"); return(NULL) }

  ctf <- intersect(sl$TF,sb$TF)
  df  <- data.frame(TF=ctf,
                    Liver=sl$md[match(ctf,sl$TF)],
                    Brain=sb$md[match(ctf,sb$TF)]) %>%
    dplyr::mutate(consistent=sign(Liver)==sign(Brain),
                  label=abs(Liver)>0.3|abs(Brain)>0.3)
  write.csv(df,file.path(od,"cross_organ_TF_delta.csv"),row.names=FALSE)

  p <- ggplot(df,aes(Liver,Brain,color=consistent))+
    geom_point(alpha=0.65,size=1.8)+
    geom_text_repel(data=dplyr::filter(df,label),aes(label=TF),
                    size=2.6,max.overlaps=30,segment.size=0.2)+
    scale_color_manual(values=c("TRUE"="#E41A1C","FALSE"="grey60"),
                       labels=c("方向不一致","肝脑一致"),name=NULL)+
    geom_hline(yintercept=0,linetype="dashed",linewidth=0.4,color="grey40")+
    geom_vline(xintercept=0,linetype="dashed",linewidth=0.4,color="grey40")+
    labs(title="肝脏间质细胞 vs 脑部 · TF差异活性（Tumor vs Control）",
         x="肝脏间质 ΔNES",y="脑部 ΔNES")+
    theme_bw(base_size=12)+
    theme(plot.title=element_text(hjust=0.5,face="bold"),legend.position="top")
  save_fig(p,file.path(od,"cross_organ_scatter"),w=9,h=8)

  cu <- df %>% dplyr::filter(consistent,Liver>0.2,Brain>0.2) %>%
    dplyr::arrange(desc(Liver+Brain))
  if (nrow(cu)>0) {
    write.csv(cu,file.path(od,"consistent_up_TF.csv"),row.names=FALSE)
    cat(sprintf("  肝脑一致上调TF: %d 个\n",nrow(cu)))
    print(head(cu,10))
  }
  cat("  ✓ 跨器官对比完成\n")
}

# ==================== 11. 执行 ====================

t0 <- proc.time()
rs2    <- run_s2_liver()
rs1    <- run_s1_liver()
rb     <- run_brain_tf()
if (!is.null(rs2) && !is.null(rb)) run_comparison(rs2, rb)

cat(sprintf("\n%s\n全部完成！耗时: %.1f 秒\n%s\n",
            strrep("=",40),(proc.time()-t0)["elapsed"],strrep("=",40)))
cat(sprintf("\n输出: %s/\n",OUTPUT_ROOT))
cat("├── Liver_S2_mouse/  策略2：小鼠间质细胞（Tumor vs Control/Blank）\n")
cat("├── Liver_S1_human/  策略1：人源肿瘤细胞（描述性）\n")
cat("├── Brain/           脑部（纯小鼠，最干净）\n")
cat("└── comparison/      肝脑跨器官TF对比\n")
