#!/usr/bin/env Rscript
source("/mnt/BioDisk/wangmingyuan/SGS10_analysis/scripts/fig_save.R")  # 三格式输出器
# ============================================================================
# paper_umap_grouped_v3_replot.R — Fig 1 主图
# v4 (客户版): 1) 左下角箭头坐标轴 UMAP1(->)/UMAP2(^), 去整条轴线
#              2) 图内细胞类型名字号 +2 (5.5->7.5)
#              3) 右侧图例照例图: 组标题+纯编号圆圈(多成员), 单成员写名字; 分2列
#              4) 两器官整图尺寸/主面板尺寸一致 -> 可拼接对齐
# ============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(ggplot2); library(ggrepel); library(patchwork)
})

OUT_DIR <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis/umap_paper_grade"
stopifnot(dir.exists(OUT_DIR))
RESTYLE_DIR <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis/figures/restyle_v1"
dir.create(RESTYLE_DIR, showWarnings = FALSE, recursive = TRUE)
cat("R:", R.version.string, " ggplot2:", as.character(packageVersion("ggplot2")), "\n\n")

# ============================================================================
# Group map / Color
# ============================================================================
LIVER_GROUP_MAP <- list(
  Hep = c("Hep_Periportal","Hep_Perivenous","Hep_Proliferating"),
  LSEC= c("LSEC","LSEC_Pericentral","LSEC_Periportal","LVEC"),
  HSC = c("HSC_Quiescent","HSC_Activated_ECM","HSC_Activated_Chemo",
          "HSC_Activated_Remodel","HSC_Proliferating"),
  KC  = c("Kupffer_KC"),
  Mac = c("Macrophage_M1","Macrophage_M2","Macrophage_MDM","Macrophage_NAM"),
  DC  = c("DC_cDC1","DC_cDC2","DC_pDC"),
  B   = c("B_Memory","B_Naive","B_Plasma"),
  T_NK= c("T_CD4","T_CD8","T_NKT","T_Treg","NK"),
  Granulocyte = c("Neutrophil","Basophil"),
  Cholangio   = c("Cholangiocyte"),
  Vasc        = c("Pericyte"),
  Erythrocyte = c("Erythrocyte"),
  LCM         = c("LCM")
)
LIVER_GROUP_COLOR <- c(
  Hep="#FFD700", LSEC="#1E90FF", HSC="#8B4513", KC="#FF6347",
  Mac="#FF1493", DC="#FF8C00", B="#9400D3", T_NK="#00CED1",
  Granulocyte="#32CD32", Cholangio="#DC143C", Vasc="#87CEEB",
  Erythrocyte="#8B0000", LCM="#FFB6C1"
)
BRAIN_GROUP_MAP <- list(
  Ex_Cortex = c("Excitatory_L23","Excitatory_L45","Excitatory_L6"),
  Inh_Cortex= c("Inhibitory_Meis2","Inhibitory_Pvalb","Inhibitory_Sst","Inhibitory_Vip"),
  Hpc       = c("Hpc_CA1","Hpc_CA3","Hpc_DG"),
  Oligo     = c("Oligo_Mature","OPCs"),
  Astro     = c("Astrocyte_Resting","Astrocyte_Gfap"),
  Microglia = c("Microglia_Homeostatic","Microglia_M1","Microglia_M2","Microglia_DAM_strict"),
  Endo      = c("Endothelial"),
  Pericyte  = c("Pericytes"),
  VLMC      = c("VLMC"),
  CPEC      = c("CPECs"),
  Ependymal = c("Ependymal")
)
BRAIN_GROUP_COLOR <- c(
  Ex_Cortex="#228B22", Inh_Cortex="#9ACD32", Hpc="#DC143C",
  Oligo="#FFD700", Astro="#8B4513", Microglia="#FF1493",
  Endo="#1E90FF", Pericyte="#87CEEB", VLMC="#9400D3",
  CPEC="#FF8C00", Ependymal="#FFB6C1"
)

# ============================================================================
# Helpers (调色)
# ============================================================================
shade_palette <- function(base_hex, n) {
  if (n == 1) return(base_hex)
  hsv_mat <- rgb2hsv(col2rgb(base_hex))
  h <- hsv_mat["h",1]; s <- hsv_mat["s",1]; v <- hsv_mat["v",1]
  v_seq <- seq(min(v*1.10,0.95), max(v*0.45,0.40), length.out=n)
  s_seq <- seq(max(s*0.95,0.85), 1.00, length.out=n)
  unname(mapply(function(vv,ss) hsv(h,ss,vv), v_seq, s_seq))
}
build_subtype_palette <- function(group_map, group_color) {
  pal <- c()
  for (g in names(group_map)) {
    cols <- shade_palette(group_color[g], length(group_map[[g]]))
    names(cols) <- group_map[[g]]; pal <- c(pal, cols)
  }
  pal
}
# 圆圈内编号文字颜色: 亮底用黑, 暗底用白
txt_on <- function(hex) {
  c <- col2rgb(hex)/255
  ifelse(0.299*c[1,]+0.587*c[2,]+0.114*c[3,] > 0.62, "black", "white")
}

# ============================================================================
# 右侧图例 (照例图: 多成员=组标题+编号圆圈; 单成员=圆圈+名字; 分 n_col 列)
# ============================================================================
build_right_legend <- function(leg_data, group_order, group_color,
                               n_col = 2L, circ_per_row = 3L, leg_span = 24,
                               title = "Cell type",
                               COL_DX = NULL, CIRC_DX = 1.15, INDENT = 0.5,
                               NAME_DX = 0.75, GAP = 0.08, CHAR_W = 0.46, COL_PAD = 0.9,
                               CIRC_SIZE = 10.5, NUM_SIZE = 5.0,
                               HDR_SIZE = 7.0, NAME_SIZE = 5.5, TITLE_SIZE = 8.5) {
  # COL_DX=NULL: 自动= 左列最长标题字符数*CHAR_W + COL_PAD (各器官适配, 标题不互压)
  # leg_span: 视觉密度(单位行映射到的总y跨度); 越大像素行距越小->越紧凑; 内容竖直居中
  sizes <- tapply(leg_data$idx, leg_data$supergroup, length)
  multi  <- group_order[sapply(group_order, function(g) sizes[[g]] > 1)]
  single <- group_order[sapply(group_order, function(g) sizes[[g]] == 1)]

  circ <- list(); hdr <- list()

  # ---- 左列(x=0): 多成员大类竖排 (标题 + 圆圈) ----
  y <- 0
  for (g in multi) {
    fillc <- unname(group_color[g])
    subs <- leg_data[leg_data$supergroup == g, ]; subs <- subs[order(subs$idx), ]
    hdr[[length(hdr)+1]] <- data.frame(x=0, y=y, label=g)
    for (j in seq_len(sizes[[g]])) {
      rr <- (j-1) %/% circ_per_row; cj <- (j-1) %% circ_per_row
      circ[[length(circ)+1]] <- data.frame(
        x=INDENT+cj*CIRC_DX, y=y-1-rr, idx=subs$idx[j],
        fillc=fillc, is_single=FALSE, name=NA_character_)
    }
    y <- y - 1 - ceiling(sizes[[g]]/circ_per_row) - GAP
  }
  left_bottom <- y

  # 列距: 左列内容最右端(标题宽 vs 圆圈行宽 取大) + 余量
  if (is.null(COL_DX)) {
    ext <- sapply(multi, function(g) {
      ncirc <- min(sizes[[g]], circ_per_row)
      max(nchar(g) * CHAR_W, INDENT + (ncirc - 1) * CIRC_DX + 0.5)
    })
    COL_DX <- (if (length(ext)) max(ext) else 3) + COL_PAD
  }

  # ---- 右列(x=COL_DX): 单成员长名字竖排 (圆圈 + 名字), 顶对齐 ----
  y <- 0
  for (g in single) {
    subs <- leg_data[leg_data$supergroup == g, ]
    circ[[length(circ)+1]] <- data.frame(
      x=COL_DX, y=y, idx=subs$idx[1], fillc=unname(group_color[g]),
      is_single=TRUE, name=subs$subtype[1])
    y <- y - 1 - GAP
  }

  circ_df <- do.call(rbind, circ); hdr_df <- do.call(rbind, hdr)
  circ_df$txtcol <- txt_on(circ_df$fillc)

  # 竖直居中(不拉满, 两器官同密度)
  y_top <- 1.6; y_bot <- min(circ_df$y) - 0.6
  content_h <- y_top - y_bot
  if (content_h > leg_span)
    warning(sprintf("[legend] 内容高 %.1f > leg_span %.1f, 调大 leg_span", content_h, leg_span))
  pad <- max(0, (leg_span - content_h) / 2)

  single_df <- circ_df[circ_df$is_single, ]
  # 右边界: 右列单成员名字末端
  name_chars <- if (nrow(single_df)) max(nchar(single_df$name)) else 0
  x_right <- COL_DX + NAME_DX + name_chars*CHAR_W + 0.4
  p_leg <- ggplot() +
    annotate("text", x=0, y=1.3, label=title, hjust=0, fontface="bold", size=TITLE_SIZE) +
    geom_point(data=circ_df, aes(x, y, fill=fillc), shape=21,
               color="white", stroke=0.7, size=CIRC_SIZE) +
    scale_fill_identity() +
    geom_text(data=circ_df, aes(x, y, label=idx, color=txtcol),
              size=NUM_SIZE, fontface="bold") +
    scale_color_identity() +
    { if (nrow(single_df)>0)
        geom_text(data=single_df, aes(x=x+NAME_DX, y=y, label=name),
                  hjust=0, size=NAME_SIZE, color="black") } +
    geom_text(data=hdr_df, aes(x, y, label=label),
              hjust=0, fontface="bold", size=HDR_SIZE, color="black") +
    scale_x_continuous(limits=c(-0.7, x_right), expand=expansion(0)) +
    scale_y_continuous(limits=c(y_bot - pad, y_top + pad), expand=expansion(0)) +
    coord_cartesian(clip="off") +
    theme_void() + theme(plot.margin = margin(6, 4, 6, 8))
  list(gg = p_leg, x_extent = x_right + 0.7)   # x_extent 供外层按内容缩放列宽
}

# ============================================================================
# 主图 (箭头轴 + 大字 + theme_void) + 右侧图例 -> 横向拼
# ============================================================================
plot_paper_umap_v4 <- function(obj, group_map, group_color, organ_name,
                               n_col = 2L, leg_span = 24) {
  pal <- build_subtype_palette(group_map, group_color)
  obj$celltype_refined <- as.character(obj$celltype_refined)
  flat <- stack(group_map); colnames(flat) <- c("subtype","supergroup")
  flat$subtype <- as.character(flat$subtype); flat$supergroup <- as.character(flat$supergroup)
  obj$supergroup <- flat$supergroup[match(obj$celltype_refined, flat$subtype)]

  emb <- as.data.frame(Embeddings(obj, "umap")); colnames(emb) <- c("UMAP1","UMAP2")
  emb$subtype <- obj$celltype_refined; emb$supergroup <- obj$supergroup
  emb <- emb[!is.na(emb$subtype) & emb$subtype %in% names(pal), ]
  emb$subtype <- factor(emb$subtype, levels=names(pal))

  centroids_g <- emb %>% group_by(supergroup) %>%
    summarise(UMAP1=median(UMAP1), UMAP2=median(UMAP2), .groups="drop")
  centroids_s <- emb %>% mutate(subtype=as.character(subtype)) %>%
    group_by(subtype, supergroup) %>%
    summarise(UMAP1=median(UMAP1), UMAP2=median(UMAP2), n=n(), .groups="drop") %>%
    arrange(supergroup, desc(n)) %>% mutate(idx=row_number())
  emb_ordered <- emb %>% group_by(subtype) %>%
    mutate(nn=n()) %>% ungroup() %>% arrange(desc(nn))

  # 箭头轴几何 (左下角)
  rx <- range(emb$UMAP1); ry <- range(emb$UMAP2)
  dx <- diff(rx); dy <- diff(ry)
  ox <- rx[1] - 0.05*dx; oy <- ry[1] - 0.05*dy
  ax <- 0.17*dx; ay <- 0.17*dy
  ar <- arrow(length=unit(0.28,"cm"), type="closed")

  p <- ggplot(emb_ordered, aes(UMAP1, UMAP2, color=subtype)) +
    geom_point(size=0.25, alpha=0.95, stroke=0) +
    scale_color_manual(values=pal, breaks=names(pal), limits=names(pal),
                       drop=FALSE, na.value="#FF00FF", guide="none") +
    scale_x_continuous(expand=expansion(mult=c(0.10, 0.03))) +
    scale_y_continuous(expand=expansion(mult=c(0.10, 0.03))) +
    ggrepel::geom_label_repel(data=centroids_s,
               aes(x=UMAP1, y=UMAP2, label=idx, fill=supergroup),
               color="white", fontface="bold", size=3.0, label.size=0,
               label.padding=unit(0.09,"lines"), label.r=unit(0.3,"lines"),
               box.padding=0.20, point.padding=0.05, force=2.2, force_pull=0.8,
               min.segment.length=0, segment.color="grey55", segment.size=0.3,
               max.overlaps=Inf, seed=7,
               show.legend=FALSE, inherit.aes=FALSE) +
    scale_fill_manual(values=group_color, guide="none") +
    ggrepel::geom_text_repel(
      data=centroids_g, aes(x=UMAP1, y=UMAP2, label=supergroup),
      color="black", size=7.5, fontface="bold",     # [v4] 5.5->7.5 (+2号)
      bg.color="white", bg.r=0.18, box.padding=0.9, point.padding=0.6,
      force=3.5, force_pull=0.3, min.segment.length=0,
      segment.color="grey60", segment.size=0.35, max.overlaps=Inf,
      seed=42, inherit.aes=FALSE) +
    # 箭头坐标轴
    annotate("segment", x=ox, xend=ox+ax, y=oy, yend=oy, arrow=ar, linewidth=0.7) +
    annotate("segment", x=ox, xend=ox, y=oy, yend=oy+ay, arrow=ar, linewidth=0.7) +
    annotate("text", x=ox+ax*0.5, y=oy-0.045*dy, label="UMAP1",
             hjust=0.5, vjust=1, size=5.0, fontface="bold") +
    annotate("text", x=ox-0.04*dx, y=oy+ay*0.5, label="UMAP2",
             angle=90, hjust=0.5, vjust=0, size=5.0, fontface="bold") +
    coord_cartesian(clip="off") +
    theme_void() +
    theme(plot.margin = margin(10, 6, 10, 10))

  leg_data <- centroids_s %>% mutate(supergroup=as.character(supergroup)) %>%
    select(idx, subtype, supergroup)
  leg_res <- build_right_legend(leg_data, names(group_map), group_color,
                                n_col=n_col, leg_span=leg_span)
  p_leg <- leg_res$gg

  # 图例列宽 = x跨度 * 单位英寸(两器官同 UNIT_IN -> 圆圈大小一致); 整图宽固定
  UNIT_IN <- 0.34; FIG_W <- 12
  leg_in  <- leg_res$x_extent * UNIT_IN
  main_in <- FIG_W - leg_in
  combined <- p | p_leg
  combined <- combined + plot_layout(widths = c(main_in, leg_in))  # [align] 同 UNIT_IN
  list(main=p, legend=p_leg, combined=combined, centroids_s=centroids_s)
}

# ============================================================================
# Run  (两器官同尺寸 14x9, 同 widths -> 主面板等高同位, 可拼接)
# ============================================================================
LEG_SPAN <- 24   # 越大行距越小; 须 > 肝左列内容高(~21)

cat("==== LIVER ====\n")
liver_obj <- readRDS(file.path(OUT_DIR, "liver_integrated.rds"))
liver_plot <- plot_paper_umap_v4(liver_obj, LIVER_GROUP_MAP, LIVER_GROUP_COLOR,
                                 "Liver", n_col=2L, leg_span=LEG_SPAN)
save_fig_all(liver_plot$combined, file.path(RESTYLE_DIR, "liver_umap_grouped_v3_9"),
             width=12, height=9)
write.csv(liver_plot$centroids_s,
          file.path(OUT_DIR, "liver_subtype_centroids_v3_9.csv"), row.names=FALSE)
rm(liver_obj, liver_plot); gc(verbose=FALSE)
cat("[Liver] done.\n\n")

cat("==== BRAIN ====\n")
brain_obj <- readRDS(file.path(OUT_DIR, "brain_integrated.rds"))
brain_plot <- plot_paper_umap_v4(brain_obj, BRAIN_GROUP_MAP, BRAIN_GROUP_COLOR,
                                 "Brain", n_col=2L, leg_span=LEG_SPAN)
save_fig_all(brain_plot$combined, file.path(RESTYLE_DIR, "brain_umap_grouped_v3_9"),
             width=12, height=9)
write.csv(brain_plot$centroids_s,
          file.path(OUT_DIR, "brain_subtype_centroids_v3_9.csv"), row.names=FALSE)
cat("[Brain] done.\n")
