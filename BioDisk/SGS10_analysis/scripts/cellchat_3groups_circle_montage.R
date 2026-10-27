#!/usr/bin/env Rscript
# ============================================================================
# cellchat_3groups_circle_montage.R   [restyle, replot-only]
# 重画 docx img15(brain)/img22(liver) 的 6-panel 聚合通信圈图.
#   布局: 2x3 -> 上排 count(Number of interactions), 下排 weight(Interaction weight)
#         列: Tumor / Blank / Control
#   只读已存的 cellchat 对象, 不重算 (秒出); 三条件统一配色; 放大标签.
#   输出: figures/restyle_v1/<organ>_3groups_circle.{pdf,tiff,png}
# ============================================================================
suppressPackageStartupMessages({
  library(CellChat)
  library(dplyr)
})
source("/mnt/BioDisk/wangmingyuan/SGS10_analysis/scripts/fig_save.R")

ROOT    <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
CC_DIR  <- file.path(ROOT, "taskB", "cellchat")
OUT_DIR <- file.path(ROOT, "figures", "restyle_v1")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# 显示名 -> rds 文件中缀
CONDS    <- c(Tumor = "tumor", Blank = "blank", Control = "control")
VLAB_CEX <- 0.8     # 节点标签字号(可调; 太挤就调小, 想更大就调大)

# ---- 画单个圈图(count 或 weight) ----
draw_one <- function(cc, mat_type, organ_cap, cond_disp, cols) {
  mat <- if (mat_type == "count") cc@net$count else cc@net$weight
  vw  <- as.numeric(table(cc@idents))   # 节点大小 = 各类型 spot 数
  ttl <- if (mat_type == "count")
    sprintf("%s %s - Number of interactions (n=%d)", organ_cap, cond_disp, sum(mat))
  else
    sprintf("%s %s - Interaction weight (total=%.2f)", organ_cap, cond_disp, sum(mat))

  netVisual_circle(
    mat,
    vertex.weight    = vw,
    weight.scale     = TRUE,
    label.edge       = FALSE,
    color.use        = cols[rownames(mat)],   # 统一配色(按名字取, 跨条件一致)
    vertex.label.cex = VLAB_CEX,
    title.name       = ttl
  )
}

# ---- 一个器官的 2x3 拼图 ----
make_montage <- function(organ, organ_cap) {
  cat(sprintf("\n==== %s ====\n", organ_cap))
  ccs <- lapply(CONDS, function(fn) {
    f <- file.path(CC_DIR, sprintf("%s_%s_cellchat.rds", organ, fn))
    if (!file.exists(f)) stop(sprintf("缺对象: %s", f))
    cat(sprintf("  读 %s\n", basename(f)))
    readRDS(f)
  })
  names(ccs) <- names(CONDS)

  # 三条件细胞类型并集 -> 统一配色(CellChat 默认调色板)
  all_ct <- sort(unique(unlist(lapply(ccs, function(cc) rownames(cc@net$count)))))
  cols   <- setNames(CellChat::scPalette(length(all_ct)), all_ct)
  cat(sprintf("  细胞类型并集: %d 种\n", length(all_ct)))

  render <- function() {
    par(mfrow = c(2, 3), mar = c(1, 1, 3, 1), xpd = TRUE)
    for (cond in names(ccs)) draw_one(ccs[[cond]], "count",  organ_cap, cond, cols)
    for (cond in names(ccs)) draw_one(ccs[[cond]], "weight", organ_cap, cond, cols)
  }

  save_fig_all(drawfun = render,
               file   = file.path(OUT_DIR, sprintf("%s_3groups_circle", organ)),
               width  = 18, height = 12, dpi = 300)
  cat(sprintf("  [OK] %s_3groups_circle -> restyle_v1 (PDF+TIFF+PNG)\n", organ))
  rm(ccs); gc(verbose = FALSE)
}

make_montage("liver", "Liver")
make_montage("brain", "Brain")

cat("\nDone.\n")
