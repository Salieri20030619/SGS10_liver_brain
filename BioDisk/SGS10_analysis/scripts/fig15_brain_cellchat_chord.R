\a/#!/usr/bin/env Rscript
# ============================================================================
# fig15_brain_cellchat_chord.R  (v2: diff_T_vs_B ranked, Excit_L45 highlighted)
# [restyle] 标签放大 0.55->0.85, 标题 1.2->1.4, 输出 restyle_v1, PNG 300dpi,
#           统一走 save_fig_all (cairo PDF + TIFF/PNG 300). 分析逻辑不变.
# Brain CellChat 3-group chord. Top neuronal LR pairs are housekeeping
# (Vip/Excit, ~7-8 weight, nearly identical across groups). To show the
# TUMOR-driven changes, we rank by diff_T_vs_B (Tumor weight - Blank weight)
# and pick top 25. This surfaces Excit_L45 as the convergence point.
# All 3 panels use the SAME 25 pairs (selected on Tumor side).
# ============================================================================
suppressPackageStartupMessages({
  library(circlize)
  library(dplyr)
})
source("/mnt/BioDisk/wangmingyuan/SGS10_analysis/scripts/paper_theme.R")
source("/mnt/BioDisk/wangmingyuan/SGS10_analysis/scripts/fig_save.R")   # [restyle]

csv_path <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis/taskB/cellchat/3groups_plots_v3/brain_3groups_top_pairs.csv"
df <- read.csv(csv_path, stringsAsFactors = FALSE)
cat(sprintf("[fig15] %d total pairs loaded\n", nrow(df)))

N_TOP <- 25

# Select SAME pair set across 3 panels: rank by diff_T_vs_B desc
selected <- df %>%
  arrange(desc(diff_T_vs_B)) %>%
  slice_head(n = N_TOP) %>%
  pull(key)

df_sel <- df %>% filter(key %in% selected)
cat(sprintf("[fig15] selected %d pairs (ranked by diff_T_vs_B)\n", nrow(df_sel)))
cat("[fig15] senders involved:\n")
print(table(df_sel$sender))
cat("[fig15] receivers involved:\n")
print(table(df_sel$receiver))

draw_chord <- function(df, weight_col, group_name, paper_key_sender) {
  d <- df[, c("sender", "receiver", weight_col)]
  colnames(d) <- c("sender", "receiver", "weight")
  d <- d %>% filter(weight > 0)
  cat(sprintf("  [%s] %d pairs\n", group_name, nrow(d)))

  all_ct <- sort(unique(c(d$sender, d$receiver)))
  n_ct <- length(all_ct)
  ct_colors <- setNames(
    colorRampPalette(c("#6A5ACD", "#7BAFD4", "#A8D8E5", "#D6BDA1",
                      "#E66D50", "#3B82C4"))(n_ct),
    all_ct
  )
  if (paper_key_sender %in% all_ct) {
    ct_colors[paper_key_sender] <- "#C0143B"
  }
  d$ribbon_col <- ct_colors[d$sender]

  circos.clear()
  circos.par(gap.after = 2.5, start.degree = 90,
             canvas.xlim = c(-1.45, 1.45), canvas.ylim = c(-1.45, 1.45),
             points.overflow.warning = FALSE)
  chordDiagram(
    d,
    grid.col = ct_colors,
    col = d$ribbon_col,
    transparency = 0.35,
    directional = 1,
    direction.type = c("arrows", "diffHeight"),
    diffHeight = -0.04,
    link.arr.type = "big.arrow",
    annotationTrack = "grid",
    preAllocateTracks = list(track.height = 0.10),
    big.gap = 12
  )
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    sector.name <- get.cell.meta.data("sector.index")
    xlim <- get.cell.meta.data("xlim")
    ylim <- get.cell.meta.data("ylim")
    is_key <- sector.name == paper_key_sender
    circos.text(mean(xlim), ylim[1] + 0.5, sector.name,
                facing = "clockwise", niceFacing = TRUE,
                adj = c(0, 0.5),
                cex = 0.85,                       # [restyle] 0.55->0.85
                font = if (is_key) 2 else 1,
                col = if (is_key) "#C0143B" else "black")
  }, bg.border = NA)

  title(group_name, cex.main = 1.4, font.main = 2,    # [restyle] 1.2->1.4
        col.main = "black", line = -0.5)
}

out_dir  <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis/figures/restyle_v1"  # [restyle]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
basename <- "fig15_brain_cellchat_chord"
w <- 13; h <- 5

paper_key <- "Excitatory_L45"

render_all <- function() {
  par(mfrow = c(1, 3), mar = c(2, 2, 3, 2))
  draw_chord(df_sel, "weight_Tumor",   "Tumor",   paper_key)
  draw_chord(df_sel, "weight_Blank",   "Blank",   paper_key)
  draw_chord(df_sel, "weight_Control", "Control", paper_key)
}

# [restyle] 统一出 PDF(cairo) + TIFF(300,lzw) + PNG(300), 白底
save_fig_all(drawfun = render_all,
             file = file.path(out_dir, basename),
             width = w, height = h, dpi = 300)

cat(sprintf("[OK] %s  %.1f x %.1f in  ->  restyle_v1 PDF+TIFF+PNG\n",
            basename, w, h))
