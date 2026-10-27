#!/usr/bin/env Rscript
# =============================================================================
# SGS10 D3 step3: paper-quality figures (重做 D3 visualization)
# -----------------------------------------------------------------------------
# 输入: D3_signature_output/{spot_scores_*.csv, animal_summary.csv,
#                            filtered_DEG/*.csv}
#
# 输出 (D3_signature_output/figures_paper/):
#   F1_signature_heatmap_paper.pdf      ⭐ ComplexHeatmap row-scaled
#   F2_stage_trend_per_signature.pdf    每个 sig 独立图,大字体,baseline reference
#   F3_per_animal_paired.pdf            paired-line 而非 scatter, 标 stage 方向
#   F4_signature_violins.pdf            spot-level 分布 (远比 mean 信息多)
#   F5_top_genes_dotplot.pdf            top 20 真 UP × samples dotplot
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(scales)
  library(RColorBrewer)
  library(ComplexHeatmap)
  library(circlize)
})

option_list <- list(
  make_option("--d3_dir", type="character",
              default="/mnt/BioDisk/wangmingyuan/SGS10_analysis/D3_signature_output"),
  make_option("--out_subdir", type="character", default="figures_paper")
)
opt <- parse_args(OptionParser(option_list=option_list))

CFG <- list(
  d3_dir = opt$d3_dir,
  out_dir = file.path(opt$d3_dir, opt$out_subdir)
)

dir.create(CFG$out_dir, showWarnings=FALSE, recursive=TRUE)

cat_log <- function(...) cat(..., sep="")
sec <- function(name) {
  cat_log("\n", strrep("-", 70), "\n  ", name, "\n",
          strrep("-", 70), "\n")
}

cat_log("=====================================================================\n")
cat_log("D3 Step 3: paper-quality figures\n")
cat_log("Start: ", format(Sys.time()), "\n")
cat_log("=====================================================================\n")

# ============ Color palettes ============
zone_colors <- c(
  "tumor_core"    = "#d7191c",
  "tumor_margin"  = "#fdae61",
  "invaded_host"  = "#ffffbf",
  "distant_host"  = "#abd9e9",
  "baseline_host" = "#2c7bb6"
)
stage_colors <- c(
  "baseline" = "#999999",
  "early"    = "#5e4fa2",
  "mid"      = "#fdae61",
  "late"     = "#d7191c"
)
sig_colors <- c(
  "A2" = "#d7191c",   # liver invaded
  "A3" = "#fdae61",   # liver distant
  "B3" = "#9ecae1",   # liver progression
  "C4" = "#2c7bb6"    # brain
)

# ============ Load inputs ============
sec("Load D3 outputs")
animal_df <- read_csv(file.path(CFG$d3_dir, "animal_summary.csv"),
                      show_col_types=FALSE)
cat_log("  animal_summary: ", nrow(animal_df), " rows\n")

spot_liver <- read_csv(file.path(CFG$d3_dir, "spot_scores_liver.csv"),
                       show_col_types=FALSE)
spot_brain <- read_csv(file.path(CFG$d3_dir, "spot_scores_brain.csv"),
                       show_col_types=FALSE)
cat_log("  spot scores liver: ", nrow(spot_liver), "\n")
cat_log("  spot scores brain: ", nrow(spot_brain), "\n")

# Stage factor levels
stage_levels <- c("baseline", "early", "mid", "late")
animal_df$stage <- factor(animal_df$stage, levels = stage_levels)
spot_liver$stage <- factor(spot_liver$stage, levels = stage_levels)
spot_brain$stage <- factor(spot_brain$stage, levels = stage_levels)

zone_levels <- c("tumor_core","tumor_margin","invaded_host",
                 "distant_host","baseline_host")
animal_df$zone <- factor(animal_df$zone, levels = zone_levels)

# ============ FIG 1: Signature heatmap (paper) ============
sec("F1: Signature heatmap (ComplexHeatmap, row-scaled per signature)")

# Build matrix: rows=signatures, cols=cell_id (zone__stage)
# But better: cols=samples, rows=signatures × zones
# We use animal_summary to get per-(animal, zone, signature, organ) means

heat_input <- animal_df %>%
  filter(!is.na(stage)) %>%
  mutate(col_id = paste0(animal_id, "\n", zone),
         row_id = paste0(organ, ":", signature)) %>%
  select(col_id, row_id, mean_score, organ, signature, zone, stage, animal_id, n_spots)

# Make wide matrix
mat_df <- heat_input %>%
  select(row_id, col_id, mean_score) %>%
  pivot_wider(names_from = col_id, values_from = mean_score, values_fill = NA)
mat <- as.matrix(mat_df[, -1])
rownames(mat) <- mat_df$row_id

# Column annotation
col_meta <- heat_input %>% distinct(col_id, organ, zone, stage, animal_id, n_spots)
col_meta <- col_meta[match(colnames(mat), col_meta$col_id), ]

# Order columns: organ, then stage, then zone
col_order <- col_meta %>%
  mutate(stage_n = as.numeric(stage),
         zone_n = as.numeric(factor(zone, levels = zone_levels))) %>%
  arrange(organ, stage_n, zone_n) %>%
  pull(col_id)
mat <- mat[, col_order]
col_meta <- col_meta[match(colnames(mat), col_meta$col_id), ]

# Row order: organ, then signature
row_order <- order(sub(":.*", "", rownames(mat)),
                    sub(".*:", "", rownames(mat)))
mat <- mat[row_order, ]

# Row-scale per signature for visibility (z-score across cols within row)
mat_scaled <- t(apply(mat, 1, function(r) {
  rs <- (r - mean(r, na.rm = TRUE)) / sd(r, na.rm = TRUE)
  if (all(is.na(rs))) rep(0, length(r)) else rs
}))
rownames(mat_scaled) <- rownames(mat)
colnames(mat_scaled) <- colnames(mat)

# Top annotation: organ + stage + zone
ha_top <- HeatmapAnnotation(
  organ = col_meta$organ,
  stage = as.character(col_meta$stage),
  zone  = as.character(col_meta$zone),
  n_spots = anno_barplot(col_meta$n_spots, gp = gpar(fill = "gray80"),
                         border = FALSE, height = unit(1.2, "cm")),
  col = list(
    organ = c("liver" = "#d8b365", "brain" = "#5ab4ac"),
    stage = stage_colors[stage_levels],
    zone = zone_colors[zone_levels]
  ),
  show_annotation_name = TRUE,
  annotation_name_side = "left",
  annotation_legend_param = list(
    organ = list(title = "organ"),
    stage = list(title = "stage"),
    zone = list(title = "zone")
  ),
  gp = gpar(col = "white")
)

# Row split by organ
row_split <- sub(":.*", "", rownames(mat_scaled))

col_fun <- colorRamp2(c(-2, 0, 2), c("#3182bd", "white", "#e34a33"))

pdf(file.path(CFG$out_dir, "F1_signature_heatmap_paper.pdf"),
    width = 14, height = 6)
ht <- Heatmap(
  mat_scaled,
  name = "z-score\n(per row)",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = row_split,
  top_annotation = ha_top,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 11),
  rect_gp = gpar(col = "white", lwd = 0.5),
  column_title = "Signature scores (z-score per signature) - paper figure",
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 9),
    labels_gp = gpar(fontsize = 8)
  )
)
draw(ht, merge_legend = TRUE)
dev.off()
cat_log("  saved F1\n")

# Also save raw-scale heatmap (without row-scaling) for completeness
pdf(file.path(CFG$out_dir, "F1b_signature_heatmap_rawscale.pdf"),
    width = 14, height = 6)
col_fun_raw <- colorRamp2(
  c(0, max(mat, na.rm = TRUE) * 0.3, max(mat, na.rm = TRUE)),
  c("white", "#fdae61", "#a50026"))
ht2 <- Heatmap(
  mat,
  name = "UCell\nscore",
  col = col_fun_raw,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = row_split,
  top_annotation = ha_top,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 11),
  rect_gp = gpar(col = "white", lwd = 0.5),
  column_title = "Signature scores (raw UCell)",
  column_title_gp = gpar(fontsize = 12, fontface = "bold")
)
draw(ht2, merge_legend = TRUE)
dev.off()
cat_log("  saved F1b (raw scale)\n")

# ============ FIG 2: Per-signature stage trend (separate panels) ============
sec("F2: Stage trend per signature (separated, with baseline reference)")

trend_plots <- list()
sig_titles <- list(
  A2 = "A2: Liver invaded_host signature\n(Fmo3, Cyp450, Cd74, ...)",
  A3 = "A3: Liver distant_host signature\n(Cyp450, Slco1a4, ...)",
  B3 = "B3: Liver progression late vs early\n(Mup, Gm pseudogenes — sample-driven)",
  C4 = "C4: Brain distant_host signature\n(Nr2f1, Apoe, Lrp1, Mef2c, ...)"
)

for (sig in names(sig_titles)) {
  sub_data <- animal_df %>%
    filter(signature == sig & !is.na(stage))
  if (nrow(sub_data) == 0) next
  
  # Compute baseline reference (per organ) for shading
  baseline_ref <- sub_data %>%
    filter(stage == "baseline") %>%
    group_by(organ) %>%
    summarise(bl_mean = mean(mean_score, na.rm = TRUE),
              bl_min = min(mean_score, na.rm = TRUE),
              bl_max = max(mean_score, na.rm = TRUE),
              .groups = "drop")
  
  # Plot per organ
  organ_plots <- list()
  for (org in unique(sub_data$organ)) {
    df_org <- sub_data %>% filter(organ == org)
    bl <- baseline_ref %>% filter(organ == org)
    
    p <- ggplot(df_org, aes(x = stage, y = mean_score)) +
      # Baseline reference band
      {if (nrow(bl) > 0)
        geom_hline(yintercept = bl$bl_mean, color = "gray60",
                   linetype = "dashed", linewidth = 0.4)} +
      {if (nrow(bl) > 0)
        annotate("rect", xmin = -Inf, xmax = Inf,
                 ymin = bl$bl_min, ymax = bl$bl_max,
                 fill = "gray", alpha = 0.15)} +
      # Per-zone lines (only relevant for liver where multiple zones exist)
      geom_line(aes(color = zone, group = zone),
                linewidth = 1.0, alpha = 0.85) +
      geom_point(aes(color = zone, size = n_spots),
                 alpha = 0.9) +
      scale_color_manual(values = zone_colors, drop = FALSE) +
      scale_size(range = c(2, 6), name = "n_spots", guide = "none") +
      labs(title = paste0(org, " - ", sig),
           x = NULL, y = "mean UCell score") +
      theme_classic(base_size = 11) +
      theme(plot.title = element_text(size = 10, face = "bold"),
            legend.position = "right",
            legend.text = element_text(size = 8),
            legend.title = element_text(size = 9))
    organ_plots[[org]] <- p
  }
  
  combined <- wrap_plots(organ_plots, ncol = length(organ_plots)) +
    plot_annotation(title = sig_titles[[sig]],
                    theme = theme(plot.title = element_text(size = 12,
                                                              face = "bold",
                                                              hjust = 0.5)))
  trend_plots[[sig]] <- combined
}

# Stack all signatures vertically
if (length(trend_plots) > 0) {
  full_trend <- wrap_plots(trend_plots, ncol = 1)
  ggsave(file.path(CFG$out_dir, "F2_stage_trend_per_signature.pdf"),
         full_trend,
         width = 12, height = 4 * length(trend_plots),
         limitsize = FALSE)
  cat_log("  saved F2\n")
}

# ============ FIG 3: Per-animal paired (replaces ugly correlation scatter) ============
sec("F3: Per-animal paired liver-brain (line plot replaces N=3 scatter)")

# For each signature, plot 3 paired animals as connected lines
# X axis: organ (liver -> brain), Y: score
# Color: stage (early/mid/late)

paired_animals <- c("animal_8", "animal_3", "animal_56")
stage_for_animal <- c(animal_8 = "early", animal_3 = "mid", animal_56 = "late")

paired_panels <- list()
for (sig in c("A2", "A3", "C4")) {  # B3 skipped (within-organ)
  # Pick the relevant zones
  liver_zone <- if (sig %in% c("A1","A2")) "invaded_host" else "distant_host"
  brain_zone <- "distant_host"
  
  ld <- animal_df %>%
    filter(organ == "liver", signature == sig,
           zone == liver_zone, animal_id %in% paired_animals) %>%
    select(animal_id, liver_score = mean_score)
  bd <- animal_df %>%
    filter(organ == "brain", signature == sig,
           zone == brain_zone, animal_id %in% paired_animals) %>%
    select(animal_id, brain_score = mean_score)
  
  paired <- inner_join(ld, bd, by = "animal_id")
  if (nrow(paired) < 2) next
  paired$stage <- stage_for_animal[paired$animal_id]
  paired$stage <- factor(paired$stage, levels = c("early","mid","late"))
  
  long <- paired %>%
    pivot_longer(c(liver_score, brain_score),
                 names_to = "organ_score",
                 values_to = "score") %>%
    mutate(organ = ifelse(organ_score == "liver_score",
                          paste0("liver\n(", liver_zone, ")"),
                          paste0("brain\n(", brain_zone, ")")),
           organ = factor(organ, levels = c(paste0("liver\n(", liver_zone, ")"),
                                              paste0("brain\n(", brain_zone, ")"))))
  
  p <- ggplot(long, aes(x = organ, y = score,
                         group = animal_id, color = stage)) +
    geom_line(linewidth = 1.2, alpha = 0.7) +
    geom_point(size = 4) +
    geom_text_repel(aes(label = animal_id), size = 3, show.legend = FALSE,
                     min.segment.length = 0) +
    scale_color_manual(values = stage_colors[c("early","mid","late")]) +
    labs(title = paste0(sig, " signature - paired liver-brain (3 animals)"),
         subtitle = paste0("liver: ", liver_zone, " · brain: ", brain_zone,
                           " · paper: hypothesis only (N=3)"),
         x = NULL, y = "mean UCell score per animal") +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(size = 9, color = "gray40"))
  paired_panels[[sig]] <- p
}

if (length(paired_panels) > 0) {
  full_paired <- wrap_plots(paired_panels, ncol = length(paired_panels))
  ggsave(file.path(CFG$out_dir, "F3_per_animal_paired.pdf"),
         full_paired,
         width = 5 * length(paired_panels), height = 5,
         limitsize = FALSE)
  cat_log("  saved F3\n")
}

# ============ FIG 4: Spot-level violin (full distribution, not just mean) ============
sec("F4: Spot-level signature score violins")

violin_panels <- list()

for (sig in c("A2","A3","B3","C4")) {
  score_col <- paste0("score_", sig)
  if (!score_col %in% colnames(spot_liver)) next
  
  # Build long df from both organs
  liver_long <- spot_liver %>%
    select(sample, zone, stage, animal_id, score = !!score_col) %>%
    mutate(organ = "liver")
  brain_long <- spot_brain %>%
    select(sample, zone, stage, animal_id, score = !!score_col) %>%
    mutate(organ = "brain")
  
  combined <- bind_rows(liver_long, brain_long) %>%
    filter(!is.na(score), !is.na(zone)) %>%
    mutate(zone = factor(zone, levels = zone_levels),
           stage = factor(stage, levels = stage_levels),
           cell_id = paste0(zone, "__", stage))
  
  # Subsample to max 5000 spots per cell_id for plotting speed
  # (slice_sample n must be constant; group_modify handles per-group size)
  set.seed(42)
  combined <- combined %>%
    group_by(organ, cell_id) %>%
    group_modify(~ {
      if (nrow(.x) > 5000) slice_sample(.x, n = 5000) else .x
    }) %>%
    ungroup()
  
  # Order x by organ then by cell sequence
  cell_order_liver <- expand.grid(stage = stage_levels, zone = zone_levels) %>%
    mutate(cell_id = paste0(zone, "__", stage)) %>%
    pull(cell_id)
  combined$cell_id <- factor(combined$cell_id,
                              levels = intersect(cell_order_liver, combined$cell_id))
  
  p <- ggplot(combined, aes(x = cell_id, y = score, fill = zone)) +
    geom_violin(scale = "width", alpha = 0.85, color = NA) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.4,
                 fill = "white") +
    scale_fill_manual(values = zone_colors, drop = FALSE) +
    facet_grid(. ~ organ, scales = "free_x", space = "free_x") +
    labs(title = paste0(sig, " - per-spot UCell score distribution"),
         x = "zone × stage cell", y = "UCell score") +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          legend.position = "bottom",
          strip.background = element_rect(fill = "gray95"),
          strip.text = element_text(face = "bold"))
  violin_panels[[sig]] <- p
}

if (length(violin_panels) > 0) {
  full_violins <- wrap_plots(violin_panels, ncol = 1)
  ggsave(file.path(CFG$out_dir, "F4_signature_violins.pdf"),
         full_violins,
         width = 14, height = 4 * length(violin_panels),
         limitsize = FALSE)
  cat_log("  saved F4\n")
}

# ============ FIG 5: Top genes × samples dotplot ============
sec("F5: Top genes per signature × samples dotplot")

# Need raw spot data for per-sample expression of top genes.
# Skip if too heavy (would need to reload raw counts). Instead use a proxy:
# Cross-reference filtered_DEG csv (which has per-gene log2FC + pct already)

dotplot_panels <- list()
for (sig in c("A2", "C4")) {  # 2 most important
  fp <- file.path(CFG$d3_dir, "filtered_DEG", paste0(sig, "_UP_filtered.csv"))
  if (!file.exists(fp)) next
  df <- read_csv(fp, show_col_types = FALSE)
  
  top_n <- min(20, nrow(df))
  top_genes <- df %>% arrange(desc(avg_log2FC)) %>% head(top_n)
  
  # Build a simple plot: log2FC vs pct.1 with gene labels
  p <- ggplot(top_genes,
              aes(x = avg_log2FC, y = reorder(gene, avg_log2FC),
                  size = pct.1, color = pct.1 - pct.2)) +
    geom_point() +
    scale_color_gradient2(low = "#3182bd", mid = "gray80", high = "#e34a33",
                          midpoint = 0,
                          name = "Δ pct\n(grp1 - grp2)") +
    scale_size(range = c(2, 7), name = "pct.1\n(% expressing)") +
    labs(title = paste0(sig, " - Top ", top_n, " UP genes (null-filtered)"),
         x = "avg_log2FC (Wilcoxon)", y = NULL) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          axis.text.y = element_text(size = 9))
  dotplot_panels[[sig]] <- p
}

if (length(dotplot_panels) > 0) {
  full_dot <- wrap_plots(dotplot_panels, ncol = length(dotplot_panels))
  ggsave(file.path(CFG$out_dir, "F5_top_genes_dotplot.pdf"),
         full_dot,
         width = 6 * length(dotplot_panels), height = 7,
         limitsize = FALSE)
  cat_log("  saved F5\n")
}

# ============ Summary ============
cat_log("\n=== DONE: ", format(Sys.time()), " ===\n")
cat_log("\nGenerated paper-quality figures in: ", CFG$out_dir, "\n")
cat_log("  F1_signature_heatmap_paper.pdf       (z-scored, ComplexHeatmap)\n")
cat_log("  F1b_signature_heatmap_rawscale.pdf   (raw UCell)\n")
cat_log("  F2_stage_trend_per_signature.pdf     (separate sigs, baseline ref)\n")
cat_log("  F3_per_animal_paired.pdf             (paired-line replaces scatter)\n")
cat_log("  F4_signature_violins.pdf             (full spot distribution)\n")
cat_log("  F5_top_genes_dotplot.pdf             (A2 & C4 top 20 genes)\n\n")
