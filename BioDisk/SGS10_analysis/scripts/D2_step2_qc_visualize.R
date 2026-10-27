#!/usr/bin/env Rscript
# =============================================================================
# SGS10 D2 step2: QC + visualization for zone × stage assignment
# -----------------------------------------------------------------------------
# Input:  D2_zone_output/spot_zone_stage_<organ>.rds  (from step1)
#         spatial_results/rctd/per_sample/*_rctd_results.rds  (for celltype)
#
# Outputs:
#   D2_zone_output/figures/
#     - zone_map_<organ>.pdf           : 5-color zone map (per sample)
#     - zone_dist_hist_<organ>.pdf     : min_dist_to_tumor histogram by zone
#     - zone_stage_trend_<organ>.pdf   : zone counts vs stage (line + bar)
#     - cell_celltype_heatmap_<organ>.pdf : 13-cell × top celltypes
#   D2_zone_output/
#     - cell_summary_<organ>.csv       : per-cell n_spots + top celltypes
#     - cell_celltype_composition_<organ>.csv : per-cell × per-celltype mean prop
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

option_list <- list(
  make_option("--organ",    type="character", default="liver"),
  make_option("--base_dir", type="character",
              default="/mnt/BioDisk/wangmingyuan/SGS10_analysis"),
  make_option("--d2_dir",   type="character",
              default="/mnt/BioDisk/wangmingyuan/SGS10_analysis/D2_zone_output"),
  make_option("--top_n_celltypes", type="integer", default=10,
              help="show top N celltypes in heatmap"),
  make_option("--rctd_intersect_min", type="integer", default=20,
              help="minimum RCTD-matched spots per cell to include in celltype analysis")
)
opt <- parse_args(OptionParser(option_list=option_list))

CFG <- list(
  organ = opt$organ,
  base_dir = opt$base_dir,
  d2_dir = opt$d2_dir,
  rctd_per_sample_dir = file.path(opt$base_dir, "spatial_results/rctd/per_sample"),
  fig_dir = file.path(opt$d2_dir, "figures"),
  top_n_celltypes = opt$top_n_celltypes,
  rctd_intersect_min = opt$rctd_intersect_min
)

dir.create(CFG$fig_dir, showWarnings=FALSE, recursive=TRUE)

cat_log <- function(...) cat(..., sep="")
sec <- function(name) {
  cat_log("\n", strrep("-", 70), "\n  ", name, "\n",
          strrep("-", 70), "\n")
}

cat_log("=====================================================================\n")
cat_log("D2 Step 2: QC + visualization\n")
cat_log("Organ: ", CFG$organ, "\n")
cat_log("Start: ", format(Sys.time()), "\n")
cat_log("=====================================================================\n")

# ============ Load D2 step1 ============
sec("Load spot_zone_stage")
zsd_path <- file.path(CFG$d2_dir,
                     paste0("spot_zone_stage_", CFG$organ, ".rds"))
if (!file.exists(zsd_path)) stop("Not found: ", zsd_path,
                                 ". Run D2_step1 first.")
zsd <- readRDS(zsd_path)
cat_log("  loaded ", nrow(zsd), " spots\n")
cat_log("  zone x stage:\n")
print(table(zsd$zone, zsd$stage_label, useNA="ifany"))

# Color scheme for zones (used everywhere)
zone_colors <- c(
  "tumor_core"    = "#d7191c",  # red
  "tumor_margin"  = "#fdae61",  # orange
  "invaded_host"  = "#ffffbf",  # pale yellow
  "distant_host"  = "#abd9e9",  # light blue
  "baseline_host" = "#2c7bb6"   # dark blue
)

zone_levels <- c("tumor_core","tumor_margin","invaded_host",
                 "distant_host","baseline_host")
zsd$zone <- factor(zsd$zone, levels=zone_levels)

stage_levels <- c("baseline","early","mid","late")
zsd$stage_label <- factor(zsd$stage_label, levels=stage_levels)

# ============ FIG 1: Zone spatial map ============
sec("FIG 1: Zone spatial map")
samples_in <- unique(as.character(zsd$sample))
# Order: baseline first, then by stage
sample_order <- c()
for (lbl in stage_levels) {
  s_in_lbl <- unique(as.character(zsd$sample[zsd$stage_label == lbl]))
  sample_order <- c(sample_order, s_in_lbl)
}
sample_order <- intersect(sample_order, samples_in)
cat_log("  sample order: ", paste(sample_order, collapse=", "), "\n")

p_zone_list <- lapply(sample_order, function(s) {
  sub <- zsd[zsd$sample == s, ]
  st <- as.character(unique(sub$stage_label))[1]
  ttl <- paste0(s, " (", st, ", n=", nrow(sub), ")")
  ggplot(sub, aes(x=x, y=y, color=zone)) +
    geom_point(size=0.3, alpha=0.85) +
    scale_color_manual(values=zone_colors, drop=FALSE,
                       name="zone") +
    coord_equal() + theme_void() +
    labs(title=ttl) +
    theme(legend.position="none",
          plot.title=element_text(size=10, hjust=0.5))
})

p_zone <- wrap_plots(p_zone_list, ncol=length(sample_order)) +
  plot_annotation(
    title=paste(CFG$organ, "zone map (D2)"),
    subtitle="red=tumor_core, orange=tumor_margin, yellow=invaded_host, light blue=distant_host, dark blue=baseline_host"
  )

# Add legend strip at bottom (no cowplot needed)
legend_strip <- ggplot(data.frame(x=zone_levels, y=1,
                                   zone=factor(zone_levels, levels=zone_levels)),
                        aes(x=x, y=y, fill=zone)) +
  geom_tile(width=0.9, height=0.5) +
  scale_fill_manual(values=zone_colors, guide="none") +
  geom_text(aes(label=x), size=3.5, vjust=2.5) +
  ylim(0, 2) +
  theme_void() +
  theme(plot.margin=margin(20,20,30,20))

p_zone_final <- p_zone / legend_strip + plot_layout(heights=c(8,1))

out1 <- file.path(CFG$fig_dir, paste0("zone_map_", CFG$organ, ".pdf"))
ggsave(out1, p_zone_final, width=4 + 3.5*length(sample_order), height=6,
       limitsize=FALSE)
cat_log("  saved: ", out1, "\n")

# ============ FIG 2: min_dist_to_tumor histogram ============
sec("FIG 2: min_dist_to_tumor distribution by zone")
host_zsd <- zsd[zsd$zone %in% c("invaded_host","distant_host") &
                !is.na(zsd$min_dist_to_tumor), ]
if (nrow(host_zsd) > 0) {
  cat_log("  ", nrow(host_zsd), " host spots with dist info\n")
  p_dist <- ggplot(host_zsd, aes(x=min_dist_to_tumor, fill=zone)) +
    geom_histogram(binwidth=20, alpha=0.7, position="identity") +
    scale_fill_manual(values=zone_colors[c("invaded_host","distant_host")]) +
    facet_wrap(~ sample, scales="free_y", ncol=3) +
    coord_cartesian(xlim=c(0, 600)) +
    geom_vline(xintercept=200, linetype="dashed", color="black", alpha=0.5) +
    labs(title=paste(CFG$organ, "min distance to tumor (px)"),
         subtitle="dashed line = 200 px threshold (invaded vs distant)",
         x="min distance to tumor neighbor (px)", y="n spots") +
    theme_minimal(base_size=10) +
    theme(legend.position="bottom")
  out2 <- file.path(CFG$fig_dir, paste0("zone_dist_hist_", CFG$organ, ".pdf"))
  ggsave(out2, p_dist, width=10, height=6)
  cat_log("  saved: ", out2, "\n")
} else {
  cat_log("  no host spots with dist info, skipping\n")
}

# ============ FIG 3: zone × stage trend ============
sec("FIG 3: zone × stage trend")
trend_df <- zsd %>%
  filter(stage_label != "baseline") %>%
  count(stage_label, zone) %>%
  group_by(stage_label) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()
# Drop baseline_host (won't appear in non-baseline stages anyway)
trend_df <- trend_df %>% filter(zone != "baseline_host")

p_trend_count <- ggplot(trend_df, aes(x=stage_label, y=n,
                                       color=zone, group=zone)) +
  geom_line(size=1.2) +
  geom_point(size=3) +
  scale_color_manual(values=zone_colors, drop=FALSE) +
  labs(title="Spot counts by zone across stages",
       x="stage", y="n spots", color="zone") +
  theme_minimal(base_size=11) +
  theme(legend.position="right")

p_trend_prop <- ggplot(trend_df, aes(x=stage_label, y=prop,
                                      fill=zone)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=zone_colors, drop=FALSE) +
  scale_y_continuous(labels=percent_format()) +
  labs(title="Zone composition by stage",
       x="stage", y="proportion", fill="zone") +
  theme_minimal(base_size=11) +
  theme(legend.position="right")

p_trend <- p_trend_count / p_trend_prop
out3 <- file.path(CFG$fig_dir, paste0("zone_stage_trend_", CFG$organ, ".pdf"))
ggsave(out3, p_trend, width=8, height=8)
cat_log("  saved: ", out3, "\n")
cat_log("\n  trend table:\n")
print(as.data.frame(trend_df))

# ============ Celltype composition per cell ============
sec("Celltype composition per cell (RCTD)")
samples <- unique(as.character(zsd$sample))
weights_list <- list()
for (s in samples) {
  fp <- file.path(CFG$rctd_per_sample_dir, paste0(s, "_rctd_results.rds"))
  if (!file.exists(fp)) {
    cat_log("  skip ", s, " (no RCTD)\n"); next
  }
  r <- readRDS(fp)
  if (!("weights" %in% names(r@results))) next
  w <- r@results$weights
  rs <- Matrix::rowSums(w); rs[rs == 0] <- 1
  w <- w / rs
  rownames(w) <- paste0(s, "_", rownames(w))
  weights_list[[s]] <- w
  cat_log("  ", s, ": ", nrow(w), " x ", ncol(w), "\n")
  rm(r); gc(verbose=FALSE)
}

if (length(weights_list) == 0) {
  cat_log("  no RCTD weights -- skipping celltype analysis\n")
} else {
  all_cts <- sort(unique(unlist(lapply(weights_list, colnames))))
  aligned <- lapply(weights_list, function(w) {
    missing <- setdiff(all_cts, colnames(w))
    if (length(missing) > 0) {
      zeros <- Matrix(0, nrow=nrow(w), ncol=length(missing), sparse=TRUE,
                      dimnames=list(rownames(w), missing))
      w <- cbind(w, zeros)
    }
    w[, all_cts, drop=FALSE]
  })
  weights <- do.call(rbind, aligned)
  rm(weights_list, aligned); gc(verbose=FALSE)

  zsd$cell_id <- paste0(as.character(zsd$zone), "__",
                        as.character(zsd$stage_label))
  common <- intersect(zsd$spot_key, rownames(weights))
  cat_log("\n  RCTD intersect: ", length(common), " (",
          round(100*length(common)/nrow(zsd),1), "%)\n")
  if (length(common) < 100) {
    cat_log("  too few intersects -- skipping celltype heatmap\n")
  } else {
    w_sub <- as.matrix(weights[common, , drop=FALSE])
    z_sub <- zsd[match(common, zsd$spot_key), ]

    df <- as.data.frame(w_sub)
    df$cell_id <- z_sub$cell_id

    cell_comp <- df %>%
      group_by(cell_id) %>%
      summarise(n_rctd_spots = n(),
                across(everything(), mean), .groups="drop")
    cat_log("\n  cells with celltype data:\n")
    cell_comp %>% select(cell_id, n_rctd_spots) %>%
      arrange(desc(n_rctd_spots)) %>% print(n=Inf)

    # filter cells with too few RCTD spots
    cell_comp_use <- cell_comp %>%
      filter(n_rctd_spots >= CFG$rctd_intersect_min)

    out_csv <- file.path(CFG$d2_dir,
                         paste0("cell_celltype_composition_",
                                CFG$organ, ".csv"))
    write_csv(cell_comp, out_csv)
    cat_log("  saved: ", out_csv, "\n")

    # cell summary
    cell_summary <- zsd %>%
      group_by(cell_id) %>%
      summarise(n_spots = n(),
                mean_human_pct = mean(human_pct, na.rm=TRUE),
                med_human_pct  = median(human_pct, na.rm=TRUE),
                .groups="drop") %>%
      left_join(cell_comp %>% select(cell_id, n_rctd_spots), by="cell_id") %>%
      arrange(desc(n_spots))
    out_summary <- file.path(CFG$d2_dir,
                            paste0("cell_summary_", CFG$organ, ".csv"))
    write_csv(cell_summary, out_summary)
    cat_log("  saved: ", out_summary, "\n")

    # ===== FIG 4: heatmap (cell × top celltype) =====
    sec("FIG 4: Celltype heatmap (cell × top celltypes)")
    # Pick top N celltypes by total weight across cells
    if (nrow(cell_comp_use) > 0) {
      mean_weights <- colMeans(cell_comp_use %>%
                                select(-cell_id, -n_rctd_spots))
      top_ct <- names(sort(mean_weights, decreasing=TRUE))[
        seq_len(min(CFG$top_n_celltypes, length(mean_weights)))]
      cat_log("  top ", CFG$top_n_celltypes, " celltypes: ",
              paste(top_ct, collapse=", "), "\n")

      hm_df <- cell_comp_use %>%
        select(cell_id, all_of(top_ct)) %>%
        pivot_longer(-cell_id, names_to="celltype", values_to="proportion")

      # Order cells by zone then stage
      cell_order <- expand.grid(stage_label=stage_levels,
                                zone=zone_levels) %>%
        mutate(cell_id = paste0(zone, "__", stage_label)) %>%
        pull(cell_id)
      cell_order <- intersect(cell_order, unique(hm_df$cell_id))
      hm_df$cell_id <- factor(hm_df$cell_id, levels=cell_order)
      hm_df$celltype <- factor(hm_df$celltype, levels=top_ct)

      p_hm <- ggplot(hm_df, aes(x=cell_id, y=celltype, fill=proportion)) +
        geom_tile() +
        scale_fill_gradientn(
          colors=c("white","#fdae61","#d7191c"),
          limits=c(0, max(0.5, max(hm_df$proportion, na.rm=TRUE)))) +
        labs(title=paste(CFG$organ, "celltype composition per cell (zone × stage)"),
             x="cell (zone__stage)", y="celltype",
             fill="mean weight") +
        theme_minimal(base_size=10) +
        theme(axis.text.x=element_text(angle=45, hjust=1, size=8),
              panel.grid=element_blank())
      out4 <- file.path(CFG$fig_dir,
                       paste0("cell_celltype_heatmap_", CFG$organ, ".pdf"))
      ggsave(out4, p_hm, width=12, height=6)
      cat_log("  saved: ", out4, "\n")
    } else {
      cat_log("  no cells pass intersect threshold\n")
    }
  }
}

# ============ Summary ============
sec("Summary")
cat_log("\nFinal cell sizes (zone × stage):\n")
final_pivot <- zsd %>%
  count(zone, stage_label) %>%
  pivot_wider(names_from=stage_label, values_from=n, values_fill=0)
print(as.data.frame(final_pivot))

cat_log("\n=== DONE: ", format(Sys.time()), " ===\n")
cat_log("Figures: ", CFG$fig_dir, "\n\n")
