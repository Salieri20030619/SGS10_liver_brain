# =============================================================
# P1_upgraded_plots.R
# Nat Comm / Cell Rep 级美学的 3 张升级图
#
# Figures:
#   P1b_v2_liver_dumbbell.png     哑铃图 (替换 boxplot)
#   P1b_v2_brain_dumbbell.png
#   P1c_v2_liver_lollipop.png     棒棒糖 + callout (替换 FC bar)
#   P1c_v2_brain_lollipop.png
#   P1e_liver3_spatial_density.png  空间密度图 (全新)
# =============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggtext)      # for markdown-enabled labels (fallback to plain)
  library(scales)
  library(png)
  library(grid)
})

T4_DIR   <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis/SGS10_deliverables_2026-04-20/08_spatial_Task_B/T4_spatial_domain"
P1_DIR   <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis/P1_tumor_altered_niche"
RCTD_DIR <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis/spatial_results/rctd/per_sample"
P0_DIR   <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis/liver_species_split_v2"

# HE image paths (from previous v4 HE overlay work)
HE_PATHS <- list(
  Liver_3 = "/mnt/BioData/bio/SGS10例/SGS25090559-8例样本数据/data/matrix/25090559_Liver_3_expression/Liver_3_aligned_HE.png"
)

OUT_DIR <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis/P1_upgraded_plots"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------
# Nat Comm / Cell Rep 色板 (Lancet / NEJM palette inspired)
# -------------------------------------------------------------
COLORS <- list(
  up         = "#E64B35",   # Lancet red
  down       = "#4DBBD5",   # Lancet blue
  neutral    = "#B4B4B4",
  accent     = "#F39B7F",   # coral accent
  tumor      = "#DC0000",   # Cell Rep red
  baseline   = "#3C5488",   # Cell Rep blue
  # domain palette (10 colors, premium quality)
  d05        = "#FDAE61",   # warm orange (Pericytes tumor-core)
  d09        = "#D73027",   # deep red (host-induced)
  d10        = "#1A9850",   # green (normal)
  d01        = "#9E9E9E",
  d02        = "#9E9E9E",
  d03        = "#9E9E9E",
  d04        = "#9E9E9E",
  d06        = "#9E9E9E",
  d07        = "#9E9E9E",
  d08        = "#9E9E9E"
)

# -------------------------------------------------------------
# Base theme (Nat Comm-ish: clean, sparse grid, Arial-like sans)
# -------------------------------------------------------------
nc_theme <- function(base_size = 11) {
  theme_classic(base_size = base_size) %+replace%
  theme(
    text             = element_text(family = "sans", color = "grey15"),
    plot.title       = element_text(size = base_size + 2, face = "bold",
                                     color = "grey5", hjust = 0,
                                     margin = margin(b = 4)),
    plot.subtitle    = element_text(size = base_size - 1, color = "grey40",
                                     hjust = 0, margin = margin(b = 10)),
    plot.caption     = element_text(size = base_size - 2, color = "grey50",
                                     hjust = 1, face = "italic"),
    axis.title       = element_text(size = base_size, color = "grey25"),
    axis.text        = element_text(size = base_size - 1, color = "grey30"),
    axis.line        = element_line(color = "grey40", linewidth = 0.4),
    axis.ticks       = element_line(color = "grey40", linewidth = 0.4),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title     = element_text(size = base_size - 1, color = "grey25", face = "bold"),
    legend.text      = element_text(size = base_size - 2, color = "grey30"),
    legend.key.size  = unit(0.9, "lines"),
    legend.background = element_blank(),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin      = margin(15, 18, 12, 15)
  )
}

log_msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

# =============================================================
# Load P1 stats (already computed)
# =============================================================
load_p1_data <- function(organ) {
  per_sample <- read.csv(file.path(T4_DIR, paste0(organ, "_per_sample_domain.csv")),
                         stringsAsFactors = FALSE)
  naming     <- read.csv(file.path(T4_DIR, paste0(organ, "_domain_naming.csv")),
                         stringsAsFactors = FALSE)
  summary    <- read.csv(file.path(P1_DIR, paste0("P1_summary_", organ, ".csv")),
                         stringsAsFactors = FALSE)
  list(per_sample = per_sample, naming = naming, summary = summary)
}

# =============================================================
# P1b v2: Dumbbell plot
# =============================================================
plot_dumbbell <- function(data, organ_label) {
  # Prepare: per domain, compute mean & SD for each group
  per_sample <- data$per_sample
  summary    <- data$summary
  
  # join sample variance
  stats <- per_sample %>%
    group_by(domain, group) %>%
    summarise(mean_pct = mean(pct), sd_pct = sd(pct),
              n = n(), .groups = "drop") %>%
    left_join(data$naming[, c("domain", "auto_name")], by = "domain")
  
  # Pivot wider for line drawing
  wide <- stats %>%
    select(domain, group, mean_pct, sd_pct, auto_name) %>%
    pivot_wider(names_from = group, values_from = c(mean_pct, sd_pct))
  
  wide$delta     <- wide$mean_pct_tumor - wide$mean_pct_baseline
  wide$direction <- ifelse(wide$delta > 0, "Up in tumor", "Down in tumor")
  wide <- wide[order(wide$delta, decreasing = TRUE), ]
  wide$auto_name <- factor(wide$auto_name, levels = rev(wide$auto_name))
  
  # Long format for points
  long <- wide %>%
    select(auto_name, mean_pct_baseline, mean_pct_tumor,
           sd_pct_baseline, sd_pct_tumor, direction) %>%
    pivot_longer(cols = c(mean_pct_baseline, mean_pct_tumor),
                 names_to = "group", values_to = "mean_pct") %>%
    mutate(
      group = sub("mean_pct_", "", group),
      sd = ifelse(group == "tumor", sd_pct_tumor, sd_pct_baseline)
    )
  
  # Background shading rows for visual rhythm
  n_dom <- nrow(wide)
  bg <- data.frame(
    y = seq_len(n_dom),
    color = rep(c("#FAFAFA", "white"), length.out = n_dom)
  )
  
  p <- ggplot() +
    # subtle alternating background
    geom_rect(data = bg,
              aes(ymin = y - 0.5, ymax = y + 0.5, fill = color),
              xmin = -Inf, xmax = Inf, color = NA, show.legend = FALSE) +
    scale_fill_identity() +
    # zero reference
    geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3, linetype = "dashed") +
    # connecting segments (thick, colored by direction)
    geom_segment(data = wide,
                 aes(x = mean_pct_baseline, xend = mean_pct_tumor,
                     y = auto_name, yend = auto_name,
                     color = direction),
                 linewidth = 1.6, lineend = "round", alpha = 0.9) +
    # SD whiskers on tumor point (per sample variance)
    geom_segment(data = wide,
                 aes(x = mean_pct_tumor - sd_pct_tumor,
                     xend = mean_pct_tumor + sd_pct_tumor,
                     y = auto_name, yend = auto_name),
                 color = COLORS$tumor, linewidth = 0.5, alpha = 0.5) +
    # baseline point (smaller, hollow)
    geom_point(data = long[long$group == "baseline", ],
               aes(x = mean_pct, y = auto_name),
               color = COLORS$baseline, fill = "white",
               shape = 21, stroke = 1.4, size = 3.5) +
    # tumor point (larger, filled)
    geom_point(data = long[long$group == "tumor", ],
               aes(x = mean_pct, y = auto_name),
               color = "white", fill = COLORS$tumor,
               shape = 21, stroke = 0.8, size = 4.8) +
    # delta label at far right
    geom_text(data = wide,
              aes(x = pmax(mean_pct_baseline, mean_pct_tumor) + 1.5,
                  y = auto_name,
                  label = sprintf("%+.1f%%", delta),
                  color = direction),
              hjust = 0, size = 3.3, fontface = "bold") +
    scale_color_manual(values = c("Up in tumor" = COLORS$up,
                                   "Down in tumor" = COLORS$down),
                       name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.18)),
                       labels = function(x) paste0(x, "%")) +
    labs(
      title    = sprintf("Domain proportion shift — %s", organ_label),
      subtitle = sprintf(
        "Each line connects mean baseline (○) to tumor (●) across %d niches. Bar = tumor ± SD across 3 tumor samples.",
        n_dom
      ),
      x = "Mean proportion of spots",
      y = NULL,
      caption = "Baseline: Blank + Control (n=2)   •   Tumor: n=3"
    ) +
    nc_theme() +
    theme(
      legend.position = "top",
      legend.justification = "left",
      legend.margin = margin(0, 0, 4, 0),
      legend.key = element_rect(fill = NA),
      axis.text.y = element_text(hjust = 1, color = "grey20")
    ) +
    guides(color = guide_legend(override.aes = list(linewidth = 2)))
  
  p
}

# =============================================================
# P1c v2: Lollipop with callout
# =============================================================
plot_lollipop <- function(data, organ_label) {
  df <- data$summary
  df <- df[order(df$delta, decreasing = TRUE), ]
  df$auto_name <- factor(df$auto_name, levels = rev(df$auto_name))
  df$direction <- ifelse(df$delta > 0, "Up in tumor", "Down in tumor")
  df$neg_log_p <- -log10(pmax(df$wilcox_p, 1e-10))
  
  # Find the top tumor-induced niche for callout
  top_up <- df[which.max(df$delta), ]
  
  # 用于标签位置的偏移量
  x_range <- max(df$delta) - min(df$delta)
  
  # 标签放在棒棒糖圆点的"外侧"，避免遮挡
  df$label_x   <- df$delta + sign(df$delta) * (x_range * 0.045) + 
                  ifelse(df$delta == 0, x_range * 0.03, 0)
  df$label_hjust <- ifelse(df$delta > 0, 0, 1)
  
  # 副标题缩短（避免被截断）
  sub_text <- sprintf(
    "Δ proportion by spatial domain.  Strongest shift: %s (+%.1f%%, fc=%.2f×).",
    sub("^(D\\d+) \\[.*", "\\1", as.character(top_up$auto_name)),
    top_up$delta, top_up$fc
  )
  
  p <- ggplot(df, aes(x = delta, y = auto_name)) +
    # zero axis
    geom_vline(xintercept = 0, color = "grey45", linewidth = 0.5) +
    # |delta| = 3 threshold lines
    geom_vline(xintercept = c(-3, 3), color = "grey80",
               linewidth = 0.3, linetype = "dotted") +
    # lollipop stems
    geom_segment(aes(x = 0, xend = delta,
                     y = auto_name, yend = auto_name,
                     color = direction),
                 linewidth = 1.1, alpha = 0.85) +
    # dots sized by -log10(p)
    geom_point(aes(color = direction, size = neg_log_p)) +
    # value labels (放在点外侧)
    geom_text(aes(label = sprintf("%+.1f%%", delta),
                  color = direction,
                  x = label_x,
                  hjust = label_hjust),
              size = 3.0, fontface = "bold") +
    scale_color_manual(values = c("Up in tumor" = COLORS$up,
                                  "Down in tumor" = COLORS$down),
                       name = NULL) +
    scale_size_continuous(
      range  = c(3, 7),
      name   = expression(-log[10]*"(p)"),
      breaks = c(0.5, 1.0, 1.5),
      limits = c(0, 2)
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.16, 0.22)),
                       labels = function(x) paste0(x, "%")) +
    labs(
      title    = sprintf("Tumor vs Baseline — %s", organ_label),
      subtitle = sub_text,
      x = expression("Δ proportion (tumor − baseline)"),
      y = NULL,
      caption = "Dot size ∝ −log₁₀(Wilcoxon p).   Dotted lines = ±3% threshold.   n=3 tumor, n=2 baseline."
    ) +
    nc_theme() +
    theme(
      legend.position = "top",
      legend.box = "horizontal",
      legend.justification = "left",
      legend.margin = margin(0, 0, 4, 0),
      legend.key = element_rect(fill = NA),
      axis.text.y = element_text(hjust = 1, color = "grey20", size = 10),
      plot.subtitle = element_text(size = 10)
    ) +
    guides(
      color = guide_legend(order = 1,
                           override.aes = list(linewidth = 0, size = 4)),
      size  = guide_legend(order = 2,
                           override.aes = list(color = "grey40"))
    )
  
  p
}

# =============================================================
# P1e: Spatial density overlay on HE
# =============================================================
assign_spot_to_domain <- function(weights, profile_df) {
  prof <- as.matrix(profile_df[, -1, drop = FALSE])
  rownames(prof) <- profile_df[[1]]
  common_ct <- intersect(colnames(weights), colnames(prof))
  w <- as.matrix(weights[, common_ct, drop = FALSE])
  p <- prof[, common_ct, drop = FALSE]
  w_norm <- w / pmax(sqrt(rowSums(w^2)), 1e-9)
  p_norm <- p / pmax(sqrt(rowSums(p^2)), 1e-9)
  sims <- w_norm %*% t(p_norm)
  rownames(prof)[max.col(sims)]
}

build_spatial_domain_df <- function(sample) {
  rctd_path <- file.path(RCTD_DIR, paste0(sample, "_rctd_results.rds"))
  rctd <- readRDS(rctd_path)
  coords  <- rctd@spatialRNA@coords
  weights <- rctd@results$weights
  
  profile_df <- read.csv(file.path(T4_DIR, "liver_domain_ct_profile.csv"),
                         stringsAsFactors = FALSE, check.names = FALSE)
  colnames(profile_df)[1] <- "domain"
  
  dom <- assign_spot_to_domain(weights, profile_df)
  df <- data.frame(
    barcode = rownames(weights),
    x = coords[rownames(weights), "x"],
    y = coords[rownames(weights), "y"],
    domain = dom,
    stringsAsFactors = FALSE
  )
  # also get size
  attr(df, "size_x") <- max(coords$x) * 1.05
  attr(df, "size_y") <- max(coords$y) * 1.05
  df
}

plot_spatial_density <- function(sample) {
  log_msg("Building spatial domain df for ", sample, "...")
  sp <- build_spatial_domain_df(sample)
  
  # Focus domains: D05 (tumor-core), D09 (host-induced), D10 (normal)
  focus_doms <- c("D05", "D09", "D10")
  focus_labels <- c(
    D05 = "D05 [Pericytes]  — tumor-core",
    D09 = "D09 [LSEC_Mid+M1]  — host-induced",
    D10 = "D10 [Hep_Periportal]  — normal"
  )
  focus_colors <- c(
    D05 = COLORS$d05,
    D09 = COLORS$d09,
    D10 = COLORS$d10
  )
  
  sp_focus <- sp[sp$domain %in% focus_doms, ]
  sp_focus$domain_label <- focus_labels[sp_focus$domain]
  sp_focus$domain_label <- factor(sp_focus$domain_label, levels = focus_labels)
  
  # Read HE image
  he_path <- HE_PATHS[[sample]]
  he_img <- NULL
  size_x <- attr(sp, "size_x"); size_y <- attr(sp, "size_y")
  
  if (!is.null(he_path) && file.exists(he_path)) {
    log_msg("  Reading HE: ", he_path)
    he_img <- png::readPNG(he_path)
    he_grob <- grid::rasterGrob(he_img, interpolate = TRUE)
  }
  
  # Aspect
  plot_width  <- 15
  plot_height <- plot_width * size_y / size_x * 2.2  # 2 rows
  
  # Build faceted plot
  # Left: HE + contour overlay; Right: pure density heatmap
  
  # we'll use facet_grid with version + domain
  # Build two dataframes
  df_left <- sp_focus; df_left$version <- "Spatial context"
  df_right <- sp_focus; df_right$version <- "Density (isolated)"
  df_all <- rbind(df_left, df_right)
  df_all$version <- factor(df_all$version, levels = c("Spatial context", "Density (isolated)"))
  
  p <- ggplot(df_all, aes(x = x, y = y))
  
  # Add HE only for "Spatial context" facet
  if (!is.null(he_img)) {
    # Trick: add HE rasters only to one facet via annotation_custom on a dummy layer
    # Simpler: add it globally then alpha via facet (not elegant)
    # Cleanest: use patchwork to compose separately. Do that instead.
    NULL
  }
  
  # Use patchwork for real composition
  library(patchwork)
  
  # Left column: 3 rows (one per domain) with HE + 填充等高线
  build_left_panel <- function(dom_code) {
    df_d <- sp_focus[sp_focus$domain == dom_code, ]
    df_other <- sp[!sp$domain %in% dom_code, ]
    color_main <- focus_colors[dom_code]
    
    g <- ggplot()
    if (!is.null(he_img)) {
      g <- g + annotation_custom(he_grob,
                                  xmin = 0, xmax = size_x,
                                  ymin = size_y, ymax = 0)
    }
    g <- g +
      # 其他 domain 的点 (浅灰背景)
      geom_point(data = df_other,
                 aes(x = x, y = y),
                 color = "grey80", size = 0.12, alpha = 0.15) +
      # 目标 domain 点 (较小,让位给等高线)
      geom_point(data = df_d,
                 aes(x = x, y = y),
                 color = color_main,
                 size = 0.35, alpha = 0.55) +
      # 填充等高线 (关键改进: 用 polygon + alpha 映射)
      stat_density_2d(data = df_d,
                      aes(x = x, y = y,
                          fill = after_stat(level),
                          alpha = after_stat(level)),
                      geom = "polygon",
                      bins = 6,
                      contour = TRUE,
                      show.legend = FALSE) +
      # 最核心 2 层等高线 (深色勾边, 简洁)
      stat_density_2d(data = df_d,
                      aes(x = x, y = y),
                      color = color_main,
                      linewidth = 0.5, bins = 3,
                      contour = TRUE) +
      scale_fill_gradient(low = color_main, high = color_main) +
      scale_alpha(range = c(0, 0.45)) +
      coord_fixed(xlim = c(0, size_x), ylim = c(size_y, 0), expand = FALSE) +
      labs(subtitle = focus_labels[dom_code]) +
      nc_theme(base_size = 10) +
      theme(
        axis.line   = element_blank(),
        axis.text   = element_blank(),
        axis.ticks  = element_blank(),
        axis.title  = element_blank(),
        plot.subtitle = element_text(face = "bold",
                                      color = color_main,
                                      size = 11, margin = margin(b = 2, t = 4)),
        plot.margin = margin(4, 6, 2, 6)
      )
    g
  }
  
  build_right_panel <- function(dom_code) {
    df_d <- sp_focus[sp_focus$domain == dom_code, ]
    color_main <- focus_colors[dom_code]
    
    g <- ggplot(df_d, aes(x = x, y = y)) +
      stat_density_2d(aes(fill = after_stat(density)),
                      geom = "raster", contour = FALSE,
                      n = 220) +
      scale_fill_gradient(low = "white", high = color_main,
                          guide = "none") +
      # 减少等高线 (4 条即可,避免拥挤)
      stat_density_2d(color = "grey25", linewidth = 0.25, bins = 4,
                      alpha = 0.5) +
      coord_fixed(xlim = c(0, size_x), ylim = c(size_y, 0), expand = FALSE) +
      labs(subtitle = paste0(dom_code, " isolated density")) +
      nc_theme(base_size = 10) +
      theme(
        axis.line   = element_blank(),
        axis.text   = element_blank(),
        axis.ticks  = element_blank(),
        axis.title  = element_blank(),
        plot.subtitle = element_text(color = "grey30", size = 10,
                                      margin = margin(b = 2, t = 4)),
        plot.margin = margin(4, 6, 2, 6),
        panel.background = element_rect(fill = "white", color = NA)
      )
    g
  }
  
  # Compose: 3 rows × 2 columns
  rows <- lapply(focus_doms, function(d) {
    build_left_panel(d) | build_right_panel(d)
  })
  
  combined <- Reduce(`/`, rows)
  combined <- combined +
    plot_annotation(
      title    = sprintf("Spatial organization of key domains — %s", sample),
      subtitle = "Left: HE with domain spots and 2D density contours. Right: isolated density for structure clarity.",
      theme = nc_theme() + theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11, color = "grey40",
                                      margin = margin(b = 8))
      )
    )
  
  list(plot = combined, height = plot_height)
}

# =============================================================
# MAIN
# =============================================================

log_msg("========== P1 Upgraded Plots ==========")

# --- P1b v2: dumbbell (both organs) ---
for (organ in c("liver", "brain")) {
  log_msg("Building dumbbell for ", organ, "...")
  data <- load_p1_data(organ)
  p <- plot_dumbbell(data, tools::toTitleCase(organ))
  out <- file.path(OUT_DIR, sprintf("P1b_v2_%s_dumbbell.png", organ))
  ggsave(out, p, width = 9, height = 5.5, dpi = 200, bg = "white")
  log_msg("  Saved: ", out)
}

# --- P1c v2: lollipop ---
for (organ in c("liver", "brain")) {
  log_msg("Building lollipop for ", organ, "...")
  data <- load_p1_data(organ)
  p <- plot_lollipop(data, tools::toTitleCase(organ))
  out <- file.path(OUT_DIR, sprintf("P1c_v2_%s_lollipop.png", organ))
  ggsave(out, p, width = 9, height = 5.2, dpi = 200, bg = "white")
  log_msg("  Saved: ", out)
}

# --- P1e: spatial density ---
log_msg("Building spatial density for Liver_3...")
sd_res <- plot_spatial_density("Liver_3")
out <- file.path(OUT_DIR, "P1e_Liver_3_spatial_density.png")
ggsave(out, sd_res$plot,
       width = 14, height = sd_res$height,
       dpi = 180, bg = "white", limitsize = FALSE)
log_msg("  Saved: ", out)

log_msg("========== DONE ==========")
for (f in list.files(OUT_DIR)) {
  info <- file.info(file.path(OUT_DIR, f))
  log_msg(sprintf("  %-50s %6.1f KB", f, info$size / 1024))
}
