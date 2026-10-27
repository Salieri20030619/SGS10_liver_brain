# =============================================================
# P1_tumor_altered_niche.R
# Purpose: 识别 tumor-altered spatial niche
# 
# Inputs: T4 已有的 per_sample_domain + tumor_vs_baseline + naming
#         P0 的 spot_meta.csv (只肝有 human_frac)
# Outputs: 
#   P1a_*_domain_proportion_three_groups.png   三组独立展示 (附图)
#   P1b_*_domain_proportion_tumor_vs_baseline.png  主图
#   P1c_*_domain_foldchange.png  FC 条形图
#   P1d_liver_domain_human_frac.png  肝端按 domain 的 human_frac boxplot
#   P1_summary_liver.csv, P1_summary_brain.csv  统计表
#   P1_significant_niches.csv  最终挑出的 altered niches
# =============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
})

# -------------------------------------------------------------
# Config
# -------------------------------------------------------------
T4_DIR  <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis/SGS10_deliverables_2026-04-20/08_spatial_Task_B/T4_spatial_domain"
P0_DIR  <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis/liver_species_split_v2"
RCTD_DIR <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis/spatial_results/rctd/per_sample"
OUT_DIR <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis/P1_tumor_altered_niche"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Sample group mapping (for 3-group split in supplementary figure)
SAMPLE_GROUP <- list(
  liver = list(
    tumor    = c("gan56", "Liver_3", "Liver_8"),
    Bla      = c("Bla_liver_1"),
    Con      = c("Con_liver_1")
  ),
  brain = list(
    tumor    = c("nao56", "Brain_3", "Brain_8"),
    Bla      = c("Bla_brain_1"),
    Con      = c("Con_brain_1")
  )
)

# Domain color palette (reuse HE overlay v4's)
DOMAIN_COLORS <- c(
  D01 = "#e41a1c", D02 = "#377eb8", D03 = "#4daf4a", D04 = "#984ea3",
  D05 = "#ff7f00", D06 = "#ffff33", D07 = "#a65628", D08 = "#f781bf",
  D09 = "#66c2a5", D10 = "#b3b3b3"
)

GROUP_COLORS <- c(tumor = "#d73027", baseline = "#4575b4",
                  Bla = "#74add1", Con = "#abd9e9")

# -------------------------------------------------------------
# Base theme (consistent with HE overlay v4 style)
# -------------------------------------------------------------
base_theme <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey40", size = 9),
    legend.position = "right"
  )

log_msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

# -------------------------------------------------------------
# STEP 1: Load T4 data for one organ
# -------------------------------------------------------------
load_organ_data <- function(organ) {
  list(
    per_sample = read.csv(file.path(T4_DIR, paste0(organ, "_per_sample_domain.csv")),
                          stringsAsFactors = FALSE),
    tvb        = read.csv(file.path(T4_DIR, paste0(organ, "_tumor_vs_baseline.csv")),
                          stringsAsFactors = FALSE),
    naming     = read.csv(file.path(T4_DIR, paste0(organ, "_domain_naming.csv")),
                          stringsAsFactors = FALSE)
  )
}

# Add sub-group (tumor / Bla / Con) to per_sample based on SAMPLE_GROUP
annotate_subgroup <- function(per_sample, organ) {
  sg <- SAMPLE_GROUP[[organ]]
  per_sample$subgroup <- ifelse(per_sample$sample %in% sg$tumor, "tumor",
                         ifelse(per_sample$sample %in% sg$Bla,   "Bla",
                         ifelse(per_sample$sample %in% sg$Con,   "Con", NA)))
  per_sample
}

# -------------------------------------------------------------
# STEP 2: Wilcoxon test per domain (tumor vs combined baseline)
# -------------------------------------------------------------
compute_wilcox <- function(per_sample) {
  per_sample %>%
    group_by(domain) %>%
    summarise(
      n_tumor    = sum(group == "tumor"),
      n_baseline = sum(group == "baseline"),
      mean_tumor = mean(pct[group == "tumor"]),
      mean_base  = mean(pct[group == "baseline"]),
      delta      = mean_tumor - mean_base,
      fc         = mean_tumor / pmax(mean_base, 0.01),
      log2fc     = log2(pmax(fc, 0.01)),
      wilcox_p   = tryCatch(
        wilcox.test(pct[group == "tumor"], pct[group == "baseline"],
                    exact = FALSE)$p.value,
        error = function(e) NA_real_
      ),
      .groups = "drop"
    ) %>%
    arrange(desc(abs(delta)))
}

# -------------------------------------------------------------
# STEP 3: Plots
# -------------------------------------------------------------

# ---- P1a: stacked bar per sample (3 groups) ----
plot_per_sample_3groups <- function(per_sample, organ, naming) {
  # order samples: tumor first, then Bla, then Con
  sg <- SAMPLE_GROUP[[organ]]
  sample_order <- c(sg$tumor, sg$Bla, sg$Con)
  ps <- per_sample
  ps$sample <- factor(ps$sample, levels = sample_order)
  
  # Add auto_name for domain labels
  ps$domain_label <- naming$auto_name[match(ps$domain, naming$domain)]
  dom_order <- naming$domain
  ps$domain <- factor(ps$domain, levels = dom_order)
  
  p <- ggplot(ps, aes(x = sample, y = pct, fill = domain)) +
    geom_bar(stat = "identity", width = 0.75) +
    scale_fill_manual(values = DOMAIN_COLORS,
                      labels = naming$auto_name,
                      name = "Spatial domain") +
    facet_grid(~ subgroup, scales = "free_x", space = "free_x") +
    labs(title = sprintf("Spatial domain composition — %s (per sample)",
                         tools::toTitleCase(organ)),
         subtitle = "Tumor (n=3) vs Blank (n=1) vs Control (n=1)",
         x = NULL, y = "Proportion (%)") +
    base_theme +
    theme(strip.background = element_rect(fill = "grey90", color = NA),
          strip.text = element_text(face = "bold"))
  p
}

# ---- P1b: tumor vs baseline comparison (main figure) ----
plot_tumor_vs_baseline_summary <- function(per_sample, tvb, organ, naming) {
  ps <- per_sample
  ps$domain_label <- naming$auto_name[match(ps$domain, naming$domain)]
  dom_order <- naming$domain
  ps$domain <- factor(ps$domain, levels = dom_order)
  
  # boxplot (each sample as a point), grouped by tumor/baseline
  p <- ggplot(ps, aes(x = domain, y = pct, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6, position = position_dodge(0.7), width = 0.6) +
    geom_jitter(aes(color = group),
                position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.7),
                size = 2, alpha = 0.9) +
    scale_fill_manual(values = c(tumor = "#d73027", baseline = "#4575b4")) +
    scale_color_manual(values = c(tumor = "#800000", baseline = "#1a3d66")) +
    scale_x_discrete(labels = setNames(naming$auto_name, naming$domain)) +
    labs(title = sprintf("Per-sample domain proportion — %s",
                         tools::toTitleCase(organ)),
         subtitle = "Each point = one sample. Tumor (n=3) vs Baseline (Bla+Con, n=2)",
         x = NULL, y = "Proportion (%)",
         fill = "Group", color = "Group") +
    base_theme
  p
}

# ---- P1c: fold change bar chart ----
plot_foldchange <- function(stats_df, organ, naming) {
  df <- stats_df
  df$auto_name <- naming$auto_name[match(df$domain, naming$domain)]
  df <- df[order(df$delta, decreasing = TRUE), ]
  df$auto_name <- factor(df$auto_name, levels = df$auto_name)
  df$direction <- ifelse(df$delta > 0, "Up in tumor", "Down in tumor")
  
  # annotate p-value when significant-ish
  df$p_label <- ifelse(df$wilcox_p < 0.05, "*",
                ifelse(df$wilcox_p < 0.1, ".", ""))
  
  p <- ggplot(df, aes(x = auto_name, y = delta, fill = direction)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = 0, color = "grey40") +
    geom_text(aes(label = p_label,
                  y = delta + sign(delta) * 0.3),
              size = 5, hjust = 0.5) +
    scale_fill_manual(values = c("Up in tumor" = "#d73027",
                                 "Down in tumor" = "#4575b4")) +
    coord_flip() +
    labs(title = sprintf("Tumor vs Baseline domain shift — %s",
                         tools::toTitleCase(organ)),
         subtitle = "Delta = mean(tumor%) − mean(baseline%). * p<0.05, . p<0.1 (Wilcoxon)",
         x = NULL, y = "Δ proportion (%)",
         fill = "Direction") +
    base_theme +
    theme(axis.text.y = element_text(angle = 0, hjust = 1))
  p
}

# -------------------------------------------------------------
# STEP 4: Liver-specific: per-domain human_frac
# -------------------------------------------------------------
# Need per-spot domain assignment — recompute via cosine similarity
# (same method HE overlay v4 used)

load_rctd_coords_and_weights <- function(sample) {
  fp <- file.path(RCTD_DIR, paste0(sample, "_rctd_results.rds"))
  if (!file.exists(fp)) return(NULL)
  rctd <- readRDS(fp)
  list(
    coords  = rctd@spatialRNA@coords,
    weights = rctd@results$weights
  )
}

# Assign each spot to nearest domain by cosine similarity with profile
assign_spot_to_domain <- function(weights, profile_df) {
  prof <- as.matrix(profile_df[, -1, drop = FALSE])
  rownames(prof) <- profile_df[[1]]
  
  # align CT columns
  common_ct <- intersect(colnames(weights), colnames(prof))
  w <- as.matrix(weights[, common_ct, drop = FALSE])
  p <- prof[, common_ct, drop = FALSE]
  
  # row-normalize both
  w_norm <- w / pmax(sqrt(rowSums(w^2)), 1e-9)
  p_norm <- p / pmax(sqrt(rowSums(p^2)), 1e-9)
  
  sims <- w_norm %*% t(p_norm)
  max_idx <- max.col(sims)
  rownames(prof)[max_idx]
}

compute_human_frac_per_domain <- function(organ) {
  if (organ != "liver") return(NULL)
  
  profile_df <- read.csv(file.path(T4_DIR, "liver_domain_ct_profile.csv"),
                         stringsAsFactors = FALSE, check.names = FALSE)
  # first column is domain (unnamed originally)
  colnames(profile_df)[1] <- "domain"
  
  results <- list()
  for (sample in c("gan56", "Liver_3", "Liver_8")) {
    log_msg("  Computing human_frac by domain for ", sample, "...")
    meta <- read.csv(file.path(P0_DIR, paste0(sample, "_spot_meta.csv")),
                     stringsAsFactors = FALSE)
    rw <- load_rctd_coords_and_weights(sample)
    if (is.null(rw)) { log_msg("   skip (no RCTD)"); next }
    
    # assign domain
    dom <- assign_spot_to_domain(rw$weights, profile_df)
    bc <- rownames(rw$weights)
    
    df <- data.frame(barcode = bc, domain = dom, stringsAsFactors = FALSE)
    df <- merge(df, meta[, c("barcode", "human_frac")], by = "barcode")
    df$sample <- sample
    results[[sample]] <- df
  }
  do.call(rbind, results)
}

plot_liver_human_frac <- function(hf_df, naming) {
  hf_df$domain_label <- naming$auto_name[match(hf_df$domain, naming$domain)]
  dom_order <- naming$domain
  hf_df$domain <- factor(hf_df$domain, levels = dom_order)
  
  p <- ggplot(hf_df, aes(x = domain, y = human_frac, fill = domain)) +
    geom_violin(scale = "width", alpha = 0.6, trim = TRUE) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.9) +
    scale_fill_manual(values = DOMAIN_COLORS, guide = "none") +
    scale_x_discrete(labels = setNames(naming$auto_name, naming$domain)) +
    facet_wrap(~ sample, ncol = 1) +
    labs(title = "Per-spot human mRNA fraction by spatial domain — Liver tumor samples",
         subtitle = "High human_frac = tumor-dominant spots. D09 [LSEC_Mid+Macrophage_M1] is tumor-induced niche.",
         x = NULL, y = "Human mRNA fraction per spot") +
    base_theme +
    theme(strip.background = element_rect(fill = "grey90", color = NA),
          strip.text = element_text(face = "bold"))
  p
}

# -------------------------------------------------------------
# STEP 5: MAIN
# -------------------------------------------------------------
log_msg("========== P1: Tumor-altered niche identification ==========")

all_stats <- list()
all_sig_niches <- list()

for (organ in c("liver", "brain")) {
  log_msg("\n----- Processing: ", organ, " -----")
  
  dat <- load_organ_data(organ)
  per_sample_annotated <- annotate_subgroup(dat$per_sample, organ)
  
  # Wilcoxon stats
  stats <- compute_wilcox(dat$per_sample)
  stats$auto_name <- dat$naming$auto_name[match(stats$domain, dat$naming$domain)]
  stats$organ <- organ
  write.csv(stats, file.path(OUT_DIR, paste0("P1_summary_", organ, ".csv")), row.names = FALSE)
  log_msg("  Wrote stats: P1_summary_", organ, ".csv")
  print(stats[, c("domain", "auto_name", "mean_tumor", "mean_base", "delta", "fc", "wilcox_p")])
  all_stats[[organ]] <- stats
  
  # Identify "altered niches": delta abs > 1.5% AND (|log2fc| > 0.3 OR p < 0.1)
  sig_niches <- stats[abs(stats$delta) > 1.5 &
                      (abs(stats$log2fc) > 0.3 | stats$wilcox_p < 0.1), ]
  sig_niches <- sig_niches[order(abs(sig_niches$delta), decreasing = TRUE), ]
  all_sig_niches[[organ]] <- sig_niches
  
  log_msg("  Significant altered niches: ", nrow(sig_niches))
  if (nrow(sig_niches) > 0) {
    cat("  Top altered:\n")
    print(sig_niches[, c("domain", "auto_name", "delta", "fc", "wilcox_p")])
  }
  
  # P1a: 3-group stacked (supplementary)
  p1a <- plot_per_sample_3groups(per_sample_annotated, organ, dat$naming)
  ggsave(file.path(OUT_DIR, paste0("P1a_", organ, "_domain_proportion_three_groups.png")),
         p1a, width = 10, height = 5.5, dpi = 150, bg = "white")
  log_msg("  Saved: P1a_", organ, "_domain_proportion_three_groups.png")
  
  # P1b: boxplot with per-sample points (MAIN)
  p1b <- plot_tumor_vs_baseline_summary(per_sample_annotated, dat$tvb, organ, dat$naming)
  ggsave(file.path(OUT_DIR, paste0("P1b_", organ, "_domain_proportion_tumor_vs_baseline.png")),
         p1b, width = 11, height = 5.5, dpi = 150, bg = "white")
  log_msg("  Saved: P1b_", organ, "_domain_proportion_tumor_vs_baseline.png")
  
  # P1c: fold change bar
  p1c <- plot_foldchange(stats, organ, dat$naming)
  ggsave(file.path(OUT_DIR, paste0("P1c_", organ, "_domain_foldchange.png")),
         p1c, width = 8, height = 5, dpi = 150, bg = "white")
  log_msg("  Saved: P1c_", organ, "_domain_foldchange.png")
}

# ---- P1d: Liver human_frac by domain ----
log_msg("\n----- Computing liver per-domain human_frac -----")
hf_df <- compute_human_frac_per_domain("liver")
if (!is.null(hf_df) && nrow(hf_df) > 0) {
  write.csv(hf_df, file.path(OUT_DIR, "P1d_liver_human_frac_by_domain.csv"), row.names = FALSE)
  
  naming_liver <- read.csv(file.path(T4_DIR, "liver_domain_naming.csv"), stringsAsFactors = FALSE)
  p1d <- plot_liver_human_frac(hf_df, naming_liver)
  ggsave(file.path(OUT_DIR, "P1d_liver_domain_human_frac.png"),
         p1d, width = 10, height = 8, dpi = 150, bg = "white")
  log_msg("  Saved: P1d_liver_domain_human_frac.png")
  
  # Summary per domain across 3 tumor samples
  hf_summary <- hf_df %>%
    group_by(domain) %>%
    summarise(
      n_spots        = n(),
      mean_human_frac = mean(human_frac, na.rm = TRUE),
      median_human_frac = median(human_frac, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_human_frac))
  hf_summary$auto_name <- naming_liver$auto_name[match(hf_summary$domain, naming_liver$domain)]
  write.csv(hf_summary, file.path(OUT_DIR, "P1d_liver_human_frac_summary.csv"), row.names = FALSE)
  log_msg("  Liver human_frac by domain summary:")
  print(hf_summary)
}

# ---- Combined sig niches ----
combined_sig <- do.call(rbind, all_sig_niches)
write.csv(combined_sig,
          file.path(OUT_DIR, "P1_significant_niches.csv"), row.names = FALSE)
log_msg("\nSaved combined significant niches: P1_significant_niches.csv")

log_msg("\n========== P1 DONE ==========")
log_msg("Output dir: ", OUT_DIR)
for (f in list.files(OUT_DIR)) {
  info <- file.info(file.path(OUT_DIR, f))
  log_msg(sprintf("  %-60s %6.2f %s",
                  f,
                  ifelse(info$size > 1e6, info$size/1e6, info$size/1e3),
                  ifelse(info$size > 1e6, "MB", "KB")))
}
