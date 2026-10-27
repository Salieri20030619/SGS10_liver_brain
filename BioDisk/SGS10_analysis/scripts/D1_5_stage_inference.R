#!/usr/bin/env Rscript
# =============================================================================
# SGS10 D1.5: Tumor stage inference for 3 tumor liver samples
# -----------------------------------------------------------------------------
# Goal: Use D1 outputs + spot-level human_pct to rank 3 tumor liver samples
#       (Liver_3, Liver_8, gan56) by infiltration stage (early/mid/late).
#
# Method: 5 independent indicators, then majority vote
#   1. mean_human_pct          : higher -> later stage
#   2. pct_tumor_spots (>50%)  : higher -> later stage
#   3. tumor_niche_proportion  : higher -> later stage (uses D1 niches)
#   4. margin_to_tumor_ratio   : LOWER  -> later stage (margin shrinks)
#   5. tumor_compactness (gini): higher -> later stage (consolidates)
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
})

option_list <- list(
  make_option("--base_dir", type="character",
              default="/mnt/BioDisk/wangmingyuan/SGS10_analysis"),
  make_option("--out_dir",  type="character",
              default="/mnt/BioDisk/wangmingyuan/SGS10_analysis/D1_niche_output"),
  make_option("--tumor_threshold", type="double", default=50,
              help="human_pct above this = tumor spot [default %default]"),
  make_option("--margin_low",  type="double", default=5,
              help="human_pct above this = at least margin [default %default]"),
  make_option("--margin_high", type="double", default=50,
              help="human_pct below this = at most margin [default %default]"),
  make_option("--tumor_samples", type="character", default="Liver_3,Liver_8,gan56"),
  make_option("--baseline_samples", type="character", default="Con_liver_1,Bla_liver_1")
)
opt <- parse_args(OptionParser(option_list=option_list))

CFG <- list(
  enriched_dir = file.path(opt$base_dir, "spatial_results/rctd/enriched"),
  out_dir      = opt$out_dir,
  niches_rds   = file.path(opt$out_dir, "niches_by_sample_liver.rds"),
  tumor_threshold = opt$tumor_threshold,
  margin_low      = opt$margin_low,
  margin_high     = opt$margin_high,
  tumor_samples    = strsplit(opt$tumor_samples, ",")[[1]],
  baseline_samples = strsplit(opt$baseline_samples, ",")[[1]]
)
all_samples <- c(CFG$baseline_samples, CFG$tumor_samples)

cat("=====================================================================\n")
cat("D1.5 Stage Inference\n")
cat("Start:", format(Sys.time()), "\n")
cat("Tumor threshold (human_pct >):", CFG$tumor_threshold, "\n")
cat("Tumor samples to rank:", paste(CFG$tumor_samples, collapse=", "), "\n")
cat("=====================================================================\n\n")

# ============ 1. Load D1 niches ============
cat("[1] Loading D1 niches RDS...\n")
if (!file.exists(CFG$niches_rds)) {
  stop("Cannot find ", CFG$niches_rds,
       "\n  -> Run D1 with --organ liver first.")
}
ni_list <- readRDS(CFG$niches_rds)
cat("  Samples in rds:", paste(names(ni_list), collapse=", "), "\n")
all_niches <- bind_rows(ni_list)
cat("  Total spots:", nrow(all_niches), "\n")
cat("  Unique niches:", length(unique(all_niches$niche)), "\n")

# ============ 2. Load human_pct + coords from enriched RDS ============
cat("\n[2] Loading meta + coords from enriched RDS...\n")

load_sample_meta <- function(s) {
  fp <- file.path(CFG$enriched_dir, paste0(s, "_seekspace_rctd.rds"))
  if (!file.exists(fp)) { warning("Missing: ", fp); return(NULL) }
  cat("  ", s, "...\n")
  so <- readRDS(fp)
  meta <- so@meta.data
  if (!"human_pct" %in% colnames(meta)) {
    warning(s, ": no human_pct"); return(NULL)
  }
  coords <- if ("spatial" %in% names(so@reductions)) {
    as.data.frame(Embeddings(so, "spatial"))
  } else {
    data.frame(x=NA, y=NA)
  }
  colnames(coords)[1:2] <- c("x","y")
  out <- data.frame(
    spot_id   = rownames(meta),
    sample    = s,
    spot_key  = paste0(s, "__", rownames(meta)),
    human_pct = meta$human_pct,
    x = coords$x, y = coords$y
  )
  rm(so); gc(verbose=FALSE)
  out
}

all_meta <- bind_rows(lapply(all_samples, load_sample_meta))
cat("  Total spots loaded:", nrow(all_meta), "\n")

# ============ 3. Indicator 1+2: per-sample human_pct stats ============
cat("\n[3] Indicators 1+2: human_pct global stats\n")

stage_stats <- all_meta %>%
  group_by(sample) %>%
  summarise(
    n_spots          = n(),
    mean_human_pct   = mean(human_pct, na.rm=TRUE),
    med_human_pct    = median(human_pct, na.rm=TRUE),
    pct_tumor_spots  = sum(human_pct > CFG$tumor_threshold, na.rm=TRUE) / n(),
    pct_margin_spots = sum(human_pct > CFG$margin_low &
                           human_pct <= CFG$margin_high, na.rm=TRUE) / n(),
    pct_host_spots   = sum(human_pct <= CFG$margin_low, na.rm=TRUE) / n(),
    .groups = "drop"
  )

cat("\n  All samples:\n")
print(as.data.frame(stage_stats))

# ============ 4. Indicator 3: tumor-enriched niche proportion ============
cat("\n[4] Indicator 3: tumor-enriched niche proportion (uses D1 niches)\n")

merged <- merge(all_niches, all_meta[, c("spot_key","human_pct")],
                by="spot_key", all.x=TRUE)

niche_global <- merged %>%
  filter(!is.na(human_pct)) %>%
  group_by(niche) %>%
  summarise(
    niche_mean_human_pct = mean(human_pct, na.rm=TRUE),
    niche_pct_tumor      = sum(human_pct > CFG$tumor_threshold, na.rm=TRUE) / n(),
    .groups="drop"
  ) %>%
  arrange(desc(niche_mean_human_pct))

cat("\n  Niche tumor enrichment (sorted):\n")
print(as.data.frame(niche_global))

tumor_niches  <- niche_global$niche[niche_global$niche_pct_tumor > 0.3]
margin_niches <- niche_global$niche[niche_global$niche_pct_tumor <= 0.3 &
                                    niche_global$niche_mean_human_pct > CFG$margin_low]
cat("\n  Tumor niches (pct_tumor > 0.3):  ", paste(tumor_niches, collapse=", "), "\n")
cat("  Margin niches (mean_pct > ", CFG$margin_low, "):", paste(margin_niches, collapse=", "), "\n")

tumor_per_sample <- merged %>%
  filter(niche %in% tumor_niches) %>%
  group_by(sample) %>%
  summarise(tumor_niche_n = n(), .groups="drop")

margin_per_sample <- merged %>%
  filter(niche %in% margin_niches) %>%
  group_by(sample) %>%
  summarise(margin_niche_n = n(), .groups="drop")

stage_stats <- stage_stats %>%
  left_join(tumor_per_sample, by="sample") %>%
  left_join(margin_per_sample, by="sample") %>%
  mutate(
    tumor_niche_n  = ifelse(is.na(tumor_niche_n),  0, tumor_niche_n),
    margin_niche_n = ifelse(is.na(margin_niche_n), 0, margin_niche_n),
    tumor_niche_proportion  = tumor_niche_n  / n_spots,
    margin_niche_proportion = margin_niche_n / n_spots
  )

# ============ 5. Indicator 4: margin/tumor ratio ============
cat("\n[5] Indicator 4: margin/tumor ratio (lower = more advanced)\n")
stage_stats <- stage_stats %>%
  mutate(
    margin_to_tumor_ratio = ifelse(tumor_niche_n > 0,
                                   margin_niche_n / tumor_niche_n, NA)
  )

# ============ 6. Indicator 5: tumor spatial compactness (Gini) ============
cat("\n[6] Indicator 5: tumor compactness Gini (higher = more focal/consolidated)\n")

gini_compactness <- function(spot_df, threshold=50, n_bins=10) {
  tumor_spots <- spot_df[spot_df$human_pct > threshold &
                         !is.na(spot_df$human_pct), ]
  if (nrow(tumor_spots) < 50) return(NA)
  x_range <- range(tumor_spots$x, na.rm=TRUE)
  y_range <- range(tumor_spots$y, na.rm=TRUE)
  if (diff(x_range) == 0 || diff(y_range) == 0) return(NA)
  x_breaks <- seq(x_range[1], x_range[2], length.out=n_bins+1)
  y_breaks <- seq(y_range[1], y_range[2], length.out=n_bins+1)
  x_bin <- cut(tumor_spots$x, breaks=x_breaks, include.lowest=TRUE)
  y_bin <- cut(tumor_spots$y, breaks=y_breaks, include.lowest=TRUE)
  tab <- as.numeric(table(paste(x_bin, y_bin)))
  if (length(tab) < 2) return(NA)
  tab <- sort(tab)
  n <- length(tab)
  numer <- sum((2 * seq_len(n) - n - 1) * tab)
  denom <- n * sum(tab)
  numer / denom
}

compactness <- all_meta %>%
  group_by(sample) %>%
  group_modify(~ data.frame(
    tumor_compactness = gini_compactness(.x,
                                          threshold = CFG$tumor_threshold,
                                          n_bins = 10)
  )) %>%
  ungroup()

stage_stats <- stage_stats %>%
  left_join(compactness, by="sample")

cat("\n=== Final stats table (all samples) ===\n")
print(as.data.frame(stage_stats))

# ============ 7. Rank tumor samples by each indicator ============
cat("\n[7] Ranking tumor samples by each indicator (1=earliest, 3=latest)\n")

tumor_only <- stage_stats %>% filter(sample %in% CFG$tumor_samples)

ranks <- data.frame(sample = tumor_only$sample)
ranks$rank_human_pct       <- rank(tumor_only$mean_human_pct,         ties.method="average")
ranks$rank_pct_tumor       <- rank(tumor_only$pct_tumor_spots,        ties.method="average")
ranks$rank_tumor_niche     <- rank(tumor_only$tumor_niche_proportion, ties.method="average")
ranks$rank_margin_to_tumor <- rank(-tumor_only$margin_to_tumor_ratio, ties.method="average")
ranks$rank_compactness     <- rank(tumor_only$tumor_compactness,      ties.method="average")

cat("\n  Per-indicator rank:\n")
print(ranks)

# Mean rank
ranks$mean_rank <- rowMeans(ranks[, grep("^rank_", colnames(ranks))], na.rm=TRUE)
ranks$inferred_stage <- factor(rank(ranks$mean_rank), levels=1:3,
                                labels=c("early","mid","late"))

cat("\n  Mean rank + inferred stage:\n")
print(ranks[, c("sample","mean_rank","inferred_stage")])

# ============ 8. Indicator agreement check ============
cat("\n[8] Indicator agreement (rank SD across 5 indicators):\n")
ranks$rank_sd <- apply(ranks[, grep("^rank_", colnames(ranks))], 1, sd, na.rm=TRUE)
print(ranks[, c("sample","rank_sd")])

agreement <- if (all(ranks$rank_sd < 0.7, na.rm=TRUE)) "STRONG" else
             if (all(ranks$rank_sd < 1.2, na.rm=TRUE)) "MODERATE" else "WEAK"

cat("\n  Overall agreement: ", agreement, "\n")
if (agreement == "STRONG") {
  cat("  -> Stage order is reliable. Use in D2-D7 stage analysis.\n")
} else if (agreement == "MODERATE") {
  cat("  -> Stage order is plausible but check diagnostic plot.\n")
} else {
  cat("  -> Indicators disagree -- treat 3 samples as independent (no stage analysis).\n")
}

# ============ 9. Output ============
cat("\n[9] Writing outputs\n")

write_csv(stage_stats, file.path(CFG$out_dir, "stage_inference_table.csv"))
cat("  -> stage_inference_table.csv\n")

final_order <- ranks[order(ranks$mean_rank), "sample"]
cat("\n  Final inferred stage order (early -> late):\n  ",
    paste(final_order, collapse=" -> "), "\n")

summary_df <- data.frame(
  sample = ranks$sample,
  inferred_stage = ranks$inferred_stage,
  mean_rank = round(ranks$mean_rank, 2),
  rank_sd = round(ranks$rank_sd, 2),
  agreement = agreement,
  recommendation = paste(final_order, collapse=" -> ")
)
write_csv(summary_df, file.path(CFG$out_dir, "stage_inference_summary.csv"))
cat("  -> stage_inference_summary.csv\n")

# ============ 10. Diagnostic plot ============
cat("\n[10] Diagnostic plot\n")

stage_long <- stage_stats %>%
  filter(sample %in% CFG$tumor_samples) %>%
  select(sample, mean_human_pct, pct_tumor_spots, tumor_niche_proportion,
         margin_to_tumor_ratio, tumor_compactness) %>%
  pivot_longer(-sample, names_to="indicator", values_to="value")

stage_long$sample <- factor(stage_long$sample, levels=final_order)

p1 <- ggplot(stage_long, aes(x=sample, y=value, group=indicator)) +
  geom_line(color="steelblue") +
  geom_point(size=3, color="firebrick") +
  facet_wrap(~indicator, scales="free_y", ncol=3) +
  labs(title="Stage inference: each indicator across samples",
       subtitle=paste("Inferred order:", paste(final_order, collapse=" -> "),
                      " | Agreement:", agreement)) +
  theme_minimal(base_size=11) +
  theme(axis.text.x = element_text(angle=30, hjust=1))

stage_with_baseline <- stage_stats %>%
  select(sample, mean_human_pct, pct_tumor_spots) %>%
  pivot_longer(-sample, names_to="indicator", values_to="value") %>%
  mutate(group = ifelse(sample %in% CFG$baseline_samples, "baseline", "tumor"))

p2 <- ggplot(stage_with_baseline, aes(x=sample, y=value, fill=group)) +
  geom_col() +
  facet_wrap(~indicator, scales="free_y") +
  scale_fill_manual(values=c("baseline"="gray70","tumor"="firebrick")) +
  labs(title="Tumor vs baseline reference") +
  theme_minimal(base_size=11) +
  theme(axis.text.x = element_text(angle=30, hjust=1))

ggsave(file.path(CFG$out_dir, "stage_inference_diag.pdf"),
       p1 / p2, width=12, height=10)
cat("  -> stage_inference_diag.pdf\n")

cat("\n=====================================================================\n")
cat("Done:", format(Sys.time()), "\n")
cat("=====================================================================\n")
cat("\nNEXT:\n")
cat("  - Check stage_inference_diag.pdf visually\n")
cat("  - If agreement >= MODERATE: use this order in D2 with --stage_groups\n")
cat("  - Inferred stage order:", paste(final_order, collapse=" -> "), "\n\n")
