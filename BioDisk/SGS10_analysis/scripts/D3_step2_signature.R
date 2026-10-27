#!/usr/bin/env Rscript
# =============================================================================
# SGS10 D3 step2: Null-filtered signatures + UCell scoring + per-animal correlation
# -----------------------------------------------------------------------------
# Pipeline:
#   1. Filter Wilcoxon DEG by removing null (Bla-vs-Con baseline) DEG
#   2. Build top-N UP signatures per contrast (A1/A2/A3/B3/C4)
#   3. UCell scoring per spot (across all spots, both organs)
#   4. Per-animal aggregation (3 paired animals: animal_8/3/56)
#   5. Cross-organ correlation: liver invaded score vs brain distant score
#   6. Per-stage trend plots
#
# Inputs:
#   - D3_DEG_output/wilcoxon/<contrast>.csv    (per-contrast Wilcoxon DEG)
#   - D3_DEG_output/wilcoxon/_NULL_*.csv       (baseline-vs-baseline null DEGs)
#   - D2 spot_zone_stage_{liver,brain}.rds     (zone+stage+animal labels)
#   - enriched RDS                              (raw counts for UCell)
#
# Outputs:
#   D3_signature_output/
#     filtered_DEG/<contrast>_UP_filtered.csv  (true tumor-driven UP, null removed)
#     signatures/<contrast>_top<N>_UP.txt      (gene list for UCell)
#     spot_scores_liver.csv                     (per-spot UCell scores)
#     spot_scores_brain.csv
#     animal_summary.csv                        (per-animal mean scores per signature)
#     figures/
#       per_animal_correlation.pdf  ← paper Figure 3 main
#       per_stage_trend.pdf
#       signature_score_heatmap.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(UCell)
})

option_list <- list(
  make_option("--base_dir", type="character",
              default="/mnt/BioDisk/wangmingyuan/SGS10_analysis"),
  make_option("--d2_dir",   type="character",
              default="/mnt/BioDisk/wangmingyuan/SGS10_analysis/D2_zone_output"),
  make_option("--d3_dir",   type="character",
              default="/mnt/BioDisk/wangmingyuan/SGS10_analysis/D3_DEG_output"),
  make_option("--out_dir",  type="character",
              default="/mnt/BioDisk/wangmingyuan/SGS10_analysis/D3_signature_output"),
  make_option("--top_n", type="integer", default=50,
              help="top N UP genes per signature"),
  make_option("--min_consistency", type="double", default=0.67,
              help="min sample_consistency for inclusion"),
  make_option("--min_log2fc", type="double", default=0.25),
  make_option("--ncores", type="integer", default=4),
  make_option("--dry_run", action="store_true", default=FALSE)
)
opt <- parse_args(OptionParser(option_list=option_list))

CFG <- list(
  base_dir = opt$base_dir,
  enriched_dir = file.path(opt$base_dir, "spatial_results/rctd/enriched"),
  d2_dir = opt$d2_dir,
  d3_dir = opt$d3_dir,
  out_dir = opt$out_dir,
  top_n = opt$top_n,
  min_consistency = opt$min_consistency,
  min_log2fc = opt$min_log2fc,
  ncores = opt$ncores,
  dry_run = opt$dry_run,
  liver_samples = c("Con_liver_1","Bla_liver_1","Liver_3","Liver_8","gan56"),
  brain_samples = c("Con_brain_1","Bla_brain_1","Brain_3","Brain_8","nao56"),
  # Contrasts and which null reference to use
  contrast_info = list(
    A1 = list(file="A1_tumor_core_vs_baseline_host.csv",   organ="liver", null="liver"),
    A2 = list(file="A2_invaded_host_vs_baseline_host.csv", organ="liver", null="liver"),
    A3 = list(file="A3_distant_host_vs_baseline_host.csv", organ="liver", null="liver"),
    B3 = list(file="B3_distant_late_vs_distant_early.csv", organ="liver", null=NA),  # no null
    C4 = list(file="C4_brain_distant_vs_brain_baseline.csv", organ="brain", null="brain")
  )
)

dir.create(CFG$out_dir, showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(CFG$out_dir, "filtered_DEG"), showWarnings=FALSE)
dir.create(file.path(CFG$out_dir, "signatures"),   showWarnings=FALSE)
dir.create(file.path(CFG$out_dir, "figures"),      showWarnings=FALSE)

cat_log <- function(...) cat(..., sep="")
sec <- function(name) {
  cat_log("\n", strrep("-", 70), "\n  ", name, "\n",
          strrep("-", 70), "\n")
}

cat_log("=====================================================================\n")
cat_log("D3 Step 2: signature scoring + per-animal correlation\n")
cat_log("Start: ", format(Sys.time()), "\n")
cat_log("Top N per signature: ", CFG$top_n, "\n")
cat_log("Min sample consistency: ", CFG$min_consistency, "\n")
cat_log("=====================================================================\n")

# ============ Step 1: Null-filter + signature extraction ============
sec("Step 1: Null-filter and build signatures")

# Load null DEGs
null_genes <- list()
for (org in c("liver","brain")) {
  fp <- file.path(CFG$d3_dir, "wilcoxon", paste0("_NULL_Bla_vs_Con_", org, ".csv"))
  if (file.exists(fp)) {
    df <- read_csv(fp, show_col_types=FALSE)
    null_genes[[org]] <- df$gene[df$p_val_adj < 0.05]
    cat_log("  null ", org, ": ", length(null_genes[[org]]), " genes (UP+DOWN)\n")
  } else {
    cat_log("  WARNING: null file missing for ", org, ": ", fp, "\n")
    null_genes[[org]] <- character(0)
  }
}

# Process each contrast
signatures_up   <- list()
signatures_down <- list()
filter_summary  <- list()

for (ct_id in names(CFG$contrast_info)) {
  ci <- CFG$contrast_info[[ct_id]]
  fp <- file.path(CFG$d3_dir, "wilcoxon", ci$file)
  if (!file.exists(fp)) {
    cat_log("  ", ct_id, ": missing Wilcoxon file, skip\n")
    next
  }
  df <- read_csv(fp, show_col_types=FALSE)
  
  # Robust filter
  up_robust <- df %>%
    filter(avg_log2FC > CFG$min_log2fc &
           p_val_adj < 0.05 &
           !is.na(sample_consistency) &
           sample_consistency >= CFG$min_consistency)
  down_robust <- df %>%
    filter(avg_log2FC < -CFG$min_log2fc &
           p_val_adj < 0.05 &
           !is.na(sample_consistency) &
           sample_consistency >= CFG$min_consistency)
  
  # Null filter
  null_g <- if (!is.na(ci$null)) null_genes[[ci$null]] else character(0)
  up_filtered   <- up_robust %>% filter(!gene %in% null_g)
  down_filtered <- down_robust %>% filter(!gene %in% null_g)
  
  cat_log("  ", ct_id, ": UP robust=", nrow(up_robust),
          " -> filtered=", nrow(up_filtered),
          " | DOWN robust=", nrow(down_robust),
          " -> filtered=", nrow(down_filtered),
          " (null=", ci$null, ", removed=", length(null_g), ")\n")
  
  # Save filtered DEG
  out_up <- file.path(CFG$out_dir, "filtered_DEG",
                     paste0(ct_id, "_UP_filtered.csv"))
  write_csv(up_filtered %>% arrange(desc(avg_log2FC)), out_up)
  out_dn <- file.path(CFG$out_dir, "filtered_DEG",
                     paste0(ct_id, "_DOWN_filtered.csv"))
  write_csv(down_filtered %>% arrange(avg_log2FC), out_dn)
  
  # Top N UP gene list for UCell
  top_up <- up_filtered %>%
    arrange(desc(avg_log2FC)) %>%
    head(CFG$top_n) %>%
    pull(gene)
  
  signatures_up[[ct_id]] <- top_up
  
  # Save gene list
  writeLines(top_up,
             file.path(CFG$out_dir, "signatures",
                       paste0(ct_id, "_top", CFG$top_n, "_UP.txt")))
  
  filter_summary[[ct_id]] <- data.frame(
    contrast = ct_id,
    organ = ci$organ,
    n_robust_UP = nrow(up_robust),
    n_filtered_UP = nrow(up_filtered),
    n_robust_DOWN = nrow(down_robust),
    n_filtered_DOWN = nrow(down_filtered),
    n_signature_UP = length(top_up),
    null_used = if (is.na(ci$null)) "none" else ci$null,
    n_null_genes = length(null_g)
  )
}

write_csv(bind_rows(filter_summary),
          file.path(CFG$out_dir, "filter_summary.csv"))

# Print signature previews
sec("Signature previews")
for (ct_id in names(signatures_up)) {
  sig <- signatures_up[[ct_id]]
  cat_log("\n[", ct_id, "] top ", min(10, length(sig)), " UP genes:\n")
  cat_log("  ", paste(head(sig, 10), collapse=", "), "\n")
}

if (CFG$dry_run) {
  cat_log("\n[DRY RUN] OK, signatures built\n")
  quit(save="no", status=0)
}

# ============ Step 2: UCell scoring per organ ============
sec("Step 2: UCell scoring (per organ)")

load_organ_for_scoring <- function(organ) {
  samples <- if (organ == "liver") CFG$liver_samples else CFG$brain_samples
  cat_log("  loading ", organ, " (mouse-only)...\n")
  
  so_list <- list()
  for (s in samples) {
    fp <- file.path(CFG$enriched_dir, paste0(s, "_seekspace_rctd.rds"))
    if (!file.exists(fp)) next
    so <- readRDS(fp)
    counts <- GetAssayData(so, layer = "counts")
    feats <- rownames(counts)
    if (any(grepl("^mm10---", feats))) {
      keep <- grep("^mm10---", feats)
      counts <- counts[keep, ]
      rownames(counts) <- sub("^mm10---", "", rownames(counts))
    }
    new_so <- CreateSeuratObject(counts = counts)
    new_so$sample_id <- s
    so_list[[s]] <- new_so
    rm(so, counts); gc(verbose=FALSE)
  }
  merged <- merge(so_list[[1]], y = so_list[-1],
                  add.cell.ids = names(so_list))
  rm(so_list); gc(verbose=FALSE)
  merged <- JoinLayers(merged)
  merged <- NormalizeData(merged, verbose = FALSE)
  
  zsd_path <- file.path(CFG$d2_dir, paste0("spot_zone_stage_", organ, ".rds"))
  zsd <- readRDS(zsd_path)
  matched <- intersect(colnames(merged), zsd$spot_key)
  merged <- merged[, matched]
  zsd_ord <- zsd[match(colnames(merged), zsd$spot_key), ]
  merged$zone        <- zsd_ord$zone
  merged$stage_label <- zsd_ord$stage_label
  merged$animal_id   <- zsd_ord$animal_id
  merged$human_pct   <- zsd_ord$human_pct
  
  cat_log("    ", organ, " merged: ", ncol(merged), " spots\n")
  merged
}

# Score each organ with all signatures
score_one_organ <- function(seurat_obj, organ_label, signatures_list) {
  cat_log("  Scoring ", organ_label, " with ",
          length(signatures_list), " signatures...\n")
  t0 <- Sys.time()
  seurat_obj <- AddModuleScore_UCell(
    obj = seurat_obj,
    features = signatures_list,
    ncores = CFG$ncores,
    name = ""
  )
  dt <- round(as.numeric(Sys.time()-t0, units="secs"), 1)
  cat_log("    UCell done. ", round(dt/60, 1), " min\n")
  seurat_obj
}

# Get organs needed
organs_needed <- unique(sapply(CFG$contrast_info, function(ci) ci$organ))

organ_objs <- list()
for (org in organs_needed) {
  organ_objs[[org]] <- load_organ_for_scoring(org)
}

# Run UCell on liver (with all signatures) and brain (same set)
# This way every spot gets every signature score, useful for cross-organ analysis
all_sigs_for_scoring <- signatures_up
# Filter empty signatures
all_sigs_for_scoring <- all_sigs_for_scoring[
  sapply(all_sigs_for_scoring, length) >= 5]
cat_log("  signatures with >=5 genes: ", length(all_sigs_for_scoring), " / ",
        length(signatures_up), "\n")

for (org in names(organ_objs)) {
  organ_objs[[org]] <- score_one_organ(organ_objs[[org]], org, all_sigs_for_scoring)
}

# Save spot-level scores
sec("Save spot-level scores")
for (org in names(organ_objs)) {
  so <- organ_objs[[org]]
  sig_cols <- names(all_sigs_for_scoring)
  
  meta_out <- data.frame(
    spot_key   = colnames(so),
    sample     = so$sample_id,
    zone       = as.character(so$zone),
    stage      = as.character(so$stage_label),
    animal_id  = so$animal_id,
    human_pct  = so$human_pct,
    stringsAsFactors = FALSE
  )
  for (sig in sig_cols) {
    meta_out[[paste0("score_", sig)]] <- so[[sig]][, 1]
  }
  
  out_path <- file.path(CFG$out_dir, paste0("spot_scores_", org, ".csv"))
  write_csv(meta_out, out_path)
  cat_log("  saved: ", out_path, " (", nrow(meta_out), " rows)\n")
}

# ============ Step 3: Per-animal aggregation ============
sec("Step 3: Per-animal score aggregation")

per_animal_summary <- list()
for (org in names(organ_objs)) {
  so <- organ_objs[[org]]
  sig_cols <- names(all_sigs_for_scoring)
  
  for (sig in sig_cols) {
    df <- data.frame(
      organ = org,
      sample = so$sample_id,
      animal_id = so$animal_id,
      zone = as.character(so$zone),
      stage = as.character(so$stage_label),
      score = so[[sig]][, 1],
      signature = sig,
      stringsAsFactors = FALSE
    )
    
    # Aggregate per (animal, zone)
    agg <- df %>%
      group_by(organ, animal_id, zone, stage, signature) %>%
      summarise(
        mean_score = mean(score, na.rm=TRUE),
        median_score = median(score, na.rm=TRUE),
        sd_score = sd(score, na.rm=TRUE),
        n_spots = n(),
        .groups = "drop"
      )
    per_animal_summary[[paste0(org, "_", sig)]] <- agg
  }
}
animal_df <- bind_rows(per_animal_summary)
write_csv(animal_df,
          file.path(CFG$out_dir, "animal_summary.csv"))
cat_log("  saved: animal_summary.csv (", nrow(animal_df), " rows)\n")

# ============ Step 4: Per-animal correlation (paper Figure 3 main) ============
sec("Step 4: Per-animal cross-organ correlation")

# For each pair (liver_signature x brain_signature), compute per-animal score
# in invaded_host (liver) and distant_host (brain), then correlate across 3 animals

paired_animals <- c("animal_8", "animal_3", "animal_56")
stage_for_animal <- c(animal_8="early", animal_3="mid", animal_56="late")

cor_results <- list()
plot_data <- list()

for (sig_name in names(all_sigs_for_scoring)) {
  # Liver-side: in invaded_host (or other relevant zone)
  liver_zone_for_sig <- "invaded_host"  # default; override per sig if needed
  if (grepl("^A1", sig_name)) liver_zone_for_sig <- "tumor_core"
  if (grepl("^A2", sig_name)) liver_zone_for_sig <- "invaded_host"
  if (grepl("^A3", sig_name)) liver_zone_for_sig <- "distant_host"
  if (grepl("^B3", sig_name)) liver_zone_for_sig <- "distant_host"
  if (grepl("^C4", sig_name)) liver_zone_for_sig <- "invaded_host"  # use invaded for cross-test
  
  # Brain-side: distant_host always
  brain_zone <- "distant_host"
  
  liver_score_animal <- animal_df %>%
    filter(organ == "liver",
           zone == liver_zone_for_sig,
           signature == sig_name,
           animal_id %in% paired_animals) %>%
    select(animal_id, liver_score = mean_score, n_spots_liver = n_spots)
  
  brain_score_animal <- animal_df %>%
    filter(organ == "brain",
           zone == brain_zone,
           signature == sig_name,
           animal_id %in% paired_animals) %>%
    select(animal_id, brain_score = mean_score, n_spots_brain = n_spots)
  
  paired <- inner_join(liver_score_animal, brain_score_animal, by="animal_id")
  
  if (nrow(paired) < 3) {
    cat_log("  ", sig_name, ": only ", nrow(paired), " paired animals, skip\n")
    next
  }
  
  paired$stage <- stage_for_animal[paired$animal_id]
  paired$signature <- sig_name
  paired$liver_zone <- liver_zone_for_sig
  paired$brain_zone <- brain_zone
  
  # Spearman correlation
  cor_test <- cor.test(paired$liver_score, paired$brain_score, method="spearman")
  
  cat_log("  ", sig_name, " (liver_", liver_zone_for_sig,
          " vs brain_", brain_zone, "): Spearman rho = ",
          round(cor_test$estimate, 3), ", p = ",
          round(cor_test$p.value, 3), "\n")
  
  cor_results[[sig_name]] <- data.frame(
    signature = sig_name,
    liver_zone = liver_zone_for_sig,
    brain_zone = brain_zone,
    n_animals = nrow(paired),
    spearman_rho = cor_test$estimate,
    spearman_p = cor_test$p.value
  )
  plot_data[[sig_name]] <- paired
}

write_csv(bind_rows(cor_results),
          file.path(CFG$out_dir, "per_animal_correlation.csv"))

# Per-animal correlation scatter plot
sec("Step 5: Plotting")
plot_df <- bind_rows(plot_data)
if (nrow(plot_df) > 0) {
  cor_df <- bind_rows(cor_results)
  p_cor <- ggplot(plot_df, aes(x=liver_score, y=brain_score)) +
    geom_point(aes(color=stage), size=4) +
    geom_text(aes(label=animal_id), vjust=-1, size=3) +
    geom_smooth(method="lm", se=FALSE, color="black",
                linetype="dashed", size=0.5) +
    scale_color_manual(values=c(early="#2c7bb6", mid="#fdae61", late="#d7191c")) +
    facet_wrap(~ signature, scales="free", ncol=2) +
    labs(title="Per-animal cross-organ signature scores",
         subtitle="3 paired animals (animal_8/3/56), liver vs brain",
         x="Liver score (mean per animal in liver target zone)",
         y="Brain score (mean per animal in brain distant_host)") +
    theme_minimal(base_size=10)
  ggsave(file.path(CFG$out_dir, "figures", "per_animal_correlation.pdf"),
         p_cor, width=10, height=10, limitsize=FALSE)
  cat_log("  saved: per_animal_correlation.pdf\n")
}

# ============ Step 5: Per-stage trend plots ============

trend_df <- animal_df %>%
  filter(stage != "baseline" & !is.na(stage)) %>%
  mutate(stage = factor(stage, levels=c("early","mid","late")))

if (nrow(trend_df) > 0) {
  p_trend <- ggplot(trend_df,
                    aes(x=stage, y=mean_score, color=zone, group=zone)) +
    geom_line(size=0.8) +
    geom_point(size=2.5) +
    facet_grid(signature ~ organ, scales="free_y") +
    labs(title="Signature scores by stage",
         x="stage", y="mean UCell score per animal-zone") +
    theme_minimal(base_size=9) +
    theme(legend.position="bottom",
          strip.text.y=element_text(size=7))
  ggsave(file.path(CFG$out_dir, "figures", "per_stage_trend.pdf"),
         p_trend, width=12, height=14, limitsize=FALSE)
  cat_log("  saved: per_stage_trend.pdf\n")
}

# ============ Step 6: Score heatmap (zone × stage × signature) ============

heat_df <- animal_df %>%
  group_by(organ, zone, stage, signature) %>%
  summarise(mean_score = mean(mean_score, na.rm=TRUE), .groups="drop") %>%
  filter(!is.na(stage))

if (nrow(heat_df) > 0) {
  heat_df$cell_id <- paste0(heat_df$zone, "__", heat_df$stage)
  
  p_heat <- ggplot(heat_df, aes(x=cell_id, y=signature, fill=mean_score)) +
    geom_tile(color="white") +
    scale_fill_gradient2(low="#2c7bb6", mid="white", high="#d7191c",
                         midpoint=median(heat_df$mean_score, na.rm=TRUE)) +
    facet_wrap(~ organ, scales="free_x") +
    labs(title="Signature score heatmap (zone × stage)",
         x="cell (zone__stage)", y="signature", fill="mean score") +
    theme_minimal(base_size=9) +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=8),
          panel.grid=element_blank())
  ggsave(file.path(CFG$out_dir, "figures", "signature_score_heatmap.pdf"),
         p_heat, width=14, height=6, limitsize=FALSE)
  cat_log("  saved: signature_score_heatmap.pdf\n")
}

cat_log("\n=== DONE: ", format(Sys.time()), " ===\n")
cat_log("Output: ", CFG$out_dir, "\n\n")
