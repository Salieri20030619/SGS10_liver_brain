#!/usr/bin/env Rscript
# =============================================================================
# SGS10 D2 step1: Spot-level zone + stage assignment
# -----------------------------------------------------------------------------
# Input: 
#   - D1 niches_by_sample_liver.rds (per-spot niche labels + coords)
#   - per-sample enriched RDS (for human_pct field, already in coords merged)
#
# Output:
#   - spot_zone_stage_liver.rds (master spot-level annotation table)
#   - spot_zone_stage_liver.csv (human-readable)
#
# Logic:
#   Step 1: Load D1 niches + add human_pct from enriched RDS
#   Step 2: Assign zone based on human_pct
#       - tumor_core:    human_pct > 50
#       - tumor_margin:  5 < human_pct <= 50
#       - host_raw:      human_pct <= 5 (further split into invaded vs distant)
#   Step 3: For each "host_raw" spot, do spatial kNN (k=12 default)
#       - if any neighbor has human_pct > 5 -> invaded_host
#       - else -> distant_host
#       - Skip baseline samples (their host is always "baseline_host")
#   Step 4: Assign stage_label
#       - Bla_liver_1, Con_liver_1 -> baseline (no stage)
#       - Liver_8 -> early   (stage_rank = 1)
#       - Liver_3 -> mid     (stage_rank = 2)
#       - gan56   -> late    (stage_rank = 3)
#   Step 5: Final 12-cell matrix + baseline
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(FNN)
})

option_list <- list(
  make_option("--organ",    type="character", default="liver"),
  make_option("--base_dir", type="character",
              default="/mnt/BioDisk/wangmingyuan/SGS10_analysis"),
  make_option("--d1_dir",   type="character",
              default="/mnt/BioDisk/wangmingyuan/SGS10_analysis/D1_niche_output"),
  make_option("--out_dir",  type="character",
              default="/mnt/BioDisk/wangmingyuan/SGS10_analysis/D2_zone_output"),
  # Zone thresholds
  make_option("--tumor_core_min",   type="double", default=50.0,
              help="human_pct > this -> tumor_core"),
  make_option("--tumor_margin_min", type="double", default=5.0,
              help="human_pct in (this, tumor_core_min] -> tumor_margin"),
  # Spatial kNN for invaded vs distant host
  make_option("--knn_k",            type="integer", default=12,
              help="spatial kNN neighbors for invaded/distant host classification"),
  make_option("--knn_invaded_min_n",type="integer", default=1,
              help="minimum number of neighbors with human_pct>tumor_margin_min to call invaded"),
  make_option("--max_invaded_dist", type="double", default=NA,
              help="optional: spatial distance ceiling (in pixel units) above which spot is forced distant_host. NA = no cap"),
  # Output
  make_option("--save_rds",         action="store_true", default=TRUE),
  make_option("--dry_run",          action="store_true", default=FALSE)
)
opt <- parse_args(OptionParser(option_list=option_list))

CFG <- list(
  organ = opt$organ,
  base_dir = opt$base_dir,
  enriched_dir = file.path(opt$base_dir, "spatial_results/rctd/enriched"),
  d1_dir = opt$d1_dir,
  out_dir = opt$out_dir,
  tumor_core_min   = opt$tumor_core_min,
  tumor_margin_min = opt$tumor_margin_min,
  knn_k = opt$knn_k,
  knn_invaded_min_n = opt$knn_invaded_min_n,
  max_invaded_dist = opt$max_invaded_dist,
  save_rds = opt$save_rds,
  dry_run = opt$dry_run,
  liver_samples = c("Con_liver_1","Bla_liver_1","Liver_3","Liver_8","gan56"),
  brain_samples = c("Con_brain_1","Bla_brain_1","Brain_3","Brain_8","nao56"),
  baseline_samples = c("Con_liver_1","Bla_liver_1","Con_brain_1","Bla_brain_1"),
  # stage assignment (liver + brain)
  stage_map = list(
    Con_liver_1 = list(label="baseline", rank=NA_integer_),
    Bla_liver_1 = list(label="baseline", rank=NA_integer_),
    Liver_8     = list(label="early",    rank=1L),
    Liver_3     = list(label="mid",      rank=2L),
    gan56       = list(label="late",     rank=3L),
    Con_brain_1 = list(label="baseline", rank=NA_integer_),
    Bla_brain_1 = list(label="baseline", rank=NA_integer_),
    Brain_8     = list(label="early",    rank=1L),
    Brain_3     = list(label="mid",      rank=2L),
    nao56       = list(label="late",     rank=3L)
  ),
  # animal pairing: 3 paired tumor animals (liver-brain), 4 unpaired baselines
  # Confirmed: Brain_8 == Liver_8, Brain_3 == Liver_3, nao56 == gan56
  # Bla / Con baselines: liver and brain are NOT same animal
  animal_map = list(
    Liver_8     = "animal_8",   Brain_8     = "animal_8",
    Liver_3     = "animal_3",   Brain_3     = "animal_3",
    gan56       = "animal_56",  nao56       = "animal_56",
    Bla_liver_1 = "ctrl_Bla_liver",
    Con_liver_1 = "ctrl_Con_liver",
    Bla_brain_1 = "ctrl_Bla_brain",
    Con_brain_1 = "ctrl_Con_brain"
  )
)

dir.create(CFG$out_dir, showWarnings=FALSE, recursive=TRUE)

cat_log <- function(...) cat(..., sep="")
sec <- function(name) {
  cat_log("\n", strrep("-", 70), "\n  ", name, "\n",
          strrep("-", 70), "\n")
}

cat_log("=====================================================================\n")
cat_log("D2 Step 1: Spot-level zone + stage assignment\n")
cat_log("Organ: ", CFG$organ, "\n")
cat_log("Start: ", format(Sys.time()), "\n")
cat_log("Zone thresholds: tumor_core>", CFG$tumor_core_min,
        "%, tumor_margin (",CFG$tumor_margin_min,"-",CFG$tumor_core_min,"]%\n")
cat_log("kNN for invaded/distant host: k=", CFG$knn_k,
        ", min_invaded_n=", CFG$knn_invaded_min_n, "\n")
cat_log("=====================================================================\n")

samples <- if (CFG$organ == "liver") CFG$liver_samples else CFG$brain_samples

# ============ Step 1: Load D1 niches ============
sec("Step 1: Load D1 niches")
d1_rds <- file.path(CFG$d1_dir, paste0("niches_by_sample_", CFG$organ, ".rds"))
if (!file.exists(d1_rds)) stop("D1 RDS not found: ", d1_rds)

ni_list <- readRDS(d1_rds)
ni_all <- bind_rows(ni_list)
cat_log("  loaded ", nrow(ni_all), " spots from ",
        length(unique(ni_all$sample)), " samples\n")
cat_log("  niches:\n")
print(table(ni_all$niche))

# ============ Step 2: Pull human_pct from enriched RDS ============
sec("Step 2: Add human_pct + extra meta from enriched RDS")
hp_list <- list()
for (s in samples) {
  fp <- file.path(CFG$enriched_dir, paste0(s, "_seekspace_rctd.rds"))
  if (!file.exists(fp)) {
    warning("missing ", fp); next
  }
  cat_log("  loading meta ", s, " ...\n")
  so <- readRDS(fp)
  meta <- so@meta.data
  # Need: human_pct + cell barcode (origname)
  meta$spot_id_orig <- rownames(meta)
  keep_cols <- intersect(c("human_pct","human_count","mouse_count",
                           "species_class","group"), colnames(meta))
  if (length(keep_cols) == 0) {
    warning(s, ": no human_pct found")
    rm(so); next
  }
  out <- meta[, c("spot_id_orig", keep_cols), drop=FALSE]
  out$sample <- s
  # Reconstruct spot_key the same way D1 step1 did: paste0(sample, "_", origname)
  out$spot_key <- paste0(s, "_", out$spot_id_orig)
  hp_list[[s]] <- out
  rm(so, meta); gc(verbose=FALSE)
  cat_log("    ", s, ": ", nrow(out), " spots, fields: ",
          paste(keep_cols, collapse=","), "\n")
}
hp_df <- bind_rows(hp_list)
cat_log("\n  total: ", nrow(hp_df), " meta rows\n")

# Merge with niches by spot_key
n_before <- nrow(ni_all)
ni_all <- ni_all %>%
  left_join(hp_df %>% select(-sample, -spot_id_orig),
            by = "spot_key")
n_after <- nrow(ni_all)
cat_log("  after merge: ", n_after, " (was ", n_before, ")\n")

n_missing <- sum(is.na(ni_all$human_pct))
if (n_missing > 0) {
  cat_log("  WARNING: ", n_missing, " spots have NA human_pct (",
          round(100*n_missing/n_after,1), "%)\n")
  cat_log("  examples of missing:\n")
  print(head(ni_all[is.na(ni_all$human_pct), c("sample","spot_key")]))
}

# ============ Step 3: Assign zone ============
sec("Step 3: Assign zone")
ni_all$zone_raw <- "host_raw"
ni_all$zone_raw[!is.na(ni_all$human_pct) &
                ni_all$human_pct > CFG$tumor_core_min] <- "tumor_core"
ni_all$zone_raw[!is.na(ni_all$human_pct) &
                ni_all$human_pct > CFG$tumor_margin_min &
                ni_all$human_pct <= CFG$tumor_core_min] <- "tumor_margin"

cat_log("  zone_raw distribution:\n")
print(table(ni_all$zone_raw, useNA="ifany"))

cat_log("\n  by sample:\n")
print(table(ni_all$sample, ni_all$zone_raw, useNA="ifany"))

# ============ Step 4: Spatial kNN for invaded vs distant host ============
sec("Step 4: Spatial kNN for invaded vs distant host (per-sample)")

ni_all$zone <- ni_all$zone_raw  # default
# Initialize new fields
ni_all$min_dist_to_tumor <- NA_real_
ni_all$n_tumor_neighbors <- NA_integer_

for (s in samples) {
  is_baseline <- s %in% CFG$baseline_samples
  spt_idx <- which(ni_all$sample == s)
  if (length(spt_idx) == 0) next
  
  if (is_baseline) {
    # baseline: all host_raw -> baseline_host
    base_idx <- spt_idx[ni_all$zone_raw[spt_idx] == "host_raw"]
    ni_all$zone[base_idx] <- "baseline_host"
    cat_log("  ", s, " (baseline): ", length(base_idx),
            " host spots -> baseline_host\n")
    next
  }
  
  # Tumor sample: split host_raw into invaded vs distant
  sub_coords <- as.matrix(ni_all[spt_idx, c("x", "y")])
  sub_human  <- ni_all$human_pct[spt_idx]
  sub_zone   <- ni_all$zone_raw[spt_idx]
  
  is_tumor_or_margin <- !is.na(sub_human) & sub_human > CFG$tumor_margin_min
  is_host_raw <- sub_zone == "host_raw"
  
  cat_log("  ", s, ": ", sum(is_host_raw), " host_raw spots, ",
          sum(is_tumor_or_margin), " tumor+margin spots; running kNN k=",
          CFG$knn_k, "\n")
  
  if (sum(is_tumor_or_margin) == 0) {
    ni_all$zone[spt_idx[is_host_raw]] <- "distant_host"
    cat_log("    no tumor/margin in sample, all host -> distant_host\n")
    next
  }
  
  # For each host_raw spot, find k nearest neighbors among ALL spots (in this sample)
  # then count how many neighbors are tumor_core or tumor_margin
  host_idx_local <- which(is_host_raw)
  if (length(host_idx_local) == 0) next
  
  k_use <- min(CFG$knn_k + 1, nrow(sub_coords))  # +1 because nearest is self
  nn <- get.knnx(sub_coords, sub_coords[host_idx_local, , drop=FALSE], k = k_use)
  # nn$nn.index: rows = host spots, cols = k neighbors (1st col = self when query==data)
  # nn$nn.dist : same shape
  # Drop self (1st col)
  nn_idx <- nn$nn.index[, -1, drop=FALSE]
  nn_dst <- nn$nn.dist[, -1, drop=FALSE]
  
  # For each row, count neighbors that are tumor_or_margin
  is_tom_vec <- is_tumor_or_margin
  n_tumor_nb <- rowSums(matrix(is_tom_vec[nn_idx], nrow=nrow(nn_idx)))
  
  # Min distance to nearest tumor/margin neighbor
  # mask out non-tumor neighbors with Inf, take rowMin
  dist_mask <- nn_dst
  dist_mask[!is_tom_vec[nn_idx]] <- Inf
  min_dist <- apply(dist_mask, 1, min)
  min_dist[is.infinite(min_dist)] <- NA_real_
  
  # Classify
  is_invaded <- n_tumor_nb >= CFG$knn_invaded_min_n
  if (!is.na(CFG$max_invaded_dist)) {
    too_far <- !is.na(min_dist) & min_dist > CFG$max_invaded_dist
    is_invaded <- is_invaded & !too_far
  }
  
  # Update global
  ni_all$zone[spt_idx[host_idx_local[is_invaded]]]  <- "invaded_host"
  ni_all$zone[spt_idx[host_idx_local[!is_invaded]]] <- "distant_host"
  ni_all$min_dist_to_tumor[spt_idx[host_idx_local]] <- min_dist
  ni_all$n_tumor_neighbors[spt_idx[host_idx_local]] <- n_tumor_nb
  
  cat_log("    invaded_host: ", sum(is_invaded),
          " (", round(100*sum(is_invaded)/length(is_invaded), 1), "%)\n")
  cat_log("    distant_host: ", sum(!is_invaded),
          " (", round(100*sum(!is_invaded)/length(is_invaded), 1), "%)\n")
}

cat_log("\n  Final zone distribution:\n")
print(table(ni_all$zone))

cat_log("\n  zone x sample:\n")
print(table(ni_all$zone, ni_all$sample))

# ============ Step 5: Assign stage ============
sec("Step 5: Assign stage_label + stage_rank")
get_stage_field <- function(s, field) {
  m <- CFG$stage_map[[s]]
  if (is.null(m)) {
    if (field == "label") {
      return(NA_character_)
    } else {
      return(NA_integer_)
    }
  }
  v <- m[[field]]
  if (is.null(v)) {
    if (field == "label") {
      return(NA_character_)
    } else {
      return(NA_integer_)
    }
  }
  v
}
ni_all$stage_label <- vapply(ni_all$sample, get_stage_field, character(1),
                              field = "label")
ni_all$stage_rank  <- vapply(ni_all$sample, get_stage_field, integer(1),
                              field = "rank")
# animal_id (for paired liver-brain analysis in D3)
get_animal_id <- function(s) {
  v <- CFG$animal_map[[s]]
  if (is.null(v)) NA_character_ else v
}
ni_all$animal_id <- vapply(ni_all$sample, get_animal_id, character(1))

cat_log("  stage_label distribution:\n")
print(table(ni_all$stage_label, useNA="ifany"))

cat_log("\n  animal_id distribution:\n")
print(table(ni_all$animal_id, useNA="ifany"))

cat_log("\n  animal_id x sample (sanity check):\n")
print(table(ni_all$animal_id, ni_all$sample, useNA="ifany"))

# ============ Step 6: Final 12-cell + baseline ============
sec("Step 6: 12-cell matrix + baseline summary")
ni_all$cell_id <- paste0(ni_all$zone, "__", ni_all$stage_label)

cat_log("  cell_id (zone__stage) distribution:\n")
cell_tab <- table(ni_all$cell_id)
cell_tab <- sort(cell_tab, decreasing=TRUE)
print(cell_tab)

# More useful: 4 zones x 4 stage groups (3 stages + baseline)
sec("  Pivot table (zone x stage)")
pivot <- ni_all %>%
  count(zone, stage_label) %>%
  pivot_wider(names_from = stage_label, values_from = n, values_fill = 0)
print(as.data.frame(pivot))

# Save pivot too
pivot_long <- ni_all %>% count(zone, stage_label, sample)

# ============ Step 7: Save outputs ============
sec("Step 7: Save outputs")

out_csv_path <- file.path(CFG$out_dir,
                          paste0("spot_zone_stage_", CFG$organ, ".csv"))
write_csv(ni_all, out_csv_path)
cat_log("  CSV: ", out_csv_path, " (",
        round(file.info(out_csv_path)$size / 1024^2, 1), " MB)\n")

if (CFG$save_rds) {
  out_rds_path <- file.path(CFG$out_dir,
                            paste0("spot_zone_stage_", CFG$organ, ".rds"))
  attr(ni_all, "config") <- CFG
  attr(ni_all, "date") <- format(Sys.time())
  saveRDS(ni_all, out_rds_path)
  cat_log("  RDS: ", out_rds_path, " (",
          round(file.info(out_rds_path)$size / 1024^2, 1), " MB)\n")
}

# Save pivot summary
pivot_csv <- file.path(CFG$out_dir,
                      paste0("zone_stage_pivot_", CFG$organ, ".csv"))
write_csv(pivot, pivot_csv)
write_csv(pivot_long,
          file.path(CFG$out_dir,
                    paste0("zone_stage_per_sample_", CFG$organ, ".csv")))
cat_log("  pivot: ", pivot_csv, "\n")

cat_log("\n=== DONE: ", format(Sys.time()), " ===\n")
cat_log("Next: Rscript D2_step2_qc_visualize.R --organ ", CFG$organ, "\n\n")
