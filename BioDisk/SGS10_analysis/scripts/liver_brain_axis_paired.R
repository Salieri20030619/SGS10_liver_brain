#!/usr/bin/env Rscript
# ============================================================
# liver_brain_axis_paired.R  --  T5: 肝脑轴 paired correlation
# ============================================================
# 
# 目的: 验证 Task A v1.1 Section 5.1-5.2 的跨器官轴假说
#   X1 TGFb 轴  : liver host Tgfb1 <-> brain Microglia Smad3 
#   X2 甲状腺轴 : liver Hep_PP Ttr  <-> brain CPECs Ttr / Endothelial Thra
#
# 数据源:
#   1. RCTD dominant celltype: spatial_results/rctd/per_sample/<sample>_dominant_celltype.csv
#      含: barcode, dominant_celltype, dominant_score, second_celltype, second_score
#   2. SeekSpace Seurat obj:   spatial_results/<sample>_seekspace.rds
#      含: counts/data (gene × spot), zone (tumor/margin/host), species_class
#
# 配对:
#   M3  = Liver_3 <-> Brain_3      (tumor, paired)
#   M8  = Liver_8 <-> Brain_8      (tumor, paired)
#   M56 = gan56   <-> nao56        (tumor, paired)
#   Mbl = Bla_liver_1 <-> Bla_brain_1  (blank, pseudo-paired, 非同一只鼠)
#   Mcon = Con_liver_1 <-> Con_brain_1 (control, pseudo-paired, 非同一只鼠)
#
# 用法:
#   cd ~/SGS10_analysis
#   screen -S t5 -dm bash -c '
#     Rscript scripts/liver_brain_axis_paired.R 2>&1 | \
#       tee scripts/logs/liver_brain_axis.log
#   '
#
# 预计运行时间: 10-20 min (主要消耗在 5 liver + 5 brain rds 的 load)
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

select    <- dplyr::select
filter    <- dplyr::filter
mutate    <- dplyr::mutate
summarise <- dplyr::summarise
group_by  <- dplyr::group_by
bind_rows <- dplyr::bind_rows

# ============================================================
# 0. 路径与参数
# ============================================================
PROJECT_ROOT  <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
SPATIAL_DIR   <- file.path(PROJECT_ROOT, "spatial_results")
RCTD_DIR      <- file.path(SPATIAL_DIR, "rctd", "per_sample")
OUT_DIR       <- file.path(PROJECT_ROOT, "taskB", "liver_brain_axis2")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
setwd(PROJECT_ROOT)

cat("============================================================\n")
cat("T5: Liver-Brain Axis Paired Analysis (X1 TGFb + X2 Thyroid)\n")
cat(sprintf("Start: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("============================================================\n")

# ============================================================
# 1. Paired mice definition
# ============================================================
PAIRED_MICE <- list(
  list(id = "M3",   liver = "Liver_3",      brain = "Brain_3",      group = "tumor"),
  list(id = "M8",   liver = "Liver_8",      brain = "Brain_8",      group = "tumor"),
  list(id = "M56",  liver = "gan56",        brain = "nao56",        group = "tumor"),
  list(id = "Mbl",  liver = "Bla_liver_1",  brain = "Bla_brain_1",  group = "blank"),
  list(id = "Mcon", liver = "Con_liver_1",  brain = "Con_brain_1",  group = "control")
)

TUMOR_MICE <- c("M3", "M8", "M56")

# ============================================================
# 2. Gene panels (X1 + X2)
# ============================================================
# Task A v1.1 核心基因:
#   X1 TGFb: liver Tgfb1/Smad3/Tgfbr2 host; brain Smad3/Tgfbr2 Microglia
#   X2 Thyroid: liver Ttr Hep_PP; brain Ttr CPECs, Thra Endothelial

LIVER_GENES <- c("Tgfb1", "Smad3", "Tgfbr2", "Ttr", "Thrb", "Thra",
                 "Acta2", "Col1a1", "Serpine1", "F2r")  # 核心 panel
BRAIN_GENES <- c("Smad3", "Tgfbr2", "Tgfbr1", "Ttr", "Thra", "Cldn5",
                 "Apoe", "Lrp1", "Sox9", "Zbtb7a")  # 核心 panel

# Cell type target
HEP_PP_CELLS   <- c("Hep_Periportal", "Hepatocyte_PP", "Hepatocyte_Periportal", "Hep_PP")
HEP_PV_CELLS   <- c("Hep_Perivenous", "Hep_Pericentral", "Hepatocyte_PV",
                    "Hepatocyte_Pericentral", "Hep_PV")
MICROGLIA_CELLS <- c("Microglia_Homeostatic", "Microglia_homeostatic", "Microglia")
CPECS_CELLS    <- c("CPECs", "CPEC", "Choroid_Plexus_Epithelial")
ENDO_CELLS     <- c("Endothelial", "Endothelial_cells", "Endothelial_Brain",
                    "Endothelial_Vein", "Endothelial_Artery")

find_ct <- function(candidates, available) {
  intersect(candidates, available)[1]
}

# ============================================================
# 3. Per-sample reader: load rds + merge RCTD
# ============================================================
load_sample <- function(sample, gene_panel) {
  rds_path   <- file.path(SPATIAL_DIR, paste0(sample, "_seekspace.rds"))
  rctd_path  <- file.path(RCTD_DIR, paste0(sample, "_dominant_celltype.csv"))
  
  if (!file.exists(rds_path)) {
    cat(sprintf("  ⚠ Missing rds: %s\n", rds_path))
    return(NULL)
  }
  if (!file.exists(rctd_path)) {
    cat(sprintf("  ⚠ Missing RCTD csv: %s\n", rctd_path))
    return(NULL)
  }
  
  cat(sprintf("  Loading %s ... ", sample))
  t0 <- Sys.time()
  obj <- readRDS(rds_path)
  
  # Extract metadata: zone, species_class
  md <- obj@meta.data
  md$barcode <- rownames(md)
  
  # --------------------------------------------------------
  # CRITICAL: Dual-species PDX gene name resolution
  # --------------------------------------------------------
  # SGS10 samples use two reference types:
  #   (a) Dual-species (PDX tumor samples with human cells):
  #         gan56, Liver_3, Liver_8
  #       Mouse genes prefixed "mm10---Tgfb1"
  #       Human genes prefixed "GRCh38-TGFB1"
  #   (b) Mouse-only (baseline + all brain):
  #         Bla_liver_1, Con_liver_1, all brain
  #       Bare mouse symbols "Tgfb1"
  # 
  # resolve_gene tries, in order:
  #   1. exact bare name (mouse-only rds)
  #   2. mm10---<gene>  (dual-species, mouse side)
  # Human-side expression not used here (we only care about host response).
  # --------------------------------------------------------
  all_feats <- rownames(obj)
  
  resolve_gene <- function(g) {
    if (g %in% all_feats) return(g)
    mm_name <- paste0("mm10---", g)
    if (mm_name %in% all_feats) return(mm_name)
    NA_character_
  }
  
  gene_map       <- sapply(gene_panel, resolve_gene)
  available_genes <- gene_panel[!is.na(gene_map)]   # panel 里能找到的基因名 (bare)
  actual_names    <- gene_map[!is.na(gene_map)]      # 对应的 rds 实际 rowname
  missing_genes   <- gene_panel[is.na(gene_map)]
  
  if (length(missing_genes) > 0) {
    cat(sprintf("\n    ⚠ Missing in rds: %s",
                paste(missing_genes, collapse = ", ")))
  }
  prefixed <- available_genes[actual_names != available_genes]
  if (length(prefixed) > 0) {
    cat(sprintf("\n    (mm10--- prefix used for %d genes: %s)",
                length(prefixed), paste(prefixed, collapse = ", ")))
  }
  
  if (length(available_genes) > 0) {
    # 用 data slot (normalized) 做比较
    expr_mat <- GetAssayData(obj, assay = DefaultAssay(obj), layer = "data")
    gene_expr <- as.matrix(t(expr_mat[actual_names, , drop = FALSE]))
    # 用 bare name 作为列名, 下游代码统一按 bare name 访问
    colnames(gene_expr) <- available_genes
    gene_df <- as.data.frame(gene_expr)
    gene_df$barcode <- rownames(gene_df)
  } else {
    gene_df <- data.frame(barcode = rownames(md))
  }
  
  # Merge RCTD
  rctd_df <- read.csv(rctd_path, stringsAsFactors = FALSE)
  rctd_df$barcode <- gsub('"', '', rctd_df$barcode)
  
  # Merge
  merged <- md %>%
    left_join(rctd_df, by = "barcode") %>%
    left_join(gene_df, by = "barcode")
  
  rm(obj); gc(verbose = FALSE)
  t_elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("\n  done (%.1fs, %d spots, %d genes available)\n",
              t_elapsed, nrow(merged), length(available_genes)))
  
  list(
    data = merged,
    genes_available = available_genes,
    n_spots = nrow(merged)
  )
}

# ============================================================
# 4. Metric extractors
# ============================================================
# 对每个 sample 提取一组 metric (mean expression by cell type / zone)

extract_liver_metrics <- function(liver_data, sample_name) {
  df <- liver_data$data
  genes <- liver_data$genes_available
  
  # zone check
  if (!"zone" %in% colnames(df)) {
    cat(sprintf("    ⚠ %s: no 'zone' column, using species_class fallback\n", sample_name))
    if ("species_class" %in% colnames(df)) {
      df$zone <- ifelse(df$species_class == "human", "tumor",
                 ifelse(df$species_class == "mix", "margin", "host"))
    } else {
      df$zone <- "all"
    }
  }
  
  # Host zone
  host_df <- df %>% filter(zone == "host")
  n_host  <- nrow(host_df)
  n_total <- nrow(df)
  
  # Hep_PP dominant spots (within host)
  hep_pp_ct <- find_ct(HEP_PP_CELLS, unique(df$dominant_celltype))
  hep_pv_ct <- find_ct(HEP_PV_CELLS, unique(df$dominant_celltype))
  
  if (!is.null(hep_pp_ct) && !is.na(hep_pp_ct)) {
    hep_pp_df <- host_df %>% filter(dominant_celltype == hep_pp_ct)
  } else {
    hep_pp_df <- host_df[0, ]
  }
  
  # Safe mean function (returns NA if empty or all NA)
  smean <- function(x) {
    if (length(x) == 0) return(NA_real_)
    v <- x[!is.na(x)]
    if (length(v) == 0) return(NA_real_)
    mean(v)
  }
  
  spct_pos <- function(x) {
    if (length(x) == 0) return(NA_real_)
    v <- x[!is.na(x)]
    if (length(v) == 0) return(NA_real_)
    100 * sum(v > 0) / length(v)
  }
  
  gfetch <- function(d, g) if (g %in% colnames(d) && g %in% genes) d[[g]] else rep(NA_real_, nrow(d))
  
  list(
    sample = sample_name,
    n_total = n_total,
    n_host = n_host,
    hep_pp_ct = ifelse(is.null(hep_pp_ct) | is.na(hep_pp_ct), "NA", hep_pp_ct),
    hep_pp_pct = 100 * nrow(hep_pp_df) / n_total,
    hep_pv_pct = if (!is.null(hep_pv_ct) && !is.na(hep_pv_ct)) {
        100 * sum(df$dominant_celltype == hep_pv_ct, na.rm = TRUE) / n_total
      } else NA_real_,
    # ---- X1 TGFb (liver host zone) ----
    liver_host_Tgfb1_mean  = smean(gfetch(host_df, "Tgfb1")),
    liver_host_Tgfb1_pct   = spct_pos(gfetch(host_df, "Tgfb1")),
    liver_host_Smad3_mean  = smean(gfetch(host_df, "Smad3")),
    liver_host_Smad3_pct   = spct_pos(gfetch(host_df, "Smad3")),
    liver_host_Tgfbr2_mean = smean(gfetch(host_df, "Tgfbr2")),
    # ---- X2 Thyroid (liver Hep_PP + global) ----
    liver_hep_pp_Ttr  = smean(gfetch(hep_pp_df, "Ttr")),
    liver_host_Ttr    = smean(gfetch(host_df, "Ttr")),
    liver_host_Thrb   = smean(gfetch(host_df, "Thrb")),
    liver_host_Thra   = smean(gfetch(host_df, "Thra")),
    # ---- L1 HSC 激活 (核心 marker) ----
    liver_host_Acta2  = smean(gfetch(host_df, "Acta2")),
    liver_host_Col1a1 = smean(gfetch(host_df, "Col1a1")),
    liver_host_Serpine1 = smean(gfetch(host_df, "Serpine1")),
    liver_host_F2r    = smean(gfetch(host_df, "F2r"))
  )
}

extract_brain_metrics <- function(brain_data, sample_name) {
  df <- brain_data$data
  genes <- brain_data$genes_available
  
  n_total <- nrow(df)
  
  # Microglia dominant spots
  mg_ct    <- find_ct(MICROGLIA_CELLS, unique(df$dominant_celltype))
  cpecs_ct <- find_ct(CPECS_CELLS,     unique(df$dominant_celltype))
  endo_ct  <- find_ct(ENDO_CELLS,      unique(df$dominant_celltype))
  
  mg_df    <- if (!is.null(mg_ct)    && !is.na(mg_ct))    df %>% filter(dominant_celltype == mg_ct)    else df[0, ]
  cpecs_df <- if (!is.null(cpecs_ct) && !is.na(cpecs_ct)) df %>% filter(dominant_celltype == cpecs_ct) else df[0, ]
  endo_df  <- if (!is.null(endo_ct)  && !is.na(endo_ct))  df %>% filter(dominant_celltype == endo_ct)  else df[0, ]
  
  smean <- function(x) {
    if (length(x) == 0) return(NA_real_)
    v <- x[!is.na(x)]
    if (length(v) == 0) return(NA_real_)
    mean(v)
  }
  gfetch <- function(d, g) if (g %in% colnames(d) && g %in% genes) d[[g]] else rep(NA_real_, nrow(d))
  
  list(
    sample = sample_name,
    n_total = n_total,
    mg_ct = ifelse(is.null(mg_ct) | is.na(mg_ct), "NA", mg_ct),
    cpecs_ct = ifelse(is.null(cpecs_ct) | is.na(cpecs_ct), "NA", cpecs_ct),
    endo_ct = ifelse(is.null(endo_ct) | is.na(endo_ct), "NA", endo_ct),
    # Cell type 占比
    brain_microglia_pct = 100 * nrow(mg_df) / n_total,
    brain_cpecs_pct     = 100 * nrow(cpecs_df) / n_total,
    brain_endo_pct      = 100 * nrow(endo_df) / n_total,
    # ---- X1 TGFb (brain Microglia) ----
    brain_mg_Smad3_mean    = smean(gfetch(mg_df, "Smad3")),
    brain_mg_Tgfbr2_mean   = smean(gfetch(mg_df, "Tgfbr2")),
    brain_mg_Tgfbr1_mean   = smean(gfetch(mg_df, "Tgfbr1")),
    brain_mg_Apoe_mean     = smean(gfetch(mg_df, "Apoe")),
    brain_mg_Lrp1_mean     = smean(gfetch(mg_df, "Lrp1")),
    brain_global_Smad3     = smean(gfetch(df, "Smad3")),
    # ---- X2 Thyroid (brain CPECs + Endothelial) ----
    brain_cpecs_Ttr        = smean(gfetch(cpecs_df, "Ttr")),
    brain_endo_Thra        = smean(gfetch(endo_df, "Thra")),
    brain_global_Ttr       = smean(gfetch(df, "Ttr")),
    brain_global_Thra      = smean(gfetch(df, "Thra")),
    # ---- BBB signature ----
    brain_endo_Cldn5       = smean(gfetch(endo_df, "Cldn5")),
    # ---- Pan-brain Zbtb7a ----
    brain_global_Zbtb7a    = smean(gfetch(df, "Zbtb7a"))
  )
}

# ============================================================
# 5. Main loop: per-mouse metric extraction
# ============================================================
cat("\n### Loading 10 samples and extracting metrics ###\n\n")

all_metrics <- list()

for (mouse in PAIRED_MICE) {
  cat(sprintf("\n>>> %s: %s (liver) + %s (brain) [%s]\n",
              mouse$id, mouse$liver, mouse$brain, mouse$group))
  
  liver_data <- load_sample(mouse$liver, LIVER_GENES)
  brain_data <- load_sample(mouse$brain, BRAIN_GENES)
  
  if (is.null(liver_data) || is.null(brain_data)) {
    cat(sprintf("  ⚠ Skipping %s due to missing data\n", mouse$id))
    next
  }
  
  lm <- extract_liver_metrics(liver_data, mouse$liver)
  bm <- extract_brain_metrics(brain_data, mouse$brain)
  
  # Combine
  row <- c(
    list(
      mouse_id = mouse$id,
      group    = mouse$group,
      liver_sample = mouse$liver,
      brain_sample = mouse$brain
    ),
    lm, bm
  )
  # 去掉重复 sample 列
  row$sample <- NULL
  
  all_metrics[[mouse$id]] <- as.data.frame(row, stringsAsFactors = FALSE)
  
  cat(sprintf("  [%s] Liver host Tgfb1 mean=%.3f, Smad3=%.3f | Brain MG Smad3=%.3f, CPECs Ttr=%.3f\n",
              mouse$id, lm$liver_host_Tgfb1_mean, lm$liver_host_Smad3_mean,
              bm$brain_mg_Smad3_mean, bm$brain_cpecs_Ttr))
}

if (length(all_metrics) == 0) {
  stop("No samples loaded successfully. Aborting.")
}

# Align columns across rows (because different rows may have different NA patterns)
all_cols <- unique(unlist(lapply(all_metrics, names)))
for (id in names(all_metrics)) {
  missing_cols <- setdiff(all_cols, names(all_metrics[[id]]))
  for (c in missing_cols) all_metrics[[id]][[c]] <- NA
  all_metrics[[id]] <- all_metrics[[id]][, all_cols, drop = FALSE]
}

metrics_df <- bind_rows(all_metrics)
cat("\n### Combined metrics (all 5 mice) ###\n")
print(metrics_df[, c("mouse_id", "group",
                     "liver_host_Tgfb1_mean", "liver_host_Smad3_mean",
                     "liver_hep_pp_Ttr",
                     "brain_mg_Smad3_mean", "brain_cpecs_Ttr", "brain_endo_Thra")])

# Save full metrics
write.csv(metrics_df, file.path(OUT_DIR, "01_paired_full_metrics.csv"), row.names = FALSE)
cat(sprintf("\n  -> 01_paired_full_metrics.csv (%d mice x %d cols)\n",
            nrow(metrics_df), ncol(metrics_df)))

# ============================================================
# 6. Tumor vs Baseline group comparison
# ============================================================
cat("\n### Tumor vs Baseline group comparison ###\n")

baseline_ids <- c("Mbl", "Mcon")
metrics_df$group_simple <- ifelse(metrics_df$mouse_id %in% TUMOR_MICE,
                                  "tumor", "baseline")

group_cmp <- metrics_df %>%
  group_by(group_simple) %>%
  summarise(
    n = n(),
    # X1
    mean_liver_Tgfb1_mean = mean(liver_host_Tgfb1_mean, na.rm = TRUE),
    mean_liver_Tgfb1_pct  = mean(liver_host_Tgfb1_pct, na.rm = TRUE),
    mean_liver_Smad3      = mean(liver_host_Smad3_mean, na.rm = TRUE),
    mean_brain_mg_Smad3   = mean(brain_mg_Smad3_mean, na.rm = TRUE),
    mean_brain_mg_Tgfbr2  = mean(brain_mg_Tgfbr2_mean, na.rm = TRUE),
    # X2
    mean_liver_hep_pp_Ttr = mean(liver_hep_pp_Ttr, na.rm = TRUE),
    mean_brain_cpecs_Ttr  = mean(brain_cpecs_Ttr, na.rm = TRUE),
    mean_brain_endo_Thra  = mean(brain_endo_Thra, na.rm = TRUE),
    # L1
    mean_liver_Acta2      = mean(liver_host_Acta2, na.rm = TRUE),
    mean_liver_Serpine1   = mean(liver_host_Serpine1, na.rm = TRUE),
    .groups = "drop"
  )

print(group_cmp)
write.csv(group_cmp, file.path(OUT_DIR, "02_tumor_vs_baseline.csv"), row.names = FALSE)

# ============================================================
# 7. Paired correlation analyses
# ============================================================
cat("\n### Paired correlation (X1 + X2) ###\n")

tumor_df <- metrics_df %>% filter(group_simple == "tumor")

if (nrow(tumor_df) < 2) {
  cat("  ⚠ Not enough tumor mice (< 2), skipping correlation\n")
} else {
  cor_pairs <- list(
    # X1 main
    list(name = "X1_liver_Tgfb1_~_brain_MG_Smad3",
         x = "liver_host_Tgfb1_mean", y = "brain_mg_Smad3_mean"),
    list(name = "X1_liver_Smad3_~_brain_MG_Smad3",
         x = "liver_host_Smad3_mean", y = "brain_mg_Smad3_mean"),
    list(name = "X1_liver_Tgfb1_~_brain_MG_Tgfbr2",
         x = "liver_host_Tgfb1_mean", y = "brain_mg_Tgfbr2_mean"),
    # X2 main
    list(name = "X2_liver_Ttr_HepPP_~_brain_CPECs_Ttr",
         x = "liver_hep_pp_Ttr", y = "brain_cpecs_Ttr"),
    list(name = "X2_liver_Ttr_HepPP_~_brain_endo_Thra",
         x = "liver_hep_pp_Ttr", y = "brain_endo_Thra"),
    list(name = "X2_liver_Ttr_global_~_brain_global_Ttr",
         x = "liver_host_Ttr", y = "brain_global_Ttr"),
    # Cross-axis sanity
    list(name = "Cross_liver_Tgfb1_~_brain_endo_Thra",
         x = "liver_host_Tgfb1_mean", y = "brain_endo_Thra"),
    list(name = "Cross_liver_Ttr_~_brain_MG_Smad3",
         x = "liver_hep_pp_Ttr", y = "brain_mg_Smad3_mean")
  )
  
  cor_results <- list()
  for (p in cor_pairs) {
    xv <- tumor_df[[p$x]]
    yv <- tumor_df[[p$y]]
    valid <- !is.na(xv) & !is.na(yv)
    
    if (sum(valid) < 2) {
      cor_results[[p$name]] <- data.frame(
        pair = p$name, n = sum(valid), r = NA, p = NA,
        note = "insufficient data"
      )
      next
    }
    
    ct <- tryCatch(
      cor.test(xv[valid], yv[valid], method = "pearson"),
      error = function(e) NULL
    )
    if (is.null(ct)) {
      cor_results[[p$name]] <- data.frame(
        pair = p$name, n = sum(valid), r = NA, p = NA,
        note = "cor.test failed"
      )
      next
    }
    cor_results[[p$name]] <- data.frame(
      pair = p$name,
      n = sum(valid),
      r = as.numeric(ct$estimate),
      p = ct$p.value,
      note = sprintf("n=%d paired tumor mice", sum(valid))
    )
    cat(sprintf("  %-45s  r = %+.3f, p = %.3f (n=%d)\n",
                p$name, ct$estimate, ct$p.value, sum(valid)))
  }
  
  cor_df <- bind_rows(cor_results)
  write.csv(cor_df, file.path(OUT_DIR, "03_paired_correlations.csv"), row.names = FALSE)
  cat(sprintf("\n  -> 03_paired_correlations.csv\n"))
}

# ============================================================
# 8. Scatter plots: X1 + X2 main axes
# ============================================================
cat("\n### Scatter plots ###\n")

make_scatter <- function(d, x_col, y_col, title_text, subtitle_text = NULL,
                         xlab = x_col, ylab = y_col) {
  plot_df <- d %>%
    select(mouse_id, group_simple, x = all_of(x_col), y = all_of(y_col)) %>%
    filter(!is.na(x) & !is.na(y))
  
  if (nrow(plot_df) == 0) return(NULL)
  
  tumor_only <- plot_df %>% filter(group_simple == "tumor")
  r_text <- ""
  if (nrow(tumor_only) >= 2) {
    ct <- tryCatch(cor.test(tumor_only$x, tumor_only$y),
                   error = function(e) NULL)
    if (!is.null(ct)) {
      r_text <- sprintf("tumor n=%d: r = %+.2f, p = %.3f",
                        nrow(tumor_only), ct$estimate, ct$p.value)
    }
  }
  sub <- if (is.null(subtitle_text)) r_text else paste(subtitle_text, r_text, sep = "\n")
  
  p <- ggplot(plot_df, aes(x = x, y = y, color = group_simple, shape = group_simple)) +
    geom_point(size = 5, alpha = 0.85) +
    geom_text_repel(aes(label = mouse_id), size = 3.8, show.legend = FALSE) +
    { if (nrow(tumor_only) >= 2) 
        geom_smooth(data = tumor_only, method = "lm", se = FALSE,
                    color = "#d62728", linetype = "dashed", size = 0.7)
    } +
    scale_color_manual(values = c(tumor = "#d62728", baseline = "#1f77b4")) +
    scale_shape_manual(values = c(tumor = 17, baseline = 16)) +
    labs(title = title_text,
         subtitle = sub,
         x = xlab, y = ylab,
         color = "group", shape = "group") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 9, color = "#555555"))
  p
}

# X1 main
p1 <- make_scatter(metrics_df,
                   "liver_host_Tgfb1_mean", "brain_mg_Smad3_mean",
                   title_text = "X1 TGFb axis (main)",
                   xlab = "Liver host Tgfb1 mean expression",
                   ylab = "Brain Microglia Smad3 mean expression")

# X1 secondary: liver Smad3 (autocrine) <-> brain MG Smad3
p2 <- make_scatter(metrics_df,
                   "liver_host_Smad3_mean", "brain_mg_Smad3_mean",
                   title_text = "X1 TGFb axis (secondary: host Smad3 parallel)",
                   xlab = "Liver host Smad3 mean",
                   ylab = "Brain Microglia Smad3 mean")

# X2 main
p3 <- make_scatter(metrics_df,
                   "liver_hep_pp_Ttr", "brain_cpecs_Ttr",
                   title_text = "X2 Thyroid axis (main)",
                   xlab = "Liver Hep_PP Ttr mean",
                   ylab = "Brain CPECs Ttr mean")

# X2 secondary: liver Ttr -> brain endothelial Thra
p4 <- make_scatter(metrics_df,
                   "liver_hep_pp_Ttr", "brain_endo_Thra",
                   title_text = "X2 Thyroid axis (secondary: liver Ttr -> brain Endo Thra)",
                   xlab = "Liver Hep_PP Ttr mean",
                   ylab = "Brain Endothelial Thra mean")

panels <- list(p1, p2, p3, p4)
panels <- panels[!sapply(panels, is.null)]

if (length(panels) > 0) {
  combined <- wrap_plots(panels, ncol = 2) +
    plot_annotation(
      title = "Liver-Brain Axis Paired Analysis (X1 TGFb + X2 Thyroid)",
      subtitle = sprintf("3 paired tumor mice (M3/M8/M56) + 2 pseudo-paired baseline (Mbl/Mcon)  |  %s",
                        format(Sys.time(), "%Y-%m-%d")),
      theme = theme(plot.title = element_text(face = "bold", size = 13))
    )
  
  ggsave(file.path(OUT_DIR, "04_scatter_x1x2_combined.pdf"),
         combined, width = 13, height = 10)
  ggsave(file.path(OUT_DIR, "04_scatter_x1x2_combined.png"),
         combined, width = 13, height = 10, dpi = 200)
  
  # Also save individual plots
  if (!is.null(p1)) ggsave(file.path(OUT_DIR, "04a_x1_main_scatter.pdf"),
                          p1, width = 7, height = 5.5)
  if (!is.null(p3)) ggsave(file.path(OUT_DIR, "04b_x2_main_scatter.pdf"),
                          p3, width = 7, height = 5.5)
  
  cat(sprintf("  -> 04_scatter_x1x2_combined.{pdf,png}\n"))
  cat(sprintf("  -> 04a_x1_main_scatter.pdf, 04b_x2_main_scatter.pdf\n"))
}

# ============================================================
# 9. Forest plot: 8 key metrics tumor vs baseline
# ============================================================
cat("\n### Forest plot: tumor vs baseline ###\n")

forest_metrics <- c(
  "liver_host_Tgfb1_mean", "liver_host_Smad3_mean", "liver_hep_pp_Ttr",
  "liver_host_Acta2", "liver_host_Serpine1",
  "brain_mg_Smad3_mean", "brain_cpecs_Ttr", "brain_endo_Thra",
  "brain_endo_Cldn5"
)

forest_df <- data.frame(metric = character(), 
                        tumor_mean = numeric(), baseline_mean = numeric(),
                        fc = numeric(), delta = numeric(),
                        stringsAsFactors = FALSE)

for (m in forest_metrics) {
  if (!m %in% colnames(metrics_df)) next
  tumor_vals <- metrics_df[metrics_df$group_simple == "tumor", m]
  base_vals  <- metrics_df[metrics_df$group_simple == "baseline", m]
  tm <- mean(tumor_vals, na.rm = TRUE)
  bm <- mean(base_vals,  na.rm = TRUE)
  if (is.nan(tm) || is.nan(bm) || is.na(tm) || is.na(bm)) next
  fc <- if (bm > 1e-10) tm / bm else NA_real_
  forest_df <- rbind(forest_df, data.frame(
    metric = m, tumor_mean = tm, baseline_mean = bm,
    fc = fc, delta = tm - bm
  ))
}

write.csv(forest_df, file.path(OUT_DIR, "05_tumor_vs_baseline_forest.csv"),
          row.names = FALSE)

if (nrow(forest_df) > 0) {
  # Bar plot: log2 FC 
  forest_df$log2fc <- log2(pmax(forest_df$fc, 0.01))
  forest_df$metric_short <- gsub("(liver|brain)_(host|hep_pp|mg|cpecs|endo|global)_", "", forest_df$metric)
  forest_df$side <- ifelse(grepl("^liver", forest_df$metric), "Liver", "Brain")
  
  fp <- ggplot(forest_df, aes(x = reorder(metric_short, log2fc),
                              y = log2fc, fill = side)) +
    geom_col(alpha = 0.85) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    geom_text(aes(label = sprintf("%.1f×", fc)),
              hjust = ifelse(forest_df$log2fc >= 0, -0.1, 1.1), size = 3) +
    coord_flip() +
    scale_fill_manual(values = c(Liver = "#e58e25", Brain = "#5b8dd6")) +
    labs(title = "Tumor vs Baseline fold-change (key metrics)",
         subtitle = "log2 FC of mean expression, 3 tumor vs 2 baseline mice",
         x = NULL, y = "log2 (tumor / baseline)") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "bottom")
  
  ggsave(file.path(OUT_DIR, "05_tumor_vs_baseline_forest.pdf"),
         fp, width = 8, height = 6)
  ggsave(file.path(OUT_DIR, "05_tumor_vs_baseline_forest.png"),
         fp, width = 8, height = 6, dpi = 200)
  cat(sprintf("  -> 05_tumor_vs_baseline_forest.{pdf,png}\n"))
}

# ============================================================
# 10. Summary markdown
# ============================================================
cat("\n### Writing summary.md ###\n")

md <- c(
  "# T5 Liver-Brain Axis Paired Analysis Summary",
  "",
  sprintf("**Generated:** %s", format(Sys.time(), "%Y-%m-%d %H:%M")),
  "",
  "## Hypothesis (Task A v1.1 Section 5.1-5.2)",
  "",
  "- **X1 TGFb axis**: Liver host Tgfb1/Smad3 (HSC + hepatocyte autocrine amplification) → circulating TGFb → brain Microglia Smad3/Tgfbr2 activation",
  "- **X2 Thyroid axis**: Liver Hep_PP Ttr → blood TTR/T4 → brain CPECs Ttr ↑ (CSF compensation) + Endothelial Thra ↑ (BBB receptor sensitization)",
  "",
  "## Mouse pairing",
  "",
  "| Mouse | Liver | Brain | Group | Note |",
  "|-------|-------|-------|-------|------|",
  "| M3   | Liver_3     | Brain_3     | tumor   | same mouse |",
  "| M8   | Liver_8     | Brain_8     | tumor   | same mouse |",
  "| M56  | gan56       | nao56       | tumor   | same mouse |",
  "| Mbl  | Bla_liver_1 | Bla_brain_1 | blank   | pseudo-paired (NOT same mouse) |",
  "| Mcon | Con_liver_1 | Con_brain_1 | control | pseudo-paired (NOT same mouse) |",
  "",
  "## Per-mouse key metrics",
  "",
  "| Mouse | Group | Liver Tgfb1 host | Liver Hep_PP Ttr | Brain MG Smad3 | Brain CPECs Ttr | Brain Endo Thra |",
  "|-------|-------|------------------:|------------------:|---------------:|----------------:|----------------:|"
)

for (i in seq_len(nrow(metrics_df))) {
  r <- metrics_df[i, ]
  fmt <- function(x) if (is.na(x)) "NA" else sprintf("%.3f", x)
  md <- c(md,
    sprintf("| %s | %s | %s | %s | %s | %s | %s |",
            r$mouse_id, r$group,
            fmt(r$liver_host_Tgfb1_mean), fmt(r$liver_hep_pp_Ttr),
            fmt(r$brain_mg_Smad3_mean), fmt(r$brain_cpecs_Ttr),
            fmt(r$brain_endo_Thra)))
}

# Correlation results section
if (exists("cor_df") && nrow(cor_df) > 0) {
  md <- c(md, "", "## Paired correlations (tumor only, n = 3)", "",
          "| Pair | n | r | p | Interpretation |",
          "|------|---|---|---|---------------|")
  for (i in seq_len(nrow(cor_df))) {
    r <- cor_df[i, ]
    interp <- if (is.na(r$r)) "insufficient data" else
              if (r$r > 0.5) "strong positive (supports axis)" else
              if (r$r > 0.2) "weak positive" else
              if (r$r < -0.5) "strong negative (contradicts axis)" else
              if (r$r < -0.2) "weak negative" else "~ no correlation"
    md <- c(md, sprintf("| %s | %d | %s | %s | %s |",
                        r$pair, r$n,
                        if (is.na(r$r)) "NA" else sprintf("%+.3f", r$r),
                        if (is.na(r$p)) "NA" else sprintf("%.3f", r$p),
                        interp))
  }
}

md <- c(md, "",
  "## Limitations",
  "",
  "1. **n=3 paired tumor mice** - statistical power limited; report direction + effect size, not p-values",
  "2. **Baseline n=2 pseudo-paired** - Bla/Con are NOT same mouse liver-brain; used as reference only",
  "3. **Spatial-level, not per-cell** - mean expression per sample, not single-cell resolution",
  "",
  "## Output files",
  "",
  "- `01_paired_full_metrics.csv` - all metrics per mouse (~30 columns)",
  "- `02_tumor_vs_baseline.csv` - group comparison",
  "- `03_paired_correlations.csv` - 8 axis correlation stats",
  "- `04_scatter_x1x2_combined.{pdf,png}` - 4-panel main figure",
  "- `04a_x1_main_scatter.pdf` - X1 single panel",
  "- `04b_x2_main_scatter.pdf` - X2 single panel",
  "- `05_tumor_vs_baseline_forest.{pdf,png}` - log2FC bar plot"
)

writeLines(md, file.path(OUT_DIR, "summary.md"))
cat(sprintf("  -> summary.md (%d lines)\n", length(md)))

# ============================================================
# Done
# ============================================================
cat("\n============================================================\n")
cat("T5 liver-brain axis paired analysis COMPLETE\n")
cat(sprintf("Finish: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(sprintf("Output: %s\n", OUT_DIR))
cat(sprintf("Files:  %d\n", length(list.files(OUT_DIR))))
cat("============================================================\n")
