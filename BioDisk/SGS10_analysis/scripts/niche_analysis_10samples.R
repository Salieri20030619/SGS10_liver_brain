#!/usr/bin/env Rscript
# ============================================================
# niche_analysis_10samples.R  --  T1 + T2 + T3
# ============================================================
# 
# T1 核心: 找 tumor-induced niche
#   对每个 anchor cell type (如 HSC_Activated_ECM, Microglia_Homeostatic),
#   算它的 spatial 邻居 cell type 富集 (用 k-NN on spatial coords + RCTD)
# 
# T2 扩展: 免疫-内皮-神经三角 (客户原话命中)
#   特别关注 Brain 里: Microglia ↔ Endothelial ↔ Neuron 三者 co-localization
# 
# T3 脑端 4-axis triangle:
#   Immune-Endo-Neuron-Glia 四类分组后看 tumor vs baseline
# 
# 输入:
#   - spatial_results/<sample>_seekspace.rds  (coords)
#   - spatial_results/rctd/per_sample/<sample>_dominant_celltype.csv
# 
# 输出:
#   taskB/niche_v2_10samples/
#     liver_niche_<anchor_ct>_neighbor_enrichment.csv     每个 anchor 的邻居富集
#     liver_niche_top_tumor_induced.csv                   tumor vs baseline 最大 shift 的 niche
#     liver_niche_<anchor>_heatmap.pdf                    可视化
#     brain_ (同上) + brain_immune_endo_neuron_triangle.{csv,pdf}
#     summary.md
# 
# 依赖:
#   FNN (k-NN)
#   dplyr, tidyr, ggplot2, pheatmap, RColorBrewer
# 
# 预计: 20-30 min
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(FNN)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
})

select    <- dplyr::select
filter    <- dplyr::filter
mutate    <- dplyr::mutate
summarise <- dplyr::summarise
group_by  <- dplyr::group_by
bind_rows <- dplyr::bind_rows

PROJECT_ROOT  <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
SPATIAL_DIR   <- file.path(PROJECT_ROOT, "spatial_results")
RCTD_DIR      <- file.path(SPATIAL_DIR, "rctd", "per_sample")
OUT_DIR       <- file.path(PROJECT_ROOT, "taskB", "niche_v2_10samples")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

set.seed(42)
K_NEIGHBORS <- 10   # k-NN

# ============================================================
# Samples per organ
# ============================================================
SAMPLES <- list(
  liver = list(
    tumor    = c("Liver_3", "Liver_8", "gan56"),
    baseline = c("Bla_liver_1", "Con_liver_1")
  ),
  brain = list(
    tumor    = c("Brain_3", "Brain_8", "nao56"),
    baseline = c("Bla_brain_1", "Con_brain_1")
  )
)

# ============================================================
# Anchor cell types per organ (for niche analysis)
# ============================================================
LIVER_ANCHORS <- c("HSC_Activated_ECM", "Hep_Periportal", "Hep_Perivenous",
                   "Kupffer", "Macrophage_M1", "Endothelial", "LSEC")
BRAIN_ANCHORS <- c("Microglia_Homeostatic", "Astrocyte_Resting", "CPECs",
                   "Endothelial", "Pericytes", "Excitatory_L23", "Hpc_CA1",
                   "Oligo_Mature")

# ============================================================
# Brain 4-axis triangle: immune / endo / neuron / glia
# ============================================================
BRAIN_AXES <- list(
  immune  = c("Microglia_Homeostatic", "Microglia_other"),
  endo    = c("Endothelial", "Pericytes", "CPECs", "VLMC", "Ependymal"),
  neuron  = c("Excitatory_L23", "Excitatory_L45", "Excitatory_L6",
              "Inhibitory_Pvalb", "Inhibitory_Sst", "Inhibitory_Vip",
              "Inhibitory_Meis2", "Hpc_CA1", "Hpc_CA3", "Hpc_DG", "Neurons"),
  glia    = c("Astrocyte_Resting", "Astrocyte_Gfap", "Oligo_Mature", "OPCs")
)

cat("============================================================\n")
cat("T1+T2+T3: Tumor-induced niche analysis\n")
cat(sprintf("Start: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("============================================================\n")

# ============================================================
# Load 1 sample: coords + dominant_celltype
# ============================================================
load_sample_for_niche <- function(sample) {
  rds_path  <- file.path(SPATIAL_DIR, paste0(sample, "_seekspace.rds"))
  rctd_path <- file.path(RCTD_DIR, paste0(sample, "_dominant_celltype.csv"))
  
  if (!file.exists(rds_path) || !file.exists(rctd_path)) {
    cat(sprintf("  ⚠ Missing data for %s\n", sample))
    return(NULL)
  }
  
  cat(sprintf("  Loading %s ... ", sample))
  t0 <- Sys.time()
  obj <- readRDS(rds_path)
  
  # Extract spatial coordinates
  # SeekSpace stores spatial coords in reductions$spatial (2D embedding)
  # Column names: spatial_1 (x), spatial_2 (y)
  coords <- tryCatch({
    # Primary: SeekSpace convention - reductions$spatial
    if ("spatial" %in% names(obj@reductions)) {
      emb <- Embeddings(obj, reduction = "spatial")
      if (ncol(emb) >= 2) {
        data.frame(barcode = rownames(emb),
                   x = emb[, 1],
                   y = emb[, 2])
      } else NULL
    } else if ("image_x" %in% colnames(obj@meta.data) &&
               "image_y" %in% colnames(obj@meta.data)) {
      # Fallback: meta.data columns
      data.frame(barcode = rownames(obj@meta.data),
                 x = obj@meta.data$image_x,
                 y = obj@meta.data$image_y)
    } else if (length(obj@images) > 0) {
      # Fallback: images slot
      img_name <- names(obj@images)[1]
      img_coords <- GetTissueCoordinates(obj, image = img_name)
      cc <- colnames(img_coords)
      if ("imagerow" %in% cc && "imagecol" %in% cc) {
        data.frame(barcode = rownames(img_coords),
                   x = img_coords$imagecol,
                   y = img_coords$imagerow)
      } else {
        num_cols <- sapply(img_coords, is.numeric)
        cc_num <- cc[num_cols]
        data.frame(barcode = rownames(img_coords),
                   x = img_coords[[cc_num[1]]],
                   y = img_coords[[cc_num[2]]])
      }
    } else {
      NULL
    }
  }, error = function(e) {
    cat(sprintf("\n    coord error: %s\n", conditionMessage(e))); NULL
  })
  
  if (is.null(coords) || nrow(coords) == 0) {
    cat(sprintf("  ⚠ No spatial coords found\n"))
    rm(obj); gc(verbose = FALSE)
    return(NULL)
  }
  
  # zone info
  md <- obj@meta.data
  md$barcode <- rownames(md)
  zone_col <- if ("zone" %in% colnames(md)) md$zone else rep("all", nrow(md))
  
  rctd_df <- read.csv(rctd_path, stringsAsFactors = FALSE)
  rctd_df$barcode <- gsub('"', '', rctd_df$barcode)
  
  merged <- coords %>%
    left_join(data.frame(barcode = md$barcode, zone = zone_col),
              by = "barcode") %>%
    left_join(rctd_df[, c("barcode", "dominant_celltype")],
              by = "barcode") %>%
    filter(!is.na(dominant_celltype) & !is.na(x) & !is.na(y))
  
  rm(obj); gc(verbose = FALSE)
  cat(sprintf("done (%.1fs, %d spots)\n",
              as.numeric(Sys.time()-t0), nrow(merged)))
  merged
}

# ============================================================
# Compute k-NN neighbor enrichment for one sample
# ============================================================
compute_niche <- function(df, k = K_NEIGHBORS) {
  # For each spot, find k nearest neighbors (spatial)
  # Return: long df with anchor_ct, neighbor_ct, pair_count
  coords <- as.matrix(df[, c("x", "y")])
  nn <- get.knn(coords, k = k)
  
  # For each spot, its k neighbors' cell types
  anchor_ct <- df$dominant_celltype[rep(seq_len(nrow(df)), each = k)]
  neighbor_idx <- as.vector(t(nn$nn.index))
  neighbor_ct <- df$dominant_celltype[neighbor_idx]
  
  pair_df <- data.frame(
    anchor = anchor_ct,
    neighbor = neighbor_ct,
    stringsAsFactors = FALSE
  ) %>%
    filter(!is.na(anchor) & !is.na(neighbor))
  
  # Count pairs
  pair_counts <- pair_df %>%
    group_by(anchor, neighbor) %>%
    summarise(n = n(), .groups = "drop")
  
  # Normalize: fraction of neighbors per anchor type
  pair_counts <- pair_counts %>%
    group_by(anchor) %>%
    mutate(pct = 100 * n / sum(n)) %>%
    ungroup()
  
  pair_counts
}

# ============================================================
# Run niche for one organ
# ============================================================
run_organ_niche <- function(organ, anchors) {
  cat(sprintf("\n\n### === %s ===\n", toupper(organ)))
  
  samples       <- c(SAMPLES[[organ]]$tumor, SAMPLES[[organ]]$baseline)
  tumor_samples <- SAMPLES[[organ]]$tumor
  base_samples  <- SAMPLES[[organ]]$baseline
  
  # --------------------------------------------------------
  # Load all samples, compute per-sample niche
  # --------------------------------------------------------
  cat(sprintf("\n[1/4] Compute per-sample niche (k=%d)\n", K_NEIGHBORS))
  per_sample_niche <- list()
  per_sample_ct_pct <- list()  # per-sample ct proportion (for background)
  
  for (s in samples) {
    df <- load_sample_for_niche(s)
    if (is.null(df)) next
    
    niche <- compute_niche(df, k = K_NEIGHBORS)
    niche$sample <- s
    niche$group <- ifelse(s %in% tumor_samples, "tumor", "baseline")
    per_sample_niche[[s]] <- niche
    
    # ct composition
    ct_pct <- df %>%
      group_by(dominant_celltype) %>%
      summarise(pct = 100 * n() / nrow(df), .groups = "drop") %>%
      mutate(sample = s, group = niche$group[1])
    per_sample_ct_pct[[s]] <- ct_pct
  }
  
  if (length(per_sample_niche) == 0) {
    cat(sprintf("  ⚠ No samples loaded for %s\n", organ))
    return(NULL)
  }
  
  all_niche <- bind_rows(per_sample_niche)
  all_ct_pct <- bind_rows(per_sample_ct_pct)
  
  write.csv(all_niche, file.path(OUT_DIR, paste0(organ, "_all_niche_pairs.csv")),
            row.names = FALSE)
  
  # --------------------------------------------------------
  # For each anchor cell type, group-level tumor vs baseline
  # --------------------------------------------------------
  cat(sprintf("\n[2/4] Group-level neighbor enrichment per anchor\n"))
  all_anchors_in_data <- unique(all_niche$anchor)
  anchors_present <- intersect(anchors, all_anchors_in_data)
  if (length(anchors_present) == 0) {
    cat(sprintf("  ⚠ None of the requested anchors found in %s data\n", organ))
    cat(sprintf("    Available: %s\n",
                paste(head(all_anchors_in_data, 10), collapse = ", ")))
    anchors_present <- head(all_anchors_in_data[!is.na(all_anchors_in_data)], 5)
  }
  
  cat(sprintf("  Running on anchors: %s\n",
              paste(anchors_present, collapse = ", ")))
  
  anchor_enrichments <- list()
  for (anc in anchors_present) {
    sub <- all_niche %>% filter(anchor == anc)
    if (nrow(sub) == 0) next
    
    group_niche <- sub %>%
      group_by(group, neighbor) %>%
      summarise(mean_pct = mean(pct), .groups = "drop") %>%
      pivot_wider(names_from = group, values_from = mean_pct,
                  values_fill = 0) %>%
      mutate(delta = tumor - baseline,
             fc = (tumor + 0.5) / (baseline + 0.5)) %>%
      arrange(desc(abs(delta)))
    
    group_niche$anchor <- anc
    anchor_enrichments[[anc]] <- group_niche
    
    write.csv(group_niche,
              file.path(OUT_DIR, sprintf("%s_niche_%s_enrichment.csv", organ,
                                          gsub("[/ ]", "_", anc))),
              row.names = FALSE)
  }
  
  # --------------------------------------------------------
  # Top tumor-induced niches overall
  # --------------------------------------------------------
  cat(sprintf("\n[3/4] Top tumor-induced niches (across all anchors)\n"))
  all_enrichments <- bind_rows(anchor_enrichments)
  
  top_tumor_induced <- all_enrichments %>%
    filter(delta > 0 & tumor > 1) %>%  # 至少占 1%
    arrange(desc(delta)) %>%
    head(20)
  
  write.csv(top_tumor_induced,
            file.path(OUT_DIR, paste0(organ, "_top_tumor_induced.csv")),
            row.names = FALSE)
  
  cat("  Top 10 tumor-induced anchor-neighbor pairs:\n")
  print(head(top_tumor_induced, 10))
  
  # --------------------------------------------------------
  # Plots: per-anchor heatmap
  # --------------------------------------------------------
  cat(sprintf("\n[4/4] Plotting per-anchor heatmaps\n"))
  
  for (anc in anchors_present) {
    df <- anchor_enrichments[[anc]]
    if (is.null(df) || nrow(df) == 0) next
    
    # Top 10 neighbors by |delta|
    top_n <- df %>% arrange(desc(abs(delta))) %>% head(10)
    
    # Bar plot of delta
    p <- ggplot(top_n, aes(x = reorder(neighbor, delta),
                           y = delta,
                           fill = delta > 0)) +
      geom_col(alpha = 0.85) +
      geom_text(aes(label = sprintf("%+.1f", delta)),
                hjust = ifelse(top_n$delta >= 0, -0.1, 1.1), size = 3) +
      geom_hline(yintercept = 0, color = "gray40") +
      coord_flip() +
      scale_fill_manual(values = c("TRUE" = "#d62728", "FALSE" = "#1f77b4"),
                        guide = "none") +
      labs(title = sprintf("%s niche: neighbor enrichment change",
                           anc),
           subtitle = sprintf("Anchor: %s (tumor vs baseline, k=%d)",
                              anc, K_NEIGHBORS),
           x = "Neighbor cell type", y = "Δ % (tumor − baseline)") +
      theme_bw(base_size = 10)
    
    fn_safe <- gsub("[/ ]", "_", anc)
    ggsave(file.path(OUT_DIR, sprintf("%s_niche_%s_bar.pdf", organ, fn_safe)),
           p, width = 7, height = 5)
  }
  
  list(organ = organ, anchors = anchors_present,
       top_induced = head(top_tumor_induced, 10))
}

# ============================================================
# T3: Brain 4-axis triangle
# ============================================================
run_brain_4axis <- function(all_niche_brain) {
  cat(sprintf("\n\n### === Brain 4-axis triangle (Immune-Endo-Neuron-Glia) ===\n"))
  
  # Recode each cell to its axis
  axis_map <- do.call(c, lapply(names(BRAIN_AXES), function(nm) {
    setNames(rep(nm, length(BRAIN_AXES[[nm]])), BRAIN_AXES[[nm]])
  }))
  
  annotated <- all_niche_brain %>%
    mutate(anchor_axis = axis_map[anchor],
           neighbor_axis = axis_map[neighbor]) %>%
    filter(!is.na(anchor_axis) & !is.na(neighbor_axis))
  
  # Aggregate: axis × axis pairs
  axis_pairs <- annotated %>%
    group_by(sample, group, anchor_axis, neighbor_axis) %>%
    summarise(total_n = sum(n), .groups = "drop") %>%
    group_by(sample, anchor_axis) %>%
    mutate(pct = 100 * total_n / sum(total_n)) %>%
    ungroup()
  
  # Group-level mean
  group_axis <- axis_pairs %>%
    group_by(group, anchor_axis, neighbor_axis) %>%
    summarise(mean_pct = mean(pct), .groups = "drop") %>%
    pivot_wider(names_from = group, values_from = mean_pct,
                values_fill = 0) %>%
    mutate(delta = tumor - baseline)
  
  write.csv(group_axis,
            file.path(OUT_DIR, "brain_immune_endo_neuron_glia_triangle.csv"),
            row.names = FALSE)
  
  cat("  Axis × axis delta (tumor − baseline):\n")
  print(group_axis %>% arrange(desc(abs(delta))))
  
  # Heatmap matrix
  mat_tumor    <- group_axis %>% select(anchor_axis, neighbor_axis, tumor)    %>%
                  pivot_wider(names_from = neighbor_axis, values_from = tumor) %>%
                  as.data.frame()
  rownames(mat_tumor) <- mat_tumor$anchor_axis; mat_tumor$anchor_axis <- NULL
  
  mat_base    <- group_axis %>% select(anchor_axis, neighbor_axis, baseline) %>%
                  pivot_wider(names_from = neighbor_axis, values_from = baseline) %>%
                  as.data.frame()
  rownames(mat_base) <- mat_base$anchor_axis; mat_base$anchor_axis <- NULL
  
  mat_delta   <- group_axis %>% select(anchor_axis, neighbor_axis, delta) %>%
                  pivot_wider(names_from = neighbor_axis, values_from = delta) %>%
                  as.data.frame()
  rownames(mat_delta) <- mat_delta$anchor_axis; mat_delta$anchor_axis <- NULL
  
  pdf(file.path(OUT_DIR, "brain_immune_endo_neuron_glia_triangle.pdf"),
      width = 12, height = 4.5)
  par(mfrow = c(1, 1))
  layout(matrix(1:3, nrow = 1))
  
  pheatmap(as.matrix(mat_base),
           color = colorRampPalette(c("white", "steelblue"))(100),
           cluster_rows = FALSE, cluster_cols = FALSE,
           display_numbers = TRUE, number_format = "%.1f",
           main = "Baseline: anchor → neighbor (%)",
           fontsize_row = 10, fontsize_col = 10)
  dev.off()
  
  # 3-panel: tumor / baseline / delta
  mats <- list(baseline = mat_base, tumor = mat_tumor, delta = mat_delta)
  pal <- list(baseline = colorRampPalette(c("white", "steelblue"))(100),
              tumor    = colorRampPalette(c("white", "darkred"))(100),
              delta    = colorRampPalette(c("steelblue", "white", "darkred"))(100))
  
  for (nm in names(mats)) {
    pdf(file.path(OUT_DIR, sprintf("brain_triangle_%s.pdf", nm)),
        width = 5, height = 4)
    pheatmap(as.matrix(mats[[nm]]),
             color = pal[[nm]],
             cluster_rows = FALSE, cluster_cols = FALSE,
             display_numbers = TRUE, number_format = "%.1f",
             main = sprintf("Brain 4-axis: %s", nm),
             fontsize_row = 10, fontsize_col = 10)
    dev.off()
  }
  cat(sprintf("  -> brain_triangle_{baseline,tumor,delta}.pdf\n"))
  
  # Immune-Endo-Neuron (客户原话) focused plot
  ien_axes <- c("immune", "endo", "neuron")
  ien_data <- group_axis %>%
    filter(anchor_axis %in% ien_axes & neighbor_axis %in% ien_axes)
  
  write.csv(ien_data,
            file.path(OUT_DIR, "brain_immune_endo_neuron_focus.csv"),
            row.names = FALSE)
  
  p_ien <- ggplot(ien_data,
                  aes(x = paste(anchor_axis, "→", neighbor_axis),
                      y = delta, fill = delta > 0)) +
    geom_col(alpha = 0.85) +
    geom_hline(yintercept = 0, color = "gray40") +
    geom_text(aes(label = sprintf("%+.1f", delta)),
              vjust = ifelse(ien_data$delta >= 0, -0.3, 1.2), size = 3) +
    scale_fill_manual(values = c("TRUE" = "#d62728", "FALSE" = "#1f77b4"),
                      guide = "none") +
    labs(title = "Brain: immune × endo × neuron cross-axis shift",
         subtitle = "Client-requested analysis: 免疫-内皮-神经三角 (tumor - baseline)",
         x = NULL, y = "Δ % neighbor") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  
  ggsave(file.path(OUT_DIR, "brain_immune_endo_neuron_focus.pdf"),
         p_ien, width = 9, height = 5)
  ggsave(file.path(OUT_DIR, "brain_immune_endo_neuron_focus.png"),
         p_ien, width = 9, height = 5, dpi = 200)
  cat(sprintf("  -> brain_immune_endo_neuron_focus.{pdf,png}\n"))
  
  group_axis
}

# ============================================================
# Run
# ============================================================
liver_result <- run_organ_niche("liver", LIVER_ANCHORS)

brain_result <- run_organ_niche("brain", BRAIN_ANCHORS)

# T3: Brain 4-axis on all pairs (not just anchors)
if (!is.null(brain_result)) {
  all_niche_brain <- read.csv(file.path(OUT_DIR, "brain_all_niche_pairs.csv"),
                              stringsAsFactors = FALSE)
  brain_4axis <- run_brain_4axis(all_niche_brain)
}

# ============================================================
# Summary.md
# ============================================================
cat("\n### Writing summary.md ###\n")

md <- c(
  "# T1+T2+T3 Tumor-Induced Niche Analysis",
  "",
  sprintf("**Generated:** %s", format(Sys.time(), "%Y-%m-%d %H:%M")),
  "",
  "## Method",
  "",
  sprintf("- k-NN (k=%d) on spatial coordinates per sample", K_NEIGHBORS),
  "- Each anchor cell type's k neighbors' cell type composition",
  "- Tumor vs baseline mean % shift",
  "- Brain 4-axis: immune / endo / neuron / glia aggregation",
  "",
  "## Liver anchors tested",
  ""
)
if (!is.null(liver_result)) {
  md <- c(md, sprintf("- %s", liver_result$anchors), "",
          "Top 10 tumor-induced liver niches:", "",
          "| Anchor | Neighbor | Tumor % | Baseline % | Δ |",
          "|--------|----------|--------:|-----------:|--:|")
  for (i in seq_len(nrow(liver_result$top_induced))) {
    r <- liver_result$top_induced[i, ]
    md <- c(md, sprintf("| %s | %s | %.1f | %.1f | %+.1f |",
                        r$anchor, r$neighbor, r$tumor, r$baseline, r$delta))
  }
}
md <- c(md, "", "## Brain anchors tested", "")
if (!is.null(brain_result)) {
  md <- c(md, sprintf("- %s", brain_result$anchors), "",
          "Top 10 tumor-induced brain niches:", "",
          "| Anchor | Neighbor | Tumor % | Baseline % | Δ |",
          "|--------|----------|--------:|-----------:|--:|")
  for (i in seq_len(nrow(brain_result$top_induced))) {
    r <- brain_result$top_induced[i, ]
    md <- c(md, sprintf("| %s | %s | %.1f | %.1f | %+.1f |",
                        r$anchor, r$neighbor, r$tumor, r$baseline, r$delta))
  }
}
md <- c(md, "",
  "## Output files",
  "",
  "### Liver",
  "- liver_all_niche_pairs.csv, liver_niche_<anchor>_*.csv/pdf",
  "- liver_top_tumor_induced.csv",
  "",
  "### Brain",
  "- brain_all_niche_pairs.csv, brain_niche_<anchor>_*.csv/pdf",
  "- brain_top_tumor_induced.csv",
  "- **brain_immune_endo_neuron_glia_triangle.{csv,pdf}**  (T3, client-requested)",
  "- **brain_immune_endo_neuron_focus.{pdf,png}**  (T2, client's 3-axis triangle)")

writeLines(md, file.path(OUT_DIR, "summary.md"))

cat("\n============================================================\n")
cat("T1+T2+T3 COMPLETE\n")
cat(sprintf("Output: %s\n", OUT_DIR))
cat(sprintf("Files:  %d\n", length(list.files(OUT_DIR))))
cat("============================================================\n")
