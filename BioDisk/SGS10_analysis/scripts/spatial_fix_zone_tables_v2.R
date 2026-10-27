#!/usr/bin/env Rscript
# v2: 自动检测 gene 命名 (mm10--- 前缀 vs 无前缀), 修复所有样本的 zone table

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

PROJECT_ROOT <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
setwd(PROJECT_ROOT)
OUT_DIR <- file.path(PROJECT_ROOT, "spatial_results")
ZONE_DIR <- file.path(OUT_DIR, "zone_tables")

# Mouse gene symbols (不含前缀)
MOUSE_SYMBOLS <- c(
  # HSC / fibrosis
  "Acta2","Col1a1","Col1a2","Col3a1","Col4a1","Col4a2",
  "Tagln","Lox","Loxl1","Loxl2","Myh11","Fn1",
  # TGFb pathway
  "Tgfb1","Tgfb2","Tgfb3",
  "Smad3","Smad2","Smad7",
  "Pmepa1","Serpine1",
  "Tgfbr1","Tgfbr2","Tgfbr3",
  # PARs
  "F2","Plg","F2r","Pard3",
  # Myeloid hub
  "Csf1","Csf1r","Gas6","Mertk","C3","Grn",
  # NAM / myeloid
  "Cd74","Cd68","Itgb2","Itgam","Itgax",
  # Cholangiocyte
  "Krt19","Krt18","Spp1","Sox9",
  # Hepatocyte
  "Alb","Glul","Cyp2e1","Pck1","Cps1","Apoa1","Apob",
  # Angio
  "Vegfa","Vegfb","Vegfc",
  # Cytokines
  "Il6","Ifng","Tnf","Il1b",
  # Brain: microglia
  "Aif1","Trem2","P2ry12","Apoe","Tmem119","Cx3cr1",
  # Brain: BBB / MHC-I
  "H2-D1","H2-K1","B2m","Elmo1",
  # Brain: astrocyte
  "Gfap","Aldh1l1","Slc1a3",
  # Brain: oligo
  "Mbp","Mobp","Olig1","Olig2","Plp1",
  # Brain: endothelial
  "Cldn5","Pecam1","Cdh5","Vwf"
)

HUMAN_SYMBOLS <- c(
  "AFP","DLK1","IGF2","ALB","APOB","TF",
  "MKI67","EPCAM",
  "TGFB1","TGFB2","TGFB3","VEGFA","CSF1",
  "F2","PLG","FN1","COL1A1",
  "IL6","IFNG","TNF","IL1B",
  "ACTA2","VIM"
)

# 所有样本 (包括已经成功的 3 个 tumor, 为了统一重生成)
ALL_SAMPLES <- c("Liver_3","Liver_8","gan56",
                 "Bla_liver_1","Con_liver_1",
                 "Brain_3","Brain_8","nao56",
                 "Bla_brain_1","Con_brain_1")

cat("========================================\n")
cat("Fix zone tables v2: auto-detect gene naming\n")
cat("========================================\n\n")

for (sname in ALL_SAMPLES) {
  rds_path <- file.path(OUT_DIR, sprintf("%s_seekspace.rds", sname))
  if (!file.exists(rds_path)) {
    cat(sprintf("  ⚠ %s: RDS not found, skip\n", sname))
    next
  }

  cat(sprintf("\n--- %s ---\n", sname))
  t0 <- proc.time()

  obj <- readRDS(rds_path)
  all_genes <- rownames(obj)
  n_cells <- ncol(obj)

  # 检测命名
  has_mm10_prefix <- any(grepl("^mm10---", all_genes))
  has_human <- any(grepl("^GRCh38-", all_genes))

  cat(sprintf("  genes: %d | mm10 prefix: %s | human: %s | spots: %d\n",
              length(all_genes), has_mm10_prefix, has_human, n_cells))

  # 确保 zone 列正确长度
  if (is.null(obj$zone) || length(obj$zone) != n_cells) {
    if ("human_pct" %in% colnames(obj@meta.data) && max(obj$human_pct) > 0) {
      h_pct <- obj$human_pct
      obj$zone <- ifelse(h_pct > 70, "tumor",
                  ifelse(h_pct >= 30, "margin",
                  ifelse(h_pct < 5, "host", "host_dom")))
    } else {
      obj$zone <- rep("host", n_cells)
    }
  }

  # 根据命名决定 mouse marker 前缀
  if (has_mm10_prefix) {
    mouse_markers <- paste0("mm10---", MOUSE_SYMBOLS)
  } else {
    mouse_markers <- MOUSE_SYMBOLS
  }

  mouse_in <- intersect(mouse_markers, all_genes)

  # Human markers (如果是双物种)
  human_in <- character(0)
  if (has_human) {
    human_markers <- paste0("GRCh38-", HUMAN_SYMBOLS)
    human_in <- intersect(human_markers, all_genes)
  }

  markers_in <- c(mouse_in, human_in)
  cat(sprintf("  markers found: %d mouse + %d human = %d\n",
              length(mouse_in), length(human_in), length(markers_in)))

  if (length(markers_in) == 0) {
    cat(sprintf("  ⚠ no markers, skip\n"))
    next
  }

  data_layer <- GetAssayData(obj, assay = "RNA", layer = "data")

  # 生成 zone table
  zone_df <- data.frame()
  for (zn in unique(obj$zone)) {
    sp <- which(obj$zone == zn)
    if (length(sp) < 5) next
    for (g in markers_in) {
      vals <- data_layer[g, sp]
      zone_df <- rbind(zone_df, data.frame(
        sample = sname,
        zone = zn,
        n_spots = length(sp),
        gene = g,
        pct = round(mean(vals > 0) * 100, 2),
        mean = round(mean(vals), 4),
        stringsAsFactors = FALSE
      ))
    }
  }

  # 加 species + gene_clean (for cross-sample matching)
  zone_df$species <- ifelse(grepl("^GRCh38-", zone_df$gene), "human", "mouse")
  zone_df$gene_clean <- gsub("^mm10---", "", zone_df$gene)
  zone_df$gene_clean <- gsub("^GRCh38-", "", zone_df$gene_clean)

  out_csv <- file.path(ZONE_DIR, sprintf("%s_zone_markers.csv", sname))
  write.csv(zone_df, out_csv, row.names = FALSE)
  cat(sprintf("  ✓ wrote %d rows (zones: %s)\n",
              nrow(zone_df),
              paste(unique(zone_df$zone), collapse=",")))

  rm(obj); gc(verbose = FALSE)
  cat(sprintf("  elapsed: %.1f sec\n", (proc.time() - t0)["elapsed"]))
}

# Merge
cat("\n========================================\n")
cat("Merging all zone tables\n")
cat("========================================\n")

zone_files <- list.files(ZONE_DIR, pattern = "_zone_markers\\.csv$",
                         full.names = TRUE)
all_zone <- do.call(rbind, lapply(zone_files, read.csv, stringsAsFactors = FALSE))

write.csv(all_zone,
          file.path(OUT_DIR, "summary", "all_samples_zone_markers.csv"),
          row.names = FALSE)

cat(sprintf("\nTotal: %d rows, %d samples\n",
            nrow(all_zone), length(unique(all_zone$sample))))
cat("\nSample × zone counts:\n")
print(table(all_zone$sample, all_zone$zone))
