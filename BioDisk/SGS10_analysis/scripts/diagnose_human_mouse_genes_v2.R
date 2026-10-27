# =============================================================
# diagnose_human_mouse_genes_v2.R
# Purpose:
#   1. Inspect gene name format in official SeekSpace tumor RDS
#   2. Inspect already-split human/mouse annotated RDS (gan56)
#   3. Check Liver_3 / Liver_8 split status
# =============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

# -------------------------------------------------------------
# Three categories of RDS to inspect
# -------------------------------------------------------------

OFFICIAL_RDS <- list(
  gan56    = "/mnt/BioData/bio/SGS10例/SGS24070264-gan56/data/rds/24070264_gan_56_expression/gan_56.rds",
  Liver_3  = "/mnt/BioData/bio/SGS10例/SGS25090559-8例样本数据/data/rds/25090559_Liver_3_expression/Liver_3.rds",
  Liver_8  = "/mnt/BioData/bio/SGS10例/SGS25090559-8例样本数据/data/rds/25090559_Liver_8_expression/Liver_8.rds"
)

# Already-split RDS (found via `find`)
SPLIT_RDS <- list(
  gan56_human = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/gan56/human/gan56_human_annotated_v4.rds",
  gan56_mouse = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/gan56/mouse/gan56_mouse_annotated_v4.rds"
)

# Check if Liver_3 / Liver_8 have the same split structure
cat("=================================================\n")
cat(" Checking for existing Liver_3 / Liver_8 human/mouse splits:\n")
cat("=================================================\n")
for (sample in c("Liver_3", "Liver_8", "liver_3", "liver_8")) {
  for (sub in c("human", "mouse")) {
    pattern <- sprintf("/mnt/BioDisk/wangmingyuan/SGS10_analysis/%s/%s", sample, sub)
    if (dir.exists(pattern)) {
      rds_files <- list.files(pattern, pattern = "\\.rds$", full.names = TRUE)
      for (f in rds_files) cat("  FOUND:", f, "\n")
    }
  }
}
cat("\n")

# -------------------------------------------------------------
# Diagnose function (reusable)
# -------------------------------------------------------------
diagnose_one <- function(label, rds_path) {
  cat("\n#################################################\n")
  cat("# [", label, "]\n", sep = "")
  cat("# Path:", rds_path, "\n")
  cat("#################################################\n")
  
  if (!file.exists(rds_path)) {
    cat("!! File does not exist. Skipping.\n\n")
    return(invisible(NULL))
  }
  
  obj <- tryCatch(readRDS(rds_path), error = function(e) {
    cat("ERROR reading RDS:", conditionMessage(e), "\n")
    NULL
  })
  if (is.null(obj)) return(invisible(NULL))
  
  cat("Object class:", paste(class(obj), collapse = ", "), "\n")
  
  # Try to extract count matrix
  counts <- NULL
  if (inherits(obj, "Seurat")) {
    assay_names <- Assays(obj)
    cat("Assays:", paste(assay_names, collapse = ", "), "\n")
    # try RNA first, then Spatial, then first available
    assay_to_use <- intersect(c("RNA", "Spatial"), assay_names)[1]
    if (is.na(assay_to_use)) assay_to_use <- assay_names[1]
    cat("Using assay:", assay_to_use, "\n")
    counts <- tryCatch(GetAssayData(obj, assay = assay_to_use, slot = "counts"),
                       error = function(e) GetAssayData(obj, assay = assay_to_use, layer = "counts"))
  } else if (inherits(obj, "dgCMatrix") || inherits(obj, "Matrix") || is.matrix(obj)) {
    counts <- obj
  } else if (is.list(obj)) {
    cat("List element names:", paste(head(names(obj), 10), collapse = ", "), "\n")
    # SeekSpace official RDS might be a list
    if ("counts" %in% names(obj))          counts <- obj$counts
    else if ("RNA" %in% names(obj))        counts <- obj$RNA
    else if ("expression" %in% names(obj)) counts <- obj$expression
    else if ("matrix" %in% names(obj))     counts <- obj$matrix
    else {
      # try the first matrix-like element
      for (nm in names(obj)) {
        if (inherits(obj[[nm]], c("dgCMatrix", "Matrix")) || is.matrix(obj[[nm]])) {
          counts <- obj[[nm]]; cat("Extracted matrix from list$", nm, "\n", sep = ""); break
        }
      }
    }
  }
  
  if (is.null(counts)) {
    cat("Could not extract count matrix. Object structure:\n")
    print(str(obj, max.level = 2))
    rm(obj); gc(verbose = FALSE); return(invisible(NULL))
  }
  
  genes <- rownames(counts)
  n_genes <- length(genes)
  n_spots <- ncol(counts)
  
  cat(sprintf("Matrix: %d genes x %d spots\n\n", n_genes, n_spots))
  
  # ---- First 15 gene names ----
  cat(">>> First 15 gene names:\n")
  print(head(genes, 15))
  cat("\n")
  
  # ---- Prefix patterns ----
  cat(">>> Prefix counts:\n")
  patterns <- list(
    "GRCh38-/GRCh38_/hg38-"          = "^(GRCh38[-_]|hg38[-_])",
    "mm10-/mm10_/GRCm38-/GRCm39-"    = "^(mm10[-_]|GRCm3[89][-_])",
    "h-/h_/H-/H_/human_"             = "^([hH][-_]|human_)",
    "m-/m_/M-/M_/mouse_"             = "^([mM][-_]|mouse_)",
    "ENSG (human Ensembl)"           = "^ENSG\\d{11}",
    "ENSMUSG (mouse Ensembl)"        = "^ENSMUSG\\d{11}"
  )
  for (p_name in names(patterns)) {
    n_match <- sum(grepl(patterns[[p_name]], genes))
    if (n_match > 0) {
      cat(sprintf("  %-40s %d (%.1f%%)\n", p_name, n_match, 100*n_match/n_genes))
    }
  }
  cat("\n")
  
  # ---- Capitalization patterns ----
  genes_stripped <- sub("^(GRCh38[-_]|hg38[-_]|mm10[-_]|GRCm3[89][-_]|[hHmM][-_]|human_|mouse_)", "", genes)
  is_allcaps <- grepl("^[A-Z][A-Z0-9-]*$", genes_stripped) & nchar(genes_stripped) >= 2
  is_titlecase <- grepl("^[A-Z][a-z]", genes_stripped)
  cat(sprintf(">>> Capitalization (after prefix strip):\n"))
  cat(sprintf("    ALLCAPS  (human-like, e.g. KRT19): %d (%.1f%%)\n",
              sum(is_allcaps), 100*sum(is_allcaps)/n_genes))
  cat(sprintf("    Titlecase (mouse-like, e.g. Krt19): %d (%.1f%%)\n",
              sum(is_titlecase), 100*sum(is_titlecase)/n_genes))
  cat("\n")
  
  # ---- Canonical marker test ----
  markers <- c("KRT19","ALB","HSPA8","EPCAM","GAPDH","ACTB")
  cat(">>> Marker presence check:\n")
  for (m in markers) {
    forms <- c(m, tools::toTitleCase(tolower(m)),
               paste0("GRCh38-", m), paste0("mm10-", tools::toTitleCase(tolower(m))),
               paste0("h-", m), paste0("m-", tools::toTitleCase(tolower(m))))
    found <- forms[forms %in% genes]
    if (length(found) > 0) cat(sprintf("    %s: found as %s\n", m, paste(found, collapse = ", ")))
  }
  cat("\n")
  
  # ---- Per-spot human/mouse fraction ----
  cat(">>> Per-spot species fraction:\n")
  
  human_idx <- grepl("^(GRCh38[-_]|hg38[-_]|[hH][-_]|human_|ENSG\\d)", genes)
  mouse_idx <- grepl("^(mm10[-_]|GRCm3[89][-_]|[mM][-_]|mouse_|ENSMUSG\\d)", genes)
  
  if (sum(human_idx) == 0 && sum(mouse_idx) == 0) {
    cat("    (no prefix -> using capitalization)\n")
    human_idx <- is_allcaps
    mouse_idx <- is_titlecase
  }
  
  if (sum(human_idx) > 0 && sum(mouse_idx) > 0) {
    cat(sprintf("    human genes: %d, mouse genes: %d\n", sum(human_idx), sum(mouse_idx)))
    h_per_spot <- Matrix::colSums(counts[human_idx, , drop = FALSE])
    m_per_spot <- Matrix::colSums(counts[mouse_idx, , drop = FALSE])
    total <- h_per_spot + m_per_spot
    h_frac <- h_per_spot / pmax(total, 1)
    
    cat(sprintf("    total human reads: %s\n", format(sum(h_per_spot), big.mark = ",")))
    cat(sprintf("    total mouse reads: %s\n", format(sum(m_per_spot), big.mark = ",")))
    cat(sprintf("    overall human fraction: %.2f%%\n\n", 100*sum(h_per_spot)/sum(total)))
    
    cat("    Per-spot human mRNA fraction quantiles:\n")
    qs <- quantile(h_frac, c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1.0), na.rm = TRUE)
    for (i in seq_along(qs)) cat(sprintf("      %5s:  %.3f\n", names(qs)[i], qs[i]))
    
    cat("\n    Spot counts by human fraction threshold:\n")
    for (thr in c(0.1, 0.3, 0.5, 0.7, 0.9)) {
      n <- sum(h_frac >= thr, na.rm = TRUE)
      cat(sprintf("      >= %.0f%%: %5d (%.1f%% of %d spots)\n",
                  100*thr, n, 100*n/n_spots, n_spots))
    }
  } else if (sum(human_idx) > 0) {
    cat(sprintf("    Pure human genome: %d genes\n", sum(human_idx)))
  } else if (sum(mouse_idx) > 0) {
    cat(sprintf("    Pure mouse genome: %d genes\n", sum(mouse_idx)))
  } else {
    cat("    Cannot classify species from gene names.\n")
  }
  
  cat("\n")
  rm(obj, counts); gc(verbose = FALSE)
  invisible(NULL)
}

# -------------------------------------------------------------
# Run diagnostics
# -------------------------------------------------------------
cat("=================================================\n")
cat(" SECTION 1: Official SeekSpace RDS (tumor liver)\n")
cat("=================================================\n")
for (nm in names(OFFICIAL_RDS)) diagnose_one(nm, OFFICIAL_RDS[[nm]])

cat("\n=================================================\n")
cat(" SECTION 2: Already-split human/mouse RDS (gan56)\n")
cat("=================================================\n")
for (nm in names(SPLIT_RDS)) diagnose_one(nm, SPLIT_RDS[[nm]])

cat("\n=================================================\n")
cat(" Diagnosis complete. Send the entire log to Salieri.\n")
cat("=================================================\n")
