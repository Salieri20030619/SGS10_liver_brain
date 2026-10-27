#!/usr/bin/env Rscript
# Step B: 应用 dam_clean_decide.R 的决策, 修改 v4 RDS, 保存为 _v4_clean.rds
# 内存友好: 单样本顺序处理

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

PROJECT_ROOT <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
setwd(PROJECT_ROOT)
source("scripts/diff_common.R")

DECISION_DIR <- file.path(PROJECT_ROOT, "dam_clean_v2")
LOG_DIR <- file.path(PROJECT_ROOT, "scripts","logs")
dir.create(LOG_DIR, showWarnings=FALSE, recursive=TRUE)

# 决策映射: decision label → 最终 celltype_refined 值
# Microglia_other 归到 Microglia_Homeostatic
# Microglia_DAM_strict 保留 (但数量极少)
# low_quality_* 保留 low_quality 标签 (后续 DEG 排除)
# 其他直接用 lineage 名 (已经是 v4_refined 里存在的 broad 类)
DECISION_TO_CELLTYPE <- c(
  "Microglia_DAM_strict"   = "Microglia_DAM_strict",
  "Microglia_other"        = "Microglia_Homeostatic",
  "Astrocyte"              = "Astrocyte_Resting",  # 默认 Resting, 后续可细分
  "Pericyte"               = "Pericytes",
  "Endothelial"            = "Endothelial",
  "Oligodendrocyte"        = "Oligo_Mature",
  "Neuron"                 = "Neurons_unassigned",
  "low_quality_high_mt"    = "low_quality",
  "low_quality_no_lineage" = "low_quality"
)

samples_brain <- c("brain3","brain8","nao56","bla_brain","con_brain")

cat("========================================\n")
cat("Apply DAM clean decisions to brain RDS\n")
cat("========================================\n\n")

audit_summary <- list()

for (s in samples_brain) {
  fp <- tryCatch(resolve_rds_path(s, "v4"), error = function(e) NULL)
  if (is.null(fp) || !file.exists(fp)) {
    cat(sprintf("  %s: RDS not found, skip\n", s))
    next
  }

  decision_csv <- file.path(DECISION_DIR, sprintf("%s_dam_decisions.csv", s))
  if (!file.exists(decision_csv)) {
    cat(sprintf("  %s: decision CSV not found, skip\n", s))
    next
  }

  cat(sprintf("\n--- %s ---\n", s))
  t0 <- proc.time()

  obj <- readRDS(fp)
  cat(sprintf("  loaded: %d cells\n", ncol(obj)))

  decisions <- read.csv(decision_csv, stringsAsFactors = FALSE)
  cat(sprintf("  decisions: %d cells\n", nrow(decisions)))

  # backup 原 celltype_refined
  obj$celltype_refined_v4orig <- obj$celltype_refined

  # 应用 decision
  md <- obj@meta.data
  ref_new <- as.character(md$celltype_refined)

  matched <- 0
  for (i in seq_len(nrow(decisions))) {
    bc <- decisions$barcode[i]
    new_dec <- decisions$new_label[i]
    new_ct <- DECISION_TO_CELLTYPE[new_dec]
    if (is.na(new_ct)) {
      new_ct <- "Unknown"
    }

    cell_idx <- match(bc, colnames(obj))
    if (!is.na(cell_idx)) {
      ref_new[cell_idx] <- new_ct
      matched <- matched + 1
    }
  }
  cat(sprintf("  matched %d / %d barcodes\n", matched, nrow(decisions)))

  obj$celltype_refined <- ref_new

  # 同步 broad 注释
  if ("celltype_broad" %in% colnames(md)) {
    BROAD_MAP <- c(
      "Astrocyte_Resting"      = "Astrocyte",
      "Pericytes"              = "Vasculature",
      "Endothelial"            = "Vasculature",
      "Oligo_Mature"           = "Oligodendrocyte",
      "Neurons_unassigned"     = "Neuron",
      "Microglia_DAM_strict"   = "Microglia",
      "Microglia_Homeostatic"  = "Microglia",
      "low_quality"            = "low_quality",
      "Unknown"                = "Unknown"
    )
    broad_new <- as.character(md$celltype_broad)
    for (i in seq_len(nrow(decisions))) {
      bc <- decisions$barcode[i]
      cell_idx <- match(bc, colnames(obj))
      new_ref <- ref_new[cell_idx]
      if (!is.na(BROAD_MAP[new_ref])) {
        broad_new[cell_idx] <- BROAD_MAP[new_ref]
      }
    }
    obj$celltype_broad <- broad_new
  }

  # 添加 metadata 标记: 哪些 cell 是被 clean 重新分类的
  obj$dam_clean_modified <- colnames(obj) %in% decisions$barcode

  # 保存
  out_rds <- sub("\\.rds$", "_clean.rds", fp)
  saveRDS(obj, out_rds)
  cat(sprintf("  ✓ saved: %s\n", basename(out_rds)))

  # Audit: 修改前后的 celltype 分布
  cat("\n  Before vs After (修改前后的 brain 细胞分布):\n")
  before_tab <- table(obj$celltype_refined_v4orig)
  after_tab <- table(obj$celltype_refined)
  changed_types <- union(names(before_tab), names(after_tab))
  
  changes <- data.frame(
    celltype = changed_types,
    before = as.integer(before_tab[changed_types]),
    after = as.integer(after_tab[changed_types])
  )
  changes[is.na(changes)] <- 0
  changes$delta <- changes$after - changes$before
  changes <- changes[abs(changes$delta) > 0, ]
  changes <- changes[order(-abs(changes$delta)), ]
  print(changes)

  audit_summary[[s]] <- changes

  rm(obj); gc(verbose = FALSE)
  cat(sprintf("  elapsed: %.1f sec\n", (proc.time()-t0)["elapsed"]))
}

# 写一个汇总 audit 表
audit_df <- do.call(rbind, lapply(names(audit_summary), function(s) {
  df <- audit_summary[[s]]
  df$sample <- s
  df
}))

audit_csv <- file.path(DECISION_DIR, "clean_audit_summary.csv")
write.csv(audit_df, audit_csv, row.names = FALSE)
cat(sprintf("\n✓ Audit summary: %s\n", audit_csv))

cat("\n\n========================================\n")
cat("Cleaning complete. New RDS saved as <sample>_v4_clean.rds\n")
cat("========================================\n")
cat("Next steps:\n")
cat("  1. Verify Microglia_Homeostatic + Microglia_DAM_strict signal\n")
cat("  2. Re-run brain Stage 1 DEG with clean RDS\n")
cat("  3. Re-run Microglia pseudotime\n")
cat("  4. Update v2.5 brain 论点\n")
