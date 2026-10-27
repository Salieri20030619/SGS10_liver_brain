#!/usr/bin/env Rscript
# Step A: 计算 Microglia_DAM 每个 cell 的重分类决策, 写 CSV
# 内存友好: 单样本顺序处理, subset DAM 后立刻释放主 obj
# 不修改 RDS

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

PROJECT_ROOT <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
setwd(PROJECT_ROOT)
source("scripts/diff_common.R")

OUT_DIR <- file.path(PROJECT_ROOT, "dam_clean_v2")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# === Lineage gates (硬性必须满足) ===
LINEAGE_GATES <- list(
  Microglia = list(
    required_any = c("Csf1r","Hexb"),       # 必须至少一个表达
    supportive   = c("Cx3cr1","P2ry12","Tmem119","Aif1","Trem2","Fcrls")
  ),
  Astrocyte = list(
    required_any = c("Slc1a2","Slc1a3","Aqp4"),
    supportive   = c("Gfap","Aldh1l1","Sparcl1","Sox9","Gja1","Ntsr2")
  ),
  Pericyte = list(
    required_any = c("Pdgfrb","Rgs5"),
    supportive   = c("Vtn","Notch3","Cspg4","Kcnj8","Abcc9","Cald1")
  ),
  Endothelial = list(
    required_any = c("Cldn5","Cdh5","Pecam1"),
    supportive   = c("Flt1","Vwf","Tek","Adgrl4")
  ),
  Oligodendrocyte = list(
    required_any = c("Mbp","Plp1"),
    supportive   = c("Mog","Mobp","Olig1","Olig2","Sox10")
  ),
  Neuron = list(
    required_any = c("Snap25","Syn1","Map2"),
    supportive   = c("Rbfox3","Nefl","Nefm","Stmn2")
  )
)

# DAM sub-signature
DAM_SIGNATURE <- list(
  required_any = c("Trem2","Gpnmb","Spp1"),  # 真 DAM 必须有 DAM signature
  supportive   = c("Cst7","Lpl","Apoe","Cd9","Itgax","Axl")
)

# QC 阈值: percent.mt 高于此为 low quality
MT_HIGH_THRESHOLD <- 3.0

# Score 阈值: 决定是否归到某类
LINEAGE_SCORE_MIN <- 0.05  # AddModuleScore 的最低值, 太低则 ambiguous

samples_brain <- c("brain3","brain8","nao56","bla_brain","con_brain")

cat("========================================\n")
cat("DAM clean: compute reassignment per cell\n")
cat("========================================\n\n")

all_decisions <- list()

for (s in samples_brain) {
  fp <- tryCatch(resolve_rds_path(s, "v4"), error = function(e) NULL)
  if (is.null(fp)) {
    cat(sprintf("  %s: RDS not found, skip\n", s))
    next
  }

  cat(sprintf("\n--- %s ---\n", s))
  t0 <- proc.time()

  obj <- readRDS(fp)
  obj <- tryCatch(JoinLayers(obj), error = function(e) obj)
  cat(sprintf("  loaded: %d cells\n", ncol(obj)))

  md <- obj@meta.data
  dam_mask <- md$celltype_refined == "Microglia_DAM"
  dam_cells <- colnames(obj)[dam_mask]
  cat(sprintf("  DAM cells: %d\n", length(dam_cells)))

  if (length(dam_cells) < 10) {
    rm(obj); gc(verbose=FALSE)
    next
  }

  # Subset 到 DAM cells, 释放主 obj
  dam_obj <- subset(obj, cells = dam_cells)

  # 取 mt 信息 (不同列名兼容)
  mt_col <- intersect(c("percent.mt","percent.mito","mito"),
                      colnames(dam_obj@meta.data))[1]
  mt_pct <- if (!is.na(mt_col)) dam_obj@meta.data[[mt_col]] else rep(NA, ncol(dam_obj))

  # nFeature
  nfeat <- dam_obj$nFeature_RNA

  rm(obj); gc(verbose=FALSE)
  cat(sprintf("  main obj released, dam_obj: %d cells\n", ncol(dam_obj)))

  if (!"data" %in% Layers(dam_obj[["RNA"]])) {
    dam_obj <- NormalizeData(dam_obj, verbose = FALSE)
  }

  # 对每个 lineage 计算 score
  cat("  computing lineage scores...\n")
  lineage_scores <- list()
  lineage_required_pass <- list()

  for (ln_name in names(LINEAGE_GATES)) {
    ln <- LINEAGE_GATES[[ln_name]]
    all_genes <- c(ln$required_any, ln$supportive)
    genes_in <- intersect(all_genes, rownames(dam_obj))
    if (length(genes_in) < 2) {
      lineage_scores[[ln_name]] <- rep(NA, ncol(dam_obj))
      lineage_required_pass[[ln_name]] <- rep(FALSE, ncol(dam_obj))
      next
    }
    # AddModuleScore 在小 subset 上不稳, 用简单 mean expression 替代
    exp_mat <- GetAssayData(dam_obj, assay="RNA", layer="data")[genes_in, , drop=FALSE]
    lineage_scores[[ln_name]] <- colMeans(as.matrix(exp_mat))

    # required_any gate
    req_in <- intersect(ln$required_any, rownames(dam_obj))
    if (length(req_in) > 0) {
      req_mat <- GetAssayData(dam_obj, assay="RNA", layer="data")[req_in, , drop=FALSE]
      lineage_required_pass[[ln_name]] <- colSums(as.matrix(req_mat) > 0) > 0
    } else {
      lineage_required_pass[[ln_name]] <- rep(FALSE, ncol(dam_obj))
    }
  }

  # DAM signature score
  dam_sig_genes <- c(DAM_SIGNATURE$required_any, DAM_SIGNATURE$supportive)
  dam_sig_in <- intersect(dam_sig_genes, rownames(dam_obj))
  dam_sig_score <- if (length(dam_sig_in) >= 2) {
    colMeans(as.matrix(GetAssayData(dam_obj, assay="RNA", layer="data")[dam_sig_in, , drop=FALSE]))
  } else {
    rep(0, ncol(dam_obj))
  }
  dam_required_pass <- {
    req_in <- intersect(DAM_SIGNATURE$required_any, rownames(dam_obj))
    if (length(req_in) > 0) {
      req_mat <- GetAssayData(dam_obj, assay="RNA", layer="data")[req_in, , drop=FALSE]
      colSums(as.matrix(req_mat) > 0) > 0
    } else { rep(FALSE, ncol(dam_obj)) }
  }

  # === 决策逻辑 (策略 C) ===
  n_cells <- ncol(dam_obj)
  decisions <- character(n_cells)
  reasons <- character(n_cells)
  best_lineage <- character(n_cells)
  best_score <- numeric(n_cells)

  for (i in seq_len(n_cells)) {
    # 1. QC filter: high mt → low_quality
    if (!is.na(mt_pct[i]) && mt_pct[i] > MT_HIGH_THRESHOLD) {
      decisions[i] <- "low_quality_high_mt"
      reasons[i] <- sprintf("mt=%.1f%%", mt_pct[i])
      best_lineage[i] <- "NA"
      best_score[i] <- NA
      next
    }

    # 2. Find lineage with highest score AMONG those passing required gate
    passing_lineages <- names(lineage_required_pass)[
      sapply(lineage_required_pass, function(x) x[i])
    ]
    if (length(passing_lineages) == 0) {
      decisions[i] <- "low_quality_no_lineage"
      reasons[i] <- "no required marker"
      best_lineage[i] <- "NA"
      best_score[i] <- NA
      next
    }

    scores_i <- sapply(lineage_scores[passing_lineages], function(x) x[i])
    best_idx <- which.max(scores_i)
    best_ln <- passing_lineages[best_idx]
    best_sc <- scores_i[best_idx]

    best_lineage[i] <- best_ln
    best_score[i] <- best_sc

    # 3. If best lineage is Microglia + DAM signature met → real DAM
    if (best_ln == "Microglia" && dam_required_pass[i]) {
      decisions[i] <- "Microglia_DAM_strict"
      reasons[i] <- sprintf("MG_score=%.2f, DAM_sig=%.2f", best_sc, dam_sig_score[i])
      next
    }
    # 4. If best is Microglia but no DAM signature → Microglia_Homeostatic-like
    if (best_ln == "Microglia") {
      decisions[i] <- "Microglia_other"
      reasons[i] <- sprintf("MG_score=%.2f, no_DAM_sig", best_sc)
      next
    }
    # 5. Other lineage with high score → reassign
    if (best_sc > LINEAGE_SCORE_MIN) {
      decisions[i] <- best_ln
      reasons[i] <- sprintf("%s_score=%.2f", best_ln, best_sc)
      next
    }
    # 6. Score too low → ambiguous low quality
    decisions[i] <- "low_quality_ambiguous"
    reasons[i] <- sprintf("max_score=%.2f", best_sc)
  }

  # 输出 per-cell decision table
  decision_df <- data.frame(
    barcode = colnames(dam_obj),
    sample = s,
    nFeature_RNA = nfeat,
    percent_mt = mt_pct,
    Microglia_score = lineage_scores$Microglia,
    Microglia_pass = lineage_required_pass$Microglia,
    Astrocyte_score = lineage_scores$Astrocyte,
    Astrocyte_pass = lineage_required_pass$Astrocyte,
    Pericyte_score = lineage_scores$Pericyte,
    Pericyte_pass = lineage_required_pass$Pericyte,
    Endothelial_score = lineage_scores$Endothelial,
    Oligo_score = lineage_scores$Oligodendrocyte,
    Neuron_score = lineage_scores$Neuron,
    DAM_sig_score = dam_sig_score,
    DAM_sig_pass = dam_required_pass,
    best_lineage = best_lineage,
    best_score = best_score,
    new_label = decisions,
    reason = reasons,
    stringsAsFactors = FALSE
  )

  out_csv <- file.path(OUT_DIR, sprintf("%s_dam_decisions.csv", s))
  write.csv(decision_df, out_csv, row.names = FALSE)

  # 统计输出
  cat(sprintf("  decisions:\n"))
  dec_tab <- sort(table(decisions), decreasing = TRUE)
  for (k in names(dec_tab)) {
    cat(sprintf("    %-30s %5d  (%.1f%%)\n",
                k, dec_tab[k], dec_tab[k]/n_cells*100))
  }

  all_decisions[[s]] <- decision_df

  rm(dam_obj, decision_df); gc(verbose = FALSE)
  cat(sprintf("  elapsed: %.1f sec\n", (proc.time()-t0)["elapsed"]))
}

# 合并 + 总统计
cat("\n\n========================================\n")
cat("Overall reassignment summary\n")
cat("========================================\n")

merged_dec <- do.call(rbind, all_decisions)
write.csv(merged_dec, file.path(OUT_DIR, "all_dam_decisions.csv"),
          row.names = FALSE)

cat(sprintf("\nTotal DAM cells: %d\n", nrow(merged_dec)))
cat("\nNew label distribution (all 5 brain samples):\n")
overall_tab <- sort(table(merged_dec$new_label), decreasing=TRUE)
for (k in names(overall_tab)) {
  cat(sprintf("  %-30s %5d  (%.1f%%)\n",
              k, overall_tab[k], overall_tab[k]/nrow(merged_dec)*100))
}

cat("\nPer-sample × new label:\n")
print(table(merged_dec$sample, merged_dec$new_label))

cat(sprintf("\nOutput: %s\n", OUT_DIR))
cat("  - <sample>_dam_decisions.csv (5 个)\n")
cat("  - all_dam_decisions.csv (合并)\n")
