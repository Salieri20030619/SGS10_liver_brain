#!/usr/bin/env Rscript
# =============================================================================
# diagnose_microenv_celltypes.R
# -----------------------------------------------------------------------------
# 快速诊断: Phase 2 输出 specificity 文件按 log2FC 排, 但有缺陷
# (cDC2 等微小 celltype 因 baseline 极低 → log2FC 高, 但绝对 mean_weight 很低)
# 
# 这个脚本按 absolute mean_weight 重新排序 + 显示 LSEC_Mid + Macrophage_M1 在
# 各微环境的位置, 验证微环境合并是否正确
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
})

PR <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
P2_DIR <- file.path(PR, "D1_v2_RCTD_composition/microenv_key_celltypes")

cat("=============================================================\n")
cat("微环境关键细胞诊断\n")
cat("=============================================================\n\n")

for (organ in c("liver", "brain")) {
  cat(sprintf("\n========== %s ==========\n\n", toupper(organ)))
  
  full <- read.csv(file.path(P2_DIR, paste0(organ, "_microenv_celltypes_full.csv")))
  
  # 1. Top 5 by absolute mean_weight per microenv
  cat("### Top 5 celltypes by absolute mean_weight (NOT log2FC) ###\n")
  top5_abs <- full %>%
    group_by(microenv_id, microenv_label, n_spots) %>%
    arrange(desc(mean_weight)) %>%
    slice_head(n = 5) %>%
    summarise(
      top_5 = paste(sprintf("%s(%.3f)", celltype, mean_weight),
                    collapse = " | "),
      .groups = "drop"
    ) %>%
    arrange(microenv_id)
  
  for (i in seq_len(nrow(top5_abs))) {
    cat(sprintf("**%s %s** (n=%d):\n  %s\n",
                top5_abs$microenv_id[i],
                top5_abs$microenv_label[i],
                top5_abs$n_spots[i],
                top5_abs$top_5[i]))
  }
  
  # 2. 关键 celltype 在各微环境的位置
  cat("\n### Key celltype tracking across microenvs ###\n")
  
  if (organ == "liver") {
    key_cts <- c("LSEC_Mid", "Macrophage_M1", "Kupffer", "B_Plasma",
                 "HSC_Quiescent", "Hep_Perivenous", "Hep_Periportal",
                 "LSEC", "Pericytes")
  } else {
    key_cts <- c("Hpc_CA1", "Hpc_CA3", "Hpc_DG",
                 "Excitatory_L23", "Excitatory_L6", "Inhibitory_Meis2",
                 "Endothelial", "Oligo_Mature", "Astrocyte_Resting", "Neurons")
  }
  
  for (ct in key_cts) {
    rows <- full %>%
      filter(celltype == ct) %>%
      arrange(desc(mean_weight)) %>%
      mutate(rank = row_number())
    
    if (nrow(rows) > 0) {
      cat(sprintf("\n%s (across all microenvs, sorted by mean_weight):\n", ct))
      for (i in seq_len(nrow(rows))) {
        marker <- if (i == 1) " ⭐" else ""
        cat(sprintf("  %s %s: mean_weight=%.3f%s\n",
                    rows$microenv_id[i], rows$microenv_label[i],
                    rows$mean_weight[i], marker))
      }
    }
  }
  
  # 3. 对每微环境的命名 vs 实际 top1 做对照
  cat("\n### Naming consistency check ###\n")
  cat("Expected: microenv_label 中的 celltype 应该是该 microenv 的 top1\n\n")
  
  top1 <- full %>%
    group_by(microenv_id, microenv_label) %>%
    arrange(desc(mean_weight)) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    select(microenv_id, microenv_label, top1_celltype = celltype, top1_weight = mean_weight)
  
  for (i in seq_len(nrow(top1))) {
    label <- top1$microenv_label[i]
    top1_ct <- top1$top1_celltype[i]
    # 看 label 里是否包含 top1_celltype 的关键词
    consistent <- grepl(top1_ct, label, fixed = TRUE) ||
      (top1_ct == "Hep_Perivenous" && grepl("zonation|Hepatocyte", label)) ||
      (top1_ct == "Hep_Periportal" && grepl("zonation|Hepatocyte", label)) ||
      (top1_ct == "Hpc_CA1" && grepl("Hippocampus|Hpc", label)) ||
      (top1_ct == "Hpc_CA3" && grepl("Hippocampus|Hpc", label)) ||
      (top1_ct == "Hpc_DG" && grepl("Hippocampus|Hpc", label)) ||
      (top1_ct == "LSEC" && grepl("Sinusoidal|LSEC", label)) ||
      (top1_ct == "Pericytes" && grepl("Vascular|Pericytes", label))
    
    flag <- if (consistent) "✓" else "⚠ MISMATCH"
    cat(sprintf("%s %s '%s' → actual top1 = %s (%.3f)\n",
                flag,
                top1$microenv_id[i],
                top1$microenv_label[i],
                top1$top1_celltype[i],
                top1$top1_weight[i]))
  }
}

cat("\n=============================================================\n")
cat("诊断完成\n")
cat("=============================================================\n")1~#!/usr/bin/env Rscript
# =============================================================================
# diagnose_microenv_celltypes.R
# -----------------------------------------------------------------------------
# 快速诊断: Phase 2 输出 specificity 文件按 log2FC 排, 但有缺陷
# (cDC2 等微小 celltype 因 baseline 极低 → log2FC 高, 但绝对 mean_weight 很低)
# 
# 这个脚本按 absolute mean_weight 重新排序 + 显示 LSEC_Mid + Macrophage_M1 在
# 各微环境的位置, 验证微环境合并是否正确
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
})

PR <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
P2_DIR <- file.path(PR, "D1_v2_RCTD_composition/microenv_key_celltypes")

cat("=============================================================\n")
cat("微环境关键细胞诊断\n")
cat("=============================================================\n\n")

for (organ in c("liver", "brain")) {
  cat(sprintf("\n========== %s ==========\n\n", toupper(organ)))
  
  full <- read.csv(file.path(P2_DIR, paste0(organ, "_microenv_celltypes_full.csv")))
  
  # 1. Top 5 by absolute mean_weight per microenv
  cat("### Top 5 celltypes by absolute mean_weight (NOT log2FC) ###\n")
  top5_abs <- full %>%
    group_by(microenv_id, microenv_label, n_spots) %>%
    arrange(desc(mean_weight)) %>%
    slice_head(n = 5) %>%
    summarise(
      top_5 = paste(sprintf("%s(%.3f)", celltype, mean_weight),
                    collapse = " | "),
      .groups = "drop"
    ) %>%
    arrange(microenv_id)
  
  for (i in seq_len(nrow(top5_abs))) {
    cat(sprintf("**%s %s** (n=%d):\n  %s\n",
                top5_abs$microenv_id[i],
                top5_abs$microenv_label[i],
                top5_abs$n_spots[i],
                top5_abs$top_5[i]))
  }
  
  # 2. 关键 celltype 在各微环境的位置
  cat("\n### Key celltype tracking across microenvs ###\n")
  
  if (organ == "liver") {
    key_cts <- c("LSEC_Mid", "Macrophage_M1", "Kupffer", "B_Plasma",
                 "HSC_Quiescent", "Hep_Perivenous", "Hep_Periportal",
                 "LSEC", "Pericytes")
  } else {
    key_cts <- c("Hpc_CA1", "Hpc_CA3", "Hpc_DG",
                 "Excitatory_L23", "Excitatory_L6", "Inhibitory_Meis2",
                 "Endothelial", "Oligo_Mature", "Astrocyte_Resting", "Neurons")
  }
  
  for (ct in key_cts) {
    rows <- full %>%
      filter(celltype == ct) %>%
      arrange(desc(mean_weight)) %>%
      mutate(rank = row_number())
    
    if (nrow(rows) > 0) {
      cat(sprintf("\n%s (across all microenvs, sorted by mean_weight):\n", ct))
      for (i in seq_len(nrow(rows))) {
        marker <- if (i == 1) " ⭐" else ""
        cat(sprintf("  %s %s: mean_weight=%.3f%s\n",
                    rows$microenv_id[i], rows$microenv_label[i],
                    rows$mean_weight[i], marker))
      }
    }
  }
  
  # 3. 对每微环境的命名 vs 实际 top1 做对照
  cat("\n### Naming consistency check ###\n")
  cat("Expected: microenv_label 中的 celltype 应该是该 microenv 的 top1\n\n")
  
  top1 <- full %>%
    group_by(microenv_id, microenv_label) %>%
    arrange(desc(mean_weight)) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    select(microenv_id, microenv_label, top1_celltype = celltype, top1_weight = mean_weight)
  
  for (i in seq_len(nrow(top1))) {
    label <- top1$microenv_label[i]
    top1_ct <- top1$top1_celltype[i]
    # 看 label 里是否包含 top1_celltype 的关键词
    consistent <- grepl(top1_ct, label, fixed = TRUE) ||
      (top1_ct == "Hep_Perivenous" && grepl("zonation|Hepatocyte", label)) ||
      (top1_ct == "Hep_Periportal" && grepl("zonation|Hepatocyte", label)) ||
      (top1_ct == "Hpc_CA1" && grepl("Hippocampus|Hpc", label)) ||
      (top1_ct == "Hpc_CA3" && grepl("Hippocampus|Hpc", label)) ||
      (top1_ct == "Hpc_DG" && grepl("Hippocampus|Hpc", label)) ||
      (top1_ct == "LSEC" && grepl("Sinusoidal|LSEC", label)) ||
      (top1_ct == "Pericytes" && grepl("Vascular|Pericytes", label))
    
    flag <- if (consistent) "✓" else "⚠ MISMATCH"
    cat(sprintf("%s %s '%s' → actual top1 = %s (%.3f)\n",
                flag,
                top1$microenv_id[i],
                top1$microenv_label[i],
                top1$top1_celltype[i],
                top1$top1_weight[i]))
  }
}

cat("\n=============================================================\n")
cat("诊断完成\n")
cat("=============================================================\n")
