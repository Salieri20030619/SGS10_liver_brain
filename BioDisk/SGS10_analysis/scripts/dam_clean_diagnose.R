#!/usr/bin/env Rscript
# DAM 清洗诊断: 算每个 DAM cell 的 true DAM score vs pericyte score
# 用 5 个 brain 样本全部, 不仅是 bla_brain

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

PROJECT_ROOT <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
setwd(PROJECT_ROOT)
source("scripts/diff_common.R")

OUT_DIR <- file.path(PROJECT_ROOT, "dam_clean_diagnose")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# True DAM panel (tissue-resident DAM / LAM signature)
DAM_GENES <- c("Trem2","Spp1","Gpnmb","Cst7","Lpl","Apoe",
               "Ctsd","Ctsb","Fabp5","Lgals3")

# Pericyte / vasculature panel
PERI_GENES <- c("Vtn","Pdgfrb","Notch3","Rgs5","Kcnj8","Abcc9",
                "Cald1","Atp13a5","Carmn","Ebf1")

# 其他 microglia marker (对照, 真 microglia 应该高表达)
MG_CORE <- c("Aif1","P2ry12","Tmem119","Cx3cr1","Csf1r")

samples_brain <- c("brain3","brain8","nao56","bla_brain","con_brain")

all_scores <- list()

for (s in samples_brain) {
  fp <- tryCatch(resolve_rds_path(s, "v4"), error = function(e) NULL)
  if (is.null(fp)) next
  cat(sprintf("\n--- %s ---\n", s))

  obj <- readRDS(fp)
  obj <- tryCatch(JoinLayers(obj), error = function(e) obj)

  # Subset Microglia_DAM
  md <- obj@meta.data
  dam_col <- if ("celltype_refined" %in% colnames(md)) "celltype_refined" else "celltype_v32"
  dam_cells <- which(md[[dam_col]] == "Microglia_DAM")
  cat(sprintf("  DAM cells: %d\n", length(dam_cells)))

  if (length(dam_cells) < 10) { rm(obj); next }

  dam_obj <- subset(obj, cells = colnames(obj)[dam_cells])

  # Calculate scores (use NormalizeData if not present)
  if (!"data" %in% Layers(dam_obj[["RNA"]])) {
    dam_obj <- NormalizeData(dam_obj, verbose = FALSE)
  }

  # AddModuleScore for 3 panels
  dam_obj <- AddModuleScore(dam_obj, features = list(DAM_GENES),
                            name = "DAM_score_", assay = "RNA")
  dam_obj <- AddModuleScore(dam_obj, features = list(PERI_GENES),
                            name = "Peri_score_", assay = "RNA")
  dam_obj <- AddModuleScore(dam_obj, features = list(MG_CORE),
                            name = "MG_core_score_", assay = "RNA")

  scores_df <- data.frame(
    barcode = colnames(dam_obj),
    sample = s,
    DAM_score = dam_obj$DAM_score_1,
    Peri_score = dam_obj$Peri_score_1,
    MG_core_score = dam_obj$MG_core_score_1
  )

  # 也加上几个关键 marker 的表达 (直接看)
  scores_df$Trem2_expr <- FetchData(dam_obj, "Trem2")[, 1]
  scores_df$Vtn_expr <- FetchData(dam_obj, "Vtn")[, 1]
  scores_df$Pdgfrb_expr <- FetchData(dam_obj, "Pdgfrb")[, 1]
  scores_df$Aif1_expr <- FetchData(dam_obj, "Aif1")[, 1]
  scores_df$P2ry12_expr <- FetchData(dam_obj, "P2ry12")[, 1]

  all_scores[[s]] <- scores_df

  cat(sprintf("  DAM_score: mean %.2f (range %.2f to %.2f)\n",
              mean(scores_df$DAM_score), min(scores_df$DAM_score),
              max(scores_df$DAM_score)))
  cat(sprintf("  Peri_score: mean %.2f (range %.2f to %.2f)\n",
              mean(scores_df$Peri_score), min(scores_df$Peri_score),
              max(scores_df$Peri_score)))
  cat(sprintf("  Aif1 pct>0: %.1f%%\n",
              mean(scores_df$Aif1_expr > 0) * 100))
  cat(sprintf("  Vtn pct>0: %.1f%%\n",
              mean(scores_df$Vtn_expr > 0) * 100))

  rm(obj, dam_obj); gc(verbose = FALSE)
}

# Merge all
merged <- do.call(rbind, all_scores)
rownames(merged) <- NULL

# 分类 heuristic (初步): Peri_score > DAM_score 记为 "suspected pericyte"
merged$tentative_label <- ifelse(
  merged$Peri_score > merged$DAM_score & merged$Peri_score > 0.1,
  "suspected_pericyte",
  ifelse(
    merged$DAM_score > merged$Peri_score & merged$DAM_score > 0.1,
    "likely_DAM",
    "ambiguous"
  )
)

# Save
write.csv(merged, file.path(OUT_DIR, "dam_scoring_all_samples.csv"),
          row.names = FALSE)

cat("\n========================================\n")
cat("Summary tentative labels:\n")
cat("========================================\n")
print(table(merged$sample, merged$tentative_label))

cat("\nScore stats per tentative label:\n")
print(merged %>% group_by(tentative_label) %>%
        summarise(
          n = n(),
          DAM_mean = round(mean(DAM_score), 2),
          Peri_mean = round(mean(Peri_score), 2),
          MG_core_mean = round(mean(MG_core_score), 2),
          Aif1_pct = round(mean(Aif1_expr > 0) * 100, 1),
          Vtn_pct = round(mean(Vtn_expr > 0) * 100, 1)
        ))

# Plot scatter
p1 <- ggplot(merged, aes(x = DAM_score, y = Peri_score,
                          color = tentative_label)) +
  geom_point(alpha = 0.6, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
  facet_wrap(~sample, ncol = 3) +
  theme_bw() +
  labs(title = "DAM score vs Pericyte score per cell",
       subtitle = "Cells labeled Microglia_DAM by v4_refined, 5 brain samples")

ggsave(file.path(OUT_DIR, "dam_peri_scatter_per_sample.pdf"),
       plot = p1, width = 12, height = 8)

p2 <- ggplot(merged, aes(x = tentative_label, y = DAM_score,
                         fill = tentative_label)) +
  geom_violin() +
  theme_bw() +
  labs(title = "DAM score by tentative label")

p3 <- ggplot(merged, aes(x = tentative_label, y = Peri_score,
                         fill = tentative_label)) +
  geom_violin() +
  theme_bw() +
  labs(title = "Pericyte score by tentative label")

library(patchwork)
ggsave(file.path(OUT_DIR, "score_distributions.pdf"),
       plot = p2 / p3, width = 8, height = 10)

cat(sprintf("\n✓ Output: %s\n", OUT_DIR))
cat("  - dam_scoring_all_samples.csv\n")
cat("  - dam_peri_scatter_per_sample.pdf\n")
cat("  - score_distributions.pdf\n")
