#!/usr/bin/env Rscript
# Step C: 用 clean RDS 重跑脑端 Stage 1 DEG
# 内存友好: 单 contrast 顺序处理, 每次只加载需要的样本

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

PROJECT_ROOT <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
setwd(PROJECT_ROOT)

OUT_DIR <- file.path(PROJECT_ROOT, "diff_stage1_brain_clean")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "deg_v4_refined"), showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "summary"), showWarnings = FALSE)

# 用 clean RDS 路径 (注意是 manual 不是 mouse)
get_clean_brain_rds <- function(sample) {
  fp <- file.path(PROJECT_ROOT, sample, "mouse",
                  sprintf("%s_manual_annotated_v4_clean.rds", sample))
  if (!file.exists(fp)) {
    stop(sprintf("Clean RDS not found: %s", fp))
  }
  fp
}

# 三个 brain 对比
CONTRASTS <- list(
  brain_Tumor_vs_Control = list(
    test = c("brain3","brain8","nao56"),
    ctrl = c("con_brain")
  ),
  brain_Tumor_vs_Blank = list(
    test = c("brain3","brain8","nao56"),
    ctrl = c("bla_brain")
  ),
  brain_Blank_vs_Control = list(
    test = c("bla_brain"),
    ctrl = c("con_brain")
  )
)

# DEG 配置
DEG_PARAMS <- list(
  min.pct = 0.1,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  max.cells.per.ident = 3000  # 加速 + 节省内存
)

# 排除 cell types: low_quality 等
EXCLUDE_CT <- c("low_quality","Unknown","Other","Microglia_M1","Microglia_M2",
                "Neurons_unassigned")

cat("========================================\n")
cat("Brain Stage 1 DEG (with cleaned RDS)\n")
cat("========================================\n\n")

# helper: 加载并合并多样本
load_merged <- function(samples, contrast_name) {
  cat(sprintf("  loading %d samples: %s\n",
              length(samples), paste(samples, collapse=",")))
  obj_list <- list()
  for (s in samples) {
    fp <- get_clean_brain_rds(s)
    o <- readRDS(fp)
    o <- tryCatch(JoinLayers(o), error = function(e) o)
    o$contrast_group <- contrast_name
    o$sample_id <- s
    # rename cells to avoid collision
    o <- RenameCells(o, add.cell.id = s)
    obj_list[[s]] <- o
  }
  if (length(obj_list) == 1) {
    merged <- obj_list[[1]]
  } else {
    merged <- merge(obj_list[[1]], y = obj_list[-1])
    merged <- tryCatch(JoinLayers(merged), error = function(e) merged)
  }
  rm(obj_list); gc(verbose=FALSE)
  merged
}

# 对每个 contrast 跑 DEG
all_summary <- list()

for (cn in names(CONTRASTS)) {
  cfg <- CONTRASTS[[cn]]
  cat(sprintf("\n\n========================================\n"))
  cat(sprintf("Contrast: %s\n", cn))
  cat(sprintf("  test: %s\n", paste(cfg$test, collapse=",")))
  cat(sprintf("  ctrl: %s\n", paste(cfg$ctrl, collapse=",")))
  cat(sprintf("========================================\n"))

  t0 <- proc.time()

  # 加载 test 组
  test_obj <- load_merged(cfg$test, "test")
  cat(sprintf("  test merged: %d cells\n", ncol(test_obj)))
  test_obj$group <- "test"

  # 加载 ctrl 组
  ctrl_obj <- load_merged(cfg$ctrl, "ctrl")
  cat(sprintf("  ctrl merged: %d cells\n", ncol(ctrl_obj)))
  ctrl_obj$group <- "ctrl"

  # 合并
  combined <- merge(test_obj, y = ctrl_obj)
  combined <- tryCatch(JoinLayers(combined), error = function(e) combined)
  rm(test_obj, ctrl_obj); gc(verbose=FALSE)
  cat(sprintf("  combined: %d cells\n", ncol(combined)))

  if (!"data" %in% Layers(combined[["RNA"]])) {
    combined <- NormalizeData(combined, verbose = FALSE)
  }

  # 决定 cell types: 包括所有 brain cell types, 但排除 low_quality
  cts <- setdiff(unique(combined$celltype_refined), EXCLUDE_CT)
  cts <- cts[!is.na(cts)]
  ct_counts <- table(combined$celltype_refined)
  cts <- cts[ct_counts[cts] >= 20]  # 至少 20 cells 才跑
  cat(sprintf("  cell types to test: %d\n", length(cts)))

  contrast_dir <- file.path(OUT_DIR, "deg_v4_refined", cn)
  dir.create(contrast_dir, showWarnings = FALSE, recursive = TRUE)

  contrast_summary <- data.frame()

  for (ct in cts) {
    # subset
    sub_cells <- which(combined$celltype_refined == ct)
    if (length(sub_cells) < 20) next

    sub <- subset(combined, cells = colnames(combined)[sub_cells])
    Idents(sub) <- sub$group
    n_test <- sum(sub$group == "test")
    n_ctrl <- sum(sub$group == "ctrl")

    if (n_test < 5 || n_ctrl < 5) {
      cat(sprintf("    %s: skip (test=%d, ctrl=%d)\n", ct, n_test, n_ctrl))
      next
    }

    deg <- tryCatch({
      FindMarkers(sub,
                  ident.1 = "test",
                  ident.2 = "ctrl",
                  min.pct = DEG_PARAMS$min.pct,
                  logfc.threshold = DEG_PARAMS$logfc.threshold,
                  test.use = DEG_PARAMS$test.use,
                  max.cells.per.ident = DEG_PARAMS$max.cells.per.ident,
                  verbose = FALSE)
    }, error = function(e) {
      cat(sprintf("    %s: FAIL - %s\n", ct, e$message))
      NULL
    })

    if (is.null(deg) || nrow(deg) == 0) {
      cat(sprintf("    %s: no DEG\n", ct))
      next
    }

    deg$gene <- rownames(deg)
    deg$cell_type <- ct
    deg$n_test <- n_test
    deg$n_ctrl <- n_ctrl
    deg <- deg %>%
      select(cell_type, gene, avg_log2FC, pct.1, pct.2,
             p_val, p_val_adj, n_test, n_ctrl) %>%
      arrange(p_val_adj, desc(abs(avg_log2FC)))

    out_csv <- file.path(contrast_dir, sprintf("%s_deg.csv",
                                                 gsub("/","_",ct)))
    write.csv(deg, out_csv, row.names = FALSE)

    n_up <- sum(deg$avg_log2FC > 0 & deg$p_val_adj < 0.05)
    n_down <- sum(deg$avg_log2FC < 0 & deg$p_val_adj < 0.05)
    cat(sprintf("    %-25s test=%5d ctrl=%5d  up=%4d down=%4d\n",
                ct, n_test, n_ctrl, n_up, n_down))

    contrast_summary <- rbind(contrast_summary, data.frame(
      contrast = cn,
      cell_type = ct,
      n_test_cells = n_test,
      n_ctrl_cells = n_ctrl,
      n_DEG_total = nrow(deg),
      n_DEG_up_sig = n_up,
      n_DEG_down_sig = n_down
    ))

    rm(sub, deg); gc(verbose = FALSE)
  }

  all_summary[[cn]] <- contrast_summary

  rm(combined); gc(verbose = FALSE)

  elapsed <- (proc.time() - t0)["elapsed"] / 60
  cat(sprintf("\n  contrast elapsed: %.1f min\n", elapsed))
}

# 写汇总表
final_summary <- do.call(rbind, all_summary)
write.csv(final_summary,
          file.path(OUT_DIR, "summary", "deg_summary.csv"),
          row.names = FALSE)

cat("\n\n========================================\n")
cat("DONE\n")
cat("========================================\n")
cat("Output dir:", OUT_DIR, "\n\n")
print(final_summary)
