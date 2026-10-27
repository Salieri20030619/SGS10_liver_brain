#!/usr/bin/env Rscript
# ============================================================
# Secondary Annotation Pipeline v4
# ============================================================
# 用法:
#   Rscript secondary_annotation.R <sample_name>
#   Rscript secondary_annotation.R con_liver
#   Rscript secondary_annotation.R liver8     # 自动处理 mouse + human 两个 RDS
#
# 输入:
#   ~/SGS10_analysis/<sample>/<sample>_*_annotated.rds       (纯鼠样本)
#   ~/SGS10_analysis/<sample>/mouse/<sample>_mouse_*.rds     (PDX 鼠源)
#   ~/SGS10_analysis/<sample>/human/<sample>_human_*.rds     (PDX 人源)
#
# 输出 (在原 RDS 同目录下):
#   <sample>_*_v4.rds                  -- 增加 3 列 metadata 的新 RDS
#   <sample>_v4_annotation_report.md   -- 二次注释结果对比报告
#   <sample>_v4_umap.pdf               -- v3.2 vs v4 UMAP 对比图
#
# 设计:
#   1) 加载对应 marker 库 (liver mouse / liver human / brain)
#   2) AddModuleScore 给每个 cell 打分所有亚型
#   3) 每个 cell 取最高 score 的亚型作为 celltype_refined
#   4) 计算 score margin (top1 - top2),结合 cluster size 标 confidence
#   5) 应用 mapping CSV 生成 celltype_broad
#   6) 保存 RDS + 生成 markdown report + UMAP PDF
#
# 置信度规则 (严格版):
#   - high:   该 refined 亚型 ≥100 cells AND score margin > median(margin)
#   - medium: 30 ≤ 细胞数 < 100  OR score margin <= median(margin)
#   - low:    细胞数 < 30  (保留标注但下游脚本可能 filter)
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

# ============================================================
# 0. 解析参数 & 路径
# ============================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("用法: Rscript secondary_annotation.R <sample_name>")
}
SAMPLE <- args[1]

PROJECT_ROOT <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
SCRIPT_DIR   <- file.path(PROJECT_ROOT, "scripts")
SAMPLE_DIR   <- file.path(PROJECT_ROOT, SAMPLE)

if (!dir.exists(SAMPLE_DIR)) {
  stop("样本目录不存在: ", SAMPLE_DIR)
}

# 自动判断是 PDX 还是纯鼠样本,自动判断器官
PDX_LIVER_SAMPLES  <- c("liver8", "gan56", "liver3")
PURE_LIVER_SAMPLES <- c("con_liver", "bla_liver")
BRAIN_SAMPLES      <- c("con_brain", "bla_brain", "nao56", "brain3", "brain8")

is_pdx_liver  <- SAMPLE %in% PDX_LIVER_SAMPLES
is_pure_liver <- SAMPLE %in% PURE_LIVER_SAMPLES
is_brain      <- SAMPLE %in% BRAIN_SAMPLES

if (!(is_pdx_liver || is_pure_liver || is_brain)) {
  stop("未知样本: ", SAMPLE)
}

cat("==========================================\n")
cat("Secondary Annotation v4\n")
cat("Sample:", SAMPLE, "\n")
cat("Type:  ",
    ifelse(is_pdx_liver,  "PDX liver (mouse + human)",
    ifelse(is_pure_liver, "Pure mouse liver", "Pure mouse brain")), "\n")
cat("Time:  ", format(Sys.time()), "\n")
cat("==========================================\n\n")

# ============================================================
# 1. 加载 marker 库 + mapping CSV
# ============================================================
cat("[1/6] 加载 marker 库和 mapping...\n")

if (is_brain) {
  source(file.path(SCRIPT_DIR, "marker_lib_brain_v4.R"))
} else {
  source(file.path(SCRIPT_DIR, "marker_lib_liver_v4.R"))
}

mapping_df <- read.csv(file.path(SCRIPT_DIR, "celltype_mapping_broad.csv"),
                       comment.char = "#", stringsAsFactors = FALSE,
                       strip.white = TRUE)
mapping_df <- mapping_df[mapping_df$refined != "" & !is.na(mapping_df$refined), ]

cat("  Mapping rules:", nrow(mapping_df), "rows\n\n")

# ============================================================
# 2. 找到所有需要处理的 RDS
# ============================================================
cat("[2/6] 寻找输入 RDS 文件...\n")

find_rds_files <- function(sample_dir, sample_name, is_pdx) {
  result <- list()

  if (is_pdx) {
    # PDX: mouse + human 两个子目录
    for (sp in c("mouse", "human")) {
      sub_dir <- file.path(sample_dir, sp)
      if (dir.exists(sub_dir)) {
        rds_files <- list.files(sub_dir,
                                pattern = paste0("^", sample_name, "_", sp,
                                                 ".*annotated\\.rds$"),
                                full.names = TRUE)
        # 过滤掉 v4 输出
        rds_files <- rds_files[!grepl("_v4\\.rds$", rds_files)]
        if (length(rds_files) > 0) {
          # 取最新的(按 mtime)
          latest <- rds_files[which.max(file.mtime(rds_files))]
          result[[sp]] <- latest
        }
      }
    }
  } else {
    # 纯样本: 直接在 sample_dir 下
    rds_files <- list.files(sample_dir,
                            pattern = paste0("^", sample_name,
                                             ".*annotated\\.rds$"),
                            full.names = TRUE)
    rds_files <- rds_files[!grepl("_v4\\.rds$", rds_files)]
    if (length(rds_files) > 0) {
      latest <- rds_files[which.max(file.mtime(rds_files))]
      species <- if (is_brain) "mouse_brain" else "mouse_liver"
      result[[species]] <- latest
    }
  }

  return(result)
}

rds_inputs <- find_rds_files(SAMPLE_DIR, SAMPLE, is_pdx_liver)

if (length(rds_inputs) == 0) {
  stop("未找到任何 RDS 文件,请检查路径: ", SAMPLE_DIR)
}

for (sp in names(rds_inputs)) {
  cat("  [", sp, "]", basename(rds_inputs[[sp]]), "\n")
}
cat("\n")

# ============================================================
# 3. 核心函数: 对一个 Seurat 对象做二次注释
# ============================================================
annotate_with_markers <- function(seu, marker_list, mapping_df,
                                  species_label = "mouse",
                                  organ_label = "liver") {

  cat("  [3a] 过滤 marker 库 (只保留 RDS 中存在的基因)...\n")
  available_genes <- rownames(seu)
  marker_filtered <- list()
  marker_skip_log <- character(0)
  for (ct in names(marker_list)) {
    genes <- intersect(marker_list[[ct]], available_genes)
    if (length(genes) >= 3) {
      marker_filtered[[ct]] <- genes
    } else {
      marker_skip_log <- c(marker_skip_log,
                           sprintf("    SKIP %s: only %d/%d markers found",
                                   ct, length(genes), length(marker_list[[ct]])))
    }
  }
  if (length(marker_skip_log) > 0) {
    cat(paste(marker_skip_log, collapse = "\n"), "\n")
  }
  cat("    Kept", length(marker_filtered), "/", length(marker_list), "cell types\n")

  if (length(marker_filtered) == 0) {
    stop("没有任何 marker 集通过过滤,marker 库可能不匹配该物种")
  }

  cat("  [3b] AddModuleScore 给每个亚型打分...\n")
  # 用临时 name 防止覆盖
  score_names <- paste0("MS_", seq_along(marker_filtered))
  seu <- AddModuleScore(seu,
                        features = marker_filtered,
                        name = "MS_",
                        seed = 42)
  # AddModuleScore 会自动加 1,2,3... 后缀
  score_cols <- paste0("MS_", seq_along(marker_filtered))
  score_mat  <- as.matrix(seu@meta.data[, score_cols])
  colnames(score_mat) <- names(marker_filtered)

  cat("  [3c] 分配 celltype_refined (取每个 cell 最高分亚型)...\n")
  refined <- colnames(score_mat)[apply(score_mat, 1, which.max)]
  names(refined) <- rownames(score_mat)  # 必须加 cell barcode 作为 names

  # 计算 score margin (top1 - top2)
  top2_score <- apply(score_mat, 1, function(x) {
    sorted <- sort(x, decreasing = TRUE)
    sorted[1] - sorted[2]
  })
  names(top2_score) <- rownames(score_mat)

  cat("  [3d] 计算置信度 (high/medium/low)...\n")
  refined_sizes <- table(refined)
  margin_threshold <- median(top2_score)

  confidence <- mapply(function(ct, margin) {
    n <- refined_sizes[ct]
    if (n >= 100 && margin > margin_threshold) {
      "high"
    } else if (n >= 30) {
      "medium"
    } else {
      "low"
    }
  }, refined, top2_score, USE.NAMES = FALSE)
  names(confidence) <- rownames(score_mat)

  cat("    margin threshold (median):", round(margin_threshold, 4), "\n")
  cat("    confidence dist: high=", sum(confidence == "high"),
      " medium=", sum(confidence == "medium"),
      " low=", sum(confidence == "low"), "\n", sep = "")

  cat("  [3e] 应用 mapping 生成 celltype_broad...\n")
  # 限定 mapping 到合适的 organ
  if (organ_label == "liver_human") {
    map_subset <- mapping_df[mapping_df$organ == "liver_human", ]
  } else if (organ_label == "liver") {
    map_subset <- mapping_df[mapping_df$organ %in% c("liver"), ]
  } else {
    map_subset <- mapping_df[mapping_df$organ == "brain", ]
  }

  # 同名 refined 取第一条
  map_subset <- map_subset[!duplicated(map_subset$refined), ]
  refined_to_broad <- setNames(map_subset$broad, map_subset$refined)

  broad <- refined_to_broad[refined]
  broad[is.na(broad)] <- "Unknown_NotMapped"
  names(broad) <- rownames(score_mat)  # 关键: refined_to_broad[refined] 会丢失 names

  unmapped <- unique(refined[is.na(refined_to_broad[refined])])
  if (length(unmapped) > 0) {
    cat("    WARNING: 以下 refined 类型未在 mapping CSV 中找到 broad 对应:\n")
    for (u in unmapped) cat("      -", u, "\n")
  }

  cat("  [3f] 写入 metadata 列...\n")
  # 删除临时 score 列
  for (sc in score_cols) seu[[sc]] <- NULL

  # 关键: 去掉 names,让 Seurat 按位置匹配 cell barcode
  # 否则 broad <- refined_to_broad[refined] 产生的命名向量会被 Seurat
  # 误用 names 去匹配 cell ID,导致 "No cell overlap" 错误
  refined    <- unname(refined)
  broad      <- unname(broad)
  confidence <- unname(confidence)

  seu$celltype_refined    <- refined
  seu$celltype_broad      <- broad
  seu$celltype_confidence <- confidence
  seu$celltype_v32        <- seu$celltype  # 备份原标注

  return(list(
    seurat = seu,
    score_matrix = score_mat,
    margin_threshold = margin_threshold,
    n_marker_used = length(marker_filtered),
    n_marker_total = length(marker_list)
  ))
}

# ============================================================
# 4. 处理每个 RDS
# ============================================================
cat("[3/6] 开始注释...\n\n")

results_per_species <- list()

for (sp in names(rds_inputs)) {
  rds_path <- rds_inputs[[sp]]
  cat("--- Processing [", sp, "] ---\n")
  cat("Loading:", basename(rds_path), "\n")

  seu <- readRDS(rds_path)
  cat("  Cells:", ncol(seu), " Features:", nrow(seu), "\n")

  # 选择 marker 库 + organ label
  if (is_brain) {
    marker_list  <- brain_markers_v4
    organ_label  <- "brain"
  } else if (sp == "human") {
    marker_list  <- liver_markers_human_v4
    organ_label  <- "liver_human"
  } else {
    marker_list  <- liver_markers_mouse_v4
    organ_label  <- "liver"
  }

  ann_result <- annotate_with_markers(seu, marker_list, mapping_df,
                                      species_label = sp,
                                      organ_label = organ_label)

  # 保存新 RDS
  out_path <- sub("\\.rds$", "_v4.rds", rds_path)
  saveRDS(ann_result$seurat, out_path)
  cat("  Saved:", basename(out_path), "\n\n")

  results_per_species[[sp]] <- list(
    seurat = ann_result$seurat,
    out_path = out_path,
    in_path = rds_path,
    margin_threshold = ann_result$margin_threshold,
    n_marker_used = ann_result$n_marker_used,
    n_marker_total = ann_result$n_marker_total
  )
}

# ============================================================
# 5. 生成 UMAP 对比图
# ============================================================
cat("[4/6] 生成 UMAP 对比图...\n")

generate_umap_comparison <- function(seu, sp, sample) {
  if (!"umap" %in% names(seu@reductions)) {
    cat("  WARNING: ", sp, " 没有 umap reduction,跳过 UMAP 图\n", sep = "")
    return(NULL)
  }

  p1 <- DimPlot(seu, group.by = "celltype_v32", label = TRUE,
                repel = TRUE, label.size = 2.5) +
    ggtitle(paste0(sample, " [", sp, "] - v3.2 (original)")) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 6),
          plot.title = element_text(size = 9))

  p2 <- DimPlot(seu, group.by = "celltype_refined", label = TRUE,
                repel = TRUE, label.size = 2.5) +
    ggtitle(paste0(sample, " [", sp, "] - v4 refined")) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 6),
          plot.title = element_text(size = 9))

  p3 <- DimPlot(seu, group.by = "celltype_broad", label = TRUE,
                repel = TRUE, label.size = 2.5) +
    ggtitle(paste0(sample, " [", sp, "] - v4 broad")) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 6),
          plot.title = element_text(size = 9))

  p4 <- DimPlot(seu, group.by = "celltype_confidence",
                cols = c("high" = "#2ca02c", "medium" = "#ff7f0e",
                         "low" = "#d62728")) +
    ggtitle(paste0(sample, " [", sp, "] - confidence")) +
    theme(legend.position = "bottom",
          plot.title = element_text(size = 9))

  return((p1 | p2) / (p3 | p4))
}

umap_pdf_path <- file.path(SAMPLE_DIR, paste0(SAMPLE, "_v4_umap.pdf"))
pdf(umap_pdf_path, width = 14, height = 12)
for (sp in names(results_per_species)) {
  p <- generate_umap_comparison(results_per_species[[sp]]$seurat, sp, SAMPLE)
  if (!is.null(p)) print(p)
}
dev.off()
cat("  Saved:", basename(umap_pdf_path), "\n\n")

# ============================================================
# 6. 生成 Markdown 报告
# ============================================================
cat("[5/6] 生成 markdown report...\n")

report_path <- file.path(SAMPLE_DIR, paste0(SAMPLE, "_v4_annotation_report.md"))
report_lines <- character(0)

add <- function(...) report_lines <<- c(report_lines, paste0(...))

add("# Secondary Annotation v4 — ", SAMPLE)
add("")
add("**Generated:** ", format(Sys.time()))
add("")
add("**Sample type:** ",
    ifelse(is_pdx_liver,  "PDX liver (mouse + human)",
    ifelse(is_pure_liver, "Pure mouse liver", "Pure mouse brain")))
add("")
add("**Strategy:** Plan 3 (dual annotation: refined + broad + confidence)")
add("")
add("**Marker library sources:**")
if (is_brain) {
  for (s in brain_marker_sources) add("- ", s)
} else {
  for (s in liver_marker_sources) add("- ", s)
}
add("")
add("---")
add("")

for (sp in names(results_per_species)) {
  res <- results_per_species[[sp]]
  seu <- res$seurat
  add("## [", sp, "]")
  add("")
  add("- **Input RDS:** `", basename(res$in_path), "`")
  add("- **Output RDS:** `", basename(res$out_path), "`")
  add("- **Cells:** ", ncol(seu))
  add("- **Markers used:** ", res$n_marker_used, " / ", res$n_marker_total,
      " cell types from library")
  add("- **Margin threshold (median):** ", round(res$margin_threshold, 4))
  add("")

  # confidence 统计
  conf_tab <- table(seu$celltype_confidence)
  add("### Confidence distribution")
  add("")
  add("| Confidence | N cells | % |")
  add("|---|---|---|")
  for (lv in c("high", "medium", "low")) {
    n <- ifelse(lv %in% names(conf_tab), conf_tab[[lv]], 0)
    pct <- round(n / ncol(seu) * 100, 2)
    add("| ", lv, " | ", n, " | ", pct, "% |")
  }
  add("")

  # v3.2 vs v4 refined 分布对比
  add("### v3.2 → v4 refined cell type distribution")
  add("")
  add("**v3.2 (original):**")
  add("")
  v32_tab <- sort(table(seu$celltype_v32), decreasing = TRUE)
  add("| Cell type | N | % |")
  add("|---|---|---|")
  for (ct in names(v32_tab)) {
    n <- v32_tab[[ct]]
    pct <- round(n / ncol(seu) * 100, 2)
    add("| ", ct, " | ", n, " | ", pct, "% |")
  }
  add("")
  add("**v4 refined (literature-based):**")
  add("")
  v4_tab <- sort(table(seu$celltype_refined), decreasing = TRUE)
  v4_conf <- table(seu$celltype_refined, seu$celltype_confidence)
  add("| Cell type | N | % | Confidence |")
  add("|---|---|---|---|")
  for (ct in names(v4_tab)) {
    n <- v4_tab[[ct]]
    pct <- round(n / ncol(seu) * 100, 2)
    # 取该亚型最常见的 confidence 等级
    conf_for_ct <- v4_conf[ct, ]
    main_conf <- names(conf_for_ct)[which.max(conf_for_ct)]
    flag <- ifelse(main_conf == "low", " ⚠️",
            ifelse(main_conf == "medium", " ◐", " ✓"))
    add("| ", ct, " | ", n, " | ", pct, "% | ", main_conf, flag, " |")
  }
  add("")

  add("**v4 broad (collapsed):**")
  add("")
  vb_tab <- sort(table(seu$celltype_broad), decreasing = TRUE)
  add("| Broad class | N | % |")
  add("|---|---|---|")
  for (ct in names(vb_tab)) {
    n <- vb_tab[[ct]]
    pct <- round(n / ncol(seu) * 100, 2)
    add("| ", ct, " | ", n, " | ", pct, "% |")
  }
  add("")
  add("---")
  add("")
}

add("## Notes")
add("")
add("- Confidence levels: ✓ high (≥100 cells + score margin > median), ",
    "◐ medium (30-99 cells or borderline), ⚠️ low (<30 cells, kept but flagged)")
add("- Low confidence cells are NOT removed; downstream analysis scripts should ",
    "decide whether to filter based on the confidence column.")
add("- Original v3.2 annotations preserved in `celltype_v32` metadata column.")
add("- This is automated assignment based on AddModuleScore. Manual ",
    "verification of biologically meaningful subtypes is recommended.")

writeLines(report_lines, report_path)
cat("  Saved:", basename(report_path), "\n\n")

# ============================================================
# 7. 完成
# ============================================================
cat("[6/6] 完成!\n\n")
cat("==========================================\n")
cat("Output files (in ", SAMPLE_DIR, "):\n", sep = "")
for (sp in names(results_per_species)) {
  cat("  - ", basename(results_per_species[[sp]]$out_path), "\n", sep = "")
}
cat("  - ", basename(umap_pdf_path), "\n", sep = "")
cat("  - ", basename(report_path), "\n", sep = "")
cat("==========================================\n")
