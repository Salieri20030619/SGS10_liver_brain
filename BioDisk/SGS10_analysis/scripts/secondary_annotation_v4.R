#!/usr/bin/env Rscript
# ============================================================
# Secondary Annotation Pipeline v4  (rewritten)
# ============================================================
# 策略: 选项 A + 选项 3 + 阈值 (c)
#   - 信任 v3.2 final_celltype,v4 做验证+精修+broad,不推倒重来
#   - 只对免疫系粗类做子类精修 (DC/T/B/脑 Microglia)
#   - score 验证阈值 top3: v3.2 名进 marker score top3 算 agree,
#     否则翻到 score top1
#
# 用法:
#   Rscript secondary_annotation_v4.R <sample_name>
#   Rscript secondary_annotation_v4.R con_liver
#   Rscript secondary_annotation_v4.R liver8     # 自动处理 mouse + human
#
# 输入 (严格只从子目录读取):
#   <ROOT>/<sample>/mouse/<sample>_*annotated.rds   (纯鼠 & PDX 鼠源)
#   <ROOT>/<sample>/human/<sample>_*annotated.rds   (PDX 人源)
#
# 输出 (写回 <sample>/<species>/ 同目录):
#   <sample>_<sp>_annotated_v4.rds         加 6 列 metadata 的新 RDS
#   <sample>_v4_annotation_report.md       样本级报告
#   <sample>_v4_umap.pdf                   UMAP 对比图
#
# 新增的 metadata 列:
#   celltype_v32                  原 final_celltype 备份
#   celltype_refined              v4 精修后的细类 (最终推荐)
#   celltype_broad                按 mapping CSV 归的粗类
#   celltype_refinement_source    alias | rescore_agree | rescore_flipped |
#                                 subtype_rescore | passthrough | no_markers
#   celltype_v32_rank             v3.2 名在 marker score 里的排名 (1 = top)
#   celltype_confidence           high | medium | low
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

# ============================================================
# 0. 参数 & 路径
# ============================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("用法: Rscript secondary_annotation_v4.R <sample_name>")
}
SAMPLE <- args[1]

PROJECT_ROOT <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
SCRIPT_DIR   <- file.path(PROJECT_ROOT, "scripts")
SAMPLE_DIR   <- file.path(PROJECT_ROOT, SAMPLE)

# 固定配置
SCORE_MARGIN_FIXED   <- 0.05   # 不再用 median 自适应
RESCORE_TOPK         <- 3      # top3 内算 agree
MIN_MARKERS_PER_SET  <- 3

if (!dir.exists(SAMPLE_DIR)) {
  stop("样本目录不存在: ", SAMPLE_DIR)
}

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
cat("Secondary Annotation v4 (rewrite)\n")
cat("Sample:", SAMPLE, "\n")
cat("Type:  ",
    ifelse(is_pdx_liver,  "PDX liver (mouse + human)",
    ifelse(is_pure_liver, "Pure mouse liver", "Pure mouse brain")), "\n")
cat("Time:  ", format(Sys.time()), "\n")
cat("==========================================\n\n")

# ============================================================
# 1. 加载 marker 库 + alias + mapping
# ============================================================
cat("[1/6] 加载 marker 库 / aliases / mapping...\n")

if (is_brain) {
  source(file.path(SCRIPT_DIR, "marker_lib_brain_v4.R"))
} else {
  source(file.path(SCRIPT_DIR, "marker_lib_liver_v4.R"))
}
source(file.path(SCRIPT_DIR, "aliases_v4.R"))

mapping_df <- read.csv(file.path(SCRIPT_DIR, "celltype_mapping_broad.csv"),
                       comment.char = "#", stringsAsFactors = FALSE,
                       strip.white = TRUE)
mapping_df <- mapping_df[mapping_df$refined != "" & !is.na(mapping_df$refined), ]

cat("  Mapping rules:", nrow(mapping_df), "\n\n")

# ============================================================
# 2. 严格从子目录找 RDS
# ============================================================
cat("[2/6] 寻找输入 RDS (仅从 mouse/ 和 human/ 子目录)...\n")

find_rds_files <- function(sample_dir, sample_name, is_pdx) {
  result <- list()
  species_list <- if (is_pdx) c("mouse", "human") else "mouse"

  for (sp in species_list) {
    sub_dir <- file.path(sample_dir, sp)
    if (!dir.exists(sub_dir)) {
      if (sp == "mouse") {
        stop("必需的子目录不存在: ", sub_dir)
      } else {
        cat("  WARN:", sp, "子目录不存在,跳过\n")
        next
      }
    }
    rds_files <- list.files(sub_dir, pattern = "annotated\\.rds$",
                            full.names = TRUE)
    rds_files <- rds_files[!grepl("_v4\\.rds$", rds_files)]
    if (length(rds_files) == 0) {
      if (sp == "mouse") {
        stop("在 ", sub_dir, " 下没找到 *annotated.rds")
      } else {
        cat("  WARN:", sub_dir, "下没 annotated.rds,跳过\n")
        next
      }
    }
    latest <- rds_files[which.max(file.mtime(rds_files))]
    result[[sp]] <- latest
  }
  return(result)
}

rds_inputs <- find_rds_files(SAMPLE_DIR, SAMPLE, is_pdx_liver)
for (sp in names(rds_inputs)) {
  cat("  [", sp, "]", basename(rds_inputs[[sp]]), "\n")
}
cat("\n")

# ============================================================
# 3. 工具函数
# ============================================================

# 找 v3.2 注释列,优先 final_celltype,要求 unique > 1
V32_CANDIDATES <- c("final_celltype", "celltype", "cell_type",
                    "manual_celltype", "Annotation", "predicted.id")

detect_v32_col <- function(seu) {
  for (col in V32_CANDIDATES) {
    if (col %in% colnames(seu@meta.data)) {
      vals <- unique(as.character(seu@meta.data[[col]]))
      vals <- vals[!is.na(vals) & vals != ""]
      if (length(vals) > 1) return(col)
    }
  }
  return(NA_character_)
}

# 清理上次跑留下的 MS_\d+ 残留列
clean_ms_residuals <- function(seu) {
  ms_cols <- grep("^MS_\\d+$", colnames(seu@meta.data), value = TRUE)
  if (length(ms_cols) > 0) {
    cat("    清理残留 score 列:", length(ms_cols), "个\n")
    for (c in ms_cols) seu[[c]] <- NULL
  }
  seu
}

# 过滤 marker 库: 只保留 RDS 里存在的基因,且每个集合 >= MIN_MARKERS_PER_SET
# 同时把 catch-all 集排除掉 (避免通吃)
filter_marker_list <- function(marker_list, available_genes,
                               catchall_names = character(0)) {
  out <- list()
  skip_log <- character(0)
  for (ct in names(marker_list)) {
    if (ct %in% catchall_names) {
      skip_log <- c(skip_log, sprintf("    SKIP %s: catch-all set", ct))
      next
    }
    genes <- intersect(marker_list[[ct]], available_genes)
    if (length(genes) >= MIN_MARKERS_PER_SET) {
      out[[ct]] <- genes
    } else {
      skip_log <- c(skip_log,
                    sprintf("    SKIP %s: only %d/%d markers",
                            ct, length(genes), length(marker_list[[ct]])))
    }
  }
  list(markers = out, skip_log = skip_log)
}

# 给一个 Seurat 对象一次性算所有 marker 集的 score,返回矩阵
compute_score_matrix <- function(seu, marker_filtered) {
  if (length(marker_filtered) == 0) {
    return(matrix(numeric(0), nrow = ncol(seu), ncol = 0,
                  dimnames = list(colnames(seu), character(0))))
  }
  seu <- AddModuleScore(seu, features = marker_filtered,
                        name = "MS_", seed = 42)
  score_cols <- paste0("MS_", seq_along(marker_filtered))
  mat <- as.matrix(seu@meta.data[, score_cols, drop = FALSE])
  colnames(mat) <- names(marker_filtered)
  # 清掉临时列
  for (c in score_cols) seu[[c]] <- NULL
  return(mat)
}

# ============================================================
# 4. 核心: 对一个 Seurat 对象做二次注释
# ============================================================
annotate_one <- function(seu, marker_list, alias_tbl, catchall_names,
                         organ_label, sp_label, sample_name,
                         refined_lineage_map = NULL,
                         v32_lineage_map = NULL) {

  # --- 4.1 检测 v3.2 列 ---
  v32_col <- detect_v32_col(seu)
  if (is.na(v32_col)) {
    stop("找不到有效的 v3.2 注释列 (候选: ",
         paste(V32_CANDIDATES, collapse = ", "), ")")
  }
  cat("  [4.1] v3.2 列: `", v32_col, "` (",
      length(unique(seu@meta.data[[v32_col]])), " types)\n", sep = "")
  seu$celltype_v32 <- as.character(seu@meta.data[[v32_col]])

  # --- 4.2 清残留 ---
  seu <- clean_ms_residuals(seu)

  # --- 4.3 过滤 marker 库 ---
  cat("  [4.2] 过滤 marker 库...\n")
  fil <- filter_marker_list(marker_list, rownames(seu), catchall_names)
  marker_filtered <- fil$markers
  if (length(fil$skip_log) > 0) {
    cat(paste(fil$skip_log, collapse = "\n"), "\n")
  }
  cat("    Kept", length(marker_filtered), "/", length(marker_list),
      "cell types (catch-all excluded:",
      length(intersect(names(marker_list), catchall_names)), ")\n")

  if (length(marker_filtered) == 0) {
    stop("marker 过滤后为空,marker 库与物种不匹配?")
  }

  # --- 4.4 一次性计算全局 score matrix (给所有细胞所有集合打分) ---
  cat("  [4.3] 计算全局 score matrix...\n")
  score_mat <- compute_score_matrix(seu, marker_filtered)
  cat("    score matrix:", nrow(score_mat), "cells x ",
      ncol(score_mat), "sets\n")

  # --- 4.5 初始化输出列 ---
  n_cells <- ncol(seu)
  cell_ids <- colnames(seu)
  refined <- character(n_cells)
  source_tag <- character(n_cells)
  v32_rank <- rep(NA_integer_, n_cells)
  names(refined) <- cell_ids
  names(source_tag) <- cell_ids
  names(v32_rank) <- cell_ids

  v32_names <- seu$celltype_v32

  # --- 4.6 对每个 v3.2 类别分路径处理 ---
  cat("  [4.4] 分路径处理 v3.2 类别...\n")

  unique_v32 <- unique(v32_names)
  marker_names_set <- names(marker_filtered)

  for (v32_name in unique_v32) {
    cell_mask <- v32_names == v32_name
    n_this <- sum(cell_mask)
    if (n_this == 0) next

    # --- 路径 A: 精修子空间 (只限 refinement_rules 里的) ---
    if (v32_name %in% names(refinement_rules)) {
      sub_targets <- intersect(refinement_rules[[v32_name]], marker_names_set)
      if (length(sub_targets) >= 2) {
        sub_mat <- score_mat[cell_mask, sub_targets, drop = FALSE]
        picked <- sub_targets[apply(sub_mat, 1, which.max)]
        refined[cell_mask] <- picked
        source_tag[cell_mask] <- "subtype_rescore"
        # rank = 在子空间里的位置 (总是 1 因为是 winner)
        v32_rank[cell_mask] <- 1L
        cat("    [rescore] ", v32_name, " (", n_this, " cells) -> ",
            paste(unique(picked), collapse = "/"), "\n", sep = "")
        next
      } else {
        cat("    WARN: refinement rule for ", v32_name,
            " has <2 available targets,fallback to passthrough\n", sep = "")
      }
    }

    # --- 路径 B: alias 改名 ---
    if (v32_name %in% names(alias_tbl)) {
      new_name <- alias_tbl[[v32_name]]
      refined[cell_mask] <- new_name
      source_tag[cell_mask] <- "alias"
      # 如果 alias 后的名字在 marker 库里,也记一下 rank
      if (new_name %in% marker_names_set) {
        sub_scores <- score_mat[cell_mask, , drop = FALSE]
        ranks <- apply(sub_scores, 1, function(x) {
          rank(-x, ties.method = "min")[new_name]
        })
        v32_rank[cell_mask] <- as.integer(ranks)
      }
      next
    }

    # --- 路径 C: passthrough (保留原名,不打分) ---
    if (v32_name %in% passthrough_names) {
      refined[cell_mask] <- v32_name
      source_tag[cell_mask] <- "passthrough"
      next
    }

    # --- 路径 D: score 验证 (v3.2 名直接在 marker 库里) ---
    if (v32_name %in% marker_names_set) {
      # ---- lineage gate (仅脑启用) ----
      # 如果 v32_name 有 lineage 映射, 把 score_mat 限制在同 lineage 兄弟类
      # 这样 rank / top1 都只在 lineage 内取, flip 不会跨谱系
      if (!is.null(v32_lineage_map) && v32_name %in% names(v32_lineage_map)) {
        my_lineage <- v32_lineage_map[[v32_name]]
        gated_cols <- names(refined_lineage_map)[
          sapply(refined_lineage_map, function(x) x == my_lineage)
        ]
        gated_cols <- intersect(gated_cols, marker_names_set)
        if (length(gated_cols) >= 2) {
          sub_scores <- score_mat[cell_mask, gated_cols, drop = FALSE]
        } else {
          # lineage 内只有 1 个或 0 个集合 -> 无法 gate, 保留 v3.2 不 flip
          refined[cell_mask] <- v32_name
          source_tag[cell_mask] <- "rescore_agree"
          cat("    [gate-only] ", v32_name, " (", n_this,
              " cells): lineage `", my_lineage,
              "` has <2 members, kept as-is\n", sep = "")
          next
        }
      } else {
        sub_scores <- score_mat[cell_mask, , drop = FALSE]
      }

      # v3.2 名的 rank (在 sub_scores 里)
      ranks <- apply(sub_scores, 1, function(x) {
        rank(-x, ties.method = "min")[v32_name]
      })
      ranks <- as.integer(ranks)
      v32_rank[cell_mask] <- ranks
      # top3 内 -> agree -> 保留 v3.2 名
      # top3 外 -> flipped -> 用 score top1 (已限定在 lineage 内)
      top1 <- colnames(sub_scores)[apply(sub_scores, 1, which.max)]
      agree_mask <- ranks <= RESCORE_TOPK
      # 写回
      cell_idx <- which(cell_mask)
      refined[cell_idx[agree_mask]]  <- v32_name
      source_tag[cell_idx[agree_mask]] <- "rescore_agree"
      refined[cell_idx[!agree_mask]] <- top1[!agree_mask]
      source_tag[cell_idx[!agree_mask]] <- "rescore_flipped"
      n_flipped <- sum(!agree_mask)
      if (n_flipped > 0) {
        cat("    [verify] ", v32_name, ": ", sum(agree_mask),
            " agree, ", n_flipped, " flipped\n", sep = "")
      }
      next
    }

    # --- 路径 E: 未知名字,保留但标记 ---
    refined[cell_mask] <- v32_name
    source_tag[cell_mask] <- "no_markers"
    cat("    [no_markers] ", v32_name, " (", n_this,
        " cells): 不在 alias/refinement/marker 库,保留原名\n", sep = "")
  }

  # --- 4.7 计算 score margin 做 confidence ---
  cat("  [4.5] 计算 confidence...\n")
  # 对所有细胞算 top1-top2 margin (基于全局 score matrix)
  margins <- apply(score_mat, 1, function(x) {
    s <- sort(x, decreasing = TRUE)
    s[1] - s[2]
  })

  refined_sizes <- table(refined)
  confidence <- character(n_cells)
  for (i in seq_len(n_cells)) {
    n <- refined_sizes[refined[i]]
    src <- source_tag[i]
    # 精修和 alias 不经过 margin 判断,用群体大小
    if (src %in% c("subtype_rescore", "alias")) {
      confidence[i] <- if (n >= 100) "high" else if (n >= 30) "medium" else "low"
    } else if (src == "passthrough" || src == "no_markers") {
      confidence[i] <- if (n >= 100) "medium" else "low"  # passthrough 天然降级
    } else {
      # rescore_agree / rescore_flipped: margin 参与判断
      m <- margins[i]
      if (n >= 100 && m > SCORE_MARGIN_FIXED) {
        confidence[i] <- "high"
      } else if (n >= 30) {
        confidence[i] <- "medium"
      } else {
        confidence[i] <- "low"
      }
    }
  }
  cat("    margin threshold (fixed):", SCORE_MARGIN_FIXED, "\n")
  cat("    confidence: high=", sum(confidence == "high"),
      " medium=", sum(confidence == "medium"),
      " low=", sum(confidence == "low"), "\n", sep = "")

  # --- 4.8 应用 mapping 生成 broad ---
  cat("  [4.6] 应用 mapping 生成 celltype_broad...\n")

  if (organ_label == "liver_human") {
    map_subset <- mapping_df[mapping_df$organ == "liver_human", ]
  } else if (organ_label == "liver") {
    map_subset <- mapping_df[mapping_df$organ == "liver", ]
  } else {
    map_subset <- mapping_df[mapping_df$organ == "brain", ]
  }
  map_subset <- map_subset[!duplicated(map_subset$refined), ]
  refined_to_broad <- setNames(map_subset$broad, map_subset$refined)

  broad <- refined_to_broad[refined]
  unmapped <- unique(refined[is.na(broad)])
  if (length(unmapped) > 0) {
    cat("    WARN: unmapped refined 名字 (", length(unmapped), "):\n", sep = "")
    for (u in unmapped) cat("      -", u, "\n")
    broad[is.na(broad)] <- "Unknown_NotMapped"
  }

  # --- 4.9 写回 metadata (unname 防 Seurat 按 name 匹配 barcode) ---
  cat("  [4.7] 写回 metadata...\n")
  seu$celltype_refined            <- unname(refined)
  seu$celltype_broad              <- unname(broad)
  seu$celltype_refinement_source  <- unname(source_tag)
  seu$celltype_v32_rank           <- unname(v32_rank)
  seu$celltype_confidence         <- unname(confidence)

  list(
    seurat = seu,
    v32_col = v32_col,
    score_matrix = score_mat,
    n_marker_used = length(marker_filtered),
    n_marker_total = length(marker_list),
    skip_log = fil$skip_log,
    unmapped = unmapped
  )
}

# ============================================================
# 5. 处理每个 RDS
# ============================================================
cat("[3/6] 开始注释...\n\n")

results_per_species <- list()

for (sp in names(rds_inputs)) {
  rds_path <- rds_inputs[[sp]]
  cat("--- Processing [", sp, "] ---\n", sep = "")
  cat("Loading:", basename(rds_path), "\n")

  seu <- readRDS(rds_path)
  cat("  Cells:", ncol(seu), " Features:", nrow(seu), "\n")

  # 选 marker 库 / alias / catch-all / organ
  refined_lineage_map <- NULL
  v32_lineage_map     <- NULL
  if (is_brain) {
    marker_list    <- brain_markers_v4
    alias_tbl      <- alias_brain
    catchall_names <- catchall_sets_brain
    organ_label    <- "brain"
    # 脑启用 lineage gate
    if (exists("brain_lineage_map"))     refined_lineage_map <- brain_lineage_map
    if (exists("brain_v32_lineage_map")) v32_lineage_map     <- brain_v32_lineage_map
  } else if (sp == "human") {
    marker_list    <- liver_markers_human_v4
    alias_tbl      <- alias_liver_human
    catchall_names <- catchall_sets_liver_human
    organ_label    <- "liver_human"
  } else {
    marker_list    <- liver_markers_mouse_v4
    alias_tbl      <- alias_liver_mouse
    catchall_names <- catchall_sets_liver_mouse
    organ_label    <- "liver"
  }

  ann <- annotate_one(seu, marker_list, alias_tbl, catchall_names,
                      organ_label, sp, SAMPLE,
                      refined_lineage_map = refined_lineage_map,
                      v32_lineage_map     = v32_lineage_map)

  # 保存新 RDS
  out_path <- sub("\\.rds$", "_v4.rds", rds_path)
  saveRDS(ann$seurat, out_path)
  cat("  Saved:", basename(out_path), "\n\n")

  results_per_species[[sp]] <- c(ann, list(
    out_path = out_path,
    in_path = rds_path
  ))
}

# ============================================================
# 6. UMAP 对比图
# ============================================================
cat("[4/6] 生成 UMAP 对比图...\n")

plot_comparison <- function(seu, sp, sample) {
  if (!"umap" %in% names(seu@reductions)) {
    cat("  WARN:", sp, "没有 umap,跳过\n")
    return(NULL)
  }
  mk <- function(grp, title) {
    DimPlot(seu, group.by = grp, label = TRUE, repel = TRUE,
            label.size = 2.5) +
      ggtitle(paste0(sample, " [", sp, "] - ", title)) +
      theme(legend.position = "bottom",
            legend.text = element_text(size = 6),
            plot.title = element_text(size = 9))
  }
  p1 <- mk("celltype_v32",     "v3.2 (original)")
  p2 <- mk("celltype_refined", "v4 refined")
  p3 <- mk("celltype_broad",   "v4 broad")
  p4 <- DimPlot(seu, group.by = "celltype_refinement_source",
                cols = c("alias" = "#1f77b4",
                         "rescore_agree" = "#2ca02c",
                         "rescore_flipped" = "#d62728",
                         "subtype_rescore" = "#9467bd",
                         "passthrough" = "#7f7f7f",
                         "no_markers" = "#e377c2")) +
    ggtitle(paste0(sample, " [", sp, "] - refinement source")) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 6),
          plot.title = element_text(size = 9))
  (p1 | p2) / (p3 | p4)
}

umap_pdf_path <- file.path(SAMPLE_DIR, paste0(SAMPLE, "_v4_umap.pdf"))
pdf(umap_pdf_path, width = 14, height = 12)
for (sp in names(results_per_species)) {
  p <- plot_comparison(results_per_species[[sp]]$seurat, sp, SAMPLE)
  if (!is.null(p)) print(p)
}
dev.off()
cat("  Saved:", basename(umap_pdf_path), "\n\n")

# ============================================================
# 7. Markdown 报告
# ============================================================
cat("[5/6] 生成 markdown report...\n")

report_path <- file.path(SAMPLE_DIR, paste0(SAMPLE, "_v4_annotation_report.md"))
lines <- character(0)
add <- function(...) lines <<- c(lines, paste0(...))

add("# Secondary Annotation v4 — ", SAMPLE)
add("")
add("**Generated:** ", format(Sys.time()))
add("")
add("**Strategy:** A + option 3 + rescore top-", RESCORE_TOPK)
add("")
add("- Trust v3.2 `final_celltype` as the base annotation")
add("- Immune coarse classes get subtype rescoring ",
    "(DC/T/B in liver, Microglia in brain)")
add("- Other v3.2 classes: alias normalization -> score verification ",
    "-> passthrough")
add("- `celltype_v32` kept in metadata as the original annotation")
add("")
add("**Marker library sources:**")
if (is_brain) {
  for (s in brain_marker_sources) add("- ", s)
} else {
  for (s in liver_marker_sources) add("- ", s)
}
add("")
add("**Fixed parameters:**")
add("- score margin threshold (high confidence): ", SCORE_MARGIN_FIXED)
add("- rescore top-K (agree if v3.2 name in score top-K): ", RESCORE_TOPK)
add("- min markers per set: ", MIN_MARKERS_PER_SET)
add("")
add("---")
add("")

for (sp in names(results_per_species)) {
  res <- results_per_species[[sp]]
  seu <- res$seurat
  add("## [", sp, "]")
  add("")
  add("- Input:  `", basename(res$in_path), "`")
  add("- Output: `", basename(res$out_path), "`")
  add("- Cells: ", ncol(seu))
  add("- v3.2 column detected: `", res$v32_col, "`")
  add("- Markers used: ", res$n_marker_used, " / ", res$n_marker_total)
  add("")

  # --- Refinement source 分布 ---
  add("### Refinement source")
  add("")
  add("| Source | N cells | % |")
  add("|---|---:|---:|")
  src_tab <- table(seu$celltype_refinement_source)
  src_order <- c("rescore_agree", "rescore_flipped", "subtype_rescore",
                 "alias", "passthrough", "no_markers")
  for (s in src_order) {
    n <- if (s %in% names(src_tab)) src_tab[[s]] else 0
    pct <- round(n / ncol(seu) * 100, 2)
    add("| ", s, " | ", n, " | ", pct, "% |")
  }
  add("")

  # --- Confidence 分布 ---
  add("### Confidence distribution")
  add("")
  add("| Confidence | N | % |")
  add("|---|---:|---:|")
  conf_tab <- table(seu$celltype_confidence)
  for (lv in c("high", "medium", "low")) {
    n <- if (lv %in% names(conf_tab)) conf_tab[[lv]] else 0
    pct <- round(n / ncol(seu) * 100, 2)
    add("| ", lv, " | ", n, " | ", pct, "% |")
  }
  add("")

  # --- v3.2 vs v4 flipped 细胞列表 (重点审查目标) ---
  add("### v3.2 → v4 flipped cells (score verify disagree)")
  add("")
  flipped_mask <- seu$celltype_refinement_source == "rescore_flipped"
  n_flip <- sum(flipped_mask)
  if (n_flip == 0) {
    add("*(no flipped cells - all v3.2 rescore classes agreed within top-",
        RESCORE_TOPK, ")*")
  } else {
    add("**Total flipped: ", n_flip, " cells (",
        round(n_flip / ncol(seu) * 100, 2), "%)**")
    add("")
    flip_tab <- table(seu$celltype_v32[flipped_mask],
                      seu$celltype_refined[flipped_mask])
    add("| v3.2 name | v4 new name | N |")
    add("|---|---|---:|")
    ft_df <- as.data.frame(flip_tab, stringsAsFactors = FALSE)
    ft_df <- ft_df[ft_df$Freq > 0, ]
    ft_df <- ft_df[order(-ft_df$Freq), ]
    for (i in seq_len(nrow(ft_df))) {
      add("| `", ft_df$Var1[i], "` | `", ft_df$Var2[i], "` | ", ft_df$Freq[i], " |")
    }
  }
  add("")

  # --- v3.2 vs v4 refined 完整分布 ---
  add("### v3.2 → v4 refined distribution")
  add("")
  add("**v3.2:**")
  add("")
  add("| Cell type | N | % |")
  add("|---|---:|---:|")
  v32_tab <- sort(table(seu$celltype_v32), decreasing = TRUE)
  for (ct in names(v32_tab)) {
    n <- v32_tab[[ct]]
    add("| ", ct, " | ", n, " | ", round(n/ncol(seu)*100, 2), "% |")
  }
  add("")
  add("**v4 refined:**")
  add("")
  add("| Cell type | N | % | Main src | Main conf |")
  add("|---|---:|---:|---|---|")
  v4_tab <- sort(table(seu$celltype_refined), decreasing = TRUE)
  src_x_ct <- table(seu$celltype_refined, seu$celltype_refinement_source)
  conf_x_ct <- table(seu$celltype_refined, seu$celltype_confidence)
  for (ct in names(v4_tab)) {
    n <- v4_tab[[ct]]
    pct <- round(n/ncol(seu)*100, 2)
    main_src  <- names(which.max(src_x_ct[ct, ]))
    main_conf <- names(which.max(conf_x_ct[ct, ]))
    flag <- ifelse(main_conf == "high", "✓",
            ifelse(main_conf == "medium", "◐", "⚠"))
    add("| ", ct, " | ", n, " | ", pct, "% | ", main_src,
        " | ", main_conf, " ", flag, " |")
  }
  add("")

  add("**v4 broad:**")
  add("")
  add("| Broad class | N | % |")
  add("|---|---:|---:|")
  vb_tab <- sort(table(seu$celltype_broad), decreasing = TRUE)
  for (ct in names(vb_tab)) {
    n <- vb_tab[[ct]]
    add("| ", ct, " | ", n, " | ", round(n/ncol(seu)*100, 2), "% |")
  }
  add("")

  # --- Skip log (marker 过滤) ---
  if (length(res$skip_log) > 0) {
    add("### Marker filter log")
    add("")
    add("```")
    for (l in res$skip_log) add(sub("^\\s+", "", l))
    add("```")
    add("")
  }

  # --- Unmapped ---
  if (length(res$unmapped) > 0) {
    add("### :warning: Unmapped refined names (broad = Unknown_NotMapped)")
    add("")
    for (u in res$unmapped) add("- `", u, "`")
    add("")
  }

  add("---")
  add("")
}

add("## Notes")
add("")
add("- `celltype_refinement_source`:")
add("  - `rescore_agree`: v3.2 name is in top-", RESCORE_TOPK,
    " of marker scores")
add("  - `rescore_flipped`: v3.2 name fell outside top-", RESCORE_TOPK,
    ", flipped to score top1 (review these!)")
add("  - `subtype_rescore`: v3.2 coarse immune class refined to subtype")
add("  - `alias`: v3.2 name normalized to v4 standard (no rescoring)")
add("  - `passthrough`: kept as-is (v4 lib doesn't have markers for this)")
add("  - `no_markers`: unexpected v3.2 name not in alias/refinement/lib")
add("- `celltype_v32_rank`: rank of v3.2 name in the cell's marker score ",
    "vector (1 = top), NA if v3.2 name not in marker lib")
add("- Flipped cells are the primary review target. ",
    "High flipped count in one v3.2 class suggests v4 lib and v3.2 disagree ",
    "on cell type definition; manual verification recommended.")

writeLines(lines, report_path)
cat("  Saved:", basename(report_path), "\n\n")

# ============================================================
# 8. 完成
# ============================================================
cat("[6/6] 完成!\n\n")
cat("==========================================\n")
cat("Output (in ", SAMPLE_DIR, "):\n", sep = "")
for (sp in names(results_per_species)) {
  cat("  - ", sp, "/",
      basename(results_per_species[[sp]]$out_path), "\n", sep = "")
}
cat("  - ", basename(umap_pdf_path), "\n", sep = "")
cat("  - ", basename(report_path), "\n", sep = "")
cat("==========================================\n")
