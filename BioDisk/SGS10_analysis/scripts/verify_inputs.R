#!/usr/bin/env Rscript
# ============================================================
# verify_inputs.R  --  二次注释 dry-run 前的输入验证
# ============================================================
# 用法:
#   Rscript verify_inputs.R
#
# 目的 (只读, 不修改任何文件):
#   1) 列出每个样本目录下的实际文件结构 (mouse/ human/ 子目录是否存在)
#   2) 找到所有候选 *annotated.rds 文件并打印路径
#   3) 对每个候选 RDS, 加载后报告:
#        - ncol / nrow
#        - meta.data 列名 (重点找 v3.2 注释列)
#        - 候选注释列的 unique 值 (前 30 个)
#        - 是否已有 reductions$umap
#        - 是否已有残留的 MS_* 列 (上次 v4 跑的痕迹)
#   4) 检查 mapping CSV 的 refined 列 vs marker 库的 cell type 名是否一致
#        - mapping 里有但 marker 库里没有 (死代码)
#        - marker 库里有但 mapping 里没有 (会落到 Unknown_NotMapped)
#   5) 检查 marker 库内的基因共享情况 (诊断 winner-takes-all 偏差)
#
# 输出:
#   stdout: 详细日志
#   ~/SGS10_analysis/scripts/verify_inputs_report.md  汇总报告
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
})

# ------------------------------------------------------------
# 0. 路径配置  (和 secondary_annotation.R 保持一致)
# ------------------------------------------------------------
PROJECT_ROOT <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
SCRIPT_DIR   <- file.path(PROJECT_ROOT, "scripts")

PDX_LIVER_SAMPLES  <- c("liver8", "gan56", "liver3")
PURE_LIVER_SAMPLES <- c("con_liver", "bla_liver")
BRAIN_SAMPLES      <- c("con_brain", "bla_brain", "nao56", "brain3", "brain8")
ALL_SAMPLES        <- c(PURE_LIVER_SAMPLES, PDX_LIVER_SAMPLES, BRAIN_SAMPLES)

REPORT_PATH <- file.path(SCRIPT_DIR, "verify_inputs_report.md")

# 收集报告内容
report <- character(0)
add <- function(...) {
  line <- paste0(...)
  report <<- c(report, line)
  cat(line, "\n", sep = "")
}

add("# Secondary Annotation v4 — Input Verification Report")
add("")
add("**Generated:** ", format(Sys.time()))
add("")
add("**Purpose:** dry-run check, no files modified.")
add("")
add("---")
add("")

# ------------------------------------------------------------
# 1. 检查每个样本的目录结构
# ------------------------------------------------------------
add("## 1. Sample directory structure")
add("")

# 候选 v3.2 注释列名 (按优先级)
V32_CANDIDATE_COLS <- c(
  "celltype", "cell_type", "Celltype", "CellType",
  "celltype_v32", "celltype_main", "celltype_final",
  "Annotation", "annotation", "predicted.id",
  "seurat_annotations", "manual_annotation"
)

# 收集每个 RDS 的诊断
all_rds_diag <- list()

for (sample in ALL_SAMPLES) {
  sample_dir <- file.path(PROJECT_ROOT, sample)
  add("### `", sample, "`")
  add("")

  if (!dir.exists(sample_dir)) {
    add("- :x: **Sample directory missing:** `", sample_dir, "`")
    add("")
    next
  }

  # 列出顶层内容
  top_entries <- list.files(sample_dir, full.names = FALSE,
                            include.dirs = TRUE, no.. = TRUE)
  add("- Top-level entries: ", paste(top_entries, collapse = ", "))

  has_mouse_dir <- dir.exists(file.path(sample_dir, "mouse"))
  has_human_dir <- dir.exists(file.path(sample_dir, "human"))
  add("- `mouse/` subdir: ", ifelse(has_mouse_dir, ":white_check_mark:", ":x:"))
  add("- `human/` subdir: ", ifelse(has_human_dir, ":white_check_mark:", ":x:"))
  add("")

  # 找候选 rds:
  #   方案 A (期望的新结构): <sample>/mouse/*.rds 和 <sample>/human/*.rds
  #   方案 B (旧结构遗留):    <sample>/*.rds
  candidates <- list()

  if (has_mouse_dir) {
    rds <- list.files(file.path(sample_dir, "mouse"),
                      pattern = "\\.rds$", full.names = TRUE)
    rds <- rds[!grepl("_v4\\.rds$", rds)]
    if (length(rds) > 0) candidates[["mouse"]] <- rds
  }
  if (has_human_dir) {
    rds <- list.files(file.path(sample_dir, "human"),
                      pattern = "\\.rds$", full.names = TRUE)
    rds <- rds[!grepl("_v4\\.rds$", rds)]
    if (length(rds) > 0) candidates[["human"]] <- rds
  }
  toplevel_rds <- list.files(sample_dir, pattern = "\\.rds$", full.names = TRUE)
  toplevel_rds <- toplevel_rds[!grepl("_v4\\.rds$", toplevel_rds)]
  if (length(toplevel_rds) > 0) candidates[["TOPLEVEL_LEGACY"]] <- toplevel_rds

  if (length(candidates) == 0) {
    add("- :warning: **No candidate .rds files found.**")
    add("")
    next
  }

  add("**Candidate RDS files:**")
  add("")
  for (loc in names(candidates)) {
    for (f in candidates[[loc]]) {
      sz_mb <- round(file.size(f) / 1024^2, 1)
      mtime <- format(file.mtime(f), "%Y-%m-%d %H:%M")
      add("- [`", loc, "`] `", basename(f), "` (", sz_mb, " MB, ", mtime, ")")
    }
  }
  add("")

  # 对每个候选 RDS 做诊断 (只取每个 location 下最新的)
  for (loc in names(candidates)) {
    files <- candidates[[loc]]
    latest <- files[which.max(file.mtime(files))]

    add("**Diagnostic for `", basename(latest), "` [", loc, "]:**")
    add("")

    diag <- tryCatch({
      seu <- readRDS(latest)
      if (!inherits(seu, "Seurat")) {
        list(error = paste0("Not a Seurat object, class = ",
                            paste(class(seu), collapse = "/")))
      } else {
        # 找 v3.2 注释列
        meta_cols <- colnames(seu@meta.data)
        v32_hit <- intersect(V32_CANDIDATE_COLS, meta_cols)
        v32_col <- if (length(v32_hit) > 0) v32_hit[1] else NA

        # 列出 unique values
        v32_vals <- if (!is.na(v32_col)) {
          uvals <- unique(as.character(seu@meta.data[[v32_col]]))
          uvals[!is.na(uvals)]
        } else character(0)

        # 残留 MS_ 列
        ms_residual <- grep("^MS_\\d+$", meta_cols, value = TRUE)

        list(
          n_cells = ncol(seu),
          n_features = nrow(seu),
          meta_cols = meta_cols,
          v32_col = v32_col,
          v32_hit_all = v32_hit,
          v32_vals = v32_vals,
          has_umap = "umap" %in% names(seu@reductions),
          reductions = names(seu@reductions),
          assays = names(seu@assays),
          default_assay = DefaultAssay(seu),
          ms_residual = ms_residual,
          error = NULL
        )
      }
    }, error = function(e) {
      list(error = conditionMessage(e))
    })

    if (!is.null(diag$error)) {
      add("- :x: **Failed to load:** ", diag$error)
      add("")
      next
    }

    add("- Cells / features: **", diag$n_cells, "** / ", diag$n_features)
    add("- Default assay: `", diag$default_assay, "`  | assays: ",
        paste0("`", diag$assays, "`", collapse = ", "))
    add("- Reductions: ", ifelse(length(diag$reductions) > 0,
                                 paste0("`", diag$reductions, "`",
                                        collapse = ", "),
                                 "*(none)*"))
    add("- UMAP available: ",
        ifelse(diag$has_umap, ":white_check_mark:", ":x:"))
    add("- meta.data columns (", length(diag$meta_cols), "): ",
        paste0("`", diag$meta_cols, "`", collapse = ", "))

    if (length(diag$v32_hit_all) == 0) {
      add("- :warning: **NO v3.2 annotation column found** ",
          "(searched: ", paste(V32_CANDIDATE_COLS, collapse = ", "), ")")
    } else {
      add("- v3.2 annotation column matches: ",
          paste0("`", diag$v32_hit_all, "`", collapse = ", "))
      add("- **Selected:** `", diag$v32_col, "`")
      add("- Unique values (", length(diag$v32_vals), "): ",
          paste0("`", head(diag$v32_vals, 30), "`", collapse = ", "),
          ifelse(length(diag$v32_vals) > 30, " ...", ""))
    }

    if (length(diag$ms_residual) > 0) {
      add("- :warning: **Residual MS_ columns from previous v4 run:** ",
          paste0("`", diag$ms_residual, "`", collapse = ", "))
    }

    add("")

    # 存起来给后面的全局检查用
    key <- paste0(sample, "::", loc)
    all_rds_diag[[key]] <- list(
      sample = sample, loc = loc, path = latest, diag = diag
    )
  }
}

add("---")
add("")

# ------------------------------------------------------------
# 2. v3.2 注释列汇总 (跨样本一致性)
# ------------------------------------------------------------
add("## 2. v3.2 annotation column consistency across samples")
add("")

v32_summary <- data.frame(
  key      = character(0),
  v32_col  = character(0),
  n_unique = integer(0),
  stringsAsFactors = FALSE
)
for (k in names(all_rds_diag)) {
  d <- all_rds_diag[[k]]$diag
  if (is.null(d$error)) {
    v32_summary <- rbind(v32_summary, data.frame(
      key      = k,
      v32_col  = ifelse(is.na(d$v32_col), "(none)", d$v32_col),
      n_unique = length(d$v32_vals),
      stringsAsFactors = FALSE
    ))
  }
}

if (nrow(v32_summary) == 0) {
  add("*(no successfully-loaded RDS to summarize)*")
} else {
  add("| Sample::location | v3.2 column | n unique types |")
  add("|---|---|---|")
  for (i in seq_len(nrow(v32_summary))) {
    add("| `", v32_summary$key[i], "` | `", v32_summary$v32_col[i],
        "` | ", v32_summary$n_unique[i], " |")
  }
  add("")
  unique_cols <- unique(v32_summary$v32_col)
  unique_cols <- setdiff(unique_cols, "(none)")
  if (length(unique_cols) > 1) {
    add(":warning: **Inconsistent v3.2 column names across samples:** ",
        paste0("`", unique_cols, "`", collapse = ", "))
    add("(secondary_annotation.R must use a fallback list, not a fixed name)")
  } else if (length(unique_cols) == 1) {
    add(":white_check_mark: All samples use the same v3.2 column: `",
        unique_cols, "`")
  }
}
add("")
add("---")
add("")

# ------------------------------------------------------------
# 3. mapping CSV vs marker 库 一致性检查
# ------------------------------------------------------------
add("## 3. mapping CSV vs marker library consistency")
add("")

mapping_path <- file.path(SCRIPT_DIR, "celltype_mapping_broad.csv")
liver_lib_path <- file.path(SCRIPT_DIR, "marker_lib_liver_v4.R")
brain_lib_path <- file.path(SCRIPT_DIR, "marker_lib_brain_v4.R")

if (!file.exists(mapping_path)) {
  add(":x: mapping CSV not found: `", mapping_path, "`")
} else if (!file.exists(liver_lib_path) || !file.exists(brain_lib_path)) {
  add(":x: marker library files not found in `", SCRIPT_DIR, "`")
} else {
  mapping_df <- read.csv(mapping_path, comment.char = "#",
                         stringsAsFactors = FALSE, strip.white = TRUE)
  mapping_df <- mapping_df[mapping_df$refined != "" &
                           !is.na(mapping_df$refined), ]

  source(liver_lib_path)
  source(brain_lib_path)

  check_consistency <- function(label, marker_names, mapping_subset) {
    add("### ", label)
    add("")
    add("- Marker library cell types: ", length(marker_names))
    add("- Mapping rows: ", nrow(mapping_subset))

    in_lib_not_in_map <- setdiff(marker_names, mapping_subset$refined)
    in_map_not_in_lib <- setdiff(mapping_subset$refined, marker_names)
    dup_in_map <- mapping_subset$refined[
      duplicated(mapping_subset$refined)]

    if (length(in_lib_not_in_map) > 0) {
      add("- :warning: **In marker lib but NOT in mapping** (will become `Unknown_NotMapped`):")
      for (x in in_lib_not_in_map) add("    - `", x, "`")
    } else {
      add("- :white_check_mark: All marker lib types are mapped")
    }

    if (length(in_map_not_in_lib) > 0) {
      add("- :information_source: In mapping but NOT in marker lib (dead rules, harmless):")
      for (x in in_map_not_in_lib) add("    - `", x, "`")
    }

    if (length(dup_in_map) > 0) {
      add("- :warning: **Duplicate `refined` keys in mapping (first wins):**")
      for (x in unique(dup_in_map)) add("    - `", x, "`")
    }
    add("")
  }

  check_consistency("Liver mouse",
                    names(liver_markers_mouse_v4),
                    mapping_df[mapping_df$organ == "liver", ])

  check_consistency("Liver human",
                    names(liver_markers_human_v4),
                    mapping_df[mapping_df$organ == "liver_human", ])

  check_consistency("Brain mouse",
                    names(brain_markers_v4),
                    mapping_df[mapping_df$organ == "brain", ])
}

add("---")
add("")

# ------------------------------------------------------------
# 4. Marker 库内基因共享情况 (winner-takes-all 偏差诊断)
# ------------------------------------------------------------
add("## 4. Marker gene sharing within each library")
add("")
add("Winner-takes-all assignment is biased toward longer / more generic ",
    "marker sets when subsets share genes. This section flags cell types ",
    "whose markers are >=70% contained in a *larger* set in the same library.")
add("")

diagnose_sharing <- function(label, marker_list) {
  add("### ", label)
  add("")
  add("| cell type | n markers | max overlap with larger set | overlap target |")
  add("|---|---:|---:|---|")
  ns <- sapply(marker_list, length)
  for (ct in names(marker_list)) {
    genes <- marker_list[[ct]]
    n <- length(genes)
    larger <- names(marker_list)[ns > n]
    if (length(larger) == 0) {
      add("| `", ct, "` | ", n, " | (largest) | — |")
      next
    }
    overlaps <- sapply(larger, function(other) {
      length(intersect(genes, marker_list[[other]])) / n
    })
    max_idx <- which.max(overlaps)
    max_ov  <- overlaps[max_idx]
    flag <- ""
    if (max_ov >= 0.7) flag <- " :warning:"
    if (max_ov >= 0.9) flag <- " :rotating_light:"
    add("| `", ct, "` | ", n, " | ", round(max_ov * 100), "%", flag,
        " | `", names(overlaps)[max_idx], "` |")
  }
  add("")
}

if (exists("liver_markers_mouse_v4")) {
  diagnose_sharing("Liver mouse", liver_markers_mouse_v4)
}
if (exists("liver_markers_human_v4")) {
  diagnose_sharing("Liver human", liver_markers_human_v4)
}
if (exists("brain_markers_v4")) {
  diagnose_sharing("Brain", brain_markers_v4)
}

add("---")
add("")
add("## 5. Summary checklist")
add("")
add("Before running `secondary_annotation.R` on any sample, the following ",
    "must be true:")
add("")
add("- [ ] Every sample has its expected RDS file under `mouse/` (and `human/` for PDX)")
add("- [ ] No sample relies on the legacy top-level `<sample>/*.rds` location")
add("- [ ] All loaded RDS expose the same v3.2 annotation column name (or ",
    "the main script uses a fallback list)")
add("- [ ] mapping CSV covers every cell type in the v4 marker libraries ",
    "(no `Unknown_NotMapped` warnings expected)")
add("- [ ] marker libraries do not contain large generic 'catch-all' sets ",
    "that subsume specific subtypes (see Section 4 :warning: rows)")
add("")

# ------------------------------------------------------------
# 写报告
# ------------------------------------------------------------
writeLines(report, REPORT_PATH)
cat("\n==========================================\n")
cat("Report written to:\n  ", REPORT_PATH, "\n", sep = "")
cat("==========================================\n")
