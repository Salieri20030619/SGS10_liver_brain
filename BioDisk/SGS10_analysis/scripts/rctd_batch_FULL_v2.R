#!/usr/bin/env Rscript
# ============================================================
# rctd_batch_FULL_v2.R — 安全版 (max_cores=2)
#
# 修复:
#   v1 (rctd_batch_10samples_FULL.R) 用 max_cores=8 (4 worker × 2 OpenMP)
#   导致 worker 死锁 (主进程等已死的 worker, 4-5 hr 无产出)
#
# v2 改动:
#   - MAX_CORES: 8 → 2 (单进程, 无 RSOCK worker, 不会死锁)
#   - SUBSAMPLE_RATE: 0.2 → 1.0 (FULL 不抽样)
#   - OUT_DIR: spatial_results/rctd → spatial_results/rctd_full
#   - DOUBLET_MODE: "full" 保留 (与原始一致, paper standard)
#
# 时间预估:
#   max_cores=2 串行 + 100% spots:
#     单 sample 1.5-3 hr (vs 多核 1.5-2 hr 但会死锁)
#     总 15-30 hr (今晚启动, 明早/中午看完成)
#
# 用法:
#   cd ~/SGS10_analysis
#   nohup nice -n 19 Rscript scripts/rctd_batch_FULL_v2.R \
#       > scripts/logs/rctd_full_v2_$(date +%Y%m%d_%H%M%S).log 2>&1 &
#   tail -f scripts/logs/rctd_full_v2_*.log
#
# 中断恢复:
#   原脚本有 skip-if-exists, 重启自动跳过已完成 sample
# ============================================================

cat("============================================================\n")
cat("RCTD FULL v2 (max_cores=2 safe mode)\n")
cat("Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("============================================================\n\n")

# ---------- 路径 ----------
PROJECT_ROOT <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
NEW_OUT_DIR <- file.path(PROJECT_ROOT, "spatial_results", "rctd_full")
ORIG_SCRIPT <- file.path(PROJECT_ROOT, "scripts", "rctd_batch_10samples.R")

# ---------- 创建输出目录 ----------
dir.create(NEW_OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(NEW_OUT_DIR, "per_sample"), showWarnings = FALSE)
dir.create(file.path(NEW_OUT_DIR, "enriched"), showWarnings = FALSE)
dir.create(file.path(NEW_OUT_DIR, "summary"), showWarnings = FALSE)

# 写 readme 说明
LOG_FILE_INIT <- file.path(NEW_OUT_DIR, "summary", "FULL_v2_README.txt")
writeLines(c(
  "=== RCTD FULL v2 (max_cores=2 safe mode) ===",
  "",
  paste0("Started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  paste0("OUT_DIR:           ", NEW_OUT_DIR),
  "SUBSAMPLE_RATE:    1.0  (vs 0.2 in original)",
  "RCTD_DOUBLET_MODE: full (same as original)",
  "MAX_CORES:         2    (vs 8 in v1, avoid 4-worker deadlock)",
  "",
  "Original 20% subsampled data preserved at:",
  "  spatial_results/rctd/",
  "",
  "Skip-if-exists logic preserved.",
  "Est total: 15-30 hours serial."
), LOG_FILE_INIT)

# ---------- 检查现有 rctd_full 已完成 sample ----------
existing <- list.files(file.path(NEW_OUT_DIR, "per_sample"),
                        pattern = "_rctd_results\\.rds$",
                        full.names = FALSE)
if (length(existing) > 0) {
  cat("Existing completed samples (will skip):\n")
  cat(paste("  ", sub("_rctd_results.rds", "", existing), collapse="\n"), "\n\n")
}

# ---------- 检查原脚本存在 ----------
if (!file.exists(ORIG_SCRIPT)) {
  stop("原脚本未找到: ", ORIG_SCRIPT)
}
cat("Reading orig script:", ORIG_SCRIPT, "\n")

orig_lines <- readLines(ORIG_SCRIPT)
cat("Orig script lines:", length(orig_lines), "\n\n")

# ---------- 改 1: OUT_DIR ----------
out_dir_idx <- grep('^OUT_DIR <- file.path\\(PROJECT_ROOT', orig_lines)
if (length(out_dir_idx) != 1) {
  stop("无法定位 OUT_DIR 行, 找到 ", length(out_dir_idx), " 处")
}
orig_lines[out_dir_idx] <- 
  'OUT_DIR <- file.path(PROJECT_ROOT, "spatial_results", "rctd_full")'
cat("[1/3] Replaced OUT_DIR at line", out_dir_idx, "\n")

# ---------- 改 2: SUBSAMPLE_RATE ----------
sub_idx <- grep('^SUBSAMPLE_RATE <- ', orig_lines)
if (length(sub_idx) != 1) {
  stop("无法定位 SUBSAMPLE_RATE 行")
}
orig_lines[sub_idx] <- 
  'SUBSAMPLE_RATE <- 1.0   # FULL run v2 (max_cores=2 safe mode)'
cat("[2/3] Replaced SUBSAMPLE_RATE at line", sub_idx, "\n")

# ---------- 改 3: MAX_CORES (关键!) ----------
# 原始可能写法多种, 找:
#   MAX_CORES <- 8
#   MAX_CORES <- ...
#   max_cores = 8
#   max_cores=8 in run.RCTD()
core_idx <- grep('^MAX_CORES <- ', orig_lines)
if (length(core_idx) >= 1) {
  for (i in core_idx) {
    orig_lines[i] <- 'MAX_CORES <- 2   # FULL v2 safe mode'
  }
  cat("[3/3] Replaced MAX_CORES at line(s):", paste(core_idx, collapse=","), "\n")
} else {
  # Fallback: 在 create.RCTD() 或 run.RCTD() 中找 max_cores=
  rctd_lines_idx <- grep('max_cores\\s*=\\s*[0-9]+', orig_lines)
  if (length(rctd_lines_idx) >= 1) {
    for (i in rctd_lines_idx) {
      orig_lines[i] <- gsub('max_cores\\s*=\\s*[0-9]+',
                            'max_cores = 2',
                            orig_lines[i])
    }
    cat("[3/3] Replaced max_cores= at line(s):", paste(rctd_lines_idx, collapse=","),
        "\n")
  } else {
    cat("⚠ 未找到 MAX_CORES 或 max_cores=, 检查 RCTD 是否有这参数\n")
    cat("  RCTD 可能默认 max_cores=8, 这会再次死锁!\n")
    cat("  脚本继续, 但风险存在\n\n")
  }
}

# ---------- 写到临时文件 + source ----------
TMP_SCRIPT <- tempfile(fileext = "_FULL_v2.R")
writeLines(orig_lines, TMP_SCRIPT)
cat("\nModified script written to:", TMP_SCRIPT, "\n")

# 写一份到 logs 目录方便调试
DEBUG_COPY <- file.path(PROJECT_ROOT, "scripts", "logs",
                        paste0("rctd_FULL_v2_modified_",
                               format(Sys.time(), "%Y%m%d_%H%M%S"), ".R"))
file.copy(TMP_SCRIPT, DEBUG_COPY, overwrite = TRUE)
cat("Debug copy:", DEBUG_COPY, "\n\n")

cat("============================================================\n")
cat("Sourcing modified RCTD script...\n")
cat("============================================================\n\n")

source(TMP_SCRIPT)

cat("\n============================================================\n")
cat("RCTD FULL v2 batch complete.\n")
cat("Output:", NEW_OUT_DIR, "\n")
cat("Done time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("============================================================\n")
