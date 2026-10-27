#!/usr/bin/env Rscript
# 只跑 brain markers (liver 4/30 已完成)
# 关键修复: future.globals.maxSize 16GB → 50GB

# 把原脚本 source 进来 (但不执行底部主循环)
# 用 trick: 只跑到 process_organ 定义结束
script <- readLines("~/SGS10_analysis/scripts/microenv_marker_genes_RCTD100.R")
# 找到 "results <- list()" 这行, 只 source 之前的部分
end_line <- grep("^results <- list\\(\\)", script)[1]
if (is.na(end_line)) stop("Cannot find 'results <- list()' line")
cat(sprintf("Sourcing lines 1-%d (definitions only)\n", end_line - 1))

# 写入临时文件
tmp <- tempfile(fileext = ".R")
writeLines(script[1:(end_line - 1)], tmp)
source(tmp)

cat("\n========== Running BRAIN ONLY ==========\n")
brain_result <- process_organ("brain")

if (is.null(brain_result)) {
  cat("⚠ Brain markers FAILED\n")
} else {
  cat(sprintf("✓ Brain markers DONE: %d microenvs, %d markers\n",
              brain_result$n_microenv, brain_result$n_markers))
}
