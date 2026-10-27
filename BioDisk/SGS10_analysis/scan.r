#!/usr/bin/env Rscript
# 30秒微扫: 确认 D1 niches_by_sample_brain.rds 和 D2 spot_zone_stage_brain.csv 结构
PR <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"

cat("===== 1. niches_by_sample_brain.rds 结构 =====\n")
fp <- file.path(PR, "D1_niche_output/niches_by_sample_brain.rds")
cat("路径:", fp, "\n")
cat("size:", round(file.info(fp)$size / 1024^2, 1), "MB\n\n")

x <- readRDS(fp)
cat("class:", class(x), "\n")
cat("length:", length(x), "\n")
cat("names:", names(x), "\n\n")

if (is.list(x) && length(x) > 0) {
  e <- x[[1]]
  cat("--- 第一个元素 (", names(x)[1], ") ---\n", sep="")
  cat("class:", class(e), "\n")
  if (is.data.frame(e)) {
    cat("dim:", dim(e), "\n")
    cat("colnames:", colnames(e), "\n\n")
    cat("head 3 rows:\n"); print(head(e, 3))
    cat("\nniche column unique values (head):\n")
    if ("niche" %in% colnames(e)) print(head(table(e$niche), 20))
  } else {
    str(e, max.level=2)
  }
}

cat("\n\n===== 2. spot_zone_stage_brain.csv 结构 =====\n")
fp2 <- file.path(PR, "D2_zone_output/spot_zone_stage_brain.csv")
cat("路径:", fp2, "\n")
df <- read.csv(fp2, nrows=5)
cat("colnames:\n")
print(colnames(df))
cat("\nhead 3:\n"); print(head(df, 3))

# 行数
all_lines <- length(readLines(fp2))
cat("\n总行数:", all_lines, "\n")

cat("\n\n===== 3. 5 个 brain samples 在 niches_by_sample 中 =====\n")
if (is.list(x)) {
  for (s in names(x)) {
    e <- x[[s]]
    if (is.data.frame(e) && "niche" %in% colnames(e)) {
      cat(s, ":", nrow(e), " spots, niche distribution:\n", sep="")
      print(table(e$niche))
      cat("\n")
    }
  }
}
