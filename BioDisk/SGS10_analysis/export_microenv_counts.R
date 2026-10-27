suppressPackageStartupMessages(library(dplyr))

DATA_DIR <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis/taskB/spatial_domain_RCTD100"

for (organ in c("liver", "brain")) {
  csv <- file.path(DATA_DIR, sprintf("%s_spot_microenv_assignment.csv", organ))
  df <- read.csv(csv, stringsAsFactors = FALSE)
  # 重映射 group
  df$group <- ifelse(grepl("^Con_", df$sample), "control",
              ifelse(grepl("^Bla_", df$sample), "blank", "tumor"))
  tab <- df |>
    count(sample, group, microenv_name) |>
    arrange(sample, microenv_name)
  cat("\n==========", organ, "==========\n")
  # 紧凑 CSV
  write.csv(tab, sprintf("/tmp/%s_counts.csv", organ), row.names = FALSE)
  # 也直接打印
  for (i in seq_len(nrow(tab))) {
    cat(sprintf("%s,%s,%s,%d\n",
                tab$sample[i], tab$group[i],
                tab$microenv_name[i], tab$n[i]))
  }
}
