# =============================================================
# check_T4_csv_content.R
# 看 T4 输出的 4 种 CSV 文件实际长什么样
# =============================================================

T4_DIR <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis/SGS10_deliverables_2026-04-20/08_spatial_Task_B/T4_spatial_domain"

show_csv <- function(name, n_head = 10) {
  fp <- file.path(T4_DIR, name)
  if (!file.exists(fp)) {
    cat("!! NOT FOUND:", fp, "\n\n"); return(invisible(NULL))
  }
  cat("============================================================\n")
  cat(" File:", name, "\n")
  cat("============================================================\n")
  df <- read.csv(fp, stringsAsFactors = FALSE, check.names = FALSE)
  cat("Dim: ", nrow(df), " rows x ", ncol(df), " cols\n", sep = "")
  cat("Columns:\n"); print(colnames(df))
  cat("\nFirst ", n_head, " rows:\n", sep = "")
  print(head(df, n_head))
  cat("\n")
  invisible(df)
}

# ---- Brain ----
cat("################ BRAIN ################\n\n")
show_csv("brain_domain_ct_profile.csv", 15)
show_csv("brain_domain_naming.csv", 15)
show_csv("brain_per_sample_domain.csv", 15)
show_csv("brain_tumor_vs_baseline.csv", 15)

# ---- Liver ----
cat("\n################ LIVER ################\n\n")
show_csv("liver_domain_ct_profile.csv", 15)
show_csv("liver_domain_naming.csv", 15)
show_csv("liver_per_sample_domain.csv", 15)
show_csv("liver_tumor_vs_baseline.csv", 15)

# ---- Check if per-spot domain exists anywhere ----
cat("\n\n============================================================\n")
cat(" Searching for per-spot domain assignment files\n")
cat("============================================================\n")
system("find /mnt/BioDisk/wangmingyuan/SGS10_analysis -name '*spot*domain*' 2>/dev/null")
system("find /mnt/BioDisk/wangmingyuan/SGS10_analysis -name '*domain*spot*' 2>/dev/null")
system("find /mnt/BioDisk/wangmingyuan/SGS10_analysis -name '*_domain_assignment*' 2>/dev/null")
system("find /mnt/BioDisk/wangmingyuan/SGS10_analysis -name 'per_spot*' 2>/dev/null")
