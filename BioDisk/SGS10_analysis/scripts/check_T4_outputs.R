# =============================================================
# check_T4_outputs.R
# 确认 T4 spatial domain 的文件结构和 per-spot 分配在哪
# =============================================================

suppressPackageStartupMessages(library(Matrix))

# Candidate T4 paths
candidates <- c(
  "~/SGS10_analysis/SGS10_deliverables_2026-04-20/08_spatial_Task_B/T4_spatial_domain",
  "~/SGS10_analysis/spatial_results/T4_spatial_domain",
  "~/SGS10_analysis/spatial_results/t4",
  "/mnt/BioDisk/wangmingyuan/SGS10_analysis/SGS10_deliverables_2026-04-20/08_spatial_Task_B/T4_spatial_domain",
  "/mnt/BioDisk/wangmingyuan/SGS10_analysis/spatial_results/T4_spatial_domain"
)

cat("=== Searching for T4 outputs ===\n")
for (p in candidates) {
  p_expanded <- path.expand(p)
  if (dir.exists(p_expanded)) {
    cat("\nFOUND:", p_expanded, "\n")
    cat("Contents (recursive, 2 levels):\n")
    files <- list.files(p_expanded, recursive = TRUE, full.names = FALSE)
    for (f in head(files, 30)) cat("  ", f, "\n")
    if (length(files) > 30) cat("  ... (", length(files), "files total)\n")
  }
}

cat("\n=== Also search by pattern anywhere ===\n")
system("find /mnt/BioDisk/wangmingyuan/SGS10_analysis -type d -iname '*T4*' 2>/dev/null | head -10")
system("find /mnt/BioDisk/wangmingyuan/SGS10_analysis -type d -iname '*domain*' 2>/dev/null | head -10")
system("find /mnt/BioDisk/wangmingyuan/SGS10_analysis -name 'domain_ct_profile*' 2>/dev/null | head -10")
system("find /mnt/BioDisk/wangmingyuan/SGS10_analysis -name '*_naming*' 2>/dev/null | head -10")

cat("\n=== Load one liver T4 file to see its structure ===\n")
# Will be filled after the find commands show actual paths
