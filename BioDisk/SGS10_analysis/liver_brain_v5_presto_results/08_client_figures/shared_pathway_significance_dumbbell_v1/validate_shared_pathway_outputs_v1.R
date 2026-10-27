options(cli.width=100L, cli.dynamic=FALSE, cli.ansi=FALSE)
Sys.setenv(R_CLI_DYNAMIC="false", NO_COLOR="1")

project_root <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
source_csv <- file.path(
  project_root, "liver_brain_v5_presto_results", "05_crosstalk",
  "shared_KEGG_pathways.csv"
)
out_dir <- file.path(
  project_root, "liver_brain_v5_presto_results", "08_client_figures",
  "shared_pathway_significance_dumbbell_v1"
)
wide_tsv <- file.path(out_dir, "shared_KEGG_pathway_significance_wide_v1.tsv")

old <- read.csv(source_csv, stringsAsFactors = FALSE, check.names = FALSE)
new <- read.delim(wide_tsv, stringsAsFactors = FALSE, check.names = FALSE)

checks <- data.frame(
  check_id = character(),
  observed = character(),
  expected = character(),
  pass = logical(),
  stringsAsFactors = FALSE
)

add_check <- function(id, observed, expected, pass) {
  checks <<- rbind(
    checks,
    data.frame(
      check_id = id,
      observed = as.character(observed),
      expected = as.character(expected),
      pass = isTRUE(pass),
      stringsAsFactors = FALSE
    )
  )
}

add_check("source_rows", nrow(old), 13L, nrow(old) == 13L)
add_check("new_rows", nrow(new), 13L, nrow(new) == 13L)
add_check(
  "pathway_id_set",
  paste(sort(new$pathway_id), collapse = ","),
  paste(sort(old$ID), collapse = ","),
  setequal(new$pathway_id, old$ID)
)

m <- merge(
  old[, c("ID", "Description", "Liver_minP", "Brain_minP")],
  new[, c("pathway_id", "pathway", "Liver_fdr", "Brain_fdr")],
  by.x = "ID", by.y = "pathway_id", all = TRUE, sort = TRUE
)

max_liver_delta <- max(abs(m$Liver_minP - m$Liver_fdr), na.rm = TRUE)
max_brain_delta <- max(abs(m$Brain_minP - m$Brain_fdr), na.rm = TRUE)
add_check("liver_fdr_max_abs_delta", format(max_liver_delta, scientific = TRUE), "0", max_liver_delta < 1e-15)
add_check("brain_fdr_max_abs_delta", format(max_brain_delta, scientific = TRUE), "0", max_brain_delta < 1e-15)
add_check("description_match", sum(m$Description == m$pathway), 13L, all(m$Description == m$pathway))
add_check("all_fdr_lt_0.05", all(c(new$Liver_fdr, new$Brain_fdr) < 0.05), TRUE,
          all(c(new$Liver_fdr, new$Brain_fdr) < 0.05))

required_files <- c(
  "shared_KEGG_pathway_significance_dumbbell_v1.pdf",
  "shared_KEGG_pathway_significance_grouped_bar_v1.pdf",
  "shared_KEGG_pathway_significance_dumbbell_v1.png",
  "shared_KEGG_pathway_significance_grouped_bar_v1.png",
  "shared_KEGG_pathway_significance_long_v1.tsv",
  "shared_KEGG_pathway_significance_wide_v1.tsv"
)
for (f in required_files) {
  fp <- file.path(out_dir, f)
  add_check(
    paste0("file_", f),
    if (file.exists(fp)) file.info(fp)$size else 0,
    ">0",
    file.exists(fp) && file.info(fp)$size > 0
  )
}

write.table(
  checks,
  file.path(out_dir, "independent_validation_checks_v1.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE, na = "NA"
)

status <- if (all(checks$pass)) "PASS" else "FAIL"
summary_lines <- c(
  paste0("status\t", status),
  paste0("n_checks\t", nrow(checks)),
  paste0("n_pass\t", sum(checks$pass)),
  paste0("n_fail\t", sum(!checks$pass)),
  paste0("liver_fdr_max_abs_delta\t", format(max_liver_delta, scientific = TRUE)),
  paste0("brain_fdr_max_abs_delta\t", format(max_brain_delta, scientific = TRUE))
)
writeLines(summary_lines, file.path(out_dir, "independent_validation_summary_v1.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")

if (status != "PASS") quit(status = 1L)
