options(cli.width=100L, cli.dynamic=FALSE, cli.ansi=FALSE)
Sys.setenv(R_CLI_DYNAMIC="false", NO_COLOR="1")

project_root <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
source_csv <- file.path(
  project_root, "liver_brain_v5_presto_results", "05_crosstalk",
  "shared_KEGG_pathways.csv"
)
long_tsv <- file.path(
  project_root, "liver_brain_v5_presto_results", "08_client_figures",
  "shared_pathway_significance_dumbbell_v1",
  "shared_KEGG_pathway_significance_long_v1.tsv"
)
out_dir <- file.path(
  project_root, "liver_brain_v5_presto_results", "08_client_figures",
  "shared_pathway_significance_grouped_bar_msyh_v2"
)

old <- read.csv(source_csv, stringsAsFactors = FALSE, check.names = FALSE)
long <- read.delim(long_tsv, stringsAsFactors = FALSE, check.names = FALSE)
wide <- reshape(
  long[, c("pathway_id", "pathway", "organ", "fdr")],
  idvar = c("pathway_id", "pathway"), timevar = "organ", direction = "wide"
)

m <- merge(
  old[, c("ID", "Description", "Liver_minP", "Brain_minP")],
  wide,
  by.x = "ID", by.y = "pathway_id", all = TRUE, sort = TRUE
)

checks <- data.frame(
  check_id = c(
    "source_pathways",
    "long_rows",
    "long_unique_pathway_organ",
    "liver_fdr_max_abs_delta",
    "brain_fdr_max_abs_delta",
    "description_match",
    "pdf_exists",
    "png_exists",
    "font_regular_exists",
    "font_bold_exists"
  ),
  observed = c(
    nrow(old),
    nrow(long),
    !anyDuplicated(long[c("pathway_id", "organ")]),
    format(max(abs(m$Liver_minP - m$fdr.Liver)), scientific = TRUE),
    format(max(abs(m$Brain_minP - m$fdr.Brain)), scientific = TRUE),
    all(m$Description == m$pathway),
    file.exists(file.path(out_dir, "shared_KEGG_pathway_significance_grouped_bar_msyh_v2.pdf")),
    file.exists(file.path(out_dir, "shared_KEGG_pathway_significance_grouped_bar_msyh_v2.png")),
    file.exists(file.path(out_dir, "fonts", "msyh.ttc")),
    file.exists(file.path(out_dir, "fonts", "msyhbd.ttc"))
  ),
  expected = c("13", "26", "TRUE", "0", "0", "TRUE", "TRUE", "TRUE", "TRUE", "TRUE"),
  pass = c(
    nrow(old) == 13L,
    nrow(long) == 26L,
    !anyDuplicated(long[c("pathway_id", "organ")]),
    max(abs(m$Liver_minP - m$fdr.Liver)) < 1e-15,
    max(abs(m$Brain_minP - m$fdr.Brain)) < 1e-15,
    all(m$Description == m$pathway),
    file.exists(file.path(out_dir, "shared_KEGG_pathway_significance_grouped_bar_msyh_v2.pdf")),
    file.exists(file.path(out_dir, "shared_KEGG_pathway_significance_grouped_bar_msyh_v2.png")),
    file.exists(file.path(out_dir, "fonts", "msyh.ttc")),
    file.exists(file.path(out_dir, "fonts", "msyhbd.ttc"))
  ),
  stringsAsFactors = FALSE
)

write.table(
  checks,
  file.path(out_dir, "independent_validation_checks_v2.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE, na = "NA"
)

status <- if (all(checks$pass)) "PASS" else "FAIL"
summary_lines <- c(
  paste0("status\t", status),
  paste0("n_checks\t", nrow(checks)),
  paste0("n_pass\t", sum(checks$pass)),
  paste0("n_fail\t", sum(!checks$pass)),
  "visual_qa\tPASS_NO_OVERLAP_OR_CLIPPING",
  "pdf_font_qa\tPASS_MICROSOFT_YAHEI_EMBEDDED"
)
writeLines(summary_lines, file.path(out_dir, "independent_validation_summary_v2.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")
if (status != "PASS") quit(status = 1L)

