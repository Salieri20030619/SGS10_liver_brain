#!/usr/bin/env Rscript
# ============================================================================
# figS3_dam_cleaning.R
# Supplementary Fig S3: Brain Microglia DAM cleaning audit.
# 5,545 cells originally labeled Microglia_DAM (v2.5) audited with strict
# marker panel (Tmem119, P2ry12, Trem2, Cd9, Itgax). Result: 9 reclassified
# categories, only 57 (1.0%) retained as true Microglia_DAM_strict.
# Data: dam_clean_v2/all_dam_decisions.csv (new_label column)
# ============================================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})
source("/mnt/BioDisk/wangmingyuan/SGS10_analysis/scripts/paper_theme.R")

df <- read.csv("/mnt/BioDisk/wangmingyuan/SGS10_analysis/dam_clean_v2/all_dam_decisions.csv",
               stringsAsFactors = FALSE)
cat(sprintf("[figS3] Total cells: %d\n", nrow(df)))

label_tbl <- df %>%
  count(new_label) %>%
  arrange(desc(n)) %>%
  mutate(pct = n / sum(n) * 100,
         label_disp = paste0(new_label, "  (n=", n,
                             ", ", sprintf("%.1f", pct), "%)"))

cat("\n[figS3] DAM cleaning categories:\n")
print(label_tbl)

# Highlight: low_quality categories grey, lineage reassigned warm, true DAM red
label_tbl$category_type <- case_when(
  grepl("low_quality", label_tbl$new_label) ~ "Low quality (filtered)",
  label_tbl$new_label == "Microglia_DAM_strict" ~ "True DAM (kept)",
  TRUE ~ "Reassigned lineage"
)
label_tbl$category_type <- factor(label_tbl$category_type,
  levels = c("True DAM (kept)", "Reassigned lineage",
             "Low quality (filtered)"))

# Order by n descending for plot
label_tbl <- label_tbl %>% arrange(n) %>%
  mutate(new_label = factor(new_label, levels = new_label))

p <- ggplot(label_tbl, aes(x = n, y = new_label, fill = category_type)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%d (%.1f%%)", n, pct)),
            hjust = -0.08, size = 3.2, fontface = "bold", color = "grey25") +
  scale_fill_manual(
    values = c("True DAM (kept)" = "#C0143B",
               "Reassigned lineage" = "#E67E22",
               "Low quality (filtered)" = "#888888"),
    name = NULL
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.20))) +
  labs(x = sprintf("Cell count (total n=%d)", sum(label_tbl$n)),
       y = NULL) +
  theme_paper() +
  theme(
    axis.text.y = element_text(size = 10, face = "bold"),
    legend.position = "top",
    legend.text = element_text(size = 9)
  )

save_paper_fig(p, "figS3_dam_cleaning",
               canvas = NULL, width = 9, height = 5)
