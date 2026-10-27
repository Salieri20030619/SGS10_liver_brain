#!/usr/bin/env Rscript
# =============================================================================
# SGS10 Task A 补图 - Patch v3: 修复 F02 / F05 / F06 / F11 / F15 / F16
# Date: 2026-04-20
# =============================================================================
# 修复点:
#   F02: legend labels 和 values 颠倒 (ns 染红色, sig 染灰色)  -> 改用 named vector
#   F05: Zbtb7a 画了 tumor_vs_control + tumor_vs_blank 40+ 行 -> 只保留 tumor_vs_control 22 个
#   F06: 负向 bar 的 padj 标签被 bar 覆盖                     -> 改 hjust 逻辑
#   F11: y 轴上限太低红虚线贴顶, sig threshold 标签被截       -> 提升 ylim 到 2
#   F15: legend 冗余 (只有 2 柱 颜色已经区分)                  -> 去 legend
#   F16: Mertk 错误地画在 HSC_ECM 里 (Mertk 是 Kupffer 表达)   -> 三分区 + Kupffer 层
# Usage:
#   cd ~/SGS10_analysis/
#   Rscript taskA_figures_v1/patch_v3_fix_six_figs.R
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr); library(tidyr)
  library(stringr); library(patchwork); library(scales); library(ggrepel)
  library(purrr)
})

# ---------- 全局配置 ----------
ROOT_DEG_LIVER <- "~/SGS10_analysis/diff_stage1_results/deg_v4_refined/liver"
ROOT_DEG_BRAIN <- "~/SGS10_analysis/SGS10_deliverables_2026-04-20/03_sc_DEG/brain_clean/brain_Tumor_vs_Control"
ROOT_TF_BRAIN  <- "~/SGS10_analysis/SGS10_deliverables_2026-04-20/05_TF_activity/brain_clean/per_contrast"
OUT_DIR        <- "~/SGS10_analysis/taskA_figures_2026-04-20_v1"
DPI            <- 300
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

LOG_SUMM <- file.path(OUT_DIR, "generation_summary_patch_v3.log")
LOG_MISS <- file.path(OUT_DIR, "missing_data_patch_v3.log")
writeLines(c(paste("# patch_v3 summary -", Sys.time()), ""), LOG_SUMM)
writeLines(c(paste("# patch_v3 missing -", Sys.time()), ""), LOG_MISS)

# ---------- theme ----------
paper_theme <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey92", size = 0.3),
        axis.text  = element_text(color = "black"),
        axis.title = element_text(color = "black", face = "bold"),
        plot.title = element_text(face = "bold", size = 12, hjust = 0),
        plot.subtitle = element_text(size = 9, color = "grey30"),
        legend.position = "right",
        strip.background = element_rect(fill = "grey95", color = NA),
        strip.text = element_text(face = "bold"))

col_up   <- "#c0392b"
col_down <- "#2980b9"
col_ns   <- "#95a5a6"
col_hl   <- "#e67e22"

# ---------- 工具 ----------
log_miss <- function(fig_id, msg) {
  cat(paste0("[MISS] ", fig_id, ": ", msg, "\n"))
  cat(paste0("[", Sys.time(), "] ", fig_id, ": ", msg, "\n"),
      file = LOG_MISS, append = TRUE)
}
log_summ <- function(line) {
  cat(line, "\n")
  cat(line, "\n", file = LOG_SUMM, append = TRUE, sep = "")
}
safe_read_csv <- function(path, fig_id) {
  if (!file.exists(path)) {
    log_miss(fig_id, paste("File not found:", path)); return(NULL)
  }
  tryCatch({
    df <- read_csv(path, show_col_types = FALSE)
    cat(paste0("[OK]   ", fig_id, ": loaded ", basename(path),
               " (", nrow(df), " rows)\n"))
    df
  }, error = function(e) {
    log_miss(fig_id, paste("Read error:", path, "|", e$message)); NULL
  })
}
pick_col <- function(df, candidates) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0) NA_character_ else hit[1]
}
save_png <- function(p, filename, width = 7, height = 5) {
  path <- file.path(OUT_DIR, filename)
  ggsave(path, p, width = width, height = height, dpi = DPI,
         units = "in", bg = "white")
  cat(paste0("[SAVE] ", filename, " (", width, "x", height, " in @ ", DPI, " DPI)\n"))
}

# =============================================================================
# F02 FIXED: legend label 对齐 values
# =============================================================================
make_F02 <- function() {
  fig_id <- "F02"
  f <- file.path(ROOT_DEG_LIVER,
                 "liver_tumor_vs_control_HSC_Activated_ECM.csv")
  hsc <- safe_read_csv(f, fig_id)
  if (is.null(hsc)) return(FALSE)

  genes <- c("Acta2", "Serpine1", "Col1a1", "F2r", "Tagln",
             "Pmepa1", "Col1a2", "Smad3", "Lox", "Loxl1")
  sub <- hsc %>% filter(gene %in% genes) %>%
    mutate(sig_flag = factor(ifelse(p_val_adj < 0.05, "sig", "ns"),
                             levels = c("sig", "ns")))  # 强制 factor level 顺序

  if (nrow(sub) == 0) {
    log_miss(fig_id, "No L1 genes matched")
    return(FALSE)
  }

  # 用 named vector 精确对齐 (key 是 factor level, 不靠字母序)
  fill_map <- c("sig" = col_up, "ns" = col_ns)
  label_map <- c("sig" = "sig (padj<0.05)", "ns" = "ns")

  p <- ggplot(sub, aes(x = reorder(gene, avg_log2FC),
                       y = avg_log2FC, fill = sig_flag)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = sprintf("padj=%.1e", p_val_adj)),
              hjust = -0.05, size = 2.8, color = "grey30") +
    coord_flip(clip = "off") +
    scale_fill_manual(values = fill_map, labels = label_map,
                      breaks = c("sig", "ns")) +  # breaks 固定 legend 顺序
    labs(title = "L1 HSC myofibroblast activation (sc DEG Top 10)",
         subtitle = "HSC_Activated_ECM: Tumor vs Control (Acta2 log2FC 8.36, Serpine1 8.83)",
         x = NULL, y = "avg_log2FC", fill = NULL) +
    paper_theme + theme(legend.position = "top",
                        plot.margin = margin(5, 40, 5, 5))

  save_png(p, "F02_L1_HSC_activation.png", width = 8, height = 5.5)
  TRUE
}

# =============================================================================
# F05 FIXED: 只保留 tumor_vs_control 的 22 个 cell type
# =============================================================================
make_F05 <- function() {
  fig_id <- "F05"
  if (!dir.exists(path.expand(ROOT_TF_BRAIN))) {
    log_miss(fig_id, paste("TF dir not found:", ROOT_TF_BRAIN))
    return(FALSE)
  }

  # 关键修复: pattern 严格限制 tumor_vs_control, 排除 tumor_vs_blank
  tf_files <- list.files(ROOT_TF_BRAIN,
                         pattern = "tumor_vs_control",
                         full.names = TRUE)
  # 保险: 手动过滤掉包含 _vs_blank 的
  tf_files <- tf_files[!grepl("_vs_blank", tf_files, ignore.case = TRUE)]
  tf_files <- tf_files[grepl("\\.csv$", tf_files)]

  if (length(tf_files) == 0) {
    log_miss(fig_id, "No tumor_vs_control TF CSV found")
    return(FALSE)
  }
  cat(paste0("[INFO] F05: scanning ", length(tf_files),
             " tumor_vs_control TF files (blank filtered)\n"))

  zbtb <- purrr::map_dfr(tf_files, function(f) {
    df <- tryCatch(read_csv(f, show_col_types = FALSE), error = function(e) NULL)
    if (is.null(df)) return(NULL)
    tf_col    <- pick_col(df, c("source", "tf", "TF", "tf_name", "regulon"))
    score_col <- pick_col(df, c("score", "activity", "t_value", "est"))
    p_col     <- pick_col(df, c("p_value", "p_adj", "padj", "pval", "p.adjust"))
    if (any(is.na(c(tf_col, score_col, p_col)))) return(NULL)

    # 清洗 cell_type 名称: 从 basename 移除 contrast + 路径噪声
    bname <- basename(f)
    ct <- bname %>%
      str_remove("\\.csv$") %>%
      str_remove("_tf_activity$") %>%
      str_remove("^brain__?") %>%
      str_remove("^liver__?") %>%
      str_remove("tumor_vs_control__?") %>%
      str_remove("__?tumor_vs_control") %>%
      str_replace_all("^_+|_+$", "")

    df %>% filter(.data[[tf_col]] == "Zbtb7a") %>%
      transmute(cell_type = ct,
                score     = .data[[score_col]],
                padj      = .data[[p_col]])
  })

  if (nrow(zbtb) == 0) {
    log_miss(fig_id, "No Zbtb7a rows found")
    return(FALSE)
  }

  # 去重 (防止意外重复)
  zbtb <- zbtb %>% distinct(cell_type, .keep_all = TRUE) %>%
    arrange(score) %>%
    mutate(cell_type = factor(cell_type, levels = cell_type),
           sig_level = case_when(
             padj < 0.05 ~ "strict (padj<0.05)",
             padj < 0.10 ~ "relaxed (padj<0.10)",
             TRUE        ~ "ns"),
           sig_level = factor(sig_level,
                              levels = c("strict (padj<0.05)",
                                         "relaxed (padj<0.10)", "ns")))

  cat(paste0("[INFO] F05: rendering ", nrow(zbtb), " cell types\n"))

  p <- ggplot(zbtb, aes(x = score, y = cell_type, color = sig_level)) +
    geom_vline(xintercept = 0, color = "grey60", size = 0.3) +
    geom_segment(aes(x = 0, xend = score, yend = cell_type), size = 0.5) +
    geom_point(size = 3) +
    geom_text(aes(label = sprintf("padj=%.3f", padj)),
              hjust = 1.15, size = 2.5, color = "grey30") +
    scale_color_manual(values = c("strict (padj<0.05)"   = col_down,
                                  "relaxed (padj<0.10)"  = "#5DADE2",
                                  "ns"                   = col_ns)) +
    labs(title = "Zbtb7a pan-brain signature (new TF signature)",
         subtitle = sprintf("%d cell types consistently down, tumor_vs_control",
                           nrow(zbtb)),
         x = "TF activity score (Tumor vs Control)", y = NULL,
         color = "Significance") +
    paper_theme

  save_png(p, "F05_Zbtb7a_pan_brain_forest.png", width = 8, height = 7)
  TRUE
}

# =============================================================================
# F06 FIXED: 负向 bar 的 padj label 推到 bar 左端外
# =============================================================================
make_F06 <- function() {
  fig_id <- "F06"
  mg <- safe_read_csv(file.path(ROOT_DEG_BRAIN,
                                "Microglia_Homeostatic_deg.csv"), fig_id)
  if (is.null(mg)) return(FALSE)

  genes_up   <- c("Smad3", "Spred2", "Lrp1", "Tgfbr2", "Apoe",
                  "Ski", "Cx3cr1", "Smurf1", "Hipk2", "Tgfb1")
  genes_down <- c("Tgfbr1", "Elmo1", "Dock4", "Mertk", "Fkbp5")

  sub <- mg %>% filter(gene %in% c(genes_up, genes_down)) %>%
    mutate(direction = ifelse(gene %in% genes_up,
                              "Up (Smad axis)", "Down (MerTK axis)"),
           # 关键修复: 正值 hjust=-0.05 (向右推); 负值 hjust=1.05 (向左推), 
           # 但 text 的 y 位置用 0 附近的小偏移, 这样 text 不会被 bar 覆盖
           label_hjust = ifelse(avg_log2FC > 0, -0.05, 1.05),
           label_text = sprintf("padj=%.1e", p_val_adj))

  p <- ggplot(sub, aes(x = reorder(gene, avg_log2FC),
                       y = avg_log2FC, fill = direction)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = 0, color = "grey40", size = 0.3) +
    # 关键: text 的 y 用 avg_log2FC (bar 末端), hjust 分正负推开
    geom_text(aes(label = label_text, hjust = label_hjust),
              size = 2.5, color = "grey30") +
    coord_flip(clip = "off") +
    scale_fill_manual(values = c("Up (Smad axis)"    = col_up,
                                 "Down (MerTK axis)" = col_down)) +
    # 扩充 y 轴范围给 label 让位 (正负两侧各多 0.5)
    scale_y_continuous(expand = expansion(mult = c(0.25, 0.25))) +
    labs(title = "N1 Microglia TGFb-driven state transition",
         subtitle = "Asymmetric regulation: Smad signaling up + MerTK clearance down",
         x = NULL, y = "avg_log2FC") +
    paper_theme + theme(plot.margin = margin(5, 20, 5, 20))

  save_png(p, "F06_N1_Microglia_state.png", width = 9, height = 5.5)
  TRUE
}

# =============================================================================
# F11 FIXED: 扩 y 轴上限给红虚线留空间
# =============================================================================
make_F11 <- function() {
  fig_id <- "F11"
  df <- tribble(
    ~cell_type,              ~feature,   ~padj_deg,  ~padj_tf,
    "Microglia_Homeostatic", "Smad3",    3.7e-13,   0.326,
    "Microglia_Homeostatic", "Apoe",     2.3e-18,   0.800,
    "Astrocyte_Resting",     "Sox9",     9.9e-36,   0.580,
    "Astrocyte_Resting",     "Cdh4",     1.5e-38,   0.604,
    "Pericytes",             "Cldn5",    2.0e-07,   0.711,
    "Pericytes",             "Mfge8",    1.5e-02,   0.842,
    "Endothelial",           "Cldn5",    8.8e-47,   0.974,
    "Endothelial",           "H2-D1",    3.1e-45,   0.974
  ) %>%
    mutate(log10_deg = -log10(padj_deg),
           log10_tf  = -log10(padj_tf))

  p <- ggplot(df, aes(x = log10_deg, y = log10_tf, color = cell_type)) +
    # 红虚线 = TF padj 0.05 阈值, 即 -log10(0.05) = 1.301
    geom_hline(yintercept = -log10(0.05), color = "red",
               linetype = "dashed", size = 0.4) +
    geom_point(size = 4) +
    geom_text_repel(aes(label = feature), size = 3, color = "grey20",
                    min.segment.length = 0, box.padding = 0.3) +
    # 在红虚线正上方放标注 (x 用数据范围中间位置)
    annotate("text", x = 25, y = -log10(0.05) + 0.15,
             label = "TF sig threshold (padj = 0.05)",
             color = "red", size = 3, hjust = 0.5, fontface = "italic") +
    scale_color_brewer(palette = "Set1") +
    coord_cartesian(ylim = c(-0.1, 2)) +  # 关键修复: 扩上限到 2
    labs(title = "TF-decoupled paradigm in non-neuronal brain cells",
         subtitle = "DEG extremely sig (right) but TF activity all NS (bottom)",
         x = "-log10(DEG padj)", y = "-log10(TF padj)",
         color = "Cell type") +
    paper_theme

  save_png(p, "F11_TFdecoupled_dichotomy.png", width = 9, height = 5.5)
  TRUE
}

# =============================================================================
# F15 FIXED: 去掉冗余 legend
# =============================================================================
make_F15 <- function() {
  fig_id <- "F15"
  df <- tribble(
    ~pair,                              ~r,     ~p,    ~conclusion,
    "Liver Tgfb1 ~ MG Smad3",          -0.904, 0.012, "Rejects ligand-driven",
    "Liver Smad3 ~ MG Smad3",          +0.891, 0.014, "Supports parallel autocrine"
  )

  p <- ggplot(df, aes(x = reorder(pair, r), y = r,
                      fill = factor(sign(r)))) +
    geom_col(width = 0.5) +
    geom_hline(yintercept = 0, color = "grey40", size = 0.3) +
    geom_text(aes(label = sprintf("r = %.3f\np = %.3f\n%s", r, p, conclusion)),
              vjust = ifelse(df$r > 0, -0.15, 1.15),
              size = 3.2, color = "grey10", lineheight = 1) +
    coord_cartesian(ylim = c(-1.35, 1.35)) +
    scale_fill_manual(values = c("-1" = col_down, "1" = col_up)) +
    labs(title = "X1 TGFb cross-organ axis - NOT ligand-driven",
         subtitle = "Task B T5 paired correlation (3 tumor mice)",
         x = NULL, y = "Pearson r") +
    paper_theme + theme(legend.position = "none")  # 关键修复

  save_png(p, "F15_X1_TGFb_axis_schematic.png", width = 7, height = 5.5)
  TRUE
}

# =============================================================================
# F16 FIXED: 三层对比 (Liver ligand / Liver receptor / Brain receptor)
# =============================================================================
make_F16 <- function() {
  fig_id <- "F16"

  # 三个数据源
  hsc <- safe_read_csv(file.path(ROOT_DEG_LIVER,
                                 "liver_tumor_vs_control_HSC_Activated_ECM.csv"),
                       fig_id)
  kup <- safe_read_csv(file.path(ROOT_DEG_LIVER,
                                 "liver_tumor_vs_control_Kupffer.csv"),
                       fig_id)
  pe  <- safe_read_csv(file.path(ROOT_DEG_BRAIN, "Pericytes_deg.csv"),
                       fig_id)

  if (is.null(hsc) || is.null(kup) || is.null(pe)) return(FALSE)

  # Liver ligand layer (HSC_ECM 作为 Gas6/C3 source)
  liver_lig <- hsc %>% filter(gene %in% c("Gas6", "C3", "Axl")) %>%
    mutate(panel = "Liver HSC_ECM (ligand source)")

  # Liver receptor layer (Kupffer 作为 Mertk receptor)
  liver_rec <- kup %>% filter(gene %in% c("Mertk", "Axl")) %>%
    mutate(panel = "Liver Kupffer (professional phagocyte)")

  # Brain non-professional phagocyte layer (Pericytes)
  brain_rec <- pe %>% filter(gene %in% c("Mfge8", "Bsg", "Calr", "Pecam1")) %>%
    mutate(panel = "Brain Pericytes (non-professional acquisition)")

  df <- bind_rows(liver_lig, liver_rec, brain_rec) %>%
    mutate(panel = factor(panel, levels = c(
      "Liver HSC_ECM (ligand source)",
      "Liver Kupffer (professional phagocyte)",
      "Brain Pericytes (non-professional acquisition)"
    )))

  if (nrow(df) == 0) {
    log_miss(fig_id, "No efferocytosis genes matched in any layer")
    return(FALSE)
  }

  p <- ggplot(df, aes(x = reorder(gene, avg_log2FC),
                      y = avg_log2FC, fill = panel)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = 0, color = "grey40", size = 0.3) +
    geom_text(aes(label = sprintf("padj=%.1e", p_val_adj)),
              hjust = -0.05, size = 2.5, color = "grey30") +
    coord_flip(clip = "off") +
    facet_wrap(~panel, scales = "free_y", ncol = 1) +
    scale_fill_manual(values = c(
      "Liver HSC_ECM (ligand source)"                    = "#8B4513",
      "Liver Kupffer (professional phagocyte)"           = "#D2691E",
      "Brain Pericytes (non-professional acquisition)"   = "#6A5ACD"
    )) +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.35))) +
    labs(title = "X3 Efferocytosis parallel response (not a causal axis)",
         subtitle = "3-layer parallel: Liver ligand -> Liver pro-phagocyte + Brain non-pro-phagocyte",
         x = NULL, y = "avg_log2FC (Tumor vs Control)") +
    paper_theme + theme(legend.position = "none",
                        plot.margin = margin(5, 40, 5, 5),
                        strip.text = element_text(size = 9))

  save_png(p, "F16_X3_Efferocytosis_parallel.png", width = 9, height = 8)
  TRUE
}

# =============================================================================
# Mini dispatcher - 只跑 6 张 FIX
# =============================================================================
log_summ(paste0("Start patch_v3: ", Sys.time()))
log_summ(paste0("ROOT_DEG_LIVER: ", path.expand(ROOT_DEG_LIVER)))
log_summ(paste0("ROOT_DEG_BRAIN: ", path.expand(ROOT_DEG_BRAIN)))
log_summ(paste0("ROOT_TF_BRAIN : ", path.expand(ROOT_TF_BRAIN)))
log_summ(paste0("OUT_DIR       : ", path.expand(OUT_DIR)))
log_summ("")

fixes <- list(
  F02 = make_F02,
  F05 = make_F05,
  F06 = make_F06,
  F11 = make_F11,
  F15 = make_F15,
  F16 = make_F16
)

results <- tibble(figure_id = character(), status = character(),
                  time_sec = numeric())

for (fid in names(fixes)) {
  log_summ(paste0("--- ", fid, " ---"))
  t0 <- Sys.time()
  ok <- tryCatch(fixes[[fid]](),
                 error = function(e) {
                   log_summ(paste0("[ERR ] ", fid, ": ", e$message))
                   FALSE
                 })
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  st <- if (isTRUE(ok)) "OK" else "FAIL"
  log_summ(paste0("[", st, "] ", fid, " (", round(dt, 1), " s)"))
  results <- bind_rows(results, tibble(figure_id = fid, status = st,
                                       time_sec = round(dt, 1)))
}

log_summ("")
log_summ("========== PATCH v3 SUMMARY ==========")
n_ok <- sum(results$status == "OK")
log_summ(paste0("Fixed: ", nrow(results),
                "  OK: ", n_ok,
                "  FAIL: ", sum(results$status == "FAIL")))
log_summ("")
for (i in seq_len(nrow(results))) {
  log_summ(sprintf("  %-6s  %-5s  %5.1f s",
                   results$figure_id[i], results$status[i],
                   results$time_sec[i]))
}
log_summ(paste0("\nEnd: ", Sys.time()))
