#!/usr/bin/env Rscript
# =============================================================================
# SGS10 Task A v1.1 -> Figures generation
# Author: Claude (handoff)
# Date: 2026-04-20
# Purpose: 为 Task A 的 13 条论点 + 3 节 TF 章节 生成 16 张 paper-grade PNG
# =============================================================================
# 使用方法:
#   cd ~/SGS10_analysis/
#   Rscript taskA_figures_v1/generate_all_figures.R
#
# 脚本容错原则:
#   - 每张图独立 tryCatch, 单张失败不影响其他
#   - 文件不存在: 打 [MISS] 警告, 记录到 missing_data.log, SKIP
#   - 列名不匹配: 打 [COLS] 警告, 列出实际列名, SKIP
#   - 最后输出 summary: X/16 成功
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(patchwork)   # 多 panel 拼图
  library(scales)
  library(ggrepel)     # volcano 基因标签
  library(purrr)       # map_dfr / imap_dfr
})

# ---------- 全局配置 ----------
ROOT       <- "~/SGS10_analysis/SGS10_deliverables_2026-04-20"
OUT_DIR    <- "~/SGS10_analysis/taskA_figures_2026-04-20_v1"
CONFIG_CSV <- "~/SGS10_analysis/taskA_figures_v1/config_figures.csv"   # 如放别处请改
DPI        <- 300
LOG_MISS   <- file.path(OUT_DIR, "missing_data.log")
LOG_SUMM   <- file.path(OUT_DIR, "generation_summary.log")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
# 清空日志
writeLines(c(paste("# missing_data.log -", Sys.time()), ""), LOG_MISS)
writeLines(c(paste("# generation_summary.log -", Sys.time()), ""), LOG_SUMM)

# ---------- paper 风格 theme ----------
paper_theme <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey92", size = 0.3),
    axis.text  = element_text(color = "black"),
    axis.title = element_text(color = "black", face = "bold"),
    plot.title = element_text(face = "bold", size = 12, hjust = 0),
    plot.subtitle = element_text(size = 9, color = "grey30"),
    legend.position = "right",
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold")
  )

# 统一配色
col_up   <- "#c0392b"   # up-regulated (红)
col_down <- "#2980b9"   # down-regulated (蓝)
col_ns   <- "#95a5a6"   # NS (灰)
col_hl   <- "#e67e22"   # highlight 橙

# ---------- 工具函数 ----------
log_miss <- function(fig_id, msg) {
  cat(paste0("[MISS] ", fig_id, ": ", msg, "\n"))
  cat(paste0("[", Sys.time(), "] ", fig_id, ": ", msg, "\n"),
      file = LOG_MISS, append = TRUE)
}
log_summ <- function(line) {
  cat(line, "\n")
  cat(line, "\n", file = LOG_SUMM, append = TRUE, sep = "")
}

# 智能 CSV 读取 (尝试容错)
safe_read_csv <- function(path, fig_id) {
  if (!file.exists(path)) {
    log_miss(fig_id, paste("File not found:", path))
    return(NULL)
  }
  tryCatch({
    df <- read_csv(path, show_col_types = FALSE)
    cat(paste0("[OK]   ", fig_id, ": loaded ", basename(path),
               " (", nrow(df), " rows, cols: ",
               paste(head(colnames(df), 6), collapse = ","),
               if (ncol(df) > 6) "..." else "", ")\n"))
    df
  }, error = function(e) {
    log_miss(fig_id, paste("Read error:", path, "|", e$message))
    NULL
  })
}

# 尝试多个候选列名
pick_col <- function(df, candidates) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

# 标准化 DEG 列名 (Seurat FindMarkers 输出一般是 avg_log2FC, p_val_adj, gene 作 rownames 或 column)
normalize_deg <- function(df, fig_id) {
  # gene 列
  g <- pick_col(df, c("gene", "Gene", "symbol", "Symbol", "...1", "X1", "rowname"))
  if (is.na(g)) {
    # 如果没有 gene 列, 检查是否 gene 是 rownames (read_csv 后 rownames 通常在第一列)
    first_col <- colnames(df)[1]
    if (!first_col %in% c("avg_log2FC", "p_val", "p_val_adj", "pct.1", "pct.2")) {
      # 第一列可能就是 gene
      colnames(df)[1] <- "gene"
    } else {
      log_miss(fig_id, paste("No gene column found. Columns:",
                             paste(colnames(df), collapse = ",")))
      return(NULL)
    }
  } else if (g != "gene") {
    df <- dplyr::rename(df, gene = all_of(g))
  }

  # log2FC 列
  fc <- pick_col(df, c("avg_log2FC", "log2FC", "log2FoldChange", "log2fc", "avg_logFC"))
  if (is.na(fc)) {
    log_miss(fig_id, paste("No log2FC column. Cols:",
                           paste(colnames(df), collapse = ",")))
    return(NULL)
  }
  if (fc != "avg_log2FC") df <- dplyr::rename(df, avg_log2FC = all_of(fc))

  # padj 列
  pa <- pick_col(df, c("p_val_adj", "padj", "p.adj", "FDR", "adj_pval"))
  if (is.na(pa)) {
    log_miss(fig_id, paste("No padj column. Cols:",
                           paste(colnames(df), collapse = ",")))
    return(NULL)
  }
  if (pa != "p_val_adj") df <- dplyr::rename(df, p_val_adj = all_of(pa))

  df
}

# 保存 PNG 统一接口
save_png <- function(p, filename, width = 7, height = 5) {
  path <- file.path(OUT_DIR, filename)
  ggsave(path, p, width = width, height = height, dpi = DPI,
         units = "in", bg = "white")
  cat(paste0("[SAVE] ", filename, " (", width, "x", height, " in @ ", DPI, " DPI)\n"))
}

# =============================================================================
# 每张图一个独立函数 make_Fxx()
# =============================================================================

# ---------- F01: X2 Thyroid axis 4-panel ----------
make_F01 <- function() {
  fig_id <- "F01"
  # panel a: Hepatocyte_PV DEG (Ttr, Thra, Thrb)
  pv_f <- file.path(ROOT, "03_sc_DEG", "liver", "Hepatocyte_PV_deg.csv")
  if (!file.exists(pv_f)) {
    # 尝试替代路径
    pv_f <- Sys.glob(file.path(ROOT, "03_sc_DEG", "liver", "*Hepatocyte*PV*deg*.csv"))[1]
  }
  pv <- safe_read_csv(pv_f, fig_id)
  if (is.null(pv)) return(FALSE)
  pv <- normalize_deg(pv, fig_id); if (is.null(pv)) return(FALSE)

  genes_liver <- c("Ttr", "Thra", "Thrb")
  pv_sub <- pv %>% filter(gene %in% genes_liver) %>%
    mutate(compartment = "Liver Hepatocyte_PV")

  # panel b: CPECs + Endothelial (brain)
  cp_f  <- file.path(ROOT, "03_sc_DEG", "brain_clean", "brain_Tumor_vs_Control",
                     "CPECs_deg.csv")
  en_f  <- file.path(ROOT, "03_sc_DEG", "brain_clean", "brain_Tumor_vs_Control",
                     "Endothelial_deg.csv")
  cp <- safe_read_csv(cp_f, fig_id); if (is.null(cp)) return(FALSE)
  en <- safe_read_csv(en_f, fig_id); if (is.null(en)) return(FALSE)
  cp <- normalize_deg(cp, fig_id); en <- normalize_deg(en, fig_id)

  brain_rows <- bind_rows(
    cp %>% filter(gene %in% c("Ttr", "Atp1a1", "Folr1")) %>%
      mutate(compartment = "Brain CPECs"),
    en %>% filter(gene == "Thra") %>% mutate(compartment = "Brain Endothelial")
  )

  # 合并 panel a+b 画双分区 bar
  df_ab <- bind_rows(pv_sub, brain_rows) %>%
    mutate(
      compartment = factor(compartment,
                           levels = c("Liver Hepatocyte_PV",
                                      "Brain CPECs",
                                      "Brain Endothelial")),
      gene_label = paste0(gene, "\npadj=", signif(p_val_adj, 2))
    )

  p_ab <- ggplot(df_ab, aes(x = reorder(gene, avg_log2FC),
                            y = avg_log2FC,
                            fill = compartment)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = 0, color = "grey40", size = 0.3) +
    geom_text(aes(label = paste0("padj=", signif(p_val_adj, 2))),
              hjust = -0.1, size = 2.8, color = "grey30") +
    coord_flip() +
    facet_wrap(~compartment, scales = "free_y", ncol = 1) +
    scale_fill_manual(values = c("#8B4513", "#6A5ACD", "#2E8B57")) +
    labs(title = "X2 Thyroid hormone cross-organ axis",
         subtitle = "Liver anchor (Ttr, Thra, Thrb) -> Brain gateway (CPECs Ttr padj 9.6e-10) -> BBB sensor (Endo Thra padj 2.5e-06)",
         x = NULL, y = "avg_log2FC (Tumor vs Control)") +
    paper_theme + theme(legend.position = "none",
                        strip.text = element_text(face = "bold", size = 10))

  # panel c: KEGG thyroid pathway enrichment 跨 cell type (可选, 若找得到)
  kegg_dir <- file.path(ROOT, "04_enrichment", "brain_clean")
  if (dir.exists(kegg_dir)) {
    # 搜索包含 "thyroid" 的 KEGG 结果
    kegg_files <- list.files(kegg_dir, pattern = "KEGG.*csv$",
                             recursive = TRUE, full.names = TRUE)
    if (length(kegg_files) > 0) {
      kegg_combined <- purrr::map_dfr(kegg_files, function(f) {
        tryCatch({
          df <- read_csv(f, show_col_types = FALSE)
          df$source_file <- basename(f)
          df
        }, error = function(e) NULL)
      })
      # 筛选 thyroid-related pathways
      if (nrow(kegg_combined) > 0) {
        desc_col <- pick_col(kegg_combined,
                             c("Description", "description", "term", "Term", "pathway"))
        padj_col <- pick_col(kegg_combined,
                             c("p.adjust", "padj", "p_val_adj", "FDR"))
        if (!is.na(desc_col) && !is.na(padj_col)) {
          thy <- kegg_combined %>%
            filter(str_detect(tolower(.data[[desc_col]]), "thyroid")) %>%
            mutate(cell_type = str_extract(source_file, "^[^_]+")) %>%
            filter(.data[[padj_col]] < 0.05)
          if (nrow(thy) > 0) {
            p_c <- ggplot(thy, aes(x = reorder(cell_type, -log10(.data[[padj_col]])),
                                   y = -log10(.data[[padj_col]]))) +
              geom_col(fill = "#2E8B57", width = 0.7) +
              coord_flip() +
              labs(title = "KEGG Thyroid pathway enrichment across cell types",
                   x = NULL, y = "-log10(p.adjust)") +
              paper_theme
          } else { p_c <- NULL }
        } else { p_c <- NULL }
      } else { p_c <- NULL }
    } else { p_c <- NULL }
  } else { p_c <- NULL }

  # 合成
  if (!is.null(p_c)) {
    final <- p_ab / p_c + plot_layout(heights = c(2, 1)) +
      plot_annotation(tag_levels = "a")
    save_png(final, "F01_X2_thyroid_axis_4panel.png", width = 9, height = 10)
  } else {
    save_png(p_ab, "F01_X2_thyroid_axis_4panel.png", width = 9, height = 7)
  }
  TRUE
}

# ---------- F02: L1 HSC activation ----------
make_F02 <- function() {
  fig_id <- "F02"
  hsc_f <- file.path(ROOT, "03_sc_DEG", "liver", "HSC_Activated_ECM_deg.csv")
  if (!file.exists(hsc_f)) {
    cand <- Sys.glob(file.path(ROOT, "03_sc_DEG", "liver",
                               "*HSC*Activated*ECM*deg*.csv"))
    if (length(cand) > 0) hsc_f <- cand[1]
  }
  hsc <- safe_read_csv(hsc_f, fig_id)
  if (is.null(hsc)) return(FALSE)
  hsc <- normalize_deg(hsc, fig_id); if (is.null(hsc)) return(FALSE)

  genes <- c("Acta2", "Serpine1", "Col1a1", "F2r", "Tagln",
             "Pmepa1", "Col1a2", "Smad3", "Lox", "Loxl1")
  sub <- hsc %>% filter(gene %in% genes) %>%
    mutate(sig = ifelse(p_val_adj < 0.05, "sig", "ns"))

  p_bar <- ggplot(sub, aes(x = reorder(gene, avg_log2FC),
                           y = avg_log2FC, fill = sig)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = paste0("padj=", signif(p_val_adj, 2))),
              hjust = -0.05, size = 2.7, color = "grey30") +
    coord_flip(clip = "off") +
    scale_fill_manual(values = c(sig = col_up, ns = col_ns)) +
    labs(title = "L1 HSC myofibroblast activation (sc DEG)",
         subtitle = "HSC_Activated_ECM: Tumor vs Control",
         x = NULL, y = "avg_log2FC") +
    paper_theme + theme(legend.position = "none",
                        plot.margin = margin(5, 30, 5, 5))

  save_png(p_bar, "F02_L1_HSC_activation.png", width = 7.5, height = 5)
  TRUE
}

# ---------- F03: L4 HSC-Immune hub ----------
make_F03 <- function() {
  fig_id <- "F03"
  hsc_f <- file.path(ROOT, "03_sc_DEG", "liver", "HSC_Activated_ECM_deg.csv")
  if (!file.exists(hsc_f)) {
    cand <- Sys.glob(file.path(ROOT, "03_sc_DEG", "liver",
                               "*HSC*Activated*ECM*deg*.csv"))
    if (length(cand) > 0) hsc_f <- cand[1]
  }
  hsc <- safe_read_csv(hsc_f, fig_id)
  if (is.null(hsc)) return(FALSE)
  hsc <- normalize_deg(hsc, fig_id); if (is.null(hsc)) return(FALSE)

  genes <- c("Cd74", "Cd68", "Csf1r", "Csf1", "Mertk", "Gas6")
  sub <- hsc %>% filter(gene %in% genes)

  if (nrow(sub) == 0) {
    log_miss(fig_id, "None of myeloid markers found in HSC_ECM DEG - likely wrong cell type; these are myeloid markers not HSC DEG targets. Consider switching to Macrophage/Kupffer DEG file.")
    # 尝试 Macrophage_M1
    m1_f <- Sys.glob(file.path(ROOT, "03_sc_DEG", "liver",
                               "*Macrophage*M1*deg*.csv"))[1]
    if (!is.na(m1_f) && file.exists(m1_f)) {
      cat("[INFO] F03: using Macrophage_M1 DEG instead\n")
      hsc <- safe_read_csv(m1_f, fig_id)
      if (is.null(hsc)) return(FALSE)
      hsc <- normalize_deg(hsc, fig_id); if (is.null(hsc)) return(FALSE)
      sub <- hsc %>% filter(gene %in% genes)
    }
  }

  if (nrow(sub) == 0) return(FALSE)

  p <- ggplot(sub, aes(x = reorder(gene, avg_log2FC), y = avg_log2FC)) +
    geom_col(fill = "#8e44ad", width = 0.7) +
    geom_text(aes(label = paste0("padj=", signif(p_val_adj, 2))),
              hjust = -0.05, size = 2.7, color = "grey30") +
    coord_flip(clip = "off") +
    labs(title = "L4 HSC-Immune hub (myeloid markers)",
         subtitle = "See also Task B Figure 2b for spatial feature plot validation",
         x = NULL, y = "avg_log2FC") +
    paper_theme + theme(plot.margin = margin(5, 30, 5, 5))

  save_png(p, "F03_L4_HSC_immune_hub.png", width = 7, height = 4)
  TRUE
}

# ---------- F04: N4 Endothelial BBB ----------
make_F04 <- function() {
  fig_id <- "F04"
  en <- safe_read_csv(file.path(ROOT, "03_sc_DEG", "brain_clean",
                                "brain_Tumor_vs_Control", "Endothelial_deg.csv"),
                      fig_id)
  if (is.null(en)) return(FALSE)
  en <- normalize_deg(en, fig_id); if (is.null(en)) return(FALSE)

  genes <- c("Cldn5", "H2-D1", "H2-K1", "Cxcl12", "Thra",
             "Ctnnb1", "Cdh5", "Flt1", "Pecam1")
  sub <- en %>% filter(gene %in% genes) %>%
    mutate(is_cldn5 = gene == "Cldn5",
           label_text = sprintf("padj=%.1e", p_val_adj))

  p <- ggplot(sub, aes(x = reorder(gene, avg_log2FC), y = avg_log2FC,
                       fill = is_cldn5)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = label_text),
              hjust = -0.05, size = 2.7,
              color = ifelse(sub$is_cldn5, col_hl, "grey30"),
              fontface = ifelse(sub$is_cldn5, "bold", "plain")) +
    coord_flip(clip = "off") +
    scale_fill_manual(values = c("TRUE" = col_hl, "FALSE" = col_up)) +
    labs(title = "N4 Endothelial BBB dual-role",
         subtitle = "Cldn5 padj 8.8e-47 - global strongest significance in project",
         x = NULL, y = "avg_log2FC") +
    paper_theme + theme(legend.position = "none",
                        plot.margin = margin(5, 40, 5, 5))

  save_png(p, "F04_N4_Endothelial_BBB.png", width = 7.5, height = 5)
  TRUE
}

# ---------- F05: Zbtb7a pan-brain forest plot ----------
make_F05 <- function() {
  fig_id <- "F05"
  tf_dir <- file.path(ROOT, "05_TF_activity", "brain_clean", "per_contrast")
  if (!dir.exists(tf_dir)) {
    log_miss(fig_id, paste("TF dir not found:", tf_dir))
    return(FALSE)
  }
  tf_files <- list.files(tf_dir, pattern = "Tumor_vs_Control.*\\.csv$",
                         full.names = TRUE)
  if (length(tf_files) == 0) {
    tf_files <- list.files(tf_dir, pattern = "\\.csv$", full.names = TRUE)
  }
  if (length(tf_files) == 0) {
    log_miss(fig_id, paste("No TF CSV in", tf_dir))
    return(FALSE)
  }

  zbtb <- purrr::map_dfr(tf_files, function(f) {
    df <- tryCatch(read_csv(f, show_col_types = FALSE), error = function(e) NULL)
    if (is.null(df)) return(NULL)
    tf_col    <- pick_col(df, c("source", "tf", "TF", "tf_name", "regulon"))
    score_col <- pick_col(df, c("score", "activity", "t_value", "est"))
    p_col     <- pick_col(df, c("p_value", "p_adj", "padj", "pval", "p.adjust"))
    if (any(is.na(c(tf_col, score_col, p_col)))) return(NULL)

    df %>% filter(.data[[tf_col]] == "Zbtb7a") %>%
      transmute(cell_type = str_remove(basename(f),
                                       "_Tumor_vs_Control.*|\\.csv$|_tf.*"),
                score = .data[[score_col]],
                padj  = .data[[p_col]])
  })

  if (nrow(zbtb) == 0) {
    log_miss(fig_id, "No Zbtb7a rows found across TF files - check TF column name in your files")
    return(FALSE)
  }

  zbtb <- zbtb %>% arrange(score) %>%
    mutate(cell_type = factor(cell_type, levels = cell_type),
           sig_level = case_when(
             padj < 0.05 ~ "strict (padj<0.05)",
             padj < 0.10 ~ "relaxed (padj<0.10)",
             TRUE        ~ "ns"),
           sig_level = factor(sig_level,
                              levels = c("strict (padj<0.05)",
                                         "relaxed (padj<0.10)", "ns")))

  p <- ggplot(zbtb, aes(x = score, y = cell_type, color = sig_level)) +
    geom_vline(xintercept = 0, color = "grey60", size = 0.3) +
    geom_segment(aes(x = 0, xend = score, yend = cell_type), size = 0.5) +
    geom_point(size = 3) +
    geom_text(aes(label = sprintf("padj=%.3f", padj)),
              hjust = 1.1, size = 2.5, color = "grey30") +
    scale_color_manual(values = c("strict (padj<0.05)"   = col_down,
                                  "relaxed (padj<0.10)" = "#5DADE2",
                                  "ns"                  = col_ns)) +
    labs(title = "Zbtb7a pan-brain signature (new signature)",
         subtitle = "22 cell types consistently down-regulated, 12 strict/relaxed sig",
         x = "TF activity score (Tumor vs Control)", y = NULL,
         color = "Significance") +
    paper_theme

  save_png(p, "F05_Zbtb7a_pan_brain_forest.png", width = 8, height = 7)
  TRUE
}

# ---------- F06: N1 Microglia Top10 (up + down 双向) ----------
make_F06 <- function() {
  fig_id <- "F06"
  mg <- safe_read_csv(file.path(ROOT, "03_sc_DEG", "brain_clean",
                                "brain_Tumor_vs_Control",
                                "Microglia_Homeostatic_deg.csv"), fig_id)
  if (is.null(mg)) return(FALSE)
  mg <- normalize_deg(mg, fig_id); if (is.null(mg)) return(FALSE)

  genes_up   <- c("Smad3", "Spred2", "Lrp1", "Tgfbr2", "Apoe",
                  "Ski", "Cx3cr1", "Smurf1", "Hipk2", "Tgfb1")
  genes_down <- c("Tgfbr1", "Elmo1", "Dock4", "Mertk", "Fkbp5")

  sub <- mg %>% filter(gene %in% c(genes_up, genes_down)) %>%
    mutate(direction = ifelse(gene %in% genes_up,
                              "Up (Smad axis)", "Down (MerTK axis)"))

  p <- ggplot(sub, aes(x = reorder(gene, avg_log2FC),
                       y = avg_log2FC, fill = direction)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = 0, color = "grey40", size = 0.3) +
    geom_text(aes(label = sprintf("padj=%.1e", p_val_adj),
                  hjust = ifelse(avg_log2FC > 0, -0.05, 1.05)),
              size = 2.5, color = "grey30") +
    coord_flip(clip = "off") +
    scale_fill_manual(values = c("Up (Smad axis)" = col_up,
                                 "Down (MerTK axis)" = col_down)) +
    labs(title = "N1 Microglia TGFb-driven state transition",
         subtitle = "Asymmetric regulation: Smad signaling up + MerTK clearance down",
         x = NULL, y = "avg_log2FC") +
    paper_theme + theme(plot.margin = margin(5, 45, 5, 45))

  save_png(p, "F06_N1_Microglia_state.png", width = 8, height = 5.5)
  TRUE
}

# ---------- F07: L6 Silent Inflammation ----------
make_F07 <- function() {
  fig_id <- "F07"
  # 先用报告内已知数字创建 df (若有 spatial CSV 可改为读取真实值)
  inline <- tibble(
    cytokine = rep(c("Il1b", "Tnf", "Il6", "Ifng"), each = 2),
    metric   = rep(c("Max pct (10 samples)", "Final report pct"), 4),
    value    = c(0.72, 0.57,  0.68, 0.26,  0.08, 0.02,  0.06, 0.06)
  )

  p <- ggplot(inline, aes(x = cytokine, y = value, fill = metric)) +
    geom_col(position = position_dodge(0.8), width = 0.7) +
    geom_hline(yintercept = 5, color = col_down, linetype = "dashed") +
    annotate("text", x = 0.8, y = 5.5,
             label = "Reference: sepsis model TNF 5-20%",
             hjust = 0, size = 3, color = col_down) +
    scale_fill_manual(values = c("Max pct (10 samples)" = "#c0392b",
                                 "Final report pct"     = "#f1948a")) +
    labs(title = "L6 Silent Inflammation (negative finding across 10 samples)",
         subtitle = "All 4 cytokines far below true inflammation response",
         x = NULL, y = "Expression pct (%)", fill = NULL) +
    paper_theme +
    coord_cartesian(ylim = c(0, 8))

  save_png(p, "F07_L6_Silent_Inflammation.png", width = 7.5, height = 5)
  TRUE
}

# ---------- F08: L2 TGFb autocrine (inline 数据) ----------
make_F08 <- function() {
  fig_id <- "F08"
  # panel c: Pearson r 近零 (3 样本)
  pearson <- tibble(
    sample = c("Liver_3", "Liver_8", "gan56"),
    r      = c(-0.054, -0.060, -0.070)
  )
  p_c <- ggplot(pearson, aes(x = sample, y = r)) +
    geom_hline(yintercept = 0, color = "grey50", linetype = "dashed") +
    geom_segment(aes(xend = sample, y = 0, yend = r), color = "grey60") +
    geom_point(size = 5, color = col_down) +
    geom_text(aes(label = sprintf("r=%.3f", r)),
              vjust = 1.8, size = 3.5, color = "grey20") +
    coord_cartesian(ylim = c(-0.15, 0.15)) +
    labs(title = "PDX TGFb ligand vs host TGFb response",
         subtitle = "All 3 samples near-zero r - rejects PDX paracrine",
         x = NULL, y = "Pearson r") +
    paper_theme

  # panel b: Smad3 分 zone (3 样本 × 3 zone)
  zone_df <- tibble(
    sample = rep(c("Liver_3", "Liver_8", "gan56"), each = 3),
    zone   = rep(c("tumor", "margin", "host"), 3),
    pct    = c(0.19, 1.60, 11.89,  0.31, 6.65, 10.03,  0.65, 10.01, 30.85)
  ) %>% mutate(zone = factor(zone, levels = c("tumor", "margin", "host")))
  p_b <- ggplot(zone_df, aes(x = zone, y = pct, fill = sample)) +
    geom_col(position = position_dodge(0.8), width = 0.7) +
    geom_text(aes(label = sprintf("%.1f", pct)),
              position = position_dodge(0.8), vjust = -0.3, size = 2.8) +
    labs(title = "Smad3 spatial: host > margin > tumor",
         x = NULL, y = "Smad3 pct (%)") +
    paper_theme + theme(legend.position = "top")

  # panel a: 7 TGFb genes spatial pct 5 样本 (inline)
  tgfb_df <- tribble(
    ~gene,    ~Con,  ~Bla,  ~Liver_3, ~Liver_8, ~gan56,
    "Smad3",   4.59,  4.77, 11.89,    10.03,    30.85,
    "Smad2",   3.95,  4.93, 11.51,    9.82,     10.78,
    "Tgfbr2",  3.58,  3.74, 20.75,    12.10,    5.35,
    "Tgfbr1",  5.07,  4.73, 14.00,    9.88,     10.09,
    "Tgfb1",   1.09,  1.36,  6.79,    3.35,     2.61,
    "Smad7",   3.05,  2.05,  7.06,    4.97,     14.81,
    "Tgfbr3",  5.69,  5.58, 13.39,    8.03,     9.44
  ) %>% pivot_longer(-gene, names_to = "sample", values_to = "pct") %>%
    mutate(sample = factor(sample, levels = c("Con", "Bla",
                                              "Liver_3", "Liver_8", "gan56")))
  p_a <- ggplot(tgfb_df, aes(x = sample, y = gene, fill = pct)) +
    geom_tile(color = "white", size = 0.3) +
    geom_text(aes(label = sprintf("%.1f", pct)), size = 2.6, color = "grey10") +
    scale_fill_gradient(low = "#fdf5e6", high = "#c0392b") +
    labs(title = "L2 Host TGFb pathway activation",
         subtitle = "7 genes x 5 samples spatial expression pct",
         x = NULL, y = NULL, fill = "pct (%)") +
    paper_theme + theme(axis.text.x = element_text(angle = 15, hjust = 1))

  final <- (p_a / (p_b | p_c)) +
    plot_layout(heights = c(1.2, 1)) +
    plot_annotation(tag_levels = "a")
  save_png(final, "F08_L2_TGFb_autocrine.png", width = 10, height = 9)
  TRUE
}

# ---------- F09: N2 Astrocyte Sox9 ----------
make_F09 <- function() {
  fig_id <- "F09"
  as <- safe_read_csv(file.path(ROOT, "03_sc_DEG", "brain_clean",
                                "brain_Tumor_vs_Control",
                                "Astrocyte_Resting_deg.csv"), fig_id)
  if (is.null(as)) return(FALSE)
  as <- normalize_deg(as, fig_id); if (is.null(as)) return(FALSE)

  # 按功能分 4 组
  gene_groups <- list(
    "K+ buffering"      = c("Kcnj10", "Atp1a2"),
    "Glutamate handling" = c("Grm3", "Slc1a2"),
    "Adhesion / structure" = c("Cdh4", "Sox9", "Mtss2", "Sorcs2"),
    "Metabolism / other" = c("Ubc", "Xylt1", "Fasn", "Gfap", "Sparcl1")
  )
  df <- purrr::imap_dfr(gene_groups, function(genes, grp) {
    as %>% filter(gene %in% genes) %>% mutate(func_group = grp)
  })

  p <- ggplot(df, aes(x = reorder(gene, avg_log2FC), y = avg_log2FC,
                      fill = func_group)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = sprintf("padj=%.1e", p_val_adj)),
              hjust = -0.05, size = 2.5, color = "grey30") +
    coord_flip(clip = "off") +
    facet_wrap(~func_group, scales = "free_y", ncol = 2) +
    scale_fill_brewer(palette = "Set2") +
    labs(title = "N2 Astrocyte Sox9-driven reprogramming",
         subtitle = "Functional upgrade rather than classic GFAP-high reactive astrogliosis",
         x = NULL, y = "avg_log2FC") +
    paper_theme + theme(legend.position = "none",
                        plot.margin = margin(5, 25, 5, 5))

  save_png(p, "F09_N2_Astrocyte_Sox9.png", width = 9, height = 6.5)
  TRUE
}

# ---------- F10: N3 Pericyte triple activation ----------
make_F10 <- function() {
  fig_id <- "F10"
  pe <- safe_read_csv(file.path(ROOT, "03_sc_DEG", "brain_clean",
                                "brain_Tumor_vs_Control",
                                "Pericytes_deg.csv"), fig_id)
  if (is.null(pe)) return(FALSE)
  pe <- normalize_deg(pe, fig_id); if (is.null(pe)) return(FALSE)

  gene_groups <- list(
    "BBB structural"     = c("Cldn5", "Klf4"),
    "Efferocytosis"      = c("Mfge8", "Bsg", "Calr", "Pecam1"),
    "IFN-gamma response" = c("Gbp4")
  )
  df <- purrr::imap_dfr(gene_groups, function(genes, grp) {
    pe %>% filter(gene %in% genes) %>% mutate(activation = grp)
  })

  p <- ggplot(df, aes(x = reorder(gene, avg_log2FC),
                      y = avg_log2FC, fill = activation)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = sprintf("padj=%.1e", p_val_adj)),
              hjust = -0.05, size = 2.6, color = "grey30") +
    coord_flip(clip = "off") +
    facet_wrap(~activation, scales = "free_y", ncol = 1) +
    scale_fill_brewer(palette = "Dark2") +
    labs(title = "N3 Pericyte triple vascular-immune activation",
         subtitle = "+ Zbtb7a TF -4.28 padj 0.003 (1/173 strict sig) - see F05",
         x = NULL, y = "avg_log2FC") +
    paper_theme + theme(legend.position = "none",
                        plot.margin = margin(5, 35, 5, 5))

  save_png(p, "F10_N3_Pericyte_triple.png", width = 7.5, height = 6.5)
  TRUE
}

# ---------- F11: TF-decoupled dichotomy scatter ----------
make_F11 <- function() {
  fig_id <- "F11"
  # 用报告内已知 top DEG padj + 典型 TF padj (inline)
  df <- tribble(
    ~cell_type,              ~feature,   ~padj_deg,  ~padj_tf,
    "Microglia_Homeostatic", "Smad3",    3.7e-13,   0.326,   # Stat1
    "Microglia_Homeostatic", "Apoe",     2.3e-18,   0.800,   # Nfkb1
    "Astrocyte_Resting",     "Sox9",     9.9e-36,   0.580,   # Stat3
    "Astrocyte_Resting",     "Cdh4",     1.5e-38,   0.604,   # Stat1
    "Pericytes",             "Cldn5",    2.0e-07,   0.711,   # Stat1
    "Pericytes",             "Mfge8",    1.5e-02,   0.842,   # Nfkb1
    "Endothelial",           "Cldn5",    8.8e-47,   0.974,   # Stat1
    "Endothelial",           "H2-D1",    3.1e-45,   0.974    # Rela
  ) %>%
    mutate(log10_deg = -log10(padj_deg),
           log10_tf  = -log10(padj_tf))

  p <- ggplot(df, aes(x = log10_deg, y = log10_tf, color = cell_type)) +
    geom_hline(yintercept = -log10(0.05), color = "red",
               linetype = "dashed", size = 0.3) +
    annotate("text", x = 30, y = -log10(0.05) + 0.1,
             label = "TF sig threshold (padj=0.05)",
             color = "red", size = 2.8, hjust = 0) +
    geom_point(size = 4) +
    geom_text_repel(aes(label = feature), size = 3, color = "grey20") +
    scale_color_brewer(palette = "Set1") +
    labs(title = "TF-decoupled paradigm in non-neuronal brain cells",
         subtitle = "DEG extremely sig (right) but TF activity all NS (bottom)",
         x = "-log10(DEG padj)", y = "-log10(TF padj)",
         color = "Cell type") +
    paper_theme

  save_png(p, "F11_TFdecoupled_dichotomy.png", width = 8.5, height = 5.5)
  TRUE
}

# ---------- F12: Excitatory_L45 TF-driven forest ----------
make_F12 <- function() {
  fig_id <- "F12"
  tf_f <- file.path(ROOT, "05_TF_activity", "brain_clean", "per_contrast",
                    "Excitatory_L45_Tumor_vs_Control_tf_activity.csv")
  if (!file.exists(tf_f)) {
    cand <- Sys.glob(file.path(ROOT, "05_TF_activity", "brain_clean",
                               "per_contrast", "*Excitatory*L45*.csv"))
    if (length(cand) > 0) tf_f <- cand[1]
  }
  tf <- safe_read_csv(tf_f, fig_id)
  if (is.null(tf)) {
    # fallback: 用报告内 inline 数据
    log_miss(fig_id, "Using inline data from report")
    tf <- tribble(
      ~source, ~score,  ~padj,
      "Tead1",   -4.72, 3.6e-04,
      "Tead4",   -3.91, 7.0e-03,
      "Pou2f2",  -3.78, 7.9e-03,
      "Pou4f2",  -3.52, 1.6e-02,
      "Mef2a",   -3.39, 2.1e-02,
      "Meis2",   -3.26, 2.7e-02,
      "Zbtb7a",  -3.17, 2.7e-02,
      "Maz",     +3.14, 2.7e-02,
      "Ikzf1",   -3.12, 2.7e-02,
      "Irf4",    -3.01, 3.2e-02
    )
  } else {
    tf_col    <- pick_col(tf, c("source", "tf", "TF", "tf_name", "regulon"))
    score_col <- pick_col(tf, c("score", "activity", "t_value"))
    p_col     <- pick_col(tf, c("p_value", "p_adj", "padj", "pval", "p.adjust"))
    if (any(is.na(c(tf_col, score_col, p_col)))) {
      log_miss(fig_id, "TF columns not matched - check column names")
      return(FALSE)
    }
    tf <- tf %>%
      transmute(source = .data[[tf_col]],
                score  = .data[[score_col]],
                padj   = .data[[p_col]]) %>%
      filter(padj < 0.05) %>%
      arrange(padj) %>% head(19)
  }

  tf <- tf %>% arrange(score) %>%
    mutate(source = factor(source, levels = source),
           dir = ifelse(score > 0, "Up", "Down"))

  p <- ggplot(tf, aes(x = score, y = source, color = dir)) +
    geom_vline(xintercept = 0, color = "grey60", size = 0.3) +
    geom_segment(aes(x = 0, xend = score, yend = source), size = 0.5) +
    geom_point(size = 3) +
    geom_text(aes(label = sprintf("padj=%.0e", padj)),
              hjust = 1.1, size = 2.5, color = "grey30") +
    scale_color_manual(values = c("Up" = col_up, "Down" = col_down)) +
    labs(title = "Excitatory_L45 classical TF-driven response",
         subtitle = "Neuronal identity TFs (Tead/Pou/Mef2) systematically down",
         x = "TF activity score", y = NULL, color = "Direction") +
    paper_theme

  save_png(p, "F12_ExcL45_TFdriven.png", width = 7.5, height = 6)
  TRUE
}

# ---------- F13: L3 Margin width vs HSC response ----------
make_F13 <- function() {
  fig_id <- "F13"
  df <- tibble(
    sample   = c("Liver_3", "Liver_8", "gan56"),
    margin   = c(12.3, 1.5, 2.4),
    acta2    = c(1.37, 0.27, 0.20)
  )

  p <- ggplot(df, aes(x = margin, y = acta2, label = sample)) +
    geom_point(size = 5, color = col_up) +
    geom_text_repel(size = 3.5, color = "grey10", seed = 42) +
    geom_smooth(method = "lm", se = FALSE, color = "grey40",
                linetype = "dashed", size = 0.5) +
    labs(title = "L3 Margin width positively correlates with HSC activation",
         subtitle = "3 tumor samples; rank-consistent across all HSC markers",
         x = "Margin pct (%)", y = "HSC Acta2 pct (%)") +
    paper_theme

  save_png(p, "F13_L3_margin_response.png", width = 6.5, height = 5)
  TRUE
}

# ---------- F14: L5 F2r sensitization ----------
make_F14 <- function() {
  fig_id <- "F14"
  df <- tribble(
    ~gene, ~cell_type,                ~condition, ~pct,
    "F2r",   "HSC_Activated_ECM",      "Tumor",    34.8,
    "F2r",   "HSC_Activated_ECM",      "Control",  1.2,
    "F2r",   "HSC (general)",          "Tumor",    3.2,
    "F2r",   "HSC (general)",          "Control",  2.5,
    "F2",    "Hepatocyte_PP",          "Tumor",    65.3,
    "F2",    "Hepatocyte_PP",          "Control",  73.2,
    "Plg",   "Hepatocyte_PP",          "Tumor",    86.8,
    "Plg",   "Hepatocyte_PP",          "Control",  91.5
  ) %>% mutate(label = paste0(gene, " in ", cell_type),
               condition = factor(condition, levels = c("Control", "Tumor")))

  p <- ggplot(df, aes(x = label, y = pct, fill = condition)) +
    geom_col(position = position_dodge(0.8), width = 0.7) +
    geom_text(aes(label = sprintf("%.1f", pct)),
              position = position_dodge(0.8), vjust = -0.3, size = 2.8) +
    coord_flip() +
    scale_fill_manual(values = c("Control" = col_ns, "Tumor" = col_up)) +
    labs(title = "L5 F2r (PAR1) 20x receptor sensitization",
         subtitle = "F2r in HSC_ECM jumps 1.2% -> 34.8% while ligand F2 in Hep stays stable",
         x = NULL, y = "pct (%)", fill = NULL) +
    paper_theme + theme(legend.position = "top")

  save_png(p, "F14_L5_F2r_sensitization.png", width = 8, height = 5)
  TRUE
}

# ---------- F15: X1 TGFb axis schematic ----------
make_F15 <- function() {
  fig_id <- "F15"
  # 用已知数字 (来自 Task B T5 summary)
  df <- tribble(
    ~pair,                              ~r,     ~p,    ~conclusion,
    "Liver Tgfb1 ~ MG Smad3",          -0.904, 0.012, "Rejects ligand-driven",
    "Liver Smad3 ~ MG Smad3",          +0.891, 0.014, "Supports parallel autocrine"
  )

  p <- ggplot(df, aes(x = reorder(pair, r), y = r,
                      fill = factor(sign(r), levels = c(-1, 1),
                                    labels = c("Negative", "Positive")))) +
    geom_col(width = 0.5) +
    geom_hline(yintercept = 0, color = "grey40", size = 0.3) +
    geom_text(aes(label = sprintf("r=%.3f\np=%.3f\n%s", r, p, conclusion)),
              vjust = ifelse(df$r > 0, -0.1, 1.1), size = 3, color = "grey10") +
    coord_cartesian(ylim = c(-1.3, 1.3)) +
    scale_fill_manual(values = c("Negative" = col_down, "Positive" = col_up)) +
    labs(title = "X1 TGFb cross-organ axis - NOT ligand-driven",
         subtitle = "Task B T5 paired correlation (3 tumor mice)",
         x = NULL, y = "Pearson r", fill = NULL) +
    paper_theme + theme(legend.position = "top")

  save_png(p, "F15_X1_TGFb_axis_schematic.png", width = 7, height = 5.5)
  TRUE
}

# ---------- F16: X3 Efferocytosis parallel ----------
make_F16 <- function() {
  fig_id <- "F16"
  # liver side
  hsc_f <- file.path(ROOT, "03_sc_DEG", "liver", "HSC_Activated_ECM_deg.csv")
  if (!file.exists(hsc_f)) {
    cand <- Sys.glob(file.path(ROOT, "03_sc_DEG", "liver",
                               "*HSC*Activated*ECM*deg*.csv"))
    if (length(cand) > 0) hsc_f <- cand[1]
  }
  hsc <- safe_read_csv(hsc_f, fig_id)
  pe  <- safe_read_csv(file.path(ROOT, "03_sc_DEG", "brain_clean",
                                 "brain_Tumor_vs_Control",
                                 "Pericytes_deg.csv"), fig_id)
  if (is.null(hsc) || is.null(pe)) return(FALSE)
  hsc <- normalize_deg(hsc, fig_id); pe <- normalize_deg(pe, fig_id)
  if (is.null(hsc) || is.null(pe)) return(FALSE)

  liver_df <- hsc %>% filter(gene %in% c("Gas6", "Mertk", "Axl", "C3")) %>%
    mutate(organ = "Liver HSC_ECM")
  brain_df <- pe %>% filter(gene %in% c("Mfge8", "Bsg", "Calr", "Pecam1")) %>%
    mutate(organ = "Brain Pericytes")

  if (nrow(liver_df) == 0 && nrow(brain_df) == 0) return(FALSE)

  df <- bind_rows(liver_df, brain_df)
  p <- ggplot(df, aes(x = reorder(gene, avg_log2FC), y = avg_log2FC, fill = organ)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = sprintf("padj=%.1e", p_val_adj)),
              hjust = -0.05, size = 2.5, color = "grey30") +
    coord_flip(clip = "off") +
    facet_wrap(~organ, scales = "free_y", ncol = 1) +
    scale_fill_manual(values = c("Liver HSC_ECM" = "#8B4513",
                                 "Brain Pericytes" = "#6A5ACD")) +
    labs(title = "X3 Efferocytosis parallel response (not causal axis)",
         subtitle = "Liver professional phagocyte + Brain non-professional acquisition",
         x = NULL, y = "avg_log2FC") +
    paper_theme + theme(legend.position = "none",
                        plot.margin = margin(5, 35, 5, 5))

  save_png(p, "F16_X3_Efferocytosis_parallel.png", width = 7.5, height = 6)
  TRUE
}

# =============================================================================
# 主调度: 读 config, 按顺序跑
# =============================================================================

log_summ(paste0("Start: ", Sys.time()))
log_summ(paste0("ROOT: ", path.expand(ROOT)))
log_summ(paste0("OUT : ", path.expand(OUT_DIR)))
log_summ("")

# 映射 figure_id -> 函数
fn_map <- list(
  F01 = make_F01, F02 = make_F02, F03 = make_F03, F04 = make_F04,
  F05 = make_F05, F06 = make_F06, F07 = make_F07, F08 = make_F08,
  F09 = make_F09, F10 = make_F10, F11 = make_F11, F12 = make_F12,
  F13 = make_F13, F14 = make_F14, F15 = make_F15, F16 = make_F16
)

# 读 config (允许 config 不存在也跑)
config <- tryCatch(read_csv(CONFIG_CSV, show_col_types = FALSE),
                   error = function(e) NULL)
if (!is.null(config)) {
  log_summ(paste0("Config loaded: ", nrow(config), " figures"))
  config <- config %>% arrange(priority)
} else {
  log_summ("Config not found, running in default order F01-F16")
  config <- tibble(figure_id = names(fn_map), priority = seq_along(fn_map))
}

results <- tibble(figure_id = character(), status = character(),
                  time_sec = numeric())

for (i in seq_len(nrow(config))) {
  fid_full <- config$figure_id[i]
  # 从 F01_X2_thyroid... 提取 F01
  fid <- str_extract(fid_full, "^F\\d+")
  if (is.na(fid) || is.null(fn_map[[fid]])) {
    log_summ(paste0("[SKIP] unknown figure_id: ", fid_full))
    results <- bind_rows(results, tibble(figure_id = fid_full,
                                         status = "skip_unknown",
                                         time_sec = 0))
    next
  }
  log_summ(paste0("--- ", fid_full, " ---"))
  t0 <- Sys.time()
  ok <- tryCatch(fn_map[[fid]](),
                 error = function(e) {
                   log_summ(paste0("[ERR ] ", fid, ": ", e$message))
                   FALSE
                 })
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  st <- if (isTRUE(ok)) "OK" else "FAIL"
  log_summ(paste0("[", st, "] ", fid_full, " (", round(dt, 1), " s)"))
  results <- bind_rows(results, tibble(figure_id = fid_full, status = st,
                                       time_sec = round(dt, 1)))
}

# 汇总
log_summ("")
log_summ("========== SUMMARY ==========")
n_ok <- sum(results$status == "OK")
log_summ(paste0("Total: ", nrow(results),
                "  OK: ", n_ok,
                "  FAIL: ", sum(results$status == "FAIL"),
                "  SKIP: ", sum(results$status == "skip_unknown")))
log_summ("")
log_summ("Per-figure status:")
for (i in seq_len(nrow(results))) {
  log_summ(sprintf("  %-35s  %-5s  %5.1f s",
                   results$figure_id[i], results$status[i],
                   results$time_sec[i]))
}
write_csv(results, file.path(OUT_DIR, "generation_results.csv"))
log_summ(paste0("\nEnd: ", Sys.time()))
log_summ(paste0("Full log: ", LOG_SUMM))
log_summ(paste0("Missing data: ", LOG_MISS))
