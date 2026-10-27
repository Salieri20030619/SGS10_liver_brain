#!/usr/bin/env Rscript
# =============================================================================
# SGS10 Task A 补图 - Patch v2: 修复 F01/F02/F03/F16
# Date: 2026-04-20
# Purpose: 上一轮 4 张图 FAIL, 根因: 脚本去 03_sc_DEG/liver/ 找细粒度 CSV
#          但 HSC_Activated_ECM / Hepatocyte_PV pooled DEG 其实在
#          ~/SGS10_analysis/diff_stage1_results/deg_v4_refined/liver/
#          (数字 100% 对齐报告: Acta2=8.36, Serpine1=8.83, F2r=3.91)
# Usage:
#   cd ~/SGS10_analysis/
#   Rscript taskA_figures_v1/patch_v2_fix_F01_F02_F03_F16.R
# =============================================================================

# ---------- 先 source 主脚本, 获得所有工具函数 + 其他 12 张已 OK 的函数 ----------
# 我们只重定义 4 个 FAIL 的函数, 不触发主脚本的 main loop
# 办法: 主脚本 main loop 在最底, 我们只 source 前面的函数定义部分
# 但主脚本结构是一条龙的, 所以最简单: 直接在这里重新 library + 定义所需工具 + 4 个函数 + 局部 dispatcher

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr); library(tidyr)
  library(stringr); library(patchwork); library(scales); library(ggrepel)
  library(purrr)
})

# ---------- 全局配置 (跟主脚本对齐) ----------
ROOT          <- "~/SGS10_analysis/SGS10_deliverables_2026-04-20"
# 新增: pooled DEG 真正来源
ROOT_DEG_LIVER <- "~/SGS10_analysis/diff_stage1_results/deg_v4_refined/liver"
ROOT_DEG_BRAIN <- "~/SGS10_analysis/SGS10_deliverables_2026-04-20/03_sc_DEG/brain_clean/brain_Tumor_vs_Control"
OUT_DIR       <- "~/SGS10_analysis/taskA_figures_2026-04-20_v1"
DPI           <- 300
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

LOG_MISS <- file.path(OUT_DIR, "missing_data_patch_v2.log")
LOG_SUMM <- file.path(OUT_DIR, "generation_summary_patch_v2.log")
writeLines(c(paste("# patch_v2 missing_data.log -", Sys.time()), ""), LOG_MISS)
writeLines(c(paste("# patch_v2 summary.log -", Sys.time()), ""), LOG_SUMM)

# ---------- paper 风格 theme (同主脚本) ----------
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
# deg_v4_refined 的实际列结构: gene, avg_log2FC, pct1, pct2, p_val, p_val_adj, sig
# 已经标准化, 不需要 normalize 函数, 直接用
save_png <- function(p, filename, width = 7, height = 5) {
  path <- file.path(OUT_DIR, filename)
  ggsave(path, p, width = width, height = height, dpi = DPI,
         units = "in", bg = "white")
  cat(paste0("[SAVE] ", filename, " (", width, "x", height, " in @ ", DPI, " DPI)\n"))
}

# =============================================================================
# F01: X2 Thyroid axis (4-panel) [FIXED]
# =============================================================================
# 改动: Hepatocyte_PV 从 deg_v4_refined/liver/ 读; CPECs/Endo 保持 brain_clean
make_F01 <- function() {
  fig_id <- "F01"

  # panel a: Hepatocyte_PV DEG (liver side)
  pv_f <- file.path(ROOT_DEG_LIVER, "liver_tumor_vs_control_Hepatocyte_PV.csv")
  pv <- safe_read_csv(pv_f, fig_id)
  if (is.null(pv)) return(FALSE)

  genes_liver <- c("Ttr", "Thra", "Thrb")
  pv_sub <- pv %>% filter(gene %in% genes_liver) %>%
    mutate(compartment = "Liver Hepatocyte_PV")

  # panel b: CPECs + Endothelial (brain)
  cp <- safe_read_csv(file.path(ROOT_DEG_BRAIN, "CPECs_deg.csv"), fig_id)
  en <- safe_read_csv(file.path(ROOT_DEG_BRAIN, "Endothelial_deg.csv"), fig_id)
  if (is.null(cp) || is.null(en)) return(FALSE)

  brain_rows <- bind_rows(
    cp %>% filter(gene %in% c("Ttr", "Atp1a1", "Folr1")) %>%
      mutate(compartment = "Brain CPECs"),
    en %>% filter(gene == "Thra") %>% mutate(compartment = "Brain Endothelial")
  )

  df_ab <- bind_rows(pv_sub, brain_rows) %>%
    mutate(compartment = factor(compartment,
                                levels = c("Liver Hepatocyte_PV",
                                           "Brain CPECs", "Brain Endothelial")))

  if (nrow(df_ab) == 0) {
    log_miss(fig_id, "No thyroid genes matched in any CSV")
    return(FALSE)
  }

  p_ab <- ggplot(df_ab, aes(x = reorder(gene, avg_log2FC),
                            y = avg_log2FC, fill = compartment)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = 0, color = "grey40", size = 0.3) +
    geom_text(aes(label = sprintf("padj=%.1e", p_val_adj)),
              hjust = -0.1, size = 2.8, color = "grey30") +
    coord_flip(clip = "off") +
    facet_wrap(~compartment, scales = "free_y", ncol = 1) +
    scale_fill_manual(values = c("#8B4513", "#6A5ACD", "#2E8B57")) +
    labs(title = "X2 Thyroid hormone cross-organ axis",
         subtitle = "Liver Hepatocyte_PV (Ttr) -> Brain CPECs (Ttr padj 9.6e-10) -> Endo (Thra padj 2.5e-06)",
         x = NULL, y = "avg_log2FC (Tumor vs Control)") +
    paper_theme + theme(legend.position = "none",
                        strip.text = element_text(face = "bold", size = 10),
                        plot.margin = margin(5, 40, 5, 5))

  save_png(p_ab, "F01_X2_thyroid_axis_4panel.png", width = 9, height = 8)
  TRUE
}

# =============================================================================
# F02: L1 HSC activation (Top 10) [FIXED]
# =============================================================================
# 改动: 路径改为 deg_v4_refined/liver/
make_F02 <- function() {
  fig_id <- "F02"
  f <- file.path(ROOT_DEG_LIVER,
                 "liver_tumor_vs_control_HSC_Activated_ECM.csv")
  hsc <- safe_read_csv(f, fig_id)
  if (is.null(hsc)) return(FALSE)

  genes <- c("Acta2", "Serpine1", "Col1a1", "F2r", "Tagln",
             "Pmepa1", "Col1a2", "Smad3", "Lox", "Loxl1")
  sub <- hsc %>% filter(gene %in% genes) %>%
    mutate(sig_flag = ifelse(p_val_adj < 0.05, "sig", "ns"))

  if (nrow(sub) == 0) {
    log_miss(fig_id, "No L1 genes matched")
    return(FALSE)
  }

  p <- ggplot(sub, aes(x = reorder(gene, avg_log2FC),
                       y = avg_log2FC, fill = sig_flag)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = sprintf("padj=%.1e", p_val_adj)),
              hjust = -0.05, size = 2.8, color = "grey30") +
    coord_flip(clip = "off") +
    scale_fill_manual(values = c(sig = col_up, ns = col_ns),
                      labels = c("sig (padj<0.05)", "ns")) +
    labs(title = "L1 HSC myofibroblast activation (sc DEG Top 10)",
         subtitle = "HSC_Activated_ECM: Tumor vs Control (Acta2 log2FC 8.36, Serpine1 8.83)",
         x = NULL, y = "avg_log2FC", fill = NULL) +
    paper_theme + theme(legend.position = "top",
                        plot.margin = margin(5, 40, 5, 5))

  save_png(p, "F02_L1_HSC_activation.png", width = 8, height = 5.5)
  TRUE
}

# =============================================================================
# F03: L4 HSC-Immune hub [FIXED - 改用 Kupffer 的 myeloid marker DEG]
# =============================================================================
# 改动: 原来找 HSC_ECM 里的 Csf1r/Mertk 等 myeloid marker 不合逻辑 (HSC 不是髓系)
#       改用 Kupffer (髓系, 有 CSV) 的 DEG, 再加一个 HSC 侧的 ligand (Csf1, Gas6)
make_F03 <- function() {
  fig_id <- "F03"

  # HSC 侧: 它是 ligand source, 看 Csf1 / Gas6 / C3 等
  hsc_f <- file.path(ROOT_DEG_LIVER,
                     "liver_tumor_vs_control_HSC_Activated_ECM.csv")
  kup_f <- file.path(ROOT_DEG_LIVER,
                     "liver_tumor_vs_control_Kupffer.csv")
  hsc <- safe_read_csv(hsc_f, fig_id)
  kup <- safe_read_csv(kup_f, fig_id)
  if (is.null(hsc) || is.null(kup)) return(FALSE)

  hsc_ligands <- hsc %>% filter(gene %in% c("Csf1", "Gas6", "C3", "Cxcl12")) %>%
    mutate(role = "Ligand in HSC_ECM (source)")
  kup_receptors <- kup %>% filter(gene %in% c("Csf1r", "Mertk", "Axl",
                                              "Cd74", "Cd68", "Mrc1",
                                              "Itgam", "Itgax")) %>%
    mutate(role = "Receptor/Marker in Kupffer (target)")

  df <- bind_rows(hsc_ligands, kup_receptors)
  if (nrow(df) == 0) {
    log_miss(fig_id, "No L4 genes matched in HSC_ECM or Kupffer")
    return(FALSE)
  }

  p <- ggplot(df, aes(x = reorder(gene, avg_log2FC),
                      y = avg_log2FC, fill = role)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = sprintf("padj=%.1e", p_val_adj)),
              hjust = -0.05, size = 2.5, color = "grey30") +
    coord_flip(clip = "off") +
    facet_wrap(~role, scales = "free_y", ncol = 1) +
    scale_fill_manual(values = c("Ligand in HSC_ECM (source)" = "#8e44ad",
                                 "Receptor/Marker in Kupffer (target)" = "#2980b9")) +
    labs(title = "L4 HSC-Immune hub: source ligand + target receptor",
         subtitle = "See Task B Figure 2b for spatial validation (Cd74/Csf1r/Mrc1/Pecam1)",
         x = NULL, y = "avg_log2FC (Tumor vs Control)") +
    paper_theme + theme(legend.position = "none",
                        plot.margin = margin(5, 40, 5, 5))

  save_png(p, "F03_L4_HSC_immune_hub.png", width = 9, height = 6.5)

  # Bonus: 顺便把 viz/ 下现成的 HSC_Activated_ECM volcano + heatmap 复制一份过来
  viz_candidates <- c(
    "~/SGS10_analysis/diff_stage1_results/viz/01_volcano/liver/liver_tumor_vs_control_HSC_Activated_ECM_volcano.png",
    "~/SGS10_analysis/diff_stage1_results/viz/02_heatmap/liver/liver_tumor_vs_control_HSC_Activated_ECM_heatmap.png"
  )
  for (src in viz_candidates) {
    src_exp <- path.expand(src)
    if (file.exists(src_exp)) {
      dst <- file.path(OUT_DIR, paste0("F03_bonus_", basename(src_exp)))
      file.copy(src_exp, dst, overwrite = TRUE)
      cat(paste0("[COPY] ", basename(src_exp), " -> ", basename(dst), "\n"))
    }
  }
  TRUE
}

# =============================================================================
# F16: X3 Efferocytosis parallel [FIXED]
# =============================================================================
# 改动: liver 侧 HSC_ECM 从 deg_v4_refined/liver/ 读; brain 侧 Pericytes 保持 brain_clean
make_F16 <- function() {
  fig_id <- "F16"
  hsc <- safe_read_csv(file.path(ROOT_DEG_LIVER,
                                 "liver_tumor_vs_control_HSC_Activated_ECM.csv"),
                       fig_id)
  pe  <- safe_read_csv(file.path(ROOT_DEG_BRAIN, "Pericytes_deg.csv"),
                       fig_id)
  if (is.null(hsc) || is.null(pe)) return(FALSE)

  liver_df <- hsc %>% filter(gene %in% c("Gas6", "Mertk", "Axl", "C3")) %>%
    mutate(organ = "Liver HSC_Activated_ECM")
  brain_df <- pe %>% filter(gene %in% c("Mfge8", "Bsg", "Calr", "Pecam1")) %>%
    mutate(organ = "Brain Pericytes")

  df <- bind_rows(liver_df, brain_df)
  if (nrow(df) == 0) {
    log_miss(fig_id, "No efferocytosis genes matched")
    return(FALSE)
  }

  p <- ggplot(df, aes(x = reorder(gene, avg_log2FC),
                      y = avg_log2FC, fill = organ)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = sprintf("padj=%.1e", p_val_adj)),
              hjust = -0.05, size = 2.5, color = "grey30") +
    coord_flip(clip = "off") +
    facet_wrap(~organ, scales = "free_y", ncol = 1) +
    scale_fill_manual(values = c("Liver HSC_Activated_ECM" = "#8B4513",
                                 "Brain Pericytes"         = "#6A5ACD")) +
    labs(title = "X3 Efferocytosis parallel response (not causal axis)",
         subtitle = "Liver professional phagocyte + Brain non-professional acquisition",
         x = NULL, y = "avg_log2FC") +
    paper_theme + theme(legend.position = "none",
                        plot.margin = margin(5, 40, 5, 5))

  save_png(p, "F16_X3_Efferocytosis_parallel.png", width = 8, height = 6.5)
  TRUE
}

# =============================================================================
# Mini dispatcher - 只跑 4 张 FIX
# =============================================================================
log_summ(paste0("Start patch_v2: ", Sys.time()))
log_summ(paste0("ROOT_DEG_LIVER: ", path.expand(ROOT_DEG_LIVER)))
log_summ(paste0("OUT_DIR:        ", path.expand(OUT_DIR)))
log_summ("")

fixes <- list(
  F01 = make_F01,
  F02 = make_F02,
  F03 = make_F03,
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
log_summ("========== PATCH v2 SUMMARY ==========")
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
