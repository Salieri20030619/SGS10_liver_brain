# ============================================================
# Replot HSC_Activated_ECM volcano — 去除 IL-6 通路强调
# 客户更新需求: Il6st/Il6ra/Jak2/Stat3 不再标橙色三角和文字
# 其他保持与原图一致
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
})

# ---- 1. 读取 DEG ----
csv_path <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis/diff_stage1_results/deg_v4_refined/liver/liver_tumor_vs_control_HSC_Activated_ECM.csv"
df <- read.csv(csv_path, stringsAsFactors = FALSE)

# 阈值（与原图一致）
fc_thr  <- 0.20
padj_thr <- 0.10

# ---- 2. 计算 sig 分类（覆盖 csv 里的 sig 列，确保和原图阈值一致）----
df$neglog10_padj <- -log10(df$p_val_adj)
df$sig <- "NS"
df$sig[df$avg_log2FC >  fc_thr & df$p_val_adj < padj_thr] <- "Up"
df$sig[df$avg_log2FC < -fc_thr & df$p_val_adj < padj_thr] <- "Down"
df$sig <- factor(df$sig, levels = c("Down", "NS", "Up"))

n_up   <- sum(df$sig == "Up")
n_down <- sum(df$sig == "Down")
cat(sprintf("Up = %d, Down = %d\n", n_up, n_down))   # 期望 173 / 38

# ---- 3. 选要标签的基因（与原图一致，IL-6 那 4 个已移除）----
label_down <- c("Gm19951","AY036118","Gm48099",
                "Camk1d","Il31ra","Ntm","Nrxn1","Ank3",
                "Cmss1","Egfem1")
label_up   <- c("Col1a1","Col3a1","Col4a1","Sparc","Slit3",
                "Ahnak","Gas6","C3","H2-D1","Hspg2")
label_set  <- c(label_down, label_up)
df$label <- ifelse(df$gene %in% label_set, df$gene, NA)

# ---- 4. y 轴上限（cap 极小 padj，与原图视觉一致）----
y_cap <- 95
df$neglog10_padj_capped <- pmin(df$neglog10_padj, y_cap)

# ---- 5. 配色（与原图一致）----
cols <- c("Down" = "#3B82C4", "NS" = "#BBBBBB", "Up" = "#E04A4A")

# ---- 6. caption（去掉 "orange ▲ = IL-6 pathway"）----
cap <- sprintf("FC>%.2f | padj<%.2f | Up=%d  Down=%d",
               fc_thr, padj_thr, n_up, n_down)

# ---- 7. 绘图 ----
p <- ggplot(df, aes(x = avg_log2FC, y = neglog10_padj_capped, color = sig)) +
  geom_hline(yintercept = -log10(padj_thr), linetype = "dashed",
             color = "black", linewidth = 0.4) +
  geom_vline(xintercept = c(-fc_thr, fc_thr), linetype = "dashed",
             color = "black", linewidth = 0.4) +
  # NS 先画（灰底）
  geom_point(data = subset(df, sig == "NS"),
             size = 1.0, alpha = 0.35, shape = 16) +
  # Up / Down 后画（盖在 NS 上）
  geom_point(data = subset(df, sig != "NS"),
             size = 1.4, alpha = 0.75, shape = 16) +
  # 文字标签（与原图相同的基因集，不含 IL-6 那 4 个）
  geom_text_repel(
    data = subset(df, !is.na(label)),
    aes(label = label),
    size = 3.2,
    max.overlaps = Inf,
    box.padding = 0.35,
    point.padding = 0.20,
    segment.size = 0.25,
    segment.alpha = 0.6,
    min.segment.length = 0.1,
    show.legend = FALSE
  ) +
  scale_color_manual(values = cols, name = NULL) +
  labs(x = "avg_log2FC",
       y = "-log10(p_val_adj)",
       title = cap) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(size = 11, hjust = 0.5,
                              color = "#444444", face = "plain"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.5),
    legend.position = "right",
    legend.key = element_blank(),
    legend.text = element_text(size = 11)
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

# ---- 8. 输出（PDF 矢量 + PNG）----
out_dir <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis/figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
pdf_out <- file.path(out_dir, "Fig_HSC_volcano_no_il6.pdf")
png_out <- file.path(out_dir, "Fig_HSC_volcano_no_il6.png")

ggsave(pdf_out, p, width = 12, height = 9, units = "in", device = cairo_pdf)
ggsave(png_out, p, width = 12, height = 9, units = "in", dpi = 300)

cat("Saved:\n  ", pdf_out, "\n  ", png_out, "\n")
