# ============================================================
# aliases_v4.R  --  v3.2 → v4 命名归一化 + 免疫系精修配置
# ============================================================
# 这个文件只定义查找表,不做任何计算。
# 主脚本 secondary_annotation_v4.R 会 source 它。
#
# 修改原则:
#   - 加新 alias 只需在下面对应的 list 加一行,不改主脚本
#   - alias 的方向: v3.2 名字 -> v4 marker 库里的标准名字
#   - 精修规则: 指定某个 v3.2 粗类应该在哪些 v4 子类之间打分
# ============================================================


# ------------------------------------------------------------
# 1) ALIAS 表: v3.2 名 -> v4 标准名
#    这是"纯改名",不会触发重新打分
# ------------------------------------------------------------

# --- 肝鼠 ---
alias_liver_mouse <- c(
  # macrophage 家族 -> 统一到 Kupffer_KC
  "KC_Homeostatic"    = "Kupffer_KC",
  "Kupffer"           = "Kupffer_KC",
  # DC 命名
  "cDC1"              = "DC_cDC1",
  "cDC2"              = "DC_cDC2",
  "pDC"               = "DC_pDC",
  # HSC 未分型 -> 静息态
  "Stellate"          = "HSC_Quiescent",
  # T 细胞命名
  "CD4_T"             = "T_CD4",
  "CD8_T"             = "T_CD8",
  "Treg"              = "T_Treg",
  "NKT"               = "T_NKT",
  # 复数->单数
  "Pericytes"         = "Pericyte",
  "Neutrophils"       = "Neutrophil",
  "Erythrocytes"      = "Erythrocyte",
  "Basophils"         = "Basophil",
  # 大血管内皮 (v3.2 的 Endothelial 非 LSEC)
  "Endothelial"       = "LVEC",
  # Monocytes 走 MDM 路径
  "Monocytes"         = "Macrophage_MDM",
  # LSEC_Mid 并到通用 LSEC
  "LSEC_Mid"          = "LSEC",
  # 复数 -> 单数
  "NK_cells"          = "NK",
  "Cholangiocytes"    = "Cholangiocyte"
  # 说明: LCM 在主脚本里不做 alias,保留原名,由 mapping CSV 归到 Macrophage_M1M2 broad
  #       因为 v4 marker 库里没有对应的子集,score 验证会跳过它
)

# --- 肝人 ---
# 注: PDX 人源样本里 v3.2 的名字多数不带 _h 后缀,
#     v4 库统一用 _h 后缀,alias 补上
alias_liver_human <- c(
  "Tumor_Hepatocyte"  = "Tumor_Hepatocyte_like",
  "Hep_Periportal"    = "Hep_Periportal_h",
  "Hep_Perivenous"    = "Hep_Perivenous_h",
  "HSC_Activated"     = "HSC_Activated_h",
  "Stellate"          = "HSC_Activated_h",
  "Endothelial"       = "Endothelial_h",
  "LSEC"              = "LSEC_h",
  "Cholangiocyte"     = "Cholangiocyte_h",
  "CD4_T"             = "T_cell_h",
  "CD8_T"             = "T_cell_h",
  "T_cells"           = "T_cell_h",
  "B_cells"           = "B_cell_h",
  "NK_cells"          = "NK_cell_h",
  "Macrophage_NAM"    = "Macrophage_h",
  "Kupffer_Macrophage"= "Macrophage_h",
  "Monocytes"         = "Macrophage_h",
  "Plasma"            = "Plasma_h",
  "B_Plasma"          = "Plasma_h",
  "Neutrophils"       = "Neutrophil_h"
)

# --- 脑 ---
# v4.1-c:
#   - Astrocytes 粗类 -> refinement_rules 精修 {Resting, Gfap} (见下方 rules)
#     不在 alias 里,因为 refinement 比 alias 优先
#   - Oligodendrocytes 粗类 alias 到 Oligo_Mature (绝大多数脑 Oligo 是成熟态)
#   - Microglia 粗类 -> refinement_rules 精修 {Homeostatic, DAM}
#   - Neurons / 神经元亚类 / OPCs / Microglia_M1/M2 -> passthrough
alias_brain <- c(
  "Pericyte"          = "Pericytes",
  "Oligodendrocytes"  = "Oligo_Mature"
)


# ------------------------------------------------------------
# 2) 免疫精修规则: v3.2 粗类 -> v4 子空间
#    用户选择的策略: 只对免疫系 (DC/T/B + 脑小胶质) 做精修
# ------------------------------------------------------------
#
# 格式: list(
#   "v3.2 粗类名" = c("v4 子类 1", "v4 子类 2", ...)
# )
# 主脚本会:
#   1) 找到 final_celltype == 粗类名 的细胞
#   2) 只用子空间里的 marker 集做 AddModuleScore
#   3) winner-takes-all,细胞的 celltype_refined 变成 top1 子类
#   4) 标记 celltype_refinement_source = "subtype_rescore"

refinement_rules <- list(

  # --- 肝鼠 ---
  "Dendritic" = c("DC_cDC1", "DC_cDC2", "DC_pDC"),
  "T_cells"   = c("T_CD4", "T_CD8", "T_Treg", "T_NKT"),
  "B_cells"   = c("B_Naive", "B_Memory", "B_Plasma"),

  # --- 脑 ---
  # v3.2 已经有 Microglia_Homeostatic / Microglia_DAM 两个细分,
  # 剩下的粗 Microglia (con_brain 241 cells, brain3 193 cells) 才精修
  "Microglia" = c("Microglia_Homeostatic", "Microglia_DAM"),

  # v3.2 的 Astrocytes 粗类 -> 区分静息 / Gfap+ 反应态
  # PDX 脑的肝癌远程影响可能诱导反应性 astrogliosis,
  # Gfap+ 比例变化是潜在的生物学指标
  "Astrocytes" = c("Astrocyte_Resting", "Astrocyte_Gfap")
)


# ------------------------------------------------------------
# 3) Score 验证豁免名单
#    这些 v3.2 名字既不在 alias 也不在 refinement,也不在 v4 marker 库,
#    但我们知道它们存在且不应该被 "翻译" 或 "重打分",
#    直接 pass-through 到 broad collapse 阶段。
# ------------------------------------------------------------

passthrough_names <- c(
  # 肝鼠: LCM 保留,让 mapping 归到 Macrophage_M1M2
  "LCM",

  # ============ 脑的 passthrough 策略 ============
  # 研究目标: 肝癌 PDX 对脑的远程影响
  # -> 核心关注免疫(Microglia)+ 血管 + 屏障
  # -> 神经元和大胶质作为背景, 完全信任 v3.2 聚类, 不做 score 验证
  #
  # v4 在脑里只精修: Microglia_Homeostatic/DAM/PerivascularMac,
  #                  Endothelial, Pericytes, VLMC, CPECs, Ependymal
  # 其他脑细胞类 -> passthrough (只做 alias + broad)

  # 神经元 - 皮层兴奋性
  "Excitatory_L23", "Excitatory_L45", "Excitatory_L6",
  # 神经元 - 海马
  "Hpc_CA1", "Hpc_CA2", "Hpc_CA3", "Hpc_DG",
  # 神经元 - 抑制性
  "Inhibitory_Pvalb", "Inhibitory_Sst", "Inhibitory_Vip",
  "Inhibitory_Meis2", "Inhibitory_Lamp5",
  # 神经元 - 兜底
  "Neurons",
  # 大胶质 (Astrocyte/Oligo 粗类走 alias, 不在这里;
  # 但如果 v3.2 以后出现更细的亚型名字没在库里, 这里兜底)
  "OPCs",

  # 脑: Microglia_M1/M2 保留 (v4 库里没有对应),
  # mapping CSV 归到 Microglia broad
  "Microglia_M1", "Microglia_M2",

  # 通用兜底
  "Other", "Unknown"
  # 说明:
  # - Astrocytes / Oligodendrocytes 走 alias 到 Astrocyte_Resting / Oligo_Mature
  # - 不把 Microglia_Homeostatic/DAM 放这里 -> 它们要走 score 验证
  # - 不把 Endothelial/Pericytes/VLMC/CPECs/Ependymal 放这里 -> 要走 score 验证
)


# ------------------------------------------------------------
# 4) Catch-all 集名单: 从 AddModuleScore 里排除掉,
#    避免长且通用的基因集通吃子分型
#    (验证报告 Section 4 诊断出的 100% overlap 问题)
# ------------------------------------------------------------

catchall_sets_brain <- character(0)  # v4.1 库已重写,全部集合都有谱系锚定,无 catch-all

catchall_sets_liver_mouse <- character(0)  # LSEC 现参与打分,做 score 验证

catchall_sets_liver_human <- character(0)  # 人源库里没有明显的 catch-all


# ------------------------------------------------------------
# 5) 打印加载信息
# ------------------------------------------------------------
cat("[aliases_v4] Loaded:\n")
cat("  - liver_mouse aliases:", length(alias_liver_mouse), "\n")
cat("  - liver_human aliases:", length(alias_liver_human), "\n")
cat("  - brain aliases:      ", length(alias_brain), "\n")
cat("  - refinement rules:   ", length(refinement_rules), "\n")
cat("  - passthrough names:  ", length(passthrough_names), "\n")
cat("  - catchall (brain):   ", length(catchall_sets_brain), "\n")
cat("  - catchall (liver_m): ", length(catchall_sets_liver_mouse), "\n")
