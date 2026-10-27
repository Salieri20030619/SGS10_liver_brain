# ==================== 脑 Marker 库 v4.1 (rewrite) ====================
# 文件: ~/SGS10_analysis/scripts/marker_lib_brain_v4.R
#
# v4.1 相比 v4.0 的核心设计变化:
#   1. 每个子类 marker 集强制包含 pan-lineage anchor 基因
#      (例: 所有 Excitatory_* 必须含 Slc17a7 或 Camk2a;
#            所有 Astrocyte_* 必须含 Aqp4/Aldh1l1/Gja1)
#   2. 删除 Mfge8/Excitatory_DL/Excitatory_Thal/Astro_Hippocampal/Astro_Thalamic
#      这些在 winner-takes-all 打分下会跨谱系通吃的"伪子类"集合
#      (诊断: v4.0 brain3 有 42% 细胞被误翻, 根因是这些集合缺乏谱系锚定)
#   3. 保留 v4 里真正可靠的亚类区分: Hpc_CA1/CA2/CA3/DG, L23/L45/L6,
#      Pvalb/Sst/Vip/Meis2/Lamp5, Microglia_Homeostatic/DAM, Oligo_Early/Pre/Mature
#
# 来源:
#   - Zeisel et al. 2015, Science (cortex/hippocampus)
#   - Saunders et al. 2018, Cell (Drop-seq mouse brain)
#   - Tasic et al. 2018, Nature (cortical interneuron taxonomy)
#   - Marques et al. 2016, Science (oligo lineage)
#   - Keren-Shaul et al. 2017, Cell (DAM signature)
# ============================================================

# pan-lineage anchor 基因 (在每个子类集合里必须出现 >= 2 个)
# 这些是 lineage gate 机制参考的锚定集,同时也混在每个子集的 marker 里
# 避免全局打分时跨谱系通吃
PAN_NEURON      <- c("Snap25", "Syp", "Stmn2", "Rbfox3")
PAN_EXCITATORY  <- c("Slc17a7", "Camk2a", "Slc17a6")
PAN_INHIBITORY  <- c("Gad1", "Gad2", "Slc32a1")
PAN_ASTROCYTE   <- c("Aqp4", "Aldh1l1", "Gja1", "Slc1a3")
PAN_OLIGO       <- c("Olig1", "Olig2", "Sox10", "Plp1")
PAN_MICROGLIA   <- c("C1qa", "C1qb", "Cx3cr1", "Csf1r", "Hexb")
PAN_ENDOTHELIAL <- c("Pecam1", "Cldn5", "Cdh5")

brain_markers_v4 <- list(

  # ==================== 神经元: 兴奋性皮层各层 ====================
  # 每集都强制含 Slc17a7 + Camk2a (pan-excitatory anchor)
  # v4.1 第二轮: 均衡每个子集的 specific 基因数到 5-6 个

  # L2/3: Cux1/Cux2 是经典 upper-layer marker
  Excitatory_L23 = c(
    "Slc17a7", "Camk2a", "Snap25",                   # anchor
    "Cux2", "Cux1", "Lamp5", "Rasgrf2", "Calb1"      # L2/3 specific
  ),

  # L4/5: Rorb 是 L4, Fezf2/Bcl11b 是 L5
  Excitatory_L45 = c(
    "Slc17a7", "Camk2a", "Snap25",                   # anchor
    "Rorb", "Fezf2", "Bcl11b", "Etv1", "Scnn1a"      # L4/5 specific
  ),

  # L6: Foxp2/Tle4/Sulf1
  Excitatory_L6 = c(
    "Slc17a7", "Camk2a", "Snap25",                   # anchor
    "Foxp2", "Tle4", "Sulf1", "Syt6", "Rprm"         # L6 specific
  ),

  # ==================== 神经元: 海马 ====================

  # CA1: Wfs1/Pou3f1/Fibcd1/Spink8
  Hpc_CA1 = c(
    "Slc17a7", "Camk2a", "Snap25",           # anchor
    "Wfs1", "Pou3f1", "Fibcd1", "Spink8",    # CA1 specific
    "Mpped1", "Dkk3"
  ),

  # CA2: 扩充到 6 个 specific (Lein 2005, Cembrowski 2016)
  Hpc_CA2 = c(
    "Slc17a7", "Camk2a", "Snap25",                   # anchor
    "Rgs14", "Amigo2", "Map3k15", "Pcp4",            # 原 CA2 specific
    "Cacng5", "Str6", "Cadps2"                       # 新增 specific
  ),

  # CA3: 扩充到 6 个 specific (Cembrowski 2016 hippocampal heterogeneity)
  # 之前只有 3 个 specific 造成短集合偷跑, 抢走 4815 个 CA1 cells
  Hpc_CA3 = c(
    "Slc17a7", "Camk2a", "Snap25",                   # anchor
    "Cpne7", "Iyd", "Necab2",                        # 原 CA3 specific
    "Chgb", "Gpr161", "Kcnh7", "Prkcd"               # 新增 specific
  ),

  # DG: 扩充到 7 个 specific (Prox1 是决定性 granule cell marker)
  Hpc_DG = c(
    "Slc17a7", "Camk2a", "Snap25",                   # anchor
    "Prox1", "C1ql2", "Stxbp6", "Tdo2",              # 原 DG specific
    "Dsp", "Fam19a2", "Rasl10a"                      # 新增 specific
  ),

  # ==================== 神经元: 抑制性 ====================
  # 每集都强制含 Gad1 + Gad2 (pan-inhibitory anchor)
  # v4.1 第二轮: 扩充 specific 到 5-6 个,避免短集合偷跑
  # (Tasic 2018 Nature cortical interneuron taxonomy)

  Inhibitory_Pvalb = c(
    "Gad1", "Gad2", "Slc32a1", "Snap25",             # anchor
    "Pvalb", "Tac1", "Pthlh", "Cox6a2", "Btbd11"     # Pvalb specific
  ),

  Inhibitory_Sst = c(
    "Gad1", "Gad2", "Slc32a1", "Snap25",             # anchor
    "Sst", "Chodl", "Nos1", "Tacr3", "Mme", "Crhbp"  # Sst specific
  ),

  Inhibitory_Vip = c(
    "Gad1", "Gad2", "Slc32a1", "Snap25",             # anchor
    "Vip", "Penk", "Crh", "Prok2", "Tac2"            # Vip specific
  ),

  # Meis2 亚型: 嗅球及腹侧 GABAergic
  Inhibitory_Meis2 = c(
    "Gad1", "Gad2", "Slc32a1", "Snap25",             # anchor
    "Meis2", "Pax6", "Sp8", "Dlx1", "Dlx5"           # Meis2 specific
  ),

  # Lamp5 + 5HT3a 中间神经元
  Inhibitory_Lamp5 = c(
    "Gad1", "Gad2", "Slc32a1", "Snap25",             # anchor
    "Lamp5", "Htr3a", "Ndnf", "Ntn1", "Reln"         # Lamp5 specific
  ),

  # ==================== 胶质细胞: Astrocyte ====================
  # 所有 astrocyte 集合必须含 Aqp4 + Aldh1l1 + Gja1
  #
  # 注意: v4.0 里的 Astrocyte_Mfge8 / Astrocyte_Hippocampal / Astrocyte_Thalamic
  # 全部删除。原因:
  #   - Mfge8 在脑里不是 astrocyte 特异,是通吃元凶
  #   - 海马/丘脑 astrocyte 没有可靠的 regional identity marker in scRNA
  # 替代: 只保留两个按激活态分的亚型

  # 通用 / 静息态 astrocyte (Zeisel Astro2 的大部分)
  Astrocyte_Resting = c(
    "Aqp4", "Aldh1l1", "Gja1", "Slc1a3",     # anchor
    "Slc1a2", "Glul", "Aldoc", "S100b"       # astrocyte specific
  ),

  # 反应性 / GFAP-high astrocyte (Zeisel Astro1, border/white matter)
  Astrocyte_Gfap = c(
    "Aqp4", "Aldh1l1", "Gja1", "Slc1a3",     # anchor
    "Gfap", "Vim", "Lcn2"                    # reactive specific
  ),

  # ==================== 胶质细胞: Oligodendrocyte ====================
  # 所有 oligo 集合必须含 Olig1 + Olig2 + Sox10 + Plp1

  # Early / pre-OL: 未开始大量髓鞘化
  Oligo_Pre = c(
    "Olig1", "Olig2", "Sox10",               # anchor
    "Bmp4", "Gpr17", "Nkx2-2", "Itpr2"       # pre-OL specific
  ),

  # Mature / myelinating: MBP/MOG
  Oligo_Mature = c(
    "Olig1", "Olig2", "Sox10", "Plp1",       # anchor
    "Mbp", "Mog", "Mag", "Mobp",             # myelin specific
    "Cnp", "Aspa", "Klk6"
  ),

  # ==================== OPC ====================
  # OPC 必须含 Pdgfra 和 Cspg4,不走 oligo anchor
  OPCs = c(
    "Pdgfra", "Cspg4", "Sox10",              # OPC specific anchor
    "Olig1", "Olig2", "Vcan", "Cspg5"
  ),

  # ==================== 免疫细胞: Microglia ====================
  # 所有 microglia 集合必须含 C1qa + C1qb + Cx3cr1 + Csf1r

  # Homeostatic: 静息态小胶质 (Keren-Shaul 2017)
  Microglia_Homeostatic = c(
    "C1qa", "C1qb", "Cx3cr1", "Csf1r", "Hexb",     # microglia anchor
    "Tmem119", "P2ry12", "Selplg", "Sall1"          # homeostatic specific
  ),

  # DAM (Disease-Associated Microglia, Keren-Shaul 2017 signature)
  # 必须同时含 homeostatic 和 DAM marker,因为 DAM 是激活态但仍是 microglia
  Microglia_DAM = c(
    "C1qa", "C1qb", "Cx3cr1", "Csf1r",       # microglia anchor
    "Trem2", "Apoe", "Cd9", "Spp1",          # DAM core signature
    "Itgax", "Lpl", "Cst7", "Axl"
  ),

  # Perivascular macrophage: Lyve1+/Mrc1+ (区别于 microglia 的 Lyve1-)
  PerivascularMac = c(
    "C1qa", "C1qb", "Csf1r",                 # macrophage anchor
    "Lyve1", "Mrc1", "Cd163", "F13a1"        # PVM specific
  ),

  # ==================== 血管/支持 ====================

  # Endothelial: 强锚定 Pecam1/Cldn5/Cdh5
  Endothelial = c(
    "Pecam1", "Cldn5", "Cdh5", "Flt1",       # endothelial anchor
    "Kdr", "Tek", "Ly6c1", "Slco1a4"
  ),

  # Pericyte: Pdgfrb + Rgs5 是核心
  Pericytes = c(
    "Pdgfrb", "Rgs5", "Vtn",                 # pericyte anchor
    "Mcam", "Higd1b", "Atp13a5"
  ),

  # VLMC (vascular leptomeningeal cells): Dcn/Lum + Col1a1
  # 注意: 不能用 Pdgfra 作 VLMC anchor,因为和 OPC 共享
  VLMC = c(
    "Dcn", "Lum", "Col1a1", "Col1a2",        # VLMC specific
    "Slc6a13", "Ogn"
  ),

  # ==================== 脉络丛/室管膜 ====================

  # CPEC: Ttr 是决定性的 transthyretin,脉络丛 epithelial 专属
  CPECs = c(
    "Ttr", "Folr1", "Kl",                    # CPEC core
    "Clic6", "Otx2", "Aqp1"
  ),

  # Ependymal: Foxj1 + 纤毛基因
  Ependymal = c(
    "Foxj1", "Rfx2",                         # master TF
    "Spag6", "Ccdc153", "Pifo", "Tekt1", "Dnah11"
  )

  # ==================== 已删除的集合 ====================
  # v4.0 有但 v4.1 删掉的:
  #   - Excitatory_L6b       (marker 和 L6 完全重叠,区分不了)
  #   - Excitatory_DL        (Foxp2/Bcl11b/Tbr1 不特异,通吃)
  #   - Excitatory_Thal      (Pvalb/Calb2/Cck 全是中间神经元 marker,完全错)
  #   - Astrocyte_Mfge8      (Mfge8 不是 astrocyte 特异,通吃元凶)
  #   - Astrocytes (通用)    (catch-all,被 filter 掉,走 passthrough)
  #   - Astrocyte_Hippocampal(无可靠 regional marker)
  #   - Astrocyte_Thalamic   (无可靠 regional marker)
  #   - Oligo_Early          (和 Oligo_Pre 重叠,合并到 Pre)
  #   - Oligodendrocytes     (catch-all,被 filter 掉,走 passthrough)
  #   - Microglia            (catch-all,被 filter 掉,走 passthrough)
  #   - Neuron_Excitatory    (catch-all)
  #   - Neuron_Inhibitory    (catch-all)
  #   - Neurons              (catch-all)
  #   - NK_brain             (太罕见,4 marker 会在 migrated T cell 上虚高)
  #   - Other                (垃圾桶)
)


# ==================== 元数据 ====================

brain_marker_min_genes <- 3

brain_marker_sources <- c(
  "Zeisel et al. 2015 Science (mouse cortex/hippocampus)",
  "Saunders et al. 2018 Cell (Drop-seq adult mouse brain)",
  "Tasic et al. 2018 Nature (cortical interneuron taxonomy)",
  "Marques et al. 2016 Science (oligodendrocyte lineage)",
  "Keren-Shaul et al. 2017 Cell (DAM signature)"
)

# ==================== Lineage gate 结构 ====================
# 每个 v4 子类所属的 broad lineage
# 主脚本在 rescore 时只允许 v3.2 名字翻到同 lineage 的兄弟类
# 这是 v4.1 的关键工程防护

brain_lineage_map <- list(
  # 兴奋性神经元
  Excitatory_L23         = "excit_neuron",
  Excitatory_L45         = "excit_neuron",
  Excitatory_L6          = "excit_neuron",
  # 海马兴奋性
  Hpc_CA1                = "excit_neuron",
  Hpc_CA2                = "excit_neuron",
  Hpc_CA3                = "excit_neuron",
  Hpc_DG                 = "excit_neuron",
  # 抑制性神经元
  Inhibitory_Pvalb       = "inhib_neuron",
  Inhibitory_Sst         = "inhib_neuron",
  Inhibitory_Vip         = "inhib_neuron",
  Inhibitory_Meis2       = "inhib_neuron",
  Inhibitory_Lamp5       = "inhib_neuron",
  # Astrocyte
  Astrocyte_Resting      = "astrocyte",
  Astrocyte_Gfap         = "astrocyte",
  # Oligo / OPC
  Oligo_Pre              = "oligo",
  Oligo_Mature           = "oligo",
  OPCs                   = "opc",
  # Microglia
  Microglia_Homeostatic  = "microglia",
  Microglia_DAM          = "microglia",
  PerivascularMac        = "pvm",
  # 血管
  Endothelial            = "endo",
  Pericytes              = "pericyte",
  VLMC                   = "vlmc",
  # CSF
  CPECs                  = "cpec",
  Ependymal              = "ependymal"
)

# v3.2 名字 → lineage (用于 gate;如果 v3.2 名字不在这里,gate 不生效,走普通 flow)
brain_v32_lineage_map <- list(
  # 兴奋性 (v3.2 命名)
  Excitatory_L23         = "excit_neuron",
  Excitatory_L45         = "excit_neuron",
  Excitatory_L6          = "excit_neuron",
  Hpc_CA1                = "excit_neuron",
  Hpc_CA2                = "excit_neuron",
  Hpc_CA3                = "excit_neuron",
  Hpc_DG                 = "excit_neuron",
  Neurons                = "excit_neuron",  # v3.2 粗 neurons 归到兴奋性(保守)
  # 抑制性
  Inhibitory_Pvalb       = "inhib_neuron",
  Inhibitory_Sst         = "inhib_neuron",
  Inhibitory_Vip         = "inhib_neuron",
  Inhibitory_Meis2       = "inhib_neuron",
  # Astrocyte (v3.2 只有 Astrocytes 粗类)
  Astrocytes             = "astrocyte",
  # Oligo
  Oligodendrocytes       = "oligo",
  OPCs                   = "opc",
  # Microglia
  Microglia_Homeostatic  = "microglia",
  Microglia_DAM          = "microglia",
  Microglia              = "microglia",
  Microglia_M1           = "microglia",
  Microglia_M2           = "microglia",
  # 血管
  Endothelial            = "endo",
  Pericytes              = "pericyte",
  VLMC                   = "vlmc",
  # CSF
  CPECs                  = "cpec",
  Ependymal              = "ependymal"
)


cat("[marker_lib_brain_v4.1] Loaded:\n")
cat("  - Cell types:", length(brain_markers_v4), "\n")
cat("  - Lineages:  ", length(unique(unlist(brain_lineage_map))), "\n")
cat("  - Sources:\n")
for (s in brain_marker_sources) cat("      *", s, "\n")
