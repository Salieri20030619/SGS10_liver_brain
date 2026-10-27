# ==================== 肝脏 Marker 库 v4 ====================
# 文件: ~/SGS10_analysis/scripts/marker_lib_liver_v4.R
# 来源:
#   - Xiong et al. 2019, Mol Cell. "Landscape of intercellular crosstalk in
#     healthy and NASH liver revealed by single-cell secretome gene analysis"
#   - Bai et al. 2023, Mil Med Res. "Single-cell transcriptomic dissection of
#     the cellular and molecular events underlying the triclosan-induced liver
#     fibrosis in mice"
# 用途: 二次注释,用 AddModuleScore 给每个 cell 打分,分配最匹配的细分亚型
# ============================================================

# ==================== Mouse 肝 marker (鼠源 PDX TME + 纯鼠对照) ====================

liver_markers_mouse_v4 <- list(

  # ----- Hepatocyte 分区 (Halpern 2017 + Xiong 2019) -----
  Hep_Periportal = c("Hal", "Sds", "Cyp2f2", "Cdh1", "Pck1", "G6pc",
                     "Cps1", "Sdhd", "Arg1", "Mup3", "Ass1", "Cyp2a4"),

  Hep_Perivenous = c("Glul", "Cyp2e1", "Cyp1a2", "Oat", "Axin2", "Lgr5",
                     "Notum", "Lhpp", "Slc1a2", "Cyp7a1", "Cyp2c37"),

  Hep_Proliferating = c("Mki67", "Top2a", "Pcna", "Ccnb1", "Mcm2",
                        "Ccnd1", "Cenpf"),

  # ----- HSC 8 亚型 (Bai 2023, Fig 4) -----
  # qHSC: 静息态,Lrat 是肝星状细胞经典 marker
  HSC_Quiescent = c("Lrat", "Angptl6", "Vipr1", "Rgs5", "Reln", "Pth1r"),

  # aHSC1: ECM 主力,胶原蛋白高表达
  HSC_Activated_ECM = c("Col1a1", "Col1a2", "Col3a1", "Col5a2", "Col6a1",
                        "Col6a2", "Col14a1", "Acta2", "Fmod", "Fbln5"),

  # aHSC2: ECM remodeling,protease 和 inhibitor
  HSC_Activated_Remodel = c("Timp1", "Timp2", "Aebp1", "Ddr2", "Bmpr1a",
                            "Sparc", "Mmp2"),

  # aHSC4: chemotactic,招募免疫细胞
  HSC_Activated_Chemo = c("Ccl6", "Cxcl10", "Ccl24", "Ccr2", "Cxcl12",
                          "Ccl2", "Cxcl1"),

  # aHSC_prof: 增殖中
  HSC_Proliferating = c("Mki67", "Mcm6", "Top2a", "Stmn1"),

  # ----- LSEC 4 亚型 + LVEC (Xiong 2019, Fig 3;沿用 v3.2 命名) -----
  # marker 来自 Xiong 2019 zonation 分析,但命名保持 LSEC_* 与 v3.2 一致
  LSEC_Periportal = c("Ednrb", "Jag1", "Lrg1", "Efnb1", "Ltbp4", "Adgrg6",
                      "Msr1", "Hes1"),

  LSEC_Pericentral = c("Wnt9b", "Rspo3", "Cdh13", "Wnt2", "Lyve1", "Kit",
                       "Thbd", "Fabp4"),

  LSEC = c("Stab2", "Stab1", "Fcgr2b", "Gpr182", "Clec4g", "Clec4f",
           "Mrc1", "Oit3"),

  LVEC = c("Vwf", "Pecam1", "Eng", "Cdh5", "Tek", "Kdr", "Flt1"),

  # ----- 其他 epithelial -----
  Cholangiocyte = c("Krt19", "Krt7", "Epcam", "Sox9", "Spp1", "Cftr",
                    "Hnf1b"),

  # ----- Macrophage 分类 (Xiong 2019, Bai 2023) ⭐ NAM 关键 -----
  # KC: 肝驻留巨噬细胞,F4/80+ Clec4f+
  Kupffer_KC = c("Adgre1", "Clec4f", "Vsig4", "Cd5l", "Marco", "Timd4",
                 "Csf1r", "C1qa", "C1qb"),

  # MDM: 单核细胞来源巨噬细胞,Cd11b+ Ccr2+
  Macrophage_MDM = c("Itgam", "Ccr2", "Ly6c2", "Ly6c1", "Cd14"),

  # NAM ⭐⭐⭐: NASH/肿瘤反应性巨噬细胞 (Xiong 2019 核心发现)
  # 关键 marker: Trem2hi, Gpnmb+, Cd9+
  Macrophage_NAM = c("Trem2", "Gpnmb", "Cd9", "Spp1", "Cd63", "Lpl",
                     "Cstb", "Fabp5", "Ms4a7", "Mmp12", "Ctsd"),

  # M2 polarization (Bai 2023)
  Macrophage_M2 = c("Mrc1", "Wnt5a", "Tgfb2", "Arg1", "Retnla", "Chil3",
                    "Cd163"),

  # M1 polarization (经典促炎极化)
  Macrophage_M1 = c("Cd86", "Cd80", "Nos2", "Tnf", "Il1b", "Il6",
                    "Cxcl9", "Cxcl10", "Fcgr1"),

  # ----- DC 亚型 -----
  DC_cDC1 = c("Xcr1", "Clec9a", "Itgae", "Irf8", "Batf3"),
  DC_cDC2 = c("Itgax", "Cd1c", "Sirpa", "Clec10a", "Irf4"),
  DC_pDC = c("Siglech", "Bst2", "Tcf4", "Ccr9", "Runx2"),

  # ----- T 亚型 -----
  T_CD4 = c("Cd3d", "Cd3e", "Cd4", "Cd40lg"),
  T_CD8 = c("Cd3d", "Cd3e", "Cd8a", "Cd8b1", "Gzmk", "Gzmb"),
  T_Treg = c("Foxp3", "Ctla4", "Il2ra", "Ikzf2", "Tnfrsf18"),
  T_NKT = c("Cd3d", "Klrb1c", "Nkg7", "Ncr1"),

  # ----- NK -----
  NK = c("Ncr1", "Nkg7", "Klrb1c", "Klre1", "Prf1", "Gzmb", "Xcl1"),

  # ----- B 亚型 (Bai 2023) -----
  B_Naive = c("Cd79a", "Cd79b", "Ms4a1", "Ighd", "Fcmr", "Ighm", "Cd19"),
  # B_Memory: 必须包含 Ms4a1/Cd19 这种 B 谱系 "必须表达" 的基因,
  # 否则只剩 Cd38/Cd27 这些通用活化 marker 会被任何增殖/炎症细胞偷分
  # (诊断: con_liver 第一轮里 B_Memory 收到 36 个跨谱系 flip, 根因是集合太短)
  B_Memory = c("Cd79a", "Cd79b", "Ms4a1", "Cd19", "Cd38", "Cd27"),
  B_Plasma = c("Jchain", "Igha", "Mzb1", "Xbp1", "Sdc1", "Prdm1"),

  # ----- 其他 -----
  Neutrophil = c("S100a8", "S100a9", "Ly6g", "Csf3r", "Lcn2", "Mpo"),
  Pericyte = c("Rgs5", "Mcam", "Pdgfrb", "Acta2", "Myh11"),
  Erythrocyte = c("Hba-a1", "Hba-a2", "Hbb-bs", "Hbb-bt", "Alas2"),
  Basophil = c("Cpa3", "Mcpt8", "Ms4a2", "Fcer1a")
)


# ==================== Human 肝 marker (PDX 人源肿瘤) ====================
# 主要来自 Xiong 2019 (人鼠对比) + 经典 HCC 文献

liver_markers_human_v4 <- list(

  # ----- HCC 主肿瘤群 -----
  # GPC3 是 HCC 经典 marker,在 PDX 模型里大部分细胞都是这个
  HCC_GPC3hi = c("GPC3", "AFP", "ALB", "APOA1", "APOA2", "TF", "FABP1",
                 "AHSG", "SERPINA1"),

  # ----- 肿瘤亚型 (基于 HCC 异质性文献) -----
  Tumor_Hepatocyte_like = c("ALB", "TF", "APOA1", "APOA2", "ASGR1", "TTR"),

  Tumor_EMT = c("VIM", "ZEB1", "ZEB2", "SNAI1", "SNAI2", "TWIST1", "CDH2",
                "FN1", "S100A4"),

  Tumor_Proliferating = c("MKI67", "TOP2A", "PCNA", "CCNB1", "CDK1", "CCNB2",
                          "CENPF", "BIRC5"),

  Tumor_Stem = c("PROM1", "EPCAM", "CD44", "ALDH1A1", "SOX9", "KRT19",
                 "NANOG", "LGR5"),

  # ----- 正常肝细胞 (拆分泄漏的少量细胞) -----
  Hep_Periportal_h = c("HAL", "SDS", "CYP2F2", "CDH1", "PCK1", "G6PC",
                       "CPS1", "ARG1"),

  Hep_Perivenous_h = c("GLUL", "CYP2E1", "CYP1A2", "OAT", "AXIN2"),

  # ----- 人源 TME (PDX 里几乎不存在,但有少量泄漏) -----
  # 注意: PDX 模型的免疫系统是鼠源,人源免疫细胞应该 ≈ 0
  HSC_Activated_h = c("ACTA2", "COL1A1", "COL1A2", "COL3A1", "FMOD", "LRAT"),

  Endothelial_h = c("PECAM1", "VWF", "CDH5", "CD34", "KDR"),

  LSEC_h = c("STAB2", "FCGR2B", "CLEC4G", "CLEC4M"),

  Cholangiocyte_h = c("KRT19", "KRT7", "EPCAM", "SOX9"),

  # 人源免疫细胞 (PDX 不应该有,如果检测到说明拆分泄漏)
  T_cell_h = c("CD3D", "CD3E", "CD3G", "CD8A", "CD4"),
  B_cell_h = c("CD79A", "CD79B", "MS4A1", "CD19"),
  NK_cell_h = c("NCAM1", "NKG7", "GNLY", "PRF1"),
  Macrophage_h = c("CD68", "CD163", "ITGAM", "AIF1", "C1QA"),
  Plasma_h = c("JCHAIN", "IGHA1", "IGHG1", "MZB1", "XBP1"),
  Neutrophil_h = c("S100A8", "S100A9", "FCGR3B", "CSF3R")
)


# ==================== 元数据:供主脚本使用 ====================

# 每个 marker 集的最少 marker 数 (低于这个数则跳过该亚型)
liver_marker_min_genes <- 3

# 文献来源标注 (写入 report)
liver_marker_sources <- c(
  "Xiong et al. 2019 Mol Cell (NASH liver scRNA-seq)",
  "Bai et al. 2023 Mil Med Res (Triclosan-induced liver fibrosis)",
  "Halpern et al. 2017 Nature (liver zonation)"
)

cat("[marker_lib_liver_v4] Loaded:\n")
cat("  - Mouse cell types:", length(liver_markers_mouse_v4), "\n")
cat("  - Human cell types:", length(liver_markers_human_v4), "\n")
cat("  - Sources:\n")
for (s in liver_marker_sources) cat("      *", s, "\n")
