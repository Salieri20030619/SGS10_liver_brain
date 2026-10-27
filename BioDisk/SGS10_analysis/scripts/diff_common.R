# ============================================================
# diff_common.R  --  阶段 1 差异分析的共用模块
# ============================================================
# 被 diff_v4_core.R 和 diff_v4_viz_enrich.R 共用
# 提供:
#   1. 样本 -> RDS 路径 resolver (支持 v4_mouse / v4_human / v3.2)
#   2. 细类名的三版本映射 (v4_refined / v4_broad / v3.2)
#   3. extract_cells() 函数,根据 annotation_version 和 cell_type 返回 cell barcodes
#   4. 对比组定义 (肝 6 个 + 脑 6 个,3 个 annotation 版本各 6 个)
#   5. 机制基因集 (9 个模块,肝脑通讯候选)
# ============================================================

# ------------------------------------------------------------
# 1. RDS 路径 resolver
# ------------------------------------------------------------
# v4 和 v3.2 RDS 都在 <sample>/mouse/ 子目录, 区别只是文件名后缀
#
# 肝:
#   v32: <sample>/mouse/<sample>_mouse_annotated.rds
#   v4:  <sample>/mouse/<sample>_mouse_annotated_v4.rds
#
# 脑:
#   v32: <sample>/mouse/<sample>_manual_annotated.rds
#   v4:  <sample>/mouse/<sample>_manual_annotated_v4.rds
#
# 注意: 顶层 <sample>/<sample>_manual_annotated.rds 是 3/25 旧版遗留,
#       不作为 v3.2 数据源
#
# PDX 肝样本 v4 分析只用 mouse RDS (鼠源 TME),
# human RDS 留给后续肿瘤异质性分析, 本阶段不用

PROJECT_ROOT <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"

# 样本分类
PDX_LIVER_SAMPLES  <- c("liver8", "gan56", "liver3")
PURE_LIVER_SAMPLES <- c("con_liver", "bla_liver")
BRAIN_SAMPLES      <- c("con_brain", "bla_brain", "nao56", "brain3", "brain8")
ALL_SAMPLES <- c(PURE_LIVER_SAMPLES, PDX_LIVER_SAMPLES, BRAIN_SAMPLES)

# 脑样本的 v4 文件名和肝不同 (manual_annotated 不是 mouse_annotated)
is_liver_sample <- function(sample) sample %in% c(PURE_LIVER_SAMPLES, PDX_LIVER_SAMPLES)
is_brain_sample <- function(sample) sample %in% BRAIN_SAMPLES

resolve_rds_path <- function(sample, version = c("v4", "v32")) {
  version <- match.arg(version)

  # v4 和 v32 都在 <sample>/mouse/ 子目录, 区别只是文件名后缀
  # v4:  <sample>_mouse_annotated_v4.rds (肝) / <sample>_manual_annotated_v4.rds (脑)
  # v32: <sample>_mouse_annotated.rds    (肝) / <sample>_manual_annotated.rds    (脑)
  #
  # 注意: 顶层的 <sample>/<sample>_manual_annotated.rds 是 3/25 旧数据, 不用

  if (is_liver_sample(sample)) {
    fname_base <- paste0(sample, "_mouse_annotated")
  } else if (is_brain_sample(sample)) {
    fname_base <- paste0(sample, "_manual_annotated")
  } else {
    stop("Unknown sample: ", sample)
  }

  fname <- if (version == "v4") paste0(fname_base, "_v4.rds") else paste0(fname_base, ".rds")
  file.path(PROJECT_ROOT, sample, "mouse", fname)
}

# 验证所有样本的 v4 和 v3.2 RDS 都存在
check_all_rds <- function(verbose = TRUE) {
  out <- data.frame(
    sample = character(0), version = character(0),
    path = character(0), exists = logical(0),
    stringsAsFactors = FALSE
  )
  for (s in ALL_SAMPLES) {
    for (v in c("v4", "v32")) {
      p <- resolve_rds_path(s, v)
      out <- rbind(out, data.frame(
        sample = s, version = v, path = p,
        exists = file.exists(p), stringsAsFactors = FALSE
      ))
    }
  }
  if (verbose) {
    n_miss <- sum(!out$exists)
    if (n_miss > 0) {
      cat("\n[WARN]", n_miss, "RDS files missing:\n")
      for (i in which(!out$exists)) {
        cat("  -", out$sample[i], "(", out$version[i], "):", out$path[i], "\n")
      }
    } else {
      cat("\n[OK] All", nrow(out), "RDS files exist\n")
    }
  }
  out
}

# ------------------------------------------------------------
# 2. 细类名的三版本映射 (virtual cell types)
# ------------------------------------------------------------
# 每个 virtual cell type 在不同 annotation 版本下对应哪些实际细类名
# 格式: list(virtual_name = list(v4_refined = c(...), v4_broad = c(...), v32 = c(...)))
#
# NULL 表示该版本下没有对应细类, 跑该 cell type 时会 skip

LIVER_CELL_TYPES <- list(
  "All_Cells" = list(
    v4_refined = "__ALL__",
    v4_broad   = "__ALL__",
    v32        = "__ALL__"
  ),
  "Hepatocyte_PP" = list(
    v4_refined = "Hep_Periportal",
    v4_broad   = "Hepatocyte_PP",
    v32        = "Hep_Periportal"
  ),
  "Hepatocyte_PV" = list(
    v4_refined = "Hep_Perivenous",
    v4_broad   = "Hepatocyte_PV",
    v32        = "Hep_Perivenous"
  ),
  "HSC" = list(
    v4_refined = c("HSC_Quiescent", "HSC_Activated_ECM", "HSC_Activated_Chemo",
                   "HSC_Activated_Remodel", "HSC_Proliferating"),
    v4_broad   = "HSC",
    v32        = c("HSC_Quiescent", "HSC_Activated_ECM", "HSC_Activated_Chemo",
                   "HSC_Proliferating", "Stellate")
  ),
  "HSC_Activated_ECM" = list(
    # v4 refined 和 v32 都有字面上的 HSC_Activated_ECM,
    # 但 v3.2 里对照组是过度分类的, 对比结果需要谨慎解读
    v4_refined = "HSC_Activated_ECM",
    v4_broad   = NULL,          # broad 合并了, 没法单独跑
    v32        = "HSC_Activated_ECM"
  ),
  "Kupffer" = list(
    v4_refined = "Kupffer_KC",
    v4_broad   = "Kupffer",
    v32        = c("KC_Homeostatic", "Kupffer")
  ),
  "Macrophage_NAM" = list(
    v4_refined = "Macrophage_NAM",
    v4_broad   = "Macrophage_NAM",
    v32        = "Macrophage_NAM"
  ),
  "DC" = list(
    v4_refined = c("DC_cDC1", "DC_cDC2", "DC_pDC"),
    v4_broad   = "DC",
    v32        = c("cDC1", "cDC2", "pDC", "Dendritic")
  )
)

BRAIN_CELL_TYPES <- list(
  "All_Cells" = list(
    v4_refined = "__ALL__",
    v4_broad   = "__ALL__",
    v32        = "__ALL__"
  ),
  "Endothelial" = list(
    v4_refined = "Endothelial",
    v4_broad   = "Endothelial",
    v32        = "Endothelial"
  ),
  "Microglia_Homeostatic" = list(
    v4_refined = "Microglia_Homeostatic",
    v4_broad   = NULL,          # broad 合并到 Microglia, 不分亚型
    v32        = "Microglia_Homeostatic"
  ),
  "Microglia_DAM" = list(
    v4_refined = "Microglia_DAM",
    v4_broad   = "Microglia_DAM",
    v32        = "Microglia_DAM"
  ),
  "Astrocyte_Gfap" = list(
    # v3.2 没有 Gfap 分类, 整个 Astrocytes 粗类也不对等
    # 只在 v4 refined 里跑
    v4_refined = "Astrocyte_Gfap",
    v4_broad   = NULL,
    v32        = NULL
  )
)

# annotation 版本 -> metadata 列名
VERSION_TO_COL <- list(
  v4_refined = "celltype_refined",
  v4_broad   = "celltype_broad",
  v32        = "final_celltype"
)

# ------------------------------------------------------------
# 3. extract_cells: 根据 virtual cell type + version 返回 cell barcodes
# ------------------------------------------------------------
extract_cells_v2 <- function(seu, virtual_ct, version,
                             cell_type_map = NULL) {
  if (is.null(cell_type_map)) {
    stop("cell_type_map required (LIVER_CELL_TYPES or BRAIN_CELL_TYPES)")
  }
  if (!virtual_ct %in% names(cell_type_map)) {
    stop("Unknown virtual cell type: ", virtual_ct)
  }
  mapping <- cell_type_map[[virtual_ct]][[version]]
  if (is.null(mapping)) return(character(0))  # 该版本不支持此细类

  meta_col <- VERSION_TO_COL[[version]]
  if (!meta_col %in% colnames(seu@meta.data)) {
    # fallback: 试其他列
    fallback_cols <- c("celltype_refined", "celltype_broad",
                       "final_celltype", "celltype")
    hit <- intersect(fallback_cols, colnames(seu@meta.data))
    if (length(hit) == 0) {
      stop("No annotation column found in this RDS")
    }
    meta_col <- hit[1]
    warning(sprintf("Column %s not found, using %s instead",
                    VERSION_TO_COL[[version]], meta_col))
  }

  if (identical(mapping, "__ALL__")) {
    return(colnames(seu))
  }
  vals <- as.character(seu@meta.data[[meta_col]])
  colnames(seu)[vals %in% mapping]
}

# ------------------------------------------------------------
# 4. 对比组定义
# ------------------------------------------------------------
# 每个 annotation 版本跑 6 个对比 (肝 3 + 脑 3)
# 同一个对比定义跨版本通用, 通过 extract_cells_v2 的 version 参数切换

LIVER_COMPARISONS <- list(
  list(comp = "liver_tumor_vs_blank",
       g1_samples = c("gan56", "liver3", "liver8"), g1_label = "Tumor",
       g2_samples = "bla_liver",                    g2_label = "Blank"),
  list(comp = "liver_tumor_vs_control",
       g1_samples = c("gan56", "liver3", "liver8"), g1_label = "Tumor",
       g2_samples = "con_liver",                    g2_label = "Control"),
  list(comp = "liver_blank_vs_control",
       g1_samples = "bla_liver",                    g1_label = "Blank",
       g2_samples = "con_liver",                    g2_label = "Control")
)

BRAIN_COMPARISONS <- list(
  list(comp = "brain_tumor_vs_blank",
       g1_samples = c("nao56", "brain3", "brain8"), g1_label = "Tumor",
       g2_samples = "bla_brain",                    g2_label = "Blank"),
  list(comp = "brain_tumor_vs_control",
       g1_samples = c("nao56", "brain3", "brain8"), g1_label = "Tumor",
       g2_samples = "con_brain",                    g2_label = "Control"),
  list(comp = "brain_blank_vs_control",
       g1_samples = "bla_brain",                    g1_label = "Blank",
       g2_samples = "con_brain",                    g2_label = "Control")
)

# ------------------------------------------------------------
# 5. 机制基因集 (9 个肝脑通讯候选)
# ------------------------------------------------------------
# 用户科学假设: 肝脏炎症因子入血 -> 脑小胶质 IL-6/JAK-STAT3 激活
# 原脚本只覆盖 IL6/TGFb/NFkB/PI3K/MAPK/Microglia_activ 5 个模块
# v5 扩展到 9 个模块, 让用户从多候选中选一条

MECHANISM_GENES <- list(
  # 原脚本 5 个模块
  IL6_JAK_STAT3 = c("Il6", "Il6ra", "Il6st", "Jak1", "Jak2", "Stat3",
                    "Socs3", "Stat1", "Pias1"),
  TGFb_SMAD     = c("Tgfb1", "Tgfb2", "Tgfbr1", "Tgfbr2",
                    "Smad2", "Smad3", "Smad4", "Smad7"),
  NFkB          = c("Nfkb1", "Nfkb2", "Rela", "Relb", "Tnf", "Il1b",
                    "Nfkbia", "Nfkbib"),
  PI3K_AKT      = c("Pik3ca", "Pik3cb", "Pik3cd", "Akt1", "Akt2",
                    "Mtor", "Pten"),
  MAPK          = c("Mapk1", "Mapk3", "Mapk8", "Mapk14",
                    "Map2k1", "Map2k2", "Fos", "Jun"),

  # v5 新增 4 个候选机制
  Complement_C3 = c("C3", "C1qa", "C1qb", "C1qc", "C3ar1", "C5ar1",
                    "Cfb", "Cfh", "Cfp"),
  CXCL_CXCR     = c("Cxcl10", "Cxcl9", "Cxcl1", "Cxcl2", "Cxcr3",
                    "Cxcr4", "Ccl2", "Ccl5", "Ccr2", "Ccr5"),
  S100_RAGE     = c("S100a8", "S100a9", "S100a4", "Ager", "Tlr4",
                    "Tlr2", "Myd88", "Hmgb1"),
  AcutePhase    = c("Saa1", "Saa2", "Saa3", "Lcn2", "Hp", "Cp",
                    "Crp", "Apcs", "Orm1", "Orm2"),

  # 原脚本的 Microglia_activ 模块保留, 改名更清晰
  Microglia_Activation = c("Cd68", "Cx3cr1", "Tmem119", "P2ry12", "Hexb",
                           "Trem2", "Apoe", "Cst7", "Cd9", "Spp1", "Itgax",
                           "Tspan18", "Basp1", "Il6st", "Cd74", "Smad3")
)

ALL_MECH_GENES <- unique(unlist(MECHANISM_GENES))

# ------------------------------------------------------------
# 6. 参数默认值 (可被主脚本覆盖)
# ------------------------------------------------------------
DEFAULT_PARAMS <- list(
  FC_THRESHOLD = 0.20,
  PADJ_CUT     = 0.10,
  MIN_PCT      = 0.10,
  MAX_CELLS    = 500,
  MIN_CELLS    = 5,
  RANDOM_SEED  = 42
)

# Il6st 专项追踪基因 (原脚本重点验证, v5 继续保留)
IL6_PATHWAY_CORE <- c("Il6", "Il6st", "Il6ra", "Jak1", "Jak2",
                      "Stat3", "Socs3")

cat("[diff_common] Loaded:\n")
cat("  - samples:", length(ALL_SAMPLES),
    "(liver:", length(c(PDX_LIVER_SAMPLES, PURE_LIVER_SAMPLES)),
    ", brain:", length(BRAIN_SAMPLES), ")\n")
cat("  - liver virtual cell types:", length(LIVER_CELL_TYPES), "\n")
cat("  - brain virtual cell types:", length(BRAIN_CELL_TYPES), "\n")
cat("  - annotation versions: v4_refined, v4_broad, v32\n")
cat("  - comparisons: liver 3, brain 3 (x 3 versions = 18 comparison-version pairs)\n")
cat("  - mechanism modules:", length(MECHANISM_GENES),
    "(", length(ALL_MECH_GENES), "unique genes)\n")
