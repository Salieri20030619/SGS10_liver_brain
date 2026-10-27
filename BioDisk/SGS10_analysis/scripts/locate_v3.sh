#!/bin/bash
# locate_v3.sh - 按每张图的核心基因/关键词, 在 4/10 之后的 CSV 里精确定位数据
# 输出每个 hit 的文件路径 + 该文件的前 3 行 header, 让我能确认列名
# 用法: bash locate_v3.sh > /tmp/fig_v3.txt 2>&1
ROOT=/mnt/BioDisk/wangmingyuan
cd $ROOT
CUTOFF="2026-04-10"

# 帮助函数: 在 CSV 里搜核心 token, 列出命中文件 + header
# args: $1=tag $2=token1 [$3=token2 ...]
hunt() {
  local tag=$1; shift
  echo ""
  echo "############################################################"
  echo "## [$tag]  tokens: $*"
  echo "############################################################"
  for token in "$@"; do
    # -l 只列文件; --include 限定后缀; 限定时间用 newermt 后过滤
    local files=$(find . -name "*.csv" -newermt "$CUTOFF" -not -path "*/miniconda3/*" \
      -not -path "*/SGS10_backup/*" 2>/dev/null \
      | xargs grep -l -- "$token" 2>/dev/null | head -8)
    if [ -n "$files" ]; then
      echo "-- token='$token' hits --"
      echo "$files"
    fi
  done
}

# 每张图的关键 token (基因/通路/分类名都来自 README + 你已确认的真实数值)

# === Fig HSC volcano ===
hunt "Fig_HSC_volcano"  "HSC_Activated_ECM"  "Serpine1"  "Acta2"

# === Fig HSC DEG lollipop ===
hunt "Fig_HSC_DEG_lollipop"  "Col1a1"  "Tagln"  "Pmepa1"  "Loxl1"

# === Fig HSC_myeloid dumbbell (Gas6/C3 ligand, Cd74/Csf1r receptor) ===
hunt "Fig_HSC_myeloid_dumbbell"  "Gas6"  "C3"  "Cd74"  "Csf1r"  "Axl"

# === Fig silent_inflammation (Ifng Il1b Il6 Tnf) ===
hunt "Fig_silent_inflammation"  "Ifng"  "Il1b"  "cytokine"

# === Fig Excit_L45 TF lollipop (Tead4/Pou2f2/Myc/Cebpa/Nr2c2/Tfap2a) ===
hunt "Fig_Excit_L45_TF"  "Excitatory_L45"  "Tead4"  "Tfap2a"  "Pou2f2"

# === Fig TF_dual_mode ===
hunt "Fig_TF_dual_mode"  "TF_activity"  "decoupler"

# === Fig thyroid cascade (Ttr / Thra / Hep_PP / CPEC / Endo) ===
hunt "Fig_thyroid"  "Ttr"  "Thra"  "Hepatocyte_PP"  "CPECs"

# === Fig efferocytosis (Gas6/Axl/Mertk + Pecam1/Bsg/Calr/Mfge8 Brain Pericyte) ===
hunt "Fig_efferocytosis"  "Mertk"  "Mfge8"  "Pericytes"  "Pecam1"

# === Fig X1 TGFb (Smad3/Tgfb1/HSC zone) ===
hunt "Fig_X1_TGFb"  "Tgfb1"  "Smad3"  "Tgfbr"

# === Fig X2 brain spatial ===
hunt "Fig_X2_brain_spatial"  "X2_brain"

# === Fig microenv_heatmap (Liver M01..M07 / Brain M01..M06) ===
hunt "Fig_microenv_heatmap"  "M01_"  "M05_Tumor"  "M04_Vascular_BBB"  "microenv_name"

# === Fig markers_heatmap (Liver/Brain microenv top markers) ===
hunt "Fig_markers_heatmap"  "microenv_markers_top"  "avg_log2FC"

# === Fig hex_bin (Liver_3 / Brain_3 spot density) ===
hunt "Fig_hex_bin"  "Liver_3"  "Brain_3"  "spot_id"  "x_pixel"

# === Fig spatial_overlay / HSC_5x5 ===
hunt "Fig_spatial_overlay"  "rctd"  "spot_x"  "spot_y"

# === Fig CellChat 3groups ===
hunt "Fig_CellChat"  "ligand_receptor"  "interaction_name"  "source.target"

# === FigS2 species_separation ===
hunt "FigS2_species"  "mouse_pct"  "human_pct"  "species_classification"

# === FigS3 DAM_cleaning ===
hunt "FigS3_DAM"  "dam_decision"  "DAM_strict"

# === Brain markers heatmap (FigS4) ===
hunt "FigS4_Brain_markers"  "Hippocampus"  "Excitatory"  "OPCs"


echo ""
echo "############################################################"
echo "## STAGE G: 列出每个核心 CSV 的列名 (前 1 行 header)"
echo "############################################################"
# 关键 CSV 全部 dump header
CSV_LIST=(
  "SGS10_analysis/diff_stage1_results/deg_v4_refined/liver/liver_tumor_vs_control_HSC_Activated_ECM.csv"
  "SGS10_analysis/D1_v2_RCTD_composition/microenv_markers_RCTD100/liver_microenv_markers_top20.csv"
  "SGS10_analysis/D1_v2_RCTD_composition/microenv_markers_RCTD100/brain_microenv_markers_top20.csv"
  "SGS10_analysis/D1_v2_RCTD_composition/microenv_markers_RCTD100/liver_M05_tumor_specific_markers.csv"
  "SGS10_analysis/D1_v2_RCTD_composition/microenv_markers_RCTD100/brain_microenv_in_tumor_markers.csv"
  "SGS10_analysis/D1_v2_RCTD_composition/microenv_key_celltypes_RCTD100/liver_microenv_celltypes_full.csv"
  "SGS10_analysis/D1_v2_RCTD_composition/microenv_key_celltypes_RCTD100/brain_microenv_celltypes_full.csv"
  "SGS10_analysis/D1_v2_RCTD_composition/microenv_key_celltypes_RCTD100/liver_microenv_specificity.csv"
  "SGS10_analysis/D1_v2_RCTD_composition/microenv_key_celltypes_RCTD100/brain_microenv_specificity.csv"
  "SGS10_analysis/taskB/spatial_domain_RCTD100/liver_microenv_naming.csv"
  "SGS10_analysis/taskB/spatial_domain_RCTD100/brain_microenv_naming.csv"
  "SGS10_analysis/taskB/spatial_domain_RCTD100/liver_domain_ct_profile.csv"
  "SGS10_analysis/taskB/spatial_domain_RCTD100/brain_domain_ct_profile.csv"
  "SGS10_analysis/taskB/cellchat/3groups_plots_v3/liver_3groups_top_pairs.csv"
  "SGS10_analysis/taskB/cellchat/3groups_plots_v3/brain_3groups_top_pairs.csv"
  "SGS10_analysis/dorothea_results_brain_clean/per_contrast/brain__tumor_vs_blank__Excitatory_L45_tf_activity.csv"
  "SGS10_analysis/dorothea_results_brain_clean/per_contrast/brain__tumor_vs_control__Excitatory_L45_tf_activity.csv"
  "SGS10_analysis/dam_clean_v2/all_dam_decisions.csv"
  "SGS10_analysis/species_assignment/liver3_species_assignment.csv"
  "SGS10_analysis/species_assignment/gan56_species_assignment.csv"
  "SGS10_analysis/D3_DEG_output/pseudobulk_edgeR/A1_tumor_core_vs_baseline_host.csv"
  "SGS10_analysis/D3_DEG_output/pseudobulk_edgeR/A2_invaded_host_vs_baseline_host.csv"
  "SGS10_analysis/D3_DEG_output/pseudobulk_edgeR/A3_distant_host_vs_baseline_host.csv"
  "SGS10_analysis/D3_DEG_output/pseudobulk_edgeR/C4_brain_distant_vs_brain_baseline.csv"
  "SGS10_analysis/SGS10_TaskC_deliverables/SGS10_spot_coordinates_full.csv"
)
for f in "${CSV_LIST[@]}"; do
  if [ -f "$f" ]; then
    nrows=$(wc -l < "$f")
    echo ""
    echo "### $f  (rows=$nrows) ###"
    head -1 "$f"
    echo "---first 2 data rows---"
    sed -n '2,3p' "$f"
  else
    echo ""
    echo "### MISSING: $f ###"
  fi
done

echo ""
echo "############################################################"
echo "## STAGE H: 关键 RDS (Seurat / RCTD / CellChat) 大小"
echo "############################################################"
find SGS10_analysis -name "*.rds" -newermt "$CUTOFF" 2>/dev/null \
  | xargs ls -lh 2>/dev/null | awk '{print $5, $9}' | sort -rh | head -30

echo ""
echo "DONE."
