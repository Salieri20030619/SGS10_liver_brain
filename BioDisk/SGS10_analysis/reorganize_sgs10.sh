#!/bin/bash
# SGS10 Deliverables Reorganize Script
# 生成时间: 2026-04-20
# 
# 功能:
#   1. 创建新的 SGS10_deliverables_2026-04-20/ 目录结构
#   2. 从原位置复制文件 (不移动, 保留原位置不变)
#   3. 排除所有 *.rds, *.h5, *.h5ad
#   4. 重组成清晰的 00-08 编号目录
#
# 运行前请先 cd ~/SGS10_analysis/
# 预计用时: 5-10 分钟
# 预计新目录大小: ~300-400 MB

set -e
set -u

OUTPUT=~/SGS10_analysis/SGS10_deliverables_2026-04-20
SRC=~/SGS10_analysis

# 用于选择性复制 (排除 rds/h5/h5ad)
RSYNC_BASE='rsync -a --exclude=*.rds --exclude=*.h5 --exclude=*.h5ad'

echo "======================================================"
echo "SGS10 Deliverables Reorganize"
echo "输出目录: $OUTPUT"
echo "======================================================"
echo ""

# === 创建目录结构 ===
echo "[1/9] 创建目录结构..."
mkdir -p "$OUTPUT"/{01_reports,\
02_sample_annotation/01_first_annotation_apr11_12,\
02_sample_annotation/02_dam_secondary_annotation_apr16,\
02_sample_annotation/03_celltype_export/per_sample_celltype,\
03_sc_DEG/liver,\
03_sc_DEG/brain_clean,\
03_sc_DEG/per_mouse/liver,\
03_sc_DEG/per_mouse/brain,\
03_sc_DEG/viz,\
03_sc_DEG/summary,\
04_enrichment,\
05_TF_activity/brain_clean,\
05_TF_activity/liver,\
06_pseudotime/liver,\
06_pseudotime/brain,\
07_Zbtb7a_analysis,\
08_spatial_Task_B/T1_niche_10sample,\
08_spatial_Task_B/T4_spatial_domain,\
08_spatial_Task_B/T5_liver_brain_paired_axis,\
08_spatial_Task_B/T6_cellchat/lr_pairs,\
08_spatial_Task_B/T6_cellchat/pathway_pairs,\
08_spatial_Task_B/T6_cellchat/pathway_plots,\
08_spatial_Task_B/spatial_figures}

# === 01_reports ===
echo "[2/9] 复制 01_reports..."
cp -v "$SRC/task_A_v1"/*.md "$SRC/task_A_v1"/*.docx "$OUTPUT/01_reports/" 2>/dev/null || true
# v1.3 定论声明 + cross-validation 如果服务器也有的话
find "$SRC" -maxdepth 3 -type f \( -name "SGS10_TaskA_v1_3*" -o -name "SGS10_TaskA_vs_TaskB*" -o -name "SGS10_TaskB_v1_Full_Report_v1_5*" \) \
  -exec cp -v {} "$OUTPUT/01_reports/" \; 2>/dev/null || true

# === 02_sample_annotation/01_first_annotation ===
echo "[3/9] 复制 02_sample_annotation/01_first_annotation..."
for sample in brain3 brain8 nao56 bla_brain con_brain bla_liver con_liver liver3 liver8; do
    if [ -d "$SRC/$sample" ]; then
        mkdir -p "$OUTPUT/02_sample_annotation/01_first_annotation_apr11_12/$sample"
        $RSYNC_BASE "$SRC/$sample/" "$OUTPUT/02_sample_annotation/01_first_annotation_apr11_12/$sample/"
    fi
done
# gan56 保留 human/mouse 子结构
if [ -d "$SRC/gan56" ]; then
    mkdir -p "$OUTPUT/02_sample_annotation/01_first_annotation_apr11_12/gan56"
    $RSYNC_BASE "$SRC/gan56/" "$OUTPUT/02_sample_annotation/01_first_annotation_apr11_12/gan56/"
fi

# === 02_sample_annotation/02_dam_secondary ===
echo "[4/9] 复制 02_sample_annotation/02_dam_secondary..."
$RSYNC_BASE "$SRC/dam_clean_v2/" "$OUTPUT/02_sample_annotation/02_dam_secondary_annotation_apr16/"

# === 02_sample_annotation/03_celltype_export ===
echo "[5/9] 复制 02_sample_annotation/03_celltype_export..."
# per-sample celltype CSV 放到子目录
mv_src="$SRC/sc_celltype_export_v3"
for f in "$mv_src"/*_celltype.csv; do
    if [ -f "$f" ] && [ "$(basename "$f")" != "all_samples_celltype.csv" ]; then
        cp "$f" "$OUTPUT/02_sample_annotation/03_celltype_export/per_sample_celltype/"
    fi
done
cp "$mv_src/all_samples_celltype.csv" "$OUTPUT/02_sample_annotation/03_celltype_export/" 2>/dev/null || true
cp "$mv_src/NOTE.md" "$OUTPUT/02_sample_annotation/03_celltype_export/" 2>/dev/null || true
# top30_markers 整个保留
$RSYNC_BASE "$mv_src/top30_markers/" "$OUTPUT/02_sample_annotation/03_celltype_export/top30_markers/"

# === 03_sc_DEG ===
echo "[6/9] 复制 03_sc_DEG..."
# Liver DEG (flat)
$RSYNC_BASE "$SRC/diff_stage1_results/deg_v4_broad/liver/" "$OUTPUT/03_sc_DEG/liver/"
# Brain DEG - 保留对照组子目录结构
$RSYNC_BASE "$SRC/diff_stage1_brain_clean/deg_v4_refined/" "$OUTPUT/03_sc_DEG/brain_clean/"
# per-mouse DEG
$RSYNC_BASE "$SRC/diff_stage3_results/per_sample_deg/liver/" "$OUTPUT/03_sc_DEG/per_mouse/liver/"
$RSYNC_BASE "$SRC/diff_stage3_results/per_sample_deg/brain/" "$OUTPUT/03_sc_DEG/per_mouse/brain/"
# viz
$RSYNC_BASE "$SRC/diff_stage1_results/viz/" "$OUTPUT/03_sc_DEG/viz/"
# summary
$RSYNC_BASE "$SRC/diff_stage3_results/summary/" "$OUTPUT/03_sc_DEG/summary/summary_stage3/"
$RSYNC_BASE "$SRC/diff_stage3_results/response_ranking/" "$OUTPUT/03_sc_DEG/summary/response_ranking/"
$RSYNC_BASE "$SRC/diff_stage3_results/fc_correlation/" "$OUTPUT/03_sc_DEG/summary/fc_correlation/"
$RSYNC_BASE "$SRC/diff_stage3_results/overlap/" "$OUTPUT/03_sc_DEG/summary/overlap/"
cp "$SRC/diff_stage1_results/summary/deg_summary.csv" "$OUTPUT/03_sc_DEG/summary/deg_summary_stage1.csv" 2>/dev/null || true
cp "$SRC/diff_stage1_brain_clean/summary/deg_summary.csv" "$OUTPUT/03_sc_DEG/summary/deg_summary_brain_clean.csv" 2>/dev/null || true

# === 04_enrichment ===
echo "[7/9] 复制 04_enrichment..."
$RSYNC_BASE "$SRC/diff_stage1_brain_clean/enrichment/" "$OUTPUT/04_enrichment/brain_clean/"

# === 05_TF_activity ===
$RSYNC_BASE "$SRC/dorothea_results_brain_clean/" "$OUTPUT/05_TF_activity/brain_clean/"
# liver TF: 从 dorothea_results 里过滤出 liver*
if [ -d "$SRC/dorothea_results/per_contrast" ]; then
    mkdir -p "$OUTPUT/05_TF_activity/liver/per_contrast"
    for f in "$SRC/dorothea_results/per_contrast"/liver*.csv; do
        if [ -f "$f" ]; then cp "$f" "$OUTPUT/05_TF_activity/liver/per_contrast/"; fi
    done
    # 其他非 per_contrast 文件
    find "$SRC/dorothea_results" -maxdepth 1 -type f ! -name "*.rds" -exec cp {} "$OUTPUT/05_TF_activity/liver/" \; 2>/dev/null || true
fi

# === 06_pseudotime ===
echo "[8/9] 复制 06_pseudotime..."
$RSYNC_BASE "$SRC/pseudotime_v4_results/" "$OUTPUT/06_pseudotime/liver/"
$RSYNC_BASE "$SRC/pseudotime_brain_slingshot/" "$OUTPUT/06_pseudotime/brain/"

# === 07_Zbtb7a ===
$RSYNC_BASE "$SRC/taskA_v1_2/zbtb7a_analysis/" "$OUTPUT/07_Zbtb7a_analysis/"

# === 08_spatial_Task_B ===
echo "[9/9] 复制 08_spatial_Task_B..."
$RSYNC_BASE "$SRC/taskB/niche_v2_10samples/" "$OUTPUT/08_spatial_Task_B/T1_niche_10sample/"
$RSYNC_BASE "$SRC/taskB/spatial_domain/" "$OUTPUT/08_spatial_Task_B/T4_spatial_domain/"
$RSYNC_BASE "$SRC/taskB/liver_brain_axis2/" "$OUTPUT/08_spatial_Task_B/T5_liver_brain_paired_axis/"

# CellChat 拆分: lr_pairs, pathway_pairs, pathway_plots
cp "$SRC/taskB/cellchat"/*_lr_pairs.csv "$OUTPUT/08_spatial_Task_B/T6_cellchat/lr_pairs/" 2>/dev/null || true
cp "$SRC/taskB/cellchat"/*_pathway_pairs.csv "$OUTPUT/08_spatial_Task_B/T6_cellchat/pathway_pairs/" 2>/dev/null || true
$RSYNC_BASE "$SRC/taskB/cellchat/pathway_plots/" "$OUTPUT/08_spatial_Task_B/T6_cellchat/pathway_plots/"
# 可能的 aggregate network PDF
find "$SRC/taskB/cellchat" -maxdepth 1 -type f -name "*.pdf" -exec cp {} "$OUTPUT/08_spatial_Task_B/T6_cellchat/" \; 2>/dev/null || true

# spatial_figures
$RSYNC_BASE "$SRC/taskB/spatial_figures/" "$OUTPUT/08_spatial_Task_B/spatial_figures/"

# === 最终统计 ===
echo ""
echo "======================================================"
echo "Reorganize 完成!"
echo "======================================================"
echo "新目录: $OUTPUT"
echo ""
echo "各模块大小:"
du -sh "$OUTPUT"/* 2>/dev/null | sort -k1 -hr
echo ""
echo "总大小:"
du -sh "$OUTPUT"
echo ""
echo "文件总数:"
find "$OUTPUT" -type f | wc -l
echo ""
echo "=== 验证: 无 rds/h5/h5ad ==="
find "$OUTPUT" -type f \( -name "*.rds" -o -name "*.h5" -o -name "*.h5ad" \) | head -5
if [ -z "$(find "$OUTPUT" -type f \( -name "*.rds" -o -name "*.h5" -o -name "*.h5ad" \))" ]; then
    echo "✓ 无 RDS/h5 泄漏"
else
    echo "⚠️ 上面文件需要手动清理!"
fi
