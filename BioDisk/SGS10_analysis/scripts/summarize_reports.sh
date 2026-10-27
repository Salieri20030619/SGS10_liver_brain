#!/bin/bash
# ==================== SGS10 汇总脚本 ====================
# 文件: ~/SGS10_analysis/scripts/summarize_reports.sh
# 用途: 扫描所有样本目录的 report.txt,提取关键指标,输出 markdown + csv
# 用法: bash summarize_reports.sh
#       输出到 SGS10_analysis/SGS10_summary_YYYYMMDD.md 和 .csv
# ========================================================

set -u
ROOT="/mnt/BioDisk/wangmingyuan/SGS10_analysis"
DATE=$(date +%Y%m%d_%H%M)
MD_OUT="$ROOT/SGS10_summary_${DATE}.md"
CSV_OUT="$ROOT/SGS10_summary_${DATE}.csv"

# 样本清单
LIVER_PURE=(con_liver bla_liver)                  # 纯鼠肝
LIVER_PDX=(liver8 gan56 liver3)                   # PDX 肝(有 mouse + human 子目录)
BRAIN=(con_brain bla_brain nao56 brain3 brain8)   # 纯鼠脑

# ---------- 解析单个 report.txt 的核心函数 ----------
# 用法: parse_report <path_to_report.txt>
# 输出: 一行,空格分隔的 KEY=VALUE 字段
parse_report() {
    local f="$1"
    if [ ! -f "$f" ]; then
        echo "MISSING"
        return
    fi
    awk '
    BEGIN { cells=0; genes=0; in_dist=0 }
    /^Cells:/ {
        # 格式: Cells: 32188   Genes: 32201
        for (i=1; i<=NF; i++) {
            if ($i == "Cells:") cells=$(i+1)
            if ($i == "Genes:") genes=$(i+1)
        }
    }
    /^Cell type distribution/ { in_dist=1; next }
    in_dist && /^  / {
        # 格式: "  CellType    N ( P%)"
        # 拿第一个非空字段做 celltype, 倒数第二个去掉数字做 count
        ct=$1
        n=$2
        # 百分比在括号里, 形如 "(42.4%)" 或 "( 0.2%)"
        pct=$NF
        gsub(/[()%]/, "", pct)
        counts[ct]=n
        pcts[ct]=pct
    }
    END {
        printf "CELLS=%s GENES=%s", cells, genes
        for (ct in counts) {
            printf " CT_%s=%s PCT_%s=%s", ct, counts[ct], ct, pcts[ct]
        }
        printf "\n"
    }
    ' "$f"
}

# ---------- 从 parse_report 输出里提取某个细胞类型的百分比 ----------
get_pct() {
    local parsed="$1"
    local ct="$2"
    echo "$parsed" | tr ' ' '\n' | awk -F= -v key="PCT_$ct" '$1==key {print $2; exit}'
}

get_count() {
    local parsed="$1"
    local ct="$2"
    echo "$parsed" | tr ' ' '\n' | awk -F= -v key="CT_$ct" '$1==key {print $2; exit}'
}

get_field() {
    local parsed="$1"
    local key="$2"
    echo "$parsed" | tr ' ' '\n' | awk -F= -v key="$key" '$1==key {print $2; exit}'
}

# 求和多个细胞类型的百分比(用于 T 细胞总和、Microglia 总和等)
sum_pcts() {
    local parsed="$1"
    shift
    local total=0
    for ct in "$@"; do
        local v=$(get_pct "$parsed" "$ct")
        if [ -n "$v" ]; then
            total=$(awk -v a="$total" -v b="$v" 'BEGIN{printf "%.2f", a+b}')
        fi
    done
    echo "$total"
}

# 比值
ratio() {
    local a="$1"
    local b="$2"
    if [ -z "$a" ] || [ -z "$b" ] || [ "$b" = "0" ]; then
        echo "-"
    else
        awk -v x="$a" -v y="$b" 'BEGIN{printf "%.2f", x/y}'
    fi
}

# 缺失值兜底
nz() { echo "${1:-0}"; }

# ==================== 写 Markdown ====================
{
echo "# SGS10 项目样本汇总"
echo ""
echo "**生成时间:** $(date '+%Y-%m-%d %H:%M:%S')"
echo "**Marker library:** liver = v3.2 (AddModuleScore, Halpern2017+Xiong2019+Bai2023);brain = current"
echo ""
echo "---"
echo ""

# ---------- 表 1: 肝脏样本(鼠源) ----------
echo "## 表 1:肝脏样本 — 鼠源细胞(TME / 正常肝)"
echo ""
echo "| 样本 | 类型 | Cells | PP% | PV% | PP:PV | Hep_Prolif% | NAM% | HSC_ECM% | T总和% | Unknown% |"
echo "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|"

# 纯鼠肝
for s in "${LIVER_PURE[@]}"; do
    f="$ROOT/$s/${s}_mouse_report.txt"
    [ ! -f "$f" ] && f="$ROOT/$s/mouse/${s}_mouse_report.txt"
    [ ! -f "$f" ] && f=$(find "$ROOT/$s" -name "*report.txt" 2>/dev/null | head -1)
    parsed=$(parse_report "$f")
    cells=$(get_field "$parsed" "CELLS")
    pp=$(get_pct "$parsed" "Hep_Periportal")
    pv=$(get_pct "$parsed" "Hep_Perivenous")
    rt=$(ratio "$(nz $pp)" "$(nz $pv)")
    hp=$(get_pct "$parsed" "Hep_Proliferating")
    nam=$(get_pct "$parsed" "Macrophage_NAM")
    hsc=$(get_pct "$parsed" "HSC_Activated_ECM")
    t_sum=$(sum_pcts "$parsed" CD4_T CD8_T Treg T_cells NKT)
    unk=$(get_pct "$parsed" "Unknown")
    echo "| $s | 纯鼠 | $(nz $cells) | $(nz $pp) | $(nz $pv) | $rt | $(nz $hp) | $(nz $nam) | $(nz $hsc) | $t_sum | $(nz $unk) |"
done

# PDX 肝鼠源
for s in "${LIVER_PDX[@]}"; do
    f="$ROOT/$s/mouse/${s}_mouse_report.txt"
    parsed=$(parse_report "$f")
    cells=$(get_field "$parsed" "CELLS")
    pp=$(get_pct "$parsed" "Hep_Periportal")
    pv=$(get_pct "$parsed" "Hep_Perivenous")
    rt=$(ratio "$(nz $pp)" "$(nz $pv)")
    hp=$(get_pct "$parsed" "Hep_Proliferating")
    nam=$(get_pct "$parsed" "Macrophage_NAM")
    hsc=$(get_pct "$parsed" "HSC_Activated_ECM")
    t_sum=$(sum_pcts "$parsed" CD4_T CD8_T Treg T_cells NKT)
    unk=$(get_pct "$parsed" "Unknown")
    echo "| $s | PDX 鼠源 | $(nz $cells) | $(nz $pp) | $(nz $pv) | $rt | $(nz $hp) | $(nz $nam) | $(nz $hsc) | $t_sum | $(nz $unk) |"
done

echo ""
echo "**字段说明**:"
echo "- PP/PV = Hep_Periportal / Hep_Perivenous"
echo "- PP:PV 比值理想 ~1:1,正常肝 con_liver 约 2.7:1"
echo "- NAM (Trem2+ Macrophage) 在肿瘤反应中富集"
echo "- HSC_ECM (Activated stellate, ECM-producing) 是纤维化标志"
echo "- T总和 = CD4_T + CD8_T + Treg + T_cells + NKT"
echo ""
echo "---"
echo ""

# ---------- 表 2: PDX 人源(肿瘤) ----------
echo "## 表 2:PDX 肝脏样本 — 人源细胞(肿瘤)"
echo ""
echo "| 样本 | Cells | HCC_GPC3hi% | Tumor亚型总和% | 非肿瘤泄漏% |"
echo "|---|---:|---:|---:|---:|"

for s in "${LIVER_PDX[@]}"; do
    f="$ROOT/$s/human/${s}_human_report.txt"
    parsed=$(parse_report "$f")
    cells=$(get_field "$parsed" "CELLS")
    hcc=$(get_pct "$parsed" "HCC_GPC3hi")
    tumor_sub=$(sum_pcts "$parsed" Tumor_Hepatocyte Tumor_EMT Tumor_Proliferating Tumor_Stem)
    leak=$(sum_pcts "$parsed" Hep_Periportal Hep_Perivenous LSEC Endothelial Stellate HSC_Activated Plasma B_cells CD4_T CD8_T NK_cells T_cells Monocytes Neutrophils Kupffer_Macrophage Macrophage_NAM Cholangiocytes Treg Dendritic)
    echo "| $s | $(nz $cells) | $(nz $hcc) | $tumor_sub | $leak |"
done

echo ""
echo "**字段说明**:"
echo "- HCC_GPC3hi 为主要肿瘤群(GPC3+ 肝癌细胞)"
echo "- Tumor亚型总和 = Tumor_Hepatocyte + Tumor_EMT + Tumor_Proliferating + Tumor_Stem"
echo "- 非肿瘤泄漏 = 人鼠拆分时混入的非肿瘤细胞,理想 < 2%"
echo ""
echo "---"
echo ""

# ---------- 表 3: 脑样本 ----------
echo "## 表 3:脑样本(纯鼠)"
echo ""
echo "| 样本 | 类型 | Cells | Hpc_CA1% | Excit总和% | Inhib总和% | Astro% | Oligo% | OPC% | Micro总和% | DAM% | Endo% | Peri% |"
echo "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"

for s in "${BRAIN[@]}"; do
    f="$ROOT/$s/${s}_report.txt"
    [ ! -f "$f" ] && f="$ROOT/$s/mouse/${s}_report.txt"
    [ ! -f "$f" ] && f=$(find "$ROOT/$s" -name "*report.txt" 2>/dev/null | head -1)

    # 类型标签
    case "$s" in
        con_brain) label="正常" ;;
        bla_brain) label="空白" ;;
        *) label="PDX" ;;
    esac

    parsed=$(parse_report "$f")
    cells=$(get_field "$parsed" "CELLS")
    ca1=$(get_pct "$parsed" "Hpc_CA1")
    ex=$(sum_pcts "$parsed" Excitatory_L23 Excitatory_L45 Excitatory_L6)
    inh=$(sum_pcts "$parsed" Inhibitory_Pvalb Inhibitory_Meis2 Inhibitory_Sst Inhibitory_Vip)
    astro=$(get_pct "$parsed" "Astrocytes")
    oligo=$(get_pct "$parsed" "Oligodendrocytes")
    opc=$(get_pct "$parsed" "OPCs")
    micro=$(sum_pcts "$parsed" Microglia Microglia_Homeostatic Microglia_DAM Microglia_M1 Microglia_M2)
    dam=$(get_pct "$parsed" "Microglia_DAM")
    endo=$(get_pct "$parsed" "Endothelial")
    peri=$(get_pct "$parsed" "Pericytes")
    echo "| $s | $label | $(nz $cells) | $(nz $ca1) | $ex | $inh | $(nz $astro) | $(nz $oligo) | $(nz $opc) | $micro | $(nz $dam) | $(nz $endo) | $(nz $peri) |"
done

echo ""
echo "**字段说明**:"
echo "- Excit总和 = Excitatory_L23 + L45 + L6"
echo "- Inhib总和 = Inhibitory_Pvalb + Meis2 + Sst + Vip"
echo "- Micro总和 = Microglia + Homeostatic + DAM + M1 + M2"
echo "- DAM = Disease-Associated Microglia(肿瘤反应小胶质)"
echo ""
echo "---"
echo ""

# ---------- 检查清单 ----------
echo "## 文件完整性检查"
echo ""
echo "| 样本 | report.txt | annotated.rds |"
echo "|---|---|---|"

check_file() { [ -f "$1" ] && echo "✓" || echo "✗"; }

for s in "${LIVER_PURE[@]}"; do
    rpt=$(find "$ROOT/$s" -name "*report.txt" 2>/dev/null | head -1)
    rds=$(find "$ROOT/$s" -name "*annotated.rds" 2>/dev/null | head -1)
    echo "| $s | $([ -n "$rpt" ] && echo "✓" || echo "✗") | $([ -n "$rds" ] && echo "✓" || echo "✗") |"
done
for s in "${LIVER_PDX[@]}"; do
    rpt_m="$ROOT/$s/mouse/${s}_mouse_report.txt"
    rpt_h="$ROOT/$s/human/${s}_human_report.txt"
    rds_m="$ROOT/$s/mouse/${s}_mouse_annotated.rds"
    rds_h="$ROOT/$s/human/${s}_human_annotated.rds"
    echo "| ${s} (mouse) | $(check_file "$rpt_m") | $(check_file "$rds_m") |"
    echo "| ${s} (human) | $(check_file "$rpt_h") | $(check_file "$rds_h") |"
done
for s in "${BRAIN[@]}"; do
    rpt=$(find "$ROOT/$s" -name "*report.txt" 2>/dev/null | head -1)
    rds=$(find "$ROOT/$s" -name "*manual_annotated.rds" -o -name "*annotated.rds" 2>/dev/null | head -1)
    echo "| $s | $([ -n "$rpt" ] && echo "✓" || echo "✗") | $([ -n "$rds" ] && echo "✓" || echo "✗") |"
done

echo ""
echo "---"
echo ""
echo "*Report 由 \`scripts/summarize_reports.sh\` 自动生成。手动二次注释和组间差异分析请参考后续脚本。*"
} > "$MD_OUT"

# ==================== 写 CSV ====================
{
# 表头(union 所有可能的列)
echo "sample,group,species,cells,Hep_Periportal_pct,Hep_Perivenous_pct,PP_PV_ratio,Hep_Proliferating_pct,Macrophage_NAM_pct,HSC_Activated_ECM_pct,T_total_pct,HCC_GPC3hi_pct,Tumor_subtypes_pct,Hpc_CA1_pct,Excit_total_pct,Inhib_total_pct,Astrocytes_pct,Oligodendrocytes_pct,OPCs_pct,Micro_total_pct,Microglia_DAM_pct,Endothelial_pct,Pericytes_pct,Unknown_pct"

# 一个通用 emit 函数
emit_row() {
    local sample="$1" group="$2" species="$3" parsed="$4"
    local cells=$(get_field "$parsed" "CELLS")
    local pp=$(get_pct "$parsed" "Hep_Periportal")
    local pv=$(get_pct "$parsed" "Hep_Perivenous")
    local rt=$(ratio "$(nz $pp)" "$(nz $pv)")
    local hp=$(get_pct "$parsed" "Hep_Proliferating")
    local nam=$(get_pct "$parsed" "Macrophage_NAM")
    local hsc=$(get_pct "$parsed" "HSC_Activated_ECM")
    local t_sum=$(sum_pcts "$parsed" CD4_T CD8_T Treg T_cells NKT)
    local hcc=$(get_pct "$parsed" "HCC_GPC3hi")
    local tum=$(sum_pcts "$parsed" Tumor_Hepatocyte Tumor_EMT Tumor_Proliferating Tumor_Stem)
    local ca1=$(get_pct "$parsed" "Hpc_CA1")
    local ex=$(sum_pcts "$parsed" Excitatory_L23 Excitatory_L45 Excitatory_L6)
    local inh=$(sum_pcts "$parsed" Inhibitory_Pvalb Inhibitory_Meis2 Inhibitory_Sst Inhibitory_Vip)
    local astro=$(get_pct "$parsed" "Astrocytes")
    local oligo=$(get_pct "$parsed" "Oligodendrocytes")
    local opc=$(get_pct "$parsed" "OPCs")
    local micro=$(sum_pcts "$parsed" Microglia Microglia_Homeostatic Microglia_DAM Microglia_M1 Microglia_M2)
    local dam=$(get_pct "$parsed" "Microglia_DAM")
    local endo=$(get_pct "$parsed" "Endothelial")
    local peri=$(get_pct "$parsed" "Pericytes")
    local unk=$(get_pct "$parsed" "Unknown")
    echo "$sample,$group,$species,$(nz $cells),$(nz $pp),$(nz $pv),$rt,$(nz $hp),$(nz $nam),$(nz $hsc),$t_sum,$(nz $hcc),$tum,$(nz $ca1),$ex,$inh,$(nz $astro),$(nz $oligo),$(nz $opc),$micro,$(nz $dam),$(nz $endo),$(nz $peri),$(nz $unk)"
}

for s in "${LIVER_PURE[@]}"; do
    f=$(find "$ROOT/$s" -name "*report.txt" 2>/dev/null | head -1)
    emit_row "$s" "liver_control" "mouse" "$(parse_report "$f")"
done
for s in "${LIVER_PDX[@]}"; do
    emit_row "${s}_mouse" "liver_PDX" "mouse" "$(parse_report "$ROOT/$s/mouse/${s}_mouse_report.txt")"
    emit_row "${s}_human" "liver_PDX" "human" "$(parse_report "$ROOT/$s/human/${s}_human_report.txt")"
done
for s in "${BRAIN[@]}"; do
    f=$(find "$ROOT/$s" -name "*report.txt" 2>/dev/null | head -1)
    case "$s" in
        con_brain) g="brain_control" ;;
        bla_brain) g="brain_blank" ;;
        *) g="brain_PDX" ;;
    esac
    emit_row "$s" "$g" "mouse" "$(parse_report "$f")"
done
} > "$CSV_OUT"

echo "✓ Markdown 汇总: $MD_OUT"
echo "✓ CSV 汇总: $CSV_OUT"
echo ""
echo "查看 markdown:"
echo "  cat $MD_OUT"
echo "下载 csv 到本地用 Excel 打开:"
echo "  scp -P 8022 wangmingyuan@g2.cnsknowall.com:$CSV_OUT ."
