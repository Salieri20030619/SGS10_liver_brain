from pathlib import Path
from datetime import datetime
import json, hashlib, shutil, subprocess
import pandas as pd
ROOT=Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis')
O=ROOT/'client_niche_figures_20260908_v1'
stamp=datetime.now().astimezone().isoformat(timespec='seconds')
checks=pd.read_csv(O/'validation_checks.csv')
ind=pd.read_csv(O/'independent_R_validation.csv')
assert checks['pass'].all() and len(ind)==6 and ind.status.eq('PASS').all()
figures=[O/'SGS10_liver_M05_three_tumor_mice_v1.pdf',O/'SGS10_brain_M04_three_tumor_mice_v1.pdf']
for f in figures:
 info=subprocess.check_output(['pdfinfo',str(f)],text=True)
 assert any(l.split(':',1)[1].strip()=='1' for l in info.splitlines() if l.startswith('Pages:'))
 (O/(f.stem+'_pdfinfo.txt')).write_text(info)
fonts={f.name:subprocess.check_output(['pdffonts',str(f)],text=True) for f in figures}
(O/'pdf_font_audit.json').write_text(json.dumps(fonts,indent=2))
summary=pd.read_csv(O/'sample_niche_summary.csv')
roi=pd.read_csv(O/'zoom_ROI_manifest.csv')
clip=checks[checks['check'].str.contains('clipping')]
method='''# Methods: SGS10 three-mouse niche figures

## Data and sample selection
The authoritative project is /mnt/BioDisk/wangmingyuan/SGS10_analysis.
All existing tumor-group samples were included: Liver_3, Liver_8, gan56 for liver and Brain_3, Brain_8, nao56 for brain. The biological sample size is three mice per organ; sample names are retained exactly as recorded. Blank and control samples were excluded from this requested display.

Existing RCTD100 spot assignments were read from taskB/spatial_domain_RCTD100. Liver M05 is D_new_01 (Tumor_immune_niche); brain M04 is D_new_10 (Vascular_BBB). No domains were reclustered or reassigned. RCTD coordinates were exported directly from each original rctd_results.rds object. Full per-sample weight matrices came from the matching weights_matrix.csv files; each matrix was checked against its original RDS after row normalization.

## Coordinate display and close-ups
The original per-sample coordinate system was retained. Brain coordinates were rotated 180 degrees using the historical canvas dimensions (56500 x 20434 for Brain_3/Brain_8; 55050 x 19906 for nao56), matching the accepted prior orientation. Equal x/y scaling was used. Every assigned spot is included in the whole-section map. The target niche is green and the remaining microenvironments are blue-grey. Plotted coordinates and assignments are saved in the six plotted_spots.csv.gz files.

One square close-up was selected per sample. Its side length equals 0.28 times the shorter dimension of the sample's coordinate bounding box. Candidate centers are all target spots, clipped to keep the square within that bounding box. The selected center maximizes the number of target spots within the square; ties follow lexicographic coordinate order. The exact crop is marked on the full map. ROI coordinates, spot counts and the selection rule are exported in zoom_ROI_manifest.csv. Close-up selection does not affect composition calculations, which use the entire target niche in each sample.

## Estimated cellular composition
For each spot, negative solver weights were set to zero and all weights were normalized to sum to one. The count and minimum of negative values and the maximum resulting change in each sample's niche composition are preserved in validation_checks.csv. All other original values are retained. The spatial domain labels remain the frozen upstream labels.

For each mouse and reference cell subtype, the arithmetic mean of normalized weights over all spots assigned to the target niche was multiplied by 100. Subtypes were aggregated into cell classes by summing their weights, using the complete organ-specific reference_to_cell_class.csv mapping. The mapping retains ambiguous labels such as LCM explicitly. NKT is grouped with T cells, not NK cells; liver macrophages include Kupffer reference types. All subtype and class values are exported.

The six largest classes were selected by their mean percentage across the three mice with equal weight per mouse. All remaining classes were summed into Other cell classes. All displayed rows sum to 100% within each mouse. Bars and numeric labels show each mouse separately; colors match sample headings in the spatial panels. No inferential hypothesis tests, P values, FDR values, pooled spot-based error bars or significance stars were computed.

## Verification
The builder validates sample identities, complete one-to-one coordinate/assignment/weight joins, weight validity, 100% composition totals, target definitions and ROI counts. A separate R script recomputes subtype and plotted-class composition directly from source CSVs and independently checks the ROI counts using exported RCTD coordinates. PDF pages were rendered with Poppler and visually inspected. Source hashes, full values, ROI rules, code and numerical validation are retained alongside the figures.
'''
(O/'METHODS.md').write_text(method,encoding='utf-8')
readme='''# SGS10 客户新增图：肝 M05、脑 M04

## 本次交付
- SGS10_liver_M05_three_tumor_mice_v1.pdf / .png：Liver_3、Liver_8、gan56。
- SGS10_brain_M04_three_tumor_mice_v1.pdf / .png：Brain_3、Brain_8、nao56。

每张组合图包含 A：三只小鼠的全景空间分布；B：与全景框位置一致的局部放大图；C：用三种颜色区分小鼠的目标微环境细胞组成。绿色标记肝 M05 或脑 M04，其余空间点作为蓝灰色背景。细胞组成按整片目标微环境计算，不限于放大框。

## 图注与结果定位
肝图：三只肿瘤组小鼠肝脏 M05（Tumor_immune_niche）的空间分布、局部放大及逐样本细胞组成。
脑图：三只肿瘤组小鼠脑组织 M04（Vascular_BBB）的空间分布、局部放大及逐样本细胞组成。
柱长为各小鼠目标微环境内平均 RCTD 权重（%），表示反卷积估计组成。展示按三鼠等权平均排序的前六个细胞大类及其余类别合计；各样本合计为 100%。完整亚型值和亚型到大类对应关系均保留。

## 当前图版与原图需求
当前交付为真实配准坐标绘制的空间图版，未叠加 H&E。旧原始目录 /mnt/BioData/bio/SGS10例 已不可用；项目内只找到部分样本的原始 H&E 与部分纯 H&E 渲染缓存，尚未补齐六样本整套原图。原始数据新路径已向用户询问。
因此，本版不标注未经核实的病理肿瘤核心/边缘，也不添加未经校准的 mm/um 比例尺。放大框标记的是目标微环境较密集的区域，并非经病理确认的肿瘤边缘。待取得完整配准 H&E 和物理尺度后，可据现有坐标、ROI 和冻结组成表另做 H&E 叠加版。

## 证据边界
每个器官有效独立样本为三只肿瘤组小鼠，结果按描述性、探索性证据呈现。空间点数不是生物学样本量；RCTD 权重不是实测细胞计数。本次不作组间检验，不计算 P/FDR，不据此推断因果或临床效用。既有 RCTD100 参考标签未因 9 月 7 日 Table S3 清洗而重新拟合。

## 复核文件
- METHODS.md：完整计算方法。
- niche_composition_all_celltypes.csv：全部参考亚型数值。
- liver/brain_niche_composition_all_classes.csv：全部细胞大类数值。
- liver/brain_reference_to_cell_class.csv：亚型到大类的完整映射。
- liver_M05/brain_M04_composition_plotted.csv：与柱形图逐值对应的数值。
- sample_niche_summary.csv、zoom_ROI_manifest.csv：样本及区域记录。
- validation_checks.csv、independent_R_validation.csv、rctd_source_validation.csv：数值核验。

'''
readme+=f'生成与检查时间：{stamp}。构图数值检查 {len(checks)}/{len(checks)} 通过；独立 R 复算六样本全部通过。\n'
(O/'DELIVERY_README.md').write_text(readme,encoding='utf-8')
# This finalizer is executed only after the final rendered PDFs are visually approved in-session.
qa={'status':'PASS_SPATIAL_VERSION','checked_at':stamp,'pdf_pages':2,'numerical_checks_passed':len(checks),'independent_R_samples_passed':6,'visual_checks':['Whole-section maps and matched square crops inspected for all six samples','Target green and background blue-grey are distinguishable','Per-mouse bar colors and sample headings match','Titles, legends, labels and footnotes are legible without clipping','No unsupported histological boundary or physical scale is displayed'],'pending':'H&E overlay version requires complete aligned H&E source location and physical calibration'}
(O/'FINAL_VALIDATION.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
build=json.loads((O/'build_summary.json').read_text());build['status']='PASS_SPATIAL_VERSION_HE_SOURCE_PENDING';build['final_validation']='FINAL_VALIDATION.json';(O/'build_summary.json').write_text(json.dumps(build,indent=2))
backup=ROOT/'00_admin/state_backups'/datetime.now().strftime('%Y%m%dT%H%M%S_niche_figure_preupdate')
backup.mkdir(parents=True,exist_ok=False)
state_names=['PROJECT_STATUS.md','CODEX_HANDOFF.md','DECISIONS.md','RESULTS_MANIFEST.tsv']
for name in state_names:
 if (ROOT/name).exists(): shutil.copy2(ROOT/name,backup/name)
status=f'''# PROJECT STATUS

Updated: {stamp}
Project: SGS10 liver-brain single-cell and spatial analysis
Authoritative root: {ROOT}

## Current module
Client-requested three-tumor-mouse liver M05 / brain M04 spatial figures.
Status: PASS_SPATIAL_VERSION_HE_SOURCE_PENDING.

## Completed deliverables
Two organ-specific combination figures each include all three tumor samples, full spatial maps, marked matched close-ups and a three-color per-mouse cellular-composition chart. Liver: Liver_3, Liver_8, gan56; brain: Brain_3, Brain_8, nao56. RCTD100 assignments are frozen. Composition is recalculated per sample across the whole target niche with explicit subtype-to-class mappings. Figure PDFs and high-resolution PNGs are in client_niche_figures_20260908_v1/.

## Validation
Builder: {len(checks)}/{len(checks)} checks pass. Independent R recomputation: 6/6 samples pass. Original RCTD RDS and normalized CSV weights agree to floating-point precision. Both final PDFs have one page and passed rendered visual QA. Complete data, ROI records, code and audit files are retained.

## Remaining source dependency
The delivered version is a spatial-coordinate map without H&E underlay. The old /mnt/BioData/bio/SGS10例 source directory is unavailable and the six-sample H&E set has not been recovered. The user has been asked for the relocated source directory. H&E overlay and calibrated physical scale require those source assets. No pathological tumor-margin contours are asserted in this version.

## Recent prior deliverables retained
- client_supplement_package_20260908_work/SGS10_client_supplement_images_index_20260908.zip: 20 image files, 5 requested tables, 3 index/guide files; prior delivery_check.json records 177 source dates checked.
- TableS3_clean_20260907_v1/: liver 34 subtypes, brain 23, descriptive markers; prior independent validation 17/17 PASS.
- Shared KEGG significance figure: liver_brain_v5_presto_results/08_client_figures/shared_pathway_significance_grouped_bar_msyh_v2/ remains the accepted typography version for that separate module.

## Statistical and interpretation boundary
Current module: n=3 independent tumor-bearing mice per organ, descriptive/exploratory. RCTD weight percentages are estimated composition, not counted cells. Full target niches determine composition; close-up selection affects visualization only. No inferential P/FDR or spot-as-replicate tests were performed. No organ-to-organ causal inference. Legacy differential analyses retain their own sample definitions; the current tumor-only display does not change them.

## Active jobs
None for the current module after completion.

## State backup
{backup.relative_to(ROOT)}
'''
(ROOT/'PROJECT_STATUS.md').write_text(status,encoding='utf-8')
handoff=f'''# CODEX HANDOFF

Updated: {stamp}

## Latest user request
Create one liver M05 figure and one brain M04 figure covering all three tumor-bearing mice per organ. Each organ figure must include full spatial maps, matched close-ups, and a cellular-composition comparison with the three mice distinguished by color. Update project handoff records.

## Completed now
The spatial-coordinate version is complete and validated. Outputs: client_niche_figures_20260908_v1/SGS10_liver_M05_three_tumor_mice_v1.pdf and SGS10_brain_M04_three_tumor_mice_v1.pdf, with PNG counterparts. Each includes three full maps, three corresponding close-ups and one three-mouse composition chart. Green identifies the target microenvironment. Sample headings and bar colors are blue/orange/purple in source sample order.

## Exact resume point
For a cosmetic revision, use client_niche_figures_20260908_v1/build_figures.py and write a non-overwriting v2. Preserve niche_composition_all_celltypes.csv and the organ-specific full class / plotted composition tables unless the user requests a methodological change. RCTD100 domain assignments remain authoritative. independent_validate.R independently verifies subtype/class composition and ROI point counts.

For an H&E overlay revision, first obtain the relocated six-sample aligned H&E files and their true image-coordinate dimensions / physical calibration. The old /mnt/BioData/bio/SGS10例 path is missing. The user was asked for its new location; no new location was received when this record was written. Do not synthesize histology or add an unverified micron scale. The current ROI is a deterministic target-dense crop, not a pathology-confirmed tumor margin.

## Core input and calculation paths
- spatial_results/rctd_full/per_sample/<sample>_rctd_results.rds and <sample>_weights_matrix.csv.
- taskB/spatial_domain_RCTD100/<organ>_spot_microenv_assignment.csv and <organ>_microenv_naming.csv.
- Liver M05 = D_new_01; brain M04 = D_new_10.
- Each of the six RCTD coordinate exports is stored in the new module.
- Per-spot negative solver values are clipped to zero and weights renormalized for composition only; numerical impact is logged in validation_checks.csv. Frozen domain assignments are not changed.
- Whole-target-niche mean normalized RCTD weights are summarized per mouse; six leading broad cell classes are ranked by the equal-mouse mean, with a residual Other cell classes category. Mapping CSVs preserve all subtype memberships.
- Brain maps preserve the accepted 180-degree orientation.

## Validation
{len(checks)} builder checks pass; all 6 samples pass a separate R recomputation. Final one-page PDFs were rendered with Poppler and visually reviewed. See FINAL_VALIDATION.json, independent_R_validation.csv, rctd_source_validation.csv, source_manifest.json and METHODS.md.

## Prior modules
The September 8 customer supplementary package and September 7 clean Table S3 are retained. The August 31 shared KEGG Microsoft YaHei v2 remains the current figure for its separate pathway module. Root records before this update are backed up at {backup.relative_to(ROOT)}.

## Active jobs / pending work
No running jobs for this module. Pending external input: relocated aligned H&E set if the client requires the H&E-underlaid reference style. Current spatial PDFs are reviewable now.
'''
(ROOT/'CODEX_HANDOFF.md').write_text(handoff,encoding='utf-8')
old=(backup/'DECISIONS.md').read_text(encoding='utf-8') if (backup/'DECISIONS.md').exists() else ''
decisions=f'''# DECISIONS

Updated: {stamp}

## Current client niche-figure module
1. Include all three tumor-bearing mice per organ and retain original sample names.
2. Produce separate liver M05 and brain M04 combination figures with full maps, matched zooms and one three-color composition comparison each.
3. Use frozen RCTD100 domain assignments and directly exported original RCTD coordinates.
4. Use green for the target niche and blue-grey for remaining spatial points. Use blue/orange/purple for sample identity in composition bars and sample headings.
5. Select a deterministic target-dense square for illustration, export its bounds and counts, and compute composition over the entire target niche rather than the selected crop.
6. Summarize RCTD weights into documented cell classes; keep every original subtype value and mapping for audit. Rank the six displayed classes by equal-mouse mean and combine the remainder.
7. Clip negative solver values and renormalize only for composition; quantify the effect and preserve frozen source files.
8. Use descriptive/exploratory reporting (n=3 mice per organ), with no inferential P/FDR or cell/spot pseudoreplication.
9. Deliver an explicitly identified spatial-coordinate version while the complete H&E source set is unavailable. A future H&E-underlaid version requires source recovery and verified calibration; do not invent histology, boundaries or scale bars.
10. Preserve previous delivery files and back up root state before updating it.

## Preserved decisions from the separate prior pathway module

'''+old
(ROOT/'DECISIONS.md').write_text(decisions,encoding='utf-8')
manifest=ROOT/'RESULTS_MANIFEST.tsv'
import csv
with manifest.open('r',encoding='utf-8') as f:
 fields=next(csv.reader(f,delimiter='\t'))
assert fields==['artifact_id','phase','created_at','type','path','status','provenance','sha256','notes']
with manifest.open('a',encoding='utf-8',newline='') as f:
 writer=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n')
 for i,p in enumerate(sorted(O.iterdir())):
  if p.is_file() and p.suffix in ['.pdf','.png','.csv','.md','.json','.py','.R']:
   writer.writerow({'artifact_id':f'SGS10_NICHE_20260908_V1_{i:03d}','phase':'THREE_TUMOR_MOUSE_NICHE_FIGURES','created_at':stamp,'type':p.suffix.lstrip('.'),'path':str(p),'status':'PASS_SPATIAL_VERSION','provenance':'Frozen RCTD100 assignments and original RCTD coordinates; build_figures.py','sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'notes':'H&E source set and physical calibration remain pending; descriptive per-mouse RCTD composition'})
records=[]
for p in sorted(O.iterdir()):
 if p.is_file() and p.name!='OUTPUT_MANIFEST.json': records.append({'file':p.name,'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
(O/'OUTPUT_MANIFEST.json').write_text(json.dumps(records,indent=2))
print(json.dumps({'status':qa['status'],'checks':len(checks),'independent_samples':len(ind),'state_backup':str(backup),'output':str(O)},indent=2))
