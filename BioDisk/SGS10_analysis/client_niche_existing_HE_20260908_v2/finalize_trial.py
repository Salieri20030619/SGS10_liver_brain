from pathlib import Path
from datetime import datetime
import json,hashlib,csv,shutil,subprocess
import numpy as np,pandas as pd
from PIL import Image
R=Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis');O=R/'client_niche_existing_HE_20260908_v2';V1=R/'client_niche_figures_20260908_v1'
t=datetime.now().astimezone().isoformat(timespec='seconds')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
c=pd.read_csv(O/'validation_checks.csv');assert c['pass'].all()
r=pd.read_csv(O/'crop_manifest.csv');assert len(r)==6
ind=[]
for row in r.to_dict('records'):
 s=row['sample'];target=row['target'];src=Image.open(row['source']).convert('RGB')
 box=tuple(int(row[x]) for x in ['zoom_left','zoom_top','zoom_right','zoom_bottom'])
 actual=Image.open(O/'native_crops'/f'{s}_{target}_zoom_native.png').convert('RGB')
 assert np.array_equal(np.asarray(src.crop(box)),np.asarray(actual))
 assert src.size==(2400,911) and actual.size==(360,360)
 # Every original annotation intersecting the crop is entirely contained.
 x0=box[0]-row['whole_panel_left'];y0=box[1]-row['whole_panel_top']
 for l,top,rr,b in json.loads(row['annotation_bounds_in_whole_panel']):
  inter=x0<rr and x0+360>l and y0<b and y0+360>top
  contains=x0<=l and x0+360>=rr and y0<=top and y0+360>=b
  assert not inter or contains
 ind.append({'sample':s,'exact_crop_pixel_match':True,'complete_intersected_annotations':True,'native_crop_px':'360x360'})
for organ,target in [('liver','M05'),('brain','M04')]:
 f=f'{organ}_{target}_composition_plotted.csv';assert sha(O/f)==sha(V1/f)
 pdf=O/f'SGS10_{organ}_{target}_existing_HE_trial_v2.pdf'
 info=subprocess.check_output(['pdfinfo',str(pdf)],text=True)
 assert any(z.split(':',1)[1].strip()=='1' for z in info.splitlines() if z.startswith('Pages:'))
 (O/f'{organ}_pdfinfo.txt').write_text(info)
pd.DataFrame(ind).to_csv(O/'independent_crop_validation.csv',index=False)
readme='''# SGS10 现有 H&E 图直接裁剪试版 v2

本版按用户“直接在现有图上画”的要求制作，供审阅显示效果。

- 肝 M05：Liver_3、Liver_8、gan56；保留原图红色 M05。
- 脑 M04：Brain_3、Brain_8、nao56；保留原图蓝色 M04。
- 每张组合图包含三个全景、三个对应放大图和三色逐鼠细胞组成比较。

## 方法
来源为 taskB_he_overlay_RCTD100/<sample>/D_microenvironment_on_HE.png，六图均为 2400 × 911 像素。截取组织面板，移除面板外的旧标题与旧图例，并保留组织内部原有点位、H&E、配色、方向、文字、箭头和星形标注。新图例颜色直接读取每张源图的图例色块，三样本间核对一致。

每个放大图为原 PNG 中 360 × 360 像素的直接裁剪。在组织面板内选取目标颜色像素最多的方框，同时要求任何进入方框的原有文字、箭头和星形标注完整保留。目标颜色来自源图图例，RGB 距离阈值为 45；平局按图像行列顺序选择。该规则只用于选取展示区域，颜色像素数不用于细胞数量或统计推断。确切像素边界、来源哈希及显示放大倍数保存在 crop_manifest.csv。

细胞组成直接使用已验证 v1 的逐鼠整片目标微环境数值；CSV 与 v1 字节及 SHA256 完全一致。采用前六个细胞大类加其他类别，每鼠合计为 100%。柱长表示目标微环境内平均归一化 RCTD 权重，样本量为每器官三只独立肿瘤组小鼠。数值不是从放大图像估计的，本次未重新运行 RCTD、空间分区或推断性检验。

## 显示范围
这是现有已合成科研图的裁剪排版试版，原图中点位覆盖 H&E 的程度保持原样。放大显示约 3–4 倍，不能新增原 PNG 未保留的组织细节。显示倍数不是显微镜物镜倍数。本版保留源图方向，脑图方向与 v1 新绘制空间图不同。未在未经物理校准的图片上添加 mm/um 比例尺。

## 文件
两张 PDF 为单页组合图；同名 PNG 为 250 dpi 输出。native_crops/ 保存每个样本的原生分辨率全景裁片和放大裁片。source_manifest.json、crop_manifest.csv、validation_checks.csv 和 independent_crop_validation.csv 记录来源、像素边界和核验结果。完整参考亚型数值和到细胞大类的映射沿用 v1。
'''
(O/'DELIVERY_README.md').write_text(readme,encoding='utf-8')
qa={'status':'TRIAL_READY_EXISTING_HE_V2','checked_at':t,'builder_checks_passed':len(c),'independent_source_crop_checks':ind,'composition_csv_byte_identical_to_v1':True,'PDF_pages_per_figure':1,'visual_QA':'Both final PDFs rendered and inspected; all six whole images and matched close-ups readable, source text retained without cutoff, composition values and legends legible.','native_source_px':[2400,911],'native_zoom_px':[360,360],'user_acceptance':'pending; requested trial completed'}
(O/'FINAL_VALIDATION.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
backup=R/'00_admin/state_backups'/datetime.now().strftime('%Y%m%dT%H%M%S_existing_HE_trial_preupdate');backup.mkdir(parents=True)
for f in ['PROJECT_STATUS.md','CODEX_HANDOFF.md','DECISIONS.md','RESULTS_MANIFEST.tsv']:shutil.copy2(R/f,backup/f)
(R/'PROJECT_STATUS.md').write_text(f'''# PROJECT STATUS

Updated: {t}
Project root: {R}

## Current module
Status: TRIAL_READY_EXISTING_HE_V2 (user aesthetic review pending).
The user authorized trying direct crops of existing H&E-plus-microenvironment figures. Two organ-specific trial combinations are complete at client_niche_existing_HE_20260908_v2/: liver M05 and brain M04, each covering all three tumor-bearing mice, whole-section panels, marked matched close-ups and the validated three-color cellular-composition comparison.

## Completed and validated
Sources are all six taskB_he_overlay_RCTD100/<sample>/D_microenvironment_on_HE.png images (2400 x 911 px). Liver M05 retains red; brain M04 retains blue. Original image colors, internal annotations and orientations are preserved. All six 360 x 360 px crops match source pixels exactly. {len(c)} builder checks pass; a separate crop validator confirms pixel identity and complete intersected annotations. Both PDFs are one page and passed rendered visual QA. Composition CSV hashes are unchanged from v1.

## Resume point
Review the two existing-HE trial PDFs/PNGs. Further cosmetic adjustments should create non-overwriting v3. Use build_existing_he_trial.py, crop_manifest.csv and the native_crops directory. Do not alter original pixels to remove old labels or invent histology. The requested trial is complete; recovery of raw H&E is not a dependency for this version.

## Retained prior deliverables
- client_niche_figures_20260908_v1/: validated spatial-coordinate version and full per-mouse composition data; retained.
- client_supplement_package_20260908_work/: September 8 client supplement images and tables; retained.
- TableS3_clean_20260907_v1/: clean marker table; retained.
- Shared KEGG Microsoft YaHei grouped-bar v2 remains the accepted version for its separate pathway module.

## Evidence and display boundary
Three independent tumor-bearing mice per organ. Composition is descriptive mean RCTD weight from the whole target niche, not measured cell counts or a pixel-derived estimate. No inferential P/FDR were recomputed. Native resolution limits zoom detail. Source-image brain orientation is preserved and differs from v1. No unverified physical scale was added.

## Active jobs
None for this module.

## State backup
{backup.relative_to(R)}
''',encoding='utf-8')
(R/'CODEX_HANDOFF.md').write_text(f'''# CODEX HANDOFF

Updated: {t}

## Latest user authorization
After the raw H&E directory could not be accessed, the user asked whether existing figures could be used, then explicitly requested a trial. That trial is now complete; no further raw-source clarification is needed to inspect or revise it.

## Current reviewable output
client_niche_existing_HE_20260908_v2/SGS10_liver_M05_existing_HE_trial_v2.pdf
client_niche_existing_HE_20260908_v2/SGS10_brain_M04_existing_HE_trial_v2.pdf
Same-name PNGs also exist. Each figure contains three full H&E overlay panels, three boxed close-ups and one per-mouse composition comparison. Status: TRIAL_READY_EXISTING_HE_V2; user aesthetic acceptance pending.

## Inputs and method
Existing flattened images: taskB_he_overlay_RCTD100/<sample>/D_microenvironment_on_HE.png (all 2400 x 911 pixels).
Samples: Liver_3, Liver_8, gan56; Brain_3, Brain_8, nao56.
Liver target M05 is red; brain target M04 is blue. Source palettes are read from original legend swatches and checked across samples. Source orientation and all internal annotations are retained. No additional brain rotation is applied.
360 x 360 native-pixel crops maximize target-color content subject to preserving complete intersected old black/yellow annotations. All crop boundaries are in crop_manifest.csv; native_crops/ provides unaltered crops. This image-based ROI selection is visualization only. Composition comes unchanged from validated v1 CSVs and covers the whole niche per mouse.

## Validation and reproducibility
build_existing_he_trial.py is the reproducible builder. {len(c)} builder checks pass. independent_crop_validation.csv separately verifies source pixel equality and complete annotations for all six samples. Composition CSVs are byte-identical to v1. Both PDF renderings were visually inspected. See FINAL_VALIDATION.json and DELIVERY_README.md.

## Next step
Show the trial to the user. For changes, create v3 without overwriting v1 or v2. Preserve composition tables unless the user requests a new analytical method. Avoid AI reconstruction, inpainting or fake physical scales in scientific panels. Raw H&E source access remains useful only if a fresh high-resolution overlay is subsequently requested.

## Prior context and active jobs
No active job remains for this trial. v1 spatial-only figures, the September 8 supplement package, clean Table S3 and the separate shared-KEGG font revision are retained. Prior root state is backed up at {backup.relative_to(R)}.
''',encoding='utf-8')
old=(backup/'DECISIONS.md').read_text(encoding='utf-8')
(R/'DECISIONS.md').write_text(f'''# DECISIONS

Updated: {t}

## Existing-H&E trial v2
1. User authorized direct use of existing H&E overlay figures; do not make raw-source availability a blocker for this trial.
2. Preserve source pixels, colors, annotations and orientations; crop and arrange only. Liver M05 stays red and brain M04 stays blue.
3. Add matched boxes to the whole-section panels and use exact 360 x 360 pixel crops. Avoid cutting through old text, arrow and star groups.
4. Use the frozen validated v1 per-mouse composition CSVs without changing values. Image pixel counts are not biological measurements.
5. Keep v1 and deliver v2 as an explicit trial pending visual acceptance.
6. Retain crop coordinates, native patches, source hashes, code and validation. No unverified micron scale or synthesized tissue detail.

## Historical decisions

'''+old,encoding='utf-8')
manifest=R/'RESULTS_MANIFEST.tsv'
with manifest.open(encoding='utf-8') as f:fields=next(csv.reader(f,delimiter='\t'))
assert len(fields)==9 and fields[0]=='artifact_id'
with manifest.open('a',encoding='utf-8',newline='') as f:
 writer=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n')
 for i,p in enumerate(sorted(O.rglob('*'))):
  if p.is_file():writer.writerow({'artifact_id':f'SGS10_EXISTING_HE_TRIAL_V2_{i:03d}','phase':'EXISTING_HE_CROP_TRIAL','created_at':t,'type':p.suffix.lstrip('.'),'path':str(p),'status':'TRIAL_READY','provenance':'Direct native-pixel crops of existing RCTD100 H&E overlay figures; frozen v1 composition','sha256':sha(p),'notes':'User-requested trial; raw histology not reconstructed'})
files=[{'file':str(p.relative_to(O)),'bytes':p.stat().st_size,'sha256':sha(p)} for p in sorted(O.rglob('*')) if p.is_file() and p.name!='OUTPUT_MANIFEST.json']
(O/'OUTPUT_MANIFEST.json').write_text(json.dumps(files,indent=2))
print(json.dumps({'status':qa['status'],'checks':len(c),'independent_crops':len(ind),'backup':str(backup),'zoom_factors':r.display_enlargement_factor.round(2).tolist()},indent=2))
