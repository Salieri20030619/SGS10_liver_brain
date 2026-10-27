from pathlib import Path
import pandas as pd,json,hashlib,shutil,datetime,zipfile
p=Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis');o=p/'client_niche_RCTD100_HE_20260908_v9_corrected';d=o/'delivery'
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
val=json.loads((o/'PROPORTION_VALIDATION.json').read_text());assert val['all_pass'] and val['labels_changed']==0
idx=pd.read_csv(d/'TIFF_index.csv');assert len(idx)==22
for r in idx.to_dict('records'):assert sha(d/r['file'])==r['sha256']
for organ,target in [('liver','M05'),('brain','M04')]:
 for tail in ['04_composition.tiff','05_microenvironment_legend.tiff','composition_plotted.csv']:
  f=Path(f'{organ}_{target}')/f'{organ}_{target}_{tail}';assert sha(d/f)==sha(p/'client_niche_blue_TIFF_20260908_v4'/f)
auditfiles=['PROPORTION_VALIDATION.json','microenvironment_proportion_before_after.csv','cell_composition_before_after.csv','original_vs_v6_v7_vs_corrected.csv','barcode_label_audit.csv']
for name in auditfiles:shutil.copy2(o/name,d/name)
note='''# SGS10 v9 corrected original-full assignment — current delivery

User requested correction of the inappropriate change of proportions during a redraw. v9 restores the ORIGINAL 100% RCTD k-means assignments in taskB/spatial_domain_RCTD100/*_spot_microenv_assignment.csv (the full-label version used by v5), without weight smoothing or reclassification. Do not use the v6/v7 reassignments or Figure2_full_v8 summaries for this requested figure.

Validation: 10 samples, 522110 full points, zero changed labels; all microenvironment spot-proportion differences exactly 0 percentage points. Six tumor samples retain all 275839 points. M05/M04 cell composition independently recomputed from full weights; largest discrepancy from original tables 1.4210854715202004e-14 percentage points (floating-point roundoff). Composition TIFFs, legends and plotted CSVs are byte-identical to the original full v4/v5 delivery. Brain M04 endothelial: Brain_3 90.10222937704941%, Brain_8 88.4565197158029%, nao56 88.07770818906658%.

Delivery: client_niche_RCTD100_HE_20260908_v9_corrected/delivery/ (22 TIFFs, 600dpi LZW), SGS10_v9_corrected_TIFF.zip, and full before/after audit tables. Original v5 HE ROI bounds restored. Continuous bounded nearest-label raster is DISPLAY ONLY. All statistics use original points/weights, never interpolated pixel area. Co-located points may have different original labels: all remain in statistics; last barcode in sorted draw order determines only the overlap display. No labels are merged or changed. Maximum support distance remains 1.5 times 95th percentile nearest-neighbor distance.

v6/v7 remain audit artifacts, superseded for this task. The v8 figure-data export is explicitly marked DO_NOT_USE. Figure2 must use v9 original assignments for all A-F panels. No inferred sham labels or unverified physical scale should be added. Current baseline labels are original Blank (Bla) / Control (Con) until experimentally confirmed. No new significance tests.

'''
(o/'CORRECTION_HANDOFF.md').write_text(note)
(p/'Figure2_full_v8/DO_NOT_USE_FOR_REQUESTED_FIGURE.md').write_text('Superseded: these summaries use v6 smoothed reassignments, an analysis change not authorized by the redraw request. Use client_niche_RCTD100_HE_20260908_v9_corrected original-full assignments and summaries instead.\n')
stamp=datetime.datetime.now().strftime('%Y%m%dT%H%M%S');b=p/'00_admin/state_backups'/f'{stamp}_corrected_v9';b.mkdir(parents=True)
for name in ['CODEX_HANDOFF.md','PROJECT_STATUS.md']:
 f=p/name;shutil.copy2(f,b/name);f.write_text(note+'## Previous state (superseded)\n\n'+f.read_text())
readme=d/'README.txt';readme.write_text(readme.read_text(encoding='utf-8-sig')+'\n同坐标存在不同原始标签时，全部点保留在统计中，显示按条码排序后的最后一点颜色处理重叠；不据此改变任何标签。\n',encoding='utf-8-sig')
with zipfile.ZipFile(o/'SGS10_v9_corrected_TIFF.zip','w',compression=zipfile.ZIP_DEFLATED,compresslevel=2) as z:
 for f in d.rglob('*'):
  if f.is_file():z.write(f,f.relative_to(d))
val.update({'TIFF_count':22,'original_composition_files_byte_identical':True,'zip_sha256':sha(o/'SGS10_v9_corrected_TIFF.zip')})
(o/'FINAL_VALIDATION.json').write_text(json.dumps(val,indent=2));print(json.dumps(val),flush=True)
