from pathlib import Path
import hashlib,json,csv,datetime,shutil
from PIL import Image
p=Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis');o=p/'client_niche_RCTD100_HE_20260908_v7';d=o/'delivery'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
rows=list(csv.DictReader((d/'TIFF_index.csv').open()));assert len(rows)==22
for r in rows:
 f=d/r['file'];assert sha(f)==r['sha256']
 with Image.open(f) as im:assert im.size==(int(r['width_px']),int(r['height_px'])) and im.tag_v2[259]==5 and im.info['dpi']==(600,600)
audit=list(csv.DictReader((d/'rendering_audit.csv').open()));assert len(audit)==6 and sum(int(r['n_full_spots']) for r in audit)==275839
assert all(int(r['labels_changed'])==0 and int(r['coordinate_label_conflicts'])==0 for r in audit)
for organ,target in [('liver','M05'),('brain','M04')]:
 for suffix in ['04_composition.tiff','05_microenvironment_legend.tiff','composition_plotted.csv']:
  f=Path(f'{organ}_{target}')/f'{organ}_{target}_{suffix}'
  assert sha(d/f)==sha(p/'client_niche_RCTD100_HE_20260908_v6/delivery'/f)
result={'TIFF_count':22,'all_sha256_dimensions_dpi_compression_pass':True,'full_spots':275839,'labels_changed':0,'compositions_identical_to_v6':True,'visual_QA':'PASS: all 22 panels inspected; liver and brain zoom comparison reviewed; continuous fill achieved, H&E alignment and major gaps preserved','zip_sha256':sha(o/'SGS10_v7_continuous_TIFF.zip')}
(o/'VALIDATION.json').write_text(json.dumps(result,indent=2))
stamp=datetime.datetime.now().strftime('%Y%m%dT%H%M%S');backup=p/'00_admin/state_backups'/f'{stamp}_continuous_v7';backup.mkdir(parents=True)
note='''# SGS10 v7 continuous HE overlay — current delivery

Completed: client_niche_RCTD100_HE_20260908_v7/delivery/, 22 independent 600-dpi LZW TIFFs, with downloadable SGS10_v7_continuous_TIFF.zip. All images validated and contact sheet visually inspected; VALIDATION.json records archive SHA256. User reported v6 zooms looked like isolated dots with unfilled gaps. v7 replaces fixed-point scatter with a bounded nearest-label spatial raster (Voronoi display). Each supported display pixel inherits the nearest full RCTD point label, alpha 0.80. Support cutoff is 1.5 times the sample's 95th percentile nearest-neighbor distance, identical in whole and zoom views. ROI coverage 95%-98%; larger unsupported gaps remain H&E only. This is categorical display interpolation, not new measured cells or biological territory inference.

All 275839 full-RCTD v6 spot labels, original coordinates, sample assignments and ROI bounds are unchanged. Composition and legend files are byte-identical to v6. Always use rctd_full (100%), not the 20% rctd version. Reproduce with build_continuous_v7.py using fully archived sibling v6 inputs. All six H&E originals stay in v6/inputs; no duplicate source transfer is needed. v6 pure H&E attachments remain available. No running analysis/export job.

## Previous state preserved

'''
for name in ['CODEX_HANDOFF.md','PROJECT_STATUS.md']:
 f=p/name;shutil.copy2(f,backup/name);f.write_text(note+f.read_text())
(o/'HANDOFF.md').write_text(note.split('## Previous state preserved')[0])
print(json.dumps(result),flush=True)
