from pathlib import Path
import hashlib,json,shutil,datetime
import pandas as pd
from PIL import Image
root=Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis')
out=root/'client_niche_RCTD100_HE_20260908_v6'
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
sourcefiles=[Path(x) for x in (out/'server_source_paths.txt').read_text().splitlines()]
sourcefiles += [root/'scripts/he_overlay_RCTD100.R',root/'scripts/rctd_batch_FULL_v2.R',root/'scripts/rctd_batch_10samples.R',root/'spatial_results/rctd_full/summary/FULL_v2_README.txt']
manifest={str(p):{'bytes':p.stat().st_size,'sha256':sha(p)} for p in sourcefiles}
(out/'server_source_hashes.json').write_text(json.dumps(manifest,indent=2))
snap=out/'source_code_snapshot';snap.mkdir(exist_ok=True)
for p in sourcefiles:
 if p.suffix=='.R' or p.name.endswith('READMe.txt'):shutil.copy2(p,snap/p.name)
idx=pd.read_csv(out/'delivery/TIFF_index.csv');assert len(idx)==34
for row in idx.to_dict('records'):
 p=out/'delivery'/row['file'].replace(chr(92),'/');assert sha(p)==row['sha256']
 with Image.open(p) as im:assert im.size==(row['width_px'],row['height_px']) and im.tag_v2[259]==5
audit=pd.read_csv(out/'full_RCTD_audit.csv');assert (audit.full_weight_rows==audit.plotted_rows).all();assert audit.plotted_rows.sum()==275839
(out/'server_delivery_validation.json').write_text(json.dumps({'files':34,'all_tiff_sha256_match_local':True,'all_full_spots_plotted':True,'source_hash_count':len(manifest)},indent=2))
stamp=datetime.datetime.now().strftime('%Y%m%dT%H%M%S')
backup=root/'00_admin/state_backups'/f'{stamp}_RCTD100_HE_v6';backup.mkdir(parents=True)
note='''# SGS10 full RCTD + supplied HE v6 — current delivery

User default: ALWAYS use 100% RCTD (`spatial_results/rctd_full/per_sample`), never the 20% `spatial_results/rctd` version for current figures.
Completed 2026-09-08: `client_niche_RCTD100_HE_20260908_v6/delivery/` contains 22 independent main TIFF panels plus 12 pure HE panels, 600 dpi LZW, for Liver_3/Liver_8/gan56 and Brain_3/Brain_8/nao56. Main targets liver M05 (red), brain M04 (blue). All 105 builder checks pass; contact sheet inspected. Server SHA256 validation passes for all 34 TIFFs. No active job.

Audit correction: previous v5 was ALREADY full RCTD by exact barcode identity; it used unsmoothed k-means assignments. Historical HE panels instead used `smooth_weights` (self+15 neighbors) and cosine-profile `assign_domain`. v6 directly reuses those exact R functions, reading the six full RDS objects and validating weight CSV agreement (<1e-10). Counts 42456/39279/25492/53570/56356/58686, total 275839; corresponding 20% counts 8491/7856/5098/10714/11271/11737. Full means all retained RCTD spots after mouse UMI>=100 filtering, not all raw spatial spots.

New target counts: 3939/863/684/575/634/504. Do not mix v5 composition with v6 maps. Composition recalculated from original normalized weights within v6 whole-sample target labels; original class categories/order retained, tiny solver negatives clipped for composition only. New ROIs use prior deterministic maximum-target square rule on v6 labels; no pathology boundary inference or new tests. Prior analytical results remain preserved.

HE images from customer's resupplied aligned_HE PNGs are retained in `inputs/`; no stain recoloring. Six-image mapping is x_img=x*W/size_x, y_img=(size_y-y)*H/size_y. Same verified v5 alignment, unchanged raw coordinates. Color differences in gan56 HE remain the source appearance.

Reproduce: R 4.5.2 `rebuild_full.R` (server full objects + historical functions), then Python `build_panels.py` (numpy/pandas/scipy/matplotlib/Pillow). Folder includes self-contained plotting inputs, original mapping references, hashes, R session information, old/new assignment and composition audits, and source code snapshots. `finalize.py` is the local packaging utility. Customer-ready archive is local C:/codex_test/SGS10_RCTD100_HE_v6_TIFF.zip. Local complete work directory C:/codex_test/SGS10_RCTD100_HE_v6.

## Previous state preserved

'''
for name in ['CODEX_HANDOFF.md','PROJECT_STATUS.md']:
 p=root/name
 if p.exists():
  shutil.copy2(p,backup/name);old=p.read_text()
 else:old=''
 p.write_text(note+old)
(out/'HANDOFF_v6.md').write_text(note.split('## Previous state preserved')[0])
print('SERVER DELIVERY VERIFIED; HANDOFF UPDATED',flush=True)
