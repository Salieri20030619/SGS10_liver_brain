from pathlib import Path
import tarfile,shutil,hashlib,json,datetime,csv
root=Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis')
out=root/'client_niche_RCTD100_HE_20260908_v6';out.mkdir(exist_ok=True)
prefix=out.name+'/delivery/'
# Only complete customer deliverables from the completed leading part of the upload.
# No partial member is accepted: exact sizes and original hashes are verified below.
with tarfile.open(root/'SGS10_RCTD100_HE_v6_server.tar','r|') as tf:
 try:
  for m in tf:
   if m.name.startswith(out.name+'/inputs'):break
   if not m.isfile() or not m.name.startswith(prefix):continue
   dest=out/'delivery'/m.name[len(prefix):]
   assert dest.resolve().is_relative_to((out/'delivery').resolve())
   dest.parent.mkdir(parents=True,exist_ok=True)
   with tf.extractfile(m) as src,dest.open('wb') as target:shutil.copyfileobj(src,target)
   assert dest.stat().st_size==m.size
 except tarfile.ReadError:pass
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
rows=list(csv.DictReader((out/'delivery/TIFF_index.csv').open(encoding='utf-8-sig')))
assert len(rows)==34
for r in rows:assert sha(out/'delivery'/r['file'].replace(chr(92),'/'))==r['sha256']
assert (out/'delivery/README.txt').is_file()
(out/'delivery/SERVER_VALIDATED.json').write_text(json.dumps({'complete_TIFFs':34,'all_original_sha256_match':True,'status':'FIGURES_READY; ancillary full source archive upload continues separately'},indent=2))
stamp=datetime.datetime.now().strftime('%Y%m%dT%H%M%S');bak=root/'00_admin/state_backups'/f'{stamp}_v6_figures_ready';bak.mkdir(parents=True)
note='''# SGS10 v6 full RCTD HE figures — READY

Current figures: client_niche_RCTD100_HE_20260908_v6/delivery/ (22 main TIFFs + 12 pure HE TIFFs, 600 dpi LZW).
All 34 server TIFFs have been verified against local original SHA256 hashes; customer figures are complete and usable.
Default MUST use spatial_results/rctd_full/per_sample, SUBSAMPLE_RATE=1.0, never the 20% rctd directory.
Previous v5 already used full barcodes but unsmoothed clustering assignments. v6 restores the exact historical HE script's self+15-neighbor weight smoothing and cosine domain assignment. Six-sample full coverage is 42456/39279/25492/53570/56356/58686, total 275839. Composition and ROI are updated to the same v6 labels. See delivery/README.txt and full_RCTD_audit.csv.

Ancillary background task: complete original HE + reproduction archive upload still running from local C:/codex_test/SGS10_RCTD100_HE_v6_server.tar. It is NOT a computation job. A single shell completion monitor waits for this upload, then extracts the full archive, runs record_remote.py, verifies again, records source hashes and replaces this note with final archive-complete handoff. Do not start duplicate upload/rebuild jobs. Local final customer package: C:/codex_test/SGS10_RCTD100_HE_v6_TIFF.zip.

## Previous state preserved

'''
for name in ['CODEX_HANDOFF.md','PROJECT_STATUS.md']:
 p=root/name;old=p.read_text() if p.exists() else ''
 if p.exists():shutil.copy2(p,bak/name)
 p.write_text(note+old)
print('FIGURES READY ON SERVER: all 34 TIFF SHA256 checks passed. Source-only archive sync continues.',flush=True)
