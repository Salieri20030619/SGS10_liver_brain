from pathlib import Path
import json,hashlib
p=Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis');o=p/'client_niche_RCTD100_HE_20260908_v9_corrected'
m=json.loads((o/'Figure2_manifest.json').read_text())
for name,h in m.items():assert hashlib.sha256((o/name).read_bytes()).hexdigest()==h
v=json.loads((o/'FINAL_VALIDATION.json').read_text());v['Figure2_files_match_local_sha256']=True;v['Figure2_PDF_visual_QA']='PASS: final PDF rendered locally with PDFium; all A-F panels inspected; no overlap or clipping.'
(o/'FINAL_VALIDATION.json').write_text(json.dumps(v,indent=2))
note='''## Figure 2 A-F completed using corrected original-full assignments

Figure2_corrected_original_full.pdf / .png / .tiff are complete; server copies match local SHA256. Final PDF render visually inspected. All ten samples use original assignments. Blank (Bla) / Control (Con) labels retained pending experimental clarification; no unverified scale bars. B/E bubble colors are equal-sample mean original weights within each microenvironment, bubble sizes are pooled original point counts. C/F retain original full M05/M04 compositions. See CORRECTION_REPORT.md for complete audit and methods. Local Figure2 builder and local HE inputs are retained in C:/codex_test/SGS10_corrected_v9 and C:/Users/31937/Downloads/SGS10_HE_rebuild_20260908_v5/inputs. No active job.

'''
for name in ['CODEX_HANDOFF.md','PROJECT_STATUS.md']:
 f=p/name;t=f.read_text();f.write_text(t.replace('\n\n','\n\n'+note,1))
print('FIGURE2 SERVER HASHES VERIFIED; HANDOFF UPDATED',flush=True)
