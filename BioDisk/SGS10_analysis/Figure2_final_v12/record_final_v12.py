from pathlib import Path
import json,hashlib,datetime,shutil
p=Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis');o=p/'Figure2_final_v12'
v=json.loads((o/'Figure2_v12_pdf_validation.json').read_text())
assert hashlib.sha256((o/'Figure2_final_v12.pdf').read_bytes()).hexdigest()==v['pdf_sha256']
assert hashlib.sha256((o/'Figure2_final_v12.png').read_bytes()).hexdigest()==v['source_png_sha256']
assert v['pdfium_extracted_image_equals_original_png'] and v['image_objects']==1
note='''# Figure 2 v12 — current display-compatible final figure

Current files: Figure2_final_v12/Figure2_final_v12.pdf, .png, .tiff. User confirmed PNG was normal but prior vector PDF displayed an E/F horizontal displacement in their viewer. V12 PDF embeds exactly one lossless RGB image with no tiles, masks, text layers or clipping groups. PDFium extraction verified exact pixel equality to the normal PNG. Final PDF rendering was visually inspected. PDF/PNG are 2880 x 3144 (240 dpi at native page size), TIFF is 600 dpi. The editable plot builder is retained separately.

Total title is now in a separate white header above panel A. Body is not rescaled or rearranged. Preserve v10 features: only Bla_brain_1 rotated 180 degrees; other brain and all liver orientations unchanged; black/yellow ROI frames and true ROI-to-zoom connectors; E/F axes aligned. All statistical data and original full-RCTD assignments remain v9, unchanged.

Prior v10 vector PDF is retained for audit, superseded for viewing/delivery. No active job.

'''
b=p/'00_admin/state_backups'/(datetime.datetime.now().strftime('%Y%m%dT%H%M%S')+'_final_v12');b.mkdir(parents=True)
for name in ['CODEX_HANDOFF.md','PROJECT_STATUS.md']:
 f=p/name;shutil.copy2(f,b/name);f.write_text(note+'## Previous state\n\n'+f.read_text())
(o/'HANDOFF.md').write_text(note)
print('FINAL PDF AND PNG VERIFIED; HANDOFF UPDATED')
