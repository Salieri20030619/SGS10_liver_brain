from pathlib import Path
import json,hashlib,datetime,shutil
p=Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis');o=p/'Figure2_layout_v10'
v=json.loads((o/'Figure2_v10_FINAL_VALIDATION.json').read_text())
for name,h in v['outputs'].items():assert hashlib.sha256((o/name).read_bytes()).hexdigest()==h
note='''# Figure 2 v10 layout and orientation — current assembled figure

Current assembled figure: Figure2_layout_v10/Figure2_layout_orientation_v10.pdf, .png, .tiff. Final rendered PDF inspected; server hashes match the inspected local copies. E/F (and B/C) plot axes have exactly matching vertical bounds and a continuous background. ROI boxes are drawn as black/yellow vector strokes at final figure scale; zoom borders match, and connection lines originate from actual ROI boxes.

User specifically identified Bla_brain_1 as opposite to the other four brain sections. ONLY Bla_brain_1 is rotated 180 degrees as a complete H&E+overlay image. Con_brain_1, Brain_3, Brain_8, nao56 remain in v9 orientation. All liver images remain in v9 orientation because no liver reversal was confirmed. No mirroring. All original source images remain intact.

Analysis source of truth remains client_niche_RCTD100_HE_20260908_v9_corrected (original 100% RCTD labels). V10 is display/layout only: original ROI bounds, counts, labels, composition and proportions unchanged; CSV inputs verified byte-identical to the v9 archive. See Figure2_v10_FINAL_VALIDATION.json. Previous standalone 22-TIFF set remains v9. No active job.

'''
stamp=datetime.datetime.now().strftime('%Y%m%dT%H%M%S');b=p/'00_admin/state_backups'/f'{stamp}_layout_v10';b.mkdir(parents=True)
for name in ['CODEX_HANDOFF.md','PROJECT_STATUS.md']:
 f=p/name;shutil.copy2(f,b/name);f.write_text(note+'## Previous state\n\n'+f.read_text())
(o/'HANDOFF.md').write_text(note)
print('SERVER VERIFIED; HANDOFF UPDATED')
