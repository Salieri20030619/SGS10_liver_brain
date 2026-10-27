from pathlib import Path
from datetime import datetime
import json,shutil
root=Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis');out=root/'client_niche_blue_TIFF_20260908_v4'
p=out/'validation_summary.json';data=json.loads(p.read_text());data['visual_QA']='PASS: final contact sheet inspected, spatial panels, boxes, zooms, legends and composition';p.write_text(json.dumps(data,indent=2))
backup=root/'00_admin/state_backups'/('blue_tiff_v4_'+datetime.now().strftime('%Y%m%dT%H%M%S'));backup.mkdir()
update='# Blue reference TIFF v4 — current delivery\n\nUser explicitly requested the blue version after supplying a reference image. Completed at client_niche_blue_TIFF_20260908_v4/: 22 independent 600-dpi LZW TIFFs for liver M05 and brain M04, with full maps, boxed maps, matched zooms, composition and legends. All 59 checks pass and final contact sheet visually inspected.\n\nUses frozen taskB/spatial_domain_RCTD100/*_spot_microenv_assignment.csv and validated v1 RCTD coordinates. No H&E background; original montage domain colors. Liver M05 red, brain M04 blue. Exact labels and transformed coordinates match v1. ROI uses validated v1 source-coordinate bounds; v3 HE pixel ROIs are superseded. Composition files unchanged. Whole/zoom point areas 9/90 pt^2, alpha 0.9, equal across all domains, barcode-sorted draw order. No new statistical analysis.\n\nReproduce with miniconda3/envs/pdac_spatial/bin/python SGS10_analysis/build_blue_tiff_20260908_v4.py from home. Old v3 retained, superseded for client spatial panels. No active jobs.\n\n## Previous state (preserved)\n\n'
for n in ['CODEX_HANDOFF.md','PROJECT_STATUS.md']:
 p=root/n;shutil.copy2(p,backup/n);p.write_text(update+p.read_text())
print('State saved',backup)
