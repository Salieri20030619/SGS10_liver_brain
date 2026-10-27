from pathlib import Path
from datetime import datetime
import shutil,json
r=Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis')
s=datetime.now().astimezone().isoformat(timespec='seconds')
b=r/'00_admin/state_backups'/('tiff_v3_'+datetime.now().strftime('%Y%m%dT%H%M%S'));b.mkdir()
for n in ['CODEX_HANDOFF.md','PROJECT_STATUS.md']:
 p=r/n;shutil.copy2(p,b/n)
 original=p.read_text()
 update='# TIFF panel export update\n\nUpdated: '+s+'\n\nLatest request: export previous combination figures as independent TIFF panels for customer assembly. Scope inferred from latest handoff: liver M05 / brain M04 existing-HE v2; broader paper figure scope remains unconfirmed.\n\nCompleted: client_niche_tiff_panels_20260908_v3/ contains 22 LZW TIFFs (whole native, whole with ROI, native zoom for six mice; composition and domain legend for each organ). Native raster panels retain original pixel dimensions at 300 dpi; composition and legends rendered at 600 dpi. Source hashes, lossless TIFF pixels, DPI and unchanged composition were checked. Local downloaded hashes match. Contact sheet visually inspected. Zooms remain 360 x 360 native pixels.\n\nReproduce: miniconda3/envs/pdac_spatial/bin/python SGS10_analysis/export_tiff_panels_20260908_v3.py (paths relative to home). No active export job. Previous figures retained.\n\n## Previous handoff (preserved)\n\n'
 p.write_text(update+original)
print('State updated with backups:',b)
