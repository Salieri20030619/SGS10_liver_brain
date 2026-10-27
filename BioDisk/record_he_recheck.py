from pathlib import Path
from datetime import datetime
import shutil,json
r=Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis')
t=datetime.now().astimezone().isoformat(timespec='seconds')
b=r/'00_admin/state_backups'/datetime.now().strftime('%Y%m%dT%H%M%S_he_path_recheck')
b.mkdir(parents=True,exist_ok=False)
note='\n\n## H&E source path recheck ('+t+')\n\nThe user confirmed /mnt/BioData/bio/SGS10'+chr(0x4f8b)+' as the original H&E directory. The current ssh sgs10 endpoint reports hostname osv01. A server-side Python check constructed the Unicode path explicitly, excluding Windows argument-encoding ambiguity; Path.exists() returned false. /mnt/BioData is mounted from /dev/sdc1. A bounded directory search across /mnt/BioData, /mnt/BioDisk and /mnt/BioHome checked 10,222 directories and found no complete original SGS10 directory. This does not establish deletion: host or mount visibility remains to be reconciled. The user was asked for hostname and ls -ld output from the terminal where the source is visible. No H&E-underlaid version has been generated. Diagnostic: client_niche_figures_20260908_v1/he_path_recheck.json.\n'
for name in ['CODEX_HANDOFF.md','PROJECT_STATUS.md']:
 p=r/name;shutil.copy2(p,b/name)
 with p.open('a',encoding='utf-8') as f:f.write(note)
print('Updated source-path recheck in CODEX_HANDOFF.md and PROJECT_STATUS.md')
