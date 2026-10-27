import os, json, socket, subprocess
from collections import deque
from pathlib import Path
print('HOST',socket.gethostname(),flush=True)
print(subprocess.run(['findmnt','-T','/mnt/BioData/bio'],capture_output=True,text=True).stdout,flush=True)
q=deque([(r,0) for r in ['/mnt/BioData','/mnt/BioDisk','/mnt/BioHome']]); n=0; hits=[]
skip={'.cache','.git','.conda','node_modules','site-packages','miniconda3','anaconda3','.codex','.local','proc','sys'}
while q:
 p,depth=q.popleft()
 try:
  with os.scandir(p) as ds:
   for d in ds:
    if d.is_dir(follow_symlinks=False):
     if 'sgs10' in d.name.lower() or 'sgs2407026' in d.name.lower() or 'sgs25090559' in d.name.lower():
      print('MATCH',d.path,flush=True);hits.append(d.path)
     if depth<5 and d.name not in skip:q.append((d.path,depth+1))
  n+=1
 except (PermissionError,FileNotFoundError,OSError):pass
print('directories_checked',n,flush=True)
Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis/client_niche_figures_20260908_v1/he_path_recheck.json').write_text(json.dumps({'host':socket.gethostname(),'requested_path':'/mnt/BioData/bio/SGS10'+chr(0x4f8b),'exists':Path('/mnt/BioData/bio/SGS10'+chr(0x4f8b)).exists(),'matching_directories':hits,'directories_checked':n},ensure_ascii=False,indent=2))
