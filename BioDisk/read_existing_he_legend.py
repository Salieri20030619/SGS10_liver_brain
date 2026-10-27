from pathlib import Path
import json,numpy as np
from PIL import Image
from scipy.ndimage import label,find_objects
r=Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis')
for s in ['Liver_3','Liver_8','gan56','Brain_3','Brain_8','nao56']:
 a=np.array(Image.open(r/'taskB_he_overlay_RCTD100'/s/'D_microenvironment_on_HE.png').convert('RGB'))
 mask=(a.max(2).astype(int)-a.min(2).astype(int)>45);mask[:200]=False;mask[:,:2060]=False
 lab,n=label(mask);res=[]
 for i,b in enumerate(find_objects(lab),1):
  if b is None:continue
  y,x=b;sub=(lab[b]==i);area=sub.sum()
  if area>40 and 7<x.stop-x.start<35 and 7<y.stop-y.start<35:
   col=np.median(a[b][sub],axis=0).astype(int);res.append({'x':x.start,'y':y.start,'area':int(area),'rgb':col.tolist()})
 print(s,json.dumps(res))
