from pathlib import Path
from PIL import Image
import numpy as np,json
from scipy import ndimage
r=Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis');o=r/'client_niche_existing_HE_20260908_v2'
for s in ['Liver_3','Liver_8','gan56','Brain_3','Brain_8','nao56']:
 a=np.asarray(Image.open(r/'taskB_he_overlay_RCTD100'/s/'D_microenvironment_on_HE.png').convert('RGB'))
 # Print runs of content by row and columns to locate the published plot panel, excluding its header and legend.
 sat=a.max(axis=2).astype(float)-a.min(axis=2)
 mask=(sat>25)&(a.min(axis=2)<210);mask[:120]=False;mask[:,2060:]=False
 yy,xx=np.where(mask)
 print(s,'color_bbox',[int(xx.min()),int(yy.min()),int(xx.max()+1),int(yy.max()+1)])
 Image.fromarray(a).resize((1440,547)).save(o/(s+'_source_preview.png'))
