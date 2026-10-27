from pathlib import Path
from PIL import Image,ImageOps,ImageDraw
import json
r=Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis')
o=r/'client_niche_existing_HE_20260908_v2';o.mkdir(exist_ok=True)
ss=['Liver_3','Liver_8','gan56','Brain_3','Brain_8','nao56']
thumb=Image.new('RGB',(1600,3*500),'white');draw=ImageDraw.Draw(thumb);out=[]
for i,s in enumerate(ss):
 p=r/'taskB_he_overlay_RCTD100'/s/'D_microenvironment_on_HE.png'
 im=Image.open(p).convert('RGB');out.append({'sample':s,'file':str(p),'size':im.size,'bytes':p.stat().st_size})
 im.thumbnail((780,460));xx=(i%2)*800;yy=(i//2)*500
 thumb.paste(im,(xx,yy+32));draw.text((xx+8,yy+8),s,fill='black')
thumb.save(o/'source_contact.jpg',quality=95)
(o/'source_image_inventory.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
