from pathlib import Path
import json,hashlib,shutil,zipfile
import numpy as np,pandas as pd
from scipy.spatial import cKDTree
from PIL import Image,ImageDraw,ImageColor
BASE=Path(__file__).resolve().parent
if not (BASE/'Liver_3_HE_mapped_spots.csv.gz').exists():
 BASE=BASE.parent/'client_niche_RCTD100_HE_20260908_v6'
OUT=BASE.parent/'client_niche_RCTD100_HE_20260908_v7';OUT.mkdir(exist_ok=True)
D=OUT/'delivery';D.mkdir(exist_ok=True)
PAL={'liver':dict(zip(['M01','M02','M03','M04','M05','M06','M07'],['#2C7FB8','#7FCDBB','#41B6C4','#A1D99B','#E31A1C','#FDAE6B','#9E9AC8'])), 'brain':dict(zip(['M01','M02','M03','M04','M05','M06'],['#E31A1C','#33A02C','#A6761D','#1F78B4','#FF7F00','#984EA3']))}
SAMPLES={'liver':['Liver_3','Liver_8','gan56'],'brain':['Brain_3','Brain_8','nao56']}
roi=pd.read_csv(BASE/'HE_ROI_and_mapping.csv').set_index('sample')
index=[];audit=[]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(im,path):
 im.save(path,compression='tiff_lzw',dpi=(600,600))
 with Image.open(path) as q:assert q.size==im.size and q.tag_v2[259]==5 and q.info['dpi']==(600,600)
 index.append({'file':path.relative_to(D).as_posix(),'sha256':sha(path),'width_px':im.width,'height_px':im.height,'dpi':600})
def render(im,tree,colors,bounds,size,sx,sy,cutoff):
 x0,y0,x1,y1=bounds;w,h=size
 pixelbounds=(x0*im.width/sx,y0*im.height/sy,x1*im.width/sx,y1*im.height/sy)
 bg=im.transform(size,Image.Transform.EXTENT,pixelbounds,resample=Image.Resampling.BICUBIC)
 rgb=np.asarray(bg).copy();support_count=0
 xs=x0+(np.arange(w)+.5)*(x1-x0)/w
 for start in range(0,h,128):
  stop=min(h,start+128);ys=y0+(np.arange(start,stop)+.5)*(y1-y0)/h
  xx,yy=np.meshgrid(xs,ys);dist,idx=tree.query(np.c_[xx.ravel(),yy.ravel()],workers=2)
  support=(dist<=cutoff).reshape(stop-start,w);support_count+=support.sum()
  field=colors[idx].reshape(stop-start,w,3)
  blend=np.rint(.80*field+.20*rgb[start:stop]).astype(np.uint8)
  rgb[start:stop]=np.where(support[:,:,None],blend,rgb[start:stop])
 return Image.fromarray(rgb),bg,float(support_count/(w*h))
for organ,samples in SAMPLES.items():
 target='M05' if organ=='liver' else 'M04';dest=D/f'{organ}_{target}';dest.mkdir(exist_ok=True)
 for sample in samples:
  dat=pd.read_csv(BASE/f'{sample}_HE_mapped_spots.csv.gz').set_index('barcode')
  xy=dat[['x','y']].to_numpy();tree=cKDTree(xy)
  nn=tree.query(xy,k=2)[0][:,1];cutoff=float(1.5*np.quantile(nn,.95))
  color=np.array([ImageColor.getrgb(PAL[organ][x]) for x in dat.microenv_id],dtype=np.uint8)
  # Exact label retrieval at observed locations; duplicate coordinates must have equal labels.
  retrieved=tree.query(xy)[1];conflicts=int((dat.microenv_id.to_numpy()[retrieved]!=dat.microenv_id.to_numpy()).sum())
  assert conflicts==0, (sample,'conflicting duplicate coordinates')
  stem={'gan56':'gan_56','nao56':'nao_56'}.get(sample,sample)
  ip=BASE/'inputs'/f'{stem}_aligned_HE.png';im=Image.open(ip).convert('RGB')
  sx,sy=(55050,19906) if sample in ['gan56','nao56'] else (56500,20434)
  r=roi.loc[sample];bounds=np.array([r.x_min,r.y_min,r.x_max,r.y_max])
  whole,bg,frac=render(im,tree,color,[0,0,sx,sy],(3900,round(3900*sy/sx)),sx,sy,cutoff)
  save(whole,dest/f'{sample}_01_HE_whole.tiff')
  boxed=whole.copy();dd=ImageDraw.Draw(boxed);box=tuple(bounds*np.array([3900/sx,whole.height/sy,3900/sx,whole.height/sy]))
  dd.rectangle(box,outline='white',width=18);dd.rectangle(box,outline='#172B40',width=8)
  save(boxed,dest/f'{sample}_02_HE_whole_with_ROI.tiff')
  zoom,zbg,zfrac=render(im,tree,color,bounds,(2100,2100),sx,sy,cutoff)
  save(zoom,dest/f'{sample}_03_HE_zoom.tiff')
  # Compact matched comparison for visual QA.
  comparison=Image.new('RGB',(1800,600),'white');draw=ImageDraw.Draw(comparison)
  for j,(label,pic) in enumerate([('v6 fixed-size scatter',Image.open(BASE/'delivery'/f'{organ}_{target}'/f'{sample}_03_HE_zoom.tiff').convert('RGB')),('v7 continuous spatial fill',zoom.copy()),('H&E reference',zbg.copy())]):
   pic.thumbnail((580,560));comparison.paste(pic,(600*j+(600-pic.width)//2,35));draw.text((600*j+10,10),sample+' | '+label,fill='black')
  comparison.save(OUT/f'{sample}_zoom_comparison.jpg',quality=92)
  audit.append({'sample':sample,'n_full_spots':len(dat),'labels_changed':0,'coordinate_label_conflicts':conflicts,'nearest_neighbor_median':float(np.median(nn)),'nearest_neighbor_p95':float(np.quantile(nn,.95)),'maximum_fill_distance':cutoff,'whole_support_fraction':frac,'ROI_support_fraction':zfrac,'assignment_sha256':sha(BASE/f'{sample}_HE_mapped_spots.csv.gz'),'HE_sha256':sha(ip)})
  print(sample,'DONE',len(dat),'full spots; ROI filled fraction',round(zfrac,4),flush=True)
 for suffix in ['04_composition.tiff','05_microenvironment_legend.tiff']:
  name=f'{organ}_{target}_{suffix}';src=BASE/'delivery'/f'{organ}_{target}'/name;p=dest/name;shutil.copy2(src,p)
  with Image.open(p) as q:index.append({'file':p.relative_to(D).as_posix(),'sha256':sha(p),'width_px':q.width,'height_px':q.height,'dpi':600})
  assert sha(p)==sha(src)
 src=BASE/'delivery'/f'{organ}_{target}'/f'{organ}_{target}_composition_plotted.csv';shutil.copy2(src,dest/src.name)
assert len(index)==22 and sum(x['n_full_spots'] for x in audit)==275839
pd.DataFrame(index).to_csv(D/'TIFF_index.csv',index=False)
pd.DataFrame(audit).to_csv(D/'rendering_audit.csv',index=False)
shutil.copy2(BASE/'HE_ROI_and_mapping.csv',D/'HE_ROI_and_mapping.csv')
files=sorted(D.glob('*/*.tiff'));canvas=Image.new('RGB',(1500,8*285),'#eeeeee');draw=ImageDraw.Draw(canvas)
for i,p in enumerate(files):
 with Image.open(p) as im:
  im=im.convert('RGB');im.thumbnail((480,245));x=i%3*500;y=i//3*285;canvas.paste(im,(x+(500-im.width)//2,y+28));draw.text((x+8,y+7),p.stem,fill='black')
canvas.save(D/'preview_contact.jpg',quality=92)
(D/'README.txt').write_text('''SGS10 v7：全量 RCTD + 连续空间色块，独立 TIFF 小图

修正 v6 局部图固定点径造成的空隙。每个有数据支持的显示像素继承最近全量 RCTD 点的 v6 微环境标签，形成有限距离的 Voronoi 连续铺色；不再绘制固定直径圆点。最大填充距离=该样本最近邻距离95百分位的1.5倍。全景与局部使用同一空间距离阈值、同一标签场和0.80透明度，不跨越超出阈值的大空白。此为显示插值，不代表测量到了新增细胞或连续组织分区。

6个样本、275839个全量点、微环境标签、坐标、ROI及组成数值均沿用已核验v6。原始HE保持原色。22张600dpi LZW TIFF：18张全景/带框/局部图，以及4张组成图/图例。v6纯HE附件仍可使用。组成图和图例文件与v6逐字节一致，未重新推断分类或计算统计结果。精确参数、源文件哈希、覆盖比例见 rendering_audit.csv。
''',encoding='utf-8-sig')
shutil.copy2(Path(__file__),OUT/Path(__file__).name)
with zipfile.ZipFile(OUT/'SGS10_v7_continuous_TIFF.zip','w',compression=zipfile.ZIP_DEFLATED,compresslevel=2) as z:
 for p in D.rglob('*'):
  if p.is_file():z.write(p,p.relative_to(D))
print('ALL COMPLETE',flush=True)
