from pathlib import Path
import json,hashlib,shutil,csv
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from PIL import Image,ImageDraw

ROOT=Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis')
V1=ROOT/'client_niche_figures_20260908_v1'
V3=ROOT/'client_niche_tiff_panels_20260908_v3'
OUT=ROOT/'client_niche_blue_TIFF_20260908_v4';OUT.mkdir(exist_ok=True)
SAMPLES={'liver':['Liver_3','Liver_8','gan56'],'brain':['Brain_3','Brain_8','nao56']}
PAL={'liver':dict(zip(['M01','M02','M03','M04','M05','M06','M07'],['#2C7FB8','#7FCDBB','#41B6C4','#A1D99B','#E31A1C','#FDAE6B','#9E9AC8'])), 'brain':dict(zip(['M01','M02','M03','M04','M05','M06'],['#E31A1C','#33A02C','#A6761D','#1F78B4','#FF7F00','#984EA3']))}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'savefig.facecolor':'white'})
index=[];checks=[];provenance={}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def record(p):provenance[str(p)]=sha(p);return p
def check(name,ok):
 checks.append({'check':name,'pass':bool(ok)})
 if not ok:raise AssertionError(name)
def register(p):
 with Image.open(p) as im:
  im.load();check(p.name+' TIFF 600dpi LZW',im.format=='TIFF' and abs(im.info['dpi'][0]-600)<.01 and im.tag_v2[259]==5)
  index.append({'file':str(p.relative_to(OUT)),'width_px':im.width,'height_px':im.height,'dpi':600,'sha256':sha(p)})
def save(fig,p):
 fig.savefig(p,dpi=600,pil_kwargs={'compression':'tiff_lzw'})
 plt.close(fig);register(p)
roi_table=pd.read_csv(record(V1/'zoom_ROI_manifest.csv'))
summary=[]
for organ,samples in SAMPLES.items():
 target='M05' if organ=='liver' else 'M04';dest=OUT/f'{organ}_{target}';dest.mkdir(exist_ok=True)
 assignments=pd.read_csv(record(ROOT/'taskB/spatial_domain_RCTD100'/f'{organ}_spot_microenv_assignment.csv'))
 for sample in samples:
  a=assignments.loc[assignments['sample'].eq(sample)].copy()
  a['barcode']=a.spot_id.str.removeprefix(sample+'__');a=a.set_index('barcode')
  co=pd.read_csv(record(V1/f'{sample}_rctd_coords.csv')).set_index('barcode')
  frozen=pd.read_csv(record(V1/f'{sample}_plotted_spots.csv.gz'),index_col=0)
  check(sample+' complete unique assignment-coordinate join',a.index.is_unique and co.index.is_unique and set(a.index)==set(co.index)==set(frozen.index))
  dat=a.join(co[['x','y']]);check(sample+' frozen composition labels unchanged',dat.microenv_id.equals(frozen.loc[dat.index,'microenv_id']))
  if organ=='brain':
   sx,sy=(55050,19906) if sample=='nao56' else (56500,20434)
   dat.x=sx-dat.x;dat.y=sy-dat.y
  check(sample+' coordinates equal validated v1',np.array_equal(dat[['x','y']].values,frozen.loc[dat.index,['x','y']].values))
  # Same barcode-sorted draw order as merge() in the original R montage.
  dat=dat.sort_index();check(sample+' all labels have original montage colors',dat.microenv_id.isin(PAL[organ]).all())
  r=roi_table.loc[roi_table['sample'].eq(sample)].iloc[0]
  bounds=np.array([r.x_min,r.y_min,r.x_max,r.y_max]);inside=dat.x.between(bounds[0],bounds[2])&dat.y.between(bounds[1],bounds[3])
  check(sample+' ROI target and total counts match',int(inside.sum())==r.all_spots_in_ROI and int((inside&dat.microenv_id.eq(target)).sum())==r.target_spots_in_ROI)
  summary.append({'organ':organ,'sample':sample,'n_spots':len(dat),'n_target_spots':int(dat.microenv_id.eq(target).sum()),'ROI_n_spots':int(inside.sum())})
  for mode in ['01_whole','02_whole_with_ROI','03_zoom']:
   zoom=mode=='03_zoom';d=dat.loc[inside] if zoom else dat
   if zoom:
    limits=bounds;w=h=3.5
   else:
    limits=np.array([dat.x.min(),dat.y.min(),dat.x.max(),dat.y.max()],float)
    w=6.5;h=w*(limits[3]-limits[1])/(limits[2]-limits[0])
   fig=plt.figure(figsize=(w,h));ax=fig.add_axes([.006,.006,.988,.988])
   # Fixed equal point size for every domain; actual domain colors, no artificial blue background.
   pointsize=90.0 if zoom else 9.0
   ax.scatter(d.x,d.y,c=d.microenv_id.map(PAL[organ]).tolist(),s=pointsize,alpha=.9,linewidths=0,edgecolors='none')
   ax.set_xlim(limits[[0,2]]);ax.set_ylim(limits[[3,1]]);ax.set_aspect('equal');ax.axis('off')
   if mode=='02_whole_with_ROI':
    for color,lw in [('white',2.4),('#172B40',1.1)]:ax.add_patch(Rectangle(bounds[:2],bounds[2]-bounds[0],bounds[3]-bounds[1],fill=False,edgecolor=color,linewidth=lw))
   save(fig,dest/f'{sample}_{mode}.tiff')
  dat[['spot_id','microenv_id','x','y']].to_csv(OUT/f'{sample}_plotted_spots.csv.gz',compression='gzip')
 # These figures use exactly the same frozen assignment table as the validated composition.
 for filename in [f'{organ}_{target}_04_composition.tiff',f'{organ}_{target}_05_microenvironment_legend.tiff']:
  source=record(V3/f'{organ}_{target}'/filename);p=dest/filename;shutil.copy2(source,p)
  check(filename+' byte identical to validated v3',sha(source)==sha(p));register(p)
 table=record(V1/f'{organ}_{target}_composition_plotted.csv');shutil.copy2(table,dest/table.name)
 check(organ+' composition sums 100',np.allclose(pd.read_csv(table,index_col=0).sum(),100,atol=1e-8))
shutil.copy2(V1/'zoom_ROI_manifest.csv',OUT/'zoom_ROI_manifest.csv')
pd.DataFrame(index).to_csv(OUT/'TIFF_index.csv',index=False,encoding='utf-8-sig')
pd.DataFrame(checks).to_csv(OUT/'validation_checks.csv',index=False)
pd.DataFrame(summary).to_csv(OUT/'sample_summary.csv',index=False)
(OUT/'source_hashes.json').write_text(json.dumps(provenance,indent=2))
(OUT/'README.txt').write_text('SGS10 蓝底参考版独立 TIFF 小图 v4\n\n本版替代 v3 的 H&E 组织图，用于客户自行拼图。肝 M05、脑 M04 各三只小鼠，每只提供全景、带 ROI 框全景、对应放大图；各器官另附组成图和独立图例，共 22 张 TIFF，全部 600 dpi、LZW 无损压缩。\n\n“蓝底”沿用客户五样本参考图的展示方式：不叠加 H&E，使用原始微环境配色及冻结 spot 标签，不是把底图染蓝。肝 M05 为红色；脑沿用原脑图配色，M04 为蓝色，并非强制整张图蓝色。组织内部无旧星号或箭头，无新增病理边界、比例尺及 A/B/C 编号。\n\n微环境来源：taskB/spatial_domain_RCTD100/<organ>_spot_microenv_assignment.csv；坐标来源：已有 RCTD 导出的坐标。逐点标签和坐标与已验证 v1 数据完全一致。脑方向与原五样本图一致（旋转180度）。按真实点位重绘，点径按独立面板可读性设置，非对旧截图放大。\n\n放大框沿用同一冻结分配表下 v1 的目标富集展示区域，替代旧 H&E 像素裁剪框。全景框与放大区域严格对应，不代表病理确认的肿瘤边界。\n\n组成图及数值保持不变：按整片目标微环境的平均 RCTD 权重统计，不限于放大框；每器官三只独立小鼠，描述性结果，无新增推断性检验。原始全景和带框全景二选一使用。\n',encoding='utf-8')
files=sorted(OUT.glob('*/*.tiff'));check('22 independent TIFF files',len(files)==22)
canvas=Image.new('RGB',(1500,((len(files)+2)//3)*270),'#eeeeee');draw=ImageDraw.Draw(canvas)
for i,p in enumerate(files):
 with Image.open(p) as im:
  im=im.convert('RGB');im.thumbnail((480,230));x=(i%3)*500;y=(i//3)*270
  canvas.paste(im,(x+(500-im.width)//2,y+28));draw.text((x+8,y+7),p.stem,fill='black')
canvas.save(OUT/'preview_contact.jpg')
(OUT/'validation_summary.json').write_text(json.dumps({'checks':len(checks),'all_pass':all(c['pass'] for c in checks),'visual_QA':'pending'},indent=2))
print(str(OUT),len(files),'TIFFs',len(checks),'checks passed')
