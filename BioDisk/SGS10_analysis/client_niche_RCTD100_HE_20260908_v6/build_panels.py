from pathlib import Path
import hashlib,json,shutil
import numpy as np
import pandas as pd
from PIL import Image,ImageDraw
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.spatial import cKDTree

ROOT=Path(__file__).resolve().parent
OLD=Path.home()/'Downloads/SGS10_HE_rebuild_20260908_v5'
V4=OLD.parent/'SGS10_blue_TIFF_20260908_v4'
if (ROOT/'reference_v5').exists():OLD=ROOT/'reference_v5'
if (ROOT/'reference_v4').exists():V4=ROOT/'reference_v4'
OUT=ROOT/'delivery';OUT.mkdir(exist_ok=True)
SAMPLES={'liver':['Liver_3','Liver_8','gan56'],'brain':['Brain_3','Brain_8','nao56']}
PAL={'liver':dict(zip(['M01','M02','M03','M04','M05','M06','M07'],['#2C7FB8','#7FCDBB','#41B6C4','#A1D99B','#E31A1C','#FDAE6B','#9E9AC8'])), 'brain':dict(zip(['M01','M02','M03','M04','M05','M06'],['#E31A1C','#33A02C','#A6761D','#1F78B4','#FF7F00','#984EA3']))}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'savefig.facecolor':'white'})
checks=[];index=[];sources={};rois=[];comparisons=[]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def record(p):sources[str(p)]=sha(p);return p
def ck(name,ok):
 checks.append({'check':name,'pass':bool(ok)})
 if not ok:raise AssertionError(name)
def register(p):
 with Image.open(p) as im:
  im.load();ck(p.name+' TIFF 600dpi LZW',im.format=='TIFF' and abs(im.info['dpi'][0]-600)<.01 and im.tag_v2[259]==5)
  index.append({'file':str(p.relative_to(OUT)),'width_px':im.width,'height_px':im.height,'dpi':600,'sha256':sha(p)})
audit=pd.read_csv(record(ROOT/'full_RCTD_audit.csv')).set_index('sample')
comp=pd.read_csv(record(ROOT/'composition_all_celltypes.csv'))
for organ,samples in SAMPLES.items():
 target='M05' if organ=='liver' else 'M04'
 dest=OUT/f'{organ}_{target}';dest.mkdir(exist_ok=True)
 pure=OUT/'HE_only'/organ;pure.mkdir(parents=True,exist_ok=True)
 for sample in samples:
  stem={'gan56':'gan_56','nao56':'nao_56'}.get(sample,sample)
  local=ROOT/'inputs'/f'{stem}_aligned_HE.png'
  if not local.exists():
   local.parent.mkdir(exist_ok=True);shutil.copy2(next((OLD/'inputs').rglob(stem+'_aligned_HE.png')),local)
  path=record(local);im=Image.open(path).convert('RGB');iw,ih=im.size
  sx,sy=(55050,19906) if sample in ['gan56','nao56'] else (56500,20434)
  ck(sample+' HE aspect',abs((iw/ih)/(sx/sy)-1)<.001)
  dat=pd.read_csv(record(ROOT/f'{sample}_full_HE_assignments.csv.gz')).set_index('barcode').sort_index()
  old=pd.read_csv(record(OLD/f'{sample}_HE_mapped_spots.csv.gz'),index_col=0).sort_index()
  ck(sample+' full barcodes preserved',dat.index.is_unique and set(dat.index)==set(old.index) and len(dat)==audit.loc[sample,'full_weight_rows']==audit.loc[sample,'plotted_rows'])
  ck(sample+' raw coordinates match independently exported full RCTD',np.array_equal(dat.x,old.raw_x) and np.array_equal(dat.y,old.raw_y))
  ck(sample+' previous labels audited',np.array_equal(dat.previous_microenv,old.microenv_id))
  dat['raw_x']=dat.x;dat['raw_y']=dat.y;dat['y']=sy-dat.raw_y
  dat['image_x_px']=dat.x*iw/sx;dat['image_y_px']=dat.y*ih/sy
  ck(sample+' point centers within HE',dat.image_x_px.between(0,iw).all() and dat.image_y_px.between(0,ih).all())
  # Reapply the prior deterministic target-density ROI rule to the restored HE labels.
  xy=dat[['x','y']].to_numpy();target_mask=dat.microenv_id.eq(target).to_numpy();txy=xy[target_mask]
  lo=xy.min(0);hi=xy.max(0);side=float(np.min(hi-lo)*.28)
  candidates=np.unique(np.clip(txy,lo+side/2,hi-side/2),axis=0)
  scores=cKDTree(txy).query_ball_point(candidates,side/2,p=np.inf,return_length=True)
  center=candidates[np.argmax(scores)];bounds=np.r_[center-side/2,center+side/2]
  inside=dat.x.between(bounds[0],bounds[2])&dat.y.between(bounds[1],bounds[3])
  ck(sample+' selected ROI count',(inside.to_numpy()&target_mask).sum()==scores.max())
  pxbox=(int(np.floor(bounds[0]/sx*iw)),int(np.floor(bounds[1]/sy*ih)),int(np.ceil(bounds[2]/sx*iw)),int(np.ceil(bounds[3]/sy*ih)))
  ck(sample+' crop inside image',0<=pxbox[0]<pxbox[2]<=iw and 0<=pxbox[1]<pxbox[3]<=ih)
  rois.append({'sample':sample,'x_min':bounds[0],'y_min':bounds[1],'x_max':bounds[2],'y_max':bounds[3],'n_ROI_spots':int(inside.sum()),'n_target_ROI_spots':int((inside.to_numpy()&target_mask).sum()),'n_total':len(dat),'n_target':int(target_mask.sum()),'source_HE':str(path),'HE_sha256':sha(path),'source_width_px':iw,'source_height_px':ih,'mapping':'image_x=raw_x*W/size_x; image_y=(size_y-raw_y)*H/size_y'})
  dat.to_csv(ROOT/f'{sample}_HE_mapped_spots.csv.gz',compression='gzip')
  for mode in ['01_HE_whole','02_HE_whole_with_ROI','03_HE_zoom','HE_only_whole','HE_only_zoom']:
   zoom=mode.endswith('zoom');pure_he=mode.startswith('HE_only')
   w=3.5 if zoom else 6.5;h=3.5 if zoom else 6.5*sy/sx
   fig=plt.figure(figsize=(w,h));ax=fig.add_axes([0,0,1,1])
   if zoom:
    bg=im.crop(pxbox);extent=[pxbox[0]/iw*sx,pxbox[2]/iw*sx,pxbox[3]/ih*sy,pxbox[1]/ih*sy];limits=bounds
   else:
    bg=im.resize((3900,round(3900*ih/iw)),Image.Resampling.LANCZOS);extent=[0,sx,sy,0];limits=[0,0,sx,sy]
   ax.imshow(bg,extent=extent,origin='upper',interpolation='none',zorder=0)
   if not pure_he:
    d=dat.loc[inside] if zoom else dat
    ck(sample+' '+mode+' all intended points',len(d)==(int(inside.sum()) if zoom else len(dat)))
    ax.scatter(d.x,d.y,c=d.microenv_id.map(PAL[organ]).tolist(),s=40 if zoom else 4,alpha=.80,linewidths=0,edgecolors='none',zorder=2)
   if mode=='02_HE_whole_with_ROI':
    for color,lw in [('white',2.4),('#172B40',1.1)]:ax.add_patch(Rectangle(bounds[:2],bounds[2]-bounds[0],bounds[3]-bounds[1],fill=False,edgecolor=color,linewidth=lw,zorder=3))
   ax.set_xlim(limits[0],limits[2]);ax.set_ylim(limits[3],limits[1]);ax.set_aspect('equal');ax.axis('off')
   p=(pure if pure_he else dest)/f'{sample}_{mode}.tiff'
   fig.savefig(p,dpi=600,pil_kwargs={'compression':'tiff_lzw'});plt.close(fig);register(p)
  print(sample,'saved full and zoom',len(dat),'spots',flush=True)
 # Use the same broad-class categories/order as the preceding delivery.
 classmap=pd.read_csv(record(ROOT/f'{organ}_reference_to_cell_class.csv')).set_index('celltype').cell_class
 oldcomp=pd.read_csv(record(V4/f'{organ}_{target}'/f'{organ}_{target}_composition_plotted.csv'),index_col=0)
 c=comp.loc[comp.organ.eq(organ)].copy();c['cell_class']=c.celltype.map(classmap)
 ck(organ+' complete class mapping',c.cell_class.notna().all())
 c.loc[~c.cell_class.isin(oldcomp.index),'cell_class']='Other cell classes'
 table=c.pivot_table(index='cell_class',columns='sample',values='mean_RCTD_weight_pct',aggfunc='sum').reindex(index=oldcomp.index,columns=samples).fillna(0)
 ck(organ+' composition sums 100',np.allclose(table.sum(),100,atol=1e-8))
 for s in samples:ck(s+' composition target matches map',comp.loc[comp['sample'].eq(s),'n_target_spots'].eq(audit.loc[s,'new_target_spots']).all())
 table.to_csv(dest/f'{organ}_{target}_composition_plotted.csv')
 delta=(table-oldcomp).stack().rename('change_percentage_points').reset_index();delta['organ']=organ;comparisons.append(delta)
 fig,ax=plt.subplots(figsize=(7.4,4.8));pos=np.arange(len(table));bw=.24
 for j,s in enumerate(samples):ax.barh(pos+(j-1)*bw,table[s],height=bw,color=['#247BA0','#E58F44','#9264AA'][j],label=s)
 ax.set_yticks(pos,table.index);ax.invert_yaxis();ax.set_xlabel('Mean RCTD weight within '+target+' (%)');ax.legend(frameon=False,ncol=3,loc='lower center',bbox_to_anchor=(.5,1.02));ax.spines[['top','right']].set_visible(False);fig.tight_layout()
 p=dest/f'{organ}_{target}_04_composition.tiff';fig.savefig(p,dpi=600,pil_kwargs={'compression':'tiff_lzw'});plt.close(fig);register(p)
 name=f'{organ}_{target}_05_microenvironment_legend.tiff';source=record(V4/f'{organ}_{target}'/name);shutil.copy2(source,dest/name);register(dest/name)
ck('34 TIFF files',len(index)==34)
pd.DataFrame(index).to_csv(OUT/'TIFF_index.csv',index=False,encoding='utf-8-sig')
pd.DataFrame(rois).to_csv(ROOT/'HE_ROI_and_mapping.csv',index=False,encoding='utf-8-sig')
pd.DataFrame(checks).to_csv(ROOT/'validation_checks.csv',index=False)
pd.concat(comparisons).to_csv(ROOT/'composition_change_vs_v5.csv',index=False)
(ROOT/'source_hashes.json').write_text(json.dumps(sources,indent=2))
files=sorted(OUT.glob('*/*.tiff'));canvas=Image.new('RGB',(1500,((len(files)+2)//3)*285),'#eeeeee');draw=ImageDraw.Draw(canvas)
for i,p in enumerate(files):
 with Image.open(p) as im:
  im=im.convert('RGB');im.thumbnail((480,245));x=(i%3)*500;y=(i//3)*285;canvas.paste(im,(x+(500-im.width)//2,y+28));draw.text((x+8,y+7),p.stem,fill='black')
canvas.save(OUT/'preview_contact.jpg')
(ROOT/'validation_summary.json').write_text(json.dumps({'checks':len(checks),'all_pass':all(c['pass'] for c in checks),'visual_QA':'pending','main_TIFFs':22,'HE_only_TIFFs':12},indent=2))
print('DONE',len(index),'TIFFs',len(checks),'checks',flush=True)
