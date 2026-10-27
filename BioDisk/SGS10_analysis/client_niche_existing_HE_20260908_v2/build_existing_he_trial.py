from pathlib import Path
import json,hashlib,shutil
import numpy as np,pandas as pd
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle,Patch
from scipy.ndimage import label,find_objects,binary_dilation
ROOT=Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis')
OLD=ROOT/'client_niche_figures_20260908_v1'
OUT=ROOT/'client_niche_existing_HE_20260908_v2';OUT.mkdir(exist_ok=True)
NATIVE=OUT/'native_crops';NATIVE.mkdir(exist_ok=True)
SAMPLES={'liver':['Liver_3','Liver_8','gan56'],'brain':['Brain_3','Brain_8','nao56']}
TARGET={'liver':'M05','brain':'M04'}
BBOX={'Liver_3':(12,125,2060,837),'Liver_8':(12,125,2060,837),'gan56':(8,125,2002,870),'Brain_3':(250,114,1735,848),'Brain_8':(270,114,1880,833),'nao56':(35,114,1570,870)}
LABELS={'liver':['M01 Hepatocyte zonation','M02 Vascular / perisinusoidal','M03 Sinusoidal endothelium','M04 Quiescent stellate','M05 Tumor-immune niche','M06 B / plasma','M07 Kupffer'],'brain':['M01 Hippocampus','M02 Excitatory cortex','M03 Inhibitory cortex','M04 Vascular / BBB','M05 Glia','M06 General neurons']}
COLORS=['#247BA0','#E58F44','#9264AA']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'pdf.fonttype':42,'savefig.facecolor':'white'})
checks=[];records=[];provenance=[]
def ck(name,ok,detail=''):
 checks.append({'check':name,'pass':bool(ok),'detail':str(detail)})
 if not ok:raise AssertionError(name+': '+str(detail))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def palette(a):
 m=(a.max(2).astype(int)-a.min(2).astype(int)>45);m[:200]=False;m[:,:2060]=False
 lab,n=label(m);res=[]
 for i,b in enumerate(find_objects(lab),1):
  if b is None:continue
  y,x=b;ms=(lab[b]==i)
  if ms.sum()>40 and 7<x.stop-x.start<35 and 7<y.stop-y.start<35:res.append((y.start,np.median(a[b][ms],axis=0).astype(int)))
 return [x[1] for x in sorted(res,key=lambda v:v[0])]
for organ,samples in SAMPLES.items():
 target=TARGET[organ];target_index=int(target[1:])-1
 plotted_source=OLD/f'{organ}_{target}_composition_plotted.csv'
 plotted_copy=OUT/plotted_source.name;shutil.copy2(plotted_source,plotted_copy)
 ck(organ+' plotted composition byte-identical to validated v1',sha(plotted_source)==sha(plotted_copy))
 data=pd.read_csv(plotted_copy,index_col=0).reindex(columns=samples)
 ck(organ+' three samples with 100% composition',len(data.columns)==3 and np.allclose(data.sum(),100,atol=1e-8))
 fig=plt.figure(figsize=(15.6,15.0))
 fig.text(.045,.959,f'{organ.capitalize()} | {target} '+('Tumor-immune niche' if organ=='liver' else 'Vascular / BBB niche'),fontsize=23,weight='bold',color='#17324D')
 fig.text(.045,.929,'Existing H&E overlays | Three tumor-bearing mice | Whole sections, matched close-ups and cellular composition',fontsize=11.3,color='#536576')
 fig.text(.045,.891,'A   Whole-section H&E + microenvironment overlays',fontsize=13.5,weight='bold')
 fig.text(.045,.673,'B   Matched close-ups',fontsize=13.5,weight='bold')
 all_palette=[];axes_records=[]
 for j,s in enumerate(samples):
  src=ROOT/'taskB_he_overlay_RCTD100'/s/'D_microenvironment_on_HE.png'
  im=Image.open(src).convert('RGB');a=np.array(im);pal=palette(a);all_palette.append(pal)
  ck(s+' expected legend colors found',len(pal)==len(LABELS[organ]),len(pal))
  box=BBOX[s];whole=im.crop(box);w=np.array(whole)
  # Select the target-color-richest square inside the tissue panel using original pixel colors only.
  dist=np.linalg.norm(w.astype(float)-pal[target_index],axis=2)
  mask=(dist<=45).astype(np.int64)
  side=360
  ii=np.pad(mask,((1,0),(1,0))).cumsum(0).cumsum(1)
  count=ii[side:,side:]-ii[:-side,side:]-ii[side:,:-side]+ii[:-side,:-side]
  # Keep any intersected pre-existing black-text / yellow-star annotation fully inside the crop.
  black=(w.max(axis=2)<75)
  yellow=(w[:,:,0]>230)&(w[:,:,1]>185)&(w[:,:,2]<130)
  black=black|yellow
  bl,nbl=label(binary_dilation(black,structure=np.ones((7,31),bool)))
  annotations=[]
  gy,gx=np.indices(count.shape)
  valid=np.ones(count.shape,bool)
  for z,bb in enumerate(find_objects(bl),1):
   if bb is None:continue
   by,bx=bb
   if np.sum(black[bb] & (bl[bb]==z))<25:continue
   l=max(0,bx.start-5);r=min(w.shape[1],bx.stop+5);top=max(0,by.start-5);bot=min(w.shape[0],by.stop+5)
   annotations.append([l,top,r,bot])
   intersects=(gx<r)&(gx+side>l)&(gy<bot)&(gy+side>top)
   contains=(gx<=l)&(gx+side>=r)&(gy<=top)&(gy+side>=bot)
   valid &= (~intersects)|contains
  ck(s+' feasible crop without cutting annotations',valid.any())
  selected=np.where(valid,count,-1)
  yy,xx=np.unravel_index(np.argmax(selected),count.shape)
  cb=(int(xx),int(yy),int(xx+side),int(yy+side))
  ab=(box[0]+cb[0],box[1]+cb[1],box[0]+cb[2],box[1]+cb[3])
  crop=im.crop(ab)
  ck(s+' crop within whole panel',cb[0]>=0 and cb[1]>=0 and cb[2]<=whole.width and cb[3]<=whole.height)
  ck(s+' target-color pixels present',int(count[yy,xx])>100,int(count[yy,xx]))
  ck(s+' crop exactly matches source pixels',np.array_equal(np.array(crop),a[ab[1]:ab[3],ab[0]:ab[2]]))
  whole.save(NATIVE/f'{s}_whole_panel.png');crop.save(NATIVE/f'{s}_{target}_zoom_native.png')
  ck(s+' saved native crop unchanged',np.array_equal(np.array(Image.open(NATIVE/f'{s}_{target}_zoom_native.png')),np.array(crop)))
  left=.075+j*.303
  fig.text(left+.133,.862,s,fontsize=13.5,weight='bold',color=COLORS[j],ha='center')
  ax=fig.add_axes([left,.714,.267,.130])
  ax.imshow(whole,interpolation='none');ax.set_axis_off()
  for color,lw in [('white',2.8),('#182F44',1.2)]:ax.add_patch(Rectangle((cb[0]-.5,cb[1]-.5),side,side,fill=False,edgecolor=color,linewidth=lw))
  az=fig.add_axes([left+.041,.439,.185,.185*15.6/15.0])
  az.imshow(crop,interpolation='none');az.set_xticks([]);az.set_yticks([])
  for sp in az.spines.values():sp.set_linewidth(.8);sp.set_color('#748393')
  az.set_title(s+' | boxed region',fontsize=11,pad=9)
  axes_records.append((s,ax,az,whole.width,side))
  records.append({'organ':organ,'sample':s,'target':target,'source':str(src),'source_sha256':sha(src),'source_width_px':im.width,'source_height_px':im.height,'whole_panel_left':box[0],'whole_panel_top':box[1],'whole_panel_right':box[2],'whole_panel_bottom':box[3],'zoom_left':ab[0],'zoom_top':ab[1],'zoom_right':ab[2],'zoom_bottom':ab[3],'native_zoom_width_px':side,'native_zoom_height_px':side,'target_color_rgb':','.join(map(str,pal[target_index])),'target_color_pixel_count':int(count[yy,xx]),'annotation_bounds_in_whole_panel':json.dumps(annotations),'selection':'Maximum target-color pixel count (RGB distance <=45) in a 360x360 square subject to keeping intersected source annotations whole; row-major tie break. Display-only selection, not biological quantification.'})
  provenance.append({'path':str(src),'sha256':sha(src),'use':'Existing flattened H&E plus RCTD100 overlay; direct crop, colors and annotations preserved'})
 ck(organ+' palette consistent across three sources',all(np.array_equal(all_palette[0],p) for p in all_palette[1:]))
 legend=fig.legend(handles=[Patch(facecolor=rgb/255,label=lab) for rgb,lab in zip(all_palette[0],LABELS[organ])],loc='upper center',bbox_to_anchor=(.53,.415),ncol=4 if organ=='liver' else 3,frameon=False,fontsize=9.4,columnspacing=1.35,handlelength=1.05,labelspacing=.65)
 fig.text(.045,.354,f'C   {target} cellular composition',fontsize=13.5,weight='bold')
 ax=fig.add_axes([.235,.073,.72,.252]);ys=np.arange(len(data));bh=.22
 for j,s in enumerate(samples):
  vals=data[s].to_numpy()
  ax.barh(ys+(j-1)*bh,vals,height=bh*.90,color=COLORS[j],label=s)
  for y,v in zip(ys+(j-1)*bh,vals):ax.text(v+.22,y,f'{v:.1f}',fontsize=8.0,va='center',color='#273747')
 ax.set_yticks(ys,data.index,fontsize=10.5);ax.invert_yaxis()
 ax.set_xlim(0,max(data.max().max()*1.16,10));ax.set_xlabel(f'Mean RCTD weight within the entire {target} niche (%)',fontsize=10.8)
 ax.tick_params(axis='y',length=0);ax.xaxis.grid(True,alpha=.17);ax.set_axisbelow(True)
 for n in ['left','top','right']:ax.spines[n].set_visible(False)
 ax.spines['bottom'].set_color('#A0ABB5')
 ax.legend(loc='lower right',bbox_to_anchor=(1,1.03),ncol=3,frameon=False,fontsize=10)
 fig.text(.045,.029,'Close-ups are direct crops of the existing raster figures; original colors and annotations are retained. Boxes identify the exact displayed regions.',fontsize=9,color='#536576')
 fig.text(.045,.013,'Composition uses the whole target niche in each mouse, not the boxed region. Validated values retained from v1; descriptive results, n = 3 mice per organ.',fontsize=9,color='#536576')
 fig.canvas.draw()
 for s,a1,a2,width,side in axes_records:
  whole_ppi=width/(a1.get_position().width*15.6);zoom_ppi=side/(a2.get_position().width*15.6)
  rec=next(r for r in records if r['sample']==s)
  rec['display_enlargement_factor']=whole_ppi/zoom_ppi;rec['zoom_effective_native_ppi']=zoom_ppi
  a2.text(.5,-.038,f'{whole_ppi/zoom_ppi:.1f}x display enlargement',transform=a2.transAxes,fontsize=8.8,ha='center',va='top',color='#536576')
 stem=f'SGS10_{organ}_{target}_existing_HE_trial_v2'
 fig.savefig(OUT/f'{stem}.pdf',dpi=300);fig.savefig(OUT/f'{stem}.png',dpi=250)
 plt.close(fig);print(stem+' saved',flush=True)
pd.DataFrame(records).to_csv(OUT/'crop_manifest.csv',index=False)
pd.DataFrame(checks).to_csv(OUT/'validation_checks.csv',index=False)
(OUT/'source_manifest.json').write_text(json.dumps(provenance,indent=2))
shutil.copy2(OLD/'niche_composition_all_celltypes.csv',OUT/'niche_composition_all_celltypes.csv')
for organ in SAMPLES:
 for suffix in ['reference_to_cell_class.csv','niche_composition_all_classes.csv']:shutil.copy2(OLD/f'{organ}_{suffix}',OUT/f'{organ}_{suffix}')
print('CHECKS',len(checks),'PASS',flush=True)
