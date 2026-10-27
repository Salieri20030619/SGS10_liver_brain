from pathlib import Path
import ast,json,hashlib
import numpy as np,pandas as pd
from PIL import Image,ImageColor
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch,ConnectionPatch,Rectangle
import matplotlib.patheffects as pe
ROOT=Path(__file__).resolve().parent;V6=ROOT.parent/'SGS10_RCTD100_HE_v6';D=ROOT/'delivery';PANELS=ROOT/'figure2_panels';PANELS.mkdir(exist_ok=True)
# Reuse the bounded continuous display only; no smoothing/reclassification is evaluated.
code=ast.parse((V6/'build_continuous_v7.py').read_text(encoding='utf-8'))
env={'np':np,'Image':Image}
exec(compile(ast.Module(body=[n for n in code.body if isinstance(n,ast.FunctionDef) and n.name=='render'],type_ignores=[]),'<render>','exec'),env)
render=env['render']
PAL={'liver':dict(zip(['M01','M02','M03','M04','M05','M06','M07'],['#2C7FB8','#7FCDBB','#41B6C4','#A1D99B','#E31A1C','#FDAE6B','#9E9AC8'])),'brain':dict(zip(['M01','M02','M03','M04','M05','M06'],['#E31A1C','#33A02C','#A6761D','#1F78B4','#FF7F00','#984EA3']))}
SAMPLES={'liver':['Bla_liver_1','Con_liver_1','Liver_3','Liver_8','gan56'],'brain':['Bla_brain_1','Con_brain_1','Brain_3','Brain_8','nao56']}
LABELS={'liver':['Hepatocyte zonation','Vascular / perisinusoidal','Sinusoidal endothelium','Quiescent stellate','Tumor-immune niche','B / plasma','Kupffer'],'brain':['Hippocampus','Excitatory cortex','Inhibitory cortex','Vascular / BBB','Glia','General neurons']}
raw=Path.home()/'Downloads/SGS10_HE_rebuild_20260908_v5/inputs'
for organ,ss in SAMPLES.items():
 for s in ss[:2]:
  dest=PANELS/f'{s}_whole.png'
  if dest.exists():continue
  dat=pd.read_csv(ROOT/f'{s}_corrected_HE_mapped_spots.csv.gz').sort_values('barcode')
  xy=dat[['x','y']].to_numpy();nn=cKDTree(xy).query(xy,k=2)[0][:,1];limit=1.5*np.quantile(nn,.95)
  visible=dat.drop_duplicates(['x','y'],keep='last');tree=cKDTree(visible[['x','y']]);colors=np.array([ImageColor.getrgb(PAL[organ][x]) for x in visible.microenv_id],dtype=np.uint8)
  image=Image.open(next(raw.rglob(s+'_aligned_HE.png'))).convert('RGB')
  im,_,_=render(image,tree,colors,[0,0,56500,20434],(2200,796),56500,20434,limit);im.save(dest)
  print(s,'baseline map ready',flush=True)
import sys
if '--prepare-only' in sys.argv:raise SystemExit(0)
props=pd.read_csv(ROOT/'microenvironment_proportion_before_after.csv');means=pd.read_csv(ROOT/'microenvironment_composition_per_sample.csv')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'ps.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
fig=plt.figure(figsize=(12,12.8),facecolor='white')
roi_table=pd.read_csv(D/'HE_ROI_and_mapping.csv').set_index('sample')
geometry_checks=[]
orientation_rows=[]
def text(x,y,s,size=10,weight='normal',**kw):fig.text(x,y,s,fontsize=size,fontweight=weight,**kw)
def box(x,y,w,h,color):fig.add_artist(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.004,rounding_size=0.008',facecolor=color,edgecolor='none',transform=fig.transFigure,zorder=-10))
box(.025,.511,.957,.457,'#fff1f1')
box(.025,.024,.957,.472,'#eff8ff')
text(.018,.984,'Figure 2. Organ-specific spatial ecosystems in liver and brain',18,'bold',va='top')
for organ,top,bottom,shade,letter in [('liver',.954,.745,'#fff1f1','A'),('brain',.482,.267,'#eff8ff','D')]:
 box(.025,bottom,.957,top-bottom+.014,shade)
 text(.012,(top+bottom)/2,organ.capitalize(),13,'bold',rotation=90,ha='center',va='center')
 text(.045,top,letter,17,'bold');text(.077,top,'Spatial ecosystems across blank, control and tumor-bearing mice',11,'bold')
 ss=SAMPLES[organ];target='M05' if organ=='liver' else 'M04';full_axes=[]
 for j,s in enumerate(ss):
  x=.075+j*.180;label=['Blank (Bla)','Control (Con)','Tumor-1','Tumor-2','Tumor-3'][j]
  text(x+.082,top-.020,label,10,'bold',ha='center');text(x+.082,top-.033,s,7.5,ha='center',color='#555555')
  path=PANELS/f'{s}_whole.png' if j<2 else D/f'{organ}_{target}'/f'{s}_01_HE_whole.tiff'
  im=Image.open(path).convert('RGB')
  rotate=s=='Bla_brain_1'
  if rotate:im=im.transpose(Image.Transpose.ROTATE_180)
  im.thumbnail((1500,700))
  sx,sy=(55050,19906) if s in ['gan56','nao56'] else (56500,20434)
  ax=fig.add_axes([x,top-.104,.165,.068]);ax.imshow(im,extent=[0,sx,sy,0]);ax.axis('off');full_axes.append(ax)
  orientation_rows.append({'sample':s,'rotation_from_v9_degrees':180 if rotate else 0,'mirror':False})
  if j>=2:
   row=roi_table.loc[s];bounds=np.array([row.x_min,row.y_min,row.x_max,row.y_max])
   displayed=np.array([sx-bounds[2],sy-bounds[3],sx-bounds[0],sy-bounds[1]]) if rotate else bounds.copy()
   assert np.allclose(displayed[2:]-displayed[:2],bounds[2:]-bounds[:2])
   if rotate:assert np.allclose([sx-displayed[2],sy-displayed[3],sx-displayed[0],sy-displayed[1]],bounds)
   im=Image.open(D/f'{organ}_{target}'/f'{s}_03_HE_zoom.tiff').convert('RGB')
   if rotate:im=im.transpose(Image.Transpose.ROTATE_180)
   im.thumbnail((1100,1100))
   az=fig.add_axes([x,top-.191,.165,.082]);az.imshow(im);az.axis('off')
   # Draw ROI in final figure units so downsampling cannot erase the stroke.
   patch=Rectangle(displayed[:2],displayed[2]-displayed[0],displayed[3]-displayed[1],fill=False,edgecolor='#151515',linewidth=1.35,zorder=10)
   patch.set_path_effects([pe.Stroke(linewidth=3.5,foreground='#FFF04A'),pe.Normal()]);ax.add_patch(patch)
   frame=Rectangle((0,0),1,1,transform=az.transAxes,fill=False,edgecolor='#151515',linewidth=1.05,clip_on=False,zorder=10)
   frame.set_path_effects([pe.Stroke(linewidth=2.8,foreground='#FFF04A'),pe.Normal()]);az.add_patch(frame)
   connection=ConnectionPatch(((displayed[0]+displayed[2])/2,displayed[3]),(.5,1),coordsA='data',coordsB='axes fraction',axesA=ax,axesB=az,color='#252525',lw=.8,zorder=8)
   connection.set_path_effects([pe.Stroke(linewidth=1.8,foreground='white'),pe.Normal()]);fig.add_artist(connection)
   text(x+.082,top-.203,target+' zoom',8,'bold',ha='center')
   geometry_checks.append({'sample':s,'ROI_size_preserved':True,'ROI_source_bounds':bounds.tolist(),'ROI_display_bounds':displayed.tolist(),'zoom_and_whole_rotation_match':True})
 for i,(m,label) in enumerate(zip(PAL[organ],LABELS[organ])):
  y=top-.118-i*.0125
  text(.088,y,'●',9,color=PAL[organ][m]);text(.104,y,m+'  '+label,8)
for organ,bottom,letter1,letter2 in [('liver',.516,'B','C'),('brain',.029,'E','F')]:
 shade='#fff7f7' if organ=='liver' else '#f5fbff';box(.025,bottom-.005,.957,.213,shade)
 text(.045,bottom+.193,letter1,17,'bold');text(.078,bottom+.193,'Cellular composition of '+organ+' ecosystems',11,'bold')
 target='M05' if organ=='liver' else 'M04';text(.678,bottom+.193,letter2,17,'bold');text(.710,bottom+.193,target+' composition in tumor-bearing mice',10,'bold')
 ss=SAMPLES[organ];pp=props[props.organ.eq(organ)].pivot(index='microenv_id',columns='sample',values='after_pct').reindex(index=list(PAL[organ]),columns=ss).fillna(0)
 ax=fig.add_axes([.074,bottom+.045,.198,.132]);acc=np.zeros(5)
 for m in pp.index:ax.bar(np.arange(5),pp.loc[m],bottom=acc,color=PAL[organ][m],width=.76,label=m);acc+=pp.loc[m].to_numpy()
 assert np.allclose(acc,100)
 ax.set_ylim(0,100);ax.set_ylabel('Proportion of spots (%)',fontsize=9);ax.set_xticks(range(5),['Bla','Con','T1','T2','T3'],fontsize=8);ax.tick_params(axis='y',labelsize=8);ax.legend(loc='center left',bbox_to_anchor=(1, .5),frameon=False,fontsize=7,handlelength=.8,labelspacing=.35)
 cm=pd.read_csv(V6/f'{organ}_reference_to_cell_class.csv').set_index('celltype').cell_class
 data=means[means.organ.eq(organ)].copy();data['class']=data.celltype.map(cm);assert data['class'].notna().all()
 if organ=='liver':
  broad={'Hepatocytes':'Hepatocytes','Endothelial cells':'Endothelial','Hepatic stellate cells':'Stellate','Macrophages / Kupffer cells':'Mac. / Kupffer','T / NKT cells':'T / NK cells','NK cells':'T / NK cells','B / plasma cells':'B / plasma'}
 else:broad={'Neurons':'Neurons','Astrocytes':'Astrocytes','Oligodendrocytes':'Oligodendrocytes','Microglia':'Microglia','Endothelial cells':'Endothelial','Pericytes':'Pericytes'}
 data['broad']=data['class'].map(broad).fillna('Other');cats=list(dict.fromkeys(broad.values()))+['Other']
 grouped=data.groupby(['sample','microenv_id','broad']).mean_weight.sum().unstack('broad').reindex(columns=cats,fill_value=0).fillna(0)
 profile=grouped.groupby('microenv_id').mean().reindex(list(PAL[organ]))
 assert np.allclose(profile.sum(axis=1),1)
 previous=pd.read_csv(ROOT/f'{organ}_bubble_equal_sample_profile.csv',index_col=0)
 assert np.allclose(profile.values,previous.values,rtol=0,atol=1e-14)
 sizes=props[props.organ.eq(organ)].groupby('microenv_id').after_n_spots.sum().reindex(profile.index)
 ab=fig.add_axes([.383,bottom+.045,.226,.132]);xs,ys=np.meshgrid(np.arange(len(cats)),np.arange(len(profile)))
 sizes_grid=np.repeat((sizes/80000*105+4).to_numpy(),len(cats))
 dots=ab.scatter(xs.ravel(),ys.ravel(),c=profile.to_numpy().ravel(),s=sizes_grid,cmap='YlOrRd',vmin=0,vmax=1,edgecolors='none')
 ab.set_xticks(range(len(cats)),cats,rotation=50,ha='right',fontsize=7);ab.set_yticks(range(len(profile)),profile.index,fontsize=8);ab.set_ylim(-.5,len(profile)-.5);ab.set_xlim(-.6,len(cats)-.4);ab.grid(axis='y',color='#eeeeee');ab.set_axisbelow(True)
 cb=fig.add_axes([.621,bottom+.116,.007,.055]);fig.colorbar(dots,cax=cb,ticks=[0,.5,1]);cb.tick_params(labelsize=6,length=2);cb.set_title('Mean\nweight',fontsize=6)
 for j,n in enumerate([10000,40000,80000]):ab.scatter([],[],s=n/80000*105+4,c='#555555',label=f'{n//1000}k')
 ab.legend(title='n spots',loc='upper left',bbox_to_anchor=(1.045,.47),frameon=False,fontsize=5.7,title_fontsize=6,handletextpad=.2,labelspacing=.5)
 comp=pd.read_csv(D/f'{organ}_{target}'/f'{organ}_{target}_composition_plotted.csv',index_col=0)
 ac=fig.add_axes([.805,bottom+.045,.160,.132]);pos=np.arange(len(comp));bw=.28
 for j,s in enumerate(ss[2:]):
  ac.barh(pos+(j-1)*bw,comp[s],height=bw*.88,color=['#247BA0','#E58F44','#9264AA'][j],label='Tumor-'+str(j+1))
  for i,v in enumerate(comp[s]):ac.text(v+.6,pos[i]+(j-1)*bw,f'{v:.1f}',fontsize=4.8,va='center')
 lab=[x.replace('Macrophages / Kupffer cells','Mac. / Kupffer').replace('Hepatic stellate cells','Stellate cells').replace('Other cell classes','Other').replace('Endothelial cells','Endothelial') for x in comp.index]
 ac.set_yticks(pos,lab,fontsize=7);ac.invert_yaxis();ac.set_xlim(0,105 if organ=='brain' else 40);ac.tick_params(axis='x',labelsize=7);ac.set_xlabel('Mean RCTD weight (%)',fontsize=8);ac.legend(frameon=False,fontsize=6.2,loc='lower center',bbox_to_anchor=(.5,1.01),ncol=3,borderaxespad=0,handlelength=1,columnspacing=.7)
 assert np.allclose(ax.get_position().bounds[1::2],ab.get_position().bounds[1::2])
 assert np.allclose(ax.get_position().bounds[1::2],ac.get_position().bounds[1::2])
text(.027,.004,'Full RCTD; original microenvironment assignments. Continuous colors are display interpolation; proportions use original spots and weights.',7,color='#555555')
fig.savefig(ROOT/'Figure2_layout_orientation_v10.pdf',dpi=300)
fig.savefig(ROOT/'Figure2_layout_orientation_v10.png',dpi=240)
fig.savefig(ROOT/'Figure2_layout_orientation_v10.tiff',dpi=600,pil_kwargs={'compression':'tiff_lzw'})
plt.close(fig)
import pypdfium2 as pdfium
pdf=pdfium.PdfDocument(ROOT/'Figure2_layout_orientation_v10.pdf');assert len(pdf)==1
page=pdf[0];page.render(scale=2).to_pil().save(ROOT/'Figure2_v10_pdf_QA.png')
(ROOT/'Figure2_v10_layout_validation.json').write_text(json.dumps({'orientation':orientation_rows,'roi_checks':geometry_checks,'statistic_axes_top_bottom_aligned':True,'original_bubble_profiles_unchanged':True,'analysis_inputs':'unchanged v9 original-full labels and proportions'},indent=2))
print('FIGURE2 V10 COMPLETE',flush=True)
