from pathlib import Path
import json, hashlib, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch, ConnectionPatch
from scipy.spatial import cKDTree

ROOT=Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis')
OUT=ROOT/'client_niche_figures_20260908_v1'
OUT.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'pdf.fonttype':42,'ps.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white'})
SAMPLES={'liver':['Liver_3','Liver_8','gan56'],'brain':['Brain_3','Brain_8','nao56']}
TARGET={'liver':'M05','brain':'M04'}
NAMES={'liver':'Tumor-immune niche','brain':'Vascular / BBB niche'}
ANIMAL_COLORS=['#247BA0','#E58F44','#9264AA']
PALETTES={'liver':{'M01':'#AFC8DD','M02':'#B5D9D5','M03':'#79B5C1','M04':'#CDB89F','M05':'#159447','M06':'#E9B9A5','M07':'#B9AED2'},'brain':{'M01':'#DDB6B8','M02':'#BDD2B1','M03':'#CCB998','M04':'#159447','M05':'#F0C7A0','M06':'#BEB6D6'}}
checks=[]; sources={}; all_comp=[]; all_roi=[]; summary=[]
def check(name,ok,detail=''):
    checks.append({'check':name,'pass':bool(ok),'detail':str(detail)})
    if not ok: raise AssertionError(name+': '+str(detail))
def record(path):
    p=Path(path)
    if str(p) not in sources:
        sources[str(p)]={'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
    return p

def load_sample(organ,s,assignment):
    co=pd.read_csv(record(OUT/f'{s}_rctd_coords.csv')).set_index('barcode')
    ap=assignment[assignment['sample']==s].copy()
    ap['barcode']=ap['spot_id'].str.removeprefix(s+'__')
    ap=ap.set_index('barcode')
    w=pd.read_csv(record(ROOT/'spatial_results/rctd_full/per_sample'/f'{s}_weights_matrix.csv'),index_col=0)
    check(s+' unique IDs',co.index.is_unique and ap.index.is_unique and w.index.is_unique)
    check(s+' all assignments are tumor',ap['group'].eq('tumor').all())
    check(s+' complete coordinate/assignment/weight join',set(ap.index)==set(w.index)==set(co.index),f'{len(ap)}/{len(w)}/{len(co)}')
    check(s+' finite weights with only small negative solver values',np.isfinite(w.values).all() and (w.values>=-0.001).all())
    w_original=w.div(w.sum(axis=1),axis=0).copy()
    negative_count=int((w.values<0).sum())
    original_min=float(w.values.min())
    w=w.clip(lower=0)
    check(s+' positive total weights',(w.sum(axis=1)>0).all())
    w=w.div(w.sum(axis=1),axis=0)
    co=co.loc[ap.index]
    w=w.loc[ap.index]
    dat=ap.join(co[['x','y']])
    # Preserve accepted brain orientation by rotating coordinates and H&E together if later supplied.
    sx,sy=(55050,19906) if s in ['gan56','nao56'] else (56500,20434)
    if organ=='brain':
        dat['x']=sx-dat['x'];dat['y']=sy-dat['y']
    check(s+' finite coordinates',np.isfinite(dat[['x','y']].values).all())
    target=dat.microenv_id.eq(TARGET[organ])
    check(s+' target present',target.any(),int(target.sum()))
    mean=w.loc[target].mean(axis=0)*100
    original_mean=w_original.loc[ap.index].loc[target].mean(axis=0)*100
    delta=float(np.max(np.abs(mean-original_mean)))
    check(s+' clipping has negligible niche composition effect',delta<0.01,f'negative_entries={negative_count}; min={original_min}; max_mean_change_pp={delta}')
    check(s+' composition sums to 100',abs(mean.sum()-100)<1e-8,mean.sum())
    comp=pd.DataFrame({'organ':organ,'sample':s,'microenv_id':TARGET[organ],'celltype':mean.index,'mean_RCTD_weight_pct':mean.values,'n_target_spots':int(target.sum())})
    all_comp.append(comp)
    xy=dat[['x','y']].to_numpy(); txy=xy[target.to_numpy()]
    lo=xy.min(axis=0);hi=xy.max(axis=0)
    side=float(np.min(hi-lo)*0.28)
    # Deterministic densest target-centered square, identical relative width for each sample.
    candidates=np.unique(np.clip(txy,lo+side/2,hi-side/2),axis=0)
    tree=cKDTree(txy)
    scores=tree.query_ball_point(candidates,side/2,p=np.inf,return_length=True)
    center=candidates[np.argmax(scores)]
    bounds=np.r_[center-side/2,center+side/2]
    inroi=(xy[:,0]>=bounds[0])&(xy[:,0]<=bounds[2])&(xy[:,1]>=bounds[1])&(xy[:,1]<=bounds[3])
    nt=int(np.sum(inroi & target.to_numpy()))
    check(s+' ROI target count matches selection',nt==int(scores.max()),nt)
    check(s+' ROI contains context',int(inroi.sum())>nt,f'{nt}/{int(inroi.sum())}')
    all_roi.append({'organ':organ,'sample':s,'microenv_id':TARGET[organ],'x_min':bounds[0],'y_min':bounds[1],'x_max':bounds[2],'y_max':bounds[3],'target_spots_in_ROI':nt,'all_spots_in_ROI':int(inroi.sum()),'target_spots_whole_sample':int(target.sum()),'rotation_deg':180 if organ=='brain' else 0,'selection':'maximum target count in target-centered square; width=0.28*shorter tissue bounding-box span; ties lexicographic','coordinate_unit':'source coordinate units; physical calibration unavailable'})
    dat.to_csv(OUT/f'{s}_plotted_spots.csv.gz',compression='gzip')
    summary.append({'organ':organ,'sample':s,'n_spots':len(dat),'n_target_spots':int(target.sum()),'target_spot_fraction_pct':float(target.mean()*100)})
    return dat,comp,bounds,lo,hi

def spatial(ax,dat,organ,bounds=None):
    if bounds is not None:
        dat=dat.loc[dat.x.between(bounds[0],bounds[2]) & dat.y.between(bounds[1],bounds[3])]
    for m in [x for x in PALETTES[organ] if x!=TARGET[organ]]+[TARGET[organ]]:
        d=dat[dat.microenv_id==m]
        ax.scatter(d.x,d.y,s=((13 if m==TARGET[organ] else 10) if bounds is not None else (1.4 if m==TARGET[organ] else .65)),color=('#159447' if m==TARGET[organ] else '#B2C5D8'),linewidths=.18 if bounds is not None else 0,edgecolors='white' if bounds is not None else 'none',rasterized=True,zorder=3 if m==TARGET[organ] else 2,alpha=1 if m==TARGET[organ] else .65)
    ax.set_aspect('equal');ax.invert_yaxis();ax.set_xticks([]);ax.set_yticks([])
    if bounds is None:
        ax.margins(.035)
    else:
        ax.set_xlim(bounds[[0,2]]);ax.set_ylim(bounds[[3,1]])
    for sp in ax.spines.values(): sp.set_visible(False)

def broad_class(ct,organ):
    if organ=='liver':
        if ct.startswith('HSC_') or ct=='Stellate': return 'Hepatic stellate cells'
        if ct.startswith('Macrophage_') or ct in ['Kupffer','KC_Homeostatic']: return 'Macrophages / Kupffer cells'
        if ct.startswith('LSEC') or ct=='Endothelial': return 'Endothelial cells'
        if ct.startswith('Hep_'): return 'Hepatocytes'
        if ct in ['cDC1','cDC2','pDC','Dendritic']: return 'Dendritic cells'
        if ct in ['B_cells','B_Plasma']: return 'B / plasma cells'
        if ct in ['CD4_T','CD8_T','T_cells','Treg','NKT']: return 'T / NKT cells'
        if ct=='NK_cells': return 'NK cells'
        if ct=='LCM': return 'LCM (reference label)'
    else:
        if ct.startswith(('Hpc_','Excitatory_','Inhibitory_')) or ct=='Neurons': return 'Neurons'
        if ct.startswith('Astrocyte_'): return 'Astrocytes'
        if ct.startswith('Microglia_'): return 'Microglia'
        if ct=='Oligo_Mature': return 'Oligodendrocytes'
        if ct=='OPCs': return 'Oligodendrocyte precursors'
        if ct=='Endothelial': return 'Endothelial cells'
        if ct=='VLMC': return 'Vascular / leptomeningeal cells'
        if ct=='CPECs': return 'Choroid plexus epithelial cells'
        if ct=='Ependymal': return 'Ependymal cells'
    return ct.replace('_',' ')

for organ,samples in SAMPLES.items():
    target=TARGET[organ]
    a=pd.read_csv(record(ROOT/'taskB/spatial_domain_RCTD100'/f'{organ}_spot_microenv_assignment.csv'))
    naming=pd.read_csv(record(ROOT/'taskB/spatial_domain_RCTD100'/f'{organ}_microenv_naming.csv'))
    check(organ+' target definition',set(naming.loc[naming.microenv_id==target,'domain'])==({'D_new_01'} if organ=='liver' else {'D_new_10'}))
    loaded=[load_sample(organ,s,a) for s in samples]
    comp=pd.concat([d[1] for d in loaded],ignore_index=True)
    comp['cell_class']=comp.celltype.map(lambda x:broad_class(x,organ))
    comp[['celltype','cell_class']].drop_duplicates().to_csv(OUT/f'{organ}_reference_to_cell_class.csv',index=False)
    broad=comp.groupby(['sample','cell_class'],sort=True)['mean_RCTD_weight_pct'].sum().reset_index()
    broad.to_csv(OUT/f'{organ}_niche_composition_all_classes.csv',index=False)
    pivot=broad.pivot(index='cell_class',columns='sample',values='mean_RCTD_weight_pct').reindex(columns=samples)
    check(organ+' broad class sums equal subtype sums',np.allclose(pivot.sum(),100,atol=1e-8))
    check(organ+' common reference cell types',not pivot.isna().any().any())
    top=pivot.mean(axis=1).sort_values(ascending=False,kind='stable').head(6).index.tolist()
    plotted=pivot.loc[top].copy()
    plotted.loc['Other cell classes']=pivot.drop(index=top).sum()
    plotted.index.name='celltype'
    plotted.to_csv(OUT/f'{organ}_{target}_composition_plotted.csv')
    check(organ+' plotted composition complete',np.allclose(plotted.sum(),100,atol=1e-8))
    fig=plt.figure(figsize=(15.6,12.4))
    gs=fig.add_gridspec(3,3,left=.10,right=.97,bottom=.085,top=.84,hspace=.40,wspace=.19,height_ratios=[.66,.90,1.24])
    fig.text(.045,.958,f'{organ.capitalize()} | {target} {NAMES[organ]}',fontsize=23,weight='bold',color='#17324D')
    fig.text(.045,.927,'Three tumor-bearing mice | Whole-section spatial maps, matched close-ups and niche composition',fontsize=11.5,color='#556575')
    fig.text(.045,.89,'A',fontsize=19,weight='bold');fig.text(.073,.893,'Whole-section spatial distribution',fontsize=13,weight='bold')
    for j,(s,(dat,_,roi,lo,hi)) in enumerate(zip(samples,loaded)):
        ax=fig.add_subplot(gs[0,j]);spatial(ax,dat,organ)
        ax.set_title(s,fontsize=13.5,weight='bold',color=ANIMAL_COLORS[j],pad=10)
        ax.add_patch(Rectangle(roi[:2],roi[2]-roi[0],roi[3]-roi[1],fill=False,edgecolor='#172B40',linewidth=1.3,zorder=5))
        ax.text(.01,-.085,f'{int(dat.microenv_id.eq(target).sum()):,} {target} spots',transform=ax.transAxes,fontsize=9,color='#17324D',bbox={'facecolor':'white','alpha':.9,'edgecolor':'none','pad':2})
        az=fig.add_subplot(gs[1,j]);spatial(az,dat,organ,roi)
        az.set_title(f'{s} | boxed region',fontsize=11.5,pad=9)
        for sp in az.spines.values():sp.set_visible(True);sp.set_color('#81909D');sp.set_linewidth(.7)
    pos=gs[1,0].get_position(fig)
    fig.text(.045,pos.y1+.036,'B',fontsize=19,weight='bold');fig.text(.073,pos.y1+.039,'Matched close-ups',fontsize=13,weight='bold')
    # Figure-level key below the spatial panels, without implying tumor pathology zones.
    fig.legend(handles=[Patch(facecolor='#159447',label=f'{target} target niche'),Patch(facecolor='#B2C5D8',label='Other microenvironments')],loc='upper center',bbox_to_anchor=(.54,pos.y0-.01),ncol=len(PALETTES[organ]),frameon=False,fontsize=9.5,handlelength=1.1,columnspacing=1.5)
    ax=fig.add_subplot(gs[2,:]);y=np.arange(len(plotted));barh=.22
    for j,s in enumerate(samples):
        values=plotted[s].to_numpy()
        ax.barh(y+(j-1)*barh,values,height=barh*.90,color=ANIMAL_COLORS[j],label=s)
        for yy,v in zip(y+(j-1)*barh,values):ax.text(v+.22,yy,f'{v:.1f}',fontsize=8.8,va='center',color='#273747')
    ax.set_yticks(y,[x.replace('_',' ') for x in plotted.index],fontsize=10.5)
    ax.invert_yaxis();ax.set_xlim(0,max(plotted.max().max()*1.16,10))
    ax.set_xlabel(f'Mean RCTD weight within the entire {target} niche (%)',fontsize=11)
    ax.xaxis.grid(True,alpha=.17);ax.set_axisbelow(True);ax.tick_params(axis='y',length=0)
    ax.spines['left'].set_visible(False);ax.spines['bottom'].set_color('#A0ABB5')
    ax.set_title(f'C   {target} cellular composition',loc='left',fontsize=13,weight='bold',pad=32)
    ax.legend(loc='lower right',bbox_to_anchor=(1,1.03),ncol=3,frameon=False,fontsize=10)
    # Longer subtype labels require left room independent of spatial map spacing.
    ap=ax.get_position();ax.set_position([.235,ap.y0,.72,ap.height])
    fig.text(.045,.035,'Composition: six leading cell classes (equal weight per mouse for ranking) + all remaining types. Each color is one biological sample.',fontsize=9,color='#556575')
    fig.text(.045,.017,'Spatial maps use the existing RCTD100 assignments. Green marks the target niche; close-ups show real mapped spots. Descriptive results, n = 3 mice per organ.',fontsize=9,color='#556575')
    stem=f'SGS10_{organ}_{target}_three_tumor_mice_v1'
    fig.savefig(OUT/f'{stem}.pdf',dpi=350)
    fig.savefig(OUT/f'{stem}.png',dpi=250)
    plt.close(fig)
    print(stem,'SAVED',flush=True)

pd.concat(all_comp,ignore_index=True).to_csv(OUT/'niche_composition_all_celltypes.csv',index=False)
pd.DataFrame(all_roi).to_csv(OUT/'zoom_ROI_manifest.csv',index=False)
pd.DataFrame(summary).to_csv(OUT/'sample_niche_summary.csv',index=False)
pd.DataFrame(checks).to_csv(OUT/'validation_checks.csv',index=False)
(OUT/'source_manifest.json').write_text(json.dumps(sources,indent=2))
(OUT/'build_summary.json').write_text(json.dumps({'status':'NUMERICAL_CHECKS_PASS_VISUAL_QA_PENDING','checks_passed':len(checks),'samples':summary,'python':sys.executable,'physical_scale':'No micrometre scale asserted; source-coordinate calibration has not been established.','HE':'Not used in spatial-only version; original external H&E directory unavailable for all six samples.'},indent=2))
print('CHECKS',len(checks),'PASS',flush=True)
