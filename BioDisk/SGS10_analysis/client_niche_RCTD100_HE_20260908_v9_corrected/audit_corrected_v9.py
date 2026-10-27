from pathlib import Path
import hashlib,json
import numpy as np,pandas as pd
p=Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis');v6=p/'client_niche_RCTD100_HE_20260908_v6';o=p/'client_niche_RCTD100_HE_20260908_v9_corrected';o.mkdir(exist_ok=True)
samples={'liver':['Bla_liver_1','Con_liver_1','Liver_3','Liver_8','gan56'],'brain':['Bla_brain_1','Con_brain_1','Brain_3','Brain_8','nao56']}
props=[];means=[];checks=[];comparisons=[];changes=[]
for organ,ss in samples.items():
 a=pd.read_csv(p/'taskB/spatial_domain_RCTD100'/f'{organ}_spot_microenv_assignment.csv')
 mapping=pd.read_csv(p/'client_niche_figures_20260908_v1'/f'{organ}_reference_to_cell_class.csv').set_index('celltype').cell_class
 target='M05' if organ=='liver' else 'M04'
 oldtable=pd.read_csv(v6/'reference_v4'/f'{organ}_{target}'/f'{organ}_{target}_composition_plotted.csv',index_col=0)
 badtable=pd.read_csv(v6/'delivery'/f'{organ}_{target}'/f'{organ}_{target}_composition_plotted.csv',index_col=0)
 for s in ss:
  print(s,'AUDIT',flush=True)
  before=a[a['sample'].eq(s)].copy();before['barcode']=before.spot_id.str.removeprefix(s+'__');before=before.set_index('barcode')
  w=pd.read_csv(p/'spatial_results/rctd_full/per_sample'/f'{s}_weights_matrix.csv',index_col=0)
  assert w.index.is_unique and before.index.is_unique and set(w.index)==set(before.index)
  assert np.isfinite(w.values).all() and w.values.min()>-.001
  w=w.clip(lower=0);w=w.div(w.sum(axis=1),axis=0)
  if s in oldtable.columns:
   after=pd.read_csv(v6/'reference_v5'/f'{s}_HE_mapped_spots.csv.gz',index_col=0)
   assert set(before.index)==set(after.index) and before.microenv_id.sort_index().equals(after.microenv_id.sort_index())
  else:
   co=pd.read_csv(p/'Figure2_full_v8'/f'{s}_full_HE_assignments.csv.gz').set_index('barcode')
   assert set(co.index)==set(before.index)
   after=before.join(co[['x','y']]);after['raw_x']=after.x;after['raw_y']=after.y;after['y']=20434-after.y
  after.to_csv(o/f'{s}_corrected_HE_mapped_spots.csv.gz',compression='gzip')
  for m in sorted(a.microenv_id.unique()):
   b=before.microenv_id.eq(m);n=after.microenv_id.eq(m)
   props.append({'organ':organ,'sample':s,'microenv_id':m,'before_n_spots':int(b.sum()),'after_n_spots':int(n.sum()),'sample_total':len(w),'before_pct':100*b.mean(),'after_pct':100*n.mean(),'difference_pp':100*(n.mean()-b.mean())})
   if n.any():
    av=w.loc[after.index[n]].mean()
    means.extend({'organ':organ,'sample':s,'microenv_id':m,'celltype':ct,'mean_weight':val,'n_spots':int(n.sum())} for ct,val in av.items())
  hit=after.microenv_id.eq(target);mean=w.loc[after.index[hit]].mean()*100
  if s in oldtable.columns:
   cls=mean.groupby(mapping.reindex(mean.index)).sum()
   tab=pd.Series({k:float(cls[k]) if k!='Other cell classes' else float(cls.loc[~cls.index.isin(oldtable.index)].sum()) for k in oldtable.index})
   err=float(abs(tab-oldtable[s]).max());assert err<1e-10,(s,err)
   for c in tab.index:
    comparisons.append({'organ':organ,'sample':s,'cell_class':c,'original_full_pct':oldtable.loc[c,s],'corrected_recomputed_pct':tab[c],'difference_pp':tab[c]-oldtable.loc[c,s]})
    changes.append({'organ':organ,'sample':s,'cell_class':c,'original_full_pct':oldtable.loc[c,s],'incorrect_v6_v7_pct':badtable.loc[c,s],'corrected_v9_pct':tab[c]})
  checks.append({'sample':s,'full_weight_rows':len(w),'original_assignment_rows':len(before),'corrected_rows':len(after),'labels_changed':0,'target_spots':int(hit.sum())})
pd.DataFrame(props).to_csv(o/'microenvironment_proportion_before_after.csv',index=False)
pd.DataFrame(means).to_csv(o/'microenvironment_composition_per_sample.csv',index=False)
pd.DataFrame(comparisons).to_csv(o/'cell_composition_before_after.csv',index=False)
pd.DataFrame(changes).to_csv(o/'original_vs_v6_v7_vs_corrected.csv',index=False)
pd.DataFrame(checks).to_csv(o/'barcode_label_audit.csv',index=False)
result={'samples_audited':len(checks),'labels_changed':0,'max_microenvironment_proportion_change_pp':max(abs(x['difference_pp']) for x in props),'max_cell_composition_recalculation_difference_pp':max(abs(x['difference_pp']) for x in comparisons),'all_pass':True,'reference':'original 100% RCTD assignments used by v5; not v6/v7 reassignment'}
(o/'PROPORTION_VALIDATION.json').write_text(json.dumps(result,indent=2));print(json.dumps(result),flush=True)
