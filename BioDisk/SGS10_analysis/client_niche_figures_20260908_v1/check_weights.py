from pathlib import Path
import pandas as pd,numpy as np
p=Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis/spatial_results/rctd_full/per_sample')
for s in ['Liver_3','Liver_8','gan56','Brain_3','Brain_8','nao56']:
 w=pd.read_csv(p/f'{s}_weights_matrix.csv',index_col=0)
 print(s,'min',np.nanmin(w),'max',np.nanmax(w),'nonfinite',int((~np.isfinite(w)).values.sum()),'negative',int((w<0).values.sum()),'minrowsum',w.sum(axis=1).min())
