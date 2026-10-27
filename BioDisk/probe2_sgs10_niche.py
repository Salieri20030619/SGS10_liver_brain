from pathlib import Path
import pandas as pd
from PIL import Image
p=Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis')
for s, f in [('Liver_3','liver_samples/25090559_Liver_3_expression/Liver_3_filtered_feature_bc_matrix'),('Brain_3','brain_samples/Brain_3_filtered_feature_bc_matrix'),('gan56','gan56/gan_56_filtered_feature_bc_matrix'),('nao56','nao_56_filtered_feature_bc_matrix')]:
    x=pd.read_csv(p/'data'/f/'cell_locations.tsv.gz',sep='\t')
    print(s, x.shape, x.head(2).to_string(index=False), sep='\n')
for root in [p/'data',Path('/mnt/BioDisk/data'),Path('/mnt/BioHome/BAK')]:
    if root.exists():
        for d in root.iterdir():
            print('DIR',d)
print('Log dimensions:')
lines=(p/'taskB_he_overlay_RCTD100/run.log').read_text().splitlines()
for i,l in enumerate(lines):
    if 'size:' in l: print(lines[i-1],l)
