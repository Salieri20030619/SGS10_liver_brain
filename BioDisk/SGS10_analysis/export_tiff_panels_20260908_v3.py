from pathlib import Path
import hashlib, json, shutil, csv
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

ROOT = Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis')
SRC = ROOT/'client_niche_existing_HE_20260908_v2'
OUT = ROOT/'client_niche_tiff_panels_20260908_v3'
OUT.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'savefig.facecolor':'white'})
COLORS=['#247BA0','#E58F44','#9264AA']
LABELS={'liver':['M01 Hepatocyte zonation','M02 Vascular / perisinusoidal','M03 Sinusoidal endothelium','M04 Quiescent stellate','M05 Tumor-immune niche','M06 B / plasma','M07 Kupffer'], 'brain':['M01 Hippocampus','M02 Excitatory cortex','M03 Inhibitory cortex','M04 Vascular / BBB','M05 Glia','M06 General neurons']}
records=[]
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save_image(im,p,dpi=300):
    im.save(p,compression='tiff_lzw',dpi=(dpi,dpi))
    with Image.open(p) as check:
        check.load()
        assert check.format=='TIFF' and np.array_equal(np.asarray(im),np.asarray(check))
        assert abs(check.info['dpi'][0]-dpi)<.01
    records.append({'file':str(p.relative_to(OUT)),'width_px':im.width,'height_px':im.height,'dpi':dpi,'sha256':sha(p)})
def save_plot(fig,p):
    fig.savefig(p,dpi=600,bbox_inches='tight',pad_inches=.14,pil_kwargs={'compression':'tiff_lzw'})
    with Image.open(p) as im:
        im.load()
        assert im.format=='TIFF' and abs(im.info['dpi'][0]-600)<.01
        records.append({'file':str(p.relative_to(OUT)),'width_px':im.width,'height_px':im.height,'dpi':600,'sha256':sha(p)})
    plt.close(fig)
manifest=pd.read_csv(SRC/'crop_manifest.csv')
for organ, group in manifest.groupby('organ',sort=False):
    target=group.iloc[0]['target']; dest=OUT/f'{organ}_{target}';dest.mkdir(exist_ok=True)
    for r in group.to_dict('records'):
        s=r['sample']; source=Path(r['source'])
        assert sha(source)==r['source_sha256']
        im=Image.open(source).convert('RGB')
        whole=im.crop(tuple(int(r['whole_panel_'+k]) for k in ['left','top','right','bottom']))
        zoom=im.crop(tuple(int(r['zoom_'+k]) for k in ['left','top','right','bottom']))
        save_image(whole,dest/f'{s}_01_whole_native.tiff')
        save_image(zoom,dest/f'{s}_03_zoom_native.tiff')
        boxed=whole.copy();draw=ImageDraw.Draw(boxed)
        x=int(r['zoom_left']-r['whole_panel_left']);y=int(r['zoom_top']-r['whole_panel_top'])
        draw.rectangle((x-2,y-2,x+361,y+361),outline='white',width=6)
        draw.rectangle((x,y,x+359,y+359),outline='#182F44',width=2)
        save_image(boxed,dest/f'{s}_02_whole_with_ROI.tiff')
    samples=group['sample'].tolist()
    table=SRC/f'{organ}_{target}_composition_plotted.csv'
    shutil.copy2(table,dest/table.name)
    data=pd.read_csv(table,index_col=0).reindex(columns=samples)
    assert np.allclose(data.sum(),100,atol=1e-8)
    fig,ax=plt.subplots(figsize=(8.2,4.6));ys=np.arange(len(data));bh=.22
    for j,s in enumerate(samples):
        vals=data[s].to_numpy();ax.barh(ys+(j-1)*bh,vals,height=bh*.9,color=COLORS[j],label=s)
        for y,v in zip(ys+(j-1)*bh,vals):ax.text(v+.22,y,f'{v:.1f}',fontsize=9,va='center')
    ax.set_yticks(ys,data.index);ax.invert_yaxis();ax.set_xlim(0,max(data.max().max()*1.16,10))
    ax.set_xlabel(f'Mean RCTD weight within the entire {target} niche (%)')
    ax.tick_params(axis='y',length=0);ax.xaxis.grid(True,alpha=.17);ax.set_axisbelow(True)
    for n in ['left','top','right']:ax.spines[n].set_visible(False)
    ax.legend(loc='lower right',bbox_to_anchor=(1,1.02),ncol=3,frameon=False)
    save_plot(fig,dest/f'{organ}_{target}_04_composition.tiff')
    # The palette for every domain is recorded as its target in the source figure;
    # recover the exact complete palette using the original builder's pure helper.
    import ast
    tree=ast.parse((SRC/'build_existing_he_trial.py').read_text())
    fn=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='palette')
    from scipy.ndimage import label,find_objects
    ns={'np':np,'label':label,'find_objects':find_objects}
    exec(compile(ast.Module(body=[fn],type_ignores=[]),'<original palette helper>','exec'),ns)
    pal=ns['palette'](np.asarray(Image.open(group.iloc[0]['source']).convert('RGB')))
    assert len(pal)==len(LABELS[organ])
    fig,ax=plt.subplots(figsize=(4.5,2.3));ax.axis('off')
    ax.legend(handles=[Patch(facecolor=rgb/255,label=lab) for rgb,lab in zip(pal,LABELS[organ])],loc='center left',frameon=False)
    save_plot(fig,dest/f'{organ}_{target}_05_microenvironment_legend.tiff')
shutil.copy2(SRC/'crop_manifest.csv',OUT/'crop_manifest.csv')
with (OUT/'TIFF_index.csv').open('w',newline='',encoding='utf-8-sig') as f:
    writer=csv.DictWriter(f,fieldnames=records[0].keys());writer.writeheader();writer.writerows(records)
(OUT/'README.txt').write_text('SGS10 独立 TIFF 面板：肝 M05、脑 M04\n\n每个器官：三只小鼠，各有原始全景、带 ROI 框全景、对应放大裁片；另附细胞组成图和独立微环境图例。按文件名自行拼版，无新增 A/B/C 编号或总标题。\n全景与放大图以原生像素无损保存，300 dpi；放大裁片为 360×360 像素（300 dpi 下 3.048 cm），放大排版不会增加真实细节。组成图与图例重新渲染为 600 dpi。全部 TIFF 使用 LZW 无损压缩，保留源图方向、颜色和组织内部标注。\n每个器官包含 11 个 TIFF，共 22 个；原始全景和带框全景二选一使用。表中 DPI 是打印尺寸元数据，TIFF_index.csv 记录实际像素。\n组成图沿用已验证 v2 数据，表示整片目标微环境的平均 RCTD 权重；每器官三只独立肿瘤组小鼠，描述性结果。\n',encoding='utf-8')
(OUT/'validation.json').write_text(json.dumps({'tiff_count':len(records),'source_hashes_verified':True,'native_tiff_pixel_equality':True,'composition_source':'unchanged v2 CSV','dpi_verified':True},indent=2))
assert len(records)==22
print(json.dumps({'output':str(OUT),'tiff_count':len(records)}))
