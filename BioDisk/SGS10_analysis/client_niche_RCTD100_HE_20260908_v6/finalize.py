from pathlib import Path
import json,shutil,zipfile,hashlib
import pandas as pd
from PIL import Image,ImageDraw
root=Path(__file__).resolve().parent;out=root/'delivery'
old=Path.home()/'Downloads/SGS10_HE_rebuild_20260908_v5';v4=old.parent/'SGS10_blue_TIFF_20260908_v4'
(root/'reference_v5').mkdir(exist_ok=True)
for p in old.glob('*_HE_mapped_spots.csv.gz'):shutil.copy2(p,root/'reference_v5'/p.name)
for organ,target in [('liver','M05'),('brain','M04')]:
 dest=root/'reference_v4'/f'{organ}_{target}';dest.mkdir(parents=True,exist_ok=True)
 for name in [f'{organ}_{target}_composition_plotted.csv',f'{organ}_{target}_05_microenvironment_legend.tiff']:shutil.copy2(v4/f'{organ}_{target}'/name,dest/name)
audit=pd.read_csv(root/'full_RCTD_audit.csv')
readme='''SGS10 原始 HE + 全量 RCTD 历史叠图流程重绘 v6

已完成六个肿瘤组样本：Liver_3、Liver_8、gan56、Brain_3、Brain_8、nao56。
主图共 22 张独立 TIFF：每样本全景、带框全景、对应局部放大，另有每器官组成图和图例。HE_only 另附 12 张对应纯 HE 图。全部为 600 dpi、LZW 无损压缩；全景宽 3900 像素、放大图 2100×2100 像素。使用全景或带框全景之一拼版。

本版从服务器 spatial_results/rctd_full/per_sample 的六个 RCTD RDS 对象重新读取全部有效点，并逐一核验同名 weights_matrix.csv。全量运行 SUBSAMPLE_RATE=1.0；全量指纳入所有通过既定 mouse UMI>=100 过滤并进入 RCTD 结果的点，不是把低 UMI 点补回。六样本合计 275839 点，无绘图抽样、无条码连接损失。

重建严格复用历史 scripts/he_overlay_RCTD100.R 的 smooth_weights 与 assign_domain：每点取自身及 15 个邻居的权重均值，重新归一化，按与全量 domain profile 的余弦相似度确定 domain，再映射至微环境。原始 HE 保持原色，使用此前核验的 x_img=x*W/size_x、y_img=(size_y-y)*H/size_y 对齐方式。肝 M05 红色、脑 M04 蓝色。全景点面积 4 pt²，局部 40 pt²，alpha=0.80，所有微环境一致。

组成图使用本版平滑分类得到的整片目标微环境点集，汇总这些点的原始 RCTD 权重（不是用平滑权重充当组成），保持之前的类别映射和顺序。求解器微小负值截为 0 后逐点归一化，最大组成差异仅 0.00000150 个百分点。每器官 3 只独立小鼠分别展示，无新增推断性检验。局部框按既定规则选取目标点最密集的固定相对大小正方形，是展示区域，不是病理肿瘤边界。

版本核查结论：v5 已使用全量点集，不是 20% 抽样版；区别是 v5 沿用了未作该邻域平滑的聚类分配表。v6 恢复历史全量 HE 图的分类流程，因此颜色分布和目标点数发生变化，组成图及局部框已同步更新。完整新旧点数和变化数见 full_RCTD_audit.csv，不应将两种分配的组成图混用。

原始 HE 为客户重新提供的 aligned_HE.png。gan56 的深紫底色来自该原图；未做染色归一化。输出分辨率不代表增加了原始组织细节。
'''
assert int(audit.plotted_rows.sum())==275839
(out/'README.txt').write_text(readme,encoding='utf-8-sig')
shutil.copy2(root/'full_RCTD_audit.csv',out/'full_RCTD_audit.csv')
summary=json.loads((root/'validation_summary.json').read_text());summary['visual_QA']='PASS: all 22 main panels inspected on contact sheet; HE/point registration confirmed against six-source alignment preview; legends and composition labels readable.'
(root/'validation_summary.json').write_text(json.dumps(summary,indent=2))
# Large paired whole-map overview for convenient review.
canvas=Image.new('RGB',(1800,3*720),'white');draw=ImageDraw.Draw(canvas)
for row,(l,b) in enumerate(zip(['Liver_3','Liver_8','gan56'],['Brain_3','Brain_8','nao56'])):
 for col,(s,folder) in enumerate([(l,'liver_M05'),(b,'brain_M04')]):
  for j,(label,p) in enumerate([('v5: direct clustering',old/'delivery'/folder/f'{s}_01_HE_whole.tiff'),('v6: full RCTD, historical HE workflow',out/folder/f'{s}_01_HE_whole.tiff')]):
   with Image.open(p) as im:
    im=im.convert('RGB');im.thumbnail((880,310));x=col*900;y=row*720+j*355;canvas.paste(im,(x+(900-im.width)//2,y+34));draw.text((x+10,y+10),s+' | '+label,fill='black')
canvas.save(out/'v5_vs_v6_whole_maps.jpg',quality=94)
with zipfile.ZipFile(root.parent/'SGS10_RCTD100_HE_v6_TIFF.zip','w',compression=zipfile.ZIP_DEFLATED,compresslevel=3) as z:
 for p in out.rglob('*'):
  if p.is_file():z.write(p,p.relative_to(out))
print('PACKAGED',summary,flush=True)
