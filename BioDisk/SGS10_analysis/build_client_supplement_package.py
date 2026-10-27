import csv,json,shutil,subprocess,hashlib,zipfile,re,html
from pathlib import Path
from collections import Counter,defaultdict
from urllib.parse import quote
root=Path("/mnt/BioDisk/wangmingyuan/SGS10_analysis")
work=root/"client_supplement_package_20260908_work"
out=work/"SGS10_client_supplement_images_index"
qa=work/"qa"
out.mkdir(parents=True,exist_ok=True);qa.mkdir(exist_ok=True)
index=[];gaps=[];merges=[]
folders={"S1":"Tables/Table_S1_animal_sample_inventory.csv","S2":"Tables/Table_S2_ST_QC_summary.csv",
"S3":"Tables/Table_S3_clean_annotation_markers.csv","S4":"Tables/Table_S4_spatial_domain_definitions.csv",
"S5":"Tables/Table_S5_exported_significant_LR.csv"}
def read(p):
 with open(p,encoding="utf-8-sig",newline="") as f:return list(csv.DictReader(f))
def write_csv(rel,rows,fields=None):
 p=out/rel;p.parent.mkdir(parents=True,exist_ok=True)
 if fields is None:fields=list(rows[0]) if rows else []
 with p.open("w",encoding="utf-8-sig",newline="") as f:
  w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
 return p
def record(item,rel,title,source="",page="",note=""):
 index.append(dict(item=item,path=rel,title=title,source=source,page=str(page),note=note))
def copy(item,src,rel,title,note=""):
 s=root/src
 if not s.is_file() or s.stat().st_size==0:
  gaps.append(f"{item}: 未找到 {src}");return
 d=out/rel;d.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(s,d)
 record(item,rel,title,src,note=note)
def pdf_pages(p):
 r=subprocess.run(["pdfinfo",str(p)],capture_output=True,text=True,check=True)
 return int(re.search(r"^Pages:\s+(\d+)",r.stdout,re.M).group(1))
def text_pdf(p):
 return subprocess.run(["pdftotext",str(p),"-"],capture_output=True,check=True).stdout
def merge(item,sources,rel,title):
 valid=[];total=0;ledger=[];original_text=[]
 for src,label in sources:
  p=root/src
  if not p.is_file():gaps.append(f"{item}: 未找到 {src}");continue
  try:n=pdf_pages(p);txt=text_pdf(p)
  except Exception as e:gaps.append(f"{item}: PDF读取失败 {src}: {e}");continue
  if not txt.strip():
   imgs=subprocess.run(["pdfimages","-list",str(p)],capture_output=True,text=True).stdout
   if len(imgs.strip().splitlines())<=2:
    gaps.append(f"{item}: 无可提取文字或嵌入图像，未收录疑似空白PDF {src}");continue
  valid.append(p);original_text.append(txt)
  ledger.append((src,label,total+1,total+n));total+=n
 if not valid:gaps.append(f"{item}: 无可合并图片");return
 dest=out/rel;dest.parent.mkdir(parents=True,exist_ok=True)
 subprocess.run(["pdfunite",*[str(p) for p in valid],str(dest)],check=True,capture_output=True)
 assert pdf_pages(dest)==total
 norm=lambda b:re.sub(rb"\s+",b"",b)
 assert norm(text_pdf(dest))==norm(b"".join(original_text)),f"Text mismatch {rel}"
 record(item,rel,title,note=f"{total} 页；保留原图页面，未重绘")
 for src,label,start,end in ledger:
  record(item,rel,label,src,page=start,note=f"第 {start} 页" if start==end else f"第 {start}–{end} 页")
 merges.append({"path":rel,"pages":total,"sources":len(valid)})
 print("MERGED",rel,total,flush=True)
# Tables: only requested deliverables, no raw cell matrices or audit payload.
paired=read(root/"taskB/liver_brain_axis2/01_paired_full_metrics.csv")
rows=[]
for r in paired:
 rows.append({"mouse_id":r["mouse_id"],"group":r["group"],"liver_sample":r["liver_sample"],"brain_sample":r["brain_sample"],
 "liver_ST":"已有分析输出","brain_ST":"已有分析输出",
 "liver_HE_overlay":"已有图像" if (root/"taskB_he_overlay_RCTD100"/r["liver_sample"]/"D_microenvironment_on_HE.png").exists() else "未定位",
 "brain_HE_overlay":"已有图像" if (root/"taskB_he_overlay_RCTD100"/r["brain_sample"]/"D_microenvironment_on_HE.png").exists() else "未定位",
 "MRI":"当前项目未发现纳入MRI分析的记录","final_inclusion_record":"待核实最终实验台账"})
write_csv(folders["S1"],rows);record("Table S1",folders["S1"],"动物与配对样本清单","taskB/liver_brain_axis2/01_paired_full_metrics.csv",note="5只动物；MRI及最终纳入记录未补造")
summary=read(root/"spatial_results/summary/sample_summary.csv")
countrows=read(root/"TableS3_clean_20260907_v1/Table_S3_sample_celltype_counts.csv")
atlas=Counter()
for r in countrows:atlas[r["sample"]]+=int(r["n_cells"])
rctd={};domain_counts={}
for organ in ("liver","brain"):
 dc=read(root/f"taskB/spatial_domain_RCTD100/{organ}_per_sample_domain.csv");domain_counts[organ]=dc
 for r in dc:rctd[r["sample"]]=sum(int(v) for k,v in r.items() if k!="sample")
rows=[{"sample":r["sample"],"organ":r["organ"],"group":r["group"],
 "ST_summary_n_spots":r["n_spots"],"ST_matrix_n_features":r["n_genes"],
 "ST_median_genes_per_spot":r["median_nFeature"],"ST_median_UMI_per_spot":r["median_nCount"],
 "RCTD100_n_spots":rctd[r["sample"]],"clean_atlas_n_cells":atlas[r["sc_sample"]]} for r in summary]
write_csv(folders["S2"],rows);record("Table S2",folders["S2"],"ST、RCTD与清洗图谱计数摘要",
"spatial_results/summary/sample_summary.csv; taskB/spatial_domain_RCTD100/*_per_sample_domain.csv; TableS3_clean_20260907_v1/Table_S3_sample_celltype_counts.csv",
note="各列是不同分析阶段/对象的计数，不作未经确认的QC前后差值解释")
copy("Table S3","TableS3_clean_20260907_v1/Table_S3_clean_top30_markers_long.csv",folders["S3"],"清洗版肝34类＋脑23类注释marker","1710条亚型—基因记录；含效应量、检出比例和样本一致性，无旧P值")
defs=[]
for organ in ("liver","brain"):
 mapping=read(root/f"taskB/spatial_domain_RCTD100/{organ}_microenv_naming.csv")
 cts=read(root/f"D1_v2_RCTD_composition/microenv_key_celltypes_RCTD100/{organ}_microenv_celltypes_top5.csv")
 markers=read(root/f"D1_v2_RCTD_composition/microenv_markers_RCTD100/{organ}_microenv_markers_top20.csv")
 counts=Counter();names={}
 with (root/f"taskB/spatial_domain_RCTD100/{organ}_spot_microenv_assignment.csv").open() as f:
  for r in csv.DictReader(f):counts[r["microenv_id"]]+=1;names[r["microenv_id"]]=r["microenv_name"]
 for mid in sorted({r["microenv_id"] for r in mapping}):
  sub=[r for r in mapping if r["microenv_id"]==mid]
  top=[r for r in cts if r["microenv_id"]==mid]
  top.sort(key=lambda r:-float(r["mean_weight"]))
  mm=[r for r in markers if r["cluster"]==mid];mm.sort(key=lambda r:-float(r["avg_log2FC"]))
  assert all(int(r["n_spots"])==counts[mid] for r in top),f"composition count mismatch {organ} {mid}"
  defs.append({"organ":organ,"microenv_id":mid,"microenv_name":names[mid],"domains":"; ".join(r["domain"] for r in sub),
   "n_spots":counts[mid],"top5_celltypes":"; ".join(r["celltype"] for r in top[:5]),
   "top20_marker_genes":"; ".join(r["gene"] for r in mm[:20]),
   "historical_label":"; ".join(sorted({r["microenv_label"] for r in sub}))})
write_csv(folders["S4"],defs);record("Table S4",folders["S4"],"RCTD100生态位定义、组成与marker",
"taskB/spatial_domain_RCTD100/*; D1_v2_RCTD_composition/microenv_key_celltypes_RCTD100/*; D1_v2_RCTD_composition/microenv_markers_RCTD100/*",
note="肝M01–M07、脑M01–M06；肝脑同编号不代表同一生态位")
lr=[];expected=defaultdict(set);lr_sources=[]
for organ in ("liver","brain"):
 for group in ("tumor","blank","control"):
  src=f"taskB/cellchat/{organ}_{group}_lr_pairs.csv";lr_sources.append(src)
  for r in read(root/src):
   lr.append(dict(organ=organ,group=group,**r))
   if group=="tumor":expected[organ].add(r["pathway_name"])
fields=list(dict.fromkeys(k for r in lr for k in r))
for r in lr:
 for k in fields:r.setdefault(k,"")
write_csv(folders["S5"],lr,fields);record("Table S5",folders["S5"],"肝脑六组已导出的显著LR结果","; ".join(lr_sources),
note="保留全部已导出显著结果；不是无筛选的全部受检LR表；CellChat未随本次S3清洗重跑")
# Figure S1 has no independent experimental QC deliverable.
gaps.append("Figure S1：未找到独立的造模/MRI/病理QC成品图；Figure S3中的H&E背景仅作取材图像参考，不能替代实验QC。")
# Figure S2.
copy("Figure S2","figures/paper_v2/PDF/figS1_qc_violin.pdf","Figures/Figure_S2/01_per_sample_QC.pdf","按样本QC小提琴图","原文件名figS1；本包按客户新顺序归入Figure S2")
sample_order=[("liver","Con_liver_1"),("liver","Bla_liver_1"),("liver","Liver_3"),("liver","Liver_8"),("liver","gan56"),
 ("brain","Con_brain_1"),("brain","Bla_brain_1"),("brain","Brain_3"),("brain","Brain_8"),("brain","nao56")]
merge("Figure S2",[(f"spatial_results/plots/{s}/01_umap.pdf",f"{o} / {s}") for o,s in sample_order],
 "Figures/Figure_S2/02_per_sample_UMAP.pdf","10个样本逐样本UMAP")
for organ in ("liver","brain"):
 copy("Figure S2",f"taskB_spatial_figs_2026-04-20_v3/{organ}/00_merged_umap_by_sample.png",
      f"Figures/Figure_S2/03_{organ}_merged_UMAP_by_sample.png",f"{organ} 合并UMAP按样本着色",
      "现有ST模块版本；图像本身不构成已排除batch effect的结论")
# Figure S3: one microenvironment-on-HE image per sample, no redundant domain/hex variants.
for organ,sample in sample_order:
 copy("Figure S3",f"taskB_he_overlay_RCTD100/{sample}/D_microenvironment_on_HE.png",
      f"Figures/Figure_S3/{sample}_microenvironment_HE.png",f"{organ} / {sample} 微环境H&E叠加","RCTD100；三只肿瘤肝样本为Liver_3、Liver_8、gan56")
# Figures S4/S5.
for organ,item in (("liver","Figure S4"),("brain","Figure S5")):
 copy(item,f"figures/microenv_stacked/{organ}_per_sample.pdf",
      f"Figures/{item.replace(' ','_')}/{organ}_per_sample_microenvironment_proportion.pdf",
      f"{organ} 逐样本微环境比例","生成脚本replot_microenv_stacked.R读取RCTD100逐spot归属；不附重复PNG或按组汇总版本")
# Figure S6: original PDFs merged, including complete available pathway views.
merge("Figure S6",[(f"taskB/cellchat/{o}_{g}_aggregate_network.pdf",f"{o} / {g}") for o in ("liver","brain") for g in ("control","blank","tumor")],
 "Figures/Figure_S6/01_six_groups_networks.pdf","肝脑六组通信网络")
merge("Figure S6",[(f"taskB/cellchat/3groups_plots_v3/{o}_3groups_paper_grade.pdf",f"{o} 三组比较") for o in ("liver","brain")],
 "Figures/Figure_S6/02_three_group_comparisons.pdf","肝脑三组通信比较")
coverage=[]
for organ in ("liver","brain"):
 paths=sorted((root/f"taskB/cellchat/pathway_plots/{organ}").glob("*.pdf"))
 merge("Figure S6",[(str(p.relative_to(root)),p.stem) for p in paths],
  f"Figures/Figure_S6/03_{organ}_pathway_figures.pdf",f"{organ} 现有全部逐通路PDF")
 stems={p.stem for p in paths}
 missing=sorted(x for x in expected[organ] if x not in stems and not any(s.startswith(x+"_") for s in stems))
 coverage.append({"organ":organ,"tumor_LR_pathways":len(expected[organ]),"source_pathway_PDFs":len(paths),"unmatched_LR_pathways":missing})
 if missing:gaps.append(f"Figure S6 / {organ}：肿瘤LR导出中以下通路未找到同名PDF：{', '.join(missing)}")
# Chinese clickable index, with page ledger embedded.
titles={
"Table S1":("动物/样本信息","配对与分组已整理；MRI和最终纳入台账待补"),
"Table S2":("ST QC摘要","已有阶段计数与genes/UMI摘要；完整QC前后及剔除原因待补"),
"Table S3":("细胞注释marker","清洗版肝34＋脑23，已完成"),
"Table S4":("生态位/空间域定义","RCTD100定义、主要细胞组成与marker已整理"),
"Table S5":("配体—受体结果","全部已导出的显著LR；无筛选全受检表仍缺"),
"Figure S1":("造模/MRI/病理QC","缺独立实验QC图；不拿其他图充数"),
"Figure S2":("逐样本QC＋UMAP","已收录QC、10样本UMAP和合并按样本着色图"),
"Figure S3":("逐动物空间生态位","已收录肝脑全部10样本微环境图"),
"Figure S4":("肝M05/微环境定量","已收录逐样本全部7个微环境比例，包括M05"),
"Figure S5":("脑空间域定量","已收录逐样本全部6个微环境比例"),
"Figure S6":("完整CellChat/LR图片","已收录六组网络、三组比较和现有逐通路PDF；覆盖缺口见下方")}
order=["Figure S1","Figure S2","Figure S3","Figure S4","Figure S5","Figure S6","Table S1","Table S2","Table S3","Table S4","Table S5"]
e=html.escape
parts=['<!doctype html><html lang="zh-CN"><meta charset="UTF-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>SGS10 客户补充材料索引</title>',
'<style>body{font-family:Arial,"Microsoft YaHei",sans-serif;max-width:1100px;margin:36px auto;padding:0 24px;color:#243447;line-height:1.65}h1{color:#153d56}h2{border-bottom:2px solid #dbe7ed;padding-top:24px}a{color:#126b9a}table{width:100%;border-collapse:collapse}td,th{padding:9px 12px;border-bottom:1px solid #ddd;text-align:left;vertical-align:top}th{background:#eef5f8}small{color:#657987}img{max-width:100%;max-height:430px;object-fit:contain;background:#fafafa}details{margin:10px 0}.note{background:#fff7e7;padding:14px 18px}.nav{display:flex;gap:12px;flex-wrap:wrap}.card{border:1px solid #ddd;padding:14px;margin:12px 0}</style>',
'<h1>SGS10 客户补充图片与表格索引</h1><p>按客户截图的 Table S1–S5、Figure S1–S6 整理。请先完整解压，再打开本索引；点击文件即可查看，合并PDF另附页码定位。</p>',
'<p>本包仅含客户清单所需的图片、5张表和索引。PDF优先；没有PDF的图保留PNG。没有逐细胞数据、全部基因统计、RDS、代码、运行日志或审计包。</p>',
'<div class="note"><b>版本口径：</b>Table S3为清洗后肝34＋脑23类；空间生态位使用RCTD100版本。现有空间图与CellChat未在本次重跑。MRI未发现用于当前项目分析的记录，Figure S1不补造。</div>',
'<nav class="nav">'+''.join(f'<a href="#{k.replace(" ","_")}">{k}</a>' for k in order)+'</nav>',
'<h2>总览</h2><table><tr><th>客户项目</th><th>内容</th><th>状态</th></tr>']
for k in order:parts.append(f'<tr><td><a href="#{k.replace(" ","_")}">{k}</a></td><td>{e(titles[k][0])}</td><td>{e(titles[k][1])}</td></tr>')
parts.append('</table>')
for k in order:
 parts.append(f'<section id="{k.replace(" ","_")}"><h2>{k} · {e(titles[k][0])}</h2><p>{e(titles[k][1])}</p>')
 main=[r for r in index if r["item"]==k and not r["page"]]
 for r in main:
  href=quote(r["path"],safe="/")
  parts.append(f'<div class="card"><a href="{href}">{e(r["title"])}</a><br><small>{e(r["note"])}</small>')
  if r["path"].endswith(".png"):parts.append(f'<details><summary>展开图片预览</summary><a href="{href}"><img loading="lazy" src="{href}" alt="{e(r["title"])}"></a></details>')
  ledger=[q for q in index if q["path"]==r["path"] and q["page"]]
  if ledger:
   parts.append('<details><summary>展开PDF页码索引</summary><table><tr><th>内容</th><th>页码</th></tr>')
   for q in ledger:parts.append(f'<tr><td>{e(q["title"])}</td><td><a href="{href}#page={q["page"]}">{e(q["note"])}</a></td></tr>')
   parts.append('</table></details>')
  sources=[q for q in index if q["path"]==r["path"] and q["source"]]
  if sources:parts.append('<details><summary>查看项目来源</summary><small>'+ '<br>'.join(e(q["source"]) for q in sources)+'</small></details>')
  parts.append('</div>')
 parts.append('</section>')
parts.append('<h2>需补齐的材料</h2><ul>'+''.join('<li>'+e(g)+'</li>' for g in gaps)+'</ul>')
parts.append('<h2>使用口径</h2><p>5只动物：肿瘤3只、blank 1只、control 1只，配对肝脑共10个组织样本。细胞/spot数不代替生物学样本量。Table S3只报告描述性marker效应、检出比例和样本一致性；肝Basophil仅6个细胞。Table S5保留历史CellChat导出的prob/pval字段，不将它们解释为独立动物层面的组间证据。</p><p>原项目：/mnt/BioDisk/wangmingyuan/SGS10_analysis。现有图内的旧图号保留；本索引按客户的新编号对应。各合并PDF已核对页数与文字内容，不改变原图页面。</p></html>')
(out/"00_客户材料索引.html").write_text("\n".join(parts),encoding="utf-8")
readme="""# SGS10 客户补充图片与表格

完整解压后，双击 00_客户材料索引.html。
索引包含11项客户要求、对应文件、合并PDF页码、来源及缺口。

Figures/：客户需要的现有图片，PDF优先；H&E等仅有PNG的保留PNG。
Tables/：仅5张客户要求的表。
没有逐细胞数据、全部基因统计、RDS、脚本、日志或审计材料。

Figure S1缺独立造模/MRI/病理QC图，不补造。Table S1最终纳入台账及Table S2完整QC前后记录仍需补齐。
Table S3采用清洗后肝34类＋脑23类；空间域为RCTD100；CellChat为既有模块输出，未因S3清洗重跑。
Table S5仅是全部“已导出的显著LR”，不是无筛选的全部受检LR。

肝Basophil仅6个细胞；Table S3中该类仅保留描述性marker。各独立动物样本支持可在主表一致性字段中查看。
"""
(out/"01_阅读说明.md").write_text(readme,encoding="utf-8")
# QA outside the client package.
(work/"build_index.json").write_text(json.dumps({"index":index,"gaps":gaps,"merges":merges,"coverage":coverage},ensure_ascii=False,indent=2),encoding="utf-8")
qa_paths=[]
for m in merges:
 for pg in sorted({1,(m["pages"]+1)//2,m["pages"]}):
  prefix=qa/(Path(m["path"]).stem+f"_page{pg}")
  subprocess.run(["pdftoppm","-f",str(pg),"-l",str(pg),"-r","60","-singlefile","-png",str(out/m["path"]),str(prefix)],check=True,capture_output=True)
  qa_paths.append(str(prefix)+".png")
# Render direct selected figures as well.
for rel in ["Figures/Figure_S2/01_per_sample_QC.pdf","Figures/Figure_S4/liver_per_sample_microenvironment_proportion.pdf","Figures/Figure_S5/brain_per_sample_microenvironment_proportion.pdf"]:
 prefix=qa/Path(rel).stem
 subprocess.run(["pdftoppm","-f","1","-l","1","-r","75","-singlefile","-png",str(out/rel),str(prefix)],check=True,capture_output=True)
 qa_paths.append(str(prefix)+".png")
try:
 from PIL import Image,ImageOps,ImageDraw
 imgs=[]
 for path in qa_paths:
  im=Image.open(path).convert("RGB");im.thumbnail((480,350))
  tile=Image.new("RGB",(500,385),"white");tile.paste(im,((500-im.width)//2,25))
  ImageDraw.Draw(tile).text((8,5),Path(path).stem,fill="black")
  imgs.append(tile)
 for batch in range(0,len(imgs),9):
  sub=imgs[batch:batch+9];sheet=Image.new("RGB",(1500,385*((len(sub)+2)//3)),(225,230,235))
  for j,im in enumerate(sub):sheet.paste(im,((j%3)*500,(j//3)*385))
  sheet.save(qa/f"contact_sheet_{batch//9+1}.jpg",quality=90)
except ImportError:pass
files=[p for p in out.rglob("*") if p.is_file()]
assert len(list((out/"Tables").glob("*.csv")))==5
assert not any(p.suffix.lower() in (".rds",".r",".py",".log",".gz") for p in files)
assert all((out/r["path"]).exists() for r in index)
# Archive will be written after the visual QA review.
print(json.dumps({"output":str(out),"files":len(files),"image_files":sum(p.suffix.lower() in (".png",".pdf") for p in files),
 "bytes":sum(p.stat().st_size for p in files),"table_rows":{p.name:sum(1 for _ in p.open())-1 for p in (out/"Tables").glob("*.csv")},
 "gaps":gaps,"merges":merges,"coverage":coverage,"qa_contact_sheets":[str(p) for p in qa.glob("contact_sheet*.jpg")]},ensure_ascii=False,indent=2))
