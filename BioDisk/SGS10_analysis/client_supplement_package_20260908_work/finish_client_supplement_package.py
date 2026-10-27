import json,csv,zipfile,hashlib,re,html
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
from urllib.parse import unquote
root=Path("/mnt/BioDisk/wangmingyuan/SGS10_analysis")
work=root/"client_supplement_package_20260908_work"
out=work/"SGS10_client_supplement_images_index"
audit=json.loads((work/"source_time_audit.json").read_text())
ledger=json.loads((work/"build_index.json").read_text())["index"]
cut=datetime(2026,4,13,tzinfo=ZoneInfo("Asia/Shanghai"))
rows=[]
for r in audit:
 mt=datetime.fromisoformat(r["mtime"]);assert mt>=cut,r["source"]
 pd=r["pdf_creation"]
 if pd:
  dt=datetime.strptime(" ".join(pd.split()[:5]),"%a %b %d %H:%M:%S %Y").replace(tzinfo=cut.tzinfo)
  assert dt>=cut,r["source"]
 related=sorted({x["item"] for x in ledger if r["source"] in [s.strip() for s in x["source"].split(";")]})
 if not related:
  if r["source"].startswith(("D1_v2_RCTD_composition","taskB/spatial_domain_RCTD100")):related=["Table S4 / RCTD100上游"]
  elif r["source"].startswith(("umap_paper_grade","dam_clean_v2")):related=["Table S3清洗图谱上游"]
 rows.append({"客户项目":"; ".join(related),"源文件":r["source"],"源文件修改时间":r["mtime"],
              "PDF内嵌创建时间":pd,"时间依据":"PDF创建日期及文件修改时间" if pd else "文件修改时间；无内嵌创建日期",
              "时间范围核对":"符合：2026-04-13当天及以后"})
with (out/"02_源文件时间索引.csv").open("w",encoding="utf-8-sig",newline="") as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
index=out/"00_客户材料索引.html"
s=index.read_text(encoding="utf-8")
insert='<div class="note"><b>有效时间范围：</b>按用户确认，采用2026年4月13日当天及以后。已核对177个源文件，均符合该时间范围；本次合并/打包时间不作为原图生成时间。PDF记录内嵌创建日期及源文件修改时间，PNG/CSV等无内嵌创建日期者明确列出修改时间作为证据。<a href="02_%E6%BA%90%E6%96%87%E4%BB%B6%E6%97%B6%E9%97%B4%E7%B4%A2%E5%BC%95.csv">打开源文件时间索引</a>。原始测序文件不受此规则限制，本包未收录原始测序数据。</div>'
s=s.replace('<nav class="nav">',insert+'<nav class="nav">')
index.write_text(s,encoding="utf-8")
with (out/"01_阅读说明.md").open("a",encoding="utf-8") as f:f.write("\n有效时间范围：2026-04-13当天及以后（用户已确认包含当天）。177个源文件均符合此范围。具体证据见02_源文件时间索引.csv；文件修改时间不冒充已独立证实的生成时间，合并/打包时间不用于判定原始材料日期。\n")
for href in re.findall(r'(?:href|src)="([^"]+)"',s):
 if not href.startswith("#"):assert (out/unquote(href.split("#")[0])).is_file(),href
files=sorted(p for p in out.rglob("*") if p.is_file())
assert len(files)==28
assert len(list((out/"Tables").glob("*.csv")))==5
assert not any(p.suffix.lower() in (".rds",".r",".py",".log",".gz",".json",".sh") for p in files)
archive=work/"SGS10_client_supplement_images_index_20260908.zip"
with zipfile.ZipFile(archive,"w",compression=zipfile.ZIP_DEFLATED,compresslevel=5) as z:
 for p in files:z.write(p,arcname=out.name+"/"+str(p.relative_to(out)))
with zipfile.ZipFile(archive) as z:
 assert z.testzip() is None
 assert len(z.namelist())==28
h=hashlib.sha256()
with archive.open("rb") as f:
 for b in iter(lambda:f.read(1024*1024),b""):h.update(b)
result={"archive":str(archive),"bytes":archive.stat().st_size,"sha256":h.hexdigest(),"files":len(files),
"image_files":20,"requested_tables":5,"index_and_guides":3,"source_dates_checked":len(rows),
"cutoff":"2026-04-13 inclusive","all_source_dates_in_range":True}
(work/"delivery_check.json").write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf-8")
print(json.dumps(result,ensure_ascii=False,indent=2))
# end
