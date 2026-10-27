import os,json,csv,subprocess,re
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
root=Path("/mnt/BioDisk/wangmingyuan/SGS10_analysis")
work=root/"client_supplement_package_20260908_work"
data=json.loads((work/"build_index.json").read_text())
srcs=set()
for row in data["index"]:
 for s in row["source"].split(";"):
  s=s.strip()
  if s and "*" not in s:srcs.add(s)
for organ in ("brain","liver"):
 for name in ("per_sample_domain","microenv_naming","spot_microenv_assignment"):
  srcs.add(f"taskB/spatial_domain_RCTD100/{organ}_{name}.csv")
 srcs.add(f"D1_v2_RCTD_composition/microenv_key_celltypes_RCTD100/{organ}_microenv_celltypes_top5.csv")
 srcs.add(f"D1_v2_RCTD_composition/microenv_markers_RCTD100/{organ}_microenv_markers_top20.csv")
 srcs.add(f"umap_paper_grade/{organ}_integrated.rds")
 srcs.add(f"umap_paper_grade/{organ}_subtype_centroids_paper.csv")
srcs.add("dam_clean_v2/all_dam_decisions.csv")
zone=ZoneInfo("Asia/Shanghai");cut=datetime(2026,4,14,tzinfo=zone)
rows=[]
for s in sorted(srcs):
 p=root/s
 if not p.exists():continue
 mt=datetime.fromtimestamp(p.stat().st_mtime,zone)
 pdf="";pdfdt=None
 if p.suffix.lower()==".pdf":
  r=subprocess.run(["pdfinfo",str(p)],capture_output=True,text=True,env=dict(os.environ,LC_ALL="C",TZ="Asia/Shanghai"))
  m=re.search(r"^CreationDate:\s+(.*)$",r.stdout,re.M)
  if m:
   pdf=m.group(1)
   try:pdfdt=datetime.strptime(" ".join(pdf.split()[:5]),"%a %b %d %H:%M:%S %Y").replace(tzinfo=zone)
   except ValueError:pass
 eligible=mt>=cut and (pdfdt is None or pdfdt>=cut)
 rows.append({"source":s,"mtime":mt.isoformat(),"pdf_creation":pdf,"eligible_after_Apr13":eligible,
              "evidence":"PDF creation and file mtime" if pdfdt else "file mtime; no embedded creation date",
              "bytes":p.stat().st_size})
(work/"source_time_audit.json").write_text(json.dumps(rows,ensure_ascii=False,indent=2),encoding="utf-8")
with (work/"source_time_audit.csv").open("w",encoding="utf-8-sig",newline="") as f:
 w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
candidates=[]
for base in ["taskB_spatial_figs_v1","taskB_spatial_figs_v2","taskB_spatial_figs_v3","taskB_spatial_figs_2026-04-20_v1","taskB_spatial_figs_2026-04-20_v2","taskB_spatial_figs_2026-04-20_v3","umap_paper_grade","taskA_figures_v1","taskA_figures_2026-04-20_v1","taskA_figures_v1","brain_export","brain_results_export","liver_download_final"]:
 for p in (root/base).rglob("*"):
  if p.is_file() and p.suffix.lower() in (".png",".pdf",".csv") and any(k in p.name.lower() for k in ["umap","qc","summary"]):
   t=datetime.fromtimestamp(p.stat().st_mtime,zone)
   if t>=cut:candidates.append({"path":str(p.relative_to(root)),"mtime":t.isoformat()})
print(json.dumps({"total_sources":len(rows),"out_of_range":[r for r in rows if not r["eligible_after_Apr13"]],
 "eligible_count":sum(r["eligible_after_Apr13"] for r in rows),"later_candidates":candidates[:100]},ensure_ascii=False,indent=2))
