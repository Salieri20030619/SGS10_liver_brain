import csv,json,zipfile,hashlib
from pathlib import Path
p=Path("/mnt/BioDisk/wangmingyuan/SGS10_analysis/TableS3_clean_20260907_v1")
def read(name):
    with (p/name).open(newline="",encoding="utf-8-sig") as f:return list(csv.DictReader(f))
validation=json.loads((p/"independent_validation.json").read_text())
assert validation["status"]=="PASS"
summary={"status":"PASS","marker_method":"descriptive effects and sample consistency; no inferential P/FDR",
         "independent_checks_passed":len(validation["checks"]),"organs":{}}
for organ,nct in [("brain",23),("liver",34)]:
    top=read(organ+"_top30_markers.csv");can=read(organ+"_all_positive_markers.csv")
    cells=read(organ+"_sample_celltype_counts.csv")
    summary["organs"][organ]={"cell_types":nct,"cells":sum(int(r["n_cells"]) for r in cells),
        "biological_samples":len({r["sample"] for r in cells}),"top30_rows":len(top),
        "all_positive_marker_rows":len(can)}
    assert len(top)==nct*30
    if organ=="brain":
        d=[r for r in top if r["cluster"]=="Microglia_DAM_strict"]
        summary["strict_DAM"]={"cells":57,"top_gene":d[0]["gene"],
            "top_gene_log2FC":float(d[0]["avg_log2FC"]),
            "top_gene_pct_target":float(d[0]["pct.1"]),
            "evaluable_samples":int(d[0]["n_evaluable_samples"]),
            "positive_samples":int(d[0]["n_positive_samples"])}
        summary["strict_DAM_top30_genes"]=[r["gene"] for r in d]
summary["limited_support"]={"liver_Basophil":{"cells":6,"evaluable_samples":0,"interpretation":"descriptive candidates only"}}
summary["warnings"]="285 named-scalar row-name warnings reproduced in table construction; all 285 reconstructed tables preserve numeric values."
(p/"RESULT_SUMMARY.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
lines=["# Table S3 clean annotation markers: completed",
       "","Liver: 34 subtypes, 136711 cells, 5 biological samples, 1020 Top-30 rows.",
       "Brain: 23 subtypes, 258990 cells, 5 biological samples, 690 Top-30 rows.",
       "Total: 57 subtypes and 1710 subtype-gene rows (not 1710 unique genes).",
       "","Independent numerical/type validation: 17/17 PASS.",
       "3130 retained DAM-cleaning decisions match the brain atlas annotations.",
       "","Primary table: Table_S3_clean_top30_markers_long.csv.",
       "See TableS3_clean_METHODS.md for the complete descriptive marker definition.",
       "Full gene statistics, positive candidates, sample support, metadata, code, logs and validation are included.",
       "","Low-support class: liver Basophil has only 6 cells in total; no sample meets the 10-target-cell consistency threshold.",
       "Strict DAM has 57 cells; its leading marker is Trem2, positive in all 4 evaluable samples.",
       "","No cell-level inferential P values or FDR were generated or copied from the older marker export.",
       "The 10-cell threshold is a measurement-support criterion, not the biological sample size.",
       "The named-scalar row-name warning has been reproduced; all 285 reconstructed tables retain identical numeric values.",
       "","The source RDS files remain in umap_paper_grade and are not duplicated in this package."]
(p/"DELIVERY_README.md").write_text("\n".join(lines)+"\n",encoding="utf-8")
files=[f for f in sorted(p.iterdir()) if f.is_file() and "preflight" not in f.name and f.suffix!=".pid"]
manifest=[]
for f in files:
    h=hashlib.sha256()
    with f.open("rb") as src:
        for chunk in iter(lambda:src.read(1024*1024),b""):h.update(chunk)
    manifest.append({"file":f.name,"bytes":f.stat().st_size,"sha256":h.hexdigest()})
(p/"OUTPUT_MANIFEST.json").write_text(json.dumps(manifest,indent=2),encoding="utf-8")
files.append(p/"OUTPUT_MANIFEST.json")
archive=p.parent/(p.name+".zip")
with zipfile.ZipFile(archive,"w",compression=zipfile.ZIP_DEFLATED,compresslevel=4,allowZip64=True) as z:
    for f in files:z.write(f,arcname=p.name+"/"+f.name,compress_type=zipfile.ZIP_STORED if f.suffix==".gz" else zipfile.ZIP_DEFLATED)
h=hashlib.sha256()
with archive.open("rb") as src:
    for chunk in iter(lambda:src.read(1024*1024),b""):h.update(chunk)
print(json.dumps({"summary":summary,"archive":str(archive),"archive_bytes":archive.stat().st_size,"archive_sha256":h.hexdigest()},indent=2))
