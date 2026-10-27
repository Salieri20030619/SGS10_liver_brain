import csv, json, math
from pathlib import Path
root=Path("/mnt/BioDisk/wangmingyuan/SGS10_analysis")
out=root/"TableS3_clean_20260907_v1"
checks=[]
def read(p):
    with open(p,encoding="utf-8-sig",newline="") as f:return list(csv.DictReader(f))
def check(n,b,d=""):
    checks.append(dict(check=n,pass_check=bool(b),detail=str(d)))
def num(x):return float(x) if x not in ("","NA",None) else float("nan")
def near(a,b):return math.isclose(a,b,rel_tol=1e-8,abs_tol=1e-8)
combined=[]
for organ,nct in [("brain",23),("liver",34)]:
    top=read(out/f"{organ}_top30_markers.csv")
    ref=read(root/"umap_paper_grade"/f"{organ}_subtype_centroids_paper.csv")
    types={r["subtype"] for r in ref}
    check(organ+"_type_set",types=={r["cluster"] for r in top} and len(types)==nct)
    check(organ+"_unique_gene_type",len({(r["cluster"],r["gene"]) for r in top})==len(top))
    counts=read(out/f"{organ}_sample_celltype_counts.csv")
    check(organ+"_five_samples",len({r["sample"] for r in counts})==5)
    check(organ+"_counts_match_ref",all(sum(int(r["n_cells"]) for r in counts if r["celltype"]==q["subtype"])==int(q["n"]) for q in ref))
    long=read(out/f"{organ}_top30_per_sample.csv")
    groups={}
    for r in long:groups.setdefault((r["cluster"],r["gene"]),[]).append(r)
    errors=[];maxdelta=0
    for r in top:
        rows=groups[(r["cluster"],r["gene"])]
        if len(rows)!=5:errors.append("sample_count");continue
        n1=sum(int(q["n_cells_target"]) for q in rows);n0=sum(int(q["n_cells_rest"]) for q in rows)
        a1=sum(num(q["mean_normalized_target"])*int(q["n_cells_target"]) for q in rows if int(q["n_cells_target"])>0)/n1
        a0=sum(num(q["mean_normalized_rest"])*int(q["n_cells_rest"]) for q in rows if int(q["n_cells_rest"])>0)/n0
        p1=sum(num(q["pct_target"])*int(q["n_cells_target"]) for q in rows if int(q["n_cells_target"])>0)/n1
        p0=sum(num(q["pct_rest"])*int(q["n_cells_rest"]) for q in rows if int(q["n_cells_rest"])>0)/n0
        lfc=math.log2(a1+1)-math.log2(a0+1)
        for a,b in [(lfc,num(r["avg_log2FC"])),(p1,num(r["pct.1"])),(p0,num(r["pct.2"]))]:
            maxdelta=max(maxdelta,abs(a-b))
            if not near(a,b):errors.append("pooled_reconstruction")
        valid=[q for q in rows if q["evaluable_for_consistency"]=="TRUE"]
        if len(valid)!=int(r["n_evaluable_samples"]):errors.append("valid_count")
        if valid:
            positive=sum(num(q["log2FC"])>0 for q in valid)
            specific=sum(num(q["log2FC"])>0 and num(q["delta_pct"])>0 for q in valid)
            if positive!=int(r["n_positive_samples"]) or specific!=int(r["n_positive_and_specific_samples"]):errors.append("sign_support")
            if not near(sum(num(q["log2FC"]) for q in valid)/len(valid),num(r["sample_mean_log2FC"])):errors.append("sample_mean")
        if num(r["avg_log2FC"])<.25 or num(r["pct.1"])<.25:errors.append("candidate_threshold")
    check(organ+"_independent_effect_reconstruction",not errors,{"errors":len(errors),"max_absolute_delta":maxdelta})
    check(organ+"_no_inferential_p_columns",not any(k.lower() in ("pval","p_val","p_val_adj","fdr") for k in top[0]))
    candidates=read(out/f"{organ}_all_positive_markers.csv")
    expected={(r["cluster"],r["gene"]) for r in candidates if int(r["rank_within_cluster"])<=30}
    check(organ+"_top30_matches_candidates",expected=={(r["cluster"],r["gene"]) for r in top})
    combined+=top
dec=read(out/"brain_DAM_clean_decision_concordance.csv")
present=[r for r in dec if r["present_in_integrated"]=="TRUE"]
check("brain_clean_decision_concordance",len(present)>0 and all(r["annotation_concordant"]=="TRUE" for r in present),len(present))
check("combined_57_clusters",len({(r["organ"],r["cluster"]) for r in combined})==57)
check("job_completed",(out/"COMPLETED.txt").exists())
res={"status":"PASS" if all(c["pass_check"] for c in checks) else "FAIL","checks":checks}
(out/"independent_validation.json").write_text(json.dumps(res,ensure_ascii=False,indent=2),encoding="utf-8")
print(json.dumps(res,ensure_ascii=False,indent=2))
if res["status"]!="PASS":raise SystemExit(1)
