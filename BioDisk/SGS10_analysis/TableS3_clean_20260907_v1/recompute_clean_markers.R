#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(Seurat);library(Matrix)})
options(stringsAsFactors=FALSE)
root <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis"
out <- file.path(root,"TableS3_clean_20260907_v1")
dir.create(out,showWarnings=FALSE)
stamp <- function(...) {cat(format(Sys.time(),"%FT%T%z"),...,"\n");flush.console()}
write_csv <- function(d,p) write.csv(as.data.frame(d),file.path(out,p),row.names=FALSE,na="")
standard_sample <- function(v,organ) {
 z <- gsub("[^a-z0-9]","",tolower(as.character(v)))
 keys <- if(organ=="brain") c("brain3","brain8","nao56","blabrain","conbrain") else c("liver3","liver8","gan56","blaliver","conliver")
 vals <- if(organ=="brain") c("brain3","brain8","nao56","bla_brain","con_brain") else c("liver3","liver8","gan56","bla_liver","con_liver")
 ans <- rep(NA_character_,length(z))
 for(i in seq_along(keys)) ans[grepl(keys[i],z,fixed=TRUE)] <- vals[i]
 ans
}
get_samples <- function(md,organ) {
 for(nm in c("Sample","sample","sample_id","orig.ident","raw_Sample")) if(nm %in% names(md)) {
  s <- standard_sample(md[[nm]],organ)
  if(all(!is.na(s)) && length(unique(s))==5L) return(s)
 }
 s <- standard_sample(rownames(md),organ)
 if(all(!is.na(s)) && length(unique(s))==5L) return(s)
 stop("Cannot identify five biological samples for ",organ)
}
bare_barcode <- function(z) {
 m <- regexpr("[ACGT]{15,}",z)
 ans <- rep(NA_character_,length(z));ok<-m>0
 ans[ok] <- regmatches(z,m)
 ans
}
all_top <- list(); all_counts <- list(); checks <- list()
record_check <- function(name,ok,detail="") {
 checks[[length(checks)+1L]] <<- data.frame(check=name,pass=ok,detail=detail)
 write_csv(do.call(rbind,checks),"validation_checks.csv")
 if(!isTRUE(ok)) stop("Validation failed: ",name," ",detail)
}
writeLines(capture.output(sessionInfo()),file.path(out,"session_info.txt"))
inputs <- c(file.path(root,"umap_paper_grade",paste0(c("brain","liver"),"_integrated.rds")),
            file.path(root,"umap_paper_grade",paste0(c("brain","liver"),"_subtype_centroids_paper.csv")),
            file.path(root,"dam_clean_v2","all_dam_decisions.csv"))
fi <- file.info(inputs)
write_csv(data.frame(path=inputs,size_bytes=fi$size,mtime=as.character(fi$mtime)),
          "input_manifest.csv")
stamp("START descriptive annotation marker recomputation; no inferential tests")
for(organ in c("brain","liver")) {
 stamp("LOAD",organ)
 refpath <- file.path(root,"umap_paper_grade",paste0(organ,"_subtype_centroids_paper.csv"))
 ref <- read.csv(refpath,check.names=FALSE)
 types <- ref$subtype[order(ref$idx)]
 record_check(paste0(organ,"_reference_types"),length(types)==if(organ=="brain")23L else 34L)
 obj <- readRDS(file.path(root,"umap_paper_grade",paste0(organ,"_integrated.rds")))
 md <- as.data.frame(obj@meta.data);md$sample_canonical<-get_samples(md,organ)
 labels <- as.character(md$celltype_refined)
 selected <- !is.na(labels) & labels %in% types
 ct <- table(factor(labels[selected],levels=ref$subtype))
 record_check(paste0(organ,"_centroid_cell_counts"),all(as.integer(ct)==ref$n),
              paste("selected",sum(selected),"input",nrow(md)))
 excluded <- as.data.frame(table(label=labels[!selected],useNA="ifany"))
 write_csv(excluded,paste0(organ,"_excluded_annotation_counts.csv"))
 metaout <- data.frame(barcode=rownames(md),sample=md$sample_canonical,organ=organ,
                       celltype_refined=labels,included=selected,
                       inclusion_reason=ifelse(selected,"paper_atlas_annotation","outside_paper_atlas_annotation"))
 for(nm in intersect(c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","percent.mito"),names(md))) metaout[[nm]]<-md[[nm]]
 write_csv(metaout,paste0(organ,"_clean_atlas_cell_metadata.csv"))
 if(organ=="brain") {
  dec <- read.csv(file.path(root,"dam_clean_v2","all_dam_decisions.csv"),check.names=FALSE)
  mapping <- c(Microglia_DAM_strict="Microglia_DAM_strict",Microglia_other="Microglia_Homeostatic",
    Astrocyte="Astrocyte_Resting",Pericyte="Pericytes",Endothelial="Endothelial",
    Oligodendrocyte="Oligo_Mature",Neuron="Neurons_unassigned",
    low_quality_high_mt="low_quality",low_quality_no_lineage="low_quality")
  key <- paste(md$sample_canonical,bare_barcode(rownames(md)),sep="|")
  record_check("brain_unique_sample_barcode",!anyDuplicated(key))
  ix <- match(paste(dec$sample,dec$barcode,sep="|"),key)
  dec$expected_annotation <- unname(mapping[dec$new_label])
  dec$present_in_integrated <- !is.na(ix)
  dec$integrated_annotation <- labels[ix]
  dec$included_in_marker_atlas <- selected[ix]
  dec$annotation_concordant <- dec$expected_annotation==dec$integrated_annotation
  write_csv(dec,"brain_DAM_clean_decision_concordance.csv")
  record_check("brain_clean_decisions_match",sum(!is.na(ix))>0 && all(dec$annotation_concordant[!is.na(ix)]),
               paste(sum(!is.na(ix)),"of",nrow(dec),"decisions present in integrated object"))
  record_check("brain_strict_DAM_57",sum(labels[selected]=="Microglia_DAM_strict")==57L)
 }
 stamp("RNA COUNTS",organ)
 cl <- Layers(obj[["RNA"]],search="^counts")
 record_check(paste0(organ,"_RNA_counts_present"),length(cl)>0)
 if(length(cl)>1L) obj[["RNA"]] <- JoinLayers(obj[["RNA"]],layers="counts",new="counts")
 x <- GetAssayData(obj,assay="RNA",layer=if(length(cl)>1L)"counts" else cl[1])
 record_check(paste0(organ,"_matrix_metadata_alignment"),identical(colnames(x),rownames(md)))
 x <- x[,selected,drop=FALSE]
 samples <- sort(unique(md$sample_canonical[selected]))
 lab <- labels[selected];smp<-md$sample_canonical[selected]
 design <- expand.grid(celltype=types,sample=samples,stringsAsFactors=FALSE)
 group <- match(paste(smp,lab,sep="|"),paste(design$sample,design$celltype,sep="|"))
 ng <- nrow(design);gn<-tabulate(group,nbins=ng);design$n_cells<-gn;design$organ<-organ
 design$evaluable_for_consistency <- gn>=10L
 write_csv(design,paste0(organ,"_sample_celltype_counts.csv"));all_counts[[organ]]<-design
 G <- sparseMatrix(i=seq_along(group),j=group,x=1,dims=c(length(group),ng))
 libsize <- Matrix::colSums(x)
 record_check(paste0(organ,"_positive_library_sizes"),all(libsize>0))
 record_check(paste0(organ,"_RNA_count_values"),all(x@x>=0) && all(is.finite(x@x)) &&
              all(abs(x@x-round(x@x))<1e-7))
 gene <- rownames(x)
 record_check(paste0(organ,"_unique_gene_ids"),!anyDuplicated(gene))
 stamp("AGGREGATE detection",organ,length(gene),"genes",ncol(x),"cells")
 det<-x;det@x[]<-1
 D<-as.matrix(det %*% G);rm(det);gc()
 stamp("AGGREGATE library-normalized expression",organ)
 x@x <- x@x*rep.int(10000/libsize,diff(x@p))
 A<-as.matrix(x %*% G)
 rm(x,obj,md,metaout,G);gc()
 totalsA <- sapply(samples,function(s)rowSums(A[,design$sample==s,drop=FALSE]))
 totalsD <- sapply(samples,function(s)rowSums(D[,design$sample==s,drop=FALSE]))
 totalsN <- sapply(samples,function(s)sum(gn[design$sample==s]))
 totalA<-rowSums(totalsA);totalD<-rowSums(totalsD);totalN<-sum(gn)
 fullcon<-gzfile(file.path(out,paste0(organ,"_all_genes_marker_statistics.csv.gz")),"wt")
 first<-TRUE;toplist<-list();candlist<-list();sampletop<-list()
 for(k in seq_along(types)) {
  typ<-types[k];idx<-match(paste(samples,typ,sep="|"),paste(design$sample,design$celltype,sep="|"))
  n1<-gn[idx];n0<-totalsN-n1;valid<-n1>=10 & n0>=10
  a1<-A[,idx,drop=FALSE];d1<-D[,idx,drop=FALSE]
  mu1<-sweep(a1,2,pmax(n1,1),"/");mu0<-sweep(totalsA-a1,2,pmax(n0,1),"/")
  p1<-sweep(d1,2,pmax(n1,1),"/");p0<-sweep(totalsD-d1,2,pmax(n0,1),"/")
  lfc<-log2(mu1+1)-log2(mu0+1);delta<-p1-p0
  lfc[,n1==0 | n0==0]<-NA_real_;delta[,n1==0 | n0==0]<-NA_real_
  ntype<-sum(n1);nrest<-totalN-ntype
  pooled1<-rowSums(a1)/ntype;pooled0<-(totalA-rowSums(a1))/nrest
  pct1<-rowSums(d1)/ntype;pct0<-(totalD-rowSums(d1))/nrest
  if(any(valid)) {
   meanfc<-rowMeans(lfc[,valid,drop=FALSE]);minfc<-apply(lfc[,valid,drop=FALSE],1,min)
   maxfc<-apply(lfc[,valid,drop=FALSE],1,max)
   npos<-rowSums(lfc[,valid,drop=FALSE]>0)
   nspecific<-rowSums(lfc[,valid,drop=FALSE]>0 & delta[,valid,drop=FALSE]>0)
  } else {meanfc<-minfc<-maxfc<-rep(NA_real_,length(gene));npos<-nspecific<-rep(NA_integer_,length(gene))}
  df<-data.frame(organ=organ,cluster=typ,gene=gene,avg_log2FC=log2(pooled1+1)-log2(pooled0+1),
    pct.1=pct1,pct.2=pct0,delta_pct=pct1-pct0,mean_normalized_target=pooled1,
    mean_normalized_rest=pooled0,n_cells_target=ntype,n_cells_rest=nrest,
    n_samples_with_target=sum(n1>0),n_evaluable_samples=sum(valid),
    sample_mean_log2FC=meanfc,sample_min_log2FC=minfc,sample_max_log2FC=maxfc,
    n_positive_samples=npos,n_positive_and_specific_samples=nspecific,
    positive_sample_fraction=if(any(valid))npos/sum(valid) else NA_real_)
  df$marker_candidate<-df$avg_log2FC>=0.25 & df$pct.1>=0.25
  write.table(df,fullcon,sep=",",row.names=FALSE,col.names=first,quote=TRUE,na="")
  first<-FALSE
  cand<-df[df$marker_candidate,]
  cand<-cand[order(-cand$avg_log2FC,-cand$delta_pct,cand$gene),]
  cand$rank_within_cluster<-seq_len(nrow(cand))
  candlist[[typ]]<-cand
  top<-head(cand,30L);toplist[[typ]]<-top
  gi<-match(top$gene,gene)
  if(length(gi)) for(j in seq_along(samples)) sampletop[[paste(typ,samples[j])]]<-data.frame(
    organ=organ,cluster=typ,gene=gene[gi],sample=samples[j],
    n_cells_target=n1[j],n_cells_rest=n0[j],evaluable_for_consistency=valid[j],
    mean_normalized_target=if(n1[j]>0)mu1[gi,j] else NA_real_,
    mean_normalized_rest=if(n0[j]>0)mu0[gi,j] else NA_real_,
    pct_target=if(n1[j]>0)p1[gi,j] else NA_real_,
    pct_rest=if(n0[j]>0)p0[gi,j] else NA_real_,
    log2FC=lfc[gi,j],delta_pct=delta[gi,j])
  stamp("MARKERS",organ,typ,"candidates",nrow(cand),"top",nrow(top),"evaluable_samples",sum(valid))
 }
 close(fullcon)
 top<-do.call(rbind,toplist);candidates<-do.call(rbind,candlist)
 write_csv(top,paste0(organ,"_top30_markers.csv"))
 write_csv(candidates,paste0(organ,"_all_positive_markers.csv"))
 write_csv(do.call(rbind,sampletop),paste0(organ,"_top30_per_sample.csv"))
 wide<-data.frame(rank=1:30)
 for(typ in types) wide[[typ]]<-c(toplist[[typ]]$gene,rep(NA_character_,30-nrow(toplist[[typ]])))
 write_csv(wide,paste0(organ,"_top30_markers_wide.csv"))
 coverage<-data.frame(organ=organ,cluster=types,n_marker_candidates=sapply(candlist,nrow),
                     n_top30=sapply(toplist,nrow))
 write_csv(coverage,paste0(organ,"_marker_coverage.csv"))
 record_check(paste0(organ,"_all_types_have_markers"),all(coverage$n_top30>0))
 record_check(paste0(organ,"_top_marker_values"),all(is.finite(top$avg_log2FC)) &&
              all(top$pct.1>=0 & top$pct.1<=1) && all(top$pct.2>=0 & top$pct.2<=1))
 all_top[[organ]]<-top
 rm(A,D,totalsA,totalsD,candlist,candidates,toplist,sampletop);gc()
}
write_csv(do.call(rbind,all_top),"Table_S3_clean_top30_markers_long.csv")
write_csv(do.call(rbind,all_counts),"Table_S3_sample_celltype_counts.csv")
record_check("combined_57_types",length(unique(paste(do.call(rbind,all_top)$organ,do.call(rbind,all_top)$cluster)))==57L)
writeLines(c("COMPLETE",format(Sys.time(),"%FT%T%z")),file.path(out,"COMPLETED.txt"))
stamp("COMPLETE")
