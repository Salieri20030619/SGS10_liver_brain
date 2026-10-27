suppressPackageStartupMessages({library(Matrix);library(RANN)})
p <- '/mnt/BioDisk/wangmingyuan/SGS10_analysis'
out <- file.path(p,'client_niche_RCTD100_HE_20260908_v6')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
# Evaluate only the two original helper definitions, preserving their exact behavior.
original <- parse(file.path(p,'scripts/he_overlay_RCTD100.R'))
for(e in original) if(is.call(e) && identical(e[[1]],as.name('<-')) && as.character(e[[2]]) %in% c('smooth_weights','assign_domain')) eval(e)
stopifnot(exists('smooth_weights'),exists('assign_domain'))
checks <- list(); comps <- list(); sources <- character()
for(s in c('Liver_3','Liver_8','gan56','Brain_3','Brain_8','nao56')) {
 cat(format(Sys.time()),s,'START\n');flush.console()
 organ <- if(s %in% c('Liver_3','Liver_8','gan56')) 'liver' else 'brain'
 target <- if(organ=='liver') 'M05' else 'M04'
 f <- file.path(p,'spatial_results/rctd_full/per_sample',paste0(s,'_rctd_results.rds'))
 wf <- file.path(dirname(f),paste0(s,'_weights_matrix.csv'))
 pf <- file.path(p,'taskB/spatial_domain_RCTD100',paste0(organ,'_domain_ct_profile.csv'))
 mf <- file.path(dirname(pf),paste0(organ,'_microenv_naming.csv'))
 af <- file.path(dirname(pf),paste0(organ,'_spot_microenv_assignment.csv'))
 sources <- unique(c(sources,f,wf,pf,mf,af))
 obj <- readRDS(f);co <- obj@spatialRNA@coords;w <- as.matrix(obj@results$weights)
 wc <- as.matrix(read.csv(wf,row.names=1,check.names=FALSE))
 a <- read.csv(af);a <- a[a$sample==s,];rownames(a)<-sub(paste0(s,'__'),'',a$spot_id,fixed=TRUE)
 stopifnot(!anyDuplicated(rownames(co)),!anyDuplicated(rownames(w)),setequal(rownames(co),rownames(w)),setequal(rownames(w),rownames(wc)),setequal(rownames(w),rownames(a)))
 co<-co[rownames(w),,drop=FALSE];wc<-wc[rownames(w),colnames(w),drop=FALSE]
 w<-w/pmax(rowSums(w),1e-9);wc<-wc/pmax(rowSums(wc),1e-9)
 err<-max(abs(w-wc));stopifnot(is.finite(err),err<1e-10)
 prof<-as.matrix(read.csv(pf,row.names=1,check.names=FALSE));me<-read.csv(mf)
 # Historical brain orientation is applied before neighbor search, including tie behavior.
 sx<-if(s %in% c('gan56','nao56')) 55050 else 56500
 sy<-if(s %in% c('gan56','nao56')) 19906 else 20434
 knnco<-co
 if(organ=='brain'){knnco$x<-sx-co$x;knnco$y<-sy-co$y}
 sm<-smooth_weights(knnco,w,k=15)
 domain<-assign_domain(sm,prof);micro<-me$microenv_id[match(domain,me$domain)]
 stopifnot(!anyNA(micro),all(is.finite(sm)),max(abs(rowSums(sm)-1))<1e-8)
 dat<-data.frame(barcode=rownames(w),sample=s,x=co$x,y=co$y,domain=domain,microenv_id=micro,previous_microenv=a[rownames(w),'microenv_id'])
 con<-gzfile(file.path(out,paste0(s,'_full_HE_assignments.csv.gz')),'wt');write.csv(dat,con,row.names=FALSE);close(con)
 hit<-micro==target;stopifnot(any(hit))
 # Match established nonnegative composition reporting; labels use original unmodified algorithm.
 wp<-pmax(w,0);wp<-wp/rowSums(wp)
 av<-colMeans(wp[hit,,drop=FALSE])*100
 comps[[s]]<-data.frame(organ=organ,sample=s,microenv_id=target,celltype=names(av),mean_RCTD_weight_pct=av,n_target_spots=sum(hit))
 sampled_f<-file.path(p,'spatial_results/rctd/per_sample',paste0(s,'_weights_matrix.csv'))
 n20<-if(file.exists(sampled_f)) nrow(read.csv(sampled_f,row.names=1)) else NA_integer_
 checks[[s]]<-data.frame(sample=s,full_weight_rows=nrow(w),full_coordinate_rows=nrow(co),previous_assignment_rows=nrow(a),sampled_version_rows=n20,plotted_rows=nrow(dat),csv_rds_max_error=err,previous_target_spots=sum(a$microenv_id==target),new_target_spots=sum(hit),changed_labels=sum(dat$microenv_id!=dat$previous_microenv),composition_clipping_max_pp=max(abs(av-colMeans(w[hit,,drop=FALSE])*100)))
 cat(s,'DONE',nrow(dat),'spots;',sum(hit),'target;',checks[[s]]$changed_labels,'labels changed\n');flush.console()
 rm(obj,co,w,wc,sm,wp);gc()
}
write.csv(do.call(rbind,checks),file.path(out,'full_RCTD_audit.csv'),row.names=FALSE)
write.csv(do.call(rbind,comps),file.path(out,'composition_all_celltypes.csv'),row.names=FALSE)
writeLines(sources,file.path(out,'server_source_paths.txt'))
writeLines(capture.output(sessionInfo()),file.path(out,'R_sessionInfo.txt'))
cat('ALL SIX COMPLETE\n')
