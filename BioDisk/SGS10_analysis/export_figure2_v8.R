suppressPackageStartupMessages({library(Matrix);library(RANN)})
p<-'/mnt/BioDisk/wangmingyuan/SGS10_analysis';out<-file.path(p,'Figure2_full_v8');dir.create(out,showWarnings=FALSE)
ex<-parse(file.path(p,'scripts/he_overlay_RCTD100.R'))
for(e in ex)if(is.call(e)&&identical(e[[1]],as.name('<-'))&&as.character(e[[2]])%in%c('smooth_weights','assign_domain'))eval(e)
counts<-list();means<-list()
samples<-c('Bla_liver_1','Con_liver_1','Liver_3','Liver_8','gan56','Bla_brain_1','Con_brain_1','Brain_3','Brain_8','nao56')
for(s in samples){
 organ<-if(s%in%samples[1:5])'liver' else 'brain'
 cat(s,'START\n');flush.console()
 w<-as.matrix(read.csv(file.path(p,'spatial_results/rctd_full/per_sample',paste0(s,'_weights_matrix.csv')),row.names=1,check.names=FALSE));w<-w/pmax(rowSums(w),1e-9)
 old<-file.path(p,'client_niche_RCTD100_HE_20260908_v6',paste0(s,'_full_HE_assignments.csv.gz'))
 if(file.exists(old)){
  a<-read.csv(old);rownames(a)<-a$barcode;stopifnot(setequal(rownames(w),rownames(a)));a<-a[rownames(w),]
 }else{
  obj<-readRDS(file.path(p,'spatial_results/rctd_full/per_sample',paste0(s,'_rctd_results.rds')));co<-obj@spatialRNA@coords
  ww<-as.matrix(obj@results$weights);ww<-ww/pmax(rowSums(ww),1e-9)
  stopifnot(setequal(rownames(w),rownames(co)),max(abs(w-ww[rownames(w),colnames(w)]))<1e-10)
  co<-co[rownames(w),];knnco<-co
  if(organ=='brain'){knnco$x<-56500-co$x;knnco$y<-20434-co$y}
  prof<-as.matrix(read.csv(file.path(p,'taskB/spatial_domain_RCTD100',paste0(organ,'_domain_ct_profile.csv')),row.names=1,check.names=FALSE));me<-read.csv(file.path(p,'taskB/spatial_domain_RCTD100',paste0(organ,'_microenv_naming.csv')))
  domain<-assign_domain(smooth_weights(knnco,w,k=15),prof)
  a<-data.frame(barcode=rownames(w),sample=s,x=co$x,y=co$y,domain=domain,microenv_id=me$microenv_id[match(domain,me$domain)])
  gz<-gzfile(file.path(out,paste0(s,'_full_HE_assignments.csv.gz')),'wt');write.csv(a,gz,row.names=FALSE);close(gz);rm(obj,ww);gc()
 }
 stopifnot(!anyNA(a$microenv_id));wp<-pmax(w,0);wp<-wp/rowSums(wp)
 for(m in unique(a$microenv_id)){
  hit<-a$microenv_id==m;key<-paste(s,m)
  counts[[key]]<-data.frame(organ=organ,sample=s,microenv_id=m,n_spots=sum(hit),sample_total=nrow(w),fraction_pct=100*mean(hit))
  av<-colMeans(wp[hit,,drop=FALSE]);means[[key]]<-data.frame(organ=organ,sample=s,microenv_id=m,celltype=colnames(w),mean_weight=av,n_spots=sum(hit))
 }
 cat(s,'DONE',nrow(w),'\n');flush.console()
}
write.csv(do.call(rbind,counts),file.path(out,'microenvironment_proportions.csv'),row.names=FALSE)
write.csv(do.call(rbind,means),file.path(out,'microenvironment_composition_per_sample.csv'),row.names=FALSE)
cat('FIGURE2 DATA READY\n')
