p <- '/mnt/BioDisk/wangmingyuan/SGS10_analysis'
o <- file.path(p,'client_niche_figures_20260908_v1')
expected <- read.csv(file.path(o,'niche_composition_all_celltypes.csv'),check.names=FALSE)
roi <- read.csv(file.path(o,'zoom_ROI_manifest.csv'),check.names=FALSE)
audit <- list()
for (organ in c('liver','brain')) {
 a <- read.csv(file.path(p,'taskB/spatial_domain_RCTD100',paste0(organ,'_spot_microenv_assignment.csv')))
 target <- if(organ=='liver') 'M05' else 'M04'
 ss <- if(organ=='liver') c('Liver_3','Liver_8','gan56') else c('Brain_3','Brain_8','nao56')
 for(s in ss) {
  w <- as.matrix(read.csv(file.path(p,'spatial_results/rctd_full/per_sample',paste0(s,'_weights_matrix.csv')),row.names=1,check.names=FALSE))
  w[w<0] <- 0;w <- w/rowSums(w)
  selected <- a$spot_id[a$sample==s & a$microenv_id==target]
  selected <- substring(selected,nchar(s)+3)
  stopifnot(all(selected %in% rownames(w)))
  calc <- colMeans(w[selected,,drop=FALSE])*100
  e <- expected[expected$sample==s,]
  err <- max(abs(calc[e$celltype]-e$mean_RCTD_weight_pct))
  stopifnot(err<1e-10,abs(sum(calc)-100)<1e-8)
  co <- read.csv(file.path(o,paste0(s,'_rctd_coords.csv')))
  rr <- roi[roi$sample==s,]
  if(organ=='brain') {
    co$x <- (if(s=='nao56') 55050 else 56500)-co$x
    co$y <- (if(s=='nao56') 19906 else 20434)-co$y
  }
  inside <- with(co,x>=rr$x_min & x<=rr$x_max & y>=rr$y_min & y<=rr$y_max)
  stopifnot(sum(inside)==rr$all_spots_in_ROI,sum(inside & co$barcode %in% selected)==rr$target_spots_in_ROI)
  plotted <- read.csv(file.path(o,paste0(organ,'_',target,'_composition_plotted.csv')),check.names=FALSE)
  mapping <- read.csv(file.path(o,paste0(organ,'_reference_to_cell_class.csv')))
  stopifnot(setequal(mapping$celltype,names(calc)),!anyDuplicated(mapping$celltype))
  broad <- tapply(calc[mapping$celltype],mapping$cell_class,sum)
  keep <- plotted$celltype[plotted$celltype!='Other cell classes']
  values <- c(broad[keep],sum(broad[setdiff(names(broad),keep)]))
  perr <- max(abs(values-plotted[[s]]))
  stopifnot(perr<1e-10)
  audit[[s]] <- data.frame(sample=s,organ=organ,n_target_spots=length(selected),composition_max_abs_error_pct=err,plotted_max_abs_error_pct=perr,roi_target_spots=sum(inside & co$barcode %in% selected),roi_all_spots=sum(inside),status='PASS')
 }
}
res <- do.call(rbind,audit)
write.csv(res,file.path(o,'independent_R_validation.csv'),row.names=FALSE)
print(as.data.frame(res));cat('INDEPENDENT R CHECKS PASS FOR ALL SIX SAMPLES\n')
