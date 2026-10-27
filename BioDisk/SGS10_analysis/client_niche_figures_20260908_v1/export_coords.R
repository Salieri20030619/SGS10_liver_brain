suppressPackageStartupMessages(library(Matrix))
p <- '/mnt/BioDisk/wangmingyuan/SGS10_analysis'
out <- file.path(p,'client_niche_figures_20260908_v1')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
samples <- c('Liver_3','Liver_8','gan56','Brain_3','Brain_8','nao56')
checks <- list()
for (s in samples) {
  cat(format(Sys.time()),s,'loading RCTD\n'); flush.console()
  f <- file.path(p,'spatial_results/rctd_full/per_sample',paste0(s,'_rctd_results.rds'))
  r <- readRDS(f)
  co <- r@spatialRNA@coords
  w <- as.matrix(r@results$weights)
  wc <- as.matrix(read.csv(file.path(dirname(f),paste0(s,'_weights_matrix.csv')),row.names=1,check.names=FALSE))
  stopifnot(!anyDuplicated(rownames(co)),!anyDuplicated(rownames(w)))
  stopifnot(setequal(rownames(w),rownames(wc)),setequal(colnames(w),colnames(wc)))
  wc <- wc[rownames(w),colnames(w),drop=FALSE]
  w <- w / pmax(rowSums(w),1e-10); wc <- wc / pmax(rowSums(wc),1e-10)
  err <- max(abs(w-wc))
  stopifnot(err<1e-10)
  write.csv(data.frame(barcode=rownames(co),co,check.names=FALSE),file.path(out,paste0(s,'_rctd_coords.csv')),row.names=FALSE)
  checks[[s]] <- data.frame(sample=s,coordinate_rows=nrow(co),weight_rows=nrow(w),max_weight_csv_rds_difference=err)
  cat(s,'EXPORTED',nrow(co),'coordinates; CSV/RDS error',err,'\n'); flush.console()
  rm(r,co,w,wc);gc()
}
write.csv(do.call(rbind,checks),file.path(out,'rctd_source_validation.csv'),row.names=FALSE)
cat('EXPORT COMPLETE\n')
