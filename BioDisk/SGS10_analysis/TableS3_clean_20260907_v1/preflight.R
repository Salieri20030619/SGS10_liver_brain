suppressPackageStartupMessages(library(Seurat))
root <- '/mnt/BioDisk/wangmingyuan/SGS10_analysis'
out <- file.path(root,'TableS3_clean_20260907_v1')
for (organ in c('brain','liver')) {
 cat('\nORGAN',organ,format(Sys.time()),'\n'); flush.console()
 obj <- readRDS(file.path(root,'umap_paper_grade',paste0(organ,'_integrated.rds')))
 cat('DIM',dim(obj),'ASSAYS',Assays(obj),'DEFAULT',DefaultAssay(obj),'\n')
 cat('LAYERS',Layers(obj[['RNA']]),'\n')
 print(head(as.data.frame(obj@meta.data),2))
 print(table(obj$celltype_refined,useNA='ifany'))
 write.csv(as.data.frame(obj@meta.data),file.path(out,paste0(organ,'_metadata_preflight.csv')),row.names=TRUE)
 cat('GENES',head(rownames(obj[['RNA']])), '\n')
 rm(obj);gc()
}
cat('PREFLIGHT_COMPLETE\n')
