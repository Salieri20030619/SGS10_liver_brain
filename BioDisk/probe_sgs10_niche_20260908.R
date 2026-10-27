suppressPackageStartupMessages(library(Seurat))
p <- '/mnt/BioDisk/wangmingyuan/SGS10_analysis'
s <- 'Liver_3'
x <- readRDS(file.path(p,'spatial_results',paste0(s,'_seekspace.rds')))
cat('MISC\n'); str(x@misc,max.level=2)
cat('IMAGES\n'); print(names(x@images)); cat('REDUCTIONS\n'); print(names(x@reductions))
cat('META\n'); print(head(as.data.frame(x@meta.data),2))
