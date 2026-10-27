
args <- commandArgs(trailingOnly = TRUE)
sample_name <- args[1]

library(Seurat)

data_paths <- list(
  brain3 = '/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/brain_samples/Brain_3_filtered_feature_bc_matrix',
  brain8 = '/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/brain_samples/Brain_8_filtered_feature_bc_matrix',
  bla_brain = '/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/brain_samples/Bla_brain_1_filtered_feature_bc_matrix',
  con_brain = '/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/brain_samples/Con_brain_1_filtered_feature_bc_matrix'
)

data_dir <- data_paths[[sample_name]]
cat('加载样本:', sample_name, '
')
cat('数据目录:', data_dir, '
')

if(dir.exists(data_dir)) {
  cat('✓ 目录存在
')
  files <- list.files(data_dir)
  cat('文件:', paste(files, collapse=', '), '
')
  
  data <- Read10X(data_dir)
  cat('数据维度:', nrow(data), 'x', ncol(data), '
')
  
  obj <- CreateSeuratObject(counts = data)
  cat('✓ Seurat对象创建成功，细胞数:', ncol(obj), '
')
} else {
  cat('✗ 目录不存在
')
}
