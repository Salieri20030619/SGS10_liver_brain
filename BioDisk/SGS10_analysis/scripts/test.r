
library(Seurat)
cat('Seurat版本:', as.character(packageVersion('Seurat')), '\n')

# 创建测试数据
test_data <- matrix(sample(0:10, 1000, replace=TRUE), nrow=100)
rownames(test_data) <- paste0('Gene', 1:100)
colnames(test_data) <- paste0('Cell', 1:10)

# 最简形式
cat('\n尝试最简形式...\n')
obj <- CreateSeuratObject(counts = test_data)
cat('✓ 成功！细胞数:', ncol(obj), '\n')

# 带项目名称
cat('\n尝试带项目名称...\n')
obj <- CreateSeuratObject(counts = test_data, project = 'test')
cat('✓ 成功！项目名:', obj@project.name, '\n')

# 查看可用参数
cat('\nCreateSeuratObject参数:\n')
args(CreateSeuratObject)
