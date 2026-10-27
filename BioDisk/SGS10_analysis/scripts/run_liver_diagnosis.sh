# ==================== 肝脏样本数据诊断和修复脚本 ====================
# 文件名: ~/SGS10_analysis/scripts/liver_samples_diagnosis_fix.R
# 用法: Rscript liver_samples_diagnosis_fix.R sample_name

args <- commandArgs(trailingOnly = TRUE)
sample_name <- args[1]

cat("========================================\n")
cat("诊断肝脏样本:", sample_name, "\n")
cat("========================================\n\n")

library(Seurat)
library(ggplot2)
library(dplyr)

sample_dir <- file.path("/mnt/BioDisk/wangmingyuan/SGS10_analysis", sample_name)
setwd(sample_dir)

# 加载原始数据（未过滤的）
cat("加载原始数据...\n")
raw_file <- paste0(sample_name, "_raw.rds")
if(file.exists(raw_file)) {
  seurat_obj <- readRDS(raw_file)
} else {
  # 如果原始数据不存在，从10X重新读取
  data_paths <- list(
    gan56 = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/gan56/gan_56_filtered_feature_bc_matrix",
    liver3 = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/liver_samples/25090559_Liver_3_expression/Liver_3_filtered_feature_bc_matrix",
    liver8 = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/liver_samples/25090559_Liver_8_expression/Liver_8_filtered_feature_bc_matrix",
    bla_liver = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/liver_samples/25090559_Bla_liver_1_expression/Bla_liver_1_filtered_feature_bc_matrix",
    con_liver = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/liver_samples/25090559_Con_liver_1_expression/Con_liver_1_filtered_feature_bc_matrix"
  )
  data_dir <- data_paths[[sample_name]]
  data <- Read10X(data_dir)
  seurat_obj <- CreateSeuratObject(counts = data, project = sample_name)
}

cat("总基因数:", nrow(seurat_obj), "\n")
cat("总细胞数:", ncol(seurat_obj), "\n\n")

# ==================== 基因名分析 ====================
cat("========== 基因名分析 ==========\n")
gene_names <- rownames(seurat_obj)

# 统计不同类型基因
has_grch38 <- sum(grepl("^GRCh38-", gene_names))
has_mm10 <- sum(grepl("^mm10-", gene_names))
has_ensembl <- sum(grepl("^ENSG", gene_names))
has_standard <- sum(!grepl("^GRCh38-|^mm10-|^ENSG", gene_names))

cat("GRCh38-前缀基因数:", has_grch38, "\n")
cat("mm10-前缀基因数:", has_mm10, "\n")
cat("ENSG前缀基因数:", has_ensembl, "\n")
cat("标准基因名数:", has_standard, "\n\n")

# 查看前50个基因名
cat("前50个基因名:\n")
print(head(gene_names, 50))

# 检查是否有蛋白编码基因
protein_coding_keywords <- c("^GRCh38-[A-Z0-9]+$", "^mm10-[A-Z0-9]+$", 
                              "ALB$|AFP$|CYP|CD[0-9]+|PECAM|VWF|KDR|PDGF")
has_protein_coding <- sum(grepl(paste(protein_coding_keywords, collapse="|"), gene_names, ignore.case=TRUE))

cat("\n可能的蛋白编码基因数:", has_protein_coding, "\n")

# ==================== 尝试提取蛋白编码基因 ====================
cat("\n========== 尝试提取蛋白编码基因 ==========\n")

# 方法1: 去除非编码RNA特征
# 常见的非编码RNA模式
noncoding_patterns <- c("MIR", "FAM", "LINC", "AL[0-9]", "AP[0-9]", "AC[0-9]", 
                         "LOC", "RNU", "SNORD", "SCARNA", "RPP", "RPL", "RPS")

# 保留可能为蛋白编码的基因
keep_pattern <- "ALB|AFP|CYP|CD[0-9]+|PECAM|VWF|KDR|PDGF|HNF4|TTR|APOA|HBA|HBB|IGH|IGK|IGL|TCR|TRB|TRA"

protein_coding_idx <- grepl(keep_pattern, gene_names, ignore.case=TRUE) & 
                      !grepl(paste(noncoding_patterns, collapse="|"), gene_names, ignore.case=TRUE)

cat("蛋白编码基因候选数:", sum(protein_coding_idx), "\n")

if(sum(protein_coding_idx) > 0) {
  protein_genes <- gene_names[protein_coding_idx]
  cat("示例蛋白编码基因:\n")
  print(head(protein_genes, 30))
}

# ==================== 尝试从原始数据中获取完整基因名 ====================
cat("\n========== 尝试从原始数据获取完整基因信息 ==========\n")

# 检查feature文件
data_paths <- list(
  gan56 = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/gan56/gan_56_filtered_feature_bc_matrix",
  liver3 = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/liver_samples/25090559_Liver_3_expression/Liver_3_filtered_feature_bc_matrix",
  liver8 = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/liver_samples/25090559_Liver_8_expression/Liver_8_filtered_feature_bc_matrix",
  bla_liver = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/liver_samples/25090559_Bla_liver_1_expression/Bla_liver_1_filtered_feature_bc_matrix",
  con_liver = "/mnt/BioDisk/wangmingyuan/SGS10_analysis/data/liver_samples/25090559_Con_liver_1_expression/Con_liver_1_filtered_feature_bc_matrix"
)

data_dir <- data_paths[[sample_name]]
features_file <- file.path(data_dir, "features.tsv.gz")

if(file.exists(features_file)) {
  cat("读取features.tsv.gz...\n")
  features <- read.table(features_file, header=FALSE, sep="\t", stringsAsFactors=FALSE)
  colnames(features) <- c("gene_id", "gene_name", "feature_type")
  
  cat("Features文件维度:", dim(features), "\n")
  cat("Features列名:", colnames(features), "\n")
  cat("\n前20个features:\n")
  print(head(features, 20))
  
  # 查看feature_type分布
  if("feature_type" %in% colnames(features)) {
    cat("\nFeature类型分布:\n")
    print(table(features$feature_type))
    
    # 筛选Gene Expression类型
    if("Gene Expression" %in% features$feature_type) {
      gene_features <- features[features$feature_type == "Gene Expression", ]
      cat("\nGene Expression基因数:", nrow(gene_features), "\n")
      
      # 使用gene_name作为标准基因名
      if("gene_name" %in% colnames(gene_features)) {
        standard_genes <- gene_features$gene_name
        # 去除空值
        standard_genes <- standard_genes[standard_genes != ""]
        cat("有效基因名数:", length(standard_genes), "\n")
        cat("示例基因名:\n")
        print(head(standard_genes, 30))
        
        # 检查是否有肝脏marker
        liver_markers <- c("ALB", "AFP", "CYP3A4", "APOA1", "PECAM1", "CDH5", 
                           "CD68", "CD163", "PDGFRA", "PDGFRB")
        found_markers <- intersect(liver_markers, standard_genes)
        cat("\n找到的肝脏marker:", paste(found_markers, collapse=", "), "\n")
        
        # 如果找到了标准基因名，重新创建Seurat对象
        if(length(standard_genes) > 0 && nrow(seurat_obj) == nrow(features)) {
          cat("\n使用标准基因名重新创建Seurat对象...\n")
          
          # 获取原始计数矩阵
          counts <- GetAssayData(seurat_obj, layer = "counts")
          
          # 只保留Gene Expression的基因
          keep_idx <- features$feature_type == "Gene Expression" & features$gene_name != ""
          if(sum(keep_idx) > 0) {
            counts_filtered <- counts[keep_idx, ]
            rownames(counts_filtered) <- features$gene_name[keep_idx]
            
            # 去重
            dup_idx <- duplicated(rownames(counts_filtered))
            counts_filtered <- counts_filtered[!dup_idx, ]
            
            cat("过滤后基因数:", nrow(counts_filtered), "\n")
            
            # 创建新对象
            seurat_obj_fixed <- CreateSeuratObject(counts = counts_filtered, project = sample_name)
            seurat_obj_fixed[["percent.mt"]] <- PercentageFeatureSet(seurat_obj_fixed, pattern = "^MT-")
            
            # 保存
            saveRDS(seurat_obj_fixed, paste0(sample_name, "_fixed_raw.rds"))
            cat("✓ 修复后的对象已保存\n")
            
            # 继续处理
            cat("\n========== 重新处理修复后的数据 ==========\n")
            
            # 质控
            seurat_obj_fixed <- subset(seurat_obj_fixed, 
                                       subset = nFeature_RNA > 200 & 
                                                nFeature_RNA < 6000 & 
                                                percent.mt < 20)
            cat("质控后细胞数:", ncol(seurat_obj_fixed), "\n")
            
            # 标准化
            seurat_obj_fixed <- NormalizeData(seurat_obj_fixed, verbose = FALSE)
            seurat_obj_fixed <- FindVariableFeatures(seurat_obj_fixed, nfeatures = 2000, verbose = FALSE)
            seurat_obj_fixed <- ScaleData(seurat_obj_fixed, verbose = FALSE)
            seurat_obj_fixed <- RunPCA(seurat_obj_fixed, npcs = 50, verbose = FALSE)
            seurat_obj_fixed <- FindNeighbors(seurat_obj_fixed, dims = 1:30, verbose = FALSE)
            seurat_obj_fixed <- FindClusters(seurat_obj_fixed, resolution = 0.8, verbose = FALSE)
            seurat_obj_fixed <- RunUMAP(seurat_obj_fixed, dims = 1:30, verbose = FALSE)
            
            # 保存处理后的对象
            saveRDS(seurat_obj_fixed, paste0(sample_name, "_fixed_processed.rds"))
            
            # 使用marker手动注释
            source_liver_markers <- function(seurat_obj) {
              liver_markers <- list(
                Hepatocytes = c("ALB", "AFP", "CYP3A4", "APOA1", "TTR", "HNF4A"),
                Endothelial = c("PECAM1", "CDH5", "VWF", "KDR", "CLEC4G"),
                HSCs = c("PDGFRA", "PDGFRB", "DES", "RGS5", "ACTA2"),
                Kupffer = c("CD68", "CD163", "VSIG4", "MARCO", "C1QA"),
                B_cells = c("CD79A", "MS4A1", "CD19"),
                T_cells = c("CD3D", "CD3E", "CD4", "CD8A"),
                NK_cells = c("NKG7", "KLRB1", "KLRD1"),
                Monocytes = c("CD14", "FCGR3A", "LYZ")
              )
              
              expr_data <- GetAssayData(seurat_obj, layer = "data")
              cluster_ids <- Idents(seurat_obj)
              cluster_levels <- levels(cluster_ids)
              
              score_matrix <- matrix(0, nrow = length(cluster_levels), 
                                     ncol = length(liver_markers))
              rownames(score_matrix) <- cluster_levels
              colnames(score_matrix) <- names(liver_markers)
              
              for(cl in cluster_levels) {
                cells <- WhichCells(seurat_obj, idents = cl)
                if(length(cells) > 100) cells <- sample(cells, min(100, length(cells)))
                
                for(celltype in names(liver_markers)) {
                  markers <- liver_markers[[celltype]]
                  markers_present <- intersect(markers, rownames(expr_data))
                  
                  if(length(markers_present) > 0) {
                    avg_exp <- mean(expr_data[markers_present, cells])
                    score_matrix[cl, celltype] <- avg_exp
                  }
                }
              }
              
              # 为每个cluster分配类型
              assignments <- apply(score_matrix, 1, function(x) {
                max_val <- max(x)
                if(max_val > 0.1) {
                  return(names(which.max(x)))
                } else {
                  return("Unknown")
                }
              })
              
              return(assignments)
            }
            
            assignments <- source_liver_markers(seurat_obj_fixed)
            cat("\n手动注释结果:\n")
            print(table(assignments))
            
            seurat_obj_fixed$celltype <- assignments[as.character(Idents(seurat_obj_fixed))]
            Idents(seurat_obj_fixed) <- seurat_obj_fixed$celltype
            
            # 保存最终结果
            saveRDS(seurat_obj_fixed, paste0(sample_name, "_complete_fixed.rds"))
            
            # 绘制UMAP
            pdf(paste0(sample_name, "_celltypes_fixed.pdf"), width = 12, height = 8)
            print(DimPlot(seurat_obj_fixed, reduction = "umap", label = TRUE, 
                          label.size = 4, repel = TRUE) + 
                    ggtitle(paste0(sample_name, " Cell Types (Fixed)")))
            dev.off()
            
            # 保存细胞类型比例
            prop_df <- data.frame(
              CellType = names(table(Idents(seurat_obj_fixed))),
              Count = as.numeric(table(Idents(seurat_obj_fixed))),
              Percentage = as.numeric(prop.table(table(Idents(seurat_obj_fixed))) * 100)
            )
            write.csv(prop_df, paste0(sample_name, "_celltype_proportions_fixed.csv"), row.names = FALSE)
            
            cat("\n✓ 修复完成！\n")
            cat("细胞类型分布:\n")
            print(prop_df)
          }
        }
      }
    }
  }
} else {
  cat("features.tsv.gz文件不存在\n")
}

cat("\n========== 诊断完成 ==========\n")
