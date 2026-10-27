# ==================== 修复后的差异分析脚本 ====================
# 文件名: ~/SGS10_analysis/scripts/liver_tumor_diff_analysis_fixed.R
# 用法: Rscript liver_tumor_diff_analysis_fixed.R

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

# 设置工作目录
setwd("/mnt/BioDisk/wangmingyuan/SGS10_analysis")

# 创建输出目录
output_dir <- "diff_analysis_results"
if(!dir.exists(output_dir)) dir.create(output_dir)

# ==================== 定义样本对应关系 ====================
comparisons <- list(
  # gan56 vs 各对照
  gan56_vs_nao56 = list(tumor = "gan56", control = "nao56", type = "brain"),
  gan56_vs_bla_liver = list(tumor = "gan56", control = "bla_liver", type = "liver"),
  gan56_vs_con_liver = list(tumor = "gan56", control = "con_liver", type = "liver"),
  
  # liver3 vs 各对照
  liver3_vs_brain3 = list(tumor = "liver3", control = "brain3", type = "brain"),
  liver3_vs_bla_liver = list(tumor = "liver3", control = "bla_liver", type = "liver"),
  liver3_vs_con_liver = list(tumor = "liver3", control = "con_liver", type = "liver"),
  
  # liver8 vs 各对照
  liver8_vs_brain8 = list(tumor = "liver8", control = "brain8", type = "brain"),
  liver8_vs_bla_liver = list(tumor = "liver8", control = "bla_liver", type = "liver"),
  liver8_vs_con_liver = list(tumor = "liver8", control = "con_liver", type = "liver")
)

# ==================== 加载所有样本 ====================
cat("\n========== 加载所有样本 ==========\n")

samples_to_load <- unique(c(
  "gan56", "liver3", "liver8",
  "nao56", "brain3", "brain8",
  "bla_liver", "con_liver"
))

seurat_objs <- list()
for(s in samples_to_load) {
  cat("加载", s, "...\n")
  file_path <- file.path(s, paste0(s, "_manual_annotated.rds"))
  if(file.exists(file_path)) {
    obj <- readRDS(file_path)
    # 关键修复：JoinLayers
    cat("  JoinLayers...\n")
    obj <- JoinLayers(obj)
    seurat_objs[[s]] <- obj
    cat("  ✓ 细胞数:", ncol(seurat_objs[[s]]), "\n")
  } else {
    cat("  ✗ 文件不存在:", file_path, "\n")
  }
}

# ==================== 区分人源和小鼠源细胞（针对肿瘤样本） ====================
human_markers <- c("ALB", "AFP", "TTR", "HNF4A", "HNF1A", "GPC3", "EPCAM")
mouse_markers <- c("Alb", "Afp", "Ttr", "Hnf4a", "Hnf1a", "Gpc3", "Epcam")

classify_cells <- function(obj, sample_name) {
  cat("\n分类", sample_name, "细胞...\n")
  
  expr <- GetAssayData(obj, layer = "data")
  
  # 计算人源标记基因平均表达
  human_genes <- intersect(human_markers, rownames(expr))
  if(length(human_genes) > 0) {
    human_score <- colMeans(expr[human_genes, , drop = FALSE])
  } else {
    human_score <- rep(0, ncol(obj))
  }
  
  # 计算小鼠源标记基因平均表达
  mouse_genes <- intersect(mouse_markers, rownames(expr))
  if(length(mouse_genes) > 0) {
    mouse_score <- colMeans(expr[mouse_genes, , drop = FALSE])
  } else {
    mouse_score <- rep(0, ncol(obj))
  }
  
  # 分类
  cell_origin <- ifelse(human_score > mouse_score, "Human", 
                         ifelse(mouse_score > human_score, "Mouse", "Unknown"))
  
  obj$cell_origin <- cell_origin
  obj$human_score <- human_score
  obj$mouse_score <- mouse_score
  
  cat("  人源细胞:", sum(cell_origin == "Human"), "\n")
  cat("  小鼠细胞:", sum(cell_origin == "Mouse"), "\n")
  cat("  未知细胞:", sum(cell_origin == "Unknown"), "\n")
  
  return(obj)
}

# 对肿瘤样本进行分类
tumor_samples <- c("gan56", "liver3", "liver8")
for(s in tumor_samples) {
  if(s %in% names(seurat_objs)) {
    seurat_objs[[s]] <- classify_cells(seurat_objs[[s]], s)
  }
}

# ==================== 函数：提取细胞类型并合并 ====================
extract_and_merge <- function(tumor_obj, control_obj, cell_type = NULL, origin = NULL) {
  # 筛选细胞
  if(!is.null(cell_type)) {
    tumor_cells <- WhichCells(tumor_obj, expression = final_celltype == cell_type)
    control_cells <- WhichCells(control_obj, expression = final_celltype == cell_type)
  } else {
    tumor_cells <- colnames(tumor_obj)
    control_cells <- colnames(control_obj)
  }
  
  # 如果是肿瘤样本，可以进一步筛选人源或小鼠源
  if(!is.null(origin) && "cell_origin" %in% colnames(tumor_obj@meta.data)) {
    tumor_cells <- intersect(tumor_cells, 
                              WhichCells(tumor_obj, expression = cell_origin == origin))
  }
  
  cat("  肿瘤细胞数:", length(tumor_cells), "\n")
  cat("  对照细胞数:", length(control_cells), "\n")
  
  if(length(tumor_cells) == 0 || length(control_cells) == 0) {
    return(NULL)
  }
  
  # 创建子集
  tumor_subset <- subset(tumor_obj, cells = tumor_cells)
  control_subset <- subset(control_obj, cells = control_cells)
  
  # 添加分组信息
  tumor_subset$group <- "Tumor"
  control_subset$group <- "Control"
  
  # 合并
  merged <- merge(tumor_subset, control_subset)
  
  # 关键修复：合并后也需要JoinLayers
  merged <- JoinLayers(merged)
  
  return(merged)
}

# ==================== 函数：差异表达分析 ====================
run_diff_analysis <- function(merged_obj, comparison_name, output_dir) {
  cat("\n==========", comparison_name, "==========\n")
  
  # 设置身份
  Idents(merged_obj) <- merged_obj$group
  
  # 检查细胞数
  tumor_count <- sum(merged_obj$group == "Tumor")
  control_count <- sum(merged_obj$group == "Control")
  cat("Tumor细胞数:", tumor_count, "\n")
  cat("Control细胞数:", control_count, "\n")
  
  if(tumor_count < 10 || control_count < 10) {
    cat("⚠️ 细胞数不足，跳过分析\n")
    return(NULL)
  }
  
  # 运行差异分析（使用default assay）
  DefaultAssay(merged_obj) <- "RNA"
  
  # 尝试运行FindMarkers
  markers <- tryCatch({
    FindMarkers(merged_obj,
                ident.1 = "Tumor",
                ident.2 = "Control",
                min.pct = 0.1,
                logfc.threshold = 0.25,
                test.use = "wilcox",
                verbose = TRUE)
  }, error = function(e) {
    cat("错误:", e$message, "\n")
    # 尝试使用另一种方法
    cat("尝试使用MAST...\n")
    tryCatch({
      FindMarkers(merged_obj,
                  ident.1 = "Tumor",
                  ident.2 = "Control",
                  min.pct = 0.1,
                  logfc.threshold = 0.25,
                  test.use = "MAST",
                  verbose = TRUE)
    }, error = function(e2) {
      cat("MAST也失败:", e2$message, "\n")
      return(NULL)
    })
  })
  
  if(is.null(markers) || nrow(markers) == 0) {
    cat("⚠️ 未找到差异基因\n")
    return(NULL)
  }
  
  # 添加基因名
  markers$gene <- rownames(markers)
  
  # 排序
  markers_up <- markers %>% filter(avg_log2FC > 0) %>% arrange(p_val_adj)
  markers_down <- markers %>% filter(avg_log2FC < 0) %>% arrange(p_val_adj)
  
  cat("上调基因数:", nrow(markers_up), "\n")
  cat("下调基因数:", nrow(markers_down), "\n")
  
  # 保存结果
  write.csv(markers, file.path(output_dir, paste0(comparison_name, "_all_genes.csv")), row.names = FALSE)
  write.csv(markers_up, file.path(output_dir, paste0(comparison_name, "_upregulated.csv")), row.names = FALSE)
  write.csv(markers_down, file.path(output_dir, paste0(comparison_name, "_downregulated.csv")), row.names = FALSE)
  
  # 火山图
  if(requireNamespace("EnhancedVolcano", quietly = TRUE)) {
    pdf(file.path(output_dir, paste0(comparison_name, "_volcano.pdf")), width = 10, height = 8)
    p <- EnhancedVolcano::EnhancedVolcano(markers,
                       lab = rownames(markers),
                       x = 'avg_log2FC',
                       y = 'p_val_adj',
                       title = comparison_name,
                       pCutoff = 0.05,
                       FCcutoff = 0.5,
                       pointSize = 2.0,
                       labSize = 3.0)
    print(p)
    dev.off()
  }
  
  # 热图
  if(nrow(markers) >= 20) {
    top_genes <- markers %>% 
      arrange(p_val_adj) %>% 
      head(20) %>% 
      pull(gene)
    
    pdf(file.path(output_dir, paste0(comparison_name, "_heatmap.pdf")), width = 12, height = 10)
    p <- DoHeatmap(merged_obj, features = top_genes, group.by = "group") +
      scale_fill_gradientn(colors = c("blue", "white", "red")) +
      ggtitle(paste0(comparison_name, " - Top 20 DEGs"))
    print(p)
    dev.off()
  }
  
  return(markers)
}

# ==================== 主分析 ====================
cat("\n========== 开始差异分析 ==========\n")

# 定义要分析的细胞类型
cell_types <- c("Hepatocytes", "Kupffer", "Endothelial", "T_cells", "B_cells", "Stellate")

# 创建结果汇总
results_summary <- data.frame()

for(comp_name in names(comparisons)) {
  comp <- comparisons[[comp_name]]
  tumor <- comp$tumor
  control <- comp$control
  comp_type <- comp$type
  
  cat("\n\n========================================\n")
  cat("分析:", tumor, "vs", control, "(", comp_type, ")\n")
  cat("========================================\n")
  
  tumor_obj <- seurat_objs[[tumor]]
  control_obj <- seurat_objs[[control]]
  
  if(is.null(tumor_obj) || is.null(control_obj)) {
    cat("⚠️ 样本不存在\n")
    next
  }
  
  # 对于肿瘤样本，分析小鼠源细胞（肿瘤微环境）
  if("cell_origin" %in% colnames(tumor_obj@meta.data)) {
    cat("\n--- 分析小鼠源细胞（肿瘤微环境）---\n")
    merged_mouse <- extract_and_merge(tumor_obj, control_obj, origin = "Mouse")
    if(!is.null(merged_mouse)) {
      results_mouse <- run_diff_analysis(merged_mouse, 
                                          paste0(comp_name, "_mouse_origin"), 
                                          output_dir)
      if(!is.null(results_mouse)) {
        results_summary <- rbind(results_summary, 
                                 data.frame(
                                   Comparison = comp_name,
                                   CellType = "All_Mouse",
                                   Up = sum(results_mouse$avg_log2FC > 0 & results_mouse$p_val_adj < 0.05, na.rm = TRUE),
                                   Down = sum(results_mouse$avg_log2FC < 0 & results_mouse$p_val_adj < 0.05, na.rm = TRUE)
                                 ))
      }
    }
    
    # 分析人源细胞（肿瘤细胞）
    cat("\n--- 分析人源细胞（肿瘤细胞）---\n")
    merged_human <- extract_and_merge(tumor_obj, control_obj, origin = "Human")
    if(!is.null(merged_human)) {
      results_human <- run_diff_analysis(merged_human, 
                                          paste0(comp_name, "_human_origin"), 
                                          output_dir)
      if(!is.null(results_human)) {
        results_summary <- rbind(results_summary, 
                                 data.frame(
                                   Comparison = comp_name,
                                   CellType = "All_Human",
                                   Up = sum(results_human$avg_log2FC > 0 & results_human$p_val_adj < 0.05, na.rm = TRUE),
                                   Down = sum(results_human$avg_log2FC < 0 & results_human$p_val_adj < 0.05, na.rm = TRUE)
                                 ))
      }
    }
  }
  
  # 按细胞类型分析
  for(ct in cell_types) {
    cat("\n--- 分析细胞类型:", ct, "---\n")
    results <- tryCatch({
      merged <- extract_and_merge(tumor_obj, control_obj, cell_type = ct)
      if(!is.null(merged)) {
        run_diff_analysis(merged, paste0(comp_name, "_", ct), output_dir)
      } else {
        NULL
      }
    }, error = function(e) {
      cat("错误:", e$message, "\n")
      NULL
    })
    
    if(!is.null(results)) {
      results_summary <- rbind(results_summary, 
                               data.frame(
                                 Comparison = comp_name,
                                 CellType = ct,
                                 Up = sum(results$avg_log2FC > 0 & results$p_val_adj < 0.05, na.rm = TRUE),
                                 Down = sum(results$avg_log2FC < 0 & results$p_val_adj < 0.05, na.rm = TRUE)
                               ))
    }
  }
}

# ==================== 保存结果汇总 ====================
cat("\n========== 保存结果汇总 ==========\n")
write.csv(results_summary, file.path(output_dir, "diff_analysis_summary.csv"), row.names = FALSE)

# 打印汇总
cat("\n\n========== 差异分析结果汇总 ==========\n")
print(results_summary)

cat("\n\n========== 分析完成 ==========\n")
cat("结果保存在:", output_dir, "\n")
