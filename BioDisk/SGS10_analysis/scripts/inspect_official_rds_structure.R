# =============================================================
# inspect_official_rds_structure.R
# Purpose: 摸清楚官方 SeekSpace RDS 的 coords 存在哪、meta.data 有哪些列
# =============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

RDS_PATH <- "/mnt/BioData/bio/SGS10例/SGS25090559-8例样本数据/data/rds/25090559_Liver_3_expression/Liver_3.rds"

cat("Loading", RDS_PATH, "...\n")
obj <- readRDS(RDS_PATH)
cat("Class:", paste(class(obj), collapse = ","), "\n\n")

# -------------------------------------------------------------
# 1. meta.data 结构
# -------------------------------------------------------------
cat("========================================\n")
cat(" 1. meta.data\n")
cat("========================================\n")
cat("nrow:", nrow(obj@meta.data), "\n")
cat("ncol:", ncol(obj@meta.data), "\n")
cat("Columns:\n")
print(colnames(obj@meta.data))
cat("\nFirst 3 rows:\n")
print(head(obj@meta.data, 3))
cat("\n")

# -------------------------------------------------------------
# 2. 对象整体 slot
# -------------------------------------------------------------
cat("========================================\n")
cat(" 2. Seurat object slots\n")
cat("========================================\n")
cat("slotNames:\n")
print(slotNames(obj))
cat("\n")

# -------------------------------------------------------------
# 3. @images
# -------------------------------------------------------------
cat("========================================\n")
cat(" 3. @images\n")
cat("========================================\n")
cat("Names:", names(obj@images), "\n")
if (length(obj@images) > 0) {
  img1 <- obj@images[[1]]
  cat("Class of @images[[1]]:", paste(class(img1), collapse = ","), "\n")
  cat("slotNames:\n")
  print(slotNames(img1))
  
  # 常见 SpatialImage 结构
  if ("coordinates" %in% slotNames(img1)) {
    cat("\n@images[[1]]@coordinates (first 3 rows):\n")
    print(head(img1@coordinates, 3))
    cat("nrow:", nrow(img1@coordinates), "\n")
  }
  if ("boundaries" %in% slotNames(img1)) {
    cat("\n@images[[1]]@boundaries (names):\n")
    print(names(img1@boundaries))
  }
}
cat("\n")

# -------------------------------------------------------------
# 4. @misc
# -------------------------------------------------------------
cat("========================================\n")
cat(" 4. @misc\n")
cat("========================================\n")
cat("Names:", names(obj@misc), "\n")

if (!is.null(obj@misc$info)) {
  cat("\n@misc$info length:", length(obj@misc$info), "\n")
  if (length(obj@misc$info) > 0) {
    info1 <- obj@misc$info[[1]]
    cat("Names of @misc$info[[1]]:\n")
    print(names(info1))
    
    # 查每个字段的类型和长度
    for (nm in names(info1)) {
      x <- info1[[nm]]
      type_str <- paste(class(x), collapse = ",")
      if (is.data.frame(x) || is.matrix(x)) {
        cat(sprintf("  %-20s %s (nrow=%d, ncol=%d)\n",
                    nm, type_str, nrow(x), ncol(x)))
        if (nm %in% c("cell_loc", "spot_loc", "coords", "location")) {
          cat("    first 3 rows:\n")
          print(head(x, 3))
        }
      } else if (is.vector(x) && !is.list(x)) {
        cat(sprintf("  %-20s %s (length=%d)\n",
                    nm, type_str, length(x)))
      } else {
        cat(sprintf("  %-20s %s\n", nm, type_str))
      }
    }
  }
}
cat("\n")

# -------------------------------------------------------------
# 5. counts 矩阵的列名（spot barcode）
# -------------------------------------------------------------
cat("========================================\n")
cat(" 5. Counts matrix column names (first 5)\n")
cat("========================================\n")
counts <- tryCatch(
  GetAssayData(obj, assay = "RNA", layer = "counts"),
  error = function(e) GetAssayData(obj, assay = "RNA", slot = "counts")
)
cat("ncol:", ncol(counts), "\n")
cat("First 5 colnames (barcodes):\n")
print(head(colnames(counts), 5))
cat("\nrownames of meta.data (first 5):\n")
print(head(rownames(obj@meta.data), 5))
cat("\nAre they identical?", identical(colnames(counts), rownames(obj@meta.data)), "\n")

cat("\nDONE.\n")
