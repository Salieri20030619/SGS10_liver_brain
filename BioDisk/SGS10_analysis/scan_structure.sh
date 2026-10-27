#!/usr/bin/env bash
ROOT=/mnt/BioDisk/wangmingyuan/26_05_NSCLC
OUT="$ROOT/docs/project_structure_$(date +%Y%m%d).txt"
{
echo "# 项目结构扫描  $(date '+%F %T')"
echo "# 根: $ROOT   总占用: $(du -sh "$ROOT" 2>/dev/null | cut -f1)"
echo
echo "## 目录树（深度3，带文件大小；深4的 data/figures 见下方专节）"
find "$ROOT" -maxdepth 3 -not -path '*/.*' -not -path '*/__pycache__/*' -print0 \
 | sort -z | while IFS= read -r -d '' p; do
    rel=${p#$ROOT}; rel=${rel#/}
    if [ -z "$rel" ]; then lvl=0; else lvl=$(( $(tr -cd '/' <<<"$rel" | wc -c) + 1 )); fi
    ind=$(printf '%*s' $((lvl*2)) '')
    if [ -d "$p" ]; then echo "${ind}$(basename "$p")/"
    else echo "${ind}$(basename "$p")  [$(du -h "$p" 2>/dev/null | cut -f1)]"; fi
done

echo; echo "## code/ 脚本清单 (.R/.py/.sh)"
find "$ROOT/code" \( -name '*.R' -o -name '*.py' -o -name '*.sh' \) 2>/dev/null | sort | sed "s|$ROOT/||"

echo; echo "## 大文件 (>50M)"
find "$ROOT" -type f -size +50M -not -path '*/.*' 2>/dev/null -exec du -h {} + | sort -rh | sed "s|$ROOT/||"

echo; echo "## results/ 各子目录文件数"
for d in "$ROOT"/results/*/ "$ROOT"/results/*/*/; do
  [ -d "$d" ] || continue
  echo "  ${d#$ROOT/} : $(find "$d" -maxdepth 1 -type f 2>/dev/null | wc -l) 文件"
done

echo; echo "## 汇总"
echo "  目录数: $(find "$ROOT" -type d -not -path '*/.*' | wc -l)"
echo "  文件数: $(find "$ROOT" -type f -not -path '*/.*' | wc -l)"
echo "  图(pdf/png/svg): $(find "$ROOT" -type f \( -iname '*.pdf' -o -iname '*.png' -o -iname '*.svg' \) -not -path '*/.*' | wc -l)"
echo "  h5ad: $(find "$ROOT" -iname '*.h5ad' | wc -l)   rds: $(find "$ROOT" -iname '*.rds' | wc -l)"
} | tee "$OUT"
echo; echo ">> 已存 $OUT"
