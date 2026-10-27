#!/usr/bin/env bash
set -uo pipefail
cd /mnt/BioDisk/wangmingyuan/SGS10_analysis
out=TableS3_clean_20260907_v1
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
nohup /mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript --vanilla "$out/recompute_clean_markers.R" > "$out/recompute.log" 2>&1 < /dev/null &
job=$!
printf '%s\n' "$job" > "$out/recompute.pid"
printf 'MARKER_JOB_STARTED pid=%s\n' "$job"
while kill -0 "$job" 2>/dev/null; do sleep 900; done
if wait "$job"; then
 python3 "$out/validate_clean_markers.py" > "$out/independent_validation.log" 2>&1
 code=$?
 cat "$out/independent_validation.log"
 tail -8 "$out/recompute.log"
 exit "$code"
else
 code=$?
 tail -100 "$out/recompute.log"
 exit "$code"
fi
