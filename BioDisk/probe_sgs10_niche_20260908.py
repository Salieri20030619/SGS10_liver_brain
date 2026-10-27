from pathlib import Path
import sys, importlib.util
print('python', sys.executable)
for p in ['numpy','pandas','matplotlib','scipy','PIL','fitz']:
    print(p, bool(importlib.util.find_spec(p)))
for root in [Path('/mnt/BioDisk/bio'), Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis')]:
    for pat in ['**/*aligned_HE*','**/*cell_locations*','**/*inspect*log']:
        for p in root.glob(pat):
            print(p, p.stat().st_size)
