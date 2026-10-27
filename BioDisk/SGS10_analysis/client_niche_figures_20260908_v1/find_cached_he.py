from pathlib import Path
p=Path('/mnt/BioDisk/wangmingyuan/SGS10_analysis')
for d in p.rglob('*'):
 if d.is_file() and d.suffix.lower() in ['.png','.jpg','.tif','.tiff','.rds']:
  n=d.name.lower()
  if ('he_only' in n or 'he_image' in n or 'h&e' in n or n.startswith('he_') or n in ['brain_8.rds','gan56.rds','gan_56.rds']): print(d)
