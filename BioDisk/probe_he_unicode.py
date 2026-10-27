from pathlib import Path
p=Path('/mnt/BioData/bio')/('SGS10'+chr(0x4f8b))
print('exact_exists',p.exists(),repr(str(p)))
for x in Path('/mnt/BioData/bio').iterdir(): print(repr(x.name),x.is_dir())
if p.exists():
 for f in p.rglob('*aligned_HE*'): print('HE',str(f),f.stat().st_size)
 for f in p.rglob('*.rds'): print('RDS',str(f),f.stat().st_size)
