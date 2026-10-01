from pathlib import Path
import urllib.request,gzip,csv,json,concurrent.futures,sys
sys.stdout.reconfigure(encoding='utf-8')
OUT=Path('<workspace>/subtypes_2026-09-08/reference/GSE58144');OUT.mkdir(parents=True,exist_ok=True)
urls={'GSE58144_family.soft.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE58nnn/GSE58144/soft/GSE58144_family.soft.gz','GSE58144_series_matrix.txt.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE58nnn/GSE58144/matrix/GSE58144_series_matrix.txt.gz'}
def get(kv):
 n,u=kv;p=OUT/n
 if not p.exists():p.write_bytes(urllib.request.urlopen(u,timeout=90).read())
 print('DOWNLOADED',n,p.stat().st_size,flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:list(ex.map(get,urls.items()))
samples={};cur=None
with gzip.open(OUT/'GSE58144_family.soft.gz','rt',encoding='utf-8',errors='replace') as f:
 for l in f:
  if l.startswith('^SAMPLE = '):cur=l.strip().split(' = ')[1];samples[cur]={}
  elif l.startswith('^'):cur=None
  elif cur and l.startswith('!Sample_'):
   k,v=l.strip().split(' = ',1);samples[cur].setdefault(k,[]).append(v)
(OUT/'sample_metadata.json').write_text(json.dumps(samples,indent=2))
for gsm in list(samples)[:2]:
 print(gsm,json.dumps({k:v for k,v in samples[gsm].items() if any(t in k for t in ['title','label','source_name','data_processing','supplementary_file','description'])},indent=2))
