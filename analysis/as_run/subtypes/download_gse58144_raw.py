from pathlib import Path
import json,urllib.request,concurrent.futures,gzip,hashlib
ROOT=Path('<workspace>/subtypes_2026-09-08/reference/GSE58144');OUT=ROOT/'raw';OUT.mkdir(exist_ok=True)
samples=json.loads((ROOT/'sample_metadata.json').read_text());urls=[u.replace('ftp://','https://') for s in samples.values() for u in s['!Sample_supplementary_file']];assert len(urls)==230
def fetch(u):
 p=OUT/u.rsplit('/',1)[1];existing=ROOT/p.name
 if not p.exists():
  if existing.exists():p.write_bytes(existing.read_bytes())
  else:
   data=urllib.request.urlopen(u,timeout=60).read();gzip.decompress(data);p.write_bytes(data)
 return dict(file=p.name,url=u,bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
records=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=6) as ex:
 for r in ex.map(fetch,urls):
  records.append(r)
  if len(records)%20==0:print('RAW_FILES',len(records),'of',len(urls),flush=True)
(ROOT/'raw_manifest.json').write_text(json.dumps(records,indent=2))
print('ALL_RAW_COMPLETE',len(records),flush=True)
