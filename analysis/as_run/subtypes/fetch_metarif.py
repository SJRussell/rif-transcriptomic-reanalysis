from pathlib import Path
import urllib.request, json, hashlib, datetime
ROOT = Path('<workspace>/subtypes_2026-09-08')
dest = ROOT / 'reference' / 'MetaRIF'
dest.mkdir(parents=True, exist_ok=True)
req = urllib.request.Request('https://api.github.com/repos/YJ-STU/MetaRIF/commits/main', headers={'User-Agent':'rif-methods-audit'})
commit = json.load(urllib.request.urlopen(req, timeout=45))['sha']
records = []
for name in ['MetaRIF.R','MetaRIF_classifier_gene_plus.rda','README.md']:
    url = f'https://raw.githubusercontent.com/YJ-STU/MetaRIF/{commit}/{name}'
    data = urllib.request.urlopen(url, timeout=90).read()
    (dest/name).write_bytes(data)
    records.append(dict(file=name, url=url, commit=commit, bytes=len(data), sha256=hashlib.sha256(data).hexdigest()))
(dest/'source_manifest.json').write_text(json.dumps({'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'files':records},indent=2),encoding='utf-8')
print(json.dumps(records,indent=2))
