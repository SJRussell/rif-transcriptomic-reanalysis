from pathlib import Path
import json,csv,gzip
P=Path('<workspace>/subtypes_2026-09-08/reference/GSE58144')
s=json.loads((P/'sample_metadata.json').read_text());rows=[]
for gsm,d in s.items():
 ch=dict(x.split(': ',1) for x in d['!Sample_characteristics_ch1'])
 assert ch['condition'] in ('control','recurrent implantation failure after IVF')
 rows.append(dict(sample=gsm,title=d['!Sample_title'][0],condition='Control' if ch['condition']=='control' else 'RIF',batch=ch['batch'],biopsy_day=ch['time of biopsy'].split()[0],age=ch.get('age'),bmi=ch.get('bmi')))
with (P/'donor_metadata.csv').open('w',newline='') as h:
 w=csv.DictWriter(h,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
with gzip.open(P/'GSE58144_rawfiles_dye_orientation.txt.gz','rt') as h:
 orientation={r['raw_data_file'].strip():{k:v.strip() for k,v in r.items()} for r in csv.DictReader(h,delimiter='\t')}
arrays=[]
for gsm,d in s.items():
 for u in d['!Sample_supplementary_file']:
  p=P/'raw'/u.rsplit('/',1)[1];key=p.name.split('_',1)[1].removesuffix('.gz')
  assert p.exists() and key in orientation
  arrays.append(dict(sample=gsm,path=p.as_posix(),**orientation[key]))
assert len(arrays)==230 and len(rows)==115
with (P/'array_targets.csv').open('w',newline='') as h:
 w=csv.DictWriter(h,fieldnames=arrays[0]);w.writeheader();w.writerows(arrays)
print('GSE58144_METADATA',len(rows),'DONORS',len(arrays),'ARRAYS')
