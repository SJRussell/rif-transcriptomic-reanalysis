from pathlib import Path
import csv, tarfile, gzip, shutil, json, hashlib
from collections import defaultdict, Counter
ROOT=Path('<workspace>/reanalysis_2026-09-08')
SOURCE=Path('<workspace>/sdrive_audit_20260908')
for part in ['inputs','audit','per_study','sensitivity','pooled','figures','logs','expression','pathways']:
    (ROOT/part).mkdir(parents=True,exist_ok=True)
hgnc=Path('<scratch>/hgnc_complete_set.txt')
shutil.copy2(hgnc,ROOT/'inputs/hgnc_complete_set.txt')
records=list(csv.DictReader(hgnc.open(encoding='utf-8-sig'),delimiter='\t'))
approved={r['symbol']:r for r in records if r['status']=='Approved'}
maps={k:defaultdict(set) for k in ['prev_symbol','alias_symbol','ensembl_gene_id','entrez_id']}
for sym,r in approved.items():
    for k,d in maps.items():
        for a in r[k].split('|'):
            if a: d[a].add(sym)
rows=[]
keys=set(approved)
for d in maps.values(): keys.update(d)
for key in sorted(keys):
    target='';status='unmapped'
    if key in approved: target=key;status='approved'
    else:
        for kind,d in maps.items():
            candidates=d.get(key,set())
            if candidates:
                if len(candidates)==1: target=next(iter(candidates));status=kind
                else: status='ambiguous_'+kind
                break
    rows.append({'input_id':key,'GeneSymbol':target,'mapping_status':status,'locus_group':approved[target]['locus_group'] if target else ''})
with (ROOT/'inputs/hgnc_lookup.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
rawdir=ROOT/'inputs/GSE71331_raw';rawdir.mkdir(exist_ok=True)
with tarfile.open(SOURCE/'MICROARRAY_DATA/GSE71331/GSE71331_RAW.tar') as tf:
    for member in tf.getmembers():
        if member.isfile() and member.name.startswith('GSM') and member.name.endswith('.txt.gz'):
            dest=rawdir/Path(member.name).name[:-3]
            if not dest.exists():
                with gzip.GzipFile(fileobj=tf.extractfile(member)) as src, dest.open('wb') as dst: shutil.copyfileobj(src,dst)
assert len(list(rawdir.glob('GSM*.txt')))==12
shutil.copy2(SOURCE/'MICROARRAY_DATA/GSE71331/GSE71331_annotation.txt',ROOT/'inputs/GSE71331_annotation.txt')
md=list(csv.DictReader(Path('<workspace>/rebuild_2026-09/inputs/master_metadata.csv').open()))
bui=list(csv.DictReader((SOURCE/'RNASEQ_DATA/Bui_Data/Bui_metadata.csv').open()))
bui_timing={r['']:r['Timing_LH'] for r in bui}
excluded=[];primary=[]
for row in md:
    reason=''
    row['biopsy_day']=''
    if row['study']=='Bui2024':
        try:row['biopsy_day']=str(float(bui_timing[row['sample']]))
        except ValueError: reason='Biopsy timing unavailable; retained in all-sample unadjusted sensitivity'
    if row['study']=='GSE223672': reason='RIF-only cohort, no within-study control contrast'
    elif row['study']=='GSE250130' and row['condition']=='Control' and row['sample'] not in ['GSE250130_con4','GSE250130_con5','GSE250130_con6']: reason='Control biopsy outside LH+7; retained for timing sensitivity'
    if reason: excluded.append(dict(row,reason=reason))
    else: primary.append(row)
for name,rr in [('master_metadata.csv',primary),('excluded_original_samples.csv',excluded),('original_metadata.csv',md)]:
    with (ROOT/'inputs'/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rr[0]);w.writeheader();w.writerows(rr)
print(json.dumps({'primary_samples':len(primary),'studies':dict(Counter(r['study'] for r in primary)),'hgnc_lookup':len(rows),'raw_arrays':12},indent=2))
