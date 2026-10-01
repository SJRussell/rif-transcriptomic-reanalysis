from pathlib import Path
import csv,gzip,tarfile,io
root=Path('<workspace>/reanalysis_2026-09-08')
out=root/'inputs/external'
# GSE207362 is on hold: the 26+12 paper subset is not identified among 67 deposited samples.
# Never construct a subset from an unverified sample list.
metadata=[]
with tarfile.open(out/'GSE243550_RAW.tar') as tf:
    members=sorted([m for m in tf if m.name.endswith('.rawcount.txt.gz')],key=lambda m:m.name)
    assert len(members)==40
    genes=None;columns=[];names=[]
    for member in members:
        name=member.name.split('_')[0]
        with gzip.GzipFile(fileobj=tf.extractfile(member)) as gz:
            lines=io.TextIOWrapper(gz)
            first=next(lines)
            if first.startswith('#'): first=next(lines)
            header=first.rstrip().split('\t')
            assert header[0]=='Geneid'
            gene=[];counts=[]
            for row in csv.reader(lines,delimiter='\t'):
                gene.append(row[0]);counts.append(row[-1])
        if genes is None:genes=gene
        else: assert genes==gene
        names.append(name);columns.append(counts)
        metadata.append({'study':'GSE243550','sample':name,'condition':'RIF' if int(name[3:])<=7790832 else 'Control','arm':'rnaseq','tech':'bulk','expr_path':(out/'GSE243550_counts.csv').as_posix(),'role':'repository_labelled_technical_sensitivity','control_definition':'Repository-labelled control; clinical eligibility of these 40 samples unresolved','timing':'Repository design: natural cycle, 5 days after ovulation (LH+7); sample field P+5'})
        print(name,'loaded',flush=True)
with (out/'GSE243550_counts.csv').open('w',newline='') as f:
    w=csv.writer(f);w.writerow(['GeneSymbol']+names)
    for i,g in enumerate(genes): w.writerow([g]+[col[i] for col in columns])
with (out/'external_metadata.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=metadata[0]);w.writeheader();w.writerows(metadata)
print('Prepared',len(metadata),'external samples')
