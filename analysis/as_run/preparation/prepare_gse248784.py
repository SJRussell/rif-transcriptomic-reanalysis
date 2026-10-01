"""Donor pseudobulk from the 21 deposited 10x matrices. No cell-level hypothesis tests."""
from pathlib import Path
import sys,csv,gzip,tarfile,json,re,hashlib
sys.path.insert(0,'<scratch>/pydeps')
import numpy as np
from scipy.io import mmread
from threadpoolctl import threadpool_limits
ROOT=Path('<workspace>/reanalysis_2026-09-08')
out=ROOT/'inputs/external'

def read_soft(path):
    samples={};key=None
    with gzip.open(path,'rt',encoding='utf8') as f:
        for line in f:
            if line.startswith('^SAMPLE = '):key=line.strip().split(' = ')[1];samples[key]={}
            elif key and line.startswith('!Sample_title = '):samples[key]['title']=line.strip().split(' = ',1)[1]
    return samples

soft=read_soft(out/'GSE248784_family.soft.gz')
assert len(soft)==21
metadata=[];qc=[];counts=[];strict_counts=[];gene_ids=None
with tarfile.open(out/'GSE248784_RAW.tar') as tf, gzip.open(ROOT/'audit/GSE248784_cell_qc.csv.gz','wt',newline='') as cellf:
    members={m.name:m for m in tf.getmembers() if m.isfile()}
    def rows(member):
        with gzip.open(tf.extractfile(member),'rt') as f:return list(csv.reader(f,delimiter='\t'))
    feature_tables={gsm:rows(next(m for n,m in members.items() if n.startswith(gsm+'_') and n.endswith('features.tsv.gz'))) for gsm in soft}
    gene_ids=sorted(set.intersection(*[{r[0] for r in ff} for ff in feature_tables.values()]))
    assert len(gene_ids)==32058
    cw=csv.writer(cellf);cw.writerow(['GSM','barcode','UMI','genes','mt_fraction','retained','retained_mt10_sensitivity'])
    for gsm in sorted(soft):
        fs={k:next(m for n,m in members.items() if n.startswith(gsm+'_') and n.endswith(k+'.tsv.gz')) for k in ['features','barcodes']}
        matrix=next(m for n,m in members.items() if n.startswith(gsm+'_') and n.endswith('matrix.mtx.gz'))
        features=feature_tables[gsm];barcodes=[r[0] for r in rows(fs['barcodes'])]
        ids=[r[0] for r in features];symbols=[r[1] for r in features]
        assert len(set(ids))==len(ids) and set(gene_ids).issubset(ids)
        positions={g:i for i,g in enumerate(ids)};reorder=np.array([positions[g] for g in gene_ids])
        with gzip.open(tf.extractfile(matrix),'rb') as f,threadpool_limits(limits=1):x=mmread(f).tocsc()
        assert x.shape==(len(ids),len(barcodes)) and np.all(x.data>=0) and np.all(x.data==np.floor(x.data))
        x.eliminate_zeros();umi=np.asarray(x.sum(axis=0)).ravel();ngenes=np.diff(x.indptr)
        mt=np.array([s.upper().startswith('MT-') for s in symbols]);assert 10<=mt.sum()<=50
        pct=np.divide(np.asarray(x[mt,:].sum(axis=0)).ravel(),umi,out=np.zeros(len(umi),dtype=float),where=umi>0)
        keep=(ngenes>=200)&(ngenes<=7000)&(umi<40000)&(pct<.30)
        strict=keep&(pct<.10)
        assert keep.sum()>100
        counts.append(np.asarray(x[:,keep].sum(axis=1)).ravel().astype(np.int64)[reorder])
        strict_counts.append(np.asarray(x[:,strict].sum(axis=1)).ravel().astype(np.int64)[reorder])
        title=soft[gsm]['title'];condition='RIF' if title.startswith('RIF') else 'Control'
        assert title.startswith(('RIF','Normal'))
        metadata.append(dict(study='GSE248784',sample=gsm,condition=condition,arm='rnaseq',tech='pseudobulk',expr_path=(out/'GSE248784_counts.csv').as_posix(),role='external_broader_clinical_cohort',control_definition='First IVF for male-factor infertility followed by ultrasound-confirmed clinical pregnancy',timing='LH+7 protocol with dydrogesterone; natural or induced ovulation',source_title=title,reference_features=len(ids)))
        qc.append(dict(GSM=gsm,title=title,condition=condition,reference_features=len(ids),cells_deposited=x.shape[1],cells_retained=int(keep.sum()),cells_mt10=int(strict.sum()),UMI_retained=int(umi[keep].sum()),median_genes=float(np.median(ngenes[keep])),median_mt=float(np.median(pct[keep])),matrix_member=matrix.name))
        cw.writerows(zip([gsm]*len(barcodes),barcodes,umi,ngenes,pct,keep.astype(int),strict.astype(int)))
        print(gsm,title,'cells',x.shape[1],'retained',keep.sum(),flush=True)
        del x
for filename,cols in [('GSE248784_counts.csv',counts),('GSE248784_mt10_counts.csv',strict_counts)]:
    with (out/filename).open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['GeneSymbol']+[m['sample'] for m in metadata]);w.writerows(zip(gene_ids,*cols))
with (ROOT/'audit/GSE248784_donor_qc.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=qc[0]);w.writeheader();w.writerows(qc)
# Reconcile the repository-only cohort to its independently verified SOFT labels.
existing=list(csv.DictReader((out/'external_metadata.csv').open()))
soft243=read_soft(out/'GSE243550_family.soft.gz')
for m in existing:
    assert m['study']=='GSE243550'
    title=soft243[m['sample']]['title'];cond='RIF' if title.lower().startswith('rif') else 'Control'
    assert title.lower().startswith(('rif','control')) and cond.lower()==m['condition'].lower()
    m.update(role='repository_labelled_technical_sensitivity',control_definition='Repository-labelled control; clinical eligibility of these 40 samples unresolved',timing='Repository design: natural cycle, 5 days after ovulation (LH+7); sample field P+5',source_title=title)
metadata=existing+metadata
with (out/'external_metadata.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=metadata[-1]);w.writeheader();w.writerows(metadata)
(ROOT/'audit/GSE248784_processing.json').write_text(json.dumps(dict(source='GSE248784',publication='PMC13159129; PMID 41758694',unit='one donor per GSM',qc='200<=genes<=7000; UMI<40000; mitochondrial UMI fraction<0.30',features='32,058 Ensembl identifiers present in all 21 feature dictionaries; absent annotation is never imputed as zero counts',doublets='Author manual lineage-coexpression doublet removal cannot be reproduced from the deposited matrices without original annotations; no new doublet removal is claimed.',secondary_qc='Primary filters plus mitochondrial fraction<0.10; separate sensitivity restricts to the 17 donors with a shared 33,538-feature reference dictionary',total_cells=sum(q['cells_deposited'] for q in qc),retained_cells=sum(q['cells_retained'] for q in qc),n_rif=11,n_control=10),indent=2))
print('EXTERNAL_PSEUDOBULK_COMPLETE',flush=True)
