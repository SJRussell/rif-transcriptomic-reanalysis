from pathlib import Path
import sys,json,gzip,csv,zipfile
sys.path.insert(0,'<scratch>/pydeps')
import pandas as pd
import numpy as np
from scipy.stats import spearmanr
R=Path('<workspace>/subtypes_2026-09-08')
hg=pd.read_csv(R/'corrected_baseline/inputs/hgnc_complete_set.txt',sep='\t',dtype=str).query("status == 'Approved'")
# Count columns are the second occurrence of the sample name, under Raw data.
md=pd.read_csv(R/'reference/GSE307510/native_biopsy_sample_manifest.csv')
with gzip.open(R/'reference/GSE307510/GSE307510_R25GdP2652_Single_Report.txt.gz','rt') as f:
 rows=csv.reader(f,delimiter='\t');next(rows);header=next(rows);data=list(rows)
keys=md.native_processed_raw_column.tolist()
idx=[[i for i,x in enumerate(header) if x==k] for k in keys]
assert all(len(x)==2 for x in idx)
assert all(17<=x[1]<27 for x in idx)
counts=pd.DataFrame({k:[int(row[ix[1]]) for row in data] for k,ix in zip(keys,idx)},index=[x[1] for x in data])
counts.index.name='source_symbol';counts.to_csv(R/'reference/GSE307510/native_raw_counts.csv.gz',compression='gzip')
print('Lee native libraries',counts.sum().to_dict())
# Repository-to-paper fingerprint. Preserve declared sample order; do not relabel samples.
d=pd.read_csv(R/'nibi_results/Xu2026_unstranded_counts.csv.gz',index_col=0)
ens=hg.dropna(subset=['ensembl_gene_id']);ens=ens[~ens.ensembl_gene_id.duplicated(False)]
mapping=ens.set_index('ensembl_gene_id').symbol
d.index=d.index.str.replace(r'\.\d+$','',regex=True).map(mapping)
d=d.loc[d.index.notna()].groupby(level=0).sum()
fp=pd.read_csv(R/'reference/Xu2026/supplement_expression_fingerprint.csv',index_col=0)
g=fp.index.intersection(d.index);a=fp.loc[g];b=d.loc[g,a.columns]
def rowz(x):
 x=np.log2(x+1);return x.sub(x.mean(axis=1),axis=0).div(x.std(axis=1).replace(0,np.nan),axis=0)
az=rowz(a);bz=rowz(b)
res=[]
for gene in g:
 res.append({'gene':gene,'raw_count_spearman':spearmanr(a.loc[gene],b.loc[gene]).statistic if a.loc[gene].std()>0 and b.loc[gene].std()>0 else np.nan,'same_integer_values':int((a.loc[gene]==b.loc[gene]).sum()),'n_samples':len(a.columns)})
pd.DataFrame(res).to_csv(R/'audit/Xu2026_supplement_gene_fingerprint.csv',index=False)
cross=pd.DataFrame({x:{y:az[x].corr(bz[y]) for y in b.columns} for x in a.columns})
cross.index.name='raw_declared_sample';cross.to_csv(R/'audit/Xu2026_supplement_sample_correlations.csv')
identity=[{'declared_sample':s,'same_sample_correlation':cross.loc[s,s],'highest_correlation_raw_sample':cross[s].idxmax(),'highest_correlation':cross[s].max()} for s in a.columns]
pd.DataFrame(identity).to_csv(R/'audit/Xu2026_supplement_sample_fingerprint.csv',index=False)
print('Xu source fingerprint',pd.DataFrame(identity).to_string(index=False))
modules=[]
for p in (R/'nibi_bundle').rglob('*fastqc.zip'):
 with zipfile.ZipFile(p) as z:
  key=next(x for x in z.namelist() if x.endswith('/summary.txt'))
  for line in z.read(key).decode().splitlines():
   status,module,sample=line.split('\t');modules.append({'archive':str(p.relative_to(R/'nibi_bundle')),'sample':sample,'module':module,'status':status})
pd.DataFrame(modules).to_csv(R/'audit/Xu2026_fastqc_modules.csv',index=False)
print('FASTQC',pd.DataFrame(modules).groupby(['module','status']).size().to_string())
