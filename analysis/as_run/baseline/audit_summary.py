from pathlib import Path
import sys,json,itertools
sys.path.insert(0,'<scratch>/pydeps')
import pandas as pd
import numpy as np
from scipy.stats import spearmanr
ROOT=Path('<workspace>/subtypes_2026-09-08/corrected_baseline')
st=pd.read_csv(ROOT/'audit/study_summary.csv').set_index('study')
primary=pd.read_csv(ROOT/'pooled/primary.csv');de=primary[primary.FDR<=.05]
ps={s:pd.read_csv(ROOT/f'per_study/{s}_DE.csv').set_index('GeneSymbol') for s in st.index}
cov={s:pd.read_csv(ROOT/f'audit/{s}_coverage.csv').set_index('GeneSymbol') for s in st.index}
allgenes=sorted(set.union(*[set(c.index) for c in cov.values()]))
coverage=pd.DataFrame({'GeneSymbol':allgenes}).set_index('GeneSymbol')
for s in st.index:
 coverage[s+'_represented']=coverage.index.isin(cov[s].index).astype(int)
 coverage[s+'_tested']=coverage.index.isin(ps[s].index).astype(int)
for kind in ['represented','tested']:
 for arm in ['microarray','rnaseq']:
  coverage[f'{arm}_{kind}']=coverage[[s+'_'+kind for s in st.index if st.loc[s,'arm']==arm]].sum(axis=1)
 coverage[f'eligible_{kind}']=(coverage[f'microarray_{kind}']>=2)&(coverage[f'rnaseq_{kind}']>=2)
 coverage[f'all13_{kind}']=(coverage[f'microarray_{kind}']==7)&(coverage[f'rnaseq_{kind}']==6)
coverage.to_csv(ROOT/'audit/gene_coverage_matrix.csv')
pairs=[]
for a,b in itertools.combinations(st.index,2):
 d=ps[a][['logFC','z']].join(ps[b][['logFC','z']],lsuffix='_a',rsuffix='_b',how='inner').dropna()
 d=d[(d.logFC_a!=0)&(d.logFC_b!=0)]
 kind='within_'+st.loc[a,'arm'] if st.loc[a,'arm']==st.loc[b,'arm'] else 'between_platforms'
 pairs.append(dict(study_a=a,study_b=b,type=kind,genes=len(d),direction_agreement=float((np.sign(d.logFC_a)==np.sign(d.logFC_b)).mean()),spearman_logFC=float(spearmanr(d.logFC_a,d.logFC_b).statistic),spearman_Z=float(spearmanr(d.z_a,d.z_b).statistic)))
pd.DataFrame(pairs).to_csv(ROOT/'audit/pairwise_concordance.csv',index=False)
stats=dict(n_studies=len(st),n_samples=int(st.N.sum()),n_rif=int(st.n_rif.sum()),n_control=int(st.n_control.sum()),primary_genes=len(primary),DEGs=len(de),up=int((de.z>0).sum()),down=int((de.z<0).sum()),strict_genes=int(coverage.all13_tested.sum()),represented_2each=int(coverage.eligible_represented.sum()),represented_all13=int(coverage.all13_represented.sum()),microarray_label_union_approved=int((coverage.microarray_represented>0).sum()),rnaseq_label_union_approved=int((coverage.rnaseq_represented>0).sum()),microarray_DE=int((primary.FDR_microarray<=.05).sum()),rnaseq_DE=int((primary.FDR_rnaseq<=.05).sum()),same_arm_direction_all=int((np.sign(primary.z_microarray)==np.sign(primary.z_rnaseq)).sum()),same_arm_direction_DE=int((np.sign(de.z_microarray)==np.sign(de.z_rnaseq)).sum()),median_I2_DE=float(de.I2.median()),median_abs_descriptive_logFC_DE=float(de.logFC_descriptive.abs().median()),DE_I2_over50=int((de.I2>50).sum()))
for label in ['primary','microarray','rnaseq','without_GSE111974','bulk_only']:
 f=pd.read_csv(ROOT/f'pathways/fgsea_{label}.csv');sig=f[f.padj<=.05]
 stats['fgsea_'+label]=dict(tested=len(f),significant=len(sig),positive=int((sig.NES>0).sum()),negative=int((sig.NES<0).sum()))
pd.DataFrame([dict(study=s,represented=len(cov[s]),tested=len(ps[s]),DE=int((ps[s]['adj.P.Val']<=.05).sum())) for s in st.index]).to_csv(ROOT/'audit/study_gene_summary.csv',index=False)
(ROOT/'audit/final_statistics.json').write_text(json.dumps(stats,indent=2))
print(json.dumps(stats,indent=2))
print('PAIRWISE',pd.DataFrame(pairs).groupby('type')[['direction_agreement','spearman_logFC']].mean().to_string())
print('TOP GENES',de.head(12)[['GeneSymbol','z','FDR','I2']].to_string(index=False))
