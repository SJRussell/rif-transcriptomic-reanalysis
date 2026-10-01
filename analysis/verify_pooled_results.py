"""Recompute pooled statistics from the supplied per-study results and check the reported tables.

Usage: python verify_pooled_results.py /path/to/results_package
The results package (Supplementary_tables/S1 to S10) is not in this repository.
No model is fitted. Requires numpy, pandas, scipy."""
from pathlib import Path
import sys,json
import numpy as np
import pandas as pd
from scipy.stats import norm,false_discovery_control

if len(sys.argv)!=2:raise SystemExit('Usage: python verify_pooled_results.py /path/to/results_package')
package=Path(sys.argv[1]).resolve()
t=package/'Supplementary_tables'
s1=t/'S1';s2=t/'S2';s4=t/'S4';s6=t/'S6';s7=t/'S7'
st=pd.read_csv(s1/'study_summary.csv').set_index('study');md=pd.read_csv(s1/'master_metadata.csv')
p=pd.read_csv(s2/'pooled/primary.csv').set_index('GeneSymbol')
assert len(st)==13 and len(md)==299 and md['sample'].nunique()==299
assert md.groupby('study').size().to_dict()==st.N.to_dict()
assert (md.condition.str.lower()=='rif').sum()==142
assert len(p)==16008 and p.index.is_unique
z=[]
for study in st.index:
 d=pd.read_csv(s2/f'per_study/{study}_DE.csv').set_index('GeneSymbol')
 assert d.index.is_unique and d.SE.gt(0).all()
 assert d['P.Value'].between(0,1).all() and np.isfinite(d.z).all()
 assert np.all(np.sign(d.z)==np.sign(d.logFC))
 # Recalculate signed Z from the per-study p-value as well as checking pooling.
 expected_z=norm.isf(np.clip(d['P.Value'].to_numpy(),1e-300,1)/2)*np.sign(d.logFC.to_numpy())
 assert np.max(np.abs(expected_z-d.z))<1e-8
 z.append(d.reindex(p.index).z.to_numpy())
z=np.array(z).T;present=np.isfinite(z);w=np.sqrt(st.N_effective.to_numpy())
expected=np.nansum(z*w,axis=1)/np.sqrt(present@(w*w))
zerror=float(np.max(np.abs(expected-p.z)));assert zerror<1e-10
fdr=false_discovery_control(2*norm.sf(np.abs(expected)))
qerror=float(np.max(np.abs(fdr-p.FDR)));assert qerror<1e-10
assert (fdr<=.05).sum()==1840
ma=st.arm.to_numpy()=='microarray'
assert np.all(present[:,ma].sum(1)>=2) and np.all(present[:,~ma].sum(1)>=2)
strict=pd.read_csv(s4/'pooled_strict_all_studies.csv')
assert len(strict)==11349 and set(strict.GeneSymbol)==set(p.index[p.n_studies==13])
members=pd.read_csv(s7/'signature_training_membership.csv')
assert members.groupby('held_out').size().eq(50).all() and not members.duplicated(['held_out','GeneSymbol']).any()
v=pd.read_csv(s7/'signature_validation.csv');assert not v.cohort.str.contains('GSE207362').any()
assert v.loc[(v.n_control<5)|(v.n_rif<5),['CI_low','CI_high']].isna().all().all()
rep=pd.read_csv(s7/'gene_replication_summary.csv')
assert set(rep.cohort)=={'GSE243550','GSE248784'} and rep.discovery_DEGs.eq(1840).all()
fry=pd.read_csv(s6/'programs_fry_and_external.csv');po=pd.read_csv(s6/'programs_meta_fry.csv').set_index('program')
recalculated=[]
for program in po.index:
 d=fry[(fry.program==program)&fry.study.isin(st.index)].set_index('study')
 ww=np.sqrt(st.loc[d.index,'N_effective'].to_numpy())
 recalculated.append(np.sum(ww*d.z.to_numpy())/np.sqrt(np.sum(ww*ww)))
program_error=float(np.max(np.abs(np.array(recalculated)-po.z)))
assert program_error<1e-10 and len(po)==57
assert np.max(np.abs(false_discovery_control(2*norm.sf(np.abs(recalculated)))-po.FDR))<1e-10
assert po.FDR.le(.05).sum()==4
quality=pd.read_csv(s1/'QC/GSE248784_donor_qc.csv')
assert len(quality)==21 and quality.GSM.nunique()==21
assert quality.cells_retained.sum()==188633 and quality.cells_retained.le(quality.cells_deposited).all()
assert not any('invalidated_GSE207362_subset' in str(x) for x in package.rglob('*'))
print(json.dumps(dict(status='PASS',discovery_studies=13,samples=299,primary_genes=16008,DEGs=1840,strict_genes=11349,fixed_programs=57,significant_programs=4,external_donors=21,max_abs_gene_Z_error=zerror,max_abs_gene_BH_error=qerror,max_abs_program_Z_error=program_error),indent=2))

pr=pd.read_csv(t/'S8/predictions.csv');pr=pr[(pr.variant=='released_scale')&(pr.role=='outside_published_development')]
assert len(pr)==126 and (pr.condition=='rif').sum()==69 and pr.study.nunique()==10
clinical=pd.read_csv(t/'S8/clinical_enrichment_by_model_output.csv');clinical=clinical[clinical.scope=='outside_published_development'].set_index('state')
for state in ['RIF1','RIF2','Either_RIF']:
 numerator=denominator=0.
 for _,d in pr.groupby('study'):
  positive=d.pred.isin(['RIF1','RIF2']) if state=='Either_RIF' else d.pred.eq(state)
  case=d.condition.eq('rif');a=sum(case&positive);b=sum(case&~positive);c=sum(~case&positive);dd=sum(~case&~positive)
  numerator+=a*dd/len(d);denominator+=b*c/len(d)
 assert abs(numerator/denominator-clinical.loc[state,'odds_ratio'])<1e-10
assert pr[pr.condition=='rif'].pred.isin(['RIF1','RIF2']).sum()==25
assert pr[pr.condition=='control'].pred.isin(['RIF1','RIF2']).sum()==22
print('PASS: frozen model denominators and independent Mantel-Haenszel odds ratios')

# Reporting checks added on 9 September; no models or testing families are altered.
lineage=pd.read_csv(t/'S9/program_tests.csv')
strict=lineage[(lineage.variant=='strict')&(lineage.FDR_global<=.05)]
assert len(strict)==4 and strict.Direction.eq('Down').all()
assert set(zip(strict.celltype,strict.program))=={('T_cell','HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION'),('Myeloid','HALLMARK_ANGIOGENESIS'),('Myeloid','HALLMARK_MYC_TARGETS_V2'),('Myeloid','HALLMARK_UNFOLDED_PROTEIN_RESPONSE')}
assert not ((lineage.variant=='primary')&(lineage.FDR_global<=.05)).any()
mapping=pd.read_csv(t/'S1/GSE207362_published_subset_mapping.csv')
assert len(mapping)==38 and mapping.count_column.nunique()==38
assert mapping.publication_group.value_counts().to_dict()=={'control':26,'rif':12}
assert mapping.role.eq('second_external_evaluation').all()
expected_control={2,3,5,6,7,9,10,11,12,13,16,17,18,20,21,23,24,26,31,32,33,35,36,38,39,40}
expected_rif={1,4,5,6,7,9,14,17,19,21,23,25}
assert set(mapping[mapping.publication_group=='control'].publication_index)==expected_control
assert set(mapping[mapping.publication_group=='rif'].publication_index)==expected_rif
for r in mapping.itertuples():assert r.count_column==('Endo_Succeed_' if r.publication_group=='control' else 'Endo_Failed_')+str(r.publication_index)
slopes=pd.read_csv(t/'S6/healthy_cycle_programs.csv');scores=pd.read_csv(t/'S6/healthy_cycle_program_scores.csv')
assert len(slopes)==7 and len(scores)>0
g=pd.read_csv(s2/'per_study/GSE111974_DE.csv')
assert len(g)==20128 and g['adj.P.Val'].le(.05).sum()==5324
arms=pd.read_csv(t/'S5/interferon_gamma_arm_comparison.csv').set_index('analysis')
assert arms.loc['primary','NES']<0 and arms.loc['microarray','NES']<0 and arms.loc['rnaseq','NES']>0
index=pd.read_csv(t/'File_index.csv')
actual={str(p.relative_to(package)).replace('\\','/') for p in t.rglob('*') if p.is_file() and p.name not in ['File_index.csv','README.md']}
assert set(index.file)==actual
print('PASS: four strict-annotation programs, published GSE207362 subset mapping, seven LH-day slopes, GSE111974 context, FGSEA arm directions and complete supplement index')

# Second external evaluation (GSE207362).
s=t/'S10'
f=pd.read_csv(s/'GSE207362_external_frozen_signature.csv').iloc[0]
g=pd.read_csv(s/'GSE207362_external_gene_replication_summary.csv').iloc[0]
p=pd.read_csv(s/'GSE207362_external_program_replication_summary.csv').iloc[0]
pred=pd.read_csv(s/'GSE207362_external_metarif_predictions.csv')
pred=pred[pred.variant.eq('released_scale')]
assert (int(f.n_rif),int(f.n_control))==(12,26)
assert abs(float(f.AUC)-0.349358974358974)<1e-12
assert (int(g.testable),int(g.same_direction),int(g.FDR))==(1778,616,0)
assert (int(p.programs),int(p.FDR_two_sided),int(p.FDR_directional))==(57,0,0)
assert len(pred)==38 and pred['sample'].is_unique
assert pred.condition.value_counts().to_dict()=={'control':26,'rif':12}
assigned=pred.pred.isin(['RIF1','RIF2'])
assert assigned[pred.condition.eq('rif')].sum()==3
assert assigned[pred.condition.eq('control')].sum()==3
primary=pd.read_csv(t/'S1/master_metadata.csv')
assert not primary.study.str.contains('GSE207362').any()

print('PASS: GSE207362 second external evaluation (38 libraries, frozen-score AUC, no directional FDR replication, excluded from discovery)')
