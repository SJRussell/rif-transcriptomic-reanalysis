"""Fixed broad-lineage marker annotation, blind to clinical labels; donor-level counts.

These are transparent marker assignments, not reproduction of the authors' cell annotations.
All marker panels and thresholds are fixed before inspecting lineage results.
"""
from pathlib import Path
import sys,csv,gzip,tarfile,json
sys.path.insert(0,'<scratch>/pydeps')
import numpy as np
from scipy.io import mmread
from threadpoolctl import threadpool_limits
ROOT=Path('<workspace>/subtypes_2026-09-08');OLD=Path('<workspace>/reanalysis_2026-09-08')
OUT=ROOT/'celltypes';OUT.mkdir(parents=True,exist_ok=True)
markers={
 'Epithelial':['EPCAM','KRT8','KRT18','KRT19','KRT7','CDH1','CLDN4','MUC1'],
 'Stromal':['DCN','LUM','COL1A1','COL1A2','COL3A1','COL6A1','PDGFRA','DPT'],
 'Endothelial':['PECAM1','VWF','EMCN','KDR','ESAM','PLVAP','RAMP2','CDH5'],
 'Perivascular':['RGS5','CSPG4','MCAM','NOTCH3','PDGFRB','MYH11','ACTA2','TAGLN'],
 'T_cell':['CD3D','CD3E','CD3G','TRAC','TRBC1','TRBC2','LCK','CD2'],
 'NK':['GNLY','NKG7','KLRD1','KLRF1','NCR1','PRF1','XCL1','XCL2'],
 'Myeloid':['LST1','LYZ','CTSS','FCER1G','TYROBP','SPI1','AIF1','C1QA'],
 'B_cell':['CD79A','CD79B','MS4A1','CD19','CD22','BANK1'],
 'Plasma':['MZB1','JCHAIN','SDC1','IGKC','IGHG1','XBP1'],
 'Mast':['TPSAB1','TPSB2','CPA3','KIT','MS4A2','HDC']}
(OUT/'annotation_protocol.json').write_text(json.dumps(dict(markers=markers,score='Mean log1p(UMI*10000/total UMI) of represented panel genes; no case labels used',primary='Highest panel score >=0.25, >=2 panel genes detected and >=0.10 gap to runner-up; otherwise Unassigned',strict='Primary assignment with >=3 panel genes detected, >=0.25 score gap and mitochondrial fraction<0.10',T_NK_rule='Cells with >=3 T-panel genes and T score>=0.25 assigned T rather than NK when NK is highest; marked for sensitivity',cell_qc='Same as preceding pooled analysis: 200<=genes<=7000; UMI<40000; MT fraction<0.30',inference='Donor pseudobulk only; at least 50 assigned cells per donor/type, and >=5 donors per clinical group for tests',limitations='No author labels or new doublet classification; unassigned cells are retained in composition denominators; strict sensitivity reports ambiguous-marker impact'),indent=2))
meta=list(csv.DictReader((OLD/'inputs/external/external_metadata.csv').open()));meta={m['sample']:m for m in meta if m['study']=='GSE248784'}
labels=list(markers)+['Unassigned'];all_genes=sorted(set(sum(markers.values(),[])))
summary=[];pseudobulk=[];strictbulk=[];columns=[];dot=[];coverage=[]
with tarfile.open(OLD/'inputs/external/GSE248784_RAW.tar') as tf,gzip.open(OUT/'cell_assignments.csv.gz','wt',newline='') as cf:
 members={m.name:m for m in tf.getmembers() if m.isfile()}
 def rows(m):
  with gzip.open(tf.extractfile(m),'rt') as f:return list(csv.reader(f,delimiter='\t'))
 ff={gsm:rows(next(m for n,m in members.items() if n.startswith(gsm+'_') and n.endswith('features.tsv.gz'))) for gsm in meta}
 common=sorted(set.intersection(*[{r[0] for r in table} for table in ff.values()]));assert len(common)==32058
 common_symbols=set.intersection(*[{r[1] for r in table} for table in ff.values()])
 unavailable={ct:[g for g in panel if g not in common_symbols] for ct,panel in markers.items()}
 (OUT/'reference_unavailable_markers.json').write_text(json.dumps(unavailable,indent=2))
 markers={ct:[g for g in panel if g in common_symbols] for ct,panel in markers.items()}
 assert min(map(len,markers.values()))>=4
 all_genes=sorted(set(sum(markers.values(),[])))
 writer=csv.writer(cf);writer.writerow(['sample','barcode','celltype','strict_retained','top_score','margin','marker_genes_detected','T_NK_override','UMI','mt_fraction'])
 for gsm in sorted(meta):
  features=ff[gsm];ids=[r[0] for r in features];symbols=[r[1] for r in features];ipos={g:i for i,g in enumerate(ids)};reorder=np.array([ipos[g] for g in common])
  matrix=next(m for n,m in members.items() if n.startswith(gsm+'_') and n.endswith('matrix.mtx.gz'))
  bc=rows(next(m for n,m in members.items() if n.startswith(gsm+'_') and n.endswith('barcodes.tsv.gz')))
  with gzip.open(tf.extractfile(matrix),'rb') as f,threadpool_limits(limits=1):x=mmread(f).tocsc()
  x.eliminate_zeros();umi=np.asarray(x.sum(axis=0)).ravel();ngenes=np.diff(x.indptr);mt=np.array([s.upper().startswith('MT-') for s in symbols]);pct=np.asarray(x[mt,:].sum(axis=0)).ravel()/np.maximum(umi,1)
  keep=(ngenes>=200)&(ngenes<=7000)&(umi<40000)&(pct<.30);ki=np.flatnonzero(keep);xx=x[:,keep]
  scores=[];detected=[]
  for ct,panel in markers.items():
   indices=np.flatnonzero(np.isin(symbols,panel));assert len(indices)>=len(panel)
   a=xx[indices,:].toarray();v=np.log1p(a*10000/umi[keep]);scores.append(v.mean(axis=0));detected.append((a>0).sum(axis=0));coverage.append(dict(sample=gsm,celltype=ct,genes_present=len(indices),panel_size=len(panel)))
  scores=np.array(scores);detected=np.array(detected);order=np.argsort(scores,axis=0);top=order[-1];j=np.arange(len(ki));topscore=scores[top,j];margin=topscore-scores[order[-2],j];nd=detected[top,j]
  accepted=(topscore>=.25)&(margin>=.10)&(nd>=2)
  types=np.array(labels,dtype=object)[np.where(accepted,top,len(labels)-1)]
  ti=labels.index('T_cell');ni=labels.index('NK');override=(types=='NK')&(detected[ti]>=3)&(scores[ti]>=.25);types[override]='T_cell'
  strict=accepted&(nd>=3)&(margin>=.25)&(pct[keep]<.10)&~override
  writer.writerows(zip([gsm]*len(ki),[bc[i][0] for i in ki],types,strict.astype(int),topscore,margin,nd,override.astype(int),umi[keep],pct[keep]))
  panel_idx=[symbols.index(g) if g in symbols else -1 for g in all_genes];assert min(panel_idx)>=0
  marker_raw=xx[panel_idx,:].toarray();marker_log=np.log1p(marker_raw*10000/umi[keep])
  for ct in labels:
   select=types==ct;ss=select&strict;count=int(select.sum());scount=int(ss.sum())
   summary.append(dict(sample=gsm,source_title=meta[gsm]['source_title'],condition=meta[gsm]['condition'],celltype=ct,cells=count,total_QC_cells=int(keep.sum()),fraction=count/keep.sum(),strict_cells=scount,strict_total=int(strict.sum()),T_NK_overrides=int((select&override).sum())))
   if count:
    for gi,g in enumerate(all_genes):dot.append(dict(sample=gsm,celltype=ct,gene=g,mean_log_expression=float(marker_log[gi,select].mean()),fraction_detected=float((marker_raw[gi,select]>0).mean()),cells=count))
   if count>=50:
    name=gsm+'__'+ct;columns.append(dict(column=name,sample=gsm,condition=meta[gsm]['condition'],celltype=ct,cells=count,strict_cells=scount))
    pseudobulk.append(np.asarray(xx[:,select].sum(axis=1)).ravel().astype(np.int64)[reorder]);strictbulk.append(np.asarray(xx[:,ss].sum(axis=1)).ravel().astype(np.int64)[reorder])
  print('ANNOTATED',gsm,'retained',len(ki),'assigned',int((types!='Unassigned').sum()),flush=True)
  del x,xx
for name,data in [('donor_cell_composition.csv',summary),('pseudobulk_metadata.csv',columns),('marker_dotplot_data.csv',dot),('marker_coverage.csv',coverage)]:
 with (OUT/name).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=data[0]);w.writeheader();w.writerows(data)
for name,data in [('celltype_counts.csv.gz',pseudobulk),('celltype_counts_strict.csv.gz',strictbulk)]:
 with gzip.open(OUT/name,'wt',newline='') as f:w=csv.writer(f);w.writerow(['ensembl_id']+[c['column'] for c in columns]);w.writerows(zip(common,*data))
print('CELLTYPE_PSEUDOBULK_COMPLETE',len(columns),flush=True)
