.libPaths(c('<r_library>',.libPaths()))
suppressPackageStartupMessages({library(edgeR);library(limma)})
options(warn=1)
ROOT<-'<workspace>/subtypes_2026-09-08';OLD<-'<workspace>/reanalysis_2026-09-08';OUT<-file.path(ROOT,'celltypes')
dir.create(file.path(OUT,'DE'),showWarnings=FALSE)
lookup<-read.csv(file.path(OLD,'inputs/hgnc_lookup.csv'));sym<-setNames(lookup$GeneSymbol,lookup$input_id)
meta<-read.csv(file.path(OUT,'pseudobulk_metadata.csv'));comp<-read.csv(file.path(OUT,'donor_cell_composition.csv'))
features<-read.csv(file.path(ROOT,'audit/model_feature_hgnc.csv'))$symbol
panels<-jsonlite::fromJSON(file.path(OUT,'annotation_protocol.json'))$markers
sets<-readRDS(file.path(OLD,'pathways/fixed_programs.rds'));excluded<-unique(c(features,unlist(panels)))
sets<-lapply(sets,setdiff,y=excluded)
results<-scoredata<-summaries<-list()
for(variant in c('primary','strict')){
 filename<-if(variant=='primary')'celltype_counts.csv.gz' else 'celltype_counts_strict.csv.gz'
 d<-read.csv(gzfile(file.path(OUT,filename)),check.names=FALSE);m<-as.matrix(d[,-1]);ids<-unname(sym[sub('[.][0-9]+$','',d[[1]])]);good<-!is.na(ids)&ids!=''
 m<-rowsum(m[good,,drop=FALSE],group=ids[good],reorder=TRUE)
 for(ct in unique(meta$celltype)){
  md<-meta[meta$celltype==ct & if(variant=='strict')meta$strict_cells>=50 else meta$cells>=50,]
  md$group<-factor(tolower(md$condition),levels=c('control','rif'));nn<-table(md$group)
  if(any(nn<5)){summaries[[paste(variant,ct)]]<-data.frame(variant=variant,celltype=ct,n_control=nn[1],n_rif=nn[2],genes=NA,DE_FDR05=NA,status='fewer than 5 donors in a group');next}
  mm<-m[,match(md$column,colnames(m)),drop=FALSE];design<-model.matrix(~group,data=md);y<-DGEList(mm);keep<-filterByExpr(y,design=design);y<-calcNormFactors(y[keep,,keep.lib.sizes=FALSE]);v<-voom(y,design,plot=FALSE)
  fit<-eBayes(lmFit(v,design),robust=TRUE);tt<-topTable(fit,coef=2,number=Inf,sort.by='none');tt$GeneSymbol<-rownames(tt)
  con<-gzfile(file.path(OUT,'DE',paste0(ct,'_',variant,'.csv.gz')),'wt');write.csv(tt,con,row.names=FALSE);close(con)
  sums<-data.frame(variant=variant,celltype=ct,n_control=nn[1],n_rif=nn[2],genes=nrow(tt),DE_FDR05=sum(tt$adj.P.Val<=.05),status='tested');summaries[[paste(variant,ct)]]<-sums
  index<-ids2indices(sets,rownames(v$E),remove.empty=TRUE);index<-index[lengths(index)>=5];f<-fry(v,index=index,design=design,contrast=2,sort=FALSE)
  results[[paste(variant,ct)]]<-data.frame(variant=variant,celltype=ct,program=rownames(f),n_control=nn[1],n_rif=nn[2],genes=f$NGenes,Direction=f$Direction,p=f$PValue,p_mixed=f$PValue.Mixed)
  z<-t(scale(t(v$E)));z[!is.finite(z)]<-0
  for(p in names(index))scoredata[[paste(variant,ct,p)]]<-data.frame(variant=variant,celltype=ct,program=p,sample=md$sample,condition=tolower(md$condition),score=colMeans(z[index[[p]],,drop=FALSE]),cells=if(variant=='strict')md$strict_cells else md$cells)
  saveRDS(list(expr=v$E,metadata=md,design=design,counts=mm),file.path(OUT,paste0(ct,'_',variant,'.rds')))
  cat('CELLTYPE_FIT',variant,ct,'n',paste(nn,collapse='/'),'DE',sums$DE_FDR05,'\n');flush.console()
 }
}
r<-do.call(rbind,results);r$FDR_global<-ave(r$p,r$variant,FUN=function(x)p.adjust(x,'BH'));r$FDR_within_celltype<-ave(r$p,interaction(r$variant,r$celltype),FUN=function(x)p.adjust(x,'BH'))
write.csv(r,file.path(OUT,'program_tests.csv'),row.names=FALSE);write.csv(do.call(rbind,scoredata),file.path(OUT,'program_scores.csv'),row.names=FALSE);write.csv(do.call(rbind,summaries),file.path(OUT,'DE_summary.csv'),row.names=FALSE)
# Composition comparisons use one fraction per donor. Cells are never independent replicates.
cc<-list();set.seed(20260908)
for(v in c('primary','strict'))for(ct in unique(comp$celltype)){
 d<-comp[comp$celltype==ct,];y<-tolower(d$condition)=='rif';value<-if(v=='primary')d$fraction else d$strict_cells/pmax(d$strict_total,1)
 if(v=='strict'&&ct=='Unassigned')next
 a<-value[y];b<-value[!y];p<-wilcox.test(a,b,exact=FALSE)$p.value
 ci<-quantile(replicate(2000,mean(sample(a,length(a),TRUE))-mean(sample(b,length(b),TRUE))),c(.025,.975))
 cc[[paste(v,ct)]]<-data.frame(variant=v,celltype=ct,n_rif=length(a),n_control=length(b),mean_RIF=mean(a),mean_control=mean(b),difference=mean(a)-mean(b),CI_low=ci[1],CI_high=ci[2],p=p)
}
c<-do.call(rbind,cc);c$FDR<-ave(c$p,c$variant,FUN=function(x)p.adjust(x,'BH'));write.csv(c,file.path(OUT,'composition_tests.csv'),row.names=FALSE)
# Association with fixed tissue axes: cell fractions + documented reference dictionary sensitivity later.
aa<-read.csv(file.path(ROOT,'results/axes_and_predictions.csv'));aa<-aa[aa$study=='GSE248784_external'&aa$variant=='exclude_model_genes',]
wide<-reshape(comp[,c('sample','celltype','fraction')],idvar='sample',timevar='celltype',direction='wide');names(wide)<-sub('fraction[.]','',names(wide));d<-merge(aa,wide,by='sample')
d$immune_fraction<-d$T_cell+d$NK+d$Myeloid+d$B_cell+d$Plasma+d$Mast
rr<-list()
for(ax in c('immune_axis','metabolic_axis')){
 y<-d[[ax]];group<-as.numeric(d$condition=='rif')
 f0<-lm(y~group);f1<-lm(y~group+immune_fraction+Epithelial,data=d)
 for(term in c('group','immune_fraction','Epithelial')){
  f<-summary(f1)$coefficients;if(!term%in%rownames(f))next
  rr[[paste(ax,term)]]<-data.frame(axis=ax,term=term,estimate=f[term,1],SE=f[term,2],p=f[term,4],adjusted_R2=summary(f1)$adj.r.squared,unadjusted_RIF_effect=coef(f0)[2],unadjusted_RIF_p=summary(f0)$coefficients[2,4])
 }
}
write.csv(do.call(rbind,rr),file.path(OUT,'tissue_axis_composition_models.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(ROOT,'logs/celltypes_sessionInfo.txt'));cat('CELLTYPE_ANALYSIS_COMPLETE\n')
