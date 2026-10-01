.libPaths(c('<r_library>',.libPaths()))
ROOT<-'<workspace>/subtypes_2026-09-08/corrected_baseline'
bundle<-readRDS(file.path(ROOT,'pooled/pooling_inputs.rds'));ps<-bundle$ps;st<-bundle$study_summary
set.seed(20260908)
train_rank<-function(keys) {
 genes<-sort(unique(unlist(lapply(ps[keys],function(d)d$GeneSymbol))))
 zz<-matrix(NA_real_,length(genes),length(keys),dimnames=list(genes,keys))
 for(s in keys)zz[match(ps[[s]]$GeneSymbol,genes),s]<-ps[[s]]$z
 w<-sqrt(st[keys,'N_effective']);present<-!is.na(zz)
 calc<-function(ix)rowSums(sweep(zz[,ix,drop=FALSE],2,w[ix],'*'),na.rm=TRUE)/sqrt(as.vector(present[,ix,drop=FALSE]%*%w[ix]^2))
 ma<-st[keys,'arm']=='microarray';rna<-!ma
 d<-data.frame(GeneSymbol=genes,z=calc(seq_along(keys)),zma=calc(which(ma)),zrna=calc(which(rna)),nma=rowSums(present[,ma,drop=FALSE]),nrna=rowSums(present[,rna,drop=FALSE]))
 d<-d[d$nma>=2&d$nrna>=2&sign(d$zma)==sign(d$zrna),]
 d[order(-abs(d$z)),]
}
auc<-function(score,case){nc<-sum(!case);nr<-sum(case);(sum(rank(score)[case])-nr*(nr+1)/2)/(nr*nc)}
performance<-function(score,case,label,program) {
 a<-auc(score,case);ic<-which(case);inonc<-which(!case)
 boot<-replicate(2000,{i<-c(sample(ic,length(ic),replace=TRUE),sample(inonc,length(inonc),replace=TRUE));auc(score[i],case[i])})
 ci<-quantile(boot,c(.025,.975),names=FALSE)
 # Percentile bootstrap intervals are uninformative with tiny, completely separated groups.
 # Keep the observed AUC, but suppress intervals when either group has fewer than five donors.
 if(min(length(ic),length(inonc))<5||diff(ci)==0)ci<-c(NA_real_,NA_real_)
 # Fixed orientation: higher scores predict RIF; do not flip poor held-out AUCs.
 p<-wilcox.test(score[case],score[!case],alternative='two.sided',exact=FALSE,correct=FALSE)$p.value
 data.frame(cohort=label,program=program,n_rif=sum(case),n_control=sum(!case),AUC=a,CI_low=ci[1],CI_high=ci[2],p=p)
}
rank_score<-function(expr,up,down){ranks<-apply(expr,2,rank,ties.method='average')/nrow(expr);up<-intersect(up,rownames(expr));down<-intersect(down,rownames(expr));stopifnot(length(up)>=5,length(down)>=5);list(score=colMeans(ranks[up,,drop=FALSE])-colMeans(ranks[down,,drop=FALSE]),up=up,down=down)}
results<-list();predictions<-list();members<-list()
tests<-c(names(ps),'GSE243550_external','GSE248784_external','GSE248784_external_mt10','GSE248784_external_shared_reference')
for(s in tests) {
 external<-grepl('_external',s);keys<-if(external)names(ps) else setdiff(names(ps),s)
 d<-train_rank(keys)
 up<-head(d$GeneSymbol[d$z>0],25);down<-head(d$GeneSymbol[d$z<0],25)
 obj<-readRDS(file.path(ROOT,'expression',paste0(s,'.rds')))
 if(!is.null(obj$classifier_expr))obj$expr<-obj$classifier_expr
 sc<-rank_score(obj$expr,up,down);case<-tolower(obj$metadata$condition)=='rif'
 res<-performance(sc$score,case,s,'training_selected_50_gene_signature')
 res$genes_up_measured<-length(sc$up);res$genes_down_measured<-length(sc$down)
 res$design<-if(external)'independent_external' else 'leave_one_cohort_out'
 res$interpretation<-if(s=='GSE58144')'Raw sample/reference rebuild; case-blind batch removal for scoring; held-out labels not used' else if(s=='GSE243550_external')'Repository-labelled technical sensitivity; clinical eligibility unresolved' else if(s=='GSE248784_external')'External cohort with broader RIF and clinical-pregnancy comparator definitions' else 'Held-out labels not used to select genes or score orientation'
 results[[s]]<-res
 predictions[[s]]<-data.frame(cohort=s,sample=obj$metadata$sample,condition=obj$metadata$condition,score=sc$score)
 members[[s]]<-rbind(data.frame(held_out=s,GeneSymbol=up,direction='Up'),data.frame(held_out=s,GeneSymbol=down,direction='Down'))
 cat(s,'AUC',res$AUC,'\n')
}
write.csv(do.call(rbind,results),file.path(ROOT,'sensitivity/signature_validation.csv'),row.names=FALSE)
write.csv(do.call(rbind,predictions),file.path(ROOT,'sensitivity/signature_predictions.csv'),row.names=FALSE)
write.csv(do.call(rbind,members),file.path(ROOT,'sensitivity/signature_training_membership.csv'),row.names=FALSE)
# Two previously published panels, fixed before external testing. Receptivity direction is reversed
# once, explicitly, to test the hypothesis that RIF reflects reduced receptivity; no AUC flipping.
compare_root<-'<results_root>'
benchmark<-list()
for(panel in c('altmae2017','chettiar2024')) {
 fn<-paste0('comparator_',panel,'_meta_signature.tsv');file.copy(file.path(compare_root,fn),file.path(ROOT,'inputs',fn),overwrite=TRUE)
 d<-read.delim(file.path(ROOT,'inputs',fn));lk<-read.csv(file.path(ROOT,'inputs/hgnc_lookup.csv'));d$GeneSymbol<-lk$GeneSymbol[match(d$GeneSymbol,lk$input_id)];d<-d[!is.na(d$GeneSymbol),]
 up<-unique(d$GeneSymbol[tolower(d$Direction)=='up']);down<-unique(d$GeneSymbol[tolower(d$Direction)=='down'])
 if(panel=='altmae2017'){tmp<-up;up<-down;down<-tmp}
 for(s in tests[grepl('_external',tests)]) {
  obj<-readRDS(file.path(ROOT,'expression',paste0(s,'.rds')));sc<-rank_score(obj$expr,up,down)
  rr<-performance(sc$score,tolower(obj$metadata$condition)=='rif',s,panel);rr$genes_up_measured<-length(sc$up);rr$genes_down_measured<-length(sc$down)
  benchmark[[paste(panel,s)]]<-rr
 }
}
write.csv(do.call(rbind,benchmark),file.path(ROOT,'sensitivity/published_panel_benchmarks.csv'),row.names=FALSE)
# External gene-level directional replication: evaluate the fixed primary discovery list without redefining it.
primary<-bundle$primary;de<-primary[primary$FDR<=.05,];rr<-list()
for(s in c('GSE243550','GSE248784')) {
 d<-read.csv(file.path(ROOT,'sensitivity',paste0(s,'_external_DE.csv')))
 m<-merge(de,d,by='GeneSymbol',suffixes=c('_discovery','_external'))
 m$replication_p_two_sided<-m$P.Value;m$replication_FDR_two_sided<-p.adjust(m$P.Value,'BH')
 m$replicated_direction<-sign(m$z_discovery)==sign(m$z_external)
 # One-sided directional p-values explicitly test the pre-specified discovery direction.
 m$replication_p_one_sided<-pnorm(-sign(m$z_discovery)*m$z_external)
 m$replication_FDR_one_sided<-p.adjust(m$replication_p_one_sided,'BH')
 write.csv(m,file.path(ROOT,'sensitivity',paste0('gene_replication_',s,'.csv')),row.names=FALSE)
 rr[[s]]<-data.frame(cohort=s,discovery_DEGs=nrow(de),testable=nrow(m),direction_agreement=sum(m$replicated_direction),nominal_directional=sum(m$replication_p_one_sided<=.05),FDR_directional=sum(m$replication_FDR_one_sided<=.05),FDR_two_sided_same_direction=sum(m$replication_FDR_two_sided<=.05&m$replicated_direction))
}
write.csv(do.call(rbind,rr),file.path(ROOT,'sensitivity/gene_replication_summary.csv'),row.names=FALSE)
print(do.call(rbind,rr));cat('VALIDATION_COMPLETE\n')
