.libPaths(c('<r_library>',.libPaths()))
suppressPackageStartupMessages({library(limma);library(edgeR)})
BASE<-'<workspace>/subtypes_2026-09-08';ROOT<-file.path(BASE,'corrected_baseline');set.seed(20260908)
for(e in parse(file.path(ROOT,'scripts/validate_signatures.R')))if(is.call(e)&&identical(e[[1]],as.name('<-'))&&as.character(e[[2]])%in%c('auc','performance','rank_score'))eval(e)
panel<-read.csv(file.path(ROOT,'sensitivity/signature_training_membership.csv'));panel<-panel[panel$held_out=='GSE248784_external',]
programs<-readRDS(file.path(ROOT,'pathways/fixed_programs.rds'));discovery<-read.csv(file.path(ROOT,'pooled/primary.csv'));de<-discovery[discovery$FDR<=.05,];po<-read.csv(file.path(ROOT,'pathways/programs_meta_fry.csv'))
summaries<-list();genes<-list();progs<-list()
cohorts<-c('Xu2026_external','Xu2026_external_age','GSE307510_external','GSE307510_external_without_synechia')
for(s in cohorts){
 o<-readRDS(file.path(ROOT,'expression',paste0(s,'.rds')))
 sc<-rank_score(o$expr,panel$GeneSymbol[panel$direction=='Up'],panel$GeneSymbol[panel$direction=='Down'])
 perf<-performance(sc$score,tolower(o$metadata$condition)=='rif',s,'frozen_corrected_discovery_50_gene_panel');perf$up_measured<-length(sc$up);perf$down_measured<-length(sc$down);summaries[[s]]<-perf
 write.csv(data.frame(sample=o$metadata$sample,condition=o$metadata$condition,score=sc$score),file.path(BASE,'results',paste0(s,'_frozen_panel_scores.csv')),row.names=FALSE)
 d<-read.csv(file.path(ROOT,'sensitivity',paste0(s,'_DE.csv')));m<-merge(de,d,by='GeneSymbol',suffixes=c('_discovery','_external'));m$p_directional<-pnorm(-sign(m$z_discovery)*m$z_external);m$FDR_directional<-p.adjust(m$p_directional,'BH');m$same_direction<-sign(m$z_discovery)==sign(m$z_external)
 write.csv(m,file.path(BASE,'results',paste0(s,'_gene_replication.csv')),row.names=FALSE)
 genes[[s]]<-data.frame(study=s,testable=nrow(m),same_direction=sum(m$same_direction),nominal=sum(m$p_directional<=.05),FDR=sum(m$FDR_directional<=.05))
 y<-calcNormFactors(DGEList(o$counts[rownames(o$expr),]));v<-voom(y,o$design,plot=FALSE);ix<-ids2indices(programs,rownames(o$expr));ix<-ix[lengths(ix)>=5]
 f<-fry(v,ix,design=o$design,contrast=2,sort=FALSE);f$program<-rownames(f);f$z_external<-qnorm(f$PValue/2,lower.tail=FALSE)*ifelse(f$Direction=='Up',1,-1)
 q<-merge(po,f,by='program',suffixes=c('_discovery','_external'));q$p_directional<-pnorm(-sign(q$z)*q$z_external);q$FDR_directional<-p.adjust(q$p_directional,'BH')
 write.csv(q,file.path(BASE,'results',paste0(s,'_program_replication.csv')),row.names=FALSE)
 progs[[s]]<-data.frame(study=s,programs=nrow(q),FDR_two_sided=sum(q$FDR_external<=.05),FDR_directional=sum(q$FDR_directional<=.05))
 cat('NEW_EXTERNAL_COMPLETE',s,'AUC',perf$AUC,'genes',sum(m$FDR_directional<=.05),'programs',sum(q$FDR_directional<=.05),'\n');flush.console()
}
write.csv(do.call(rbind,summaries),file.path(BASE,'results/new_cohorts_frozen_signature.csv'),row.names=FALSE)
write.csv(do.call(rbind,genes),file.path(BASE,'results/new_cohorts_gene_replication_summary.csv'),row.names=FALSE)
write.csv(do.call(rbind,progs),file.path(BASE,'results/new_cohorts_program_replication_summary.csv'),row.names=FALSE)
