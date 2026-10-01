# Adapted 2026-09-18 from as_run/subtypes/evaluate_new_cohorts.R (frozen 50-gene panel, gene and program replication) plus the comparator-panel block of as_run/baseline/validate_signatures.R (lines 54-63, Altmae/Chettiar).
# Changes: paths; one cohort (GSE207362_external); panel read from inputs/signature_training_membership.csv (held_out == GSE248784_external = panel trained on all discovery cohorts); discovery DEGs from inputs/discovery_DEGs.csv; the frozen discovery program statistics (programs_meta_fry.csv, not shipped) are taken from the identical discovery columns of the frozen S10 program_replication table with FDR_discovery renamed back to FDR; comparator TSVs read from reference_inputs. Seed 20260908, thresholds, bootstrap count and test families unchanged.
suppressPackageStartupMessages({library(limma);library(edgeR)})
WD<-'<cluster_scratch>/gse207362_external_20260918';ROOT<-file.path(WD,'work');REF<-file.path(WD,'inputs/frozen/Analysis_code/reference_inputs');OUT<-file.path(WD,'outputs');set.seed(20260908)
for(e in parse(file.path(WD,'scripts/validate_signatures.R')))if(is.call(e)&&identical(e[[1]],as.name('<-'))&&as.character(e[[2]])%in%c('auc','performance','rank_score'))eval(e)
panel<-read.csv(file.path(WD,'inputs/signature_training_membership.csv'));panel<-panel[panel$held_out=='GSE248784_external',];stopifnot(nrow(panel)==50)
programs<-readRDS(file.path(REF,'fixed_programs.rds'));discovery<-read.csv(file.path(WD,'inputs/discovery_DEGs.csv'));de<-discovery[discovery$FDR<=.05,];stopifnot(nrow(de)==1840,length(programs)==57)
tmpl<-read.csv(file.path(WD,'inputs/frozen/Supplementary_tables/S10/GSE307510_external_program_replication.csv'))
po<-tmpl[,c('program','studies','z','p','same_direction','z_without_GSE111974','p_without_GSE111974','FDR_discovery','FDR_without_GSE111974')];names(po)[names(po)=='FDR_discovery']<-'FDR';stopifnot(nrow(po)==57,setequal(po$program,names(programs)))
summaries<-list();genes<-list();progs<-list()
cohorts<-c('GSE207362_external')
for(s in cohorts){
 o<-readRDS(file.path(ROOT,'expression',paste0(s,'.rds')))
 sc<-rank_score(o$expr,panel$GeneSymbol[panel$direction=='Up'],panel$GeneSymbol[panel$direction=='Down'])
 perf<-performance(sc$score,tolower(o$metadata$condition)=='rif',s,'frozen_corrected_discovery_50_gene_panel');perf$up_measured<-length(sc$up);perf$down_measured<-length(sc$down);summaries[[s]]<-perf
 write.csv(data.frame(sample=o$metadata$sample,condition=o$metadata$condition,score=sc$score),file.path(OUT,paste0(s,'_frozen_panel_scores.csv')),row.names=FALSE)
 d<-read.csv(file.path(ROOT,'sensitivity',paste0(s,'_DE.csv')));m<-merge(de,d,by='GeneSymbol',suffixes=c('_discovery','_external'));m$p_directional<-pnorm(-sign(m$z_discovery)*m$z_external);m$FDR_directional<-p.adjust(m$p_directional,'BH');m$same_direction<-sign(m$z_discovery)==sign(m$z_external)
 write.csv(m,file.path(OUT,paste0(s,'_gene_replication.csv')),row.names=FALSE)
 genes[[s]]<-data.frame(study=s,testable=nrow(m),same_direction=sum(m$same_direction),nominal=sum(m$p_directional<=.05),FDR=sum(m$FDR_directional<=.05))
 y<-calcNormFactors(DGEList(o$counts[rownames(o$expr),]));v<-voom(y,o$design,plot=FALSE);ix<-ids2indices(programs,rownames(o$expr));ix<-ix[lengths(ix)>=5]
 f<-fry(v,ix,design=o$design,contrast=2,sort=FALSE);f$program<-rownames(f);f$z_external<-qnorm(f$PValue/2,lower.tail=FALSE)*ifelse(f$Direction=='Up',1,-1)
 q<-merge(po,f,by='program',suffixes=c('_discovery','_external'));q$p_directional<-pnorm(-sign(q$z)*q$z_external);q$FDR_directional<-p.adjust(q$p_directional,'BH')
 write.csv(q,file.path(OUT,paste0(s,'_program_replication.csv')),row.names=FALSE)
 progs[[s]]<-data.frame(study=s,programs=nrow(q),FDR_two_sided=sum(q$FDR_external<=.05),FDR_directional=sum(q$FDR_directional<=.05))
 cat('NEW_EXTERNAL_COMPLETE',s,'AUC',perf$AUC,'genes',sum(m$FDR_directional<=.05),'programs',sum(q$FDR_directional<=.05),'\n');flush.console()
 # Two previously published panels, fixed before external testing (validate_signatures.R lines 54-63; altmae2017 direction reversed once, explicitly).
 lk<-read.csv(file.path(REF,'hgnc_lookup.csv'));benchmark<-list()
 for(pn in c('altmae2017','chettiar2024')) {
  dd<-read.delim(file.path(REF,paste0('comparator_',pn,'_meta_signature.tsv')));dd$GeneSymbol<-lk$GeneSymbol[match(dd$GeneSymbol,lk$input_id)];dd<-dd[!is.na(dd$GeneSymbol),]
  up<-unique(dd$GeneSymbol[tolower(dd$Direction)=='up']);down<-unique(dd$GeneSymbol[tolower(dd$Direction)=='down'])
  if(pn=='altmae2017'){tmp<-up;up<-down;down<-tmp}
  scp<-rank_score(o$expr,up,down);rr<-performance(scp$score,tolower(o$metadata$condition)=='rif',s,pn);rr$genes_up_measured<-length(scp$up);rr$genes_down_measured<-length(scp$down);benchmark[[pn]]<-rr
  write.csv(data.frame(sample=o$metadata$sample,condition=o$metadata$condition,score=scp$score),file.path(OUT,paste0(s,'_',pn,'_panel_scores.csv')),row.names=FALSE)
  cat('PANEL',pn,'AUC',rr$AUC,'\n')
 }
 write.csv(do.call(rbind,benchmark),file.path(OUT,paste0(s,'_published_panel_benchmarks.csv')),row.names=FALSE)
}
write.csv(do.call(rbind,summaries),file.path(OUT,'GSE207362_external_frozen_signature.csv'),row.names=FALSE)
write.csv(do.call(rbind,genes),file.path(OUT,'GSE207362_external_gene_replication_summary.csv'),row.names=FALSE)
write.csv(do.call(rbind,progs),file.path(OUT,'GSE207362_external_program_replication_summary.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(ROOT,'logs/evaluate_gse207362_sessionInfo.txt'))
cat('GSE207362_EVALUATION_COMPLETE\n')
