Sys.setenv(LOAD_FUNCTIONS_ONLY='1');source('<workspace>/subtypes_2026-09-08/corrected_baseline/scripts/fit_studies.R')
o<-readRDS(file.path(ROOT,'inputs/GSE71331_normalized_probes.rds'))
keep<-!is.na(o$probe_symbol)
if('ControlType' %in% colnames(o$genes))keep<-keep&o$genes$ControlType==0
m<-suppressWarnings(log2(o$expression[keep,,drop=FALSE]+1));rownames(m)<-o$probe_symbol[keep]
rawfiles<-sort(list.files(file.path(ROOT,'inputs/GSE71331_raw'),pattern='^GSM.*txt$',full.names=TRUE));colnames(m)<-sub('.*(GSM[0-9]+).*','\\1',basename(rawfiles))
m<-harmonize(m,'GSE71331_extra_log',FALSE);sm<-md[md$study=='GSE71331',]
extra<-fit_one('GSE71331',m,sm,FALSE,label='GSE71331_extra_log')
base<-read.csv(file.path(ROOT,'per_study/GSE71331_DE.csv'));d<-merge(base,extra,by='GeneSymbol',suffixes=c('_correct','_extra_log'))
write.csv(data.frame(genes=nrow(d),spearman_Z=cor(d$z_correct,d$z_extra_log,method='spearman'),same_direction=sum(sign(d$z_correct)==sign(d$z_extra_log)),median_abs_logFC_correct=median(abs(d$logFC_correct)),median_abs_logFC_extra_log=median(abs(d$logFC_extra_log))),file.path(ROOT,'audit/GSE71331_transform_diagnostic.csv'),row.names=FALSE)
cat('TRANSFORM_DIAGNOSTIC_COMPLETE\n')
