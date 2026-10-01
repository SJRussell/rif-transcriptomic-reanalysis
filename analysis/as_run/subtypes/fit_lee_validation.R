.libPaths(c('<r_library>',.libPaths()))
suppressPackageStartupMessages({library(edgeR);library(limma)})
BASE<-'<workspace>/subtypes_2026-09-08';ROOT<-file.path(BASE,'corrected_baseline')
d<-read.csv(gzfile(file.path(BASE,'reference/GSE307510/native_raw_counts.csv.gz')),check.names=FALSE)
hg<-read.delim(file.path(ROOT,'inputs/hgnc_complete_set.txt'),quote='',colClasses='character',na.strings='');hg<-hg[hg$status=='Approved',]
lookup<-setNames(hg$symbol,hg$symbol)
# Approved symbols take precedence; resolve only unique previous/alias symbols.
aa<-do.call(rbind,lapply(c('prev_symbol','alias_symbol'),function(col){i<-which(!is.na(hg[[col]]));ss<-strsplit(hg[[col]][i],'|',fixed=TRUE);data.frame(alias=unlist(ss),symbol=rep(hg$symbol[i],lengths(ss)))}))
aa<-unique(aa);aa<-aa[!duplicated(aa$alias)&!duplicated(aa$alias,fromLast=TRUE)&!aa$alias%in%names(lookup),]
lookup<-c(lookup,setNames(aa$symbol,aa$alias))
symbols<-unname(lookup[d[[1]]]);good<-!is.na(symbols)
write.csv(data.frame(source_symbol=d[[1]],approved_symbol=symbols,retained=good),file.path(BASE,'audit/GSE307510_gene_mapping.csv'),row.names=FALSE)
m<-rowsum(as.matrix(d[good,-1,drop=FALSE]),symbols[good],reorder=TRUE)
md<-read.csv(file.path(BASE,'reference/GSE307510/native_biopsy_sample_manifest.csv'));md<-md[match(colnames(m),md$native_processed_raw_column),]
stopifnot(identical(colnames(m),md$native_processed_raw_column))
sm<-data.frame(study='GSE307510',sample=colnames(m),condition=ifelse(md$geo_group=='RIF','RIF','Control'),arm='rnaseq',tech='bulk',diagnosis=md$additional_diagnosis_in_supplement,timing='Post-ovulation day 5',role='small_broader_external_exploratory')
for(e in parse(file.path(ROOT,'scripts/fit_studies.R')))if(is.call(e)&&identical(e[[1]],as.name('<-'))&&identical(e[[2]],as.name('fit_one')))eval(e)
fit_one('GSE307510',m,sm,TRUE,label='GSE307510_external')
fit_one('GSE307510',m,sm[sm$diagnosis!='Synechia',],TRUE,label='GSE307510_external_without_synechia')
write.csv(sm,file.path(BASE,'audit/GSE307510_final_metadata.csv'),row.names=FALSE)
cat('LEE_VALIDATION_FIT_COMPLETE\n')
