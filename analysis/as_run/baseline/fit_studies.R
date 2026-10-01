.libPaths(c('<r_library>',.libPaths()))
if(Sys.getenv('LOAD_FUNCTIONS_ONLY')!='1')stop('This preserved helper must be loaded with LOAD_FUNCTIONS_ONLY=1. Rebuild GSE58144 with ../../scripts/rebuild_gse58144.R; do not reuse its superseded processed CSV.')
suppressPackageStartupMessages({library(limma);library(edgeR)})
ROOT <- '<workspace>/subtypes_2026-09-08/corrected_baseline'
SRC <- '<workspace>/sdrive_audit_20260908'
args <- commandArgs(trailingOnly=TRUE)
md <- read.csv(file.path(ROOT,'inputs/master_metadata.csv'),check.names=FALSE)
lookup <- read.csv(file.path(ROOT,'inputs/hgnc_lookup.csv'),check.names=FALSE,na.strings='')
namesym <- setNames(lookup$GeneSymbol,lookup$input_id)
namestatus <- setNames(lookup$mapping_status,lookup$input_id)
harmonize <- function(m,study,counts=FALSE) {
 ids <- trimws(rownames(m)); key<-sub('^(ENSG[0-9]+)\\.[0-9]+$','\\1',ids)
 mapped<-unname(namesym[key]); status<-unname(namestatus[key]); status[is.na(status)]<-'unmapped'
 audit<-data.frame(input_id=ids,lookup_id=key,GeneSymbol=mapped,status=status)
 write.csv(audit,file.path(ROOT,'audit',paste0(study,'_identifier_mapping.csv')),row.names=FALSE)
 keep<-!is.na(mapped)&mapped!=''&apply(m,1,function(x)all(is.finite(x)))
 m<-m[keep,,drop=FALSE]; g<-mapped[keep]
 if(counts) m<-rowsum(m,g,reorder=TRUE) else m<-avereps(m,ID=g)
 stopifnot(!anyDuplicated(rownames(m)),all(is.finite(m)))
 attr(m,'source_rows')<-length(ids); m
}
read_matrix <- function(path) {
 d<-read.csv(path,check.names=FALSE)
 m<-as.matrix(d[,-1,drop=FALSE]);mode(m)<-'numeric';rownames(m)<-d[[1]]
 colnames(m)<-sub('[.]CEL$','',gsub('-','_',colnames(m),fixed=TRUE));m
}
fit_one <- function(study,m,sm,counts,label=study,covariates=NULL,quality_weights=FALSE) {
 j<-match(sm$sample,colnames(m));stopifnot(!anyNA(j),!anyDuplicated(sm$sample))
 m<-m[,j,drop=FALSE]
 grp<-factor(tolower(sm$condition),levels=c('control','rif'))
 design<-model.matrix(~grp)
 if(!is.null(covariates)) design<-cbind(design,covariates)
 stopifnot(qr(design)$rank==ncol(design),all(table(grp)>=2),nrow(design)>ncol(design))
 measured<-rownames(m)
 if(counts) {
  stopifnot(all(m>=0),max(abs(m-round(m)))<1e-6)
  y<-DGEList(m); keep<-filterByExpr(y,design=design)
  y<-calcNormFactors(y[keep,,keep.lib.sizes=FALSE],method='TMM')
  pdf(file.path(ROOT,'figures',paste0('QC_voom_',label,'.pdf')),width=7,height=6)
  v<-if(quality_weights)voomWithQualityWeights(y,design,plot=TRUE) else voom(y,design,plot=TRUE);dev.off()
  fit<-eBayes(lmFit(v,design),robust=TRUE)
  expr<-v$E
  write.csv(data.frame(sample=sm$sample,library_size=y$samples$lib.size,norm_factor=y$samples$norm.factors),file.path(ROOT,'audit',paste0(label,'_library_qc.csv')),row.names=FALSE)
  if(quality_weights)write.csv(data.frame(sample=sm$sample,sample_weight=v$targets$sample.weights),file.path(ROOT,'audit',paste0(label,'_sample_weights.csv')),row.names=FALSE)
 } else {
  keep<-apply(m,1,sd)>0
  expr<-m[keep,,drop=FALSE];fit<-eBayes(lmFit(expr,design),robust=TRUE)
 }
 tt<-topTable(fit,coef=2,number=Inf,sort.by='none')
 se<-fit$stdev.unscaled[,2]*sqrt(fit$s2.post)
 df<-fit$df.total
 logp<-log(2)+pt(-abs(fit$t[,2]),df=df,log.p=TRUE)
 z<-qnorm(logp-log(2),lower.tail=FALSE,log.p=TRUE)*sign(tt$logFC)
 res<-data.frame(GeneSymbol=rownames(tt),logFC=tt$logFC,SE=se,t=tt$t,df=df,P.Value=exp(logp),log_p=logp,z=z,adj.P.Val=p.adjust(exp(logp),'BH'),study=study,arm=sm$arm[1],tech=sm$tech[1])
 stopifnot(all(is.finite(res$z)),all(is.finite(res$SE)),all(res$SE>0))
 write.csv(res,file.path(ROOT,if(label==study)'per_study' else 'sensitivity',paste0(label,'_DE.csv')),row.names=FALSE)
 write.csv(data.frame(GeneSymbol=measured,tested=measured %in% rownames(expr)),file.path(ROOT,'audit',paste0(label,'_coverage.csv')),row.names=FALSE)
 saveRDS(list(expr=expr,metadata=sm,design=design,counts=if(counts)m else NULL),file.path(ROOT,'expression',paste0(label,'.rds')))
 print(data.frame(study=label,n_control=sum(grp=='control'),n_rif=sum(grp=='rif'),represented=length(measured),tested=nrow(res),DE=sum(res$adj.P.Val<=.05)))
 invisible(res)
}
studies<-unique(md$study)
if(Sys.getenv('LOAD_FUNCTIONS_ONLY')=='1') studies<-character()
if(length(args)) studies<-intersect(studies,args)
for(study in studies) {
 sm<-md[md$study==study,];counts<-sm$arm[1]=='rnaseq'
 if(study=='PRJNA576073') {
  f<-file.path(ROOT,'inputs/PRJNA576073/gene_counts_s0.txt')
  if(!file.exists(f)) {cat('WAITING_FOR_RECOUNT',study,'\n');next}
  d<-read.delim(f,skip=1,check.names=FALSE)
  m<-as.matrix(d[,-(1:6)]);rownames(m)<-d$Geneid;colnames(m)<-sub('_.*','',basename(colnames(m)))
 } else if(study=='GSE71331') {
  rawfiles<-sort(list.files(file.path(ROOT,'inputs/GSE71331_raw'),pattern='^GSM.*txt$',full.names=TRUE))
  cached<-file.path(ROOT,'inputs/GSE71331_normalized_probes.rds')
  if(file.exists(cached)) {o<-readRDS(cached);raw<-list(genes=o$genes);norm<-list(E=o$expression)} else {
   raw<-read.maimages(rawfiles,source='agilent',green.only=TRUE)
   norm<-normalizeBetweenArrays(backgroundCorrect(raw,method='normexp',offset=0),method='quantile')
  }
  ann<-read.delim(file.path(ROOT,'inputs/GSE71331_annotation.txt'),check.names=FALSE)
  platform<-read.delim(file.path(SRC,'MICROARRAY_DATA/GSE71331/GPL19072_052909_D_GEO_20130704.txt'),check.names=FALSE,quote='',comment.char='')
  # Annotation ProbeID is the platform feature-row ID, not Agilent's internal ProbeUID.
  # Bridge through the actual probe name. Repeated physical spots share one probe sequence/name.
  ann$ProbeName<-platform$NAME[match(as.character(ann$ProbeID),as.character(platform$ID))]
  stopifnot(!anyNA(ann$ProbeName),all(raw$genes$ProbeName %in% platform$NAME))
  candidates<-split(ann$GeneSymbol,as.character(ann$ProbeName))
  probe<-as.character(raw$genes$ProbeName)
  unique_symbols<-vapply(candidates,function(x){x<-unique(x);x<-x[!is.na(x)&x!=''];if(length(x)==1)x else NA_character_},character(1))
  syms<-unname(unique_symbols[probe])
  write.csv(data.frame(ProbeUID=raw$genes$ProbeUID,ProbeName=probe,GeneSymbol=syms),file.path(ROOT,'audit/GSE71331_corrected_probe_map.csv'),row.names=FALSE)
  keep<-!is.na(syms)
  if('ControlType' %in% colnames(raw$genes)) keep<-keep&raw$genes$ControlType==0
  m<-norm$E[keep,,drop=FALSE];rownames(m)<-syms[keep];colnames(m)<-sub('.*(GSM[0-9]+).*','\\1',basename(rawfiles))
  saveRDS(list(genes=raw$genes,expression=norm$E,probe_symbol=syms),file.path(ROOT,'inputs/GSE71331_normalized_probes.rds'))
  cat('GSE71331 normalized probe quantiles:',quantile(norm$E,c(.01,.25,.5,.75,.99)),'\n')
 } else m<-read_matrix(sm$expr_path[1])
 m<-harmonize(m,study,counts)
 covariates<-NULL
 if(study=='Bui2024') {
  covariates<-cbind(biopsy_day=as.numeric(sm$biopsy_day)-7)
  stopifnot(!anyNA(covariates))
 }
 if(study=='GSE58144') {
  ph<-read.csv(file.path(SRC,'MICROARRAY_DATA/meta data/metadata_mar27/phenotype_data_GSE58144.csv'),check.names=FALSE)
  time<-as.numeric(ph$biopsy_timing[match(sm$sample,ph$geo_accession)])
  stopifnot(!anyNA(time));covariates<-cbind(biopsy_day=time-7)
 }
 fit_one(study,m,sm,counts,covariates=covariates)
 if(study=='GSE58144') fit_one(study,m,sm,counts,label=paste0(study,'_without_timing'))
 if(study=='Bui2024') {
  old<-read.csv(file.path(ROOT,'inputs/original_metadata.csv'));old<-old[old$study==study,]
  fit_one(study,m,old,TRUE,label='Bui2024_unadjusted_all')
  fit_one(study,m,sm,TRUE,label='Bui2024_unadjusted_complete')
 }
 if(study=='GSE250130') {
  old<-read.csv(file.path(ROOT,'inputs/original_metadata.csv'));old<-old[old$study==study,]
  fit_one(study,m,old,TRUE,label='GSE250130_mixed_time')
 }
 if(study=='PRJNA576073') {
  d<-read.delim(file.path(ROOT,'inputs/PRJNA576073/gene_counts_s2.txt'),skip=1,check.names=FALSE)
  mm<-as.matrix(d[,-(1:6)]);rownames(mm)<-d$Geneid;colnames(mm)<-sub('_.*','',basename(colnames(mm)))
  mm<-harmonize(mm,paste0(study,'_reverse'),TRUE)
  fit_one(study,mm,sm,TRUE,label=paste0(study,'_reverse'))
 }
}
capture.output(sessionInfo(),file=file.path(ROOT,'logs/fit_studies_sessionInfo.txt'))
cat('STUDY_FITS_COMPLETE\n')
