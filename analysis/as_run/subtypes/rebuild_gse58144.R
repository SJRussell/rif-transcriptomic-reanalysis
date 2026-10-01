.libPaths(c('<r_library>',.libPaths()))
suppressPackageStartupMessages(library(limma))
options(warn=1)
ROOT<-'<workspace>/subtypes_2026-09-08'
OLD<-'<workspace>/reanalysis_2026-09-08'
REF<-file.path(ROOT,'reference/GSE58144');DEST<-file.path(ROOT,'corrected_baseline')
for(p in c('expression','per_study','sensitivity','audit','inputs','figures','logs','pooled','pathways','validation'))dir.create(file.path(DEST,p),recursive=TRUE,showWarnings=FALSE)
tg<-read.csv(file.path(REF,'array_targets.csv'),check.names=FALSE)
tg$Cy3<-trimws(tg$Cy3);tg$Cy5<-trimws(tg$Cy5)
stopifnot(all(tg$Cy3%in%c('sample','reference')),all(tg$Cy5%in%c('sample','reference')),all(tg$Cy3!=tg$Cy5))
md<-read.csv(file.path(REF,'donor_metadata.csv'),na.strings=c('NA',''),check.names=FALSE)
stopifnot(nrow(tg)==230,nrow(md)==115,all(table(tg$sample)==2),all(table(tg$sample,tg$Cy3)==1))
cached<-file.path(REF,'raw_normalized_donor_probes_v2_trimmed_targets.rds')
if(!file.exists(cached)) {
 R<-G<-NULL;probes<-NULL
 for(i in seq_len(nrow(tg))) {
  d<-read.delim(gzfile(tg$path[i]),check.names=FALSE,quote='',colClasses=c('character',rep('NULL',5),'numeric',rep('NULL',5),'numeric',rep('NULL',7)))
  if(i==1){probes<-d[[1]];R<-G<-matrix(NA_real_,nrow(d),nrow(tg),dimnames=list(NULL,basename(tg$path)))}
  stopifnot(identical(probes,d[[1]]),all(is.finite(d[[2]])),all(is.finite(d[[3]])))
  # Reproduce the imported foreground processing; Ch1 is treated as Cy5,
  # Ch2 as Cy3, anchored to GEO's declared sample/reference values below.
  R[,i]<-if(tg$Cy5[i]=='sample')d[[2]] else d[[3]]
  G[,i]<-if(tg$Cy5[i]=='sample')d[[3]] else d[[2]]
  stopifnot(any(R[,i]!=G[,i]))
  if(i%%20==0)cat('READ',i,'of',nrow(tg),'\n')
 }
 stopifnot(all(R>=0),all(G>=0))
 rg<-new('RGList',list(R=R,G=G,genes=data.frame(ProbeName=probes)))
 bg<-backgroundCorrect(rg,method='normexp',offset=0)
 norm<-normalizeWithinArrays(bg,method='loess')
 stopifnot(all(apply(norm$M,2,sd)>0))
 donor<-avearrays(norm$M,ID=tg$sample)
 pair<-lapply(unique(tg$sample),function(s){j<-which(tg$sample==s);data.frame(sample=s,dye_pair_correlation=cor(norm$M[,j[1]],norm$M[,j[2]]),median_absolute_pair_difference=median(abs(norm$M[,j[1]]-norm$M[,j[2]])))})
 write.csv(do.call(rbind,pair),file.path(DEST,'audit/GSE58144_dye_pair_qc.csv'),row.names=FALSE)
 saveRDS(list(probes=probes,donor=donor,normalization='foreground normexp offset=0; within-array loess; mean dye-swap log2 ratios; Ch1=Cy5 convention anchored to deposited sample/reference values'),cached)
}
o<-readRDS(cached)
gpl<-read.delim(file.path(REF,'GPL15789.tsv'),quote='',check.names=FALSE,na.strings='',fill=TRUE)
lk<-read.csv(file.path(OLD,'inputs/hgnc_lookup.csv'),na.strings=c('NA',''))
gpl_approved<-lk$GeneSymbol[match(gpl$GENE_SYMBOL,lk$input_id)]
ann<-split(gpl_approved,toupper(gpl$SPOT_ID))
sy<-vapply(ann,function(x){y<-unique(na.omit(x));if(length(y)==1)y else NA_character_},character(1))
g<-unname(sy[toupper(o$probes)]);keep<-!is.na(g)&g!=''&apply(o$donor,1,function(x)all(is.finite(x)))
m<-avereps(o$donor[keep,,drop=FALSE],ID=g[keep]);m<-m[,md$sample,drop=FALSE]
write.csv(data.frame(ProbeName=o$probes,approved_symbol=g,retained=keep),file.path(DEST,'audit/GSE58144_raw_probe_mapping.csv'),row.names=FALSE)
geo<-readRDS(file.path(REF,'GSE58144_GEO_approved_expression.rds'));genes<-intersect(rownames(m),rownames(geo));samples<-intersect(colnames(m),colnames(geo))
rows_cor<-function(a,b){a<-a-rowMeans(a);b<-b-rowMeans(b);rowSums(a*b)/sqrt(rowSums(a*a)*rowSums(b*b))}
geo_shared<-geo[genes,samples,drop=FALSE]
raw_r<-rows_cor(m[genes,samples,drop=FALSE],geo_shared)
md$batch<-factor(md$batch);md$condition<-factor(md$condition,levels=c('Control','RIF'))
covs<-model.matrix(~batch+biopsy_day,data=md)[,-1,drop=FALSE]
# Case-blind batch adjustment is only a classifier input; DE uses the raw normalized matrix and its full design.
caseblind<-removeBatchEffect(m,batch=md$batch)
adj_r<-rows_cor(caseblind[genes,samples,drop=FALSE],geo_shared)
write.csv(data.frame(GeneSymbol=genes,raw_to_GEO=raw_r,caseblind_batch_to_GEO=adj_r),file.path(DEST,'audit/GSE58144_raw_GEO_correlations.csv'),row.names=FALSE)
cat('ORIENTATION_CHECK',length(genes),'RAW_MEDIAN_R',median(raw_r),'BATCH_MEDIAN_R',median(adj_r),'POSITIVE_FRACTION',mean(adj_r>0),'\n')
stopifnot(median(raw_r)>.8,mean(raw_r>0)>.95,mean(adj_r>0)>.95)
saveRDS(list(expr=m,classifier_expr=caseblind,metadata=md,normalization=o$normalization),file.path(REF,'rebuilt_all115.rds'))
# Reuse the verified common study-fit writer with a separate output root.
parsed<-parse(file.path(OLD,'scripts/fit_studies.R'))
for(e in parsed)if(is.call(e)&&identical(e[[1]],as.name('<-'))&&identical(e[[2]],as.name('fit_one')))eval(e)
stopifnot(exists('fit_one'));ROOT<-DEST
master<-read.csv(file.path(OLD,'inputs/master_metadata.csv'),check.names=FALSE)
sm<-master[master$study=='GSE58144',];j<-match(sm$sample,md$sample);stopifnot(!anyNA(j),nrow(sm)==113)
stopifnot(all(tolower(sm$condition)==tolower(as.character(md$condition[j]))),sum(md$condition=='RIF')==43,sum(md$condition=='Control')==72)
sm$batch<-as.character(md$batch[j]);sm$biopsy_day<-md$biopsy_day[j]
fit_one('GSE58144',m,sm,FALSE,covariates=covs[j,,drop=FALSE])
fit_one('GSE58144',m,sm,FALSE,label='GSE58144_without_timing',covariates=covs[j,grepl('^batch',colnames(covs)),drop=FALSE])
allsm<-sm[rep(1,nrow(md)),];rownames(allsm)<-NULL;allsm$sample<-md$sample;allsm$condition<-as.character(md$condition);allsm$batch<-as.character(md$batch);allsm$biopsy_day<-md$biopsy_day
fit_one('GSE58144',m,allsm,FALSE,label='GSE58144_all115',covariates=covs)
for(label in c('GSE58144','GSE58144_without_timing','GSE58144_all115')) {
 p<-file.path(DEST,'expression',paste0(label,'.rds'));x<-readRDS(p);x$classifier_expr<-caseblind[rownames(x$expr),x$metadata$sample,drop=FALSE];x$normalization<-o$normalization;saveRDS(x,p)
}
capture.output(sessionInfo(),file=file.path(DEST,'logs/GSE58144_rebuild_sessionInfo.txt'))
cat('GSE58144_REBUILD_COMPLETE\n')
