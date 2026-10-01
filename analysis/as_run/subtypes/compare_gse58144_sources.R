.libPaths(c('<r_library>',.libPaths()))
suppressPackageStartupMessages(library(limma))
ROOT<-'<workspace>/subtypes_2026-09-08';REF<-file.path(ROOT,'reference/GSE58144');OLD<-'<workspace>/reanalysis_2026-09-08'
con<-gzfile(file.path(REF,'GSE58144_series_matrix.txt.gz'),'rt');lines<-readLines(con);close(con);begin<-which(lines=='!series_matrix_table_begin');end<-which(lines=='!series_matrix_table_end');stopifnot(length(begin)==1,length(end)==1)
d<-read.delim(text=paste(lines[(begin+1):(end-1)],collapse='\n'),check.names=FALSE);m<-as.matrix(d[,-1]);rownames(m)<-as.character(d[[1]])
gpl<-read.delim(file.path(REF,'GPL15789.tsv'),quote='',check.names=FALSE,na.strings='',fill=TRUE)
lookup<-read.csv(file.path(OLD,'inputs/hgnc_lookup.csv'));sy<-gpl$GENE_SYMBOL[match(rownames(m),gpl$ID)];approved<-lookup$GeneSymbol[match(sy,lookup$input_id)];good<-!is.na(approved)&apply(m,1,function(x)all(is.finite(x)));m<-avereps(m[good,,drop=FALSE],ID=approved[good])
saveRDS(m,file.path(REF,'GSE58144_GEO_approved_expression.rds'))
write.csv(data.frame(probe=d[[1]],input_symbol=sy,approved_symbol=approved),file.path(REF,'GEO_probe_mapping.csv'),row.names=FALSE)
old<-readRDS(file.path(OLD,'expression/GSE58144.rds'));shared<-intersect(colnames(m),colnames(old$expr));genes<-intersect(rownames(m),rownames(old$expr));a<-m[genes,shared];b<-old$expr[genes,shared]
r<-vapply(seq_len(nrow(a)),function(i)cor(a[i,],b[i,]),numeric(1));q<-data.frame(GeneSymbol=genes,correlation_GEO_original=r);write.csv(q,file.path(ROOT,'audit/GSE58144_GEO_original_correlations.csv'),row.names=FALSE)
cat('GEO_ROWS',nrow(m),'GEO_SAMPLES',ncol(m),'COMMON_GENES',length(genes),'MEDIAN_GENE_R',median(r,na.rm=TRUE),'NEGATIVE_FRACTION',mean(r<0,na.rm=TRUE),'\n')
e<-new.env();load(file.path(ROOT,'reference/MetaRIF/MetaRIF_classifier_gene_plus.rda'),envir=e)
hg<-read.delim(file.path(OLD,'inputs/hgnc_complete_set.txt'),colClasses='character',quote='');ei<-setNames(hg$entrez_id,hg$symbol)
ids<-ei[rownames(m)];kk<-!is.na(ids)&ids!='';mm<-avereps(m[kk,,drop=FALSE],ID=ids[kk]);mm<-t(scale(t(mm)))
for(mod in c('best_model_RIF11','best_model_RIF22')){
 td<-e[[mod]]$trainingData;s<-intersect(rownames(td),colnames(mm));g<-intersect(sub('gene_','',colnames(td)[-1]),rownames(mm));x<-t(mm[g,s]);y<-as.matrix(td[s,paste0('gene_',g)]);r<-vapply(seq_len(ncol(x)),function(j)cor(x[,j],y[,j]),numeric(1));cat(mod,'N',length(s),'MEDIAN_GEO_MODEL_R',median(r),'\n')
}
cat('GEO samples absent from the original subset:',setdiff(colnames(m),old$metadata$sample),'\n')
