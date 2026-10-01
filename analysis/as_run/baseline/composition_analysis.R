Sys.setenv(LOAD_FUNCTIONS_ONLY='1')
source('<workspace>/subtypes_2026-09-08/corrected_baseline/scripts/fit_studies.R')
comp<-read.csv(file.path(ROOT,'inputs/Cao2025_donor_composition.csv'),check.names=FALSE)
mapping<-read.csv('<workspace>/sdrive_audit_20260908/metadata_audit/GSE250130_barcode_mapping.csv')
mapping$donor<-sub('^GSM[0-9]+_(E[0-9]+)_.*','\\1',mapping$deposited_folder)
comp$archive_sample<-mapping$archive_sample[match(comp$donor,mapping$donor)]
celltypes<-c('Stromal','Unciliated','Ciliated','Endothelial','NK/T','Myeloid','B','Mast')
rows<-list()
for(label in c('matched_LH7','LH7_including_PGT','mixed_time')) {
 keep<-comp$condition=='RIF'|switch(label,matched_LH7=comp$group=='LH7',LH7_including_PGT=comp$group %in% c('LH7','LH7_PGT'),mixed_time=comp$group %in% c('LH5','LH7','LH9','LH11'))
 d<-comp[keep,];control<-d$condition=='Control';nc<-sum(control);nr<-sum(!control)
 y<-as.matrix(d[,celltypes]);clr<-log(y+1e-4)-rowMeans(log(y+1e-4))
 if(choose(nrow(d),nc)<10000) perm<-combn(nrow(d),nc) else {set.seed(20260908);perm<-replicate(10000,sample(nrow(d),nc))}
 tests<-lapply(seq_along(celltypes),function(j){obs<-mean(clr[!control,j])-mean(clr[control,j]);null<-apply(perm,2,function(ix)mean(clr[-ix,j])-mean(clr[ix,j]));p<-mean(abs(null)>=abs(obs)-1e-12);data.frame(analysis=label,cell_type=celltypes[j],n_control=nc,n_rif=nr,control_median=median(y[control,j]),rif_median=median(y[!control,j]),CLR_difference=obs,p=p,permutations=ncol(perm))})
 tt<-do.call(rbind,tests);tt$FDR<-p.adjust(tt$p,'BH');rows[[label]]<-tt
}
write.csv(do.call(rbind,rows),file.path(ROOT,'sensitivity/Cao2025_cell_composition_tests.csv'),row.names=FALSE)
write.csv(comp,file.path(ROOT,'audit/Cao2025_composition_sample_mapping.csv'),row.names=FALSE)
obj<-readRDS(file.path(ROOT,'expression/GSE250130.rds'))
d<-comp[match(obj$metadata$sample,comp$archive_sample),]
stopifnot(!anyNA(d$donor),all(d$condition==obj$metadata$condition))
y<-as.matrix(d[,celltypes]);clr<-log(y+1e-4)-rowMeans(log(y+1e-4))
pc<-prcomp(clr,center=TRUE,scale.=FALSE)
cov<-pc$x[,1:2,drop=FALSE];colnames(cov)<-c('composition_PC1','composition_PC2')
fit_one('GSE250130',obj$counts,obj$metadata,TRUE,label='GSE250130_composition_adjusted',covariates=cov)
write.csv(data.frame(sample=obj$metadata$sample,cov),file.path(ROOT,'audit/composition_covariates.csv'),row.names=FALSE)
cat('Composition PC variance:',pc$sdev[1:2]^2/sum(pc$sdev^2),'\n')
print(do.call(rbind,rows));cat('COMPOSITION_COMPLETE\n')
