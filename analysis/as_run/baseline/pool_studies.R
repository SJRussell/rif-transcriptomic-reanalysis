.libPaths(c('<r_library>',.libPaths()))
ROOT<-'<workspace>/subtypes_2026-09-08/corrected_baseline'
md<-read.csv(file.path(ROOT,'inputs/master_metadata.csv'))
files<-list.files(file.path(ROOT,'per_study'),pattern='_DE.csv$',full.names=TRUE)
ps<-lapply(files,read.csv);names(ps)<-vapply(ps,function(x)x$study[1],character(1))
stopifnot(setequal(names(ps),unique(md$study)))
st<-do.call(rbind,lapply(names(ps),function(s){m<-md[md$study==s,];nc<-sum(tolower(m$condition)=='control');nr<-sum(tolower(m$condition)=='rif');data.frame(study=s,arm=m$arm[1],tech=m$tech[1],n_control=nc,n_rif=nr,N=nc+nr,N_effective=4*nc*nr/(nc+nr),genes_tested=nrow(ps[[s]]))}))
rownames(st)<-st$study;write.csv(st,file.path(ROOT,'audit/study_summary.csv'),row.names=FALSE)
pool<-function(keys=names(ps),weights='effective',frames=ps) {
 genes<-sort(unique(unlist(lapply(frames[keys],function(x)x$GeneSymbol))))
 nr<-length(genes); k<-length(keys)
 z<-fc<-se<-matrix(NA_real_,nr,k,dimnames=list(genes,keys))
 for(s in keys) {d<-frames[[s]];i<-match(d$GeneSymbol,genes);z[i,s]<-d$z;fc[i,s]<-d$logFC;se[i,s]<-d$SE}
 n<-switch(weights,effective=st[keys,'N_effective'],total=st[keys,'N'],equal=rep(1,k));w<-sqrt(n)
 present<-!is.na(z)
 combine<-function(cols){if(!length(cols))return(rep(NA_real_,nr));rowSums(sweep(z[,cols,drop=FALSE],2,w[cols],'*'),na.rm=TRUE)/sqrt(as.vector(present[,cols,drop=FALSE]%*%n[cols]))}
 za<-combine(which(st[keys,'arm']=='microarray'));zr<-combine(which(st[keys,'arm']=='rnaseq'));zz<-combine(seq_len(k))
 nma<-rowSums(present[,st[keys,'arm']=='microarray',drop=FALSE]);nrna<-rowSums(present[,st[keys,'arm']=='rnaseq',drop=FALSE])
 lfc<-rowSums(sweep(fc,2,n,'*'),na.rm=TRUE)/as.vector(present%*%n)
 out<-data.frame(GeneSymbol=genes,z=zz,p=2*pnorm(-abs(zz)),logFC_descriptive=lfc,n_microarray=nma,n_rnaseq=nrna,n_studies=rowSums(present),z_microarray=za,p_microarray=2*pnorm(-abs(za)),z_rnaseq=zr,p_rnaseq=2*pnorm(-abs(zr)))
 out$direction_agreement_n<-rowSums(sign(z)==sign(zz),na.rm=TRUE)
 # DerSimonian-Laird heterogeneity and random-effects logFC are descriptive sensitivities.
 wi<-1/se^2;sumw<-rowSums(wi,na.rm=TRUE);mu<-rowSums(wi*fc,na.rm=TRUE)/sumw
 Q<-rowSums(wi*(fc-mu)^2,na.rm=TRUE);df<-out$n_studies-1
 C<-sumw-rowSums(wi^2,na.rm=TRUE)/sumw;tau2<-pmax(0,(Q-df)/C)
 rw<-1/(se^2+tau2);out$RE_logFC<-rowSums(rw*fc,na.rm=TRUE)/rowSums(rw,na.rm=TRUE)
 out$RE_SE<-sqrt(1/rowSums(rw,na.rm=TRUE));out$RE_p<-2*pnorm(-abs(out$RE_logFC/out$RE_SE))
 out$Q<-Q;out$I2<-ifelse(Q>0,pmax(0,(Q-df)/Q)*100,0);out$tau2<-tau2
 out
}
adj<-function(d){d$FDR<-p.adjust(d$p,'BH');d$FDR_microarray<-p.adjust(d$p_microarray,'BH');d$FDR_rnaseq<-p.adjust(d$p_rnaseq,'BH');d$RE_FDR<-p.adjust(d$RE_p,'BH');d$Direction<-ifelse(d$z>0,'Up','Down');d[order(d$p),]}
all<-pool();write.csv(all,file.path(ROOT,'pooled/all_tested_genes.csv'),row.names=FALSE)
primary<-adj(all[all$n_microarray>=2&all$n_rnaseq>=2,]);de<-primary[primary$FDR<=.05,]
write.csv(primary,file.path(ROOT,'pooled/primary.csv'),row.names=FALSE);write.csv(de,file.path(ROOT,'pooled/DEGs.csv'),row.names=FALSE)
write.csv(primary[primary$FDR_microarray<=.05&primary$FDR_rnaseq<=.05&sign(primary$z_microarray)==sign(primary$z_rnaseq),],file.path(ROOT,'pooled/replicated_both_platforms.csv'),row.names=FALSE)
summaries<-list();add<-function(d,name){d<-adj(d);write.csv(d,file.path(ROOT,'sensitivity',paste0('pooled_',name,'.csv')),row.names=FALSE);summaries[[name]]<<-data.frame(analysis=name,genes=nrow(d),DE=sum(d$FDR<=.05),up=sum(d$FDR<=.05&d$z>0),down=sum(d$FDR<=.05&d$z<0));invisible(d)}
add(primary,'primary');add(all[all$n_microarray==sum(st$arm=='microarray')&all$n_rnaseq==sum(st$arm=='rnaseq'),],'strict_all_studies')
for(weight in c('total','equal')){d<-pool(weights=weight);add(d[d$GeneSymbol %in% primary$GeneSymbol,],paste0(weight,'_weights'))}
loo<-list()
for(s in c(names(ps),'all_pseudobulk','all_microarray','all_bulk_rnaseq')) {
 omit<-switch(s,all_pseudobulk=st$study[st$tech=='pseudobulk'],all_microarray=st$study[st$arm=='microarray'],all_bulk_rnaseq=st$study[st$tech=='bulk'],s)
 d<-pool(setdiff(names(ps),omit));d<-d[match(primary$GeneSymbol,d$GeneSymbol),];unest<-!is.finite(d$p);d$p[unest]<-1;d$GeneSymbol<-primary$GeneSymbol;d<-adj(d)
 write.csv(d,file.path(ROOT,'sensitivity',paste0('leaveout_',s,'.csv')),row.names=FALSE)
 m<-match(de$GeneSymbol,d$GeneSymbol)
 loo[[s]]<-data.frame(omitted=s,family=nrow(primary),estimable=sum(!unest),DE=sum(d$FDR<=.05),baseline_retained=sum(d$FDR[m]<=.05),baseline_same_direction=sum(sign(d$z[m])==sign(de$z),na.rm=TRUE))
}
write.csv(do.call(rbind,loo),file.path(ROOT,'sensitivity/leave_one_out_summary.csv'),row.names=FALSE)
for(label in c('GSE250130_mixed_time','PRJNA576073_reverse','GSE58144_without_timing','GSE58144_all115','Bui2024_unadjusted_all','Bui2024_unadjusted_complete','GSE250130_composition_adjusted')) {
 file<-file.path(ROOT,'sensitivity',paste0(label,'_DE.csv'));if(!file.exists(file))next
 fp<-ps;d<-read.csv(file);fp[[d$study[1]]]<-d
 # Count changes for the deliberately unmatched GSE250130 contrast must change its weight as well.
 oldst<-st
 if(label=='GSE250130_mixed_time'){st['GSE250130','N']<-22;st['GSE250130','N_effective']<-4*12*10/22}
 if(label=='Bui2024_unadjusted_all'){st['Bui2024','N']<-46;st['Bui2024','N_effective']<-46}
 if(label=='GSE58144_all115'){st['GSE58144','N']<-115;st['GSE58144','N_effective']<-4*43*72/115}
 x<-pool(frames=fp);st<-oldst
 add(x[x$GeneSymbol %in% primary$GeneSymbol,],label)
}
write.csv(do.call(rbind,summaries),file.path(ROOT,'sensitivity/analysis_summary.csv'),row.names=FALSE)
print(summaries);print(head(de[,c('GeneSymbol','z','FDR','logFC_descriptive','I2')],15));print(do.call(rbind,loo))
saveRDS(list(ps=ps,study_summary=st,primary=primary),file.path(ROOT,'pooled/pooling_inputs.rds'))
cat('POOLING_COMPLETE\n')
