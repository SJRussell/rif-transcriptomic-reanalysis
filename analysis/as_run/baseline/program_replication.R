.libPaths(c('<r_library>',.libPaths()))
suppressPackageStartupMessages({library(limma);library(edgeR)})
ROOT<-'<workspace>/subtypes_2026-09-08/corrected_baseline'
programs<-readRDS(file.path(ROOT,'pathways/fixed_programs.rds'));st<-read.csv(file.path(ROOT,'audit/study_summary.csv'))
results<-list()
cohorts<-c(st$study,'GSE248784_external','GSE248784_external_mt10','GSE248784_external_shared_reference','GSE243550_external')
for(s in cohorts) {
 f<-file.path(ROOT,'expression',paste0(s,'.rds'));if(!file.exists(f))next
 o<-readRDS(f);inds<-ids2indices(programs,rownames(o$expr));inds<-inds[lengths(inds)>=5]
 if(!is.null(o$counts)){y<-calcNormFactors(DGEList(o$counts[rownames(o$expr),]));em<-voom(y,o$design,plot=FALSE)}else em<-o$expr
 fryres<-fry(em,index=inds,design=o$design,contrast=2,sort=FALSE)
 camerares<-camera(em,index=inds,design=o$design,contrast=2,inter.gene.cor=NA,sort=FALSE)
 d<-data.frame(study=s,program=rownames(fryres),genes=fryres$NGenes,direction=fryres$Direction,p=fryres$PValue,FDR=fryres$FDR,camera_direction=camerares$Direction[match(rownames(fryres),rownames(camerares))],camera_p=camerares$PValue[match(rownames(fryres),rownames(camerares))])
 d$camera_FDR<-p.adjust(d$camera_p,'BH');d$z<-qnorm(d$p/2,lower.tail=FALSE)*ifelse(d$direction=='Up',1,-1)
 results[[s]]<-d;cat(s,'FRY FDR',sum(d$FDR<=.05),'CAMERA FDR',sum(d$camera_FDR<=.05),'\n')
}
ca<-do.call(rbind,results);write.csv(ca,file.path(ROOT,'pathways/programs_fry_and_external.csv'),row.names=FALSE)
ww<-setNames(sqrt(st$N_effective),st$study)
po<-do.call(rbind,lapply(split(ca[ca$study %in% st$study,],ca$program[ca$study %in% st$study]),function(d){z<-sum(ww[d$study]*d$z)/sqrt(sum(ww[d$study]^2));i<-d$study!='GSE111974';zloo<-sum(ww[d$study[i]]*d$z[i])/sqrt(sum(ww[d$study[i]]^2));data.frame(program=d$program[1],studies=nrow(d),z=z,p=2*pnorm(-abs(z)),same_direction=sum(sign(d$z)==sign(z)),z_without_GSE111974=zloo,p_without_GSE111974=2*pnorm(-abs(zloo)))}))
po$FDR<-p.adjust(po$p,'BH');po$FDR_without_GSE111974<-p.adjust(po$p_without_GSE111974,'BH');po<-po[order(po$p),]
write.csv(po,file.path(ROOT,'pathways/programs_meta_fry.csv'),row.names=FALSE)
for(s in c('GSE248784_external','GSE248784_external_mt10','GSE248784_external_shared_reference','GSE243550_external')) {
 ext<-ca[ca$study==s,];if(!nrow(ext))next
 d<-merge(po,ext,by='program',suffixes=c('_discovery','_external'))
 d$same_direction_external<-sign(d$z_discovery)==sign(d$z_external)
 d$p_directional_external<-pnorm(-sign(d$z_discovery)*d$z_external)
 d$FDR_directional_external<-p.adjust(d$p_directional_external,'BH')
 write.csv(d,file.path(ROOT,'pathways',paste0('program_replication_',s,'.csv')),row.names=FALSE)
}
print(head(po,12));cat('PROGRAM_REPLICATION_COMPLETE\n')
