Sys.setenv(LOAD_FUNCTIONS_ONLY='1')
source('<workspace>/subtypes_2026-09-08/corrected_baseline/scripts/fit_studies.R')
b<-readRDS(file.path(ROOT,'pooled/pooling_inputs.rds'));ps<-b$ps;st<-b$study_summary;primary<-b$primary
for(s in if('--summarize' %in% commandArgs(trailingOnly=TRUE))character() else st$study[st$arm=='rnaseq']) {
 o<-readRDS(file.path(ROOT,'expression',paste0(s,'.rds')))
 cv<-if(ncol(o$design)>2)o$design[,-(1:2),drop=FALSE] else NULL
 fit_one(s,o$counts,o$metadata,TRUE,label=paste0(s,'_quality_weights'),covariates=cv,quality_weights=TRUE)
}
keys<-names(ps);g<-primary$GeneSymbol
zmat<-function(frames){z<-matrix(NA_real_,length(g),length(keys),dimnames=list(g,keys));for(s in keys)z[,s]<-frames[[s]]$z[match(g,frames[[s]]$GeneSymbol)];z}
combine<-function(z,n){rowSums(sweep(z,2,sqrt(n),'*'),na.rm=TRUE)/sqrt(as.vector((!is.na(z))%*%n))}
z<-zmat(ps);ma<-st[keys,'arm']=='microarray';rna<-!ma
# Historical two-stage weighting, on the corrected study models and fixed gene family.
za<-combine(z[,ma,drop=FALSE],st[keys[ma],'N']);zr<-combine(z[,rna,drop=FALSE],st[keys[rna],'N'])
Nma<-sum(st[keys[ma],'N']);Nrna<-sum(st[keys[rna],'N'])
two<-(sqrt(Nma)*za+sqrt(Nrna)*zr)/sqrt(Nma+Nrna)
qps<-ps;for(s in keys[rna])qps[[s]]<-read.csv(file.path(ROOT,'sensitivity',paste0(s,'_quality_weights_DE.csv')))
qz<-combine(zmat(qps),st[keys,'N_effective'])
out<-list();for(n in c('previous_two_stage','RNA_quality_weights')) {
 zz<-if(n=='previous_two_stage')two else qz
 d<-data.frame(GeneSymbol=g,z=zz,p=2*pnorm(-abs(zz)));d$FDR<-p.adjust(d$p,'BH')
 write.csv(d,file.path(ROOT,'sensitivity',paste0('pooled_',n,'.csv')),row.names=FALSE)
 out[[n]]<-data.frame(analysis=n,genes=nrow(d),DE=sum(d$FDR<=.05),up=sum(d$FDR<=.05&d$z>0),down=sum(d$FDR<=.05&d$z<0),baseline_retained=sum(d$FDR<=.05&primary$FDR<=.05))
}
# Modified Hartung-Knapp interval with DL tau2, using moderated study variances as a sensitivity.
fc<-se<-z*NA_real_;for(s in keys){ii<-match(g,ps[[s]]$GeneSymbol);fc[,s]<-ps[[s]]$logFC[ii];se[,s]<-ps[[s]]$SE[ii]}
rw<-1/(se^2+primary$tau2);mu<-rowSums(rw*fc,na.rm=TRUE)/rowSums(rw,na.rm=TRUE);df<-rowSums(is.finite(fc))-1
scale<-pmax(1,rowSums(rw*(fc-mu)^2,na.rm=TRUE)/df);ses<-sqrt(scale/rowSums(rw,na.rm=TRUE))
d<-data.frame(GeneSymbol=g,logFC=mu,SE=ses,df=df,CI_low=mu-qt(.975,df)*ses,CI_high=mu+qt(.975,df)*ses,p=2*pt(-abs(mu/ses),df));d$FDR<-p.adjust(d$p,'BH');d$I2<-primary$I2
write.csv(d,file.path(ROOT,'sensitivity/random_effects_modified_HK.csv'),row.names=FALSE)
write.csv(do.call(rbind,out),file.path(ROOT,'sensitivity/extended_analysis_summary.csv'),row.names=FALSE)
print(do.call(rbind,out));cat('Modified HK DE',sum(d$FDR<=.05),'\nFINAL_SENSITIVITIES_COMPLETE\n')
