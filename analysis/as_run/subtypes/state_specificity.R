.libPaths(c('<r_library>',.libPaths()))
suppressPackageStartupMessages(library(metafor))
ROOT<-'<workspace>/subtypes_2026-09-08'
p<-read.csv(file.path(ROOT,'results/predictions.csv'));p<-p[p$variant=='released_scale',]
rows<-list()
for(scope in c('outside_published_development','published_development','external_previously_examined'))for(state in c('RIF1','RIF2','Either_RIF')){
 d<-p[p$role==scope,];if(!nrow(d))next
 tables<-lapply(split(d,d$study),function(q){y<-q$condition=='rif';z<-if(state=='Either_RIF')q$pred%in%c('RIF1','RIF2') else q$pred==state;matrix(c(sum(y&z),sum(!y&z),sum(y&!z),sum(!y&!z)),2,2)})
 ar<-array(unlist(tables),c(2,2,length(tables)));f<-tryCatch(mantelhaen.test(ar,correct=FALSE),error=function(e)NULL)
 if(is.null(f))next
 rows[[paste(scope,state)]]<-data.frame(scope=scope,state=state,studies=length(tables),n=nrow(d),odds_ratio=unname(f$estimate),CI_low=f$conf.int[1],CI_high=f$conf.int[2],p=f$p.value)
}
r<-do.call(rbind,rows);r$FDR<-ave(r$p,r$scope,FUN=function(x)p.adjust(x,'BH'))
write.csv(r,file.path(ROOT,'results/clinical_enrichment_by_model_output.csv'),row.names=FALSE)
# Modified Hartung-Knapp sensitivity for all fixed score contrasts, with explicit BH families.
d<-read.csv(file.path(ROOT,'results/program_contrasts_per_study.csv'))
po<-read.csv(file.path(ROOT,'results/program_contrasts_meta.csv'));po$modified_HK_effect<-po$modified_HK_SE<-po$modified_HK_p<-NA_real_
for(i in seq_len(nrow(po))){
 q<-d[d$variant==po$variant[i]&d$contrast==po$contrast[i]&d$program==po$program[i],]
 q<-q[if(po$scope[i]=='original_13')!grepl('_external',q$study) else if(po$scope[i]=='outside_MetaRIF_development')q$role=='outside_published_development' else q$study=='GSE248784_external',]
 if(nrow(q)<3)next
 f<-tryCatch(rma.uni(yi=q$effect,sei=q$SE,method='REML',test='z'),error=function(e)NULL);if(is.null(f))next
 w<-1/(q$SE^2+f$tau2);mu<-as.numeric(f$b);df<-nrow(q)-1;scale<-max(1,sum(w*(q$effect-mu)^2)/df);se<-sqrt(scale/sum(w))
 po$modified_HK_effect[i]<-mu;po$modified_HK_SE[i]<-se;po$modified_HK_p[i]<-2*pt(-abs(mu/se),df)
}
po$RE_FDR<-ave(po$RE_p,interaction(po$scope,po$variant,po$contrast),FUN=function(x)p.adjust(x,'BH'))
po$modified_HK_FDR<-ave(po$modified_HK_p,interaction(po$scope,po$variant,po$contrast),FUN=function(x)p.adjust(x,'BH'))
write.csv(po,file.path(ROOT,'results/program_contrasts_meta_robust.csv'),row.names=FALSE)
# Compare RIF cases with controls assigned to the SAME model output.
# This remains exploratory: expression-derived groups are not independently observed subtypes.
objs<-readRDS(file.path(ROOT,'expression/cohorts_for_subtyping.rds'))
scores<-read.csv(file.path(ROOT,'results/fixed_program_scores.csv'));scores<-scores[scores$variant=='exclude_model_genes',]
rs<-list()
for(s in names(objs))for(state in c('RIF1','RIF2','Normal')){
 o<-objs[[s]];pp<-p[p$study==s,];pp<-pp[match(o$metadata$sample,pp$sample),];k<-pp$pred==state;grp<-as.numeric(pp$condition[k]=='rif')
 if(sum(grp)<3||sum(grp==0)<3)next
 X<-cbind(Intercept=1,group=grp)
 if(ncol(o$design)>2)for(j in 3:ncol(o$design)){cv<-o$design[k,j];if(sd(cv)>0&&qr(cbind(X,cv))$rank>ncol(X))X<-cbind(X,cv)}
 if(nrow(X)<=ncol(X)+1)next
 sc<-scores[scores$study==s,]
 for(prog in unique(sc$program)){
  q<-sc[sc$program==prog,];q<-q[match(pp$sample,q$sample),];if(anyNA(q$score)||min(q$coverage)<.8)next
  f<-lm.fit(X,q$score[k]);df<-nrow(X)-ncol(X);se<-sqrt(diag(solve(crossprod(X)))*sum(f$residuals^2)/df)[2]
  if(!is.finite(se)||se<=0)next
  b<-f$coefficients[2];rs[[paste(s,state,prog)]]<-data.frame(study=s,role=pp$role[1],state=state,program=prog,n_case=sum(grp),n_control=sum(grp==0),effect=unname(b),SE=unname(se),p=2*pt(-abs(b/se),df))
 }
}
rr<-do.call(rbind,rs);rr$FDR<-ave(rr$p,interaction(rr$study,rr$state),FUN=function(x)p.adjust(x,'BH'))
write.csv(rr,file.path(ROOT,'results/program_contrasts_same_output_controls.csv'),row.names=FALSE)
print(r);cat('STATE_SPECIFICITY_COMPLETE\n')
