.libPaths(c('<r_library>',.libPaths()))
suppressPackageStartupMessages({library(mclust);library(diptest);library(metafor)})
ROOT<-'<workspace>/subtypes_2026-09-08';set.seed(20260908)
objs<-readRDS(file.path(ROOT,'expression/cohorts_for_subtyping.rds'))
pr<-read.csv(file.path(ROOT,'results/predictions.csv'));pr<-pr[pr$variant=='released_scale',]
ps<-read.csv(file.path(ROOT,'results/fixed_program_scores.csv'))
immune<-c('HALLMARK_INFLAMMATORY_RESPONSE','HALLMARK_INTERFERON_GAMMA_RESPONSE','HALLMARK_TNFA_SIGNALING_VIA_NFKB')
metabolic<-c('HALLMARK_OXIDATIVE_PHOSPHORYLATION','HALLMARK_FATTY_ACID_METABOLISM','HALLMARK_GLYCOLYSIS')
axes<-list()
for(v in unique(ps$variant)){
 d<-ps[ps$variant==v,];r<-reshape(d[,c('study','sample','condition','program','score')],idvar=c('study','sample','condition'),timevar='program',direction='wide')
 names(r)<-sub('^score[.]','',names(r));r$immune_axis<-rowMeans(r[,immune]);r$metabolic_axis<-rowMeans(r[,metabolic]);r$variant<-v
 axes[[v]]<-merge(pr[,c('study','sample','pred','p_RIFI','p_RIFM','RIF_prob','role')],r,by=c('study','sample'),all.x=TRUE)
}
aa<-do.call(rbind,axes);write.csv(aa,file.path(ROOT,'results/axes_and_predictions.csv'),row.names=FALSE)
# Independently documented covariates are retained for within-study comparisons.
fit_contrast<-function(y,group,design){
 X<-cbind(Intercept=1,group=group)
 if(ncol(design)>2)for(j in 3:ncol(design)){z<-design[,j];if(sd(z)>0&&qr(cbind(X,z))$rank>ncol(X))X<-cbind(X,z)}
 if(nrow(X)<=ncol(X)+1||sd(y)==0)return(NULL)
 f<-lm.fit(X,y);df<-length(y)-f$rank;s2<-sum(f$residuals^2)/df;se<-sqrt(diag(solve(crossprod(X)))*s2)[2]
 if(!is.finite(se)||se<=0)return(NULL)
 b<-f$coefficients[2];p<-2*pt(-abs(b/se),df);z<-qnorm(p/2,lower.tail=FALSE)*sign(b)
 data.frame(effect=b,SE=se,t=b/se,df=df,p=p,z=z,n_case=sum(group==1),n_control=sum(group==0))
}
res<-list();axisres<-list()
for(s in names(objs)){
 obj<-objs[[s]];p<-pr[pr$study==s,];p<-p[match(obj$metadata$sample,p$sample),];stopifnot(!anyNA(p$sample));y<-p$condition=='rif'
 for(v in unique(ps$variant)){
  d<-ps[ps$study==s&ps$variant==v,]
  for(ct in c('all_RIF','RIF_I','RIF_M','RIF_called_Normal')){
   selected<-if(ct=='all_RIF')y else y&p$pred==switch(ct,RIF_I='RIF1',RIF_M='RIF2',RIF_called_Normal='Normal')
   k<-!y|selected;g<-as.numeric(selected[k]);if(sum(g)<3||sum(g==0)<3)next
   for(prog in unique(d$program)){
    q<-d[d$program==prog,];q<-q[match(p$sample,q$sample),];if(anyNA(q$score)||min(q$coverage)<.8)next
    f<-fit_contrast(q$score[k],g,obj$design[k,,drop=FALSE]);if(is.null(f))next
    res[[paste(s,v,ct,prog)]]<-cbind(data.frame(study=s,role=p$role[1],variant=v,contrast=ct,program=prog),f)
   }
  }
 }
}
rr<-do.call(rbind,res);rr$FDR_within_study_contrast<-ave(rr$p,interaction(rr$study,rr$variant,rr$contrast),FUN=function(x)p.adjust(x,'BH'))
write.csv(rr,file.path(ROOT,'results/program_contrasts_per_study.csv'),row.names=FALSE)
pool<-list()
for(scope in c('original_13','outside_MetaRIF_development','external_GSE248784')){
 ds<-rr[if(scope=='original_13')!grepl('_external',rr$study) else if(scope=='outside_MetaRIF_development')rr$role=='outside_published_development' else rr$study=='GSE248784_external',]
 for(d in split(ds,interaction(ds$variant,ds$contrast,ds$program,drop=TRUE))){
  if(nrow(d)<2&&scope!='external_GSE248784')next
  w<-sqrt(4*d$n_case*d$n_control/(d$n_case+d$n_control));z<-sum(w*d$z)/sqrt(sum(w^2));p<-2*pnorm(-abs(z))
  f<-if(nrow(d)>=3)tryCatch(rma.uni(yi=d$effect,sei=d$SE,method='REML',test='knha'),error=function(e)NULL) else NULL
  pool[[paste(scope,d$variant[1],d$contrast[1],d$program[1])]]<-data.frame(scope=scope,variant=d$variant[1],contrast=d$contrast[1],program=d$program[1],studies=nrow(d),n_case=sum(d$n_case),n_control=sum(d$n_control),z=z,p=p,same_direction=sum(sign(d$effect)==sign(z)),RE_effect=if(is.null(f))NA else as.numeric(f$b),RE_CI_low=if(is.null(f))NA else f$ci.lb,RE_CI_high=if(is.null(f))NA else f$ci.ub,RE_p=if(is.null(f))NA else f$pval,I2=if(is.null(f))NA else f$I2)
 }
}
po<-do.call(rbind,pool);po$FDR<-ave(po$p,interaction(po$scope,po$variant,po$contrast),FUN=function(x)p.adjust(x,'BH'))
write.csv(po,file.path(ROOT,'results/program_contrasts_meta.csv'),row.names=FALSE)
# Discreteness checks on fixed axes. These are exploratory distribution diagnostics.
mix<-dips<-cv<-list();base<-aa[aa$variant=='exclude_model_genes',]
for(scope in c('original_13','outside_MetaRIF_development','external_GSE248784'))for(condition in c('rif','control')){
 d<-base[base$condition==condition & if(scope=='original_13')!grepl('_external',base$study) else if(scope=='outside_MetaRIF_development')base$role=='outside_published_development' else base$study=='GSE248784_external',]
 if(nrow(d)<10)next
 x<-as.matrix(d[,c('immune_axis','metabolic_axis')]);x<-scale(x)
 for(k in 1:3){fit<-tryCatch(Mclust(x,G=k,modelNames='VVV',verbose=FALSE),error=function(e)NULL)
  if(is.null(fit))next
  mix[[paste(scope,condition,k)]]<-data.frame(scope=scope,condition=condition,n=nrow(x),components=k,BIC=fit$bic,loglik=fit$loglik,smallest_component=min(tabulate(fit$classification,nbins=k)))
 }
 for(ax in c('immune_axis','metabolic_axis')){f<-dip.test(d[[ax]]);dips[[paste(scope,condition,ax)]]<-data.frame(scope=scope,condition=condition,axis=ax,n=nrow(d),dip=unname(f$statistic),p=f$p.value)}
 if(condition=='rif'&&scope=='outside_MetaRIF_development'){
  for(s in unique(d$study)){
   train<-d[d$study!=s,];test<-d[d$study==s,];mu<-colMeans(train[,c('immune_axis','metabolic_axis')]);sig<-apply(train[,c('immune_axis','metabolic_axis')],2,sd)
   tx<-scale(as.matrix(train[,c('immune_axis','metabolic_axis')]),center=mu,scale=sig);vx<-scale(as.matrix(test[,c('immune_axis','metabolic_axis')]),center=mu,scale=sig)
   for(k in 1:2){f<-tryCatch(Mclust(tx,G=k,modelNames='VVV',verbose=FALSE),error=function(e)NULL);if(is.null(f))next
    ld<-mclust::dens(modelName=f$modelName,data=vx,parameters=f$parameters,logarithm=TRUE)
    cv[[paste(s,k)]]<-data.frame(held_out_study=s,components=k,n_train=nrow(train),n_test=nrow(test),sum_log_density=sum(ld),mean_log_density=mean(ld))
   }
  }
 }
}
write.csv(do.call(rbind,mix),file.path(ROOT,'results/axis_mixture_BIC.csv'),row.names=FALSE)
write.csv(do.call(rbind,dips),file.path(ROOT,'results/axis_unimodality_tests.csv'),row.names=FALSE)
write.csv(do.call(rbind,cv),file.path(ROOT,'results/axis_mixture_cross_cohort.csv'),row.names=FALSE)
# Timing correlations are confined to recorded variable timing; no pseudo-clock is inferred.
tim<-list()
for(s in names(objs)){
 o<-objs[[s]];if(ncol(o$design)<3||!any(grepl('biopsy_day',colnames(o$design))))next
 d<-base[base$study==s,];d<-d[match(o$metadata$sample,d$sample),];day<-o$design[,'biopsy_day']
 for(ax in c('immune_axis','metabolic_axis','HALLMARK_ESTROGEN_RESPONSE_EARLY','CILIARY_MOTILITY'))for(group in c('all','control','rif')){
  k<-if(group=='all')rep(TRUE,nrow(d)) else d$condition==group;if(sum(k)<5||length(unique(day[k]))<2)next
  f<-cor.test(day[k],d[[ax]][k],method='spearman',exact=FALSE)
  tim[[paste(s,ax,group)]]<-data.frame(study=s,group=group,program=ax,n=sum(k),rho=unname(f$estimate),p=f$p.value)
 }
}
write.csv(do.call(rbind,tim),file.path(ROOT,'results/recorded_timing_correlations.csv'),row.names=FALSE)
cat('HETEROGENEITY_COMPLETE\n');print(do.call(rbind,mix));print(do.call(rbind,cv))
