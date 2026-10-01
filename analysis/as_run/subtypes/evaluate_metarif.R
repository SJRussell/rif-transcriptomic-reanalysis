.libPaths(c('<r_library>',.libPaths()))
suppressPackageStartupMessages({library(caret);library(nnet);library(randomForest);library(edgeR);library(limma)})
ROOT<-'<workspace>/subtypes_2026-09-08';OLD<-file.path(ROOT,'corrected_baseline')
for(d in c('audit','results','expression','logs','figures'))dir.create(file.path(ROOT,d),showWarnings=FALSE,recursive=TRUE)
set.seed(20260908)
e<-new.env(parent=globalenv());load(file.path(ROOT,'reference/MetaRIF/MetaRIF_classifier_gene_plus.rda'),envir=e)
code<-parse(file.path(ROOT,'reference/MetaRIF/MetaRIF.R'))
for(z in code)if(is.call(z)&&as.character(z[[1]])%in%c('<-','=')&&as.character(z[[2]])%in%c('MSig_pred_factor','MSig_pred_prob','MetaRIF'))eval(z,envir=e)
hg<-read.delim(file.path(OLD,'inputs/hgnc_complete_set.txt'),check.names=FALSE,colClasses='character',quote='',na.strings='')
hg<-hg[hg$status=='Approved'&!is.na(hg$entrez_id)&!duplicated(hg$symbol),];entrez<-setNames(hg$entrez_id,hg$symbol)
feature_ids<-unique(c(e$gene1$ENTREZID,e$gene2$ENTREZID));feature_symbols<-hg$symbol[hg$entrez_id%in%feature_ids]
write.csv(hg[hg$entrez_id%in%feature_ids,c('symbol','entrez_id','ensembl_gene_id')],file.path(ROOT,'audit/model_feature_hgnc.csv'),row.names=FALSE)
training<-unique(c(rownames(e$best_model_RIF11$trainingData),rownames(e$best_model_RIF22$trainingData)))
td<-merge(data.frame(sample=rownames(e$best_model_RIF11$trainingData),model1_class=e$best_model_RIF11$trainingData$.outcome),data.frame(sample=rownames(e$best_model_RIF22$trainingData),model2_class=e$best_model_RIF22$trainingData$.outcome),all=TRUE)
write.csv(td,file.path(ROOT,'audit/published_training_sample_ids.csv'),row.names=FALSE)
published<-c('GSE111974','GSE58144','GSE106602','GSE71331','HRA010876')
to_entrez<-function(m){ids<-unname(entrez[rownames(m)]);keep<-!is.na(ids)&ids!='';avereps(m[keep,,drop=FALSE],ID=ids[keep])}
predict_scaled<-function(z) {
 make_x<-function(ids){x<-matrix(0,nrow=ncol(z),ncol=length(ids),dimnames=list(colnames(z),paste0('gene_',ids)));i<-match(ids,rownames(z));good<-!is.na(i);x[,good]<-t(z[i[good],,drop=FALSE]);x[!is.finite(x)]<-0;as.data.frame(x)}
 p1<-predict(e$best_model_RIF11,newdata=make_x(e$gene1$ENTREZID),type='prob')[,'cluster_RIF1']
 p2<-predict(e$best_model_RIF22,newdata=make_x(e$gene2$ENTREZID),type='prob')[,'cluster_RIF2']
 pred<-ifelse(p1>0.5&p1>p2+.1,'RIF1',ifelse(p2>0.5&p2>p1+.1,'RIF2',ifelse(p1<.5&p2<.5,'Normal','Uncertain')))
 data.frame(sample=colnames(z),p_RIFI=p1,p_RIFM=p2,pred=pred,RIF_prob=ifelse(pred%in%c('Normal','RIF1'),p1,p2),max_RIF_prob=pmax(p1,p2))
}
# Preserve published two-pass row scaling.
scale_twice<-function(m){z<-t(scale(t(m)));t(scale(t(z)))}
pred_fun<-function(m)predict_scaled(scale_twice(m))
auc<-function(y,p){r<-rank(p);a<-sum(y);b<-sum(!y);if(!a||!b)return(NA_real_);(sum(r[y])-a*(a+1)/2)/(a*b)}
ci_auc<-function(y,p){if(min(table(y))<5)return(c(NA,NA));a<-which(y);b<-which(!y);ci<-quantile(replicate(2000,{j<-c(sample(a,length(a),TRUE),sample(b,length(b),TRUE));auc(y[j],p[j])}),c(.025,.975),na.rm=TRUE);if(diff(ci)==0)c(NA_real_,NA_real_) else ci}
md<-read.csv(file.path(OLD,'inputs/master_metadata.csv'));studies<-c(unique(md$study),'GSE248784_external','GSE243550_external')
if(file.exists(file.path(OLD,'expression/Xu2026_external.rds')))studies<-c(studies,'Xu2026_external')
if(file.exists(file.path(OLD,'expression/GSE307510_external.rds')))studies<-c(studies,'GSE307510_external')
sets<-readRDS(file.path(OLD,'pathways/fixed_programs.rds'))
predictions<-coverage<-contexts<-original_tests<-summaries<-scores<-list();object_cache<-list()
for(s in studies){
 obj<-readRDS(file.path(OLD,'expression',paste0(s,'.rds')));m<-obj$expr;sm<-obj$metadata
 stopifnot(identical(as.character(sm$sample),colnames(m)))
 m0<-if(!is.null(obj$counts))cpm(obj$counts,log=FALSE,normalized.lib.sizes=FALSE) else if(!is.null(obj$classifier_expr))obj$classifier_expr else m
 em<-to_entrez(m0);em<-em[intersect(feature_ids,rownames(em)),,drop=FALSE]
 avail<-intersect(feature_ids,rownames(em));missing<-setdiff(feature_ids,avail)
 role<-if(sub('_external','',s)%in%published)'published_development' else if(s=='GSE248784_external')'external_previously_examined' else if(s=='GSE243550_external')'technical_clinical_link_unresolved' else if(s=='Xu2026_external')'repository_linked_validation_provenance_unresolved' else if(s=='GSE307510_external')'small_broader_external_exploratory' else 'outside_published_development'
 n1<-sum(e$gene1$ENTREZID%in%avail);n2<-sum(e$gene2$ENTREZID%in%avail)
 eligible<-n1/64>=.9&n2/19>=.9
 coverage[[s]]<-data.frame(study=s,role=role,n=ncol(m),n_control=sum(tolower(sm$condition)=='control'),n_rif=sum(tolower(sm$condition)=='rif'),model1_present=n1,model1_total=64,model2_present=n2,model2_total=19,coverage_pass=eligible,missing_entrez=paste(missing,collapse=';'),exact_training_id_overlap=sum(sm$sample%in%training),input=if(is.null(obj$counts))'source log expression' else 'unlogged CPM; gene-length cancels under row scaling',preprocessing_note=if(s=='GSE58144')'raw sample/reference rebuild; case-blind batch removal for classifier; batch and day in DE design' else '')
 pr<-pred_fun(em);pr$study<-s;pr$condition<-tolower(sm$condition);pr$role<-role;pr$coverage_pass<-eligible;pr$variant<-'released_scale';predictions[[s]]<-pr
 pub<-tryCatch(list(pred=e$MSig_pred_factor(scale_twice(em)),prob=e$MSig_pred_prob(scale_twice(em))),error=function(er)er)
 if(inherits(pub,'error'))original_tests[[s]]<-data.frame(study=s,status='published_helpers_error',max_prob_difference=NA,error=conditionMessage(pub)) else original_tests[[s]]<-data.frame(study=s,status=if(identical(unname(pub$pred),unname(pr$pred))&&max(abs(pub$prob$RIF_prob-pr$RIF_prob))<1e-12)'identical' else 'DIFFERENT',max_prob_difference=max(abs(pub$prob$RIF_prob-pr$RIF_prob)),error='')
 if(!is.null(obj$counts)){
  alt<-pred_fun(to_entrez(m));alt$study<-s;alt$condition<-pr$condition;alt$role<-role;alt$coverage_pass<-eligible;alt$variant<-'logCPM_sensitivity';predictions[[paste0(s,'_log')]]<-alt
 }
 y<-pr$condition=='rif';ci<-ci_auc(y,pr$RIF_prob)
 summaries[[s]]<-data.frame(study=s,role=role,coverage_pass=eligible,n_control=sum(!y),n_rif=sum(y),AUC=auc(y,pr$RIF_prob),CI_low=ci[1],CI_high=ci[2],AUC_max_score=auc(y,pr$max_RIF_prob),sensitivity=sum(y&pr$pred%in%c('RIF1','RIF2'))/sum(y),specificity=sum(!y&pr$pred=='Normal')/sum(!y),RIF_I=sum(y&pr$pred=='RIF1'),RIF_M=sum(y&pr$pred=='RIF2'),RIF_Normal=sum(y&pr$pred=='Normal'),RIF_Uncertain=sum(y&pr$pred=='Uncertain'),control_I=sum(!y&pr$pred=='RIF1'),control_M=sum(!y&pr$pred=='RIF2'))
 # Drop one other subject: predict remaining subjects with the released normalization.
 if(ncol(em)>=6){omit<-if(ncol(em)<=30)seq_len(ncol(em)) else sort(sample(seq_len(ncol(em)),30));chg<-tot<-rep(0,ncol(em));pd<-rep(0,ncol(em))
  for(j in omit){k<-setdiff(seq_len(ncol(em)),j);q<-pred_fun(em[,k,drop=FALSE]);chg[k]<-chg[k]+(q$pred!=pr$pred[k]);tot[k]<-tot[k]+1;pd[k]<-pmax(pd[k],abs(q$RIF_prob-pr$RIF_prob[k]))}
  contexts[[s]]<-data.frame(study=s,sample=sm$sample,condition=pr$condition,baseline=pr$pred,contexts=tot,changed=chg,change_fraction=chg/tot,max_probability_change=pd)
 }
 score_m<-if(!is.null(obj$classifier_expr))obj$classifier_expr else m
 ze<-t(scale(t(score_m)));ze[!is.finite(ze)]<-0
 for(v in c('all_genes','exclude_model_genes'))for(p in names(sets)){
  ids<-sets[[p]];if(v=='exclude_model_genes')ids<-setdiff(ids,feature_symbols)
  found<-intersect(ids,rownames(ze));if(length(found)<5)next
  scores[[paste(s,p,v)]]<-data.frame(study=s,sample=sm$sample,condition=pr$condition,program=p,variant=v,score=colMeans(ze[found,,drop=FALSE]),genes_present=length(found),genes_total=length(ids),coverage=length(found)/length(ids))
 }
 obj$entrez_model_input<-em;obj$prediction<-pr;object_cache[[s]]<-obj
 cat('PREDICTED',s,ncol(m),'coverage',n1,'/',n2,'AUC',round(auc(y,pr$RIF_prob),3),'\n');flush.console()
}
write.csv(do.call(rbind,coverage),file.path(ROOT,'audit/cohort_model_coverage.csv'),row.names=FALSE)
write.csv(do.call(rbind,original_tests),file.path(ROOT,'audit/published_helper_equivalence.csv'),row.names=FALSE)
write.csv(do.call(rbind,predictions),file.path(ROOT,'results/predictions.csv'),row.names=FALSE)
write.csv(do.call(rbind,summaries),file.path(ROOT,'results/classifier_by_cohort.csv'),row.names=FALSE)
write.csv(do.call(rbind,contexts),file.path(ROOT,'results/cohort_context_sensitivity.csv'),row.names=FALSE)
write.csv(do.call(rbind,scores),file.path(ROOT,'results/fixed_program_scores.csv'),row.names=FALSE)
saveRDS(object_cache,file.path(ROOT,'expression/cohorts_for_subtyping.rds'))
capture.output(sessionInfo(),file=file.path(ROOT,'logs/metarif_sessionInfo.txt'))
cat('METARIF_COMPLETE\n')
