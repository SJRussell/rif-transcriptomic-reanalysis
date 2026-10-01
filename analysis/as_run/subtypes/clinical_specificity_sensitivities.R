.libPaths(c('<r_library>',.libPaths()))
ROOT<-'<workspace>/subtypes_2026-09-08'
p<-read.csv(file.path(ROOT,'results/predictions.csv'));p<-p[p$variant=='released_scale'&p$role=='outside_published_development',]
omissions<-list(without_Bui='Bui2024',without_Bui_and_oocyte_controls=c('Bui2024','GSE92324'))
rr<-list()
for(n in names(omissions))for(state in c('RIF1','RIF2','Either_RIF')){
 d<-p[!p$study%in%omissions[[n]],]
 tab<-lapply(split(d,d$study),function(q){y<-q$condition=='rif';z<-if(state=='Either_RIF')q$pred%in%c('RIF1','RIF2') else q$pred==state;matrix(c(sum(y&z),sum(!y&z),sum(y&!z),sum(!y&!z)),2,2)})
 f<-mantelhaen.test(array(unlist(tab),c(2,2,length(tab))),correct=FALSE)
 rr[[paste(n,state)]]<-data.frame(scope=n,state=state,studies=length(tab),n_rif=sum(d$condition=='rif'),n_control=sum(d$condition=='control'),odds_ratio=unname(f$estimate),CI_low=f$conf.int[1],CI_high=f$conf.int[2],p=f$p.value)
}
r<-do.call(rbind,rr);r$FDR<-ave(r$p,r$scope,FUN=function(x)p.adjust(x,'BH'))
write.csv(r,file.path(ROOT,'results/clinical_enrichment_sensitivities.csv'),row.names=FALSE)
print(r)
