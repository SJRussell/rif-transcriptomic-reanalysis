.libPaths(c('<r_library>',.libPaths()))
ROOT<-'<workspace>/subtypes_2026-09-08/corrected_baseline'
programs<-readRDS(file.path(ROOT,'pathways/fixed_programs.rds'))
o<-readRDS(file.path(ROOT,'expression/GSE250130_mixed_time.rds'))
mapping<-read.csv('<workspace>/sdrive_audit_20260908/metadata_audit/GSE250130_barcode_mapping.csv')
group<-mapping$GEO_group[match(o$metadata$sample,mapping$archive_sample)]
print(table(group));stopifnot(!anyNA(group))
day<-as.numeric(gsub('[^0-9]','',group));day[o$metadata$condition=='RIF']<-7
expr<-o$expr;ze<-t(scale(t(expr)));ze[!is.finite(ze)]<-0
scores<-list();tests<-list()
for(p in c('CILIARY_MOTILITY','NK_CYTOTOXIC','MYELOID_LINEAGE','EPITHELIAL_LINEAGE','STROMAL_LINEAGE','DECIDUALIZATION','HALLMARK_APICAL_JUNCTION')) {
 ids<-intersect(programs[[p]],rownames(expr));sc<-colMeans(ze[ids,,drop=FALSE]);control<-o$metadata$condition=='Control'
 # Healthy-cycle trajectory: controls only. No disease labels are used to estimate this slope.
 fit<-lm(sc[control]~day[control]);cc<-summary(fit)$coefficients
 tests[[p]]<-data.frame(program=p,n_control=sum(control),genes=length(ids),slope_per_LHday=cc[2,1],SE=cc[2,2],p=cc[2,4])
 scores[[p]]<-data.frame(sample=o$metadata$sample,condition=o$metadata$condition,group=group,LHday=day,program=p,score=sc,genes=length(ids))
}
tt<-do.call(rbind,tests);tt$FDR<-p.adjust(tt$p,'BH');write.csv(tt,file.path(ROOT,'sensitivity/healthy_cycle_programs.csv'),row.names=FALSE)
write.csv(do.call(rbind,scores),file.path(ROOT,'sensitivity/healthy_cycle_program_scores.csv'),row.names=FALSE)
print(tt);cat('TIMING_PROGRAMS_COMPLETE\n')
