Sys.setenv(LOAD_FUNCTIONS_ONLY='1')
source('<workspace>/subtypes_2026-09-08/corrected_baseline/scripts/fit_studies.R')
emd<-read.csv(file.path(ROOT,'inputs/external/external_metadata.csv'))
for(study in unique(emd$study)) {
 sm<-emd[emd$study==study,];m<-read_matrix(sm$expr_path[1]);m<-harmonize(m,study,TRUE)
 fit_one(study,m,sm,TRUE,label=paste0(study,'_external'))
 if(study=='GSE248784') {
  same<-sm[sm$reference_features==33538,]
  fit_one(study,m,same,TRUE,label='GSE248784_external_shared_reference')
  mm<-read_matrix(file.path(ROOT,'inputs/external/GSE248784_mt10_counts.csv'));mm<-harmonize(mm,'GSE248784_mt10',TRUE)
  fit_one(study,mm,sm,TRUE,label='GSE248784_external_mt10')
 }
}
cat('EXTERNAL_FITS_COMPLETE\n')
