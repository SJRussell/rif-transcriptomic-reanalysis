.libPaths(c('<r_library>',.libPaths()))
suppressPackageStartupMessages({library(limma);library(edgeR);library(fgsea);library(BiocParallel)})
ROOT<-'<workspace>/subtypes_2026-09-08/corrected_baseline'
oldsets<-'<results_package_earlier>/analysis/downstream/msigdb_sets_used.rds'
sets<-readRDS(oldsets);saveRDS(sets,file.path(ROOT,'inputs/msigdb_sets_used.rds'))
lookup<-read.csv(file.path(ROOT,'inputs/hgnc_lookup.csv'));mapped<-lookup$GeneSymbol[match(sets$gene_symbol,lookup$input_id)]
sets$approved_symbol<-mapped;sets<-sets[!is.na(mapped)&mapped!='',]
pathways<-lapply(split(sets$approved_symbol,sets$gs_name),unique)
saveRDS(pathways,file.path(ROOT,'pathways/approved_gene_sets.rds'))
run_fg<-function(d,label,zcol='z') {
 ranks<-setNames(d[[zcol]],d$GeneSymbol);ranks<-sort(ranks[is.finite(ranks)],decreasing=TRUE)
 set.seed(20260908)
 fg<-fgseaMultilevel(pathways=pathways,stats=ranks,minSize=15,maxSize=500,nPermSimple=10000,eps=1e-10,nproc=1,BPPARAM=SerialParam())
 fg<-as.data.frame(fg);fg$leadingEdge<-vapply(fg$leadingEdge,paste,collapse=';',character(1));fg<-fg[order(fg$padj),]
 write.csv(fg,file.path(ROOT,'pathways',paste0('fgsea_',label,'.csv')),row.names=FALSE)
 cat('FGSEA',label,'tested',nrow(fg),'significant',sum(fg$padj<=.05,na.rm=TRUE),'\n')
}
primary<-read.csv(file.path(ROOT,'pooled/primary.csv'))
for(label in c('primary','microarray','rnaseq','without_GSE111974','bulk_only')) {
 if(label %in% c('primary','microarray','rnaseq')) d<-primary else d<-read.csv(file.path(ROOT,'sensitivity',if(label=='bulk_only')'leaveout_all_pseudobulk.csv' else 'leaveout_GSE111974.csv'))
 run_fg(d,label,if(label %in% c('microarray','rnaseq'))paste0('z_',label) else 'z')
}
# Fixed biological programs: Hallmarks plus cilia and cell-lineage markers; defined before this rerun's enrichment inspection.
programs<-pathways[grepl('^HALLMARK_',names(pathways))]
programs[['CILIARY_MOTILITY']]<-c('FOXJ1','PIFO','TPPP3','RFX2','RFX3','DNAH5','DNAH9','DNAH11','CCDC39','CCDC40','SPAG6','SPAG17','HYDIN','CFAP43','CFAP44','RSPH1','RSPH4A','RSPH9','DNAAF1','DNAAF2')
programs[['NK_CYTOTOXIC']]<-c('NKG7','GNLY','KLRD1','KLRF1','NCR1','PRF1','GZMB','GZMH','CTSW','XCL1','XCL2','IL2RB')
programs[['MYELOID_LINEAGE']]<-c('LST1','AIF1','TYROBP','FCER1G','CTSS','SPI1','CD68','LYZ','CSF1R','C1QA','C1QB','C1QC')
programs[['T_CELL_LINEAGE']]<-c('CD3D','CD3E','CD3G','TRAC','TRBC1','TRBC2','LCK','CD2','CD247','IL7R')
programs[['EPITHELIAL_LINEAGE']]<-c('EPCAM','KRT8','KRT18','KRT19','CDH1','CLDN3','CLDN4','CLDN7','MUC1','TACSTD2')
programs[['STROMAL_LINEAGE']]<-c('DCN','LUM','COL1A1','COL1A2','COL3A1','PDGFRA','COL6A1','COL6A2','COL6A3','DPT')
programs[['DECIDUALIZATION']]<-c('IGFBP1','PRL','FOXO1','HAND2','WNT4','BMP2','CEBPB','IL15','DKK1','PROK1')
saveRDS(programs,file.path(ROOT,'pathways/fixed_programs.rds'))
write.csv(do.call(rbind,lapply(names(programs),function(n)data.frame(program=n,GeneSymbol=programs[[n]]))),file.path(ROOT,'pathways/fixed_program_membership.csv'),row.names=FALSE)
md<-read.csv(file.path(ROOT,'inputs/master_metadata.csv'));results<-list();scores<-list()
for(s in unique(md$study)) {
 obj<-readRDS(file.path(ROOT,'expression',paste0(s,'.rds')));expr<-obj$expr;design<-obj$design
 if(!is.null(obj$counts)) {y<-DGEList(obj$counts[rownames(expr),]);y<-calcNormFactors(y);v<-voom(y,design,plot=FALSE);em<-v} else em<-expr
 inds<-ids2indices(programs,rownames(expr),remove.empty=TRUE);inds<-inds[lengths(inds)>=5]
 ca<-camera(em,index=inds,design=design,contrast=2,inter.gene.cor=NA,sort=FALSE)
 ca$program<-rownames(ca);ca$study<-s;ca$z<-qnorm(ca$PValue/2,lower.tail=FALSE)*ifelse(ca$Direction=='Up',1,-1)
 results[[s]]<-ca
 # Mean standardized expression scores describe programs; these are not cell fractions.
 ze<-t(scale(t(expr)));ze[!is.finite(ze)]<-0
 for(p in names(inds))scores[[paste(s,p)]]<-data.frame(study=s,sample=obj$metadata$sample,condition=obj$metadata$condition,program=p,score=colMeans(ze[inds[[p]],,drop=FALSE]),genes=length(inds[[p]]))
 cat('CAMERA',s,nrow(ca),'programs\n')
}
ca<-do.call(rbind,results);write.csv(ca,file.path(ROOT,'pathways/programs_per_study_camera.csv'),row.names=FALSE)
write.csv(do.call(rbind,scores),file.path(ROOT,'pathways/program_scores_per_sample.csv'),row.names=FALSE)
st<-read.csv(file.path(ROOT,'audit/study_summary.csv'));ww<-setNames(sqrt(st$N_effective),st$study)
po<-do.call(rbind,lapply(split(ca,ca$program),function(d){z<-sum(ww[d$study]*d$z)/sqrt(sum(ww[d$study]^2));i<-d$study!='GSE111974';zloo<-sum(ww[d$study[i]]*d$z[i])/sqrt(sum(ww[d$study[i]]^2));data.frame(program=d$program[1],studies=nrow(d),z=z,p=2*pnorm(-abs(z)),same_direction=sum(sign(d$z)==sign(z)),z_without_GSE111974=zloo,p_without_GSE111974=2*pnorm(-abs(zloo)))}))
po$FDR<-p.adjust(po$p,'BH');po$FDR_without_GSE111974<-p.adjust(po$p_without_GSE111974,'BH');po<-po[order(po$p),]
write.csv(po,file.path(ROOT,'pathways/programs_meta_camera.csv'),row.names=FALSE)
print(head(po,15));cat('PATHWAYS_COMPLETE\n')
