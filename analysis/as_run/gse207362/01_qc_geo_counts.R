#!/usr/bin/env Rscript
# GSE207362 external evaluation — step 1: verify mapped columns, species sanity, library QC.
# Input: GEO processed FeatureCounts table (human), published subset mapping.
# Output: outputs/GSE207362_QC_table.csv, outputs/GSE207362_counts_38.csv.gz, logs/01_qc_sessionInfo.txt
suppressPackageStartupMessages({library(data.table); library(edgeR); library(jsonlite)})
wd <- "<cluster_scratch>/gse207362_external_20260918"
cnt <- fread(cmd = paste("zcat", file.path(wd, "geo/GSE207362_Human_FeatureCounts.txt.gz")))
setnames(cnt, 1, "GeneSymbol")
map <- fread(file.path(wd, "inputs/GSE207362_published_subset_mapping.csv"))
cat("count table:", nrow(cnt), "genes x", ncol(cnt)-1, "libraries\n")
# 1. all mapped columns exist and positions agree with the 1-based expression column
stopifnot(all(map$count_column %in% names(cnt)))
pos <- match(map$count_column, names(cnt)) - 1L
stopifnot(all(pos == map$expression_column_1based))
cat("all", nrow(map), "mapped columns present; positions agree with mapping\n")
stopifnot(sum(map$publication_group=="rif")==12, sum(map$publication_group=="control")==26)
# 2. species sanity: human symbols are upper-case; mouse would be title-case
sym <- cnt$GeneSymbol
frac_upper <- mean(grepl("^[A-Z0-9][A-Z0-9.-]*$", sym))
n_titlecase <- sum(grepl("^[A-Z][a-z]", sym))
cat(sprintf("symbols: %d; fraction all-upper %.3f; title-case (mouse-like) %d\n", length(sym), frac_upper, n_titlecase))
human_markers <- c("GAPDH","ACTB","PAEP","LIF","HOXA10","EZH2","XIST","RPS4Y1")
cat("human marker symbols present:", paste(intersect(human_markers, sym), collapse=","), "\n")
cat("duplicated symbols:", sum(duplicated(sym)), "\n")
# 2b. GEO deposit defect: rows with NA counts cannot be used; record and drop (gene rows only, never libraries)
na_rows <- sym[!complete.cases(cnt[, -1])]
cat("gene rows with NA counts (dropped):", paste(na_rows, collapse=","), "\n")
writeLines(na_rows, file.path(wd, "outputs/GSE207362_dropped_NA_gene_rows.txt"))
cnt <- cnt[complete.cases(cnt[, -1])]; sym <- cnt$GeneSymbol
# 3. library QC on the 38 analysed libraries (and all 67 for reference)
M_all <- as.matrix(cnt[, -1]); rownames(M_all) <- sym
qc_all <- data.table(library = colnames(M_all),
                     library_size = colSums(M_all),
                     detected_genes_gt0 = colSums(M_all > 0),
                     detected_genes_ge10 = colSums(M_all >= 10))
qc_all <- merge(qc_all, map[, .(library = count_column, publication_group, publication_index)],
                by = "library", all.x = TRUE)
qc_all[, analysed_in_paper := !is.na(publication_group)]
qc_all[, column_position := match(library, colnames(M_all))]
setorder(qc_all, column_position)
M38 <- M_all[, map$count_column]
grp <- factor(map$publication_group, levels = c("control","rif"))
keep <- filterByExpr(M38, group = grp)
cat("filterByExpr keeps", sum(keep), "of", nrow(M38), "genes (38-library design)\n")
qc_all[analysed_in_paper == TRUE, genes_kept_filterByExpr := sum(keep)]
dge <- calcNormFactors(DGEList(M38, group = grp), method = "TMM")
qc_all[analysed_in_paper == TRUE, TMM_norm_factor := dge$samples$norm.factors[match(library, colnames(M38))]]
qc_all[analysed_in_paper == TRUE, frac_top_gene := apply(M38, 2, function(x) max(x)/sum(x))[match(library, colnames(M38))]]
print(qc_all[analysed_in_paper == TRUE, .(n = .N, min_libsize = min(library_size), median_libsize = median(library_size),
         max_libsize = max(library_size), min_detected = min(detected_genes_gt0), median_detected = median(detected_genes_gt0)),
     by = publication_group])
dir.create(file.path(wd, "outputs"), showWarnings = FALSE)
fwrite(qc_all, file.path(wd, "outputs/GSE207362_QC_table.csv"))
out38 <- data.table(GeneSymbol = sym, as.data.table(M38))
fwrite(out38, file.path(wd, "outputs/GSE207362_counts_38_raw_symbols.csv"))
writeLines(capture.output(sessionInfo()), file.path(wd, "logs/01_qc_sessionInfo.txt"))
cat("done\n")
