# Fixed raw-validation rules

Specified before inspecting aligned expression or case/control results, 8 September 2026. All 22 BioSample title mappings are retained (13 CON and 9 RIF).

Choose count direction from each library's fraction `forward_assigned / (forward_assigned + reverse_assigned)`. Use forward if every library is at least 0.80; reverse if every library is at most 0.20; otherwise use unstranded, with any inconsistent library investigated as a technical issue. Do not choose direction from DEG counts or classifier performance.

Inspect read retention, FastQC, alignment, gene assignment, count distributions and sample correlations. Flag fewer than one million retained reads or less than 50% unique alignment for investigation; these flags are not automatic exclusions. There is no outcome-based outlier removal. Primary model is RIF versus control, with age adjustment as a documented sensitivity. Count filtering uses edgeR filterByExpr and normalization uses TMM, followed by voom/limma with robust empirical Bayes. Map version-stripped GENCODE Ensembl identifiers to unambiguous approved HGNC symbols and sum counts mapping to the same approved gene.

This cohort is excluded from pooled discovery, gene selection, threshold selection and subtype fitting. Test the released MetaRIF model, fixed biological programs and the corrected discovery signature without flipping poor AUCs or adapting thresholds. Retain normal and uncertain predictions in all denominators. Report feature coverage and all 22 sample QC records. Apparent complete separation in small samples does not justify a zero-width confidence interval.

The article/deposit protocol and collection-date discrepancies are unresolved source annotations. Report them with the raw-data processing provenance; do not claim that the deposited reads reproduce the article's stated library preparation.
