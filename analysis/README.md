# Analysis scripts

`as_run/` holds the scripts that produced the reported results, grouped by folder. They are archived for reading. Machine-specific paths and a few identifier names were replaced with placeholders (`<workspace>`, `<original_archive>`, `<cluster_scratch>`, `<results_package>` and similar), so none of them runs unchanged.

A script that reads `scripts/<name>` refers to the copy of that file in the analysis working folder: `scripts/fit_studies.R` is `baseline/fit_studies.R`, `scripts/validate_signatures.R` is `baseline/validate_signatures.R`, and `scripts/rebuild_gse58144.R` is `subtypes/rebuild_gse58144.R`. Paths under `inputs/frozen/` refer to the frozen results package and code bundle (`Supplementary_tables/` and `Analysis_code/as_run/`). The results package is not included. `verify_pooled_results.py` is a separate check on the pooled results (see the README).

The scripts are listed in the order of the analysis. Folder `baseline/` is the discovery analysis, `subtypes/` adds the corrected GSE58144 fit, the newer external cohorts, MetaRIF and the donor-lineage analysis, and `gse207362/` is the second external cohort.

## Inputs and preparation

| Script | What it does |
|---|---|
| `preparation/prepare_inputs.py` | Builds the HGNC identifier lookup and the discovery sample metadata files. |
| `preparation/prepare_external.py` | Builds the 40 GSE243550 count columns. |
| `preparation/prepare_gse248784.py` | Quality control and donor pseudobulk of the 21 GSE248784 10x matrices. |
| `preparation/recount_prjna.py` | Recounts the six PRJNA576073 BAMs with featureCounts (unstranded and reverse-strand). |
| `preparation/extract_cao_composition.py` | Extracts donor cell-type fractions from the Cao source workbook. |
| `subtypes/fetch_gse58144_audit.py`, `subtypes/download_gse58144_raw.py` | Download the GSE58144 GEO records and the 230 raw arrays. |
| `subtypes/prepare_gse58144_rebuild.py` | Maps the 230 dye-swap arrays to 115 donors using the recorded dye key. |
| `subtypes/compare_gse58144_sources.R` | Compares the original processed matrix with the GEO deposited matrix (orientation check, see `docs/sources/`). |
| `subtypes/fetch_metarif.py` | Retrieves the MetaRIF files at the pinned commit. |
| `subtypes/nibi_xu_pipeline.py` | Cutadapt, FastQC, STAR and counting for the Xu raw reads, run under Slurm. |
| `subtypes/audit_new_validation.py` | Extracts the GSE307510 native count columns and checks the Xu supplementary table against the recounted libraries for sample order. |

## Study models and pooling

| Script | What it does |
|---|---|
| `baseline/fit_studies.R` | Defines `harmonize()`, `read_matrix()` and `fit_one()` (lines 11 to 61) and holds the study loop (from line 62). The file stops unless `LOAD_FUNCTIONS_ONLY=1` is set, which also empties the study list, so the superseded GSE58144 matrix cannot be reused. The 12 contrasts other than GSE58144 were fitted by running the loop with the study list set to those 12 studies. Study names can also be passed as arguments. |
| `subtypes/rebuild_gse58144.R` | Rebuilds GSE58144 from raw ratios and fits the 113-donor primary and 115-donor sensitivity models. |
| `baseline/transform_diagnostic.R` | Compares the GSE71331 fit with and without an extra log2 transform. |
| `baseline/composition_analysis.R` | Timing and tissue-composition sensitivities. |
| `baseline/fit_external.R` | Fits GSE248784 and GSE243550 and the quality-control sensitivities. |
| `baseline/pool_studies.R` | Pools signed study statistics with fixed effective-sample-size weights. |
| `baseline/final_sensitivities.R` | RNA quality weights, the earlier two-stage weighting and the modified Hartung-Knapp sensitivity. |
| `baseline/audit_summary.py` | Gene coverage, study agreement and influence summaries. |

## Pathways and programs

| Script | What it does |
|---|---|
| `baseline/pathways_and_programs.R` | FGSEA and CAMERA tests of the 57 fixed programs. |
| `baseline/program_replication.R` | FRY pooling and directional program replication in external cohorts. |
| `baseline/timing_programs.R` | Seven exploratory program slopes against LH day among ordinary controls. |

## Scores and external cohorts

| Script | What it does |
|---|---|
| `baseline/validate_signatures.R` | Leave-one-study-out and frozen 50-gene scores, published comparator panels, gene replication. |
| `subtypes/fit_lee_validation.R` | Fits the five GSE307510 (Lee 2025) native biopsies and the synechia sensitivity. |
| `subtypes/fit_xu_validation.R` | Fits the Xu recounts. |
| `subtypes/evaluate_new_cohorts.R` | Frozen scores, genes and programs in GSE307510 and Xu. |
| `gse207362/01_qc_geo_counts.R` | Checks the published 12 + 26 column mapping and library quality. |
| `gse207362/fit_gse207362_validation.R`, `gse207362/fit_studies.R` | External-only fit with the discovery harmonization and model. |
| `gse207362/evaluate_gse207362_cohort.R`, `gse207362/validate_signatures.R` | Frozen panel, gene and program replication. |
| `gse207362/evaluate_metarif_gse207362.R` | Released MetaRIF predictions and the cohort assignment odds ratio. |

## MetaRIF and subtypes

| Script | What it does |
|---|---|
| `subtypes/evaluate_metarif.R` | Applies the released models without refitting; coverage and helper-equivalence checks. |
| `subtypes/analyze_heterogeneity.R` | Fixed immune and metabolic axes, distributions and gene heterogeneity. |
| `subtypes/state_specificity.R` | Clinical enrichment and subtype-selected program comparisons. |
| `subtypes/clinical_specificity_sensitivities.R` | Cohort-omission sensitivity of the clinical enrichment. |

## Donor-lineage analysis (GSE248784)

| Script | What it does |
|---|---|
| `subtypes/annotate_external_cells.py` | Outcome-blind marker annotation and donor-by-lineage aggregation. |
| `subtypes/analyze_celltypes.R` | Primary and strict lineage tests with defined FDR families. |
| `subtypes/celltype_shared_reference.R` | The 17-donor shared-reference sensitivity. |

## Historical preprocessing

`as_run/historical_preprocessing/` holds the original analysis notebooks behind the archived input matrices:

- microarray series: GSE103465, GSE26787, GSE92324, and one notebook for GSE111974 with GSE188409 (Agilent arrays);
- RNA-seq: the collaborator (Bui) data;
- single-cell: GSE183837, PRJNA1141235 (the notebook file name has an extra digit), DoubletFinder doublet removal, and the merged donor aggregation.

The GSE58144 and GSE71331 notebooks are the original processing. They are superseded by `subtypes/rebuild_gse58144.R` and `baseline/fit_studies.R`.
