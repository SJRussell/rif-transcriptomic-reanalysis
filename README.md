# RIF endometrial transcriptomic reanalysis

Methods and analysis code for a systematic review, study-level meta-analysis and external evaluation of endometrial gene expression in recurrent implantation failure (RIF). The repository holds what is needed to judge the methods. It holds no raw data, expression matrices or manuscript text.

## Manuscript (under review)

*A systematic review and meta-analysis of endometrial transcriptomic reproducibility and subtype specificity in recurrent implantation failure.*

Authors: C. Keshen, C. Duval, G. Steba, S. Russell (corresponding author), C. Librach. S. Russell and C. Librach are co-senior authors.

## Question

Do endometrial expression differences between RIF and control biopsies reproduce across studies, and do published immune and metabolic RIF subtypes separate RIF from controls? Study influence, platform disagreement, biopsy timing and tissue composition are part of that question. An expression score or a subtype call is not evidence that a biopsy test improves patient outcomes.

## Design

- **Discovery:** 13 contrasts (microarray, bulk RNA-seq and single-cell donor pseudobulk), 299 samples (142 RIF, 157 controls).
- **Study models:** each study is fitted separately, RIF against control. Arrays use limma. Counts use edgeR `filterByExpr`, TMM and limma-voom. Both use robust empirical Bayes moderation.
- **Pooling:** signed Z statistics are combined with weights equal to the square root of each study's effective sample size. The primary family is 16,008 genes tested in at least two array and two RNA-seq studies. 1,840 genes pass BH FDR 0.05.
- **Sensitivity:** leave-one-study-out, platform-specific, timing, tissue-composition and weighting analyses.
- **Programs:** FGSEA pathways and 57 fixed programs (FRY and CAMERA).
- **Scores and subtypes:** a 50-gene score (selected leave-one-study-out within discovery, then frozen for external cohorts), two published comparator panels, and the released MetaRIF models applied without refitting. External cohorts are GSE248784 and GSE207362, with technical comparisons in GSE243550, GSE307510 and PRJNA1335587.
- **Seed:** 20260908 for permutations and bootstraps.

## Where to look

| Path | Contents |
|---|---|
| `docs/methods.md` | Full methods, section by section. |
| `docs/cohorts.md` | Each cohort's role and sample decisions, and the datasets considered and not used. |
| `docs/sources/` | GSE58144 channel-orientation audit and the sources for the fixed program membership. |
| `docs/source_checks/` | Sample mappings for GSE207362, GSE307510 and Xu; the pinned MetaRIF record. |
| `analysis/README.md` | Script index in analysis order. |
| `analysis/as_run/` | The scripts that produced the results. |
| `analysis/verify_pooled_results.py` | Recomputes pooled statistics from per-study outputs. |
| `provenance/` | SHA-256 of the 45 non-array inputs (`inputs.csv`) and the recount manifests. The 242 raw GSE58144 and GSE71331 array files are identified by accession. |
| `environments/` | R package versions and `sessionInfo()` output. |
| `results/reference/` | Aggregate reference tables and headline numbers (`expected.json`). |

## Using the code

This is a code repository, not a run package. The scripts are archived as run. Machine-specific paths and a few identifier names were replaced with placeholders, so none of the scripts runs unchanged. Inputs come from public downloads, restricted data and archived matrices, and a reanalysis has to supply its own.

`analysis/verify_pooled_results.py` recomputes signed Z, pooled statistics, BH adjustment, program statistics and odds ratios from the per-study outputs. It needs the full results package, which is not included. It passes against the package behind the reference tables: 13 studies, 299 samples, 16,008 genes, 1,840 DEGs, 57 programs. It does not refit any model.

## Limits of the record

- No end-to-end rerun from raw data was done. GSE58144 and GSE71331 are rebuilt from raw arrays in the scripts. The other five microarray series and the three single-cell discovery cohorts enter as archived processed matrices. Their notebooks are in `analysis/as_run/historical_preprocessing/`, but the full path from raw files to matrix was not replayed.
- PRJNA576073 counts come from recounting existing BAM files. The upstream alignment settings cannot be recovered from the recount.
- The Bui counts and metadata are restricted and are not in the repository. The Bui and GSE205398 counts enter as archived count matrices, and GSE205398 has no preprocessing notebook.
- The release of the MSigDB gene sets is not recorded here. `provenance/inputs.csv` gives the hash of the file used.
- The full database search strategies and some record-level screening details are not documented here.
- Control groups follow each source study's definition (see `docs/cohorts.md`).
