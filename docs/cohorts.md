# Cohorts and sample decisions

Sample lists come from source metadata and fixed manifests. They are never chosen or changed by looking at expression. `results/reference/Main_tables/Table1_manuscript.csv` records the primary denominators.

## Inclusion rule

Candidate datasets were judged against a written rule: age 40 or under, assisted reproduction, a mid-secretory biopsy at LH+5 to LH+7, RIF after at least three failed embryo transfers, controls with regular menstruation and at least one live birth, and a list of uterine, endometrial, chromosomal, immune and metabolic exclusions. The protocol files give different age, transfer-count and timing thresholds, and no dataset was promoted by choosing the more permissive wording.

In the discovery set, each source study's own definition of the comparison group is kept. A live-birth requirement for controls would keep only 6 of the 13 discovery cohorts. GSE92324 uses healthy oocyte donors without a documented live birth, and GSE58144 reports live birth after IVF or ICSI.

Studies are placed in one of three roles. Discovery cohorts are pooled. External cohorts are scored with the frozen discovery panel and programs and never enter selection or weighting. Technical comparisons are repository-labelled data whose link to a clinical cohort is incomplete.

## Discovery cohorts (13 contrasts, 299 samples)

| Cohort | RIF/control | Input and fixed decision |
|---|---:|---|
| GSE111974 (Bastu 2019) | 24/24 | Archived microarray matrix. 5,324 within-study DEGs; its strong influence is reported and its weight is not changed. |
| GSE188409 (Zhao 2022) | 5/5 | Archived microarray matrix. |
| GSE103465 (Guo 2018) | 3/3 | Archived microarray matrix. |
| GSE92324 (Pathare 2017) | 10/8 | Archived microarray matrix. GSE71835 overlaps it; control_OD21 and control_OD22 have different timing in the two deposits. |
| GSE26787 (Lédée 2011) | 5/5 | Archived microarray matrix. |
| GSE58144 (Koot 2016) | 42/71 | Rebuilt from the 230 raw dye-swap arrays. The primary set drops GSM1402357 (RIF, LH+8) and GSM1402432 (control, LH+6), the two samples excluded in the original PCA; the original numerical threshold was not recorded. All 115 donors (43/72) are a sensitivity. GEO lists 115 donors, MetaRIF training used 110, the primary set has 113. |
| GSE71331 (Shi 2018) | 7/5 | Rebuilt from raw Agilent files with probe-name mapping. |
| Bui 2024, shared by the source investigators | 22/22 | Archived counts. Two samples lack LH-relative timing, so the primary set is 44 and the 46-sample set is an unadjusted sensitivity. The restricted source metadata confirms the two unknowns. |
| GSE205398 (Wang 2024) | 3/3 | Archived bulk counts. |
| PRJNA576073 (Zhou 2019) | 3/3 | Unstranded recount of the existing BAMs. SRA SRP224538. Zhang 2023 reused the same samples. |
| PRJNA1141235 (Fu 2024) | 2/2 | Mid-secretory normal-endometrium RIF only: SRR30016494 and SRR30016501, with controls SRR30016499 and SRR30016504. Thin-endometrium RIF and proliferative runs are excluded. |
| GSE250130 (Cao 2025) | 10/3 | The three ordinary LH+7 controls (con4 to con6). GEO has 28 biopsies and the archive 22 (con1 to con12, rif1 to rif10). The LH+3 trio is absent from the archive. The three PGT controls are not primary: two of the three donors had structural karyotype variants and the public supplement does not say which control is the normal-karyotype donor. |
| GSE183837 (Lai 2022) | 6/3 | Archived donor pseudobulk. |

GSE223672 has 3 RIF samples and no controls, so it cannot enter a study contrast.

## External and technical cohorts

| Cohort | RIF/control | Role and caveats |
|---|---:|---|
| GSE248784 | 11/10 | External single-cell cohort named in its publication. Donor labels follow the GEO titles. RIF enrollment allowed a two-attempt route for some participants, and controls reached clinical pregnancy rather than documented live birth. No public accession was found for the paper's separate bulk cohort (48 RIF, 19 controls). |
| GSE207362 (Fukui 2023) | 12/26 | External bulk cohort. Biopsies were taken in hormone replacement cycles on day 7 of progesterone (P+7); no discovery cohort was sampled under hormone replacement. The 12 + 26 samples are the paper's published subset of the 67 deposited libraries, mapped to count-column suffixes (`docs/source_checks/GSE207362/published_subset_mapping.csv`). The mapping assumes the paper's within-group numbers equal the column suffixes; GEO metadata cannot confirm it. Controls reached clinical pregnancy. It is not a MetaRIF development cohort. |
| GSE243550 | 20/20 | Technical comparison. GEO verifies 20 RIF and 20 control P+5 libraries. The related paper describes a 20 + 20 clinical cohort, reports RT-qPCR on 10 + 10, and does not name the accession, so the library-to-patient link is unresolved. Its own discovery analysis reuses GSE58144, GSE103465 and GSE111974, which are not recounted here. |
| GSE307510 (Lee 2025) | 3/2 | Technical comparison with five donors. Native biopsies only, all on post-ovulation day 5. The five matched chip libraries come from the same donors and are excluded. The source gives no operational RIF definition or control pregnancy outcome, and RIF-2 has synechia, so a sensitivity omits that donor. Sample manifest: `docs/source_checks/GSE307510/native_biopsy_sample_manifest.csv`. |
| PRJNA1335587 (Xu 2026) | 9/13 | Technical comparison from raw reads. Sample identity (CON1 to CON13, RIF1 to RIF9) comes from BioSample titles (`docs/source_checks/Xu2026/verified_sample_mapping.csv`). The paper and the deposit disagree on library protocol, instrument, read length, collection dates and ages, and the paper's supplementary table is not a usable count matrix. Processing rules are in `docs/source_checks/Xu2026/RAW_VALIDATION_RULES.md`. |

## Datasets considered and not used

| Dataset | Reason |
|---|---|
| GSE71835 | Nine of its 12 samples are the same as GSE92324. Counted once. |
| GSE106602 (PRJNA417762) | 19 `IVF_fri` and 16 healthy samples at LH+7. The linked paper does not define the RIF cohort, the comparison group or the exclusions, and no exclusion reason was recorded. It is also a MetaRIF training cohort, so it could not test that model independently. |
| GSE108966 | miRNA-seq, outside the mRNA platforms analysed. |
| GSE252280 and GSE287072 | RIF-only. The two accessions share 29 exact duplicate profiles, 19 of them with different declared days, and the public labels give no reliable donor or day key. |
| Mo et al. 2 + 2 single-cell cohort | No public accession; available on request. The 2026 report appears to describe the same small cohort as the 2024 report. |
| GSE305811 (PRJNA1307877) | Repository labels show 25 RIF and 3 controls inside a larger neoplasia series. No linked paper or comparison-group definition. |
| GSE287278, GSE284236 | Spatial transcriptomics. |
| GSE297873, GSE295682, GSE279514, GSE269408 | Cell lines, organoids or cultured treatments, with no uncultured patient biopsy contrast. |
| PRJNA963093 | scATAC-seq, not RNA. |
| PRJNA1051672 | Duplicate of GSE250130. |
| GSE311010 | Private until 1 December 2026 at the 8 September 2026 check. No sample list was inferred. |
| HRA010876 | Not released at the 8 September 2026 check. It reuses four series already considered here. |
