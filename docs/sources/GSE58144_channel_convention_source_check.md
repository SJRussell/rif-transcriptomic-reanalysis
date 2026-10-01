# GSE58144 two-colour channel convention source check

**Checked:** 2026-09-08  
**Scope:** source and file-format provenance only; no phenotype, differential-expression, or model-performance evidence was used to choose the sign.

## Decision

For the deposited GSE58144 raw files, reconstruct the red and green intensities as **Ch1 = Cy5/red** and **Ch2 = Cy3/green**, then reverse the two intensities on arrays where the study dye key says **Cy3 = sample**. This makes the resulting limma-style `M = log2(R/G)` equal to `log2(sample/reference)` for both members of every dye-swap pair.

This conclusion is **strongly supported for GSE58144 but is not a universal ImaGene rule**. The raw text files contain no scanner wavelength, dye name, or image filename in their header; their only relevant labels are `Ch1 Signal Mean` and `Ch2 Signal Mean`. The study repositories map each raw array to sample/reference dyes, while the deposited processed values define the intended biological direction. A BioDiscovery-authored ImaGene manual is consistent with red-first/green-second ordering, but also shows why channel numbers should not be interpreted from software name alone: ImaGene treats channels as loaded images, assigns the first two display colors red then green by default, and permits users to change those colors without changing the quantified values.

## Primary-source chain

1. **Study design and intended processed direction.** The [GSE58144 GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE58144) says there are two dye-swapped profiles per sample and describes each deposited sample-table `VALUE` as “Quantile normalized and batch corrected log2 ratio (sample/reference).” The [primary publication](https://pmc.ncbi.nlm.nih.gov/articles/PMC4726345/) likewise says that each RNA sample produced two dye-swapped profiles, both compared with Universal Human Reference RNA, and that the two fold-change profiles were averaged after normalization. Neither source explicitly equates an ImaGene numeric channel with a dye.

2. **Raw-file-to-dye mapping.** The [E-MTAB-2591 SDRF](https://www.ebi.ac.uk/biostudies/files/E-MTAB-2591/E-MTAB-2591.sdrf.txt) assigns a dye and biomaterial to every raw filename. The same mappings are represented in GEO's `GSE58144_rawfiles_dye_orientation.txt.gz`. For example:

   | donor | raw file | sample dye | reference dye |
   |---|---|---|---|
   | 10C | `13223_raw.txt` | Cy5 | Cy3 |
   | 10C | `13232_raw.txt` | Cy3 | Cy5 |
   | 6P | `13088_raw.txt` | Cy3 | Cy5 |
   | 69C | `13108_raw.txt` | Cy5 | Cy3 |

   This establishes the biological contents of Cy3 and Cy5 on each array, but not the numeric `Ch1`/`Ch2` labels inside the raw file.

   GEO's per-GSM `label_ch1 = Cy3,Cy5` and `label_ch2 = Cy5,Cy3` fields are also **not** a raw-column key. In GEO, channel 1 is the endometrial sample and channel 2 is Universal Human Reference RNA; because each GSM combines the two physical dye-swap arrays, each biological channel lists both dyes in opposite order.

3. **ImaGene documentation.** The BioDiscovery-authored [ImaGene User Manual v5.5](https://manualzilla.com/doc/5685471/imagene-user-guide) states that loaded images are listed in the Images panel and that the first two images default to red and green, although the display colors can be changed without affecting quantified values. It also names the within-program ratio `Signal Mean Ratio ch1/ch2` and says that ImaGene saves quantified data from each image and leaves ratio generation to downstream software. A later [BioDiscovery ImaGene guide for Agilent arrays](https://cdn2.hubspot.net/hubfs/644585/Bio_-_InstallationActivation_Guides/imagene_guidance_brochure.pdf) similarly treats the two channels as separate images/files and does not prescribe Cy3 as channel 1. The v5.5 UI description is consistent with Ch1=red and Ch2=green, but it is older than the study's ImaGene 8.0 workflow and does not explicitly define the columns in this deposit.

4. **Same-laboratory analysis convention.** The UMC Utrecht-authored [Bioconductor `dyebias` manual](https://bioconductor.org/packages/devel/bioc/manuals/dyebias/man/dyebias.pdf) defines `R = Cy5`, `G = Cy3`, and `M = log2(R/G) = log2(Cy5/Cy3)` (manual p. 9). It warns that dye-swaps must be handled according to the identities occupying Cy3 and Cy5. This fixes the meaning of `R`, `G`, and `M`, but does not independently map GSE58144's raw numeric `Ch1` and `Ch2` columns.

5. **Deposit-level sign check.** The orientation audit records that the original reconstruction using `R = Ch2`, `G = Ch1`, followed by a swap when Cy3 contained the sample, is almost an exact negative of GEO's deposited `sample/reference` matrix (median gene-wise correlation `-0.9917559`; 99.99471% negative correlations across 18,919 approved genes and 113 shared samples). See [GSE58144_ORIENTATION_AUDIT.md](GSE58144_ORIENTATION_AUDIT.md). This is a comparison to the submitter-deposited processing direction, not to RIF status or classifier performance.

## Algebra of the correction

With `Ch1 = Cy5` and `Ch2 = Cy3`:

- sample in Cy5: `log2(Ch1/Ch2) = log2(sample/reference)`;
- sample in Cy3: swap the two channel intensities first, then `log2(R/G) = log2(sample/reference)`.

The original mapping starts with the inverse (`R = Ch2`, `G = Ch1`). Its subsequent dye-orientation swap therefore leaves both array orientations as `log2(reference/sample)`, explaining the observed near-perfect sign reversal relative to GEO.

## Residual uncertainty and safe use

- No GSE58144 raw header or accessible P-UMCU-42 protocol line explicitly says “Ch1 is Cy5 and Ch2 is Cy3.” The conclusion depends on the BioDiscovery red-first ordering plus the source-deposited `sample/reference` direction and the observed inversion.
- Public documents sometimes use the opposite generic shorthand, `channel 1 = Cy3` and `channel 2 = Cy5`. Those conventions describe other export/import pipelines and should not override this deposit's own processed direction.
- The rebuild should retain a direct, case-blind assertion that its reconstructed donor matrix agrees positively with GEO's `sample/reference` values. This verifies orientation without selecting a sign from clinical outcomes.

## Local files inspected

- `GSE58144_family.soft.gz`
- `GSE58144_series_matrix.txt.gz`
- `GSE58144_rawfiles_dye_orientation.txt.gz`
- `raw/GSM1402321_13223_raw.txt.gz`
- `raw/GSM1402321_13232_raw.txt.gz`
- `GSE58144_analysis_mar12.Rmd`
