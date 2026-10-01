# Software environments

The main analysis used R 4.6.1 and Bioconductor 3.23. Direct R requirements: limma 3.68.4, edgeR 4.10.1, fgsea 1.38.0, BiocParallel 1.46.0, metafor 5.0.1, mclust 6.1.3, diptest 0.77.2, caret 7.0-1, nnet 7.3-20, randomForest 4.7-1.2, jsonlite 2.0.0 and data.table 1.18.6.1, plus statmod for robust moderation.

| File | Contents |
|---|---|
| `R_packages.csv` | Every R package installed in the analysis environment, with version and dependencies. It is an observed inventory. |
| `observed/` | `sessionInfo()` output from the runs. The main analysis (R 4.6.1) is `celltypes_sessionInfo.txt` and `metarif_sessionInfo.txt`. The GSE207362 cluster runs used R 4.4.3 for the fit, scores and replication and R 4.5.0 for MetaRIF. |
| `requirements.txt` | Python versions for the preparation scripts and the verifier. |
