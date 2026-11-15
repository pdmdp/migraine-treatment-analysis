# Migraine treatment response

R coursework on treatment discontinuation, response at the end of three treatment cycles, and repeated measurements of monthly migraine days. The analysis combines multiple imputation, pooled logistic regression and linear mixed-effects models.

The clinical dataset is confidential. This repository contains the analysis code and an execution check using entirely fictional inputs. It includes no clinical records, fitted results, patient plots or copies of the submitted reports.

## Analysis

| Question | Approach |
| :--- | :--- |
| Which baseline variables are associated with discontinuation? | Logistic regression across 20 MICE imputations, with pooled coefficients. |
| Which variables are associated with a response of at least 50%? | Separate logistic models at month 12 of each treatment cycle. |
| How do monthly migraine days vary over time? | Random-intercept models comparing continuous time with cycle-by-month effects. |

The longitudinal models use the first completed baseline dataset. The original workflow's selection rules and model specifications are retained; [methods and limitations](docs/methods.md) explain how these choices affect interpretation.

## Run the synthetic check

From the repository root, with R installed:

```sh
Rscript --vanilla scripts/install_dependencies.R
Rscript --vanilla scripts/check.R
```

The check generates fictional baseline and longitudinal tables in a temporary directory, runs the imputation and modelling workflow, verifies the outcome exclusions in the imputation matrices, and removes the temporary files. It does not download data or reproduce the clinical findings.

## Use authorised local data

Keep the two input CSVs in a directory outside this repository. The analysis retains the source column names and semicolon separator, including the original filename spelling `MigraineLongitudianlVars.csv`.

```r
source("R/coursework.R")
analysis <- run_coursework("/path/to/authorised/private-data")
```

The function returns the fitted models and imputation objects in memory. Those objects contain information from the input records and should remain private. It writes no model files, tables or figures.

## Source

The public script was extracted from `Final.Rmd`, the recovered coursework file whose random-intercept specification agrees with the final report. The November 2026 preparation removed exploratory prose, standalone displays, summaries and plotting calls, added an external input directory and packaged the code as a function. The 180 retained computation expressions were compared with the source; only the two input-path expressions were adapted.

The original files remain separate. [Publication notes](docs/publication.json) record the source checksum and verification scope.

## Contribution

The group report, *Migraine Treatment Dynamics: A 3-Year Clinical Data Project*, assigns Task 2 to **Matteo Pietro Di Pilato**: discontinuation and responder models, and longitudinal mixed-effects analysis. The analysis builds on preprocessing and imputation developed by **Matteo De Vecchi** for Tasks 0–1. **Marco Sabbadin** carried out the study's separate phenotyping task, which is not included here.

Medical Application and Healthcare Course, 2025–2026 · BSc in Artificial Intelligence.
