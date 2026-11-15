# Migraine treatment response

Which baseline variables are associated with treatment discontinuation or response, and how do monthly migraine days change across repeated treatment cycles? This R project addresses those questions with multiple imputation, pooled logistic regression and longitudinal mixed-effects models.

**The clinical data are confidential.** The public repository contains the analysis code and a reproducible execution check on fictional inputs. It contains no patient records, clinical results, patient plots or submitted reports.

## Analysis workflow

| Question | Approach |
| :--- | :--- |
| Which baseline variables are associated with discontinuation? | Logistic regression across 20 MICE imputations, followed by coefficient pooling. |
| Which variables are associated with a response of at least 50%? | Separate logistic models at month 12 of each treatment cycle. |
| How do monthly migraine days vary over time? | Random-intercept models comparing continuous time with cycle-by-month effects. |

Multiple imputation represents uncertainty about missing baseline values. Pooling the logistic models combines information across the completed datasets. The longitudinal models account for repeated observations within a person through a random intercept; they use the first completed baseline dataset.

The implementation exposes the analysis as a function that accepts an external data directory and returns fitted objects in memory. [Methods and limitations](docs/methods.md) describe the selection rules, model specifications and treatment of missing values.

## What can be reproduced publicly

The synthetic check generates fictional baseline and longitudinal tables, runs imputation and modelling, and verifies that the specified outcomes are excluded from the imputation predictor matrices. It then removes the temporary inputs.

This verifies execution and selected workflow properties. Synthetic outputs do not reproduce clinical findings or validate treatment effects.

## Contribution

**Matteo Pietro Di Pilato** carried out Task 2: discontinuation and responder models, and longitudinal mixed-effects analysis. The workflow builds on preprocessing and imputation developed by **Matteo De Vecchi** for Tasks 0–1. **Marco Sabbadin** carried out the separate phenotyping task, which is not included here.

Medical Application and Healthcare, 2025–2026 · BSc in Artificial Intelligence.

## Limits and next steps

The models describe associations in the study data; this public code does not establish causal treatment effects. The longitudinal analysis uses one completed baseline dataset rather than pooling across imputations, and the original selection rules affect interpretation.

A methodological follow-up would assess sensitivity to imputation and model specification on authorised data. Any resulting patient-level outputs or fitted objects must remain in the private analysis environment.

## Run the synthetic check

From the repository root, with R installed:

```sh
Rscript --vanilla scripts/install_dependencies.R
Rscript --vanilla scripts/check.R
```

No clinical dataset is needed or downloaded.

For authorised local data, keep both CSVs outside the repository:

```r
source("R/coursework.R")
analysis <- run_coursework("/path/to/authorised/private-data")
```

The inputs retain the original semicolon separator, column names and filename spelling `MigraineLongitudianlVars.csv`. Returned model and imputation objects contain input-record information and should remain private. The function writes no model files, tables or figures.
