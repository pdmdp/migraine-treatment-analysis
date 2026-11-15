# Methods and scope

The group report credits Matteo De Vecchi with Tasks 0–1 (preprocessing and imputation), Matteo Pietro Di Pilato with Task 2 (response, discontinuation and longitudinal models), and Marco Sabbadin with Task 3 (phenotyping). The recovered Task 2 script includes the shared preparation steps needed to fit its models. The phenotyping analysis is outside this repository.

## Baseline and longitudinal preparation

The coursework reconciles height, weight and BMI where two measurements are available, recodes categorical variables and counts recorded preventive-drug failures. Missing entries in a partially observed drug history contribute zero to that count; a completely missing drug history remains missing. This assumption can undercount unrecorded failures.

Baseline missingness is handled with MICE: 20 completed datasets, 30 iterations and seed 123. Discontinuation, treatment duration and reason for discontinuation are excluded as imputation predictors. Baseline migraine frequency is also excluded as a predictor and is not imputed. These choices follow the recovered script.

Longitudinal records are filtered using active treatment months, `(CYCLE - 1) * 12 + MONTH`. The script compares last observation carried forward, next observation carried backward and interpolation. For Task 2, it carries selected covariates forward within each subject and cycle; it does not fill the MMD outcome. These imputation methods impose assumptions about unobserved changes.

## Discontinuation and responder models

Discontinuation is modelled from baseline covariates with logistic regression. Coefficients are pooled across imputations. The script also computes ROC curves on the same observations used to fit each model. These are apparent, in-sample discrimination estimates; they are not held-out or external validation.

A responder has at least a 50% reduction in MMDs relative to baseline. The endpoint tables are joined across all three cycles before fitting separate responder models. This restricts eligibility to subjects with endpoint rows in every cycle, with further omissions when an outcome is missing. It can favour subjects who remain in follow-up and should not be described as a deployment-time prediction study.

The recovered script treats antibody codes numerically and uses an ordered diagnosis factor. The public version preserves those encodings; interpretation depends on the original codebook. No clinical performance numbers are republished here.

## Mixed-effects models

The longitudinal analysis uses only the first completed baseline dataset. It therefore does not propagate between-imputation uncertainty into the mixed-model estimates.

The continuous model uses time and baseline migraine frequency as fixed effects, with a subject-specific random intercept. The categorical model uses the interaction of cycle and within-cycle month, with baseline frequency and the same random-intercept structure. Both are fitted by maximum likelihood (`REML = FALSE`). Calendar time adds a three-month gap between cycles, unlike the active-treatment scale used in the filtering step.

The final report describes random slopes as considered but not retained. `Final.Rmd` matches the reported random-intercept choice; the separate `Final.Rmd.py` copy contains a random-slope experiment. The publication follows `Final.Rmd`.

Associations and fitted trajectories from this observational coursework do not establish treatment effects or clinical usefulness. No clinical data were used in the publication checks.

## Publication check

The executable check uses fictional values generated independently of the clinical cohort and of the report's fitted estimates. It checks input guards, completion of the computation, imputation exclusions, the number of fitted logistic models and the mixed-model random-effects structure. It does not establish that the original cohort's results can be reproduced, assess clinical validity, or test every possible missing-data pattern.
