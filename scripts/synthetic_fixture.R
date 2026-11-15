# Fictional inputs for execution checks. No cohort data or fitted parameters are used.
make_synthetic_inputs <- function(directory, n = 300L, seed = 901L) {
  set.seed(seed)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  coin <- function(probability = 0.4) rbinom(n, 1L, probability)
  height <- runif(n, 1.5, 1.95)
  bmi <- runif(n, 19, 33)
  months <- sample(c(12L, 24L, 36L), n, replace = TRUE, prob = c(.15, .15, .7))
  baseline <- data.frame(
    SUBJECT_ID = sprintf("synthetic_%04d", seq_len(n)),
    AGE = runif(n, 20, 75), SEX = sample(1:2, n, replace = TRUE),
    HEIGTH = height, WEIGTH = bmi * height^2, BMI = bmi,
    DIAGNOSIS = sample(1:3, n, replace = TRUE),
    SIDE = sample(1:3, n, replace = TRUE),
    T0_SYMPT_TREATMENT = sample(1:3, n, replace = TRUE),
    Psycopathological = coin(), ANTIBODY = sample(1:3, n, replace = TRUE),
    INT_T0 = runif(n, 3, 10), GGCEF_T0 = runif(n, 16, 28),
    Suspension = as.integer(months < 36L), MONTHS_OF_TREAT = months,
    TREATMENT_DISC = ifelse(months < 36L, sample(1:3, n, replace = TRUE), 0L)
  )
  for (name in c("Aura", "Hypertension", "Sleep_Disorders", "PREV_T0", "DETOXPRE",
                 "PULSATING", "PAIN_MOVMENT", "FAMILIARITY", "Bbloc", "Caant",
                 "Tricyclic", "Antiepil", "SSRISNRI", "Antiipnt", "Botulin", "Pizotifene")) {
    baseline[[name]] <- coin()
  }
  longitudinal <- expand.grid(SUBJECT_ID = baseline$SUBJECT_ID,
                              MONTH = c(0L, 3L, 6L, 9L, 12L), CYCLE = 1:3)
  index <- match(longitudinal$SUBJECT_ID, baseline$SUBJECT_ID)
  random_intercept <- rnorm(n, 0, 3.5)
  time <- (longitudinal$CYCLE - 1L) * 15L + longitudinal$MONTH
  rows <- nrow(longitudinal)
  longitudinal$MMDs <- pmax(0, pmin(30, baseline$GGCEF_T0[index] * .6 -
                                   .06 * time + random_intercept[index] + rnorm(rows, 0, 3)))
  for (name in c("DOSE", "DaysWithCefalea", "GGFAR", "INT", "HADSA", "HADSD", "HIT6", "MIDAS")) {
    longitudinal[[name]] <- runif(rows, 1, 20)
  }
  longitudinal$RED_MMD_VST01 <- 100 * (baseline$GGCEF_T0[index] - longitudinal$MMDs) /
    baseline$GGCEF_T0[index]
  # Exercise imputation while leaving the outcome and endpoint rows observable.
  baseline$AGE[sample.int(n, 25L)] <- NA_real_
  missing_bmi <- sample.int(n, 25L)
  baseline[missing_bmi, c("BMI", "HEIGTH", "WEIGTH")] <- NA_real_
  longitudinal$INT[sample.int(rows, 60L)] <- NA_real_
  write.table(baseline, file.path(directory, "MigraineBaselineVars.csv"),
              sep = ";", row.names = FALSE)
  write.table(longitudinal, file.path(directory, "MigraineLongitudianlVars.csv"),
              sep = ";", row.names = FALSE)
  invisible(directory)
}
