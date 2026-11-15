# Coursework analysis extracted from Final.Rmd (February 2026).
# Data displays and plots are omitted; model specifications are retained.
# See docs/methods.md for the scope of the original analysis.
run_coursework <- function(data_dir) {
  if (missing(data_dir) || !dir.exists(data_dir)) {
    stop('Provide a directory containing the two authorised input CSVs.')
  }
  data_dir <- normalizePath(data_dir, mustWork = TRUE)
  expected <- c('MigraineBaselineVars.csv', 'MigraineLongitudianlVars.csv')
  if (!all(file.exists(file.path(data_dir, expected)))) {
    stop('The input directory must contain both required CSV files.')
  }
  library(dplyr)

  library(mice)

  library(zoo)

  library(ggplot2)


  # Baseline preparation
  df <- read.csv(file.path(data_dir, "MigraineBaselineVars.csv"), sep = ";", stringsAsFactors = FALSE)

  df <- df %>% mutate(HEIGTH = ifelse(HEIGTH > 3, HEIGTH/100, HEIGTH)) %>% mutate(BMI_calc = WEIGTH/(HEIGTH^2), 
      BMI_final = ifelse(is.na(BMI), round(BMI_calc, 2), BMI))

  df <- df %>% mutate(SEX = case_when(SEX == 1 ~ 0, SEX == 2 ~ 1, TRUE ~ NA_real_))

  df <- df %>% mutate(Psycopathological = case_when(Psycopathological == 0 ~ 0, Psycopathological > 
      0 ~ 1, TRUE ~ NA_real_))

  drug_cols <- c("Bbloc", "Caant", "Tricyclic", "Antiepil", "SSRISNRI", "Antiipnt", "Botulin", 
      "Pizotifene")

  calc_failures <- function(x) {
      if (all(is.na(x))) {
          return(NA)
      }
      else {
          return(sum(x, na.rm = TRUE))
      }
  }

  df$Num_Failed_Drugs <- apply(df[, drug_cols], 1, calc_failures)

  df <- df %>% select(-all_of(drug_cols))

  df <- df %>% mutate(BMI = case_when(is.na(BMI) & !is.na(WEIGTH) & !is.na(HEIGTH) ~ round(WEIGTH/(HEIGTH^2), 
      2), TRUE ~ BMI), WEIGTH = case_when(is.na(WEIGTH) & !is.na(BMI) & !is.na(HEIGTH) ~ round(BMI * 
      (HEIGTH^2), 1), TRUE ~ WEIGTH), HEIGTH = case_when(is.na(HEIGTH) & !is.na(WEIGTH) & !is.na(BMI) ~ 
      round(sqrt(WEIGTH/BMI), 2), TRUE ~ HEIGTH))

  df <- df %>% mutate(DIAGNOSIS_Severity = case_when(DIAGNOSIS == 3 ~ 1, DIAGNOSIS == 1 ~ 2, 
      DIAGNOSIS == 2 ~ 3, TRUE ~ NA_real_)) %>% mutate(DIAGNOSIS = DIAGNOSIS_Severity) %>% 
      select(-DIAGNOSIS_Severity)

  binary_vars <- c("SEX", "Aura", "Hypertension", "Sleep_Disorders", "PREV_T0", "DETOXPRE", 
      "PULSATING", "PAIN_MOVMENT", "FAMILIARITY")

  ordinal_vars <- c("DIAGNOSIS", "SIDE", "T0_SYMPT_TREATMENT")

  df[binary_vars] <- lapply(df[binary_vars], factor)

  df[ordinal_vars] <- lapply(df[ordinal_vars], ordered)

  df$na_row <- rowMeans(is.na(df))

  vars_future <- c("Suspension", "MONTHS_OF_TREAT", "TREATMENT_DISC")

  vars_mice <- c("AGE", "SEX", "DIAGNOSIS", "BMI", "Num_Failed_Drugs", "FAMILIARITY", "Hypertension", 
      "Sleep_Disorders", "PREV_T0", "DETOXPRE", "Psycopathological", "T0_SYMPT_TREATMENT", 
      "ANTIBODY", "INT_T0", "GGCEF_T0", vars_future)

  df_mice <- df[, vars_mice]

  pred <- make.predictorMatrix(df_mice)

  diag(pred) <- 0

  ini <- mice(df_mice, maxit = 0, print = FALSE)

  meth <- ini$method

  pred <- ini$predictorMatrix

  meth[vars_future] <- ""

  pred[, vars_future] <- 0

  meth["GGCEF_T0"] <- ""

  pred[, "GGCEF_T0"] <- 0

  meth[c("AGE", "BMI", "INT_T0", "Num_Failed_Drugs")] <- "pmm"

  meth[c("SEX", "Hypertension", "Sleep_Disorders", "PREV_T0", "DETOXPRE", "FAMILIARITY")] <- "logreg"

  meth[c("DIAGNOSIS", "T0_SYMPT_TREATMENT")] <- "polr"

  meth[c("ANTIBODY")] <- ""

  set.seed(123)


  # Baseline multiple imputation
  imp <- mice(df_mice, m = 20, maxit = 30, method = meth, predictorMatrix = pred, printFlag = TRUE)

  df_view <- complete(imp, 1)

  fit <- with(imp, lm(GGCEF_T0 ~ AGE + SEX + BMI + DIAGNOSIS + Psycopathological + Num_Failed_Drugs))

  pooled <- pool(fit)


  # Longitudinal preparation
  longit_data <- read.csv(file.path(data_dir, "MigraineLongitudianlVars.csv"), sep = ";", stringsAsFactors = FALSE)

  longit_joined <- longit_data %>% left_join(df %>% select(SUBJECT_ID, MONTHS_OF_TREAT), by = "SUBJECT_ID")

  longit_clean <- longit_joined %>% mutate(abs_month = (CYCLE - 1) * 12 + MONTH) %>% filter(abs_month <= 
      MONTHS_OF_TREAT)

  longit_clean <- longit_joined %>% mutate(abs_month = (CYCLE - 1) * 12 + MONTH) %>% filter(abs_month <= 
      MONTHS_OF_TREAT)

  longit_clean <- dplyr::arrange(longit_clean, CYCLE, MONTH)

  longit_clean <- longit_clean %>% select(-MONTHS_OF_TREAT, -abs_month)

  vars_imp <- c("MMDs", "DOSE", "DaysWithCefalea", "GGFAR", "INT", "HIT6", "MIDAS")

  vars_imp <- intersect(vars_imp, names(longit_clean))

  vars_no_imp <- c("HADSA", "HADSD", "RED_MMD_VST01")

  long_obs <- longit_clean

  long_locf <- longit_clean %>% arrange(SUBJECT_ID, CYCLE, MONTH) %>% group_by(SUBJECT_ID, 
      CYCLE) %>% mutate(across(all_of(vars_imp), ~na.locf(., na.rm = FALSE))) %>% ungroup()

  long_nocb <- longit_clean %>% arrange(SUBJECT_ID, CYCLE, MONTH) %>% group_by(SUBJECT_ID, 
      CYCLE) %>% mutate(across(all_of(vars_imp), ~na.locf(., fromLast = TRUE, na.rm = FALSE))) %>% 
      ungroup()

  long_interp <- longit_clean %>% arrange(SUBJECT_ID, CYCLE, MONTH) %>% group_by(SUBJECT_ID, 
      CYCLE) %>% mutate(across(all_of(vars_imp), ~na.approx(., x = MONTH, na.rm = FALSE))) %>% 
      ungroup()

  na_count <- function(dat, method_name) {
      tibble(Method = method_name, Variable = vars_imp, N_missing = sapply(vars_imp, function(v) sum(is.na(dat[[v]]))))
  }

  na_report <- bind_rows(na_count(long_obs, "Observed"), na_count(long_locf, "LOCF"), na_count(long_nocb, 
      "NOCB"), na_count(long_interp, "Interpolation"))

  summ_traj <- function(dat, method) {
      dat %>% group_by(CYCLE, MONTH) %>% summarise(mean_MMDs = mean(MMDs, na.rm = TRUE), mean_INT = mean(INT, 
          na.rm = TRUE), .groups = "drop") %>% mutate(Method = method)
  }

  traj_comp <- bind_rows(summ_traj(long_obs, "Observed"), summ_traj(long_locf, "LOCF"), summ_traj(long_nocb, 
      "NOCB"), summ_traj(long_interp, "Interpolation"))

  get_slopes <- function(dat, method) {
      dat %>% group_by(SUBJECT_ID, CYCLE) %>% summarise(slope_mmd = ifelse(sum(!is.na(MMDs)) >= 
          2, coef(lm(MMDs ~ MONTH))[2], NA_real_), .groups = "drop") %>% mutate(Method = method)
  }

  slopes_comp <- bind_rows(get_slopes(long_obs, "Observed"), get_slopes(long_locf, "LOCF"), 
      get_slopes(long_nocb, "NOCB"), get_slopes(long_interp, "Interpolation"))

  vars_imp_2 <- c("DOSE", "DaysWithCefalea", "GGFAR", "INT", "HIT6", "MIDAS")

  vars_imp_2 <- intersect(vars_imp_2, names(longit_clean))

  vars_no_imp <- c("HADSA", "HADSD", "RED_MMD_VST01", "MMDs")

  long_obs_t2 <- longit_clean

  long_locf_t2 <- longit_clean %>% arrange(SUBJECT_ID, CYCLE, MONTH) %>% group_by(SUBJECT_ID, 
      CYCLE) %>% mutate(across(all_of(vars_imp_2), ~na.locf(., na.rm = FALSE))) %>% ungroup()

  long_nocb_t2 <- longit_clean %>% arrange(SUBJECT_ID, CYCLE, MONTH) %>% group_by(SUBJECT_ID, 
      CYCLE) %>% mutate(across(all_of(vars_imp_2), ~na.locf(., fromLast = TRUE, na.rm = FALSE))) %>% 
      ungroup()

  long_interp_t2 <- longit_clean %>% arrange(SUBJECT_ID, CYCLE, MONTH) %>% group_by(SUBJECT_ID, 
      CYCLE) %>% mutate(across(all_of(vars_imp_2), ~na.approx(., x = MONTH, na.rm = FALSE))) %>% 
      ungroup()

  na_count <- function(dat, method_name) {
      tibble(Method = method_name, Variable = vars_imp_2, N_missing = sapply(vars_imp_2, function(v) sum(is.na(dat[[v]]))))
  }

  na_report <- bind_rows(na_count(long_obs_t2, "Observed"), na_count(long_locf_t2, "LOCF"), 
      na_count(long_nocb_t2, "NOCB"), na_count(long_interp_t2, "Interpolation"))

  baseline_list <- complete(imp, action = "all")

  baseline_list <- lapply(baseline_list, function(dat) {
      dat %>% mutate(SUBJECT_ID = df$SUBJECT_ID) %>% relocate(SUBJECT_ID)
  })

  long_task2_list <- lapply(baseline_list, function(baseline_i) {
      long_locf_t2 %>% left_join(baseline_i, by = "SUBJECT_ID")
  })

  vars_long <- c("CYCLE", "MONTH", "MMDs", "RED_MMD_VST01", "DOSE", "GGFAR", "DaysWithCefalea", 
      "INT", "HADSA", "HADSD", "HIT6", "MIDAS")


  # Discontinuation model
  model_disc_imp <- with(imp, glm(Suspension ~ AGE + SEX + DIAGNOSIS + BMI + ANTIBODY + Num_Failed_Drugs + 
      DETOXPRE + Psycopathological + Sleep_Disorders + PREV_T0 + GGCEF_T0 + INT_T0, family = binomial))

  pooled_disc <- pool(model_disc_imp)

  completed_list <- complete(imp, action = "all")

  library(pROC)

  roc_list <- vector("list", length(completed_list))

  auc_values <- numeric(length(completed_list))

  for (i in seq_along(completed_list)) {
      d <- completed_list[[i]]
      fit <- model_disc_imp$analyses[[i]]
      probs <- predict(fit, newdata = d, type = "response")
      roc_list[[i]] <- roc(d$Suspension, probs, quiet = TRUE)
      auc_values[i] <- auc(roc_list[[i]])
  }

  mean_auc <- mean(auc_values)

  sd_auc <- sd(auc_values)

  fpr_grid <- seq(0, 1, length.out = 100)

  tpr_mat <- sapply(roc_list, function(r) {
      approx(1 - r$specificities, r$sensitivities, xout = fpr_grid)$y
  })

  mean_tpr <- rowMeans(tpr_mat, na.rm = TRUE)

  library(car)

  vif_values <- lapply(model_disc_imp$analyses, vif)

  vif_mean <- Reduce("+", vif_values)/length(vif_values)

  long_locf_t2 <- long_locf_t2 %>% arrange(SUBJECT_ID, CYCLE, MONTH)


  # Response at the end of each cycle
  mmd_c1m12 <- long_locf_t2 %>% filter(CYCLE == 1, MONTH == 12) %>% select(SUBJECT_ID, MMDs) %>% 
      rename(mmd_c1m12 = MMDs)

  mmd_c2m12 <- long_locf_t2 %>% filter(CYCLE == 2, MONTH == 12) %>% select(SUBJECT_ID, MMDs) %>% 
      rename(mmd_c2m12 = MMDs)

  mmd_c3m12 <- long_locf_t2 %>% filter(CYCLE == 3, MONTH == 12) %>% select(SUBJECT_ID, MMDs) %>% 
      rename(mmd_c3m12 = MMDs)

  data_resp <- df %>% left_join(mmd_c1m12, by = "SUBJECT_ID") %>% filter(!is.na(mmd_c1m12))

  data_resp <- df %>% left_join(mmd_c2m12, by = "SUBJECT_ID") %>% filter(!is.na(mmd_c2m12))

  data_resp <- df %>% left_join(mmd_c3m12, by = "SUBJECT_ID") %>% filter(!is.na(mmd_c3m12))

  mmd_all <- mmd_c1m12 %>% inner_join(mmd_c2m12, by = "SUBJECT_ID") %>% inner_join(mmd_c3m12, 
      by = "SUBJECT_ID")

  data_resp <- df %>% inner_join(mmd_all, by = "SUBJECT_ID") %>% mutate(red_pct_c1m12 = 100 * 
      (GGCEF_T0 - mmd_c1m12)/GGCEF_T0, responder_c1m12 = as.integer(red_pct_c1m12 >= 50), red_pct_c2m12 = 100 * 
      (GGCEF_T0 - mmd_c2m12)/GGCEF_T0, responder_c2m12 = as.integer(red_pct_c2m12 >= 50), red_pct_c3m12 = 100 * 
      (GGCEF_T0 - mmd_c3m12)/GGCEF_T0, responder_c3m12 = as.integer(red_pct_c3m12 >= 50))

  vars_future <- c("Suspension", "MONTHS_OF_TREAT", "TREATMENT_DISC")

  vars_mice <- c("SUBJECT_ID", "red_pct_c1m12", "responder_c1m12", "red_pct_c2m12", "responder_c2m12", 
      "red_pct_c3m12", "responder_c3m12", "mmd_c1m12", "mmd_c2m12", "mmd_c3m12", "AGE", "SEX", 
      "DIAGNOSIS", "BMI", "Num_Failed_Drugs", "FAMILIARITY", "Hypertension", "Sleep_Disorders", 
      "PREV_T0", "DETOXPRE", "Psycopathological", "T0_SYMPT_TREATMENT", "ANTIBODY", "INT_T0", 
      "GGCEF_T0", vars_future)

  df_mice <- data_resp[, vars_mice]

  pred <- make.predictorMatrix(df_mice)

  diag(pred) <- 0

  ini <- mice(df_mice, maxit = 0, print = FALSE)

  meth <- ini$method

  pred <- ini$predictorMatrix

  meth[vars_future] <- ""

  pred[, vars_future] <- 0

  meth["GGCEF_T0"] <- ""

  meth["SUBJECT_ID"] <- ""

  pred[, "SUBJECT_ID"] <- 0

  pred[, "GGCEF_T0"] <- 0

  meth["red_pct_c1m12"] <- ""

  pred[, "red_pct_c1m12"] <- 0

  meth["responder_c1m12"] <- ""

  pred[, "responder_c1m12"] <- 0

  meth["mmd_c1m12"] <- ""

  pred[, "mmd_c1m12"] <- 0

  meth["red_pct_c2m12"] <- ""

  pred[, "red_pct_c2m12"] <- 0

  meth["responder_c2m12"] <- ""

  pred[, "responder_c2m12"] <- 0

  meth["mmd_c2m12"] <- ""

  pred[, "mmd_c2m12"] <- 0

  meth["red_pct_c3m12"] <- ""

  pred[, "red_pct_c3m12"] <- 0

  meth["responder_c3m12"] <- ""

  pred[, "responder_c3m12"] <- 0

  meth["mmd_c3m12"] <- ""

  pred[, "mmd_c3m12"] <- 0

  meth[c("AGE", "BMI", "INT_T0", "Num_Failed_Drugs")] <- "pmm"

  meth[c("SEX", "Hypertension", "Sleep_Disorders", "PREV_T0", "DETOXPRE", "FAMILIARITY")] <- "logreg"

  meth[c("DIAGNOSIS", "T0_SYMPT_TREATMENT")] <- "polr"

  meth[c("ANTIBODY")] <- ""

  set.seed(123)


  # Responder multiple imputation
  imp1 <- mice(df_mice, m = 20, maxit = 30, method = meth, predictorMatrix = pred, printFlag = TRUE)

  df_view <- complete(imp1, 1)


  # Responder models
  model_resp_c1m12 <- with(imp1, glm(responder_c1m12 ~ AGE + SEX + DIAGNOSIS + FAMILIARITY + 
      BMI + ANTIBODY + Num_Failed_Drugs + DETOXPRE + Psycopathological + Sleep_Disorders + 
      PREV_T0 + GGCEF_T0 + INT_T0, family = binomial))

  pooled_disc <- pool(model_resp_c1m12)

  model_resp_c2m12 <- with(imp1, glm(responder_c2m12 ~ AGE + SEX + DIAGNOSIS + FAMILIARITY + 
      BMI + ANTIBODY + Num_Failed_Drugs + DETOXPRE + Psycopathological + Sleep_Disorders + 
      PREV_T0 + GGCEF_T0 + INT_T0, family = binomial))

  pooled_disc <- pool(model_resp_c2m12)

  model_resp_c3m12 <- with(imp1, glm(responder_c3m12 ~ AGE + SEX + DIAGNOSIS + FAMILIARITY + 
      BMI + ANTIBODY + Num_Failed_Drugs + DETOXPRE + Psycopathological + Sleep_Disorders + 
      PREV_T0 + GGCEF_T0 + INT_T0, family = binomial))

  pooled_disc <- pool(model_resp_c3m12)

  completed_list <- complete(imp1, action = "all")

  library(pROC)

  roc_list <- vector("list", length(completed_list))

  auc_values <- numeric(length(completed_list))

  for (i in seq_along(completed_list)) {
      d <- completed_list[[i]]
      fit <- model_resp_c3m12$analyses[[i]]
      probs <- predict(fit, newdata = d, type = "response")
      roc_list[[i]] <- roc(d$responder_c3m12, probs, quiet = TRUE)
      auc_values[i] <- auc(roc_list[[i]])
  }

  mean_auc <- mean(auc_values)

  sd_auc <- sd(auc_values)

  fpr_grid <- seq(0, 1, length.out = 100)

  tpr_mat <- sapply(roc_list, function(r) {
      approx(1 - r$specificities, r$sensitivities, xout = fpr_grid)$y
  })

  mean_tpr <- rowMeans(tpr_mat, na.rm = TRUE)

  library(car)

  vif_values <- lapply(model_resp_c3m12$analyses, vif)

  vif_mean <- Reduce("+", vif_values)/length(vif_values)

  long_task2_list[[1]] <- long_task2_list[[1]] %>% mutate(TOTAL_MONTHS = case_when(CYCLE == 
      1 ~ MONTH, CYCLE == 2 ~ MONTH + 15, CYCLE == 3 ~ MONTH + 30))

  library(lme4)

  library(lmerTest)


  # Longitudinal mixed models
  null_model <- lmer(MMDs ~ 1 + (1 | SUBJECT_ID), data = long_task2_list[[1]], REML = FALSE)

  full_model <- lmer(MMDs ~ TOTAL_MONTHS + GGCEF_T0 + (1 | SUBJECT_ID), data = long_task2_list[[1]], 
      REML = FALSE)

  model_linear <- lmer(MMDs ~ TOTAL_MONTHS + GGCEF_T0 + (1 | SUBJECT_ID), data = long_task2_list[[1]], 
      REML = FALSE)

  model_categorical <- lmer(MMDs ~ factor(CYCLE) * factor(MONTH) + GGCEF_T0 + (1 | SUBJECT_ID), 
      data = long_task2_list[[1]], REML = FALSE)

  library(emmeans)

  mmd_estimates <- emmeans(model_linear, ~TOTAL_MONTHS, at = list(TOTAL_MONTHS = c(1, 42)))

  reduction_results <- pairs(mmd_estimates, reverse = TRUE)

  library(emmeans)

  library(ggplot2)

  emm_cat <- emmeans(model_categorical, ~CYCLE * MONTH, at = list(GGCEF_T0 = mean(long_task2_list[[1]]$GGCEF_T0, 
      na.rm = TRUE)))

  emm_df <- as.data.frame(emm_cat)

  library(lme4)

  library(ggplot2)

  pred_time <- data.frame(TOTAL_MONTHS = seq(min(long_task2_list[[1]]$TOTAL_MONTHS, na.rm = TRUE), 
      max(long_task2_list[[1]]$TOTAL_MONTHS, na.rm = TRUE), length.out = 200), GGCEF_T0 = mean(long_task2_list[[1]]$GGCEF_T0, 
      na.rm = TRUE))

  beta <- fixef(model_linear)

  V <- vcov(model_linear)

  tt <- delete.response(terms(model_linear))

  X <- model.matrix(tt, pred_time)

  X <- X[, names(beta), drop = FALSE]

  pred_time$fit <- as.numeric(X %*% beta)

  pred_time$se <- sqrt(diag(X %*% V %*% t(X)))

  pred_time$lwr <- pred_time$fit - 1.96 * pred_time$se

  pred_time$upr <- pred_time$fit + 1.96 * pred_time$se

  invisible(list(
    discontinuation = model_disc_imp,
    response = list(cycle1 = model_resp_c1m12, cycle2 = model_resp_c2m12, cycle3 = model_resp_c3m12),
    longitudinal = list(null = null_model, linear = model_linear, categorical = model_categorical),
    baseline_imputation = imp, response_imputation = imp1
  ))
}
