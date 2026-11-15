# Run from the repository root: Rscript --vanilla scripts/check.R
source("R/coursework.R")
source("scripts/synthetic_fixture.R")

check <- function() {
  directory <- tempfile("migraine-synthetic-")
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  stopifnot(inherits(try(run_coursework(directory), silent = TRUE), "try-error"))
  dir.create(directory)
  stopifnot(inherits(try(run_coursework(directory), silent = TRUE), "try-error"))
  make_synthetic_inputs(directory)
  # The original MICE calls print iteration progress; keep the check output short.
  warning_messages <- character()
  progress <- capture.output(result <- withCallingHandlers(
    run_coursework(directory), warning = function(w) {
      warning_messages <<- c(warning_messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  ))
  stopifnot(length(result$discontinuation$analyses) == 20L)
  stopifnot(all(vapply(result$response, function(x) length(x$analyses) == 20L, logical(1))))
  future <- c("Suspension", "MONTHS_OF_TREAT", "TREATMENT_DISC")
  stopifnot(all(result$baseline_imputation$predictorMatrix[, future] == 0))
  stopifnot(all(result$baseline_imputation$method[future] == ""))
  stopifnot(result$baseline_imputation$method["GGCEF_T0"] == "")
  excluded <- c("SUBJECT_ID", "responder_c1m12", "responder_c2m12", "responder_c3m12")
  stopifnot(all(result$response_imputation$predictorMatrix[, excluded] == 0))
  stopifnot(all(result$response_imputation$method[excluded] == ""))
  stopifnot(all(vapply(result$longitudinal, inherits, logical(1), what = "lmerMod")))
  stopifnot(length(lme4::findbars(formula(result$longitudinal$linear))) == 1L)
  stopifnot(identical(lme4::findbars(formula(result$longitudinal$linear))[[1]][[2]], 1))
  stopifnot(all(is.finite(lme4::fixef(result$longitudinal$linear))))
  stopifnot(all(vapply(result$longitudinal, function(x) {
    length(x@optinfo$conv$lme4$messages) == 0L
  }, logical(1))))
  stopifnot(all(is.finite(mice::pool(result$discontinuation)$pooled$estimate)))
  for (model in result$response) {
    stopifnot(all(is.finite(mice::pool(model)$pooled$estimate)))
  }
  cat("PASS: missing-input guards, synthetic execution, imputation exclusions,\n")
  cat("      pooled logistic models and random-intercept longitudinal models.\n")
  cat("Synthetic inputs were created in a temporary directory and will be removed.\n")
  if (length(warning_messages)) print(table(warning_messages))
  for (name in c("baseline_imputation", "response_imputation")) {
    events <- result[[name]]$loggedEvents
    if (!is.null(events)) {
      cat(name, "imputation events:\n")
      print(events)
    }
  }
}

check()
