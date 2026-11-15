packages <- c("dplyr", "mice", "zoo", "ggplot2", "pROC", "car", "lme4", "lmerTest", "emmeans")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
