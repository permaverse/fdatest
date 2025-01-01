stat_lm_glob <- function(anova) {
  stats::summary.lm(anova)$f[1]
}

stat_aov_part <- function(anova) {
  out <- summary(anova)[[1]][, 4]
  out <- out[-length(out)]
  out
}

extract_residuals <- function(x) {
  x$residuals
}
extract_fitted <- function(x) {
  x$fitted
}
