# Tests for utility functions in R/utils.R
library(fdatest)

data("NASAtemp", package = "fdatest")
d1 <- NASAtemp$milan[1:4, 1:8]
d2 <- NASAtemp$paris[1:4, 1:8]

# ---------------------------------------------------------------------------
# AVAILABLE_* constant helpers
# ---------------------------------------------------------------------------
expect_identical(
  fdatest:::AVAILABLE_ALTERNATIVES(),
  c("two.sided", "less", "greater")
)

expect_identical(
  fdatest:::AVAILABLE_METHODS(),
  c("IWT", "TWT", "PCT", "Global", "FDR")
)

expect_identical(
  fdatest:::AVAILABLE_STATISTICS(),
  c("Integral", "Max", "Integral_std", "Max_std")
)

# ---------------------------------------------------------------------------
# stat_lm_glob
# ---------------------------------------------------------------------------
fit_lm <- stats::lm(d1[, 1] ~ seq_len(nrow(d1)))
f_stat <- fdatest:::stat_lm_glob(fit_lm)
expect_true(is.numeric(f_stat))
expect_true(f_stat >= 0)

# ---------------------------------------------------------------------------
# stat_aov_part
# ---------------------------------------------------------------------------
groups <- c(0, 0, 1, 1)
fit_aov <- stats::aov(d1[, 1] ~ factor(groups))
aov_stats <- fdatest:::stat_aov_part(fit_aov)
expect_true(is.numeric(aov_stats))
expect_true(length(aov_stats) == 1L) # one factor

# ---------------------------------------------------------------------------
# extract_residuals / extract_fitted
# ---------------------------------------------------------------------------
expect_equal(fdatest:::extract_residuals(fit_lm), stats::residuals(fit_lm))
expect_equal(fdatest:::extract_fitted(fit_lm), stats::fitted(fit_lm))

# ---------------------------------------------------------------------------
# pval_correct — monotone correction
# ---------------------------------------------------------------------------
set.seed(1)
pm <- matrix(stats::runif(16), nrow = 4, ncol = 4)
pm[4, ] <- stats::runif(4, 0.1, 0.9) # bottom row (diag = 1)
corrected <- fdatest:::pval_correct(pm)
expect_equal(dim(corrected), c(4L, 4L))
expect_true(all(corrected >= 0 & corrected <= 1, na.rm = TRUE))

# Corrected values form a valid p-value matrix
expect_true(!any(corrected < 0, na.rm = TRUE))

# ---------------------------------------------------------------------------
# onesample2coeffs — matrix input
# ---------------------------------------------------------------------------
out <- fdatest:::onesample2coeffs(d1, mu = 0)
expect_identical(out$coeff, d1)
expect_equal(out$mu, 0)

# Vector mu
mu_vec <- colMeans(d1)
out2 <- fdatest:::onesample2coeffs(d1, mu = mu_vec)
expect_equal(out2$mu, mu_vec)

# Invalid data type → error
expect_error(fdatest:::onesample2coeffs(list(1, 2), mu = 0))
# Invalid mu type → error  (matrix is not a vector, so it triggers the error)
expect_error(fdatest:::onesample2coeffs(d1, mu = matrix(1, 1, 1)))

# ---------------------------------------------------------------------------
# twosamples2coeffs — matrix input
# ---------------------------------------------------------------------------
out4 <- fdatest:::twosamples2coeffs(d1, d2, mu = 0)
expect_identical(out4$coeff1, d1)
expect_identical(out4$coeff2, d2)
expect_equal(out4$mu, 0)

# Vector mu
out5 <- fdatest:::twosamples2coeffs(d1, d2, mu = mu_vec)
expect_equal(out5$mu, mu_vec)

# Invalid data1 type → error
expect_error(fdatest:::twosamples2coeffs(list(), d2, mu = 0))
# Invalid data2 type → error
expect_error(fdatest:::twosamples2coeffs(d1, list(), mu = 0))
# Invalid mu type → error  (matrix is not a vector)
expect_error(fdatest:::twosamples2coeffs(d1, d2, mu = matrix(1, 1, 1)))

# ---------------------------------------------------------------------------
# formula2coeff — extracts left-hand side matrix
# ---------------------------------------------------------------------------
groups8 <- c(rep(0L, 4L), rep(1L, 4L))
temperature <- rbind(d1, d2)

coeff <- fdatest:::formula2coeff(temperature ~ groups8)
expect_identical(coeff, temperature)

# Error on non-matrix, non-fd object
not_matrix <- as.data.frame(d1)
expect_error(fdatest:::formula2coeff(not_matrix ~ groups8))

# ---------------------------------------------------------------------------
# formula2design_matrix — builds design matrix
# ---------------------------------------------------------------------------
dm <- fdatest:::formula2design_matrix(temperature ~ groups8, coeff)
expect_true(is.matrix(dm))
expect_equal(nrow(dm), nrow(temperature))
expect_equal(ncol(dm), 2L) # intercept + groups8
expect_true("(Intercept)" %in% colnames(dm))
