# Tests for utility functions in R/utils.R

d1 <- NASAtemp$milan[1:4, 1:8]
d2 <- NASAtemp$paris[1:4, 1:8]

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
# extract_residuals, extract_fitted
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

corrected <- fdatest:::pval_correct_cpp(pm)
expect_equal(dim(corrected), c(4L, 4L))
expect_true(all(corrected >= 0 & corrected <= 1, na.rm = TRUE))

# Corrected values form a valid p-value matrix
expect_true(!any(corrected < 0, na.rm = TRUE))

# for the C++ version: Expect matrix, not array
array_pm <- array(stats::runif(16), dim = c(4, 4, 1))
expect_error(fdatest:::pval_correct_cpp(array_pm))

# for the C++ version: Expect square matrix
rect_pm <- matrix(stats::runif(12), nrow = 4, ncol = 3)
expect_error(fdatest:::pval_correct_cpp(rect_pm))

# for the C++ version: add example with NA values (NA in a non-last row)
pm_with_na <- matrix(stats::runif(16), nrow = 4, ncol = 4)
pm_with_na[1, 1] <- NA
corrected_with_na <- fdatest:::pval_correct_cpp(pm_with_na)
expect_equal(dim(corrected_with_na), c(4L, 4L))
expect_true(all(corrected_with_na >= 0 & corrected_with_na <= 1, na.rm = TRUE))

# for the C++ version: NA in the last row makes pval_var start as NaN,
# exercising the `pval_var = seg_max` branch (std::isnan(pval_var) == true)
pm_with_na_lastrow <- matrix(stats::runif(16), nrow = 4, ncol = 4)
pm_with_na_lastrow[4, 4] <- NA
corrected_na_last <- fdatest:::pval_correct_cpp(pm_with_na_lastrow)
expect_equal(dim(corrected_na_last), c(4L, 4L))
expect_true(all(corrected_na_last >= 0 & corrected_na_last <= 1, na.rm = TRUE))

# ---------------------------------------------------------------------------
# os_to_coeffs — matrix input
# ---------------------------------------------------------------------------
out <- fdatest:::os_to_coeffs(d1, mu = 0)
expect_identical(out$coeff, d1)
expect_equal(out$mu, 0)

# Vector mu
mu_vec <- colMeans(d1)
out2 <- fdatest:::os_to_coeffs(d1, mu = mu_vec)
expect_equal(out2$mu, mu_vec)

# Invalid data type → error
expect_error(fdatest:::os_to_coeffs(list(1, 2), mu = 0))
# Invalid mu type → error  (matrix is not a vector, so it triggers the error)
expect_error(fdatest:::os_to_coeffs(d1, mu = matrix(1, 1, 1)))

# ---------------------------------------------------------------------------
# ts_to_coeffs — matrix input
# ---------------------------------------------------------------------------
out4 <- fdatest:::ts_to_coeffs(d1, d2, mu = 0)
expect_identical(out4$coeff1, d1)
expect_identical(out4$coeff2, d2)
expect_equal(out4$mu, 0)

# Vector mu
out5 <- fdatest:::ts_to_coeffs(d1, d2, mu = mu_vec)
expect_equal(out5$mu, mu_vec)

# Invalid data1 type → error
expect_error(fdatest:::ts_to_coeffs(list(), d2, mu = 0))
# Invalid data2 type → error
expect_error(fdatest:::ts_to_coeffs(d1, list(), mu = 0))
# Invalid mu type → error  (matrix is not a vector)
expect_error(fdatest:::ts_to_coeffs(d1, d2, mu = matrix(1, 1, 1)))

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
expect_equal(ncol(dm), 2L) # intercept and groups8
expect_true("(Intercept)" %in% colnames(dm))

set.seed(NULL)
