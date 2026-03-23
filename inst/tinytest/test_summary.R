# Tests for summary methods: R/summary.fanova.R, R/summary.flm.R
library(fdatest)

data("NASAtemp", package = "fdatest")
d1 <- NASAtemp$milan[1:4, 1:8]
d2 <- NASAtemp$paris[1:4, 1:8]
temperature <- rbind(d1, d2)
groups <- c(rep(0L, 4L), rep(1L, 4L))

# ===========================================================================
# summary.fanova — structure checks with a real IWTaov result
# ===========================================================================
set.seed(42)
res_aov <- IWTaov(temperature ~ groups, B = 5L)
s_aov <- summary(res_aov)

expect_true(is.list(s_aov))
expect_true(!is.null(s_aov$call))
# factors data.frame: rows = factor names, cols include min p-value + signif
expect_true(is.data.frame(s_aov$factors))
expect_true("Minimum p-value" %in% names(s_aov$factors))
expect_equal(ncol(s_aov$factors), 2L) # min p-value + signif stars column
# R2 range: 2-element vector [min, max]
expect_equal(length(s_aov$R2), 2L)
expect_true(s_aov$R2[1] <= s_aov$R2[2])
# ftest: data.frame with min p-value + signif
expect_true(is.data.frame(s_aov$ftest))
expect_true("Minimum p-value" %in% names(s_aov$ftest))
expect_equal(ncol(s_aov$ftest), 2L) # min p-value + signif stars column

# ===========================================================================
# summary.fanova — significance code levels via mock objects
# The 5 levels: *** <0.001, ** <0.01, * <0.05, . <0.1, (space) >=0.1
# ===========================================================================

make_fanova <- function(p_factors, p_f) {
  # Minimal fanova mock to exercise summary.fanova
  n_pts <- length(p_f)
  obj <- list(
    call = quote(IWTaov(temperature ~ groups, B = 5L)),
    unadjusted_pval_factors = matrix(
      p_factors,
      nrow = 1L,
      ncol = n_pts,
      dimnames = list("groups", NULL)
    ),
    adjusted_pval_factors = matrix(
      p_factors,
      nrow = 1L,
      ncol = n_pts,
      dimnames = list("groups", NULL)
    ),
    unadjusted_pval_F = p_f,
    adjusted_pval_F = p_f,
    R2_eval = seq(0.2, 0.8, length.out = n_pts)
  )
  class(obj) <- "fanova"
  obj
}

codes <- c("***", "**", "*", ".", "")
p_vals <- c(0.0001, 0.005, 0.03, 0.07, 0.5)

for (i in seq_along(codes)) {
  pv <- rep(p_vals[i], 4L)
  mock <- make_fanova(p_factors = pv, p_f = pv)
  s_mock <- summary(mock)
  expect_true(
    s_mock$ftest[["Minimum p-value"]] <= p_vals[i] + 1e-9,
    info = paste("F-test min pval for pval =", p_vals[i])
  )
  expect_true(
    s_mock$ftest[[2]] == codes[i],
    info = paste("F-test signif code for pval =", p_vals[i])
  )
  expect_true(
    s_mock$factors[[2]] == codes[i],
    info = paste("factor signif code for pval =", p_vals[i])
  )
}

# ===========================================================================
# summary.flm — structure checks with a real IWTlm result
# ===========================================================================
set.seed(42)
res_lm <- IWTlm(temperature ~ groups, B = 5L)
s_lm <- summary(res_lm)

expect_true(is.list(s_lm))
expect_true(!is.null(s_lm$call))
# ttest: data.frame with rows for each coefficient
expect_true(is.data.frame(s_lm$ttest))
expect_true("Minimum p-value" %in% names(s_lm$ttest))
expect_equal(ncol(s_lm$ttest), 2L) # min p-value + signif stars column
# R2
expect_equal(length(s_lm$R2), 2L)
# ftest
expect_true(is.data.frame(s_lm$ftest))

# ===========================================================================
# summary.flm — significance code levels via mock objects
# ===========================================================================

make_flm <- function(p_part, p_f) {
  n_pts <- length(p_f)
  obj <- list(
    call = quote(IWTlm(temperature ~ groups, B = 5L)),
    unadjusted_pval_part = matrix(
      p_part,
      nrow = 2L,
      ncol = n_pts,
      dimnames = list(
        c("(Intercept)", "groups"),
        NULL
      )
    ),
    adjusted_pval_part = matrix(
      p_part,
      nrow = 2L,
      ncol = n_pts,
      dimnames = list(
        c("(Intercept)", "groups"),
        NULL
      )
    ),
    unadjusted_pval_F = p_f,
    adjusted_pval_F = p_f,
    R2_eval = seq(0.2, 0.8, length.out = n_pts)
  )
  class(obj) <- "flm"
  obj
}

for (i in seq_along(codes)) {
  pv <- rep(p_vals[i], 4L)
  mock <- make_flm(p_part = pv, p_f = pv)
  s_mock <- summary(mock)
  expect_true(
    s_mock$ftest[["Minimum p-value"]] <= p_vals[i] + 1e-9,
    info = paste("F-test min pval for pval =", p_vals[i])
  )
  expect_true(
    s_mock$ftest[[2]] == codes[i],
    info = paste("F-test signif code for pval =", p_vals[i])
  )
  # ttest has 2 rows; check at least the first
  expect_true(
    s_mock$ttest[[2]][1] == codes[i],
    info = paste("t-test signif code for pval =", p_vals[i])
  )
}
