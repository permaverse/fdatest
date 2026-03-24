# Tests for functional ANOVA:
# R/aov-iwt.R, R/aov-twt.R, R/aov-global.R,
# R/functional-anova-test.R, R/ITPaovbspline.R
library(fdatest)

data("NASAtemp", package = "fdatest")
d1 <- NASAtemp$milan[1:4, 1:8]
d2 <- NASAtemp$paris[1:4, 1:8]
temperature <- rbind(d1, d2) # 8 obs × 8 time pts
groups <- c(rep(0L, 4L), rep(1L, 4L))
p <- ncol(d1)

# ===========================================================================
# IWTaov — method = "residuals"
# ===========================================================================
set.seed(42)
res_iwt_r <- IWTaov(temperature ~ groups, B = 5L, method = "residuals")

expect_inherits(res_iwt_r, "fanova")
expect_equal(length(res_iwt_r$adjusted_pval_F), p)
expect_equal(length(res_iwt_r$unadjusted_pval_F), p)
expect_true(all(
  res_iwt_r$adjusted_pval_F >= 0 & res_iwt_r$adjusted_pval_F <= 1
))
expect_equal(dim(res_iwt_r$pval_matrix_F), c(p, p))
# factors
expect_equal(ncol(res_iwt_r$unadjusted_pval_factors), p)
expect_equal(ncol(res_iwt_r$adjusted_pval_factors), p)
# data, fitted, residuals, R2
expect_equal(dim(res_iwt_r$data_eval), dim(temperature))
expect_equal(dim(res_iwt_r$fitted_eval), dim(temperature))
expect_equal(dim(res_iwt_r$residuals_eval), dim(temperature))
expect_equal(length(res_iwt_r$R2_eval), p)

# ===========================================================================
# IWTaov — method = "responses"
# ===========================================================================
set.seed(42)
res_iwt_resp <- IWTaov(temperature ~ groups, B = 5L, method = "responses")
expect_inherits(res_iwt_resp, "fanova")
expect_equal(length(res_iwt_resp$adjusted_pval_F), p)

# ===========================================================================
# IWTaov — recycle = FALSE
# ===========================================================================
set.seed(42)
res_iwt_nr <- IWTaov(temperature ~ groups, B = 5L, recycle = FALSE)
expect_inherits(res_iwt_nr, "fanova")
expect_true(is.na(res_iwt_nr$pval_matrix_F[1, p]))

# ===========================================================================
# TWTaov — method = "residuals"
# ===========================================================================
set.seed(42)
res_twt_r <- TWTaov(temperature ~ groups, B = 5L, method = "residuals")
expect_inherits(res_twt_r, "fanova")
expect_equal(length(res_twt_r$adjusted_pval_F), p)
expect_equal(ncol(res_twt_r$adjusted_pval_factors), p)
# No pval_matrix fields in TWT
expect_null(res_twt_r$pval_matrix_F)

# ===========================================================================
# TWTaov — method = "responses"
# ===========================================================================
set.seed(42)
res_twt_resp <- TWTaov(temperature ~ groups, B = 5L, method = "responses")
expect_inherits(res_twt_resp, "fanova")

# ===========================================================================
# Globalaov — stat = "Integral", method = "residuals"
# ===========================================================================
set.seed(42)
res_g_r <- Globalaov(
  temperature ~ groups,
  B = 5L,
  stat = "Integral",
  method = "residuals"
)
expect_inherits(res_g_r, "fanova")
expect_equal(length(res_g_r$adjusted_pval_F), p)
expect_true(!is.null(res_g_r$Global_pval_F))
expect_true(res_g_r$Global_pval_F >= 0 && res_g_r$Global_pval_F <= 1)
expect_true(!is.null(res_g_r$Global_pval_factors))

# ===========================================================================
# Globalaov — stat = "Max"
# ===========================================================================
set.seed(42)
res_g_max <- Globalaov(temperature ~ groups, B = 5L, stat = "Max")
expect_inherits(res_g_max, "fanova")
expect_true(!is.null(res_g_max$Global_pval_F))

# ===========================================================================
# Globalaov — method = "responses"
# ===========================================================================
set.seed(42)
res_g_resp <- Globalaov(temperature ~ groups, B = 5L, method = "responses")
expect_inherits(res_g_resp, "fanova")

# ===========================================================================
# functional_anova_test — all corrections
# ===========================================================================
set.seed(42)
res_ft_iwt <- functional_anova_test(
  temperature ~ groups,
  correction = "IWT",
  B = 5L
)
expect_inherits(res_ft_iwt, "fanova")
expect_equal(res_ft_iwt$correction, "IWT")

set.seed(42)
res_ft_twt <- functional_anova_test(
  temperature ~ groups,
  correction = "TWT",
  B = 5L
)
expect_inherits(res_ft_twt, "fanova")
expect_equal(res_ft_twt$correction, "TWT")

set.seed(42)
res_ft_g <- functional_anova_test(
  temperature ~ groups,
  correction = "Global",
  B = 5L
)
expect_inherits(res_ft_g, "fanova")
expect_equal(res_ft_g$correction, "Global")

# ===========================================================================
# IWTaov — nvar > 1, factor() in formula, method = "residuals"
# Exercises lines 317–319 of aov_permtest (factor renaming block)
# ===========================================================================
grpA <- c(0L, 0L, 1L, 1L, 0L, 0L, 1L, 1L)
grpB <- c(0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L)

set.seed(42)
res_iwt_rf <- IWTaov(
  temperature ~ grpA + factor(grpB),
  B = 5L,
  method = "residuals"
)
expect_inherits(res_iwt_rf, "fanova")
expect_equal(length(res_iwt_rf$adjusted_pval_F), p)

# ===========================================================================
# IWTaov — nvar > 1, interaction term, method = "residuals"
# Exercises lines 331–337 of aov_permtest (the ":" interaction sub-branch)
# ===========================================================================
set.seed(42)
res_iwt_ia <- IWTaov(
  temperature ~ grpA * grpB,
  B = 5L,
  method = "residuals"
)
expect_inherits(res_iwt_ia, "fanova")
expect_equal(length(res_iwt_ia$adjusted_pval_F), p)

# ===========================================================================
# IWTaov — intercept-only (nvar = 0)
# Exercises lines 309–310 (nvar == 0 t0_glob branch) and
# lines 374–376 (sign-permutation branch in the permutation loop)
# ===========================================================================
set.seed(42)
res_iwt_int0 <- suppressWarnings(IWTaov(temperature ~ 1, B = 5L))
expect_inherits(res_iwt_int0, "fanova")

# ===========================================================================
# Deprecated wrapper: ITPaovbspline
# ===========================================================================
set.seed(42)
res_itp_aov <- suppressWarnings(ITPaovbspline(temperature ~ groups, B = 5L))
expect_inherits(res_itp_aov, "fanova")
