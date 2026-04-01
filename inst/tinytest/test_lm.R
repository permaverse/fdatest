# Tests for functional linear models:
# R/lm-iwt.R, R/lm-twt.R, R/lm-global.R,
# R/functional-lm-test.R, R/ITPlmbspline.R

d1 <- NASAtemp$milan[1:4, 1:8]
d2 <- NASAtemp$paris[1:4, 1:8]
temperature <- rbind(d1, d2) # 8 obs × 8 time pts
groups <- c(rep(0L, 4L), rep(1L, 4L))
p <- ncol(d1)

# ===========================================================================
# IWTlm — method = "residuals" (default)
# ===========================================================================
set.seed(42)
res_iwt_r <- IWTlm(temperature ~ groups, B = 5L, method = "residuals")

expect_inherits(res_iwt_r, "flm")
expect_equal(length(res_iwt_r$adjusted_pval_F), p)
expect_equal(length(res_iwt_r$unadjusted_pval_F), p)
expect_true(all(
  res_iwt_r$adjusted_pval_F >= 0 & res_iwt_r$adjusted_pval_F <= 1
))
expect_equal(dim(res_iwt_r$pval_matrix_F), c(p, p))
# partial t-tests
expect_equal(ncol(res_iwt_r$unadjusted_pval_part), p)
expect_equal(ncol(res_iwt_r$adjusted_pval_part), p)
expect_equal(dim(res_iwt_r$pval_matrix_part), c(2L, p, p)) # intercept and groups
# data, fitted, residuals, R2
expect_equal(dim(res_iwt_r$data_eval), dim(temperature))
expect_equal(dim(res_iwt_r$fitted_eval), dim(temperature))
expect_equal(dim(res_iwt_r$residuals_eval), dim(temperature))
expect_equal(length(res_iwt_r$R2_eval), p)

# ===========================================================================
# IWTlm — method = "residuals", nvar > 1 (multi-predictor)
# Exercises the paste(variables_reduced, ...) reduced-formula branch
# ===========================================================================
grp_a <- c(0L, 0L, 1L, 1L, 0L, 0L, 1L, 1L)
grp_b <- c(0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L)

set.seed(42)
res_iwt_r2 <- IWTlm(temperature ~ grp_a + grp_b, B = 5L, method = "residuals")
expect_inherits(res_iwt_r2, "flm")
expect_equal(length(res_iwt_r2$adjusted_pval_F), p)
# pval_matrix_part has one slice per predictor + intercept
expect_equal(dim(res_iwt_r2$pval_matrix_part), c(3L, p, p))

# IWTlm — method = "residuals", factor() in formula
# (exercises column-renaming branch; intercept formula uses renamed columns)
set.seed(42)
res_iwt_rf <- IWTlm(
  temperature ~ factor(groups),
  B = 5L,
  method = "residuals"
)
expect_inherits(res_iwt_rf, "flm")
expect_equal(length(res_iwt_rf$adjusted_pval_F), p)

# ===========================================================================
# IWTlm — method = "responses"
# ===========================================================================
set.seed(42)
res_iwt_resp <- IWTlm(temperature ~ groups, B = 5L, method = "responses")
expect_inherits(res_iwt_resp, "flm")
expect_equal(length(res_iwt_resp$adjusted_pval_F), p)

# ===========================================================================
# IWTlm — recycle = FALSE
# ===========================================================================
set.seed(42)
res_iwt_nr <- IWTlm(temperature ~ groups, B = 5L, recycle = FALSE)
expect_inherits(res_iwt_nr, "flm")
expect_true(is.na(res_iwt_nr$pval_matrix_F[1, p]))

# ===========================================================================
# IWTlm — intercept-only model (nvar = 0 path, forces method = "responses")
# ===========================================================================
set.seed(42)
res_iwt_int <- IWTlm(temperature ~ 1, B = 5L)
expect_inherits(res_iwt_int, "flm")
expect_equal(length(res_iwt_int$adjusted_pval_F), p)

# ===========================================================================
# TWTlm — method = "residuals"
# ===========================================================================
set.seed(42)
res_twt_r <- TWTlm(temperature ~ groups, B = 5L, method = "residuals")
expect_inherits(res_twt_r, "flm")
expect_equal(length(res_twt_r$adjusted_pval_F), p)
expect_equal(ncol(res_twt_r$adjusted_pval_part), p)
# No pval_matrix fields in TWT
expect_null(res_twt_r$pval_matrix_F)

# TWTlm — method = "residuals", nvar > 1
set.seed(42)
res_twt_r2 <- TWTlm(temperature ~ grp_a + grp_b, B = 5L, method = "residuals")
expect_inherits(res_twt_r2, "flm")
expect_equal(length(res_twt_r2$adjusted_pval_F), p)

# TWTlm — method = "residuals", factor() in formula
set.seed(42)
res_twt_rf <- TWTlm(
  temperature ~ factor(groups),
  B = 5L,
  method = "residuals"
)
expect_inherits(res_twt_rf, "flm")
expect_equal(length(res_twt_rf$adjusted_pval_F), p)

# ===========================================================================
# TWTlm — method = "responses"
# ===========================================================================
set.seed(42)
res_twt_resp <- TWTlm(temperature ~ groups, B = 5L, method = "responses")
expect_inherits(res_twt_resp, "flm")

# ===========================================================================
# Globallm — stat = "Integral", method = "residuals"
# ===========================================================================
set.seed(42)
res_g_r <- Globallm(
  temperature ~ groups,
  B = 5L,
  stat = "Integral",
  method = "residuals"
)
expect_inherits(res_g_r, "flm")
expect_equal(length(res_g_r$adjusted_pval_F), p)
expect_true(!is.null(res_g_r$Global_pval_F))
expect_true(res_g_r$Global_pval_F >= 0 && res_g_r$Global_pval_F <= 1)
expect_true(!is.null(res_g_r$Global_pval_part))

# Globallm — method = "residuals", nvar > 1
set.seed(42)
res_g_r2 <- Globallm(
  temperature ~ grp_a + grp_b,
  B = 5L,
  stat = "Integral",
  method = "residuals"
)
expect_inherits(res_g_r2, "flm")
expect_equal(length(res_g_r2$adjusted_pval_F), p)

# Globallm — method = "residuals", factor() in formula
set.seed(42)
res_g_rf <- Globallm(
  temperature ~ factor(groups),
  B = 5L,
  method = "residuals"
)
expect_inherits(res_g_rf, "flm")
expect_equal(length(res_g_rf$adjusted_pval_F), p)

# ===========================================================================
# Globallm — stat = "Max"
# ===========================================================================
set.seed(42)
res_g_max <- Globallm(temperature ~ groups, B = 5L, stat = "Max")
expect_inherits(res_g_max, "flm")
expect_true(!is.null(res_g_max$Global_pval_F))

# ===========================================================================
# Globallm — method = "responses"
# ===========================================================================
set.seed(42)
res_g_resp <- Globallm(temperature ~ groups, B = 5L, method = "responses")
expect_inherits(res_g_resp, "flm")

# ===========================================================================
# functional_lm_test — all corrections
# ===========================================================================
set.seed(42)
res_ft_iwt <- functional_lm_test(
  temperature ~ groups,
  correction = "IWT",
  B = 5L
)
expect_inherits(res_ft_iwt, "flm")
expect_equal(res_ft_iwt$correction, "IWT")

set.seed(42)
res_ft_twt <- functional_lm_test(
  temperature ~ groups,
  correction = "TWT",
  B = 5L
)
expect_inherits(res_ft_twt, "flm")
expect_equal(res_ft_twt$correction, "TWT")

set.seed(42)
res_ft_g <- functional_lm_test(
  temperature ~ groups,
  correction = "Global",
  B = 5L
)
expect_inherits(res_ft_g, "flm")
expect_equal(res_ft_g$correction, "Global")

# ===========================================================================
# Deprecated wrapper: ITPlmbspline
# ===========================================================================
set.seed(42)
res_itp_lm <- suppressWarnings(ITPlmbspline(temperature ~ groups, B = 5L))
expect_inherits(res_itp_lm, "flm")

set.seed(NULL)
