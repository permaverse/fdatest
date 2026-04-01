# Tests for mirai parallel integration.
#
# Every if (mirai::daemons_set()) branch is exercised by setting up a single
# local daemon before each function call and tearing it down afterwards:
#
#   R/one-iwt.R   – iwt1  point-wise permutation path
#   R/utils-ts.R  – twosample_alt_permtest point-wise path
#                   (shared by IWT2, TWT2, FDR2, PCT2, Global2)
#   R/two-iwt.R   – iwt2  interval-wise row-computation path
#   R/utils-aov.R – aov_permtest point-wise path
#                   (shared by IWTaov, TWTaov, Globalaov)
#   R/aov-iwt.R   – iwt_aov interval-wise row-computation path
#   R/utils-lm.R  – lm_permtest point-wise path
#                   (shared by IWTlm, TWTlm, Globallm)
#   R/lm-iwt.R    – iwt_lm interval-wise row-computation path

d1 <- NASAtemp$milan[1:4, 1:8]
d2 <- NASAtemp$paris[1:4, 1:8]
temperature <- rbind(d1, d2) # 8 obs × 8 time pts
groups <- c(rep(0L, 4L), rep(1L, 4L))
grp_a <- c(0L, 0L, 1L, 1L, 0L, 0L, 1L, 1L)
grp_b <- c(0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L)
partition <- c(rep(1L, 4L), rep(2L, 4L))
p <- ncol(d1)

# Sanity-check that daemons_set() reflects the current state correctly.
mirai::daemons(1L, seed = 42L, dispatcher = FALSE)
expect_true(mirai::daemons_set())
mirai::daemons(0)
expect_false(mirai::daemons_set())

# ===========================================================================
# R/one-iwt.R — iwt1 point-wise permutation (mirai_map over seq_len(n_perm))
# ===========================================================================

# recycle = TRUE (default)
mirai::daemons(1L, seed = 42L, dispatcher = FALSE)
res_iwt1 <- IWT1(data = d1, mu = 0, B = 5L)
mirai::daemons(0)

expect_inherits(res_iwt1, "IWT1")
expect_equal(res_iwt1$test, "1pop")
expect_equal(length(res_iwt1$adjusted_pval), p)
expect_equal(length(res_iwt1$unadjusted_pval), p)
expect_true(all(res_iwt1$adjusted_pval >= 0 & res_iwt1$adjusted_pval <= 1))
expect_true(all(res_iwt1$unadjusted_pval >= 0 & res_iwt1$unadjusted_pval <= 1))
expect_equal(dim(res_iwt1$pval_matrix), c(p, p))
expect_true(all(res_iwt1$adjusted_pval >= res_iwt1$unadjusted_pval))

# recycle = FALSE — exercises the non-recycled interval loop with mirai
mirai::daemons(1L, seed = 42L, dispatcher = FALSE)
res_iwt1_nr <- IWT1(data = d1, mu = 0, B = 5L, recycle = FALSE)
mirai::daemons(0)

expect_inherits(res_iwt1_nr, "IWT1")
expect_equal(length(res_iwt1_nr$adjusted_pval), p)
expect_true(is.na(res_iwt1_nr$pval_matrix[1, p]))

# vector mu — same mirai path, different pre-processing
mirai::daemons(1L, seed = 42L, dispatcher = FALSE)
res_iwt1_mu <- IWT1(data = d1, mu = colMeans(d1), B = 5L)
mirai::daemons(0)

expect_inherits(res_iwt1_mu, "IWT1")
expect_equal(length(res_iwt1_mu$mu), p)

# ===========================================================================
# R/utils-ts.R — twosample_alt_permtest point-wise path (mirai_map over perms)
# R/two-iwt.R  — iwt2 interval-wise row-computation path (mirai_map over rows)
# ===========================================================================

# IWT2 exercises BOTH utils-ts.R and two-iwt.R mirai branches.
mirai::daemons(1L, seed = 42L, dispatcher = FALSE)
res_iwt2 <- IWT2(data1 = d1, data2 = d2, mu = 0, B = 5L, verbose = FALSE)
mirai::daemons(0)

expect_inherits(res_iwt2, "ftwosample")
expect_equal(length(res_iwt2$adjusted_pvalues), p)
expect_equal(length(res_iwt2$unadjusted_pvalues), p)
expect_true(all(
  res_iwt2$adjusted_pvalues >= 0 & res_iwt2$adjusted_pvalues <= 1
))
expect_true(all(
  res_iwt2$unadjusted_pvalues >= 0 & res_iwt2$unadjusted_pvalues <= 1
))
expect_equal(dim(res_iwt2$pvalue_matrix), c(p, p))
expect_true(all(res_iwt2$adjusted_pvalues >= res_iwt2$unadjusted_pvalues))

# recycle = FALSE — exercises the non-recycled row loop in iwt2 via mirai
mirai::daemons(1L, seed = 42L, dispatcher = FALSE)
res_iwt2_nr <- IWT2(
  data1 = d1,
  data2 = d2,
  mu = 0,
  B = 5L,
  recycle = FALSE,
  verbose = FALSE
)
mirai::daemons(0)

expect_inherits(res_iwt2_nr, "ftwosample")
expect_true(is.na(res_iwt2_nr$pvalue_matrix[1, p]))

# TWT2 — exercises only the utils-ts.R point-wise mirai branch
# (no interval-wise row loop in TWT)
mirai::daemons(1L, seed = 42L, dispatcher = FALSE)
res_twt2 <- TWT2(data1 = d1, data2 = d2, mu = 0, B = 5L, verbose = FALSE)
mirai::daemons(0)

expect_inherits(res_twt2, "ftwosample")
expect_equal(length(res_twt2$adjusted_pvalues), p)

# FDR2 — same utils-ts.R mirai path via a different entry point
mirai::daemons(1L, seed = 42L, dispatcher = FALSE)
res_fdr2 <- FDR2(data1 = d1, data2 = d2, mu = 0, B = 5L)
mirai::daemons(0)

expect_inherits(res_fdr2, "ftwosample")
expect_equal(length(res_fdr2$adjusted_pvalues), p)

# alternative = "greater" — exercises the alternative switch inside the
# parallel lambda (utils-ts.R)
mirai::daemons(1L, seed = 42L, dispatcher = FALSE)
res_iwt2_gr <- IWT2(
  data1 = d1,
  data2 = d2,
  B = 5L,
  alternative = "greater",
  verbose = FALSE
)
mirai::daemons(0)

expect_inherits(res_iwt2_gr, "ftwosample")

# alternative = "less"
mirai::daemons(1L, seed = 42L, dispatcher = FALSE)
res_iwt2_lt <- IWT2(
  data1 = d1,
  data2 = d2,
  B = 5L,
  alternative = "less",
  verbose = FALSE
)
mirai::daemons(0)

expect_inherits(res_iwt2_lt, "ftwosample")

# paired = TRUE — exercises the paired-permutation branch inside the parallel
# lambda (utils-ts.R line 7–14)
mirai::daemons(1L, seed = 42L, dispatcher = FALSE)
res_iwt2_p <- IWT2(
  data1 = d1,
  data2 = d2,
  B = 5L,
  paired = TRUE,
  verbose = FALSE
)
mirai::daemons(0)

expect_inherits(res_iwt2_p, "ftwosample")

# ===========================================================================
# R/utils-aov.R — aov_permtest point-wise path
# R/aov-iwt.R   — iwt_aov interval-wise row-computation path
# ===========================================================================

# IWTaov, method = "residuals" — exercises both aov mirai branches
mirai::daemons(1L, seed = 42L, dispatcher = FALSE)
res_iwt_aov <- IWTaov(temperature ~ groups, B = 5L, method = "residuals")
mirai::daemons(0)

expect_inherits(res_iwt_aov, "fanova")
expect_equal(length(res_iwt_aov$adjusted_pval_F), p)
expect_equal(length(res_iwt_aov$unadjusted_pval_F), p)
expect_true(all(
  res_iwt_aov$adjusted_pval_F >= 0 & res_iwt_aov$adjusted_pval_F <= 1
))
expect_equal(dim(res_iwt_aov$pval_matrix_F), c(p, p))
expect_equal(ncol(res_iwt_aov$adjusted_pval_factors), p)

# method = "responses" — exercises the "responses" branch inside
# .aov_one_perm when running in parallel (utils-aov.R lines 43–50)
mirai::daemons(1L, seed = 42L, dispatcher = FALSE)
res_iwt_aov_resp <- IWTaov(
  temperature ~ groups,
  B = 5L,
  method = "responses"
)
mirai::daemons(0)

expect_inherits(res_iwt_aov_resp, "fanova")
expect_equal(length(res_iwt_aov_resp$adjusted_pval_F), p)

# recycle = FALSE — exercises the non-recycled row loop in iwt_aov via mirai
mirai::daemons(1L, seed = 42L, dispatcher = FALSE)
res_iwt_aov_nr <- IWTaov(
  temperature ~ groups,
  B = 5L,
  recycle = FALSE
)
mirai::daemons(0)

expect_inherits(res_iwt_aov_nr, "fanova")
expect_true(is.na(res_iwt_aov_nr$pval_matrix_F[1, p]))

# nvar > 1 (multi-factor) — exercises more iterations of the parallel row loop
# in aov-iwt.R and the multi-variable inner loop in utils-aov.R
mirai::daemons(1L, seed = 42L, dispatcher = FALSE)
res_iwt_aov_2f <- IWTaov(
  temperature ~ grp_a + grp_b,
  B = 5L,
  method = "residuals"
)
mirai::daemons(0)

expect_inherits(res_iwt_aov_2f, "fanova")
expect_equal(length(res_iwt_aov_2f$adjusted_pval_F), p)

# nvar > 1 + factor() in formula — exercises the factor-renaming branch inside
# aov_permtest together with the mirai point-wise path
mirai::daemons(1L, seed = 42L, dispatcher = FALSE)
res_iwt_aov_fac <- IWTaov(
  temperature ~ grp_a + factor(grp_b),
  B = 5L,
  method = "residuals"
)
mirai::daemons(0)

expect_inherits(res_iwt_aov_fac, "fanova")

# TWTaov — exercises only the utils-aov.R point-wise mirai branch
mirai::daemons(1L, seed = 42L, dispatcher = FALSE)
res_twt_aov <- TWTaov(temperature ~ groups, B = 5L, method = "residuals")
mirai::daemons(0)

expect_inherits(res_twt_aov, "fanova")
expect_equal(length(res_twt_aov$adjusted_pval_F), p)

# TWTaov method = "responses"
mirai::daemons(1L, seed = 42L, dispatcher = FALSE)
res_twt_aov_resp <- TWTaov(temperature ~ groups, B = 5L, method = "responses")
mirai::daemons(0)

expect_inherits(res_twt_aov_resp, "fanova")

# ===========================================================================
# R/utils-lm.R — lm_permtest point-wise path
# R/lm-iwt.R   — iwt_lm interval-wise row-computation path
# ===========================================================================

# IWTlm, method = "residuals" — exercises both lm mirai branches
mirai::daemons(1L, seed = 42L, dispatcher = FALSE)
res_iwt_lm <- IWTlm(temperature ~ groups, B = 5L, method = "residuals")
mirai::daemons(0)

expect_inherits(res_iwt_lm, "flm")
expect_equal(length(res_iwt_lm$adjusted_pval_F), p)
expect_equal(length(res_iwt_lm$unadjusted_pval_F), p)
expect_true(all(
  res_iwt_lm$adjusted_pval_F >= 0 & res_iwt_lm$adjusted_pval_F <= 1
))
expect_equal(dim(res_iwt_lm$pval_matrix_F), c(p, p))
expect_equal(ncol(res_iwt_lm$adjusted_pval_part), p)

# method = "responses" — exercises the "responses" branch inside
# .lm_one_perm when running in parallel (utils-lm.R lines 32–37)
mirai::daemons(1L, seed = 42L, dispatcher = FALSE)
res_iwt_lm_resp <- IWTlm(
  temperature ~ groups,
  B = 5L,
  method = "responses"
)
mirai::daemons(0)

expect_inherits(res_iwt_lm_resp, "flm")
expect_equal(length(res_iwt_lm_resp$adjusted_pval_F), p)

# recycle = FALSE — exercises the non-recycled row loop in iwt_lm via mirai
mirai::daemons(1L, seed = 42L, dispatcher = FALSE)
res_iwt_lm_nr <- IWTlm(
  temperature ~ groups,
  B = 5L,
  recycle = FALSE
)
mirai::daemons(0)

expect_inherits(res_iwt_lm_nr, "flm")
expect_true(is.na(res_iwt_lm_nr$pval_matrix_F[1, p]))

# nvar > 1 — exercises the multi-predictor inner loop in lm-iwt.R and
# the extra residuals loop (ii from 2 to nvar+1) in utils-lm.R
mirai::daemons(1L, seed = 42L, dispatcher = FALSE)
res_iwt_lm_2v <- IWTlm(
  temperature ~ grp_a + grp_b,
  B = 5L,
  method = "residuals"
)
mirai::daemons(0)

expect_inherits(res_iwt_lm_2v, "flm")
expect_equal(dim(res_iwt_lm_2v$pval_matrix_part), c(3L, p, p))

# factor() in formula — exercises column-renaming branch in lm_permtest
# together with the mirai point-wise path
mirai::daemons(1L, seed = 42L, dispatcher = FALSE)
res_iwt_lm_fac <- IWTlm(
  temperature ~ factor(groups),
  B = 5L,
  method = "residuals"
)
mirai::daemons(0)

expect_inherits(res_iwt_lm_fac, "flm")

# intercept-only model — forces method = "responses" internally; exercises
# the nvar == 0 branch of lm_permtest together with the mirai path
mirai::daemons(1L, seed = 42L, dispatcher = FALSE)
res_iwt_lm_int <- IWTlm(temperature ~ 1, B = 5L)
mirai::daemons(0)

expect_inherits(res_iwt_lm_int, "flm")
expect_equal(length(res_iwt_lm_int$adjusted_pval_F), p)

# TWTlm — exercises only the utils-lm.R point-wise mirai branch
mirai::daemons(1L, seed = 42L, dispatcher = FALSE)
res_twt_lm <- TWTlm(temperature ~ groups, B = 5L, method = "residuals")
mirai::daemons(0)

expect_inherits(res_twt_lm, "flm")
expect_equal(length(res_twt_lm$adjusted_pval_F), p)

# TWTlm method = "responses"
mirai::daemons(1L, seed = 42L, dispatcher = FALSE)
res_twt_lm_resp <- TWTlm(temperature ~ groups, B = 5L, method = "responses")
mirai::daemons(0)

expect_inherits(res_twt_lm_resp, "flm")

Sys.sleep(1)
