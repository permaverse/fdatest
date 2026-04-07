# Tests for one-sample IWT: R/one-iwt.R, R/ITP1bspline.R, R/ITP1fourier.R

d1 <- NASAtemp$milan[1:6, 1:8]
p <- ncol(d1)

# ---------------------------------------------------------------------------
# IWT1 — basic call (scalar mu, recycle = TRUE)
# ---------------------------------------------------------------------------
set.seed(42)
res <- IWT1(data = d1, mu = 0, B = 5L)

expect_inherits(res, "IWT1")
expect_equal(res$test, "1pop")
expect_equal(length(res$adjusted_pval), p)
expect_equal(length(res$unadjusted_pval), p)
expect_true(all(res$adjusted_pval >= 0 & res$adjusted_pval <= 1))
expect_true(all(res$unadjusted_pval >= 0 & res$unadjusted_pval <= 1))
expect_equal(dim(res$pval_matrix), c(p, p))
expect_equal(dim(res$data_eval), dim(d1))
expect_equal(res$mu, 0)

# Adjusted p-values are monotone corrections of unadjusted (>= unadjusted)
expect_true(all(res$adjusted_pval >= res$unadjusted_pval))

# ---------------------------------------------------------------------------
# IWT1 — vector mu (length = p)
# ---------------------------------------------------------------------------
set.seed(42)
mu_vec <- colMeans(d1)
res_mu <- IWT1(data = d1, mu = mu_vec, B = 5L)
expect_inherits(res_mu, "IWT1")
expect_equal(length(res_mu$mu), p)

# ---------------------------------------------------------------------------
# IWT1 — recycle = FALSE (non-recycled interval testing)
# ---------------------------------------------------------------------------
set.seed(42)
res_nr <- IWT1(data = d1, mu = 0, B = 5L, recycle = FALSE)
expect_inherits(res_nr, "IWT1")
expect_equal(length(res_nr$adjusted_pval), p)
# Upper triangle of pval_matrix should be NA (non-recycled)
expect_true(is.na(res_nr$pval_matrix[1, p]))

# ---------------------------------------------------------------------------
# iwt1() — new API returning fos class
# ---------------------------------------------------------------------------
set.seed(42)
res_fos <- iwt1(data = d1, mu = 0, n_perm = 5L)

expect_inherits(res_fos, "fos")
expect_equal(length(res_fos$adjusted_pvalues), p)
expect_equal(length(res_fos$unadjusted_pvalues), p)
expect_true(all(res_fos$adjusted_pvalues >= 0 & res_fos$adjusted_pvalues <= 1))
expect_true(all(
  res_fos$unadjusted_pvalues >= 0 & res_fos$unadjusted_pvalues <= 1
))
expect_equal(dim(res_fos$pvalue_matrix), c(p, p))
expect_equal(dim(res_fos$data), dim(d1))
expect_equal(res_fos$mu, 0)
expect_equal(res_fos$correction_method, "IWT")
expect_true(all(res_fos$adjusted_pvalues >= res_fos$unadjusted_pvalues))

# ---------------------------------------------------------------------------
# functional_one_sample_test() — interface returning fos class
# ---------------------------------------------------------------------------
set.seed(42)
res_fos2 <- functional_one_sample_test(data = d1, mu = 0, n_perm = 5L)
expect_inherits(res_fos2, "fos")
expect_equal(length(res_fos2$adjusted_pvalues), p)

# verbose = TRUE (exercises both cli_h1 progress paths)
set.seed(42)
res_fos2_v <- functional_one_sample_test(
  data = d1,
  mu = 0,
  n_perm = 5L,
  verbose = TRUE
)
expect_inherits(res_fos2_v, "fos")
expect_equal(length(res_fos2_v$adjusted_pvalues), p)

# recycle = FALSE
set.seed(42)
res_fos_nr <- iwt1(data = d1, mu = 0, n_perm = 5L, recycle = FALSE)
expect_inherits(res_fos_nr, "fos")
expect_true(is.na(res_fos_nr$pvalue_matrix[1, p]))

# ---------------------------------------------------------------------------
# Deprecated wrapper: ITP1bspline
# ---------------------------------------------------------------------------
set.seed(42)
res_bsp <- suppressWarnings(ITP1bspline(data = d1, mu = 0, B = 5L))
expect_inherits(res_bsp, "IWT1")

# ---------------------------------------------------------------------------
# Deprecated wrapper: ITP1fourier
# ---------------------------------------------------------------------------
set.seed(42)
res_fou <- suppressWarnings(ITP1fourier(data = d1, mu = 0, B = 5L))
expect_inherits(res_fou, "IWT1")

set.seed(NULL)
