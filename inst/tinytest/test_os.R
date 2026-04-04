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
