# Tests for fd-object paths in R/utils.R:
#   onesample2coeffs(), twosamples2coeffs(), formula2coeff()
# fda is in Imports so it is always available; no skip needed.
library(fdatest)

# ---------------------------------------------------------------------------
# Shared fd fixtures
# ---------------------------------------------------------------------------
basis01 <- fda::create.bspline.basis(rangeval = c(0, 1), nbasis = 6L)
basis02 <- fda::create.bspline.basis(rangeval = c(0, 2), nbasis = 6L)

set.seed(1L)
n1 <- 5L
n2 <- 4L
fd1 <- fda::fd(matrix(stats::rnorm(6L * n1), nrow = 6L, ncol = n1), basis01)
fd2 <- fda::fd(matrix(stats::rnorm(6L * n2), nrow = 6L, ncol = n2), basis01)

# fd object on a different range — used for mismatch error tests
fd_bad <- fda::fd(
  matrix(stats::rnorm(6L * n1), nrow = 6L, ncol = n1),
  basis02
)

# scalar fd for mu (1 curve on same basis)
fd_mu <- fda::fd(matrix(stats::rnorm(6L), nrow = 6L, ncol = 1L), basis01)

# ---------------------------------------------------------------------------
# onesample2coeffs — fd data, scalar mu
# ---------------------------------------------------------------------------
out <- fdatest:::onesample2coeffs(fd1, mu = 0)
expect_true(is.matrix(out$coeff))
# p columns: (rangeval[2] - rangeval[1]) / dx + 1 with default dx = 0.01
expect_equal(nrow(out$coeff), n1)
expect_equal(out$mu, 0)

# explicit dx shrinks grid
out_dx <- fdatest:::onesample2coeffs(fd1, mu = 0, dx = 0.1)
expect_true(is.matrix(out_dx$coeff))
expect_equal(nrow(out_dx$coeff), n1)
expect_true(ncol(out_dx$coeff) < ncol(out$coeff))

# ---------------------------------------------------------------------------
# onesample2coeffs — fd data, fd mu (same rangeval)
# ---------------------------------------------------------------------------
out_fdmu <- fdatest:::onesample2coeffs(fd1, mu = fd_mu)
expect_true(is.matrix(out_fdmu$coeff))
expect_true(is.matrix(out_fdmu$mu))
expect_equal(nrow(out_fdmu$coeff), n1)
expect_equal(nrow(out_fdmu$mu), 1L) # 1 mu curve

# ---------------------------------------------------------------------------
# onesample2coeffs — fd mu with mismatched rangeval → error
# ---------------------------------------------------------------------------
expect_error(
  fdatest:::onesample2coeffs(fd1, mu = fd_bad),
  "range"
)

# ---------------------------------------------------------------------------
# twosamples2coeffs — both fd, scalar mu
# ---------------------------------------------------------------------------
out2 <- fdatest:::twosamples2coeffs(fd1, fd2, mu = 0)
expect_true(is.matrix(out2$coeff1))
expect_true(is.matrix(out2$coeff2))
expect_equal(nrow(out2$coeff1), n1)
expect_equal(nrow(out2$coeff2), n2)
expect_equal(out2$mu, 0)

# ---------------------------------------------------------------------------
# twosamples2coeffs — both fd, fd mu (same rangeval)
# ---------------------------------------------------------------------------
out2_fdmu <- fdatest:::twosamples2coeffs(fd1, fd2, mu = fd_mu)
expect_true(is.matrix(out2_fdmu$mu))

# ---------------------------------------------------------------------------
# twosamples2coeffs — mismatched data rangeval → error
# ---------------------------------------------------------------------------
expect_error(
  fdatest:::twosamples2coeffs(fd1, fd_bad, mu = 0),
  "range"
)

# ---------------------------------------------------------------------------
# twosamples2coeffs — fd mu with mismatched rangeval → error
# ---------------------------------------------------------------------------
fd_mu_bad <- fda::fd(
  matrix(stats::rnorm(6L), nrow = 6L, ncol = 1L),
  basis02
)
expect_error(
  fdatest:::twosamples2coeffs(fd1, fd2, mu = fd_mu_bad),
  "range"
)

# ---------------------------------------------------------------------------
# formula2coeff — fd on LHS of formula
# ---------------------------------------------------------------------------
groups <- seq_len(n1)
coeff_fd <- fdatest:::formula2coeff(fd1 ~ groups)
expect_true(is.matrix(coeff_fd))
expect_equal(nrow(coeff_fd), n1)

# explicit dx
coeff_fd_dx <- fdatest:::formula2coeff(fd1 ~ groups, dx = 0.1)
expect_true(is.matrix(coeff_fd_dx))
expect_true(ncol(coeff_fd_dx) < ncol(coeff_fd))

# ---------------------------------------------------------------------------
# End-to-end: IWT1 and IWT2 accept fd objects and produce plottable results
# ---------------------------------------------------------------------------
set.seed(2L)
res1 <- IWT1(fd1, mu = 0, B = 5L)
expect_true(inherits(res1, "IWT1"))
# IWT1 uses base graphics via plot.IWT1 — redirect to a temp device
grDevices::pdf(tempfile(fileext = ".pdf"))
plot(res1)
grDevices::dev.off()

set.seed(3L)
res2 <- IWT2(fd1, fd2, B = 5L, verbose = FALSE)
expect_true(inherits(res2, "ftwosample"))
p2 <- autoplot(res2)
expect_true(inherits(p2, "gg") || inherits(p2, "patchwork"))
