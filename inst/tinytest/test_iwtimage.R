# Tests for IWTimage (R/IWTimage.R) and deprecated ITPimage (R/ITPimage.R)
#
# IWTimage dispatches on class:
#   * "IWT1" or "IWT2" → branch 1 (only "IWT1" is reachable from the public API)
#   * "IWTaov"         → branch 2 (dead code: IWTaov() now returns "fanova")
#
# To achieve 100% line coverage we create minimal fake objects with the old
# class names where the current public API can no longer produce them.
library(fdatest)

data("NASAtemp", package = "fdatest")
d1 <- NASAtemp$milan[1:4, 1:8]
d2 <- NASAtemp$paris[1:4, 1:8]
p <- ncol(d1)

# Helper: redirect base-graphics output to a temp PDF, clean up after
with_null_device <- function(expr) {
  tmp <- tempfile(fileext = ".pdf")
  grDevices::pdf(tmp, width = 10, height = 8)
  tryCatch(force(expr), finally = grDevices::dev.off())
  unlink(tmp)
  invisible(NULL)
}

# ===========================================================================
# IWTimage — IWT1 class (real object, 1-pop branch, scalar mu)
# ===========================================================================
set.seed(42)
res_iwt1 <- IWT1(d1, mu = 0, B = 5L)

with_null_device({
  IWTimage(res_iwt1, abscissa_range = c(0, 1))
})
expect_true(TRUE) # reached without error

# ===========================================================================
# IWTimage — IWT1 class, plot_unadjusted = TRUE
# ===========================================================================
with_null_device({
  IWTimage(res_iwt1, abscissa_range = c(0, 1), plot_unadjusted = TRUE)
})
expect_true(TRUE)

# ===========================================================================
# IWTimage — fake "IWT2" object (exercises the 2-pop sub-branch of branch 1)
#
# This branch is dead code in the current codebase (IWT2() now sets class
# "ftwosample", not "IWT2"), but the lines exist in IWTimage and must be
# covered.  We create a minimal compatible list with the old class.
# ===========================================================================
fake_iwt2 <- list(
  test = "2pop",
  mu = 0,
  adjusted_pval = rep(0.3, p),
  unadjusted_pval = rep(0.3, p),
  pval_matrix = matrix(0.3, nrow = p, ncol = p),
  data_eval = rbind(d1, d2),
  ord_labels = c(rep(1L, nrow(d1)), rep(2L, nrow(d2)))
)
class(fake_iwt2) <- "IWT2"

with_null_device({
  IWTimage(fake_iwt2, abscissa_range = c(0, 1))
})
expect_true(TRUE)

# plot_unadjusted = TRUE for IWT2 fake (exercises that lines branch too)
with_null_device({
  IWTimage(fake_iwt2, abscissa_range = c(0, 1), plot_unadjusted = TRUE)
})
expect_true(TRUE)

# IWT1 where some adjusted p-values are < alpha (exercises the rect-drawing loop)
res_iwt1_low <- res_iwt1
res_iwt1_low$adjusted_pval <- rep(0.01, p) # force all < 0.05
with_null_device({
  IWTimage(res_iwt1_low, alpha = 0.05, abscissa_range = c(0, 1))
})
expect_true(TRUE)

# ===========================================================================
# IWTimage — fake "IWTaov" object (exercises branch 2 — IWTaov dead code)
#
# IWTimage branch 2 uses old field names: pval_matrix_F, pval_matrix_factor,
# adjusted_pval_factor, unadjusted_pval_factor.  IWTaov() now stores
# "factors" (with 's') and returns class "fanova".  We craft a fake object
# with the exact fields that IWTimage expects.
# ===========================================================================
nvar <- 1L # one factor: "groups"
fake_iwt_aov <- list(
  adjusted_pval_F = rep(0.3, p),
  unadjusted_pval_F = rep(0.4, p),
  pval_matrix_F = matrix(0.3, nrow = p, ncol = p),
  adjusted_pval_factor = matrix(
    0.3,
    nrow = nvar,
    ncol = p,
    dimnames = list("groups", NULL)
  ),
  unadjusted_pval_factor = matrix(
    0.4,
    nrow = nvar,
    ncol = p,
    dimnames = list("groups", NULL)
  ),
  pval_matrix_factor = array(0.3, dim = c(nvar, p, p)),
  data_eval = rbind(d1, d2),
  design_matrix = cbind(
    `(Intercept)` = rep(1L, nrow(d1) + nrow(d2)),
    groups = c(rep(0L, nrow(d1)), rep(1L, nrow(d2)))
  )
)
class(fake_iwt_aov) <- "IWTaov"

with_null_device({
  IWTimage(fake_iwt_aov, abscissa_range = c(0, 1))
})
expect_true(TRUE)

# plot_unadjusted = TRUE with IWTaov fake
with_null_device({
  IWTimage(fake_iwt_aov, abscissa_range = c(0, 1), plot_unadjusted = TRUE)
})
expect_true(TRUE)

# IWTaov fake where some p-values < alpha (exercises rect-drawing loops)
fake_iwt_aov_low <- fake_iwt_aov
fake_iwt_aov_low$adjusted_pval_F <- rep(0.01, p)
fake_iwt_aov_low$adjusted_pval_factor <- matrix(
  0.01,
  nrow = nvar,
  ncol = p,
  dimnames = list("groups", NULL)
)
with_null_device({
  IWTimage(fake_iwt_aov_low, alpha = 0.05, abscissa_range = c(0, 1))
})
expect_true(TRUE)

# ===========================================================================
# ITPimage — deprecated wrapper (just calls IWTimage with IWT1 result)
# ===========================================================================
with_null_device({
  suppressWarnings(ITPimage(res_iwt1, abscissa_range = c(0, 1)))
})
expect_true(TRUE)
