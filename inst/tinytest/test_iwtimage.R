# Tests for IWTimage (R/plot-iwt-image.R) and deprecated ITPimage (R/ITPimage.R)
#
# IWTimage dispatches on class:
#   * "fos"    → new one-sample IWT API (iwt1 / functional_one_sample_test)
#   * "fts"    → new two-sample IWT API (iwt2 / functional_two_sample_test)
#   * "IWT1"   → legacy one-sample class (IWT1() backward-compat wrapper)
#   * "IWTaov" → legacy ANOVA class (IWTaov() now returns "faov")
#
# IWT2 is NOT a live class: IWT2() returns class "fts", not "IWT2".
# There is therefore no test for a hypothetical "IWT2" object.
library(ggplot2)
library(tinysnapshot)
using(tinysnapshot)

options(tinysnapshot_device = "png")
options(tinysnapshot_os = "Darwin")

d1 <- NASAtemp$milan[1:4, 1:8]
d2 <- NASAtemp$paris[1:4, 1:8]
p <- ncol(d1)

# ===========================================================================
# IWTimage — fos class (new one-sample API)
# ===========================================================================
set.seed(42)
res_fos <- functional_one_sample_test(d1, mu = 0, n_perm = 5L)

# Default: significance ribbon drawn (all adjusted p-values == 0)
p_fos <- IWTimage(res_fos, abscissa_range = c(0, 1))
expect_true(inherits(p_fos, "patchwork"))
expect_snapshot_plot(p_fos, "iwtimage_fos")

# plot_unadjusted = TRUE exercises the dashed unadjusted-line layer
p_fos_unadj <- IWTimage(
  res_fos,
  abscissa_range = c(0, 1),
  plot_unadjusted = TRUE
)
expect_true(inherits(p_fos_unadj, "patchwork"))

# Vector mu exercises the `else mu` branch in .iwt_image_panels()
set.seed(42)
res_fos_vmu <- functional_one_sample_test(
  d1,
  mu = seq(0, 7, length.out = p),
  n_perm = 5L
)
p_fos_vmu <- IWTimage(res_fos_vmu, abscissa_range = c(0, 1))
expect_true(inherits(p_fos_vmu, "patchwork"))

# fos without pvalue_matrix must error (correction other than "IWT")
fake_fos_no_matrix <- structure(
  list(
    data = d1,
    mu = 0,
    adjusted_pvalues = rep(0.3, p),
    unadjusted_pvalues = rep(0.3, p),
    correction_method = "BH"
  ),
  class = "fos"
)
expect_error(IWTimage(fake_fos_no_matrix), "IWT")

# ===========================================================================
# IWTimage — fts class (new two-sample API)
# ===========================================================================
set.seed(42)
res_fts <- iwt2(d1, d2, n_perm = 5L)

# All adjusted p-values > 0.05: no significance ribbon (covers that branch)
p_fts <- IWTimage(res_fts, abscissa_range = c(0, 1))
expect_true(inherits(p_fts, "patchwork"))
expect_snapshot_plot(p_fts, "iwtimage_fts")

# fts without pvalue_matrix must error (e.g. functional_two_sample_test + TWT)
fake_fts_no_matrix <- structure(
  list(
    data = rbind(d1, d2),
    mu = rep(0, p),
    group_labels = c(rep(1L, nrow(d1)), rep(2L, nrow(d2))),
    adjusted_pvalues = rep(0.3, p),
    unadjusted_pvalues = rep(0.3, p),
    correction_method = "TWT"
  ),
  class = "fts"
)
expect_error(IWTimage(fake_fts_no_matrix), "IWT")

# ===========================================================================
# IWTimage — unknown class must error
# ===========================================================================
expect_error(IWTimage(structure(list(), class = "mystery")), "Unsupported")

# ===========================================================================
# IWTimage — IWT1 class (legacy one-sample wrapper)
# ===========================================================================
set.seed(42)
res_iwt1 <- IWT1(d1, mu = 0, B = 5L)

p_iwt1 <- IWTimage(res_iwt1, abscissa_range = c(0, 1))
expect_true(inherits(p_iwt1, "patchwork"))
expect_snapshot_plot(p_iwt1, "iwtimage_iwt1")

# ===========================================================================
# ITPimage — deprecated wrapper (forwards to IWTimage with an IWT1 object)
# ===========================================================================
p_itp <- suppressWarnings(ITPimage(res_iwt1, abscissa_range = c(0, 1)))
expect_true(inherits(p_itp, "patchwork"))

# ===========================================================================
# IWTimage — IWTaov: one factor, non-interaction, single design-matrix column
# (baseline aov snapshot; non-interaction + single-col paths)
# ===========================================================================
nvar <- 1L
n_obs <- nrow(d1) + nrow(d2)

fake_iwt_aov <- structure(
  list(
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
      `(Intercept)` = rep(1L, n_obs),
      groups = c(rep(0L, nrow(d1)), rep(1L, nrow(d2)))
    )
  ),
  class = "IWTaov"
)

p_aov <- IWTimage(fake_iwt_aov, abscissa_range = c(0, 1))
expect_true(inherits(p_aov, "patchwork"))
expect_snapshot_plot(p_aov, "iwtimage_aov")

# ===========================================================================
# IWTimage — IWTaov: interaction + multi-column design-matrix match
#
# This single fixture covers all four branch combinations in .iwt_image_aov():
#   non-interaction + single-col  (factors "g1", "g2")
#   interaction + multi-col       (factor "g1:g2")
#
# For "g1:g2": grep("g1", colnames) ∩ grep("g2", colnames) gives two columns
# {"g1:g2a", "g1:g2b"} → dim(colors) not NULL → apply(paste) branch taken.
# ===========================================================================
nvar_mc <- 3L
dm_mc <- cbind(
  `(Intercept)` = rep(1L, n_obs),
  g1 = c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L),
  g2 = c(0L, 0L, 1L, 1L, 0L, 0L, 1L, 1L),
  `g1:g2a` = c(0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L),
  `g1:g2b` = c(0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L)
)
fake_iwt_aov_mc <- structure(
  list(
    adjusted_pval_F = rep(0.3, p),
    unadjusted_pval_F = rep(0.4, p),
    pval_matrix_F = matrix(0.3, nrow = p, ncol = p),
    adjusted_pval_factor = matrix(
      0.3,
      nrow = nvar_mc,
      ncol = p,
      dimnames = list(c("g1", "g2", "g1:g2"), NULL)
    ),
    unadjusted_pval_factor = matrix(
      0.4,
      nrow = nvar_mc,
      ncol = p,
      dimnames = list(c("g1", "g2", "g1:g2"), NULL)
    ),
    pval_matrix_factor = array(0.3, dim = c(nvar_mc, p, p)),
    data_eval = rbind(d1, d2),
    design_matrix = dm_mc
  ),
  class = "IWTaov"
)

p_mc <- IWTimage(fake_iwt_aov_mc, abscissa_range = c(0, 1))
expect_true(inherits(p_mc, "patchwork"))

set.seed(NULL)
