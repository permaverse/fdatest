# Tests for IWTimage (R/plot-iwt-image.R) and deprecated ITPimage (R/ITPimage.R)
#
# IWTimage dispatches on class:
#   * "fos"            → new one-sample IWT API (iwt1 / functional_one_sample_test)
#   * "fts"            → new two-sample IWT API (iwt2 / functional_two_sample_test)
#   * "IWT1" or "IWT2" → legacy classes kept for backward compatibility
#   * "IWTaov"         → legacy ANOVA class (IWTaov() now returns "faov")
#
# To achieve 100% line coverage we create minimal fake objects with the old
# class names where the current public API can no longer produce them.
library(ggplot2)
library(tinysnapshot)
using(tinysnapshot)

options(tinysnapshot_device = "png")
options(tinysnapshot_os = "Darwin")

d1 <- NASAtemp$milan[1:4, 1:8]
d2 <- NASAtemp$paris[1:4, 1:8]
p <- ncol(d1)

# ===========================================================================
# Shared legacy IWT1 result (exercises the IWT1 branch and is reused below)
# ===========================================================================
set.seed(42)
res_iwt1 <- IWT1(d1, mu = 0, B = 5L)

# ===========================================================================
# IWTimage — fos class (new one-sample API)
# ===========================================================================
set.seed(42)
res_fos <- functional_one_sample_test(d1, mu = 0, n_perm = 5L)
expect_true(inherits(res_fos, "fos"))

p_fos <- IWTimage(res_fos, abscissa_range = c(0, 1))
expect_true(inherits(p_fos, "patchwork"))
expect_snapshot_plot(p_fos, "iwtimage_fos")

p_fos_unadj <- IWTimage(
  res_fos,
  abscissa_range = c(0, 1),
  plot_unadjusted = TRUE
)
expect_true(inherits(p_fos_unadj, "patchwork"))

# fos without pvalue_matrix must error
fake_fos_no_matrix <- structure(
  list(
    data = d1,
    mu = rep(0, p),
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
expect_true(inherits(res_fts, "fts"))

p_fts <- IWTimage(res_fts, abscissa_range = c(0, 1))
expect_true(inherits(p_fts, "patchwork"))
expect_snapshot_plot(p_fts, "iwtimage_fts")

# fts without pvalue_matrix must error
fake_fts_no_matrix <- structure(
  list(
    data = list(data1 = d1, data2 = d2),
    mu = rep(0, p),
    group_labels = c(rep("1", 4L), rep("2", 4L)),
    adjusted_pvalues = rep(0.3, p),
    unadjusted_pvalues = rep(0.3, p),
    correction_method = "BH"
  ),
  class = "fts"
)
expect_error(IWTimage(fake_fts_no_matrix), "IWT")

# ===========================================================================
# IWTimage — unknown class must error
# ===========================================================================
expect_error(IWTimage(structure(list(), class = "mystery")), "Unsupported")

# ===========================================================================
# IWTimage — IWT1 class (legacy one-sample, default settings)
p_iwt1 <- IWTimage(res_iwt1, abscissa_range = c(0, 1))
expect_true(inherits(p_iwt1, "patchwork"))
expect_snapshot_plot(p_iwt1, "iwtimage_iwt1")

# plot_unadjusted = TRUE
p_iwt1_unadj <- IWTimage(
  res_iwt1,
  abscissa_range = c(0, 1),
  plot_unadjusted = TRUE
)
expect_true(inherits(p_iwt1_unadj, "patchwork"))

# Force all adjusted p-values below alpha (exercises the significance ribbon)
res_iwt1_low <- res_iwt1
res_iwt1_low$adjusted_pval <- rep(0.01, p)
p_iwt1_low <- IWTimage(res_iwt1_low, alpha = 0.05, abscissa_range = c(0, 1))
expect_true(inherits(p_iwt1_low, "patchwork"))

# ===========================================================================
# IWTimage — fake "IWT2" object (exercises the 2-pop sub-branch of IWT1/IWT2)
#
# IWT2() now sets class "fts", not "IWT2", but the legacy branch still exists.
# ===========================================================================
fake_iwt2 <- structure(
  list(
    test = "2pop",
    mu = 0,
    adjusted_pval = rep(0.3, p),
    unadjusted_pval = rep(0.3, p),
    pval_matrix = matrix(0.3, nrow = p, ncol = p),
    data_eval = rbind(d1, d2),
    ord_labels = c(rep(1L, nrow(d1)), rep(2L, nrow(d2)))
  ),
  class = "IWT2"
)

p_iwt2 <- IWTimage(fake_iwt2, abscissa_range = c(0, 1))
expect_true(inherits(p_iwt2, "patchwork"))

p_iwt2_unadj <- IWTimage(
  fake_iwt2,
  abscissa_range = c(0, 1),
  plot_unadjusted = TRUE
)
expect_true(inherits(p_iwt2_unadj, "patchwork"))

# ===========================================================================
# IWTimage — fake "IWTaov" object (one factor, no interactions)
#
# IWTaov() now returns class "faov", but the legacy IWTaov branch still exists.
# ===========================================================================
nvar <- 1L

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
      `(Intercept)` = rep(1L, nrow(d1) + nrow(d2)),
      groups = c(rep(0L, nrow(d1)), rep(1L, nrow(d2)))
    )
  ),
  class = "IWTaov"
)

p_aov <- IWTimage(fake_iwt_aov, abscissa_range = c(0, 1))
expect_true(inherits(p_aov, "patchwork"))
expect_snapshot_plot(p_aov, "iwtimage_aov")

p_aov_unadj <- IWTimage(
  fake_iwt_aov,
  abscissa_range = c(0, 1),
  plot_unadjusted = TRUE
)
expect_true(inherits(p_aov_unadj, "patchwork"))

# Force p-values below alpha (exercises significance ribbon in all panels)
fake_iwt_aov_low <- fake_iwt_aov
fake_iwt_aov_low$adjusted_pval_F <- rep(0.01, p)
fake_iwt_aov_low$adjusted_pval_factor <- matrix(
  0.01,
  nrow = nvar,
  ncol = p,
  dimnames = list("groups", NULL)
)
p_aov_low <- IWTimage(fake_iwt_aov_low, alpha = 0.05, abscissa_range = c(0, 1))
expect_true(inherits(p_aov_low, "patchwork"))

# ===========================================================================
# ITPimage — deprecated wrapper (just calls IWTimage with an IWT1 result)
# ===========================================================================
p_itp <- suppressWarnings(ITPimage(res_iwt1, abscissa_range = c(0, 1)))
expect_true(inherits(p_itp, "patchwork"))

# ===========================================================================
# IWTimage — fake "IWTaov" with two main factors + interaction
# ===========================================================================
nvar_ia <- 3L
n_obs <- nrow(d1) + nrow(d2)
grp_a <- c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L)
grp_b <- c(0L, 0L, 1L, 1L, 0L, 0L, 1L, 1L)

fake_iwt_aov_ia <- structure(
  list(
    adjusted_pval_F = rep(0.3, p),
    unadjusted_pval_F = rep(0.4, p),
    pval_matrix_F = matrix(0.3, nrow = p, ncol = p),
    adjusted_pval_factor = matrix(
      0.3,
      nrow = nvar_ia,
      ncol = p,
      dimnames = list(c("grpA", "grpB", "grpA:grpB"), NULL)
    ),
    unadjusted_pval_factor = matrix(
      0.4,
      nrow = nvar_ia,
      ncol = p,
      dimnames = list(c("grpA", "grpB", "grpA:grpB"), NULL)
    ),
    pval_matrix_factor = array(0.3, dim = c(nvar_ia, p, p)),
    data_eval = rbind(d1, d2),
    design_matrix = cbind(
      `(Intercept)` = rep(1L, n_obs),
      grpA = grp_a,
      grpB = grp_b,
      `grpA:grpB` = grp_a * grp_b
    )
  ),
  class = "IWTaov"
)

p_ia <- IWTimage(fake_iwt_aov_ia, abscissa_range = c(0, 1))
expect_true(inherits(p_ia, "patchwork"))

p_ia_unadj <- IWTimage(
  fake_iwt_aov_ia,
  abscissa_range = c(0, 1),
  plot_unadjusted = TRUE
)
expect_true(inherits(p_ia_unadj, "patchwork"))

# Significance ribbons in the interaction case
fake_iwt_aov_ia_low <- fake_iwt_aov_ia
fake_iwt_aov_ia_low$adjusted_pval_F <- rep(0.01, p)
fake_iwt_aov_ia_low$adjusted_pval_factor <- matrix(
  0.01,
  nrow = nvar_ia,
  ncol = p,
  dimnames = list(c("grpA", "grpB", "grpA:grpB"), NULL)
)
p_ia_low <- IWTimage(
  fake_iwt_aov_ia_low,
  alpha = 0.05,
  abscissa_range = c(0, 1)
)
expect_true(inherits(p_ia_low, "patchwork"))

# ===========================================================================
# IWTimage — interaction factor matching multiple design-matrix columns
#
# intersect(grep("g1", colnames), grep("g2", colnames)) = {g1:g2a, g1:g2b}
# → 2 columns → matrix → apply(paste) branch
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

# ===========================================================================
# IWTimage — non-interaction factor matching multiple design-matrix columns
#
# grep("grp", c("(Intercept)", "grpA", "grpB")) → {grpA, grpB}
# → 2 columns → apply(paste) branch
# ===========================================================================
nvar_mf <- 1L
dm_mf <- cbind(
  `(Intercept)` = rep(1L, n_obs),
  grpA = c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L),
  grpB = c(0L, 0L, 1L, 1L, 0L, 0L, 1L, 1L)
)
fake_iwt_aov_mf <- structure(
  list(
    adjusted_pval_F = rep(0.3, p),
    unadjusted_pval_F = rep(0.4, p),
    pval_matrix_F = matrix(0.3, nrow = p, ncol = p),
    adjusted_pval_factor = matrix(
      0.3,
      nrow = nvar_mf,
      ncol = p,
      dimnames = list("grp", NULL)
    ),
    unadjusted_pval_factor = matrix(
      0.4,
      nrow = nvar_mf,
      ncol = p,
      dimnames = list("grp", NULL)
    ),
    pval_matrix_factor = array(0.3, dim = c(nvar_mf, p, p)),
    data_eval = rbind(d1, d2),
    design_matrix = dm_mf
  ),
  class = "IWTaov"
)

p_mf <- IWTimage(fake_iwt_aov_mf, abscissa_range = c(0, 1))
expect_true(inherits(p_mf, "patchwork"))

set.seed(NULL)
