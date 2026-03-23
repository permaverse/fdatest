# Tests for two-sample tests:
# R/two-iwt.R, R/two-twt.R, R/two-fdr.R, R/two-pct.R, R/two-global.R,
# R/functional-two-sample-test.R
# R/ITP2bspline.R, R/ITP2fourier.R, R/ITP2pafourier.R
library(fdatest)

data("NASAtemp", package = "fdatest")
d1 <- NASAtemp$milan[1:4, 1:8]
d2 <- NASAtemp$paris[1:4, 1:8]
p <- ncol(d1)
partition <- c(rep(1L, 4L), rep(2L, 4L)) # 2 partitions over 8 grid pts

# ===========================================================================
# IWT2
# ===========================================================================
set.seed(42)
res_iwt <- IWT2(data1 = d1, data2 = d2, mu = 0, B = 5L, verbose = FALSE)

expect_inherits(res_iwt, "ftwosample")
expect_equal(length(res_iwt$adjusted_pvalues), p)
expect_equal(length(res_iwt$unadjusted_pvalues), p)
expect_true(all(res_iwt$adjusted_pvalues >= 0 & res_iwt$adjusted_pvalues <= 1))
expect_true(all(
  res_iwt$unadjusted_pvalues >= 0 & res_iwt$unadjusted_pvalues <= 1
))
expect_equal(dim(res_iwt$pvalue_matrix), c(p, p))
expect_true(all(res_iwt$adjusted_pvalues >= res_iwt$unadjusted_pvalues))

# verbose = TRUE (exercises CLI progress path)
set.seed(42)
res_iwt_v <- IWT2(data1 = d1, data2 = d2, mu = 0, B = 5L, verbose = TRUE)
expect_inherits(res_iwt_v, "ftwosample")

# recycle = FALSE
set.seed(42)
res_iwt_nr <- IWT2(
  data1 = d1,
  data2 = d2,
  mu = 0,
  B = 5L,
  recycle = FALSE,
  verbose = FALSE
)
expect_inherits(res_iwt_nr, "ftwosample")
expect_true(is.na(res_iwt_nr$pvalue_matrix[1, p])) # upper-right is NA

# paired = TRUE
set.seed(42)
res_iwt_p <- IWT2(
  data1 = d1,
  data2 = d2,
  mu = 0,
  B = 5L,
  paired = TRUE,
  verbose = FALSE
)
expect_inherits(res_iwt_p, "ftwosample")

# ===========================================================================
# Global2 — all four statistic options
# ===========================================================================
for (stat in c("Integral", "Max", "Integral_std", "Max_std")) {
  set.seed(42)
  res_g <- Global2(data1 = d1, data2 = d2, mu = 0, B = 5L, statistic = stat)
  expect_inherits(
    res_g,
    "ftwosample",
    info = paste("Global2 statistic =", stat)
  )
  expect_true(
    !is.null(res_g$global_pvalue),
    info = paste("Global2 global_pvalue for", stat)
  )
  expect_true(
    res_g$global_pvalue >= 0 & res_g$global_pvalue <= 1,
    info = paste("Global2 pvalue range for", stat)
  )
}

# Global2 paired
set.seed(42)
res_g_p <- Global2(
  data1 = d1,
  data2 = d2,
  mu = 0,
  B = 5L,
  statistic = "Integral",
  paired = TRUE
)
expect_inherits(res_g_p, "ftwosample")

# ===========================================================================
# TWT2 — all three alternative options
# ===========================================================================
for (alt in c("two.sided", "less", "greater")) {
  set.seed(42)
  res_t <- TWT2(
    data1 = d1,
    data2 = d2,
    mu = 0,
    B = 5L,
    alternative = alt,
    verbose = FALSE
  )
  expect_inherits(res_t, "ftwosample", info = paste("TWT2 alternative =", alt))
  expect_equal(
    length(res_t$adjusted_pvalues),
    p,
    info = paste("TWT2 length for", alt)
  )
}

# TWT2 paired
set.seed(42)
res_t_p <- TWT2(
  data1 = d1,
  data2 = d2,
  mu = 0,
  B = 5L,
  paired = TRUE,
  verbose = FALSE
)
expect_inherits(res_t_p, "ftwosample")

# TWT2 verbose = TRUE
set.seed(42)
res_t_v <- TWT2(data1 = d1, data2 = d2, mu = 0, B = 5L, verbose = TRUE)
expect_inherits(res_t_v, "ftwosample")

# ===========================================================================
# FDR2 — two.sided + alternative options
# ===========================================================================
set.seed(42)
res_fdr <- FDR2(data1 = d1, data2 = d2, mu = 0, B = 5L)
expect_inherits(res_fdr, "ftwosample")
expect_equal(length(res_fdr$adjusted_pvalues), p)

set.seed(42)
res_fdr_less <- FDR2(
  data1 = d1,
  data2 = d2,
  mu = 0,
  B = 5L,
  alternative = "less"
)
expect_inherits(res_fdr_less, "ftwosample")

set.seed(42)
res_fdr_gr <- FDR2(
  data1 = d1,
  data2 = d2,
  mu = 0,
  B = 5L,
  alternative = "greater"
)
expect_inherits(res_fdr_gr, "ftwosample")

# FDR2 paired
set.seed(42)
res_fdr_p <- FDR2(data1 = d1, data2 = d2, mu = 0, B = 5L, paired = TRUE)
expect_inherits(res_fdr_p, "ftwosample")

# ===========================================================================
# PCT2
# ===========================================================================
set.seed(42)
res_pct <- PCT2(data1 = d1, data2 = d2, mu = 0, B = 5L, partition = partition)
expect_inherits(res_pct, "ftwosample")
expect_equal(length(res_pct$adjusted_pvalues), p)

# ===========================================================================
# functional_two_sample_test — all corrections
# ===========================================================================
set.seed(42)
res_ft_iwt <- functional_two_sample_test(
  data1 = d1,
  data2 = d2,
  correction = "IWT",
  B = 5L,
  verbose = FALSE
)
expect_inherits(res_ft_iwt, "ftwosample")
expect_equal(res_ft_iwt$correction, "IWT")

set.seed(42)
res_ft_twt <- functional_two_sample_test(
  data1 = d1,
  data2 = d2,
  correction = "TWT",
  B = 5L,
  verbose = FALSE
)
expect_inherits(res_ft_twt, "ftwosample")
expect_equal(res_ft_twt$correction, "TWT")

set.seed(42)
res_ft_g <- functional_two_sample_test(
  data1 = d1,
  data2 = d2,
  correction = "Global",
  B = 5L
)
expect_inherits(res_ft_g, "ftwosample")
expect_equal(res_ft_g$correction, "Global")

set.seed(42)
res_ft_fdr <- functional_two_sample_test(
  data1 = d1,
  data2 = d2,
  correction = "FDR",
  B = 5L
)
expect_inherits(res_ft_fdr, "ftwosample")
expect_equal(res_ft_fdr$correction, "FDR")

# PCT without partition → error
expect_error(
  functional_two_sample_test(
    data1 = d1,
    data2 = d2,
    correction = "PCT",
    B = 5L
  )
)

# PCT with partition
set.seed(42)
res_ft_pct <- functional_two_sample_test(
  data1 = d1,
  data2 = d2,
  correction = "PCT",
  B = 5L,
  partition = partition
)
expect_inherits(res_ft_pct, "ftwosample")
expect_equal(res_ft_pct$correction, "PCT")

# ===========================================================================
# Deprecated wrappers
# ===========================================================================
set.seed(42)
res_itp2b <- suppressWarnings(ITP2bspline(
  data1 = d1,
  data2 = d2,
  mu = 0,
  B = 5L
))
expect_inherits(res_itp2b, "ftwosample")

set.seed(42)
res_itp2f <- suppressWarnings(ITP2fourier(
  data1 = d1,
  data2 = d2,
  mu = 0,
  B = 5L
))
expect_inherits(res_itp2f, "ftwosample")

set.seed(42)
res_itp2pa <- suppressWarnings(ITP2pafourier(
  data1 = d1,
  data2 = d2,
  mu = 0,
  B = 5L
))
expect_inherits(res_itp2pa, "ftwosample")
