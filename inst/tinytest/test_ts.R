# Tests for two-sample tests:
# two-iwt.R, two-twt.R, two-fdr.R, two-pct.R, two-global.R
# functional-two-sample-test.R
# ITP2bspline.R, ITP2fourier.R, ITP2pafourier.R

d1 <- NASAtemp$milan[1:4, 1:8]
d2 <- NASAtemp$paris[1:4, 1:8]
p <- ncol(d1)
partition <- c(rep(1L, 4L), rep(2L, 4L)) # 2 partitions over 8 grid pts

# ===========================================================================
# IWT2
# ===========================================================================
set.seed(42)
res_iwt <- IWT2(data1 = d1, data2 = d2, mu = 0, B = 5L, verbose = FALSE)

expect_inherits(res_iwt, "fts")
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
expect_inherits(res_iwt_v, "fts")

# recycle disabled
set.seed(42)
res_iwt_nr <- IWT2(
  data1 = d1,
  data2 = d2,
  mu = 0,
  B = 5L,
  recycle = FALSE,
  verbose = FALSE
)
expect_inherits(res_iwt_nr, "fts")
expect_true(is.na(res_iwt_nr$pvalue_matrix[1, p])) # upper-right is NA

# recycle = FALSE + verbose = TRUE (exercises the cli progress line in the else branch)
set.seed(42)
res_iwt_nr_v <- IWT2(
  data1 = d1,
  data2 = d2,
  mu = 0,
  B = 5L,
  recycle = FALSE,
  verbose = TRUE
)
expect_inherits(res_iwt_nr_v, "fts")

# paired samples
set.seed(42)
res_iwt_p <- IWT2(
  data1 = d1,
  data2 = d2,
  mu = 0,
  B = 5L,
  paired = TRUE,
  verbose = FALSE
)
expect_inherits(res_iwt_p, "fts")

# aggregation_strategy = "max" (tested via new API)
set.seed(42)
res_iwt_max <- iwt2(
  data1 = d1,
  data2 = d2,
  mu = 0,
  n_perm = 5L,
  aggregation_strategy = "max"
)
expect_inherits(res_iwt_max, "fts")
expect_equal(length(res_iwt_max$adjusted_pvalues), p)

# ===========================================================================
# global2 — all four combinations of standardize / aggregation_strategy
# ===========================================================================
for (std in c(FALSE, TRUE)) {
  for (agg in c("integral", "max")) {
    set.seed(42)
    res_g <- global2(
      data1 = d1,
      data2 = d2,
      mu = 0,
      n_perm = 5L,
      standardize = std,
      aggregation_strategy = agg
    )
    info_label <- paste0(
      "global2 standardize=",
      std,
      " aggregation_strategy=",
      agg
    )
    expect_inherits(res_g, "fts", info = info_label)
    expect_true(
      !is.null(res_g$global_pvalue),
      info = info_label
    )
    expect_true(
      res_g$global_pvalue >= 0 && res_g$global_pvalue <= 1,
      info = info_label
    )
  }
}

# Global2 paired
set.seed(42)
res_g_p <- Global2(
  data1 = d1,
  data2 = d2,
  mu = 0,
  B = 5L,
  paired = TRUE
)
expect_inherits(res_g_p, "fts")

# Global2 verbose = TRUE
set.seed(42)
res_g_v <- Global2(data1 = d1, data2 = d2, mu = 0, B = 5L, verbose = TRUE)
expect_inherits(res_g_v, "fts")

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
  expect_inherits(res_t, "fts", info = paste("TWT2 alternative =", alt))
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
expect_inherits(res_t_p, "fts")

# TWT2 verbose = TRUE
set.seed(42)
res_t_v <- TWT2(data1 = d1, data2 = d2, mu = 0, B = 5L, verbose = TRUE)
expect_inherits(res_t_v, "fts")

# aggregation_strategy = "max" (tested via new API)
set.seed(42)
res_t_max <- twt2(
  data1 = d1,
  data2 = d2,
  mu = 0,
  n_perm = 5L,
  aggregation_strategy = "max"
)
expect_inherits(res_t_max, "fts")
expect_equal(length(res_t_max$adjusted_pvalues), p)

# ===========================================================================
# FDR2 — two.sided + alternative options
# ===========================================================================
set.seed(42)
res_fdr <- FDR2(data1 = d1, data2 = d2, mu = 0, B = 5L)
expect_inherits(res_fdr, "fts")
expect_equal(length(res_fdr$adjusted_pvalues), p)

set.seed(42)
res_fdr_less <- FDR2(
  data1 = d1,
  data2 = d2,
  mu = 0,
  B = 5L,
  alternative = "less"
)
expect_inherits(res_fdr_less, "fts")

set.seed(42)
res_fdr_gr <- FDR2(
  data1 = d1,
  data2 = d2,
  mu = 0,
  B = 5L,
  alternative = "greater"
)
expect_inherits(res_fdr_gr, "fts")

# FDR2 paired
set.seed(42)
res_fdr_p <- FDR2(data1 = d1, data2 = d2, mu = 0, B = 5L, paired = TRUE)
expect_inherits(res_fdr_p, "fts")

# FDR2 verbose = TRUE
set.seed(42)
res_fdr_v <- FDR2(data1 = d1, data2 = d2, mu = 0, B = 5L, verbose = TRUE)
expect_inherits(res_fdr_v, "fts")

# ===========================================================================
# PCT2
# ===========================================================================
set.seed(42)
res_pct <- PCT2(data1 = d1, data2 = d2, mu = 0, B = 5L, partition = partition)
expect_inherits(res_pct, "fts")
expect_equal(length(res_pct$adjusted_pvalues), p)

# PCT2 verbose = TRUE
set.seed(42)
res_pct_v <- PCT2(
  data1 = d1,
  data2 = d2,
  mu = 0,
  B = 5L,
  partition = partition,
  verbose = TRUE
)
expect_inherits(res_pct_v, "fts")

# aggregation_strategy = "max" (tested via new API)
set.seed(42)
res_pct_max <- pct2(
  data1 = d1,
  data2 = d2,
  mu = 0,
  n_perm = 5L,
  partition = partition,
  aggregation_strategy = "max"
)
expect_inherits(res_pct_max, "fts")
expect_equal(length(res_pct_max$adjusted_pvalues), p)

# ===========================================================================
# functional_two_sample_test — all corrections
# ===========================================================================
set.seed(42)
res_ft_iwt <- functional_two_sample_test(
  data1 = d1,
  data2 = d2,
  correction = "IWT",
  n_perm = 5L,
  verbose = FALSE
)
expect_inherits(res_ft_iwt, "fts")
expect_equal(res_ft_iwt$correction, "IWT")

set.seed(42)
res_ft_twt <- functional_two_sample_test(
  data1 = d1,
  data2 = d2,
  correction = "TWT",
  n_perm = 5L,
  verbose = FALSE
)
expect_inherits(res_ft_twt, "fts")
expect_equal(res_ft_twt$correction, "TWT")

set.seed(42)
res_ft_g <- functional_two_sample_test(
  data1 = d1,
  data2 = d2,
  correction = "Global",
  n_perm = 5L
)
expect_inherits(res_ft_g, "fts")
expect_equal(res_ft_g$correction, "Global")

set.seed(42)
res_ft_fdr <- functional_two_sample_test(
  data1 = d1,
  data2 = d2,
  correction = "FDR",
  n_perm = 5L
)
expect_inherits(res_ft_fdr, "fts")
expect_equal(res_ft_fdr$correction, "FDR")

# PCT without partition → error
expect_error(
  functional_two_sample_test(
    data1 = d1,
    data2 = d2,
    correction = "PCT",
    n_perm = 5L
  )
)

# PCT with partition
set.seed(42)
res_ft_pct <- functional_two_sample_test(
  data1 = d1,
  data2 = d2,
  correction = "PCT",
  n_perm = 5L,
  partition = partition
)
expect_inherits(res_ft_pct, "fts")
expect_equal(res_ft_pct$correction, "PCT")

# ===========================================================================
# alternative != "two.sided" combined with standardize = TRUE
# (exercises the one-sided + standardised path in ts_prepare_data, shared
#  across all methods; iwt2 used as a representative)
# ===========================================================================
for (alt in c("less", "greater")) {
  set.seed(42)
  res_std <- iwt2(
    data1 = d1,
    data2 = d2,
    mu = 0,
    n_perm = 5L,
    alternative = alt,
    standardize = TRUE
  )
  expect_inherits(
    res_std,
    "fts",
    info = paste("iwt2 alternative =", alt, "standardize = TRUE")
  )
  expect_equal(
    length(res_std$adjusted_pvalues),
    p,
    info = paste("iwt2 alternative =", alt, "standardize = TRUE")
  )
}

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
expect_inherits(res_itp2b, "fts")

set.seed(42)
res_itp2f <- suppressWarnings(ITP2fourier(
  data1 = d1,
  data2 = d2,
  mu = 0,
  B = 5L
))
expect_inherits(res_itp2f, "fts")

set.seed(42)
res_itp2pa <- suppressWarnings(ITP2pafourier(
  data1 = d1,
  data2 = d2,
  mu = 0,
  B = 5L
))
expect_inherits(res_itp2pa, "fts")

set.seed(NULL)
