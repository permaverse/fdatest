ts_pointwise_stat <- function(coeff, n1, alternative, standardize) {
  n <- nrow(coeff)

  meandiff <- colMeans(coeff[seq_len(n1), , drop = FALSE], na.rm = TRUE) -
    colMeans(coeff[seq(n1 + 1L, n), , drop = FALSE], na.rm = TRUE)

  if (standardize) {
    s1 <- stats::cov(coeff[seq_len(n1), , drop = FALSE])
    s2 <- stats::cov(coeff[seq(n1 + 1L, n), , drop = FALSE])
    sp <- ((n1 - 1L) * s1 + (n - n1 - 1L) * s2) / (n - 2L)
  }

  if (alternative == "two.sided") {
    t0 <- meandiff^2
    if (standardize) {
      t0 <- t0 / diag(sp)
    }
    return(t0)
  }

  sign_diff <- sign(meandiff)
  sign_diff[sign_diff == -1L] <- 0L
  t0 <- switch(
    alternative,
    greater = (meandiff * sign_diff)^2,
    less = (meandiff * (sign_diff - 1L))^2
  )

  if (standardize) {
    t0 <- t0 / diag(sp)
  }

  t0
}

# Internal helper: one permutation iteration for ts_permtest.
# Returns a numeric vector of length p (the permuted test statistic row).
.ts_one_perm <- function(coeff, n1, alternative, paired, standardize) {
  n <- nrow(coeff)
  if (paired) {
    if_perm <- stats::rbinom(n1, 1, 0.5)
    coeff_perm <- coeff
    for (couple in seq_len(n1)) {
      if (if_perm[couple] == 1L) {
        coeff_perm[c(couple, n1 + couple), ] <- coeff[c(n1 + couple, couple), ]
      }
    }
  } else {
    coeff_perm <- coeff[sample(n), ]
  }

  ts_pointwise_stat(
    coeff = coeff_perm,
    n1 = n1,
    alternative = alternative,
    standardize = standardize
  )
}

# Internal helper: shared pointwise permutation test for two-sample functions
# (iwt2, twt2, fdr2, pct2). Returns list(t0, t_coeff, pval).
ts_permtest <- function(coeff, n1, n_perm, alternative, paired, standardize) {
  p <- ncol(coeff)

  t0 <- ts_pointwise_stat(coeff, n1, alternative, standardize)

  # Run permutations in parallel via mirai_map().
  # Each task returns one row of t_coeff (length p).
  perm_args <- list(
    coeff = coeff,
    n1 = n1,
    alternative = alternative,
    paired = paired,
    standardize = standardize
  )

  if (mirai::daemons_set()) {
    perm_tasks <- mirai::mirai_map(seq_len(n_perm), function(.x) {
      rlang::inject(.ts_one_perm(!!!perm_args))
    })
    perm_results <- perm_tasks[.progress]
  } else {
    perm_results <- lapply(seq_len(n_perm), function(.x) {
      rlang::inject(.ts_one_perm(!!!perm_args))
    })
  }
  t_coeff <- do.call(rbind, perm_results)

  pval <- colSums(
    t_coeff >= matrix(t0, nrow = n_perm, ncol = p, byrow = TRUE)
  ) /
    n_perm

  list(t0 = t0, t_coeff = t_coeff, pval = pval)
}

ts_prepare_data <- function(
  data1,
  data2,
  mu,
  dx,
  n_perm,
  paired,
  alternative,
  standardize
) {
  inputs <- twosamples2coeffs(data1, data2, mu, dx = dx)
  coeff1 <- inputs$coeff1
  coeff2 <- inputs$coeff2
  mu_eval <- inputs$mu

  n1 <- dim(coeff1)[1]
  n2 <- dim(coeff2)[1]
  p <- dim(coeff1)[2]
  if (ncol(coeff2) != p) {
    cli::cli_abort(
      "{.strong 'data1'} and {.strong 'data2'} must have the same number of columns (evaluation points)."
    )
  }

  data_eval <- rbind(coeff1, coeff2)

  # Re-center data by mu_eval for permutation test so that mean of sample 1
  # matches mean of sample 2 under null hypothesis.
  coeff1 <- coeff1 - matrix(data = mu_eval, nrow = n1, ncol = p)
  coeff <- rbind(coeff1, coeff2)
  group_labels <- c(rep(1, n1), rep(2, n2))

  # Pointwise permutation test.
  out <- ts_permtest(coeff, n1, n_perm, alternative, paired, standardize)
  out$data <- data_eval
  out$mu <- mu_eval
  out$group_labels <- group_labels
  out$p <- p
  out
}

compute_row_ts <- function(
  i,
  t0_2x,
  t_coeff_2x,
  n_perm,
  p,
  nvar,
  recycle,
  aggregation_strategy
) {
  js <- if (recycle) seq_len(p) else seq_len(i)
  row_vals <- numeric(if (recycle) p else i)
  for (k in seq_along(js)) {
    j <- js[k]
    inf <- j
    sup <- (p - i) + j
    t0_temp <- if (aggregation_strategy == "integral") {
      sum(t0_2x[inf:sup], na.rm = TRUE)
    } else {
      max(t0_2x[inf:sup], na.rm = TRUE)
    }
    t_temp <- if (aggregation_strategy == "integral") {
      rowSums(t_coeff_2x[, inf:sup, drop = FALSE], na.rm = TRUE)
    } else {
      apply(t_coeff_2x[, inf:sup, drop = FALSE], 1, max, na.rm = TRUE)
    }
    row_vals[k] <- sum(t_temp >= t0_temp) / n_perm
  }
  row_vals
}

compute_row_pair_ts <- function(
  i,
  t0_2x,
  t_coeff_2x,
  n_perm,
  p,
  nvar,
  recycle,
  aggregation_strategy
) {
  if (i == p - i) {
    return(list(compute_row_ts(
      i = i,
      t0_2x = t0_2x,
      t_coeff_2x = t_coeff_2x,
      n_perm = n_perm,
      p = p,
      nvar = nvar,
      recycle = recycle,
      aggregation_strategy = aggregation_strategy
    )))
  }
  list(
    compute_row_ts(
      i = i,
      t0_2x = t0_2x,
      t_coeff_2x = t_coeff_2x,
      n_perm = n_perm,
      p = p,
      nvar = nvar,
      recycle = recycle,
      aggregation_strategy = aggregation_strategy
    ),
    compute_row_ts(
      i = p - i,
      t0_2x = t0_2x,
      t_coeff_2x = t_coeff_2x,
      n_perm = n_perm,
      p = p,
      nvar = nvar,
      recycle = recycle,
      aggregation_strategy = aggregation_strategy
    )
  )
}
