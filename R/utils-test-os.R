os_to_coeffs <- function(data, mu, dx = NULL) {
  if (fda::is.fd(data)) {
    # data is a functional data object
    rangeval <- data$basis$rangeval
    if (is.null(dx)) {
      dx <- (rangeval[2] - rangeval[1]) * 0.01
    }
    abscissa <- seq(rangeval[1], rangeval[2], by = dx)
    coeff <- t(fda::eval.fd(fdobj = data, evalarg = abscissa))
  } else if (is.matrix(data)) {
    coeff <- data
  } else {
    cli::cli_abort(
      "The {.arg data} argument must be either a functional data object of class
      {.cls fd} or a matrix."
    )
  }

  if (fda::is.fd(mu)) {
    # mu is a functional data
    rangeval_mu <- mu$basis$rangeval
    if (sum(rangeval_mu == rangeval) != 2) {
      cli::cli_abort(
        "The range of values of {.arg mu} must be the same as the range of
        values of {.arg data}."
      )
    }
    abscissa <- seq(rangeval_mu[1], rangeval_mu[2], by = dx)
    mu_eval <- t(fda::eval.fd(fdobj = mu, evalarg = abscissa))
  } else if (is.vector(mu)) {
    mu_eval <- mu
  } else {
    cli::cli_abort(
      "The {.arg mu} argument must be either a functional data object of class
      {.cls fd} or a numeric vector."
    )
  }

  list(coeff = coeff, mu = mu_eval)
}

# Internal helper: compute the pointwise test statistic for one-sample test.
# Returns a numeric vector of length p (squared absolute column means).
os_pointwise_stat <- function(coeff) {
  abs(colMeans(coeff))^2
}

# Internal helper: one permutation iteration for os_permtest.
# Returns a numeric vector of length p (the permuted test statistic row).
os_single_perm <- function(coeff, n) {
  signs <- stats::rbinom(n, 1, 0.5) * 2 - 1
  os_pointwise_stat(coeff * signs)
}

# Internal helper: shared pointwise permutation test for one-sample functions.
# Returns list(t0, t_coeff, pval).
os_permtest <- function(coeff, n_perm) {
  n <- nrow(coeff)
  p <- ncol(coeff)

  t0 <- os_pointwise_stat(coeff)

  # Run permutations in parallel via mirai_map().
  # Each task returns one row of t_coeff (length p).
  perm_args <- list(coeff = coeff, n = n)

  if (mirai::daemons_set()) {
    perm_tasks <- mirai::mirai_map(seq_len(n_perm), function(.x) {
      rlang::inject(os_single_perm(!!!perm_args))
    })
    perm_results <- perm_tasks[.progress]
  } else {
    perm_results <- lapply(seq_len(n_perm), function(.x) {
      rlang::inject(os_single_perm(!!!perm_args))
    })
  }
  t_coeff <- do.call(rbind, perm_results)

  pval <- colSums(
    t_coeff >= matrix(t0, nrow = n_perm, ncol = p, byrow = TRUE)
  ) /
    n_perm

  list(t0 = t0, t_coeff = t_coeff, pval = pval)
}

# Internal helper: data preparation + pointwise permutation test for
# one-sample functions. Converts data to coefficients, centres by mu, runs
# os_permtest(), and attaches metadata. Returns a list with fields t0,
# t_coeff, pval, data, mu, p.
os_prepare_data <- function(data, mu, dx, n_perm) {
  inputs <- os_to_coeffs(data, mu, dx = dx)
  coeff <- inputs$coeff
  mu_eval <- inputs$mu

  n <- dim(coeff)[1]
  p <- dim(coeff)[2]

  data_eval <- coeff
  coeff <- coeff - matrix(data = mu_eval, nrow = n, ncol = p, byrow = TRUE)

  out <- os_permtest(coeff, n_perm)
  out$data <- data_eval
  out$mu <- mu_eval
  out$p <- p
  out
}
