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
