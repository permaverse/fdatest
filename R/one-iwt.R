#' One population Interval Wise Testing procedure
#'
#' The function implements the Interval Wise Testing procedure for testing the
#' center of symmetry of a functional population. Functional data are tested
#' locally and  unadjusted  and adjusted p-value functions are provided. The
#' unadjusted p-value function controls the point-wise error rate. The adjusted
#' p-value function controls the interval-wise error rate.
#'
#' @param data Either pointwise evaluations of the functional data set on a
#'   uniform grid, or an \code{\link[fda]{fd}}. If pointwise evaluations are
#'   provided, `data` is a matrix of dimensions `c(n, J)`, with `J` evaluations
#'   on columns and `n` units on rows.
#' @param mu The center of symmetry under the null hypothesis. Three
#'   possibilities are available for `mu`:
#'
#' - a constant (in this case, a constant function is used);
#' - a `J`-dimensional vector containing the evaluations on the same grid which
#' `data` are evaluated;
#' - a \code{\link[fda]{fd}} object containing one function.
#' Defaults to `0`.
#' @param B The number of iterations of the MC algorithm to evaluate the
#'   p-values of the permutation tests. Defaults to `1000L`.
#' @param dx Used only if an \code{\link[fda]{fd}} object is provided. In this
#'   case, `dx` is the size of the discretization step of the gridused to
#'   evaluate functional data. If set to `NULL`, a grid of size `100L` is used.
#'   Defaults to `NULL`.
#' @param recycle Flag used to decide whether the recycled version of the IWT
#'   should be used (see Pini and Vantini, 2017 for details). Defaults to
#'   `TRUE`.
#'
#' @return An object of class \code{\link{IWT1}}, which is a list containing at
#'   least the following components:
#'
#' - `test`: String vector indicating the type of test performed. In this case
#' equal to `"1pop"`.
#' - `mu`: Evaluation on a grid of the center of symmetry under the null
#' hypothesis (as entered by the user).
#' - `unadjusted_pval`: Evaluation on a grid of the unadjusted p-value function.
#' - `pval_matrix`: Matrix of dimensions `c(p, p)` of the p-values of the
#' multivariate tests. The element `(i,j)` of matrix `pval_matrix` contains the
#' p-value of the joint NPC test of the components `(j,j+1,...,j+(p-i))`.
#' - `adjusted_pval`: Evaluation on a grid of the adjusted p-value function.
#' - `data_eval`: Evaluation on a grid of the functional data.
#'
#' @seealso See also \code{\link{plot.IWT1}} and \code{\link{IWTimage}} for
#'   plotting the results.
#'
#' @references
#' A. Pini and S. Vantini (2017). The Interval Testing Procedure: Inference
#' for Functional Data Controlling the Family Wise Error Rate on Intervals.
#' *Biometrics*, 73(3): 835–845.
#'
#' A. Pini and S. Vantini (2017). Interval-wise testing for functional data.
#' *Journal of Nonparametric Statistics*, 29(2), 407-424.
#'
#' @export
#' @examples
#' # Performing the IWT for one population
#' IWT_result <- IWT1(NASAtemp$paris, mu = 4, B = 10L)
#'
#' # Plotting the results of the IWT
#' plot(IWT_result, xrange = c(0, 12), main = 'Paris temperatures')
#'
#' # Plotting the p-value heatmap
#' IWTimage(IWT_result, abscissa_range = c(0, 12))
#'
#' # Selecting the significant components at 5% level
#' which(IWT_result$adjusted_pval < 0.05)
IWT1 <- # nolint: object_name_linter.
  function(
    data,
    mu = 0,
    B = 1000L, # nolint: object_name_linter.
    dx = NULL,
    recycle = TRUE
  ) {
    iwt1(data = data, mu = mu, n_perm = B, dx = dx, recycle = recycle)
  }

#' @param n_perm An integer value specifying the number of permutations for the
#'   permutation tests. Defaults to `1000L`.
#' @rdname IWT1
#' @export
iwt1 <- function(
  data,
  mu = 0,
  n_perm = 1000L,
  dx = NULL,
  recycle = TRUE
) {
  inputs <- onesample2coeffs(data, mu, dx = dx)
  coeff <- inputs$coeff
  mu_eval <- inputs$mu

  n <- dim(coeff)[1]
  p <- dim(coeff)[2]
  data_eval <- coeff <- coeff -
    matrix(data = mu_eval, nrow = n, ncol = p, byrow = TRUE)

  cli::cli_h1("Point-wise tests")

  t0 <- abs(colMeans(coeff))^2

  # Run permutations in parallel via mirai_map().
  # Each task produces one numeric vector of length p.
  perm_tasks <- mirai::mirai_map(
    seq_len(n_perm),
    \(.x) {
      signs <- stats::rbinom(n, 1, 0.5) * 2 - 1
      abs(colMeans(coeff * signs))^2
    },
    .args = list(coeff = coeff, n = n)
  )
  t_coeff <- do.call(rbind, perm_tasks[])

  pval <- colSums(
    t_coeff >= matrix(t0, nrow = n_perm, ncol = p, byrow = TRUE)
  ) /
    n_perm

  cli::cli_h1("Interval-wise tests")

  matrice_pval_asymm <- matrix(nrow = p, ncol = p)
  matrice_pval_asymm[p, ] <- pval[seq_len(p)]
  t0_2x <- c(t0, t0)
  t_coeff_2x <- cbind(t_coeff, t_coeff)

  maxrow <- 1L
  if (recycle) {
    for (i in (p - 1L):maxrow) {
      for (j in seq_len(p)) {
        inf <- j
        sup <- (p - i) + j
        t0_temp <- sum(t0_2x[inf:sup])
        t_temp <- rowSums(t_coeff_2x[, inf:sup, drop = FALSE])
        matrice_pval_asymm[i, j] <- sum(t_temp >= t0_temp) / n_perm
      }
      cli::cli_h1(
        "Creating the p-value matrix: end of row {p - i + 1} out of {p}"
      )
    }
  } else {
    for (i in (p - 1L):maxrow) {
      for (j in seq_len(i)) {
        inf <- j
        sup <- (p - i) + j
        t0_temp <- sum(t0_2x[inf:sup])
        t_temp <- rowSums(t_coeff_2x[, inf:sup, drop = FALSE])
        matrice_pval_asymm[i, j] <- sum(t_temp >= t0_temp) / n_perm
      }
      cli::cli_h1(
        "Creating the p-value matrix: end of row {p - i + 1} out of {p}"
      )
    }
  }

  corrected_pval_matrix <- pval_correct(matrice_pval_asymm)
  corrected_pval <- corrected_pval_matrix[1, ]

  cli::cli_h1("Interval-Wise Testing completed")

  out <- list(
    test = "1pop",
    mu = mu_eval,
    adjusted_pval = corrected_pval,
    unadjusted_pval = pval,
    pval_matrix = matrice_pval_asymm,
    data_eval = data_eval
  )
  class(out) <- "IWT1"
  out
}
