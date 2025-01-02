#' Two population local testing procedures
#'
#' The function implements local testing procedures for testing mean differences
#' between two functional populations. Functional data are tested locally and
#' unadjusted and adjusted p-value functions are provided. The unadjusted
#' p-value function controls the point-wise error rate. The adjusted p-value
#' function can be computed according to the following methods:
#' 
#' - global testing (controlling the FWER weakly)
#' - interval-wise testing (controlling the interval-wise error rate)
#' - threshold-wise testing (controlling the FWER asymptotically)
#' - partition closed testing (controlling the FWER on a partition)
#' - functional Benjamini Hochberg (controlling the FDR)
#'
#' @param data1 First population's data. Either pointwise evaluations of the
#'   functional data set on a uniform grid, or a \code{fd} object from the
#'   package \code{fda}. If pointwise evaluations are provided, \code{data2} is
#'   a matrix of dimensions \code{c(n1,J)}, with \code{J} evaluations on columns
#'   and \code{n1} units on rows.
#' @param data2 Second population's data. Either pointwise evaluations of the
#'   functional data set on a uniform grid, or a \code{fd} object from the
#'   package \code{fda}. If pointwise evaluations are provided, \code{data2} is
#'   a matrix of dimensions \code{c(n1,J)}, with \code{J} evaluations on columns
#'   and \code{n2} units on rows.
#' @param method A character string specifying the method chosen to adjust the
#'   p-value. Should be one of the following: "\code{Global}", "\code{IWT}",
#'   "\code{TWT}", "\code{PCT}", "\code{FDR}".
#' @inheritParams IWT2
#' @param partition Used only if \code{method}="\code{PCT}". The partition to be
#'   used for PCT procedure. Default is \code{NULL}.
#'
#' @return An object of class `fdatest2` containing at least the following
#'   components:
#' 
#' - `test`: a string vector indicating the type of test performed. In this case
#'  equal to `"2pop"`.
#' - `mu`: evaluation on a grid of the functional mean difference under the null
#' hypothesis (as entered by the user).
#' - `unadjusted_pval`: evaluation on a grid of the unadjusted p-value function.
#' - `adjusted_pval`: evaluation on a grid of the adjusted p-value function.
#' - `data.eval`: evaluation on a grid of the functional data.
#' - `ord_labels`: vector of labels indicating the group membership of
#' `data.eval`.
#'
#' @seealso See also [`plot.fdatest2()`] for plotting the results.
#'
#' @references
#' Abramowicz, K., Pini, A., Schelin, L., Stamm, A., & Vantini, S. (2022).
#' “Domain selection and familywise error rate for functional data: A unified
#' framework. \emph{Biometrics} 79(2), 1119-1132.
#'
#' Pini, A., & Vantini, S. (2017). Interval-wise testing for functional data.
#' \emph{Journal of Nonparametric Statistics}, 29(2), 407-424
#'
#' A. Pini and S. Vantini (2017). The Interval Testing Procedure: Inference for
#' Functional Data Controlling the Family Wise Error Rate on Intervals.
#' Biometrics 73(3): 835–845.
#'
#' Lundtorp Olsen, N., Pini, A., & Vantini, S. (2021). False discovery rate for
#' functional data \emph{TEST} 30, 784–809.
#'
#' @export
#' @examples
#' # Importing the NASA temperatures data set
#' data(NASAtemp)
#'
#' # Performing the TWT for two populations
#' TWT.result <- fdatest2(
#'   NASAtemp$paris, NASAtemp$milan, 
#'   method = "TWT", B = 10L
#' )
#'
#' # Plotting the results of the TWT
#' plot(
#'   TWT.result, 
#'   xrange = c(0, 12), 
#'   main = 'TWT results for testing mean differences'
#' )
#'
#' # Selecting the significant components at 5% level
#' which(TWT.result$adjusted_pval < 0.05)
#' 
#' # Performing the IWT for two populations
#' IWT.result <- fdatest2(
#'   NASAtemp$paris, NASAtemp$milan, 
#'   method = "IWT", B = 10L
#' )
#'
#' # Plotting the results of the IWT
#' plot(
#'   IWT.result, 
#'   xrange = c(0, 12), 
#'   main = 'IWT results for testing mean differences'
#' )
#'
#' # Selecting the significant components at 5% level
#' which(IWT.result$adjusted_pval < 0.05)
fdatest2 <- function(data1, data2, method, 
                     mu = 0, 
                     B = 1000L, 
                     paired = FALSE, 
                     dx = NULL, 
                     alternative = "two.sided", 
                     recycle = TRUE, 
                     partition = NULL, 
                     verbose = TRUE) {
  alternative <- rlang::arg_match(alternative, values = AVAILABLE_ALTERNATIVES())
  method <- rlang::arg_match(method, values = AVAILABLE_METHODS())
  
  if (method == "PCT" && is.null(partition)) {
    cli::cli_abort(
      'When the {.arg method} argument is set to {.code "PCT"}, the {.arg
      partition} argument should be explicitly provided.'
    )
  }
  
  out <- switch(
    method,
    IWT = IWT2(
      data1 = data1,
      data2 = data2,
      mu = mu,
      B = B,
      paired = paired,
      dx = dx,
      recycle = recycle,
      alternative = alternative,
      verbose = verbose
    ), 
    TWT = TWT2(
      data1 = data1,
      data2 = data2,
      mu = mu,
      B = B,
      paired = paired,
      dx = dx,
      alternative = alternative
    ), 
    PCT = PCT2(
      data1 = data1,
      data2 = data2,
      mu = mu,
      B = B,
      paired = paired,
      dx = dx,
      alternative = alternative,
      partition = partition
    ), 
    Global = Global2(
      data1 = data1,
      data2 = data2,
      mu = mu,
      B = B,
      paired = paired,
      dx = dx
    ), 
    FDR = FDR2(
      data1 = data1,
      data2 = data2,
      mu = mu,
      B = B,
      paired = paired,
      dx = dx,
      alternative = alternative
    )
  )
  out$method <- method
  out
}

