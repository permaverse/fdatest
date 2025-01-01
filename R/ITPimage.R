#' Heatmap plot of the Interval Testing Procedure results
#'
#' Plotting function creating a graphical output of the ITP: the p-value
#' heat-map, the plot of the corrected p-values, and the plot of the functional
#' data.
#'
#' @param ITP.result Results of the ITP, as created by
#'   \code{\link{ITP1bspline}}, \code{\link{ITP1fourier}},
#'   \code{\link{ITP2bspline}}, \code{\link{ITP2fourier}}, and
#'   \code{\link{ITP2pafourier}}.
#' @param alpha Threshold for the interval-wise error rate used for the
#'   hypothesis test. The default is \code{alpha}=0.05.
#' @param abscissa.range Range of the plot abscissa. The default is
#'   \code{c(0,1)}.
#' @param nlevel Number of desired color levels for the p-value heatmap. The
#'   default is \code{nlevel=20}.
#'
#' @return NULL
#'
#' @seealso See \code{\link{plot.ITP1}}, \code{\link{plot.ITP2}},
#'   \code{\link{plot.ITPlm}}, and \code{\link{plot.ITPaov}} for the plot method
#'   applied to the ITP results of one- and two-population tests, linear models,
#'   and ANOVA, respectively. See also \code{\link{ITP1bspline}},
#'   \code{\link{ITP1fourier}}, \code{\link{ITP2bspline}},
#'   \code{\link{ITP2fourier}}, and \code{\link{ITP2pafourier}} for applying the
#'   ITP.
#'
#' @references
#' A. Pini and S. Vantini (2017). The Interval Testing Procedure: Inference for
#' Functional Data Controlling the Family Wise Error Rate on Intervals.
#' Biometrics 73(3): 835–845.
#'
#' @export
#' @examples
#' # Performing the ITP for two populations with the B-spline basis
#' ITP.result <- ITP2bspline(
#'   NASAtemp$milan, NASAtemp$paris, 
#'   nknots = 20, 
#'   B = 10L
#' )
#'
#' # Plotting the results of the ITP
#' ITPimage(ITP.result,abscissa.range=c(0,12))
#'
#' # Selecting the significant components for the radius at 5% level
#' which(ITP.result$corrected.pval < 0.05)
ITPimage <- function(ITP.result,
                     alpha = 0.05,
                     abscissa.range = c(0, 1),
                     nlevel = 20) {
  lifecycle::deprecate_soft(
    when = "2.2.0", 
    what = "ITPimage()", 
    with = "IWTimage()"
  )
  
  IWTimage(
    IWT_result = ITP.result, 
    alpha = alpha, 
    abscissa_range = abscissa.range, 
    nlevel = nlevel
  )
}
