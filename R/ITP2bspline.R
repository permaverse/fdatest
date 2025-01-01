#' @rdname IWT2
#' 
#' @param order Order of the B-spline basis expansion. Defaults to `2L`.
#' @param nknots Number of knots of the B-spline basis expansion. Defaults to
#'   `dim(data1)[2]`.
#'
#' @export
ITP2bspline <- function(data1, data2, 
                        mu = 0, 
                        B = 1000, 
                        paired = FALSE, 
                        order = 2, 
                        nknots = dim(data1)[2]) {
  lifecycle::deprecate_soft(
    when = "2.2.0", 
    what = "ITP2bspline()", 
    with = "IWT2()"
  )
  
  IWT2(
    data1 = data1,
    data2 = data2,
    mu = mu,
    B = B,
    paired = paired,
    verbose = FALSE
  )
}
