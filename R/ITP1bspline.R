#' @rdname IWT1
#'
#' @param order Order of the B-spline basis expansion. Defaults to `2L`.
#' @param nknots Number of knots of the B-spline basis expansion. Defaults to
#'   `dim(data)[2]`.
#' 
#' @export
ITP1bspline <- function(data, mu = 0, B = 1000, order = 2, nknots = dim(data)[2]) {
  lifecycle::deprecate_soft(
    when = "2.2.0", 
    what = "ITP1bspline()", 
    with = "IWT1()"
  )
  
  IWT1(data = data, mu = mu, B = B)
}
