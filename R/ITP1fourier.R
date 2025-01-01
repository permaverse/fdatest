#' @rdname IWT1
#' 
#' @param maxfrequency The maximum frequency to be used in the Fourier basis
#'   expansion of data. Defaults to `floor(dim(data)[2] / 2)`, leading to an
#'   interpolating expansion.
#' 
#' @export
ITP1fourier <- function(data, 
                        mu = 0, 
                        B = 1000, 
                        maxfrequency = floor(dim(data)[2] / 2)) {
  lifecycle::deprecate_soft(
    when = "2.2.0", 
    what = "ITP1fourier()", 
    with = "IWT1()"
  )
  
  IWT1(data = data, mu = mu, B = B)
}
