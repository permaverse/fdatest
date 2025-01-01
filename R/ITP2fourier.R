#' @rdname IWT2
#' 
#' @param maxfrequency The maximum frequency to be used in the Fourier basis
#'   expansion of data. Defaults to `floor(dim(data1)[2] / 2)`, leading to an
#'   interpolating expansion.
#' 
#' @export
ITP2fourier <- function(data1, data2, 
                        mu = 0, 
                        B = 1000, 
                        paired = FALSE, 
                        maxfrequency = floor(dim(data1)[2] / 2)) {
  lifecycle::deprecate_soft(
    when = "2.2.0", 
    what = "ITP2fourier()", 
    with = "IWT2()"
  )
  
  IWT2(
    data1 = data1,
    data2 = data2,
    mu = mu,
    B = B,
    paired = paired,
    verbose = TRUE
  )
}
