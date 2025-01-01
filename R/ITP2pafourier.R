#' @rdname IWT2
#' 
#' @inheritParams ITP2fourier

#' @export
ITP2pafourier <- function(data1, data2, 
                          mu = 0, 
                          B = 1000, 
                          paired = FALSE, 
                          maxfrequency = floor(dim(data1)[2] / 2)) {
  lifecycle::deprecate_soft(
    when = "2.2.0", 
    what = "ITP2pafourier()", 
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
