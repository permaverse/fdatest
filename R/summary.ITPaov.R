#' S3 specialization of the `summary()` function for ITPaov objects
#' 
#' This function is a specialization of the `summary()` function for objects of
#' class `ITPaov`. It returns a list with the minimum p-values for the factors,
#' the range of functional R-squared values, and the minimum p-value for the
#' F-test.
#'
#' @param object An object of class `ITPaov`.
#' @inheritParams base::summary
#'
#' @return A list with the minimum p-values for the factors, the range of
#' functional R-squared values, and the minimum p-value for the F-test.
#' 
#' @export
summary.ITPaov <- function(object, ...) {
  printresult = vector('list')
  printresult$call = object$call
  printresult$factors = matrix(data=apply(object$adjusted.pval.factors,1,min),ncol=1)
  var.names = rownames(object$adjusted.pval.factors)
  rownames(printresult$factors) = var.names
  printresult$factors = as.data.frame(printresult$factors)
  signif = rep('',length(var.names))
  signif[which(printresult$factors[,1] <0.001)] = '***'
  signif[which(printresult$factors[,1] <0.01 & printresult$factors[,1] >= 0.001)] = '**'
  signif[which(printresult$factors[,1] <0.05 & printresult$factors[,1] >= 0.01)] = '*'
  signif[which(printresult$factors[,1] <0.1 & printresult$factors[,1] >= 0.05)] = '.'
  printresult$factors[,2] = signif
  colnames(printresult$factors) = c('Minimum p-value','')
  
  printresult$R2 = as.matrix(range(object$R2.eval))
  colnames(printresult$R2) = 'Range of functional R-squared'
  rownames(printresult$R2) = c('Min R-squared', 'Max R-squared')
  
  printresult$ftest = as.matrix(min(object$adjusted.pval.F))
  printresult$ftest = as.data.frame(printresult$ftest)
  signif.f = ''
  signif.f[which(printresult$ftest[,1] <0.001)] = '***'
  signif.f[which(printresult$ftest[,1] <0.01 & printresult$ftest[,1] >= 0.001)] = '**'
  signif.f[which(printresult$ftest[,1] <0.05 & printresult$ftest[,1] >= 0.01)] = '*'
  signif.f[which(printresult$ftest[,1] <0.1 & printresult$ftest[,1] >= 0.05)] = '.'
  printresult$ftest[,2] = signif.f
  colnames(printresult$ftest) = c('Minimum p-value','')
  printresult
  
}
