#' @title One population Interval Testing Procedure with B-spline basis
#'
#' @description The function implements the Interval Testing Procedure for testing the center of symmetry of a
#' functional population evaluated on a uniform grid. Data are represented by means of the B-spline expansion and the
#' significance of each basis coefficient is tested with an interval-wise control of the Family Wise Error Rate.
#' The default parameters of the basis expansion lead to the piece-wise interpolating function.
#'
#' @param data Pointwise evaluations of the functional data set on a uniform grid.
#' It is a matrix of dimensions \code{c(n,J)}, with \code{J} evaluations on columns and \code{n} units on rows.
#'
#' @param mu The center of symmetry under the null hypothesis: either a constant
#' (in this case, a constant function is used) or a \code{J}-dimensional vector containing the evaluations on the
#' same grid which \code{data} are evaluated. The default is \code{mu=0}.
#'
#' @param order Order of the B-spline basis expansion. The default is \code{order=2}.
#'
#' @param nknots Number of knots of the B-spline basis expansion. The default is \code{nknots=dim(data)[2]}.
#'
#' @param B The number of iterations of the MC algorithm to evaluate the p-values of the permutation tests. The defualt is \code{B=1000}.
#'
#' @return \code{ITP1bspline} returns an object of \code{\link{class}} "\code{ITP1}".
#' An object of class "\code{ITP1}" is a list containing at least the following components:
#' \item{basis}{String vector indicating the basis used for the first phase of the algorithm. In this case equal to \code{"B-spline"}.}
#' \item{test}{String vector indicating the type of test performed. In this case equal to \code{"1pop"}.}
#' \item{mu}{Center of symmetry under the null hypothesis (as entered by the user).}
#' \item{coeff}{Matrix of dimensions \code{c(n,p)} of the \code{p} coefficients of the B-spline basis expansion. Rows are associated to units and columns to the basis index.}
#' \item{pval}{Unadjusted p-values for each basis coefficient.}
#' \item{pval.matrix}{Matrix of dimensions \code{c(p,p)} of the p-values of the multivariate tests. The element \code{(i,j)} of matrix \code{pval.matrix} contains the p-value of the joint NPC test of the components \code{(j,j+1,...,j+(p-i))}.}
#' \item{adjusted.pval}{Adjusted p-values for each basis coefficient.}
#' \item{labels}{Labels indicating the population membership of each data (in this case always equal to \code{1}).}
#' \item{data.eval}{Evaluation on a fine uniform grid of the functional data obtained through the basis expansion.}
#' \item{heatmap.matrix}{Heatmap matrix of p-values (used only for plots).}
#'
#' @seealso See also \code{\link{plot.ITP1}} and \code{\link{ITPimage}} for plotting the results,
#'  \code{\link{ITP1fourier}} for ITP based on Fourier basis, \code{\link{IWT1}} for a one-sample test that is not based on
#'  an a-priori selected basis expansion.
#'
#' @examples
#' # Performing the ITP for one population with the B-spline basis
#' ITP.result <- ITP1bspline(NASAtemp$paris,mu=4,nknots=50,B=1000)
#'
#' # Plotting the results of the ITP
#' plot(ITP.result,xrange=c(0,12),main='Paris temperatures')
#'
#' # Plotting the p-value heatmap
#' ITPimage(ITP.result,abscissa.range=c(0,12))
#'
#' # Selecting the significant components at 5% level
#' which(ITP.result$adjusted.pval < 0.05)
#'
#' @references A. Pini and S. Vantini (2017).
#' The Interval Testing Procedure: Inference for Functional Data Controlling the Family Wise Error Rate on Intervals. Biometrics 73(3): 835–845.
#'
#' @export

ITP1bspline <-
function(data,mu=0,order=2,nknots=dim(data)[2],B=1000){
  fisher_cf_L <- function(L){ #fisher on rows of the matrix L
    return(-2*rowSums(log(L)))
  }
  fisher_cf <- function(lambda){ #fisher on vector lambda
    return(-2*sum(log(lambda)))
  }
  pval.correct <- function(pval.matrix){
    matrice_pval_2_2x <- cbind(pval.matrix,pval.matrix)
    p <- dim(pval.matrix)[2]
    matrice_pval_2_2x <- matrice_pval_2_2x[,(2*p):1]
    adjusted.pval <- numeric(p)
    for(var in 1:p){
      pval_var <- matrice_pval_2_2x[p,var]
      inizio <- var
      fine <- var #inizio fisso, fine aumenta salendo nelle righe
      for(riga in (p-1):1){
        fine <- fine + 1
        pval_cono <- matrice_pval_2_2x[riga,inizio:fine]
        pval_var <- max(pval_var,pval_cono)
      }
      adjusted.pval[var] <- pval_var
    }
    adjusted.pval <- adjusted.pval[p:1]
    return(adjusted.pval)
  }
  data <- as.matrix(data)

  n <- dim(data)[1]
  J <- dim(data)[2]
  labels <- rep(1,n)
  data <- data - matrix(data=mu,nrow=n,ncol=J,byrow=TRUE)

  print('First step: basis expansion')
  #splines coefficients:
  bspl.basis <- fda::create.bspline.basis(c(1,J),norder=order,breaks=seq(1,J,length.out=nknots))
  ascissa <- seq(1,J,1)

  data.fd <- fda::Data2fd(t(data),ascissa,bspl.basis)
  coeff <- t(data.fd$coef)

  p <- dim(coeff)[2]

  #functional data
  npt <- 1000
  ascissa.2 <- seq(1,J,length.out=npt)
  bspl.eval.smooth <- fda::eval.basis(ascissa.2,bspl.basis)
  data.eval <- t(bspl.eval.smooth %*% t(coeff) )
  data.eval <- data.eval + matrix(data=mu,nrow=n,ncol=npt)

  #univariate permutations
  print('Second step: joint univariate tests')
  T0 <- abs(colMeans(coeff))  #sample mean
  T_coeff <- matrix(ncol=p,nrow=B)
  for (perm in 1:B){
    signs <- stats::rbinom(n,1,0.5)*2 - 1
    coeff_perm <- coeff*signs
    T_coeff[perm,] <- abs(colMeans(coeff_perm))
  }
  pval <- numeric(p)
  for(i in 1:p){
    pval[i] <- sum(T_coeff[,i]>=T0[i])/B
  }

  #combination
  print('Third step: interval-wise combination and correction')
  q <- numeric(B)
  L <- matrix(nrow=B,ncol=p)
  for(j in 1:p){
    ordine <- sort.int(T_coeff[,j],index.return=T)$ix
    q[ordine] <- (B:1)/(B)
    L[,j] <- q
  }

  #asymmetric combination matrix:
  matrice_pval_asymm <- matrix(nrow=p,ncol=p)
  matrice_pval_asymm[p,] <- pval[1:p]
  pval_2x <- c(pval,pval)
  L_2x <- cbind(L,L)
  for(i in (p-1):1){
    for(j in 1:p){
      inf <- j
      sup <- (p-i)+j
      T0_temp <- fisher_cf(pval_2x[inf:sup])
      T_temp <- fisher_cf_L(L_2x[,inf:sup])
      pval_temp <- sum(T_temp>=T0_temp)/B
      matrice_pval_asymm[i,j] <- pval_temp
    }
    print(paste('creating the p-value matrix: end of row ',as.character(p-i+1),' out of ',as.character(p),sep=''))
  }

  #symmetric combination matrix
  matrice_pval_symm <- matrix(nrow=p,ncol=4*p)
  for(i in 0:(p-1)){
    for(j in 1:(2*p)){
      matrice_pval_symm[p-i,j+i+p] <- matrice_pval_asymm[p-i,(j+1)%/%2]
      if(j+i>2*p-i){
        matrice_pval_symm[p-i,j+i-p] <- matrice_pval_asymm[p-i,(j+1)%/%2]
      }
    }
  }

  adjusted.pval <- pval.correct(matrice_pval_asymm)
  print('Interval Testing Procedure completed')
  ITP.result <- list(basis='B-spline',test='1pop',mu=mu,coeff=coeff,pval=pval,pval.matrix=matrice_pval_asymm,adjusted.pval=adjusted.pval,labels=labels,data.eval=data.eval,heatmap.matrix=matrice_pval_symm)
  class(ITP.result) = 'ITP1'
  return(ITP.result)
}
