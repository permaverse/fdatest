#' @title Two populations Interval Testing Procedure with B-spline basis
#'
#' @description The function implements the Interval Testing Procedure for testing the difference between two functional populations evaluated on a uniform grid. Data are represented by means of the B-spline basis and the significance of each basis coefficient is tested with an interval-wise control of the Family Wise Error Rate. The default parameters of the basis expansion lead to the piece-wise interpolating function.
#'
#' @param data1 Pointwise evaluations of the first population's functional data set on a uniform grid. \code{data1} is a matrix of dimensions \code{c(n1,J)}, with \code{J} evaluations on columns and \code{n1} units on rows.
#'
#' @param data2 Pointwise evaluations of the second population's functional data set on a uniform grid. \code{data2} is a matrix of dimensions \code{c(n2,J)}, with \code{J} evaluations on columns and \code{n2} units on rows.
#'
#' @param mu The difference between the first functional population and the second functional population under the null hypothesis. Either a constant (in this case, a constant function is used) or a \code{J}-dimensional vector containing the evaluations on the same grid which \code{data} are evaluated. The default is \code{mu=0}.
#'
#' @param order Order of the B-spline basis expansion. The default is \code{order=2}.
#'
#' @param nknots Number of knots of the B-spline basis expansion. The default is \code{nknots=dim(data1)[2]}.
#'
#' @param B The number of iterations of the MC algorithm to evaluate the p-values of the permutation tests. The defualt is \code{B=1000}.
#'
#' @param paired A logical indicating whether the test is paired. The default is \code{FALSE}.
#'
#' @return \code{ITP2bspline} returns an object of \code{\link{class}} "\code{ITP2}".
#' An object of class "\code{ITP2}" is a list containing at least the following components:
#' \item{basis}{String vector indicating the basis used for the first phase of the algorithm. In this case equal to \code{"B-spline"}.}
#' \item{test}{String vector indicating the type of test performed. In this case equal to \code{"2pop"}.}
#' \item{mu}{Difference between the first functional population and the second functional population under the null hypothesis
#' (as entered by the user).}
#' \item{paired}{Logical indicating whether the test is paired (as entered by the user).}
#' \item{coeff}{Matrix of dimensions \code{c(n,p)} of the \code{p} coefficients of the B-spline basis expansion,
#' with \code{n=n1+n2}. Rows are associated to units and columns to the basis index.
#' The first \code{n1} rows report the coefficients of the first population units and the following \code{n2} rows report the coefficients
#' of the second population units}
#' \item{pval}{Unadjusted p-values for each basis coefficient.}
#' \item{pval.matrix}{Matrix of dimensions \code{c(p,p)} of the p-values of the multivariate tests. The element \code{(i,j)} of matrix \code{pval.matrix} contains the p-value of the joint NPC test of the components \code{(j,j+1,...,j+(p-i))}.}
#' \item{adjusted.pval}{Adjusted p-values for each basis coefficient.}
#' \item{labels}{Labels indicating the population membership of each data.}
#' \item{data.eval}{Evaluation on a fine uniform grid of the functional data obtained through the basis expansion.}
#' \item{heatmap.matrix}{Heatmap matrix of p-values (used only for plots).}
#'
#' @seealso See also \code{\link{plot.ITP2}} and \code{\link{ITPimage}} for plotting the results,
#'  \code{\link{ITP2fourier}} for ITP based on Fourier basis, \code{\link{IWT2}} for a two-sample test that is not based on
#'  an a-priori selected basis expansion.
#'
#' @examples
#' # Performing the ITP
#' ITP.result <- ITP2bspline(NASAtemp$milan,NASAtemp$paris,nknots=50,B=1000)
#'
#' # Plotting the results of the ITP
#' plot(ITP.result,main='NASA data',xrange=c(1,365),xlab='Day')
#'
#' # Plotting the p-values heatmap
#' ITPimage(ITP.result,abscissa.range=c(0,12))
#'
#' # Selecting the significant components at 5% level
#' which(ITP.result$adjusted.pval < 0.05)
#'
#' @references A. Pini and S. Vantini (2017).
#' The Interval Testing Procedure: Inference for Functional Data Controlling the Family Wise Error Rate on Intervals. Biometrics 73(3): 835–845.
#'
#' @export


ITP2bspline <-
function(data1,data2,mu=0,order=2,nknots=dim(data1)[2],B=10000,paired=FALSE){
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

  data1 <- as.matrix(data1)
  data2 <- as.matrix(data2)

  n1 <- dim(data1)[1]
  n2 <- dim(data2)[1]
  J <- dim(data1)[2]
  n <- n1+n2
  etichetta_ord <- c(rep(1,n1),rep(2,n2))
  data2 <- data2 - matrix(data=mu,nrow=n2,ncol=J)

  print('First step: basis expansion')
  #splines coefficients:
  eval <- rbind(data1,data2)
  bspl.basis <- fda::create.bspline.basis(c(1,J),norder=order,breaks=seq(1,J,length.out=nknots))
  ascissa <- seq(1,J,1)

  data.fd <- fda::Data2fd(t(eval),ascissa,bspl.basis)
  coeff <- t(data.fd$coef)
  p <- dim(coeff)[2]

  #functional data
  npt <- 1000
  ascissa.2 <- seq(1,J,length.out=npt)
  bspl.eval.smooth <- fda::eval.basis(ascissa.2,bspl.basis)
  data.eval <- t(bspl.eval.smooth %*% t(coeff))
  data.eval[(n1+1):n,] <- data.eval[(n1+1):n,] + matrix(data=mu,nrow=n2,ncol=npt)

  print('Second step: joint univariate tests')
  #univariate permutations
  T0 <- abs(colMeans(coeff[1:n1,,drop=FALSE]) - colMeans(coeff[(n1+1):n,,drop=FALSE])) #sample mean difference
  T_coeff <- matrix(ncol=p,nrow=B)
  for (perm in 1:B){
    if(paired==TRUE){
      if.perm <- stats::rbinom(n1,1,0.5)
      coeff_perm <- coeff
      for(couple in 1:n1){
        if(if.perm[couple]==1){
          coeff_perm[c(couple,n1+couple),] <- coeff[c(n1+couple,couple),]
        }
      }
    }else if(paired==FALSE){
      permutazioni <- sample(n)
      coeff_perm <- coeff[permutazioni,]
    }
    T_coeff[perm,] <- abs(colMeans(coeff_perm[1:n1,,drop=FALSE]) - colMeans(coeff_perm[(n1+1):n,,drop=FALSE]))
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
  ITPresult <- list(basis='B-spline',test='2pop',mu=mu,paired=as.character(paired),coeff=coeff,pval=pval,pval.matrix=matrice_pval_asymm,adjusted.pval=adjusted.pval,labels=etichetta_ord,data.eval=data.eval,heatmap.matrix=matrice_pval_symm)
  class(ITPresult) = 'ITP2'
  return(ITPresult)
}
