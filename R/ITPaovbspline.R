#' @title Interval Testing Procedure for testing unctional analysis of variance with B-spline basis
#'
#' @description The function implements the Interval Testing Procedure for testing for significant differences between several functional population evaluated on a uniform grid, in a functional analysis of variance setting. Data are represented by means of the B-spline basis and the significance of each basis coefficient is tested with an interval-wise control of the Family Wise Error Rate. The default parameters of the basis expansion lead to the piece-wise interpolating function.
#'
#' @param formula An object of class "\code{\link{formula}}" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#'
#' @param order Order of the B-spline basis expansion. The default is \code{order=2}.
#'
#' @param nknots Number of knots of the B-spline basis expansion. The default is \code{nknots=dim(data1)[2]}.
#'
#' @param B The number of iterations of the MC algorithm to evaluate the p-values of the permutation tests. The defualt is \code{B=1000}.
#'
#' @param method Permutation method used to calculate the p-value of permutation tests. Choose "\code{residuals}" for the permutations of residuals under the reduced model, according to the Freedman and Lane scheme, and "\code{responses}" for the permutation of the responses, according to the Manly scheme.
#'
#' @return \code{ITPaovbspline} returns an object of \code{\link{class}} "\code{ITPaov}". The function \code{summary} is used to obtain and print a summary of the results.
#' An object of class "\code{ITPaov}" is a list containing at least the following components:
#' \item{call}{The matched call.}
#' \item{design.matrix}{The design matrix of the functional-on-scalar linear model.}
#' \item{basis}{String vector indicating the basis used for the first phase of the algorithm. In this case equal to \code{"B-spline"}.}
#' \item{coeff}{Matrix of dimensions \code{c(n,p)} of the \code{p} coefficients of the B-spline basis expansion. Rows are associated to units and columns to the basis index.}
#' \item{coeff.regr}{Matrix of dimensions \code{c(L+1,p)} of the \code{p} coefficients of the B-spline basis expansion of the intercept (first row) and the \code{L} effects of the covariates specified in \code{formula}. Columns are associated to the basis index.}
#' \item{pval.F}{Unadjusted p-values of the functional F-test for each basis coefficient.}
#' \item{pval.matrix.F}{Matrix of dimensions \code{c(p,p)} of the p-values of the multivariate F-tests. The element \code{(i,j)} of matrix \code{pval.matrix} contains the p-value of the joint NPC test of the components \code{(j,j+1,...,j+(p-i))}.}
#' \item{adjusted.pval.F}{Adjusted p-values of the functional F-test for each basis coefficient.}
#' \item{pval.factors}{Unadjusted p-values of the functional F-tests on each factor of the analysis of variance, separately (rows) and each basis coefficient (columns).}
#' \item{pval.matrix.factors}{Array of dimensions \code{c(L+1,p,p)} of the p-values of the multivariate F-tests on factors. The element \code{(l,i,j)} of array \code{pval.matrix} contains the p-value of the joint NPC test on factor \code{l} of the components \code{(j,j+1,...,j+(p-i))}.}
#' \item{adjusted.pval.factors}{adjusted p-values of the functional F-tests on each factor of the analysis of variance (rows) and each basis coefficient (columns).}
#' \item{data.eval}{Evaluation on a fine uniform grid of the functional data obtained through the basis expansion.}
#' \item{coeff.regr.eval}{Evaluation on a fine uniform grid of the functional regression coefficients.}
#' \item{fitted.eval}{Evaluation on a fine uniform grid of the fitted values of the functional regression.}
#' \item{residuals.eval}{Evaluation on a fine uniform grid of the residuals of the functional regression.}
#' \item{R2.eval}{Evaluation on a fine uniform grid of the functional R-squared of the regression.}
#' \item{heatmap.matrix.F}{Heatmap matrix of p-values of functional F-test (used only for plots).}
#' \item{heatmap.matrix.factors}{Heatmap matrix of p-values of functional F-tests on each factor of the analysis of variance (used only for plots).}
#'
#' @seealso See \code{\link{summary.ITPaov}} for summaries and \code{\link{plot.ITPaov}} for plotting the results.
#' See \code{\link{IWTaov}} for a functional analysis of variance test that is not based on an a-priori selected basis expansion.
#' See also \code{\link{ITPlmbspline}} to fit and test a functional-on-scalar linear model applying the ITP, and \code{\link{ITP1bspline}}, \code{\link{ITP2bspline}}, \code{\link{ITP2fourier}}, \code{\link{ITP2pafourier}} for one-population and two-population tests.
#'
#'
#' @examples
#' temperature <- rbind(NASAtemp$milan,NASAtemp$paris)
#' groups <- c(rep(0,22),rep(1,22))
#'
#' # Performing the ITP
#' ITP.result <- ITPaovbspline(temperature ~ groups,B=100,nknots=20,order=3)
#'
#' # Summary of the ITP results
#' summary(ITP.result)
#'
#' # Plot of the ITP results
#' graphics::layout(1)
#' plot(ITP.result)
#'
#' # All graphics on the same device
#' graphics::layout(matrix(1:4,nrow=2,byrow=FALSE))
#' plot(ITP.result,main='NASA data', plot.adjpval = TRUE,xlab='Day',xrange=c(1,365))
#'
#' @references
#' A. Pini and S. Vantini (2017).
#' The Interval Testing Procedure: Inference for Functional Data Controlling the Family Wise Error Rate on Intervals. Biometrics 73(3): 835–845.
#' 
#' A. Pini and S. Vantini (2017).
#' The Interval Testing Procedure: Inference for Functional Data Controlling the Family Wise Error Rate on Intervals. Biometrics 73(3): 835–845.
#'
#' Pini, A., Vantini, S., Colosimo, B. M., & Grasso, M. (2018). Domain‐selective functional analysis of variance for supervised statistical profile monitoring of signal data. \emph{Journal of the Royal Statistical Society: Series C (Applied Statistics)} 67(1), 55-81.
#'
#' Abramowicz, K., Hager, C. K., Pini, A., Schelin, L., Sjostedt de Luna, S., & Vantini, S. (2018).
#' Nonparametric inference for functional‐on‐scalar linear models applied to knee kinematic hop data after injury of the anterior cruciate ligament. \emph{Scandinavian Journal of Statistics} 45(4), 1036-1061.
#'
#' D. Freedman and D. Lane (1983). A Nonstochastic Interpretation of Reported Significance Levels. \emph{Journal of Business & Economic Statistics} 1.4, 292-298.
#'
#' B. F. J. Manly (2006). Randomization, \emph{Bootstrap and Monte Carlo Methods in Biology}. Vol. 70. CRC Press.
#'
#' @export


ITPaovbspline <-
function(formula,order=2,nknots=dim(stats::model.response(stats::model.frame(formula)))[2],B=1000,method='residuals'){
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
      fine <- var
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
  extract.pval <- function(anova){
    result <- summary(anova)[[1]][,5]
    result <- result[-length(result)]
    return(result)
  }


  variables = all.vars(formula)
  y.name = variables[1]
  #data.all = stats::model.frame(formula)
  cl <- match.call()
  design.matrix = stats::model.matrix(formula)
  mf = stats::model.frame(formula)
  data = stats::model.response(mf)

  dummynames.all <- colnames(design.matrix)[-1]

  #var.names = variables[-1]
  #nvar = length(var.names)

  n <- dim(data)[1]
  J <- dim(data)[2]


  print('First step: basis expansion')
  #splines coefficients:
  eval <- data
  bspl.basis <- fda::create.bspline.basis(c(1,J),norder=order,breaks=seq(1,J,length.out=nknots))
  ascissa <- seq(1,J,1)

  data.fd <- fda::Data2fd(t(data),ascissa,bspl.basis)
  coeff <- t(data.fd$coef)
  p <- dim(coeff)[2]

  #functional data
  npt <- 1000
  ascissa.2 <- seq(1,J,length.out=npt)
  bspl.eval.smooth <- fda::eval.basis(ascissa.2,bspl.basis)
  data.eval <- t(bspl.eval.smooth %*% t(coeff))

  print('Second step: joint univariate tests')
  #univariate permutations
  formula.const <- deparse(formula[[3]],width.cutoff = 500L) #extracting the part after ~ on formula. this will not work if the formula is longer than 500 char
  coeffnames <- paste('coeff[,',as.character(1:p),']',sep='')
  formula.coeff <- paste(coeffnames,'~',formula.const)
  formula.coeff <- sapply(formula.coeff,stats::as.formula)

  formula.temp = coeff ~ design.matrix
  mf.temp = cbind(stats::model.frame(formula.temp),as.data.frame(design.matrix[,-1]))
  aovcoeff1 <- stats::aov(formula.coeff[[1]],data=mf.temp)
  var.names <- rownames(summary(aovcoeff1)[[1]])
  df.vars <- summary(aovcoeff1)[[1]][,1]
  df.residuals <- df.vars[length(df.vars)]
  var.names <- var.names[-length(var.names)]
  nvar = length(var.names)
  for(ii in 1:nvar){
    var.names[ii] <- gsub(' ' , '',var.names[ii])
  }

  index.vars <- cbind(c(2,(cumsum(df.vars)+2)[-length(df.vars)]),cumsum(df.vars)+1)
  regr0 = stats::lm.fit(design.matrix,coeff)
  #pval_parametric <- sapply(aov0,extract.pval)
  MS0 <- matrix(nrow=nvar+1,ncol=p)
  for(var in 1:(nvar+1)){
    MS0[var,] <- colSums(rbind(regr0$effects[index.vars[var,1]:index.vars[var,2],]^2))/df.vars[var]
  }
  # test statistic:
  T0_part <- MS0[1:nvar,] / matrix(MS0[nvar+1,],nrow=nvar,ncol=p,byrow=TRUE)
  Sigma <- chol2inv(regr0$qr$qr)
  resvar <- colSums(regr0$residuals^2)/regr0$df.residual


  if(nvar >1){
    T0_glob <- colSums((regr0$fitted - matrix(colMeans(regr0$fitted),nrow=n,ncol=p,byrow=TRUE))^2)/ ((nvar)*resvar)
  }else if(nvar==1){ #only one factor -> the permutation of the residuals is equivalent to the one of responses
    method <- 'responses'
    T0_glob <- colSums((regr0$fitted - matrix(colMeans(regr0$fitted),nrow=n,ncol=p,byrow=TRUE))^2)/ ((nvar)*resvar)
  }else if(nvar==0){
    method = 'responses' # model with only intercept -> the permutation of the residuals is equivalent to the one of responses
    T0_glob = numeric(p)
  }

  #calculate residuals
  if(method=='residuals'){
    #n residuals for each coefficient of basis expansion (1:p)
    #and for each partial test + global test (nvar+1)
    #saved in array of dim (nvar+1,n,p)
    design.matrix.names2 = design.matrix
    var.names2 = var.names
    if(length(grep('factor',formula.const))>0){
      index.factor = grep('factor',var.names)
      replace.names = paste('group',(1:length(index.factor)),sep='')
      var.names2[index.factor] = replace.names
      colnames(design.matrix.names2) = var.names2
    }

    residui = array(dim=c(nvar,n,p))
    fitted_part = array(dim=c(nvar,n,p)) # fitted values of the reduced model (different for each test)
    formula.coeff_part = vector('list',nvar)
    regr0_part = vector('list',nvar)
    dummy.interaz <- grep(':',dummynames.all)
    #coeff.perm_part = array(dim=c(nvar+1,n,p))
    for(ii in 1:(nvar)){ #no test on intercept
      var.ii = var.names2[ii]
      variables.reduced = var.names2[-which(var.names2==var.ii)] #removing the current variable to test

      if(length(grep(':',var.ii))>0){ # testing interaction
        #print('interaz')
        var12 <- strsplit(var.ii,':')
        var1 <- var12[[1]][1]
        var2 <- var12[[1]][2]
        dummy.test1 <- grep(var1,dummynames.all)
        dummy.test2 <- grep(var2,dummynames.all)
        dummy.test <- intersect(dummy.test1,dummy.test2)
        dummynames.reduced <- dummynames.all[-dummy.test]
      }else{
        #print('nointeraz')
        dummy.test <- grep(var.ii,dummynames.all)
        dummy.test <- setdiff(dummy.test,dummy.interaz)
        dummynames.reduced <- dummynames.all[-dummy.test]
      }


      if(nvar>1){
        formula.temp = paste(dummynames.reduced,collapse=' + ')
      }else{
        formula.temp = '1' #removing the only variable -> reduced model only has intercept term
      }

      formula.temp2 = coeff ~ design.matrix.names2
      mf.temp2 = cbind(stats::model.frame(formula.temp2)[-((p+1):(p+nvar+1))],as.data.frame(design.matrix.names2[,-1]))

      formula.coeff.temp <- paste(coeffnames,'~',formula.temp)
      formula.coeff_part[[ii]] <- sapply(formula.coeff.temp,stats::as.formula)
      regr0_part[[ii]] = lapply(formula.coeff_part[[ii]],stats::lm,data=mf.temp2)

      residui[ii,,] = simplify2array(lapply(regr0_part[[ii]],extract_residuals))
      fitted_part[ii,,] = simplify2array(lapply(regr0_part[[ii]],extract_fitted))
    }

  }


  T_glob <- matrix(ncol=p,nrow=B)
  T_part = array(dim=c(B,nvar,p))

  for (perm in 1:B){
    # the F test is the same for both methods
    if(nvar >0){
      permutazioni <- sample(n)
      coeff_perm <- coeff[permutazioni,]
    }else{ # testing intercept -> permute signs
      signs <- stats::rbinom(n,1,0.5)*2 - 1
      coeff_perm <- coeff*signs
    }

    regr_perm = stats::lm.fit(design.matrix,coeff_perm)
    Sigma <- chol2inv(regr_perm$qr$qr)
    resvar <- colSums(regr_perm$residuals^2)/regr0$df.residual

    if(nvar > 0)
      T_glob[perm,] <- colSums((regr_perm$fitted - matrix(colMeans(regr_perm$fitted),nrow=n,ncol=p,byrow=TRUE))^2)/ ((nvar)*resvar)

    # partial tests: differ depending on the method
    if(method=='responses'){
      MSperm <- matrix(nrow=nvar+1,ncol=p)
      for(var in 1:(nvar+1)){
        MSperm[var,] <- colSums(rbind(regr_perm$effects[index.vars[var,1]:index.vars[var,2],]^2))/df.vars[var]
      }
      # test statistic:
      T_part[perm,,] <- MSperm[1:nvar,] / matrix(MSperm[nvar+1,],nrow=nvar,ncol=p,byrow=TRUE)

    }else if(method=='residuals'){
      residui_perm = residui[,permutazioni,]
      aov_perm_part = vector('list',nvar)
      for(ii in 1:(nvar)){
        coeff_perm = fitted_part[ii,,] + residui_perm[ii,,]
        regr_perm = stats::lm.fit(design.matrix,coeff_perm)
        MSperm <- matrix(nrow=nvar+1,ncol=p)
        for(var in 1:(nvar+1)){
          MSperm[var,] <- colSums(rbind(regr_perm$effects[index.vars[var,1]:index.vars[var,2],]^2))/df.vars[var]
        }
        # test statistic:
        T_part[perm,ii,] <- (MSperm[1:nvar,] / matrix(MSperm[nvar+1,],nrow=nvar,ncol=p,byrow=TRUE))[ii,]
      }
    }
  }

  pval_glob <- numeric(p)
  pval_part = matrix(nrow=nvar,ncol=p)
  for(i in 1:p){
    pval_glob[i] <- sum(T_glob[,i]>=T0_glob[i])/B
    pval_part[,i] = colSums(T_part[,,i]>=matrix(T0_part[,i],nrow=B,ncol=nvar,byrow=TRUE))/B
  }

  #combination
  print('Third step: interval-wise combination and correction')
  q <- numeric(B)
  L_glob <- matrix(nrow=B,ncol=p)
  for(j in 1:p){
    ordine <- sort.int(T_glob[,j],index.return=T)$ix
    q[ordine] <- (B:1)/(B)
    L_glob[,j] <- q
  }

  L_part <- array(dim=c(B,nvar,p))
  for(j in 1:p){
    for(i in 1:(nvar)){
      ordine <- sort.int(T_part[,i,j],index.return=T)$ix
      q[ordine] <- (B:1)/(B)
      L_part[,i,j] <- q
    }
  }

  #asymmetric combination matrix:
  matrice_pval_asymm_glob <- matrix(nrow=p,ncol=p)
  matrice_pval_asymm_glob[p,] <- pval_glob[1:p]
  pval_2x_glob <- c(pval_glob,pval_glob)
  L_2x_glob <- cbind(L_glob,L_glob)

  matrice_pval_asymm_part <- array(dim=c(nvar,p,p))
  pval_2x_part <- cbind(pval_part,pval_part)
  L_2x_part = array(dim=c(B,nvar,p*2))
  for(ii in 1:(nvar)){
    matrice_pval_asymm_part[ii,p,] <- pval_part[ii,1:p]
    L_2x_part[,ii,] <- cbind(L_part[,ii,],L_part[,ii,])
  }

  for(i in (p-1):1){
    for(j in 1:p){
      inf <- j
      sup <- (p-i)+j
      T0_temp <- fisher_cf(pval_2x_glob[inf:sup])
      T_temp <- fisher_cf_L(L_2x_glob[,inf:sup])
      pval_temp <- sum(T_temp>=T0_temp)/B
      matrice_pval_asymm_glob[i,j] <- pval_temp
      for(ii in 1:(nvar )){
        T0_temp <- fisher_cf(pval_2x_part[ii,inf:sup])
        T_temp <- fisher_cf_L(L_2x_part[,ii,inf:sup])
        pval_temp <- sum(T_temp>=T0_temp)/B
        matrice_pval_asymm_part[ii,i,j] <- pval_temp
      }

    }
    print(paste('creating the p-value matrix: end of row ',as.character(p-i+1),' out of ',as.character(p),sep=''))
  }

  #symmetric combination matrix
  matrice_pval_symm_glob <- matrix(nrow=p,ncol=4*p)
  matrice_pval_symm_part <- array(dim=c(nvar,p,4*p))

  for(i in 0:(p-1)){
    for(j in 1:(2*p)){
      matrice_pval_symm_glob[p-i,j+i+p] <- matrice_pval_asymm_glob[p-i,(j+1)%/%2]
      if(j+i>2*p-i){
        matrice_pval_symm_glob[p-i,j+i-p] <- matrice_pval_asymm_glob[p-i,(j+1)%/%2]
      }
      for(ii in 1:(nvar)){
        matrice_pval_symm_part[ii,p-i,j+i+p] <- matrice_pval_asymm_part[ii,p-i,(j+1)%/%2]
        if(j+i>2*p-i){
          matrice_pval_symm_part[ii,p-i,j+i-p] <- matrice_pval_asymm_part[ii,p-i,(j+1)%/%2]
        }
      }
    }
  }

  adjusted.pval_glob <- pval.correct(matrice_pval_asymm_glob)
  adjusted.pval_part = matrix(nrow=nvar,ncol=p)
  for(ii in 1:(nvar)){
    adjusted.pval_part[ii,] = pval.correct(matrice_pval_asymm_part[ii,,])
  }

  coeff.regr = regr0$coeff
  coeff.t <- t(bspl.eval.smooth %*% t(coeff.regr))

  fitted.regr = regr0$fitted.values
  fitted.t <- t(bspl.eval.smooth %*% t(fitted.regr))

  rownames(adjusted.pval_part) = var.names
  rownames(coeff.t) = colnames(design.matrix)
  rownames(coeff.regr) = colnames(design.matrix)
  rownames(pval_part) = var.names

  residuals.t = data.eval - fitted.t
  ybar.t = colMeans(data.eval)
  R2.t = colSums((fitted.t - matrix(data=ybar.t,nrow=n,ncol=npt,byrow=TRUE))^2)/colSums((data.eval - matrix(data=ybar.t,nrow=n,ncol=npt,byrow=TRUE))^2)

  print('Interval Testing Procedure completed')

  ITPresult <- list(call=cl,design.matrix=design.matrix,basis='B-spline',coeff=coeff,coeff.regr=coeff.regr,
                    pval.F=pval_glob,pval.matrix.F=matrice_pval_asymm_glob,adjusted.pval.F=adjusted.pval_glob,
                    pval.factors=pval_part,pval.matrix.factors=matrice_pval_asymm_part,adjusted.pval.factors=adjusted.pval_part,
                    data.eval=data.eval,coeff.regr.eval=coeff.t,fitted.eval=fitted.t,residuals.eval=residuals.t,
                    R2.eval=R2.t,
                    heatmap.matrix.F=matrice_pval_symm_glob,
                    heatmap.matrix.factors=matrice_pval_symm_part
                    #pval_parametric=pval_parametric
                    )
  class(ITPresult) = 'ITPaov'
  return(ITPresult)
}
