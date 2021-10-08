# Partition closed testing function
# arguments:
# data1 and data2: functional data sets of the two groups. They can be either data of class fd or matrices with pointwise evaluations. 
# mu: mean difference under the null hypothesis
# B: number of permutations
# paired: whether the test is paired or not
# dx: optional. If data1 and data2 are fd objects, dx is used as grid step to evaluate them
# partition: vector defining the partition of the domain
PCT2 <- function(data1,data2,mu=0,B=1000,paired=FALSE,dx=NULL,partition){
  if(fda::is.fd(data1)){ # data1 is a functional data object
    rangeval1 <- data1$basis$rangeval
    rangeval2 <- data2$basis$rangeval
    if(is.null(dx)){
      dx <- (rangeval1[2]-rangeval1[1])*0.01
    }
    if(sum(rangeval1 == rangeval2)!=2){
      stop("rangeval of data1 and data2 must coincide.")
    }
    abscissa <- seq(rangeval1[1],rangeval1[2],by=dx)
    coeff1 <- t(fda::eval.fd(fdobj=data1,evalarg=abscissa))
    coeff2 <- t(fda::eval.fd(fdobj=data2,evalarg=abscissa))
    
  }else if(is.matrix(data1)){
    coeff1 <- data1
    coeff2 <- data2
  }else{
    stop("First argument must be either a functional data object or a matrix.")
  }
  
  if (fda::is.fd(mu)){ # mu is a functional data
    rangeval.mu <- mu$basis$rangeval
    if(sum(rangeval.mu == rangeval1)!=2){
      stop("rangeval of mu must be the same as rangeval of data.")
    }
    if(is.null(dx)){
      dx <- (rangeval.mu[2]-rangeval.mu[1])*0.01
    }
    abscissa <- seq(rangeval.mu[1],rangeval.mu[2],by=dx)
    mu.eval <- t(fda::eval.fd(fdobj=mu,evalarg=abscissa))
  }else if(is.vector(mu)){
    mu.eval <- mu
  }else{
    stop("Second argument must be either a functional data object or a numeric vector.")
  }
  
  
  n1 <- dim(coeff1)[1]
  n2 <- dim(coeff2)[1]
  J <- dim(coeff1)[2]
  n <- n1+n2
  etichetta_ord <- c(rep(1,n1),rep(2,n2))
  #coeff1 <- coeff1 - matrix(data=mu,nrow=n1,ncol=J)
  
  #print('First step: basis expansion')
  #splines coefficients:
  eval <- coeff <- rbind(coeff1,coeff2)
  p <- dim(coeff)[2]
  
  data.eval <- eval
  #data.eval[1:n1,] <- data.eval[1:n1,] + matrix(data=mu,nrow=n1,ncol=J)
  
  #print('Point-wise tests')
  #univariate permutations
  T0 <- (colMeans(coeff[1:n1,,drop=F]) - colMeans(coeff[(n1+1):n,,drop=F]))^2 #sample mean difference
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
    T_coeff[perm,] <- (colMeans(coeff_perm[1:n1,,drop=F]) - colMeans(coeff_perm[(n1+1):n,,drop=F]))^2
  }
  pval <- numeric(p)
  for(i in 1:p){
    pval[i] <- sum(T_coeff[,i]>=T0[i])/B
  }
  
  #combination
  #print('Partition tests')
  partition = factor(partition)
  nintervals = length(levels(partition))
  ntests = 2^nintervals - 1
  all_combs <- matrix(nrow=ntests,ncol=p)
  labels = levels(partition)
  tt = 1
  for(nint in 1:nintervals){
    combinations = utils::combn(labels,nint)
    n.comb = dim(combinations)[2]
    for(comb in 1:n.comb){
      index = rep(0,p)
      for(ii in 1:dim(combinations)[1]){
        index = index + as.numeric(partition==combinations[ii,comb])
      }
      all_combs[tt,] = index
      tt = tt + 1
    }
  }
  
  #interval-wise tests
  adjusted.pval <- numeric(p)
  responsible.test <- matrix(nrow=p,ncol=p)
  for(test in 1:ntests){
    T0_comb <- sum(T0[which(all_combs[test,]==1)])
    T_comb <- (rowSums(T_coeff[,which(all_combs[test,]==1),drop=FALSE]))
    pval.temp <- mean(T_comb>=T0_comb)
    indexes <- which(all_combs[test,]==1)
    max <- (apply(rbind(adjusted.pval[indexes],pval.temp),2,which.max) )
    adjusted.pval[indexes] <- apply(rbind(adjusted.pval[indexes],pval.temp),2,max)
    if(2 %in% max){
      responsible.test[indexes[which(max==2)] , ] <- matrix(data=all_combs[test,],nrow=sum((max==2)), ncol=p,byrow=TRUE)
    }
  }
  
  
  result = list(
    test = '2pop', mu = mu.eval,
    adjusted_pval = adjusted.pval,
    unadjusted_pval = pval,
    data.eval=data.eval,
    ord_labels = etichetta_ord
  )
  return(result)
}
