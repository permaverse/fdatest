# Threshold wise testing function
# arguments:
# data1 and data2: functional data sets of the two groups. They can be either data of class fd or matrices with pointwise evaluations. 
# mu: mean difference under the null hypothesis
# B: number of permutations
# paired: whether the test is paired or not
# dx: optional. If data1 and data2 are fd objects, dx is used as grid step to evaluate them
TWT2 <- function(data1,data2,mu=0,B=1000,paired=FALSE,dx=NULL){
  
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
  # at the end you have two matrices coeff1 and coeff2 with the point wise evaluation of the functional data pver a grid
  
  n1 <- dim(coeff1)[1]
  n2 <- dim(coeff2)[1]
  n <- n1+n2
  coeff <- rbind(coeff1,coeff2)
  p <- dim(coeff)[2]
  etichetta_ord <- c(rep(1,n1),rep(2,n2))
  #print('Point-wise tests')
  # First part:
  # univariate permutation test for each point
  # this is the computation that needs to be done for each voxel
  T0 <- (colMeans(coeff[1:n1,,drop=F],na.rm=TRUE) - colMeans(coeff[(n1+1):n,,drop=F],na.rm=TRUE))^2 #sample mean difference
  T_coeff <- matrix(ncol=p,nrow=B)
  for (perm in 1:B){ # loop on random permutations
    if(paired==TRUE){ # paired test (for brain data we will not need it)
      if.perm <- stats::rbinom(n1,1,0.5) 
      coeff_perm <- coeff
      for(couple in 1:n1){
        if(if.perm[couple]==1){
          coeff_perm[c(couple,n1+couple),] <- coeff[c(n1+couple,couple),]
        }
      }
    }else if(paired==FALSE){ # unpaired test
      permutazioni <- sample(n)
      coeff_perm <- coeff[permutazioni,]
    }
    T_coeff[perm,] <- (colMeans(coeff_perm[1:n1,,drop=F],na.rm=TRUE) - colMeans(coeff_perm[(n1+1):n,,drop=F],na.rm=TRUE))^2
  }
  
  # p-value computation
  pval <- numeric(p)
  for(i in 1:p){ 
    pval[i] <- sum(T_coeff[,i]>=T0[i])/B
  }
  
  # Second part:
  # combination into subsets
  #print('Interval-wise tests')
  
  thresholds = c(0,sort(unique(pval)),1)
  adjusted.pval <- pval # we initialize the adjusted p-value as unadjusted one
  pval.tmp <- rep(0,p) # inizialize p-value vector resulting from combined test
  for(test in 1:length(thresholds)){
    #print(paste(test,length(thresholds)))
    # test below threshold
    points.1 = which(pval <= thresholds[test])
    T0_comb = sum(T0[points.1],na.rm=TRUE) # combined test statistic
    T_comb <- (rowSums(T_coeff[,points.1,drop=FALSE],na.rm=TRUE))
    pval.test <- mean(T_comb>=T0_comb)
    pval.tmp[points.1] <- pval.test
    # compute maximum
    adjusted.pval = apply(rbind(adjusted.pval,pval.tmp),2,max) 
    
    # test above threshold
    points.2 = which(pval > thresholds[test])
    T0_comb = sum(T0[points.2]) # combined test statistic
    T_comb <- (rowSums(T_coeff[,points.2,drop=FALSE],na.rm=TRUE))
    pval.test <- mean(T_comb>=T0_comb)
    pval.tmp[points.2] <- pval.test
    # compute maximum
    adjusted.pval = apply(rbind(adjusted.pval,pval.tmp),2,max) 
    
  }
  
  result = list(
    test = '2pop', 
    mu = mu.eval,
    adjusted_pval = adjusted.pval,
    unadjusted_pval = pval,
    data.eval = coeff,
    ord_labels = etichetta_ord
  )
  return(result)
}

