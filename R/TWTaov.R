#' Threshold Wise Testing procedure for testing functional analysis of variance
#'
#' The function implements the Threshold Wise Testing procedure for testing mean
#' differences between several functional populations in a one-way or multi-way
#' functional analysis of variance framework. Functional data are tested locally
#' and unadjusted and adjusted p-value functions are provided. The unadjusted
#' p-value function controls the point-wise error rate. The adjusted p-value
#' function controls the threshold-wise error rate.
#'
#' @param formula An object of class "\code{\link{formula}}" (or one that can be
#'   coerced to that class): a symbolic description of the model to be fitted.
#'   The output variable of the formula can be either a matrix of dimension
#'   \code{c(n,J)} collecting the pointwise evaluations of \code{n} functional
#'   data on the same grid of \code{J} points, or a \code{fd} object from the
#'   package \code{fda}.
#' @param B The number of iterations of the MC algorithm to evaluate the
#'   p-values of the permutation tests. The defualt is \code{B=1000}.
#' @param method Permutation method used to calculate the p-value of permutation
#'   tests. Choose "\code{residuals}" for the permutations of residuals under
#'   the reduced model, according to the Freedman and Lane scheme, and
#'   "\code{responses}" for the permutation of the responses, according to the
#'   Manly scheme.
#' @param dx Used only if a \code{fd} object is provided. In this case,
#'   \code{dx} is the size of the discretization step of the grid  used to
#'   evaluate functional data. If set to \code{NULL}, a grid of size 100 is
#'   used. Default is \code{NULL}.
#'
#' @returns An object of class `IWTaov`. The function \code{summary} is used to
#'   obtain and print a summary of the results. An object of class
#'   "\code{IWTaov}" is a list containing at least the following components:
#'   
#'   - `call`: The matched call.
#'   - `design_matrix`: The design matrix of the functional-on-scalar linear
#'   model.
#'   - `unadjusted_pval_F`: Evaluation on a grid of the unadjusted p-value
#'   function of the functional F-test.
#'   - `adjusted_pval_F`: Evaluation on a grid of the adjusted p-value function
#'   of the functional F-test.
#'   - `unadjusted_pval_factors`: Evaluation on a grid of the unadjusted p-value
#'   function of the functional F-tests on each factor of the analysis of
#'   variance (rows).
#'   - `adjusted_pval_factors`: Adjusted p-values of the functional F-tests on
#'   each factor of the analysis of variance (rows) and each basis coefficient
#'   (columns).
#'   - `data.eval`: Evaluation on a fine uniform grid of the functional data
#'   obtained through the basis expansion.
#'   - `coeff.regr.eval`: Evaluation on a fine uniform grid of the functional
#'   regression coefficients.
#'   - `fitted.eval`: Evaluation on a fine uniform grid of the fitted values of
#'   the functional regression.
#'   - `residuals.eval`: Evaluation on a fine uniform grid of the residuals of
#'   the functional regression.
#'   - `R2.eval`: Evaluation on a fine uniform grid of the functional R-squared
#'   of the regression.
#'
#' @seealso See \code{\link{summary.IWTaov}} for summaries and
#'   \code{\link{plot.IWTaov}} for plotting the results. See
#'   \code{\link{ITPaovbspline}} for a functional analysis of variance test
#'   based on B-spline basis expansion. See also \code{\link{IWTlm}} to fit and
#'   test a functional-on-scalar linear model applying the IWT, and
#'   \code{\link{IWT1}}, \code{\link{IWT2}}  for one-population and
#'   two-population tests.
#'
#' @references
#' Abramowicz, K., Pini, A., Schelin, L., Stamm, A., & Vantini, S. (2022).
#' “Domain selection and familywise error rate for functional data: A unified
#' framework. \emph{Biometrics} 79(2), 1119-1132.
#'
#' D. Freedman and D. Lane (1983). A Nonstochastic Interpretation of Reported
#' Significance Levels. \emph{Journal of Business & Economic Statistics} 1.4,
#' 292-298.
#'
#' B. F. J. Manly (2006). Randomization, \emph{Bootstrap and Monte Carlo Methods
#' in Biology}. Vol. 70. CRC Press.
#'
#' @export
#' @examples
#' # Importing the NASA temperatures data set
#' data(NASAtemp)
#' temperature <- rbind(NASAtemp$milan, NASAtemp$paris)
#' groups <- c(rep(0, 22), rep(1, 22))
#'
#' # Performing the TWT
#' TWT.result <- TWTaov(temperature ~ groups, B = 100L)
#'
#' # Summary of the ITP results
#' summary(TWT.result)
#'
#' # Plot of the TWT results
#' layout(1)
#' plot(TWT.result)
#'
#' # All graphics on the same device
#' layout(matrix(1:4, nrow = 2, byrow = FALSE))
#' plot(
#'   TWT.result, 
#'   main = 'NASA data', 
#'   plot_adjpval = TRUE, 
#'   xlab = 'Day', 
#'   xrange = c(1, 365)
#' )
TWTaov <- function(formula,
                   B = 1000L,
                   method = "residuals",
                   dx = NULL) {
  cl <- match.call()
  coeff <- formula2coeff(formula, dx = dx)
  
  dummynames.all <- colnames(attr(stats::terms(formula), "factors"))
  formula.const <- deparse(formula[[3]], width.cutoff = 500L) #extracting the part after ~ on formula. this will not work if the formula is longer than 500 char
  
  formula.discrete <- stats::as.formula(
    paste('coeff ~', formula.const), 
    env = environment()
  )
  design.matrix <- stats::model.matrix(formula.discrete)
  mf <- stats::model.frame(formula.discrete)
  
  n <- dim(coeff)[1]
  J <- dim(coeff)[2]
  
  p <- dim(coeff)[2]
  npt <- J
  
  print('Point-wise tests')
  #univariate permutations
  coeffnames <- paste('coeff[,', as.character(1:p), ']', sep = '')
  formula.coeff <- paste(coeffnames, '~', formula.const)
  formula.coeff <- sapply(formula.coeff, stats::as.formula, env = environment())
  
  aovcoeff1 <- stats::aov(formula.coeff[[1]], data = mf)
  var.names <- rownames(summary(aovcoeff1)[[1]])
  df.vars <- summary(aovcoeff1)[[1]][, 1]
  df.residuals <- df.vars[length(df.vars)]
  var.names <- var.names[-length(var.names)]
  nvar <- length(var.names)
  for (ii in 1:nvar) {
    var.names[ii] <- gsub(' ' , '', var.names[ii])
  }
  
  index.vars <- cbind(
    c(2, (cumsum(df.vars) + 2)[-length(df.vars)]), 
    cumsum(df.vars) + 1
  )
  regr0 <- stats::lm.fit(design.matrix, coeff)
  
  MS0 <- matrix(nrow = nvar + 1, ncol = p)
  for (var in 1:(nvar + 1)) {
    MS0[var, ] <- colSums(rbind(
      regr0$effects[index.vars[var, 1]:index.vars[var, 2], ]^2
    )) / df.vars[var]
  }
  
  # test statistic:
  T0_part <- MS0[1:nvar, ] / matrix(
    MS0[nvar + 1, ],
    nrow = nvar,
    ncol = p,
    byrow = TRUE
  )
  Sigma <- chol2inv(regr0$qr$qr)
  resvar <- colSums(regr0$residuals^2) / regr0$df.residual
  
  if (nvar > 1) {
    T0_glob <- colSums((regr0$fitted - matrix(
      colMeans(regr0$fitted),
      nrow = n,
      ncol = p,
      byrow = TRUE
    ))^2) / ((nvar) * resvar)
  } else if (nvar == 1) { #only one factor -> the permutation of the residuals is equivalent to the one of responses
    method <- "responses"
    T0_glob <- colSums((regr0$fitted - matrix(
      colMeans(regr0$fitted),
      nrow = n,
      ncol = p,
      byrow = TRUE
    ))^2) / ((nvar) * resvar)
  } else if (nvar == 0) {
    method <- "responses" # model with only intercept -> the permutation of the residuals is equivalent to the one of responses
    T0_glob <- numeric(p)
  }
  
  #calculate residuals
  if (method == "residuals") {
    #n residuals for each coefficient of basis expansion (1:p)
    #and for each partial test + global test (nvar+1)
    #saved in array of dim (nvar+1,n,p)
    design.matrix.names2 <- design.matrix
    var.names2 <- var.names
    if (length(grep('factor', formula.const)) > 0) {
      index.factor <- grep('factor', var.names)
      replace.names <- paste('group', (1:length(index.factor)), sep = '')
      var.names2[index.factor] <- replace.names
      colnames(design.matrix.names2) <- var.names2
    }
    
    residui <- array(dim = c(nvar, n, p))
    fitted_part <- array(dim = c(nvar, n, p)) # fitted values of the reduced model (different for each test)
    formula.coeff_part <- vector('list', nvar)
    regr0_part <- vector('list', nvar)
    dummy.interaz <- grep(':', dummynames.all)
    for (ii in 1:nvar) { #no test on intercept
      var.ii <- var.names2[ii]
      variables.reduced <- var.names2[-which(var.names2 == var.ii)] #removing the current variable to test
      
      if (length(grep(':', var.ii)) > 0) { # testing interaction
        var12 <- strsplit(var.ii, ':')
        var1 <- var12[[1]][1]
        var2 <- var12[[1]][2]
        dummy.test1 <- grep(var1, dummynames.all)
        dummy.test2 <- grep(var2, dummynames.all)
        dummy.test <- intersect(dummy.test1, dummy.test2)
        dummynames.reduced <- dummynames.all[-dummy.test]
      } else {
        dummy.test <- grep(var.ii, dummynames.all)
        dummy.test <- setdiff(dummy.test, dummy.interaz)
        dummynames.reduced <- dummynames.all[-dummy.test]
      }
      
      if (nvar > 1) {
        formula.temp <- paste(dummynames.reduced, collapse = ' + ')
      } else {
        formula.temp <- '1' #removing the only variable -> reduced model only has intercept term
      }
      
      formula.coeff.temp <- paste(coeffnames, '~', formula.temp)
      formula.coeff_part[[ii]] <- sapply(
        formula.coeff.temp,
        stats::as.formula,
        env = environment()
      )
      regr0_part[[ii]] <- lapply(formula.coeff_part[[ii]], stats::lm)
      
      residui[ii, , ] <- simplify2array(lapply(regr0_part[[ii]], extract_residuals))
      fitted_part[ii, , ] <- simplify2array(lapply(regr0_part[[ii]], extract_fitted))
    }
  }
  
  T_glob <- matrix(ncol = p, nrow = B)
  T_part <- array(dim = c(B, nvar, p))
  
  for (perm in 1:B) {
    # the F test is the same for both methods
    if (nvar > 0) {
      permutazioni <- sample(n)
      coeff_perm <- coeff[permutazioni, ]
    } else { # testing intercept -> permute signs
      signs <- stats::rbinom(n, 1, 0.5) * 2 - 1
      coeff_perm <- coeff * signs
    }
    
    regr_perm <- stats::lm.fit(design.matrix, coeff_perm)
    Sigma <- chol2inv(regr_perm$qr$qr)
    resvar <- colSums(regr_perm$residuals^2) / regr0$df.residual
    
    if (nvar > 0)
      T_glob[perm, ] <- colSums((
        regr_perm$fitted - matrix(
          colMeans(regr_perm$fitted),
          nrow = n,
          ncol = p,
          byrow = TRUE
        )
      )^2) / ((nvar) * resvar)
    
    # partial tests: differ depending on the method
    if (method == "responses") {
      MSperm <- matrix(nrow = nvar + 1, ncol = p)
      for (var in 1:(nvar + 1)) {
        MSperm[var, ] <- colSums(rbind(
          regr_perm$effects[index.vars[var, 1]:index.vars[var, 2], ]^2
        )) / df.vars[var]
      }
      # test statistic:
      T_part[perm, , ] <- MSperm[1:nvar, ] / matrix(
        MSperm[nvar + 1, ],
        nrow = nvar,
        ncol = p,
        byrow = TRUE
      )
    } else if (method == "residuals") {
      residui_perm <- residui[, permutazioni, ]
      aov_perm_part <- vector('list', nvar)
      for (ii in 1:nvar) {
        coeff_perm <- fitted_part[ii, , ] + residui_perm[ii, , ]
        regr_perm <- stats::lm.fit(design.matrix, coeff_perm)
        MSperm <- matrix(nrow = nvar + 1, ncol = p)
        for (var in 1:(nvar + 1)) {
          MSperm[var, ] <- colSums(rbind(
            regr_perm$effects[index.vars[var, 1]:index.vars[var, 2], ]^2
          )) / df.vars[var]
        }
        # test statistic:
        T_part[perm, ii, ] <- (MSperm[1:nvar, ] / matrix(
          MSperm[nvar + 1, ],
          nrow = nvar,
          ncol = p,
          byrow = TRUE
        ))[ii, ]
      }
    }
  }
  
  pval_glob <- numeric(p)
  pval_part <- matrix(nrow = nvar, ncol = p)
  for (i in 1:p) {
    pval_glob[i] <- sum(T_glob[, i] >= T0_glob[i]) / B
    pval_part[, i] <- colSums(T_part[, , i] >= matrix(
      T0_part[, i],
      nrow = B,
      ncol = nvar,
      byrow = TRUE
    )) / B
  }
  
  #combination
  print('Threshold-wise tests')
  
  # F-test
  thresholds <- c(0, sort(unique(pval_glob)), 1)
  adjusted.pval_glob <- pval_glob # we initialize the adjusted p-value as unadjusted one
  pval.tmp <- rep(0, p) # inizialize p-value vector resulting from combined test
  for (test in 1:length(thresholds)) {
    # test below threshold
    points.1 <- which(pval_glob <= thresholds[test])
    T0_comb <- sum(T0_glob[points.1], na.rm = TRUE) # combined test statistic
    T_comb <- (rowSums(T_glob[, points.1, drop = FALSE], na.rm = TRUE))
    pval.test <- mean(T_comb >= T0_comb)
    pval.tmp[points.1] <- pval.test
    # compute maximum
    adjusted.pval_glob <- apply(rbind(adjusted.pval_glob, pval.tmp), 2, max)
    
    # test above threshold
    points.2 <- which(pval_glob > thresholds[test])
    T0_comb <- sum(T0_glob[points.2]) # combined test statistic
    T_comb <- (rowSums(T_glob[, points.2, drop = FALSE], na.rm = TRUE))
    pval.test <- mean(T_comb >= T0_comb)
    pval.tmp[points.2] <- pval.test
    # compute maximum
    adjusted.pval_glob <- apply(rbind(adjusted.pval_glob, pval.tmp), 2, max)
  }
  
  # F-tests on single factors
  thresholds <- c(0, sort(unique(as.numeric(pval_part))), 1)
  adjusted.pval_part <- pval_part # we initialize the adjusted p-value as unadjusted one
  
  for (ii in 1:nvar) {
    pval.tmp <- rep(0, p)
    for (test in 1:length(thresholds)) {
      # test below threshold
      points.1 <- which(pval_part[ii, ] <= thresholds[test])
      T0_comb <- sum(T0_part[ii, points.1], na.rm = TRUE) # combined test statistic
      T_comb <- rowSums(T_part[, ii, points.1, drop = FALSE], na.rm = TRUE)
      pval.test <- mean(T_comb >= T0_comb)
      pval.tmp[points.1] <- pval.test
      # compute maximum
      adjusted.pval_part[ii, ] <- apply(rbind(adjusted.pval_part[ii, ], pval.tmp), 2, max)
      
      # test above threshold
      points.2 <- which(pval_part[ii, ] > thresholds[test])
      T0_comb <- sum(T0_part[ii, points.2]) # combined test statistic
      T_comb <- rowSums(T_part[, ii, points.2, drop = FALSE], na.rm = TRUE)
      pval.test <- mean(T_comb >= T0_comb)
      pval.tmp[points.2] <- pval.test
      # compute maximum
      adjusted.pval_part[ii, ] <- apply(rbind(adjusted.pval_part[ii, ], pval.tmp), 2, max) 
    }
  }
  
  coeff.regr <- regr0$coeff
  coeff.t <- coeff.regr
  
  fitted.regr <- regr0$fitted.values
  fitted.t <- fitted.regr
  
  rownames(adjusted.pval_part) <- var.names
  rownames(coeff.t) <- colnames(design.matrix)
  rownames(coeff.regr) <- colnames(design.matrix)
  rownames(pval_part) <- var.names
  
  residuals.t <- coeff - fitted.t
  ybar.t <- colMeans(coeff)
  R2.t <- colSums((fitted.t - matrix(
    data = ybar.t,
    nrow = n,
    ncol = npt,
    byrow = TRUE
  ))^2) / colSums((coeff - matrix(
    data = ybar.t,
    nrow = n,
    ncol = npt,
    byrow = TRUE
  ))^2)
  
  print('Threshold-Wise Testing completed')
  
  out <- list(
    call = cl,
    design_matrix = design.matrix,
    unadjusted_pval_F = pval_glob,
    adjusted_pval_F = adjusted.pval_glob,
    unadjusted_pval_factors = pval_part,
    adjusted_pval_factors = adjusted.pval_part,
    data.eval = coeff,
    coeff.regr.eval = coeff.t,
    fitted.eval = fitted.t,
    residuals.eval = residuals.t,
    R2.eval = R2.t
  )
  class(out) <- "IWTaov"
  out
}
