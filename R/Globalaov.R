#' Global testing procedure for testing functional analysis of variance
#'
#' The function implements the Global Testing procedure for testing mean
#' differences between several functional populations in a one-way or multi-way
#' functional analysis of variance framework. Functional data are tested
#' globally and unadjusted and adjusted p-value functions are provided. The
#' unadjusted p-value function controls the point-wise error rate. The adjusted
#' p-value function controls the family-wise error rate weakly. Since this is a
#' global test, the adjusted p-value function is constant.
#'
#' @param formula An object of class "\code{\link{formula}}" (or one that can be
#'   coerced to that class): a symbolic description of the model to be fitted.
#'   The output variable of the formula can be either a matrix of dimension
#'   \code{c(n,J)} collecting the pointwise evaluations of \code{n} functional
#'   data on the same grid of \code{J} points, or a \code{fd} object from the
#'   package \code{fda}.
#' @inheritParams IWT1
#' @param method Permutation method used to calculate the p-value of permutation
#'   tests. Choose "\code{residuals}" for the permutations of residuals under
#'   the reduced model, according to the Freedman and Lane scheme, and
#'   "\code{responses}" for the permutation of the responses, according to the
#'   Manly scheme.
#' @param stat Type of test statistic used for the global test. Possible values
#'   are: \code{'Integral'} (default) for the integral over the domain of the
#'   F-test statistic; \code{'Max'} for max over the domain of the F-test
#'   statistic.
#'
#' @return An object of class `IWTaov`. The function \code{summary} is used to
#'   obtain and print a summary of the results. This object is a list containing
#'   the following components:
#'   
#'   - `call`: The matched call.
#'   - `design_matrix`: The design matrix of the functional-on-scalar linear
#'   model.
#'   - `unadjusted_pval_F`: Evaluation on a grid of the unadjusted p-value
#'   function of the functional F-test.
#'   - `adjusted_pval_F`: Evaluation on a grid of the adjusted p-value function
#'   of the functional F-test.
#'   - `unadjusted_pval_factors`: Evaluation on a grid of the unadjusted
#'   p-value function of the functional F-tests on each factor of the analysis
#'   of variance (rows).
#'   - `adjusted.pval.factors`: Adjusted p-values of the functional F-tests on
#'   each factor of the analysis of variance (rows) and each basis coefficient
#'   (columns).
#'   - `Global_pval_F`: Global p-value of the overall test F.
#'   - `Global_pval_factors`: Global p-value of test F involving each factor
#'   separately.
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
#' - Abramowicz, K., Pini, A., Schelin, L., Stamm, A., & Vantini, S. (2022).
#' “Domain selection and familywise error rate for functional data: A unified
#' framework. \emph{Biometrics} 79(2), 1119-1132.
#' - D. Freedman and D. Lane (1983). A Nonstochastic Interpretation of Reported
#' Significance Levels. \emph{Journal of Business & Economic Statistics} 1.4,
#' 292-298.
#' - B. F. J. Manly (2006). Randomization, \emph{Bootstrap and Monte Carlo
#' Methods in Biology}. Vol. 70. CRC Press.
#'
#' @export
#' @examples
#' # Importing the NASA temperatures data set
#' data(NASAtemp)
#' temperature <- rbind(NASAtemp$milan, NASAtemp$paris)
#' groups <- c(rep(0, 22), rep(1, 22))
#'
#' # Performing the test
#' Global.result <- Globalaov(temperature ~ groups, B = 1000)
#'
#' # Summary of the test results
#' summary(Global.result)
#'
#' # Plot of the results
#' layout(1)
#' plot(Global.result)
#'
#' # All graphics on the same device
#' layout(matrix(1:4, nrow = 2, byrow = FALSE))
#' plot(
#'   Global.result, 
#'   main = 'NASA data', 
#'   plot.adjpval = TRUE, 
#'   xlab = 'Day', 
#'   xrange = c(1, 365)
#' )
Globalaov <- function(formula, 
                      B = 1000, 
                      dx = NULL, 
                      recycle = TRUE, 
                      method = 'residuals', 
                      stat = 'Integral') {
  env <- environment(formula)
  variables <- all.vars(formula)
  y.name <- variables[1]
  covariates.names <- colnames(attr(stats::terms(formula), "factors"))
  cl <- match.call()
  data <- get(y.name, envir = env)
  if (fda::is.fd(data)) { # data is a functional data object
    rangeval <- data$basis$rangeval
    if (is.null(dx)) {
      dx <- (rangeval[2] - rangeval[1]) * 0.01
    }
    abscissa <- seq(rangeval[1], rangeval[2], by = dx)
    coeff <- t(fda::eval.fd(fdobj = data, evalarg = abscissa))
  } else if (is.matrix(data)) {
    coeff <- data
  } else {
    stop("First argument of the formula must be either a functional data object or a matrix.")
  }
  
  possible_statistics <- c("Integral", "Max")
  if (!(stat %in% possible_statistics)) {
    stop(paste0(
      'Possible statistics are ',
      paste0(possible_statistics, collapse = ', ')
    ))
  }
  
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
    method <- 'responses'
    T0_glob <- colSums((regr0$fitted - matrix(
      colMeans(regr0$fitted),
      nrow = n,
      ncol = p,
      byrow = TRUE
    ))^2) / ((nvar) * resvar)
  } else if (nvar == 0) {
    method <- 'responses' # model with only intercept -> the permutation of the residuals is equivalent to the one of responses
    T0_glob <- numeric(p)
  }
  
  #calculate residuals
  if (method == 'residuals') {
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
    for (ii in 1:(nvar)) { #no test on intercept
      var.ii <- var.names2[ii]
      variables.reduced <- var.names2[-which(var.names2 == var.ii)] #removing the current variable to test
      
      if (length(grep(':',var.ii)) > 0) { # testing interaction
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
        FUN = stats::as.formula, 
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
    if (nvar >0) {
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
    if (method == 'responses') {
      MSperm <- matrix(nrow = nvar + 1, ncol = p)
      for (var in 1:(nvar+1)) {
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
    } else if (method == 'residuals') {
      residui_perm <- residui[, permutazioni, ]
      aov_perm_part <- vector('list', nvar)
      for (ii in 1:nvar) {
        coeff_perm <- fitted_part[ii, , ] + residui_perm[ii, , ]
        regr_perm <- stats::lm.fit(design.matrix, coeff_perm)
        MSperm <- matrix(nrow = nvar + 1, ncol = p)
        for (var in 1:(nvar+1)) {
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
  print('Global test')
  if (stat == 'Integral') {
    T0_temp <- sum(T0_glob[1:p])
    T_temp <- rowSums(T_glob[, 1:p])
    Global_pval_F <- sum(T_temp >= T0_temp) / B
    
    Global_pval_factors <- numeric(nvar)
    for (ii in 1:nvar) {
      T0_temp <- sum(T0_part[ii, 1:p])
      T_temp <- rowSums(T_part[, ii, 1:p])
      Global_pval_factors[ii] <- sum(T_temp >= T0_temp) / B
    }
    
    corrected.pval_glob <- rep(Global_pval_F, p)
    
    corrected.pval_part <- matrix(nrow = nvar, ncol = p)
    for (ii in 1:nvar) {
      corrected.pval_part[ii, ] <- rep(Global_pval_factors[ii], p)
    }
  } else if (stat=='Max') {
    T0_temp <- max(T0_glob)
    T_temp <- apply(T_glob, 1, max)
    Global_pval_F <- sum(T_temp >= T0_temp) / B
    
    Global_pval_factors <- numeric(nvar)
    for (ii in 1:nvar) {
      T0_temp <- max(T0_part[ii, ])
      T_temp <- apply(T_part[, ii, ], 1, max)
      Global_pval_factors[ii] <- sum(T_temp >= T0_temp) / B
    }
    
    corrected.pval_glob <- rep(Global_pval_F, p)
    
    corrected.pval_part <- matrix(nrow = nvar, ncol = p)
    for (ii in 1:nvar) {
      corrected.pval_part[ii, ] <- rep(Global_pval_factors[ii], p)
    }
  }
  
  coeff.regr <- regr0$coeff
  coeff.t <- coeff.regr
  
  fitted.regr <- regr0$fitted.values
  fitted.t <- fitted.regr
  
  rownames(corrected.pval_part) <- var.names
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
  
  print('Interval-Wise Testing completed')
  
  out <- list(
    call = cl,
    design_matrix = design.matrix,
    unadjusted_pval_F = pval_glob,
    adjusted_pval_F = corrected.pval_glob,
    unadjusted_pval_factors = pval_part,
    adjusted_pval_factors = corrected.pval_part,
    Global_pval_F = Global_pval_F,
    Global_pval_factors = Global_pval_factors,
    data.eval = coeff,
    coeff.regr.eval = coeff.t,
    fitted.eval = fitted.t,
    residuals.eval = residuals.t,
    R2.eval = R2.t
  )
  class(out) <- 'IWTaov'
  out
}
