#' Global testing procedure for testing functional-on-scalar linear models
#'
#' The function is used to fit and test functional linear models. It can be used
#' to carry out regression, and analysis of variance. It implements the global
#' testing procedure for testing the significance of the effects of scalar
#' covariates on a functional population.
#'
#' @inheritParams Globalaov
#'
#' @returns An object of class `IWTlm`. The function \code{summary} is used to
#'   obtain and print a summary of the results. This object is a list containing
#'   the following components:
#'   
#'   - `call`: Call of the function.
#'   - `design_matrix`: Design matrix of the linear model.
#'   - `unadjusted_pval_F`: Unadjusted p-value function of the F test.
#'   - `adjusted_pval_F`: Adjusted p-value function of the F test.
#'   - `unadjusted_pval_part`: Unadjusted p-value functions of the functional
#'   t-tests on each covariate, separately (rows) on each domain point
#'   (columns).
#'   - `adjusted_pval_part`: Adjusted p-values of the functional t-tests on each
#'   covariate (rows) on each domain point (columns).
#'   - `Global_pval_F`: Global p-value of the overall test F.
#'   - `Global_pval_part`: Global p-value of t-test involving each covariate
#'   separately.
#'   - `data.eval`: Evaluation of functional data.
#'   - `coeff.regr.eval`: Evaluation of the regression coefficients.
#'   - `fitted.eval`: Evaluation of the fitted values.
#'   - `residuals.eval`: Evaluation of the residuals.
#'   - `R2.eval`: Evaluation of the functional R-suared.
#'
#' @seealso See \code{\link{summary.IWTlm}} for summaries and
#'   \code{\link{plot.IWTlm}} for plotting the results. See
#'   \code{\link{ITPlmbspline}} for a functional linear model test based on an
#'   a-priori selected B-spline basis expansion. See also \code{\link{IWTaov}}
#'   to fit and test a functional analysis of variance applying the IWT, and
#'   \code{\link{IWT1}}, \code{\link{IWT2}} for one-population and
#'   two-population tests.
#'
#' @references
#' Abramowicz, K., Pini, A., Schelin, L., Stamm, A., & Vantini, S. (2022).
#' “Domain selection and familywise error rate for functional data: A unified
#' framework. \emph{Biometrics} 79(2), 1119-1132.
#'
#' D. Freedman and D. Lane (1983). A Nonstochastic Interpretation of Reported
#' Significance Levels. \emph{Journal of Business & Economic Statistics} 1(4),
#' 292-298.
#'
#' B. F. J. Manly (2006). Randomization, \emph{Bootstrap and Monte Carlo Methods
#' in Biology}. Vol. 70. CRC Press.
#'
#' @export
#' @examples
#' # Defining the covariates
#' temperature <- rbind(NASAtemp$milan, NASAtemp$paris)
#' groups <- c(rep(0, 22), rep(1, 22))
#'
#' # Performing the IWT
#' Global.result <- Globallm(temperature ~ groups, B = 1000)
#' # Summary of the IWT results
#' summary(Global.result)
#'
#' # Plot of the IWT results
#' layout(1)
#' plot(
#'   Global.result, 
#'   main = 'NASA data', 
#'   plot_adjpval = TRUE, 
#'   xlab = 'Day', 
#'   xrange = c(1, 365)
#' )
#'
#' # All graphics on the same device
#' layout(matrix(1:6, nrow = 3, byrow = FALSE))
#' plot(
#'   Global.result, 
#'   main = 'NASA data', 
#'   plot_adjpval = TRUE, 
#'   xlab = 'Day', 
#'   xrange = c(1, 365)
#' )
Globallm <- function(formula,
                     dx = NULL,
                     B = 1000L,
                     method = c("residuals", "responses"),
                     stat = c("Integral", "Max")) {
  method <- rlang::arg_match(method)
  stat <- rlang::arg_match(stat)
  cl <- match.call()
  coeff <- formula2coeff(formula, dx = dx)
  design_matrix <- formula2design_matrix(formula, coeff)
  
  nvar <- dim(design_matrix)[2] - 1
  var_names <- colnames(design_matrix)
  p <- dim(coeff)[2]
  n <- dim(coeff)[1]
  # Univariate permutations
  regr0 <- stats::lm.fit(design_matrix, coeff)
  # Test statistics
  Sigma <- chol2inv(regr0$qr$qr)
  resvar <- colSums(regr0$residuals^2) / regr0$df.residual
  se <- sqrt(
    matrix(
      diag(Sigma),
      nrow = nvar + 1,
      ncol = p,
      byrow = FALSE
    ) * matrix(
      resvar,
      nrow = nvar + 1,
      ncol = p,
      byrow = TRUE
    )
  )
  T0_part <- abs(regr0$coeff / se)^2
  if (nvar > 0) {
    T0_glob <- colSums((regr0$fitted - matrix(
      colMeans(regr0$fitted),
      nrow = n,
      ncol = p,
      byrow = TRUE
    ))^2) / (nvar * resvar)
  } else {
    method <- 'responses'
    T0_glob <- numeric(p)
    T0_part <- t(as.matrix(T0_part))
  }
  # Compute residuals
  if (method == 'residuals') {
    # n residuals for each coefficient of basis expansion (1:p) 
    # and for each partial test + global test (nvar+1) 
    # Saved in array of dim (nvar+1,n,p)
    # Extracting the part after ~ on formula. 
    # This will not work if the formula 
    # is longer than 500 char
    formula_const <- deparse(formula[[3]], width.cutoff = 500L)
    design_matrix_names2 <- design_matrix
    var_names2 <- var_names
    coeffnames <- paste('coeff[,', as.character(1:p), ']', sep = '')
    formula_temp <- coeff ~ design_matrix
    mf_temp <- cbind(
      stats::model.frame(formula_temp)[-((p + 1):(p + nvar + 1))], 
      as.data.frame(design_matrix[, -1])
    )
    if (length(grep('factor', formula_const)) > 0) {
      index_factor <- grep('factor', var_names)
      replace_names <- paste('group', (1:length(index_factor)), sep = '')
      var_names2[index_factor] <- replace_names
      colnames(design_matrix_names2) <- var_names2
    }
    residui <- array(dim = c(nvar + 1, n, p))
    fitted_part <- array(dim = c(nvar + 1, n, p)) 
    formula_coeff_part <- vector('list', nvar + 1)
    regr0_part <- vector('list', nvar + 1)
    # The first one is the intercept. Treated as special case after loop
    for (ii in 2:(nvar + 1)) { 
      var_ii <- var_names2[ii]
      variables_reduced <- var_names2[-c(1, which(var_names2 == var_ii))]
      if (nvar > 1) {
        formula_temp <- paste(variables_reduced, collapse = ' + ')
      } else {
        # Removing the unique variable -> reduced model only has intercept ter
        formula_temp <- '1' 
      }
      formula_temp2 <- coeff ~ design_matrix_names2
      mf_temp2 <- cbind(
        stats::model.frame(formula_temp2)[-((p + 1):(p + nvar + 1))], 
        as.data.frame(design_matrix_names2[, -1])
      )
      formula_coeff_temp <- paste(coeffnames, '~', formula_temp)
      formula_coeff_part[[ii]] <- sapply(formula_coeff_temp, stats::as.formula)
      regr0_part[[ii]] <- lapply(formula_coeff_part[[ii]], stats::lm, data = mf_temp2)
      residui[ii, , ] <- simplify2array(lapply(regr0_part[[ii]], extract_residuals))
      fitted_part[ii, , ] <- simplify2array(lapply(regr0_part[[ii]], extract_fitted))
    }
    ii <- 1 # intercept
    formula_temp <- paste(formula_const, ' -1', sep = '')
    formula_coeff_temp <- paste(coeffnames, '~', formula_temp)
    formula_coeff_part[[ii]] <- sapply(formula_coeff_temp, stats::as.formula)
    regr0_part[[ii]] <- lapply(formula_coeff_part[[ii]], stats::lm, data = mf_temp)
    residui[ii, , ] <- simplify2array(lapply(regr0_part[[ii]], extract_residuals))
    fitted_part[ii, , ] <- simplify2array(lapply(regr0_part[[ii]], extract_fitted))
  }
  
  cli::cli_h1("Point-wise tests")
  
  # CMC algorithm
  T_glob <- matrix(ncol = p, nrow = B)
  T_part <- array(dim = c(B, nvar + 1, p))
  for (perm in 1:B) {
    # the F test is the same for both methods
    if (nvar > 0) {
      permutazioni <- sample(n)
      coeff_perm <- coeff[permutazioni, ]
    } else { # Test on intercept permuting signs
      signs <- stats::rbinom(n, 1, 0.5) * 2 - 1
      coeff_perm <- coeff * signs
    }
    regr_perm <- stats::lm.fit(design_matrix, coeff_perm)
    Sigma <- chol2inv(regr_perm$qr$qr)
    resvar <- colSums(regr_perm$residuals^2) / regr_perm$df.residual
    if (nvar > 0) {
      T_glob[perm, ] <- colSums((
        regr_perm$fitted - matrix(
          colMeans(regr_perm$fitted),
          nrow = n,
          ncol = p,
          byrow = TRUE
        )
      )^2) / (nvar * resvar)
    }
    # Partial tests: differ depending on the method
    if (method == 'responses') {
      se <- sqrt(
        matrix(
          diag(Sigma),
          nrow = nvar + 1,
          ncol = p,
          byrow = FALSE
        ) * matrix(
          resvar,
          nrow = nvar + 1,
          ncol = p,
          byrow = TRUE
        )
      )
      T_part[perm, , ] <- abs(regr0$coeff / se)^2
    } else if (method == 'residuals') {
      residui_perm <- residui[, permutazioni, ]
      regr_perm_part <- vector('list', nvar + 1)
      for (ii in 1:(nvar + 1)) {
        coeff_perm <- fitted_part[ii, , ] + residui_perm[ii, , ]
        regr_perm <- stats::lm.fit(design_matrix, coeff_perm)
        Sigma <- chol2inv(regr_perm$qr$qr)
        resvar <- colSums(regr_perm$residuals^2) / regr_perm$df.residual
        se <- sqrt(
          matrix(
            diag(Sigma),
            nrow = nvar + 1 ,
            ncol = p,
            byrow = FALSE
          ) * matrix(
            resvar,
            nrow = nvar + 1,
            ncol = p,
            byrow = TRUE
          )
        )
        T_part[perm, ii, ] <- abs(regr_perm$coeff / se)[ii, ]^2
      }
    }
  }
  pval_glob <- numeric(p)
  pval_part <- matrix(nrow = nvar + 1, ncol = p)
  for (i in 1:p) {
    pval_glob[i] <- sum(T_glob[, i] >= T0_glob[i]) / B
    pval_part[, i] <- colSums(T_part[, , i] >= matrix(
      T0_part[, i],
      nrow = B,
      ncol = nvar + 1,
      byrow = TRUE
    )) / B
  }
  
  cli::cli_h1("Global test")
  
  if (stat == "Integral") {
    T0_temp <- sum(T0_glob[1:p])
    T_temp <- rowSums(T_glob[, 1:p])
    Global_pval_F <- sum(T_temp >= T0_temp) / B
    
    Global_pval_part <- numeric(nvar + 1)
    for (ii in 1:(nvar + 1)) {
      T0_temp <- sum(T0_part[ii, 1:p])
      T_temp <- rowSums(T_part[, ii, 1:p])
      Global_pval_part[ii] <- sum(T_temp >= T0_temp) / B
    }
    
    corrected.pval_glob <- rep(Global_pval_F, p)
    
    corrected.pval_part <- matrix(nrow = nvar + 1, ncol = p)
    for (ii in 1:(nvar + 1)) {
      corrected.pval_part[ii, ] <- rep(Global_pval_part[ii], p)
    }
  } else if (stat == 'Max') {
    T0_temp <- max(T0_glob)
    T_temp <- apply(T_glob, 1, max)
    Global_pval_F <- sum(T_temp >= T0_temp) / B
    
    Global_pval_part <- numeric(nvar + 1)
    for (ii in 1:(nvar + 1)) {
      T0_temp <- max(T0_part[ii, ])
      T_temp <- apply(T_part[, ii, ], 1, max)
      Global_pval_part[ii] <- sum(T_temp >= T0_temp) / B
    }
    
    corrected.pval_glob <- rep(Global_pval_F, p)
    corrected.pval_part <- matrix(nrow = nvar + 1, ncol = p)
    for (ii in 1:(nvar + 1)) {
      corrected.pval_part[ii, ] <- rep(Global_pval_part[ii], p)
    }
  }
  
  coeff.regr <- regr0$coeff
  coeff.t <- coeff.regr
  
  fitted.regr <- regr0$fitted
  fitted.t <- fitted.regr
  
  rownames(corrected.pval_part) <- var_names
  rownames(coeff.t) <- var_names
  rownames(coeff.regr) <- var_names
  rownames(pval_part) <- var_names
  
  data.eval <- coeff
  residuals.t <- data.eval - fitted.t
  ybar.t <- colMeans(data.eval)
  npt <- p
  R2.t <- colSums((fitted.t - matrix(
    data = ybar.t,
    nrow = n,
    ncol = npt,
    byrow = TRUE
  ))^2) / colSums((data.eval - matrix(
    data = ybar.t,
    nrow = n,
    ncol = npt,
    byrow = TRUE
  ))^2)
  
  cli::cli_h1("Global Testing completed")
  
  out <- list(
    call = cl,
    design_matrix = design_matrix,
    unadjusted_pval_F = pval_glob,
    adjusted_pval_F = corrected.pval_glob,
    unadjusted_pval_part = pval_part,
    adjusted_pval_part = corrected.pval_part,
    Global_pval_F = Global_pval_F,
    Global_pval_part = Global_pval_part,
    data.eval = coeff,
    coeff.regr.eval = coeff.t,
    fitted.eval = fitted.t,
    residuals.eval = residuals.t,
    R2.eval = R2.t
  )
  class(out) <- 'IWTlm'
  out
}
