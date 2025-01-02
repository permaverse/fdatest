AVAILABLE_ALTERNATIVES <- function() {
  c("two.sided", "less", "greater")
}

AVAILABLE_METHODS <- function() {
  c("IWT", "TWT", "PCT", "Global", "FDR")
}

stat_lm_glob <- function(anova) {
  stats::summary.lm(anova)$f[1]
}

stat_aov_part <- function(anova) {
  out <- summary(anova)[[1]][, 4]
  out <- out[-length(out)]
  out
}

extract_residuals <- function(x) {
  x$residuals
}

extract_fitted <- function(x) {
  x$fitted
}

pval_correct <- function(pval.matrix) {
  matrice_pval_2_2x <- cbind(pval.matrix, pval.matrix)
  p <- dim(pval.matrix)[2]
  matrice_pval_2_2x <- matrice_pval_2_2x[, (2 * p):1]
  corrected.pval.matrix <- matrix(nrow = p, ncol = p)
  corrected.pval.matrix[p, ] <- pval.matrix[p, p:1]
  for (var in 1:p) {
    pval_var <- matrice_pval_2_2x[p, var]
    inizio <- var
    fine <- var #inizio fisso, fine aumenta salendo nelle righe
    for (riga in (p - 1):1) {
      fine <- fine + 1
      pval_cono <- matrice_pval_2_2x[riga, inizio:fine]
      pval_var <- max(pval_var, pval_cono, na.rm = TRUE)
      corrected.pval.matrix[riga, var] <- pval_var
    }
  }
  corrected.pval.matrix[, p:1]
}

onesample2coeffs <- function(data, mu, dx = NULL) {
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
    stop("data argument must be either a functional data object or a matrix.")
  }
  
  if (fda::is.fd(mu)) { # mu is a functional data
    rangeval.mu <- mu$basis$rangeval
    if (sum(rangeval.mu == rangeval) != 2) {
      stop("rangeval of mu must be the same as rangeval of data.")
    }
    if (is.null(dx)) {
      dx <- (rangeval.mu[2] - rangeval.mu[1]) * 0.01
    }
    abscissa <- seq(rangeval.mu[1], rangeval.mu[2], by = dx)
    mu.eval <- t(fda::eval.fd(fdobj = mu, evalarg = abscissa))
  } else if (is.vector(mu)) {
    mu.eval <- mu
  } else {
    stop("mu argument must be either a functional data object or a numeric vector.")
  }
  
  list(coeff = coeff, mu = mu.eval)
}

twosamples2coeffs <- function(data1, data2, mu, dx = NULL) {
  if (fda::is.fd(data1) && fda::is.fd(data2)) {
    rangeval1 <- data1$basis$rangeval
    rangeval2 <- data2$basis$rangeval
    if (sum(rangeval1 == rangeval2) != 2) {
      stop("rangeval of data1 must be the same as rangeval of data2.")
    }
    if (is.null(dx)) {
      dx <- (rangeval1[2] - rangeval1[1]) * 0.01
    }
    abscissa <- seq(rangeval1[1], rangeval1[2], by = dx)
    coeff1 <- t(fda::eval.fd(fdobj = data1, evalarg = abscissa))
    coeff2 <- t(fda::eval.fd(fdobj = data2, evalarg = abscissa))
  } else if (is.matrix(data1) && is.matrix(data2)) {
    coeff1 <- data1
    coeff2 <- data2
  } else {
    stop("Both arguments must be either functional data objects or matrices.")
  }
  
  if (fda::is.fd(mu)) { # mu is a functional data
    rangeval.mu <- mu$basis$rangeval
    if (sum(rangeval.mu == rangeval1) != 2) {
      stop("rangeval of mu must be the same as rangeval of data.")
    }
    if (is.null(dx)) {
      dx <- (rangeval.mu[2] - rangeval.mu[1]) * 0.01
    }
    abscissa <- seq(rangeval.mu[1], rangeval.mu[2], by = dx)
    mu.eval <- t(fda::eval.fd(fdobj = mu, evalarg = abscissa))
  } else if (is.vector(mu)) {
    mu.eval <- mu
  } else {
    stop("mu must be either a functional data object or a numeric vector.")
  }
  
  list(coeff1 = coeff1, coeff2 = coeff2, mu = mu.eval)
}

formula2coeff <- function(formula, dx = NULL) {
  env <- environment(formula)
  variables <- all.vars(formula)
  y.name <- variables[1]
  covariates.names <- colnames(attr(stats::terms(formula), "factors"))
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
  
  coeff
}

formula2design_matrix <- function(formula, coeff) {
  # extracting the part after ~ on formula. this will not work if the formula is
  # longer than 500 char
  formula.const <- deparse(formula[[3]], width.cutoff = 500L) 
  formula.discrete <- stats::as.formula(
    paste('coeff ~', formula.const), 
    env = environment()
  )
  stats::model.matrix(formula.discrete)
}
