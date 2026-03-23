#' Threshold Wise Testing procedure for testing functional analysis of variance
#'
#' The function implements the Threshold Wise Testing procedure for testing mean
#' differences between several functional populations in a one-way or multi-way
#' functional analysis of variance framework. Functional data are tested locally
#' and unadjusted and adjusted p-value functions are provided. The unadjusted
#' p-value function controls the point-wise error rate. The adjusted p-value
#' function controls the threshold-wise error rate.
#'
#' @inherit functional_anova_test params return seealso
#'
#' @references
#' Abramowicz, K., Pini, A., Schelin, L., Stamm, A., & Vantini, S. (2022).
#' “Domain selection and familywise error rate for functional data: A unified
#' framework. *Biometrics* 79(2), 1119-1132.
#'
#' D. Freedman and D. Lane (1983). A Nonstochastic Interpretation of Reported
#' Significance Levels. *Journal of Business & Economic Statistics* 1.4,
#' 292-298.
#'
#' B. F. J. Manly (2006). Randomization, *Bootstrap and Monte Carlo Methods
#' in Biology*. Vol. 70. CRC Press.
#'
#' @export
#' @examples
#' temperature <- rbind(NASAtemp$milan, NASAtemp$paris)
#' groups <- c(rep(0, 22), rep(1, 22))
#'
#' # Performing the TWT
#' TWT_result <- TWTaov(temperature ~ groups, B = 100L)
#'
#' # Summary of the TWT results
#' summary(TWT_result)
#'
#' # Plot of the TWT results
#' layout(1)
#' plot(TWT_result)
#'
#' # All graphics on the same device
#' layout(matrix(1:4, nrow = 2, byrow = FALSE))
#' plot(
#'   TWT_result,
#'   main = 'NASA data',
#'   plot_adjpval = TRUE,
#'   xlab = 'Day',
#'   xrange = c(1, 365)
#' )
TWTaov <- function(
  formula,
  dx = NULL,
  B = 1000L,
  method = c("residuals", "responses")
) {
  method <- rlang::arg_match(method)
  cl <- match.call()
  coeff <- formula2coeff(formula, dx = dx)

  dummynames_all <- colnames(attr(stats::terms(formula), "factors"))
  formula_const <- deparse(formula[[3]], width.cutoff = 500L) #extracting the part after ~ on formula. this will not work if the formula is longer than 500 char

  formula_discrete <- stats::as.formula(
    paste('coeff ~', formula_const),
    env = environment()
  )
  design_matrix <- stats::model.matrix(formula_discrete)
  mf <- stats::model.frame(formula_discrete)

  n <- dim(coeff)[1]
  J <- dim(coeff)[2]

  p <- dim(coeff)[2]
  npt <- J

  cli::cli_h1("Point-wise tests")

  #univariate permutations
  coeffnames <- paste('coeff[,', as.character(1:p), ']', sep = '')
  formula_coeff <- paste(coeffnames, '~', formula_const)
  formula_coeff <- sapply(formula_coeff, stats::as.formula, env = environment())

  aovcoeff1 <- stats::aov(formula_coeff[[1]], data = mf)
  var_names <- rownames(summary(aovcoeff1)[[1]])
  df_vars <- summary(aovcoeff1)[[1]][, 1]
  df_residuals <- df_vars[length(df_vars)]
  var_names <- var_names[-length(var_names)]
  nvar <- length(var_names)
  for (ii in 1:nvar) {
    var_names[ii] <- gsub(' ', '', var_names[ii])
  }

  index_vars <- cbind(
    c(2, (cumsum(df_vars) + 2)[-length(df_vars)]),
    cumsum(df_vars) + 1
  )
  regr0 <- stats::lm.fit(design_matrix, coeff)

  MS0 <- matrix(nrow = nvar + 1, ncol = p)
  for (var in 1:(nvar + 1)) {
    MS0[var, ] <- colSums(rbind(
      regr0$effects[index_vars[var, 1]:index_vars[var, 2], ]^2
    )) /
      df_vars[var]
  }

  # test statistic:
  T0_part <- MS0[1:nvar, ] /
    matrix(
      MS0[nvar + 1, ],
      nrow = nvar,
      ncol = p,
      byrow = TRUE
    )
  Sigma <- chol2inv(regr0$qr$qr)
  resvar <- colSums(regr0$residuals^2) / regr0$df.residual

  if (nvar > 1) {
    T0_glob <- colSums(
      (regr0$fitted -
        matrix(
          colMeans(regr0$fitted),
          nrow = n,
          ncol = p,
          byrow = TRUE
        ))^2
    ) /
      ((nvar) * resvar)
  } else if (nvar == 1) {
    #only one factor -> the permutation of the residuals is equivalent to the one of responses
    method <- "responses"
    T0_glob <- colSums(
      (regr0$fitted -
        matrix(
          colMeans(regr0$fitted),
          nrow = n,
          ncol = p,
          byrow = TRUE
        ))^2
    ) /
      ((nvar) * resvar)
  } else if (nvar == 0) {
    method <- "responses" # model with only intercept -> the permutation of the residuals is equivalent to the one of responses
    T0_glob <- numeric(p)
  }

  #calculate residuals
  if (method == "residuals") {
    #n residuals for each coefficient of basis expansion (1:p)
    #and for each partial test + global test (nvar+1)
    #saved in array of dim (nvar+1,n,p)
    design_matrix_names2 <- design_matrix
    var_names2 <- var_names
    if (length(grep('factor', formula_const)) > 0) {
      index_factor <- grep('factor', var_names)
      replace_names <- paste('group', (1:length(index_factor)), sep = '')
      var_names2[index_factor] <- replace_names
      colnames(design_matrix_names2) <- var_names2
    }

    residui <- array(dim = c(nvar, n, p))
    fitted_part <- array(dim = c(nvar, n, p)) # fitted values of the reduced model (different for each test)
    formula_coeff_part <- vector('list', nvar)
    regr0_part <- vector('list', nvar)
    dummy_interaz <- grep(':', dummynames_all)
    for (ii in 1:nvar) {
      #no test on intercept
      var_ii <- var_names2[ii]
      variables_reduced <- var_names2[-which(var_names2 == var_ii)] #removing the current variable to test

      if (length(grep(':', var_ii)) > 0) {
        # testing interaction
        var12 <- strsplit(var_ii, ':')
        var1 <- var12[[1]][1]
        var2 <- var12[[1]][2]
        dummy_test1 <- grep(var1, dummynames_all)
        dummy_test2 <- grep(var2, dummynames_all)
        dummy_test <- intersect(dummy_test1, dummy_test2)
        dummynames_reduced <- dummynames_all[-dummy_test]
      } else {
        dummy_test <- grep(var_ii, dummynames_all)
        dummy_test <- setdiff(dummy_test, dummy_interaz)
        dummynames_reduced <- dummynames_all[-dummy_test]
      }

      if (nvar > 1) {
        formula_temp <- paste(dummynames_reduced, collapse = ' + ')
      } else {
        formula_temp <- '1' #removing the only variable -> reduced model only has intercept term
      }

      formula_coeff_temp <- paste(coeffnames, '~', formula_temp)
      formula_coeff_part[[ii]] <- sapply(
        formula_coeff_temp,
        stats::as.formula,
        env = environment()
      )
      regr0_part[[ii]] <- lapply(formula_coeff_part[[ii]], stats::lm)

      residui[ii, , ] <- simplify2array(lapply(
        regr0_part[[ii]],
        extract_residuals
      ))
      fitted_part[ii, , ] <- simplify2array(lapply(
        regr0_part[[ii]],
        extract_fitted
      ))
    }
  }

  T_glob <- matrix(ncol = p, nrow = B)
  T_part <- array(dim = c(B, nvar, p))

  for (perm in 1:B) {
    # the F test is the same for both methods
    if (nvar > 0) {
      permutazioni <- sample(n)
      coeff_perm <- coeff[permutazioni, ]
    } else {
      # testing intercept -> permute signs
      signs <- stats::rbinom(n, 1, 0.5) * 2 - 1
      coeff_perm <- coeff * signs
    }

    regr_perm <- stats::lm.fit(design_matrix, coeff_perm)
    Sigma <- chol2inv(regr_perm$qr$qr)
    resvar <- colSums(regr_perm$residuals^2) / regr0$df.residual

    if (nvar > 0) {
      T_glob[perm, ] <- colSums(
        (regr_perm$fitted -
          matrix(
            colMeans(regr_perm$fitted),
            nrow = n,
            ncol = p,
            byrow = TRUE
          ))^2
      ) /
        ((nvar) * resvar)
    }

    # partial tests: differ depending on the method
    if (method == "responses") {
      MSperm <- matrix(nrow = nvar + 1, ncol = p)
      for (var in 1:(nvar + 1)) {
        MSperm[var, ] <- colSums(rbind(
          regr_perm$effects[index_vars[var, 1]:index_vars[var, 2], ]^2
        )) /
          df_vars[var]
      }
      # test statistic:
      T_part[perm, , ] <- MSperm[1:nvar, ] /
        matrix(
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
        regr_perm <- stats::lm.fit(design_matrix, coeff_perm)
        MSperm <- matrix(nrow = nvar + 1, ncol = p)
        for (var in 1:(nvar + 1)) {
          MSperm[var, ] <- colSums(rbind(
            regr_perm$effects[index_vars[var, 1]:index_vars[var, 2], ]^2
          )) /
            df_vars[var]
        }
        # test statistic:
        T_part[perm, ii, ] <- (MSperm[1:nvar, ] /
          matrix(
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
    pval_part[, i] <- colSums(
      T_part[,, i] >=
        matrix(
          T0_part[, i],
          nrow = B,
          ncol = nvar,
          byrow = TRUE
        )
    ) /
      B
  }

  #combination
  cli::cli_h1("Threshold-wise tests")

  # F-test
  thresholds <- c(0, sort(unique(pval_glob)), 1)
  adjusted_pval_glob <- pval_glob # we initialize the adjusted p-value as unadjusted one
  pval_tmp <- rep(0, p) # inizialize p-value vector resulting from combined test
  for (test in 1:length(thresholds)) {
    # test below threshold
    points_1 <- which(pval_glob <= thresholds[test])
    T0_comb <- sum(T0_glob[points_1], na.rm = TRUE) # combined test statistic
    T_comb <- (rowSums(T_glob[, points_1, drop = FALSE], na.rm = TRUE))
    pval_test <- mean(T_comb >= T0_comb)
    pval_tmp[points_1] <- pval_test
    # compute maximum
    adjusted_pval_glob <- apply(rbind(adjusted_pval_glob, pval_tmp), 2, max)

    # test above threshold
    points_2 <- which(pval_glob > thresholds[test])
    T0_comb <- sum(T0_glob[points_2]) # combined test statistic
    T_comb <- (rowSums(T_glob[, points_2, drop = FALSE], na.rm = TRUE))
    pval_test <- mean(T_comb >= T0_comb)
    pval_tmp[points_2] <- pval_test
    # compute maximum
    adjusted_pval_glob <- apply(rbind(adjusted_pval_glob, pval_tmp), 2, max)
  }

  # F-tests on single factors
  thresholds <- c(0, sort(unique(as.numeric(pval_part))), 1)
  adjusted_pval_part <- pval_part # we initialize the adjusted p-value as unadjusted one

  for (ii in 1:nvar) {
    pval_tmp <- rep(0, p)
    for (test in 1:length(thresholds)) {
      # test below threshold
      points_1 <- which(pval_part[ii, ] <= thresholds[test])
      T0_comb <- sum(T0_part[ii, points_1], na.rm = TRUE) # combined test statistic
      T_comb <- rowSums(T_part[, ii, points_1, drop = FALSE], na.rm = TRUE)
      pval_test <- mean(T_comb >= T0_comb)
      pval_tmp[points_1] <- pval_test
      # compute maximum
      adjusted_pval_part[ii, ] <- apply(
        rbind(adjusted_pval_part[ii, ], pval_tmp),
        2,
        max
      )

      # test above threshold
      points_2 <- which(pval_part[ii, ] > thresholds[test])
      T0_comb <- sum(T0_part[ii, points_2]) # combined test statistic
      T_comb <- rowSums(T_part[, ii, points_2, drop = FALSE], na.rm = TRUE)
      pval_test <- mean(T_comb >= T0_comb)
      pval_tmp[points_2] <- pval_test
      # compute maximum
      adjusted_pval_part[ii, ] <- apply(
        rbind(adjusted_pval_part[ii, ], pval_tmp),
        2,
        max
      )
    }
  }

  coeff_regr <- regr0$coeff
  coeff_t <- coeff_regr

  fitted_regr <- regr0$fitted.values
  fitted_t <- fitted_regr

  rownames(adjusted_pval_part) <- var_names
  rownames(coeff_t) <- colnames(design_matrix)
  rownames(coeff_regr) <- colnames(design_matrix)
  rownames(pval_part) <- var_names

  residuals_t <- coeff - fitted_t
  ybar_t <- colMeans(coeff)
  R2_t <- colSums(
    (fitted_t -
      matrix(
        data = ybar_t,
        nrow = n,
        ncol = npt,
        byrow = TRUE
      ))^2
  ) /
    colSums(
      (coeff -
        matrix(
          data = ybar_t,
          nrow = n,
          ncol = npt,
          byrow = TRUE
        ))^2
    )

  cli::cli_h1("Threshold-Wise Testing completed")

  out <- list(
    call = cl,
    design_matrix = design_matrix,
    unadjusted_pval_F = pval_glob,
    adjusted_pval_F = adjusted_pval_glob,
    unadjusted_pval_factors = pval_part,
    adjusted_pval_factors = adjusted_pval_part,
    data_eval = coeff,
    coeff_regr_eval = coeff_t,
    fitted_eval = fitted_t,
    residuals_eval = residuals_t,
    R2_eval = R2_t
  )
  class(out) <- "fanova"
  out
}
