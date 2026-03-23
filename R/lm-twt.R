#' Threshold-wise testing procedure for testing functional-on-scalar linear
#' models
#'
#' The function is used to fit and test functional linear models. It can be used
#' to carry out regression, and analysis of variance. It implements the
#' Threshold-wise testing procedure (TWT) for testing the significance of the
#' effects of scalar covariates on a functional population.
#'
#' @inherit functional_lm_test params return seealso
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
#' # Performing the TWT
#' TWT_result <- TWTlm(temperature ~ groups, B = 100L)
#' # Summary of the TWT results
#' summary(TWT_result)
#'
#' # Plot of the TWT results
#' plot(
#'   TWT_result,
#'   main = 'NASA data',
#'   plot_adjpval = TRUE,
#'   xlab = 'Day',
#'   xrange = c(1, 365)
#' )
#'
#' plot(
#'   TWT_result,
#'   main = 'NASA data',
#'   plot_adjpval = TRUE,
#'   xlab = 'Day',
#'   xrange = c(1, 365)
#' )
TWTlm <- function(
  formula,
  dx = NULL,
  B = 1000L,
  method = c("residuals", "responses")
) {
  method <- rlang::arg_match(method)
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
    ) *
      matrix(
        resvar,
        nrow = nvar + 1,
        ncol = p,
        byrow = TRUE
      )
  )
  T0_part <- abs(regr0$coeff / se)^2
  if (nvar > 0) {
    T0_glob <- colSums(
      (regr0$fitted -
        matrix(
          colMeans(regr0$fitted),
          nrow = n,
          ncol = p,
          byrow = TRUE
        ))^2
    ) /
      (nvar * resvar)
  } else {
    method <- "responses"
    T0_glob <- numeric(p)
    # Ensure T0_part is (nvar+1) x p = 1 x p for intercept-only
    T0_part <- matrix(T0_part, nrow = 1L, ncol = p)
  }
  # Compute residuals
  if (method == "residuals") {
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
    mf_temp_raw <- stats::model.frame(formula_temp)[-((p + 1):(p + nvar + 1))]
    mf_temp_cov <- as.data.frame(design_matrix[, -1, drop = FALSE])
    colnames(mf_temp_cov) <- var_names[-1]
    mf_temp <- cbind(mf_temp_raw, mf_temp_cov)
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
      mf_temp2_raw <- stats::model.frame(formula_temp2)[
        -((p + 1):(p + nvar + 1))
      ]
      mf_temp2_cov <- as.data.frame(design_matrix_names2[, -1, drop = FALSE])
      colnames(mf_temp2_cov) <- var_names2[-1]
      mf_temp2 <- cbind(mf_temp2_raw, mf_temp2_cov)
      formula_coeff_temp <- paste(coeffnames, '~', formula_temp)
      formula_coeff_part[[ii]] <- sapply(formula_coeff_temp, stats::as.formula)
      regr0_part[[ii]] <- lapply(
        formula_coeff_part[[ii]],
        stats::lm,
        data = mf_temp2
      )
      residui[ii, , ] <- simplify2array(lapply(
        regr0_part[[ii]],
        extract_residuals
      ))
      fitted_part[ii, , ] <- simplify2array(lapply(
        regr0_part[[ii]],
        extract_fitted
      ))
    }
    ii <- 1 # intercept
    formula_temp <- paste(formula_const, ' -1', sep = '')
    formula_coeff_temp <- paste(coeffnames, '~', formula_temp)
    formula_coeff_part[[ii]] <- sapply(formula_coeff_temp, stats::as.formula)
    regr0_part[[ii]] <- lapply(
      formula_coeff_part[[ii]],
      stats::lm,
      data = mf_temp
    )
    residui[ii, , ] <- simplify2array(lapply(
      regr0_part[[ii]],
      extract_residuals
    ))
    fitted_part[ii, , ] <- simplify2array(lapply(
      regr0_part[[ii]],
      extract_fitted
    ))
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
    } else {
      # Test on intercept permuting signs
      signs <- stats::rbinom(n, 1, 0.5) * 2 - 1
      coeff_perm <- coeff * signs
    }
    regr_perm <- stats::lm.fit(design_matrix, coeff_perm)
    Sigma <- chol2inv(regr_perm$qr$qr)
    resvar <- colSums(regr_perm$residuals^2) / regr_perm$df.residual
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
        (nvar * resvar)
    }
    # Partial tests: differ depending on the method
    if (method == "responses") {
      se <- sqrt(
        matrix(
          diag(Sigma),
          nrow = nvar + 1,
          ncol = p,
          byrow = FALSE
        ) *
          matrix(
            resvar,
            nrow = nvar + 1,
            ncol = p,
            byrow = TRUE
          )
      )
      T_part[perm, , ] <- abs(regr0$coeff / se)^2
    } else if (method == "residuals") {
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
            nrow = nvar + 1,
            ncol = p,
            byrow = FALSE
          ) *
            matrix(
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
    pval_part[, i] <- colSums(
      T_part[,, i] >=
        matrix(
          T0_part[, i],
          nrow = B,
          ncol = nvar + 1,
          byrow = TRUE
        )
    ) /
      B
  }

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
    T_comb <- rowSums(T_glob[, points_2, drop = FALSE], na.rm = TRUE)
    pval_test <- mean(T_comb >= T0_comb)
    pval_tmp[points_2] <- pval_test
    # compute maximum
    adjusted_pval_glob <- apply(rbind(adjusted_pval_glob, pval_tmp), 2, max)
  }

  # F-tests on single factors
  thresholds <- c(0, sort(unique(as.numeric(pval_part))), 1)
  adjusted_pval_part <- pval_part # we initialize the adjusted p-value as unadjusted one

  for (ii in 1:(nvar + 1)) {
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

  fitted_regr <- regr0$fitted
  fitted_t <- fitted_regr

  rownames(adjusted_pval_part) <- var_names
  rownames(coeff_t) <- var_names
  rownames(coeff_regr) <- var_names
  rownames(pval_part) <- var_names

  data_eval <- coeff
  residuals_t <- data_eval - fitted_t
  ybar_t <- colMeans(data_eval)
  npt <- p
  R2_t = colSums(
    (fitted_t -
      matrix(
        data = ybar_t,
        nrow = n,
        ncol = npt,
        byrow = TRUE
      ))^2
  ) /
    colSums(
      (data_eval -
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
    unadjusted_pval_part = pval_part,
    adjusted_pval_part = adjusted_pval_part,
    data_eval = coeff,
    coeff_regr_eval = coeff_t,
    fitted_eval = fitted_t,
    residuals_eval = residuals_t,
    R2_eval = R2_t
  )
  class(out) <- "flm"
  out
}
