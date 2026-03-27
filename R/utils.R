available_alternatives <- function() {
  c("two.sided", "less", "greater")
}

available_methods <- function() {
  c("IWT", "TWT", "PCT", "Global", "FDR")
}

available_statistics <- function() {
  c("Integral", "Max", "Integral_std", "Max_std")
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

pval_correct <- function(pval_matrix) {
  p <- dim(pval_matrix)[2]
  matrice_pval_2_2x <- cbind(pval_matrix, pval_matrix)[, (2 * p):1]
  corrected_pval_matrix <- matrix(nrow = p, ncol = p)
  corrected_pval_matrix[p, ] <- pval_matrix[p, p:1]
  for (var in seq_len(p)) {
    pval_var <- matrice_pval_2_2x[p, var]
    fine <- var
    for (riga in (p - 1L):1L) {
      fine <- fine + 1L
      pval_var <- max(pval_var, matrice_pval_2_2x[riga, var:fine], na.rm = TRUE)
      corrected_pval_matrix[riga, var] <- pval_var
    }
  }
  corrected_pval_matrix[, p:1]
}

onesample2coeffs <- function(data, mu, dx = NULL) {
  if (fda::is.fd(data)) {
    # data is a functional data object
    rangeval <- data$basis$rangeval
    if (is.null(dx)) {
      dx <- (rangeval[2] - rangeval[1]) * 0.01
    }
    abscissa <- seq(rangeval[1], rangeval[2], by = dx)
    coeff <- t(fda::eval.fd(fdobj = data, evalarg = abscissa))
  } else if (is.matrix(data)) {
    coeff <- data
  } else {
    cli::cli_abort(
      "The {.arg data} argument must be either a functional data object of class
      {.cls fd} or a matrix."
    )
  }

  if (fda::is.fd(mu)) {
    # mu is a functional data
    rangeval_mu <- mu$basis$rangeval
    if (sum(rangeval_mu == rangeval) != 2) {
      cli::cli_abort(
        "The range of values of {.arg mu} must be the same as the range of
        values of {.arg data}."
      )
    }
    abscissa <- seq(rangeval_mu[1], rangeval_mu[2], by = dx)
    mu_eval <- t(fda::eval.fd(fdobj = mu, evalarg = abscissa))
  } else if (is.vector(mu)) {
    mu_eval <- mu
  } else {
    cli::cli_abort(
      "The {.arg mu} argument must be either a functional data object of class
      {.cls fd} or a numeric vector."
    )
  }

  list(coeff = coeff, mu = mu_eval)
}

twosamples2coeffs <- function(data1, data2, mu, dx = NULL) {
  if (fda::is.fd(data1) && fda::is.fd(data2)) {
    rangeval1 <- data1$basis$rangeval
    rangeval2 <- data2$basis$rangeval
    if (sum(rangeval1 == rangeval2) != 2) {
      cli::cli_abort(
        "The range of values of {.arg data1} must be the same as the range of
        values of {.arg data2}."
      )
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
    cli::cli_abort(
      "Both {.arg data1} and {.arg data2} must be either functional data objects
      of class {.cls fd} or matrices."
    )
  }

  if (fda::is.fd(mu)) {
    # mu is a functional data
    rangeval_mu <- mu$basis$rangeval
    if (sum(rangeval_mu == rangeval1) != 2) {
      cli::cli_abort(
        "The range of values of {.arg mu} must be the same as the range of
        values of {.arg data1}."
      )
    }
    abscissa <- seq(rangeval_mu[1], rangeval_mu[2], by = dx)
    mu_eval <- t(fda::eval.fd(fdobj = mu, evalarg = abscissa))
  } else if (is.vector(mu)) {
    mu_eval <- mu
  } else {
    cli::cli_abort(
      "The {.arg mu} argument must be either a functional dataobject of class
      {.cls fd} or a numeric vector."
    )
  }

  list(coeff1 = coeff1, coeff2 = coeff2, mu = mu_eval)
}

formula2coeff <- function(formula, dx = NULL) {
  env <- environment(formula)
  variables <- all.vars(formula)
  y_name <- variables[1]
  data <- get(y_name, envir = env)
  if (fda::is.fd(data)) {
    # data is a functional data object
    rangeval <- data$basis$rangeval
    if (is.null(dx)) {
      dx <- (rangeval[2] - rangeval[1]) * 0.01
    }
    abscissa <- seq(rangeval[1], rangeval[2], by = dx)
    coeff <- t(fda::eval.fd(fdobj = data, evalarg = abscissa))
  } else if (is.matrix(data)) {
    coeff <- data
  } else {
    cli::cli_abort(
      "The first argument of the formula must be either a functional data object
      of class {.cls fd} or a matrix."
    )
  }

  coeff
}

formula2design_matrix <- function(formula, coeff) {
  # extracting the part after ~ on formula. this will not work if the formula is
  # longer than 500 char
  formula_const <- deparse(formula[[3]], width.cutoff = 500L)
  # Parent the evaluation env to the original formula's env so that covariates
  # defined in the caller (e.g. `groups`) can be found by model.matrix().
  env <- new.env(parent = environment(formula))
  env$coeff <- coeff
  formula_discrete <- stats::as.formula(
    paste("coeff ~", formula_const),
    env = env
  )
  stats::model.matrix(formula_discrete)
}

# Internal helper: one permutation iteration for twosample_alt_permtest.
# Returns a numeric vector of length p (the permuted test statistic row).
.twosample_one_perm <- function(coeff, n1, alternative, paired) {
  n <- nrow(coeff)
  if (paired) {
    if_perm <- stats::rbinom(n1, 1, 0.5)
    coeff_perm <- coeff
    for (couple in seq_len(n1)) {
      if (if_perm[couple] == 1L) {
        coeff_perm[c(couple, n1 + couple), ] <- coeff[c(n1 + couple, couple), ]
      }
    }
  } else {
    coeff_perm <- coeff[sample(n), ]
  }
  meandiff_perm <- colMeans(coeff_perm[seq_len(n1), , drop = FALSE]) -
    colMeans(coeff_perm[seq(n1 + 1L, n), , drop = FALSE])
  sign_diff_perm <- sign(meandiff_perm)
  sign_diff_perm[sign_diff_perm == -1L] <- 0L
  switch(
    alternative,
    two.sided = meandiff_perm^2,
    greater = (meandiff_perm * sign_diff_perm)^2,
    less = (meandiff_perm * (sign_diff_perm - 1L))^2
  )
}

# Internal helper: shared pointwise permutation test for two-sample functions
# (iwt2, twt2, fdr2, pct2). Returns list(t0, t_coeff, pval).
twosample_alt_permtest <- function(coeff, n1, n_perm, alternative, paired) {
  n <- nrow(coeff)
  p <- ncol(coeff)

  meandiff <- colMeans(coeff[seq_len(n1), , drop = FALSE], na.rm = TRUE) -
    colMeans(coeff[seq(n1 + 1L, n), , drop = FALSE], na.rm = TRUE)
  sign_diff <- sign(meandiff)
  sign_diff[sign_diff == -1L] <- 0L
  t0 <- switch(
    alternative,
    two.sided = meandiff^2,
    greater = (meandiff * sign_diff)^2,
    less = (meandiff * (sign_diff - 1L))^2
  )

  # Run permutations in parallel via mirai_map().
  # Each task returns one row of t_coeff (length p).
  if (mirai::daemons_set()) {
    perm_args <- list(
      coeff = coeff,
      n1 = n1,
      alternative = alternative,
      paired = paired
    )
    perm_tasks <- mirai::mirai_map(
      seq_len(n_perm),
      function(.x) {
        rlang::inject(.twosample_one_perm(!!!perm_args))
      }
    )
    t_coeff <- do.call(rbind, perm_tasks[])
  } else {
    t_coeff <- matrix(nrow = n_perm, ncol = p)
    for (i in seq_len(n_perm)) {
      t_coeff[i, ] <- .twosample_one_perm(coeff, n1, alternative, paired)
    }
  }

  pval <- colSums(
    t_coeff >= matrix(t0, nrow = n_perm, ncol = p, byrow = TRUE)
  ) /
    n_perm

  list(t0 = t0, t_coeff = t_coeff, pval = pval)
}

# Internal helper: one permutation for aov_permtest.
# Returns list(t_glob_row, t_part_row) where t_part_row is a nvar-length
# numeric vector (one value per variable), already as a (nvar x p) matrix.
.aov_one_perm <- function(
  coeff,
  n,
  p,
  nvar,
  design_matrix,
  index_vars,
  df_vars,
  regr0_df_residual,
  method,
  residui,
  fitted_part
) {
  permutazioni <- sample(n)
  coeff_perm <- coeff[permutazioni, ]

  regr_perm <- stats::lm.fit(design_matrix, coeff_perm)
  resvar <- colSums(regr_perm$residuals^2) / regr0_df_residual

  t_glob_row <- if (nvar > 0) {
    colSums(
      (regr_perm$fitted -
        matrix(
          colMeans(regr_perm$fitted),
          nrow = n,
          ncol = p,
          byrow = TRUE
        ))^2
    ) /
      (nvar * resvar)
  } else {
    numeric(p)
  }

  t_part_row <- matrix(nrow = nvar, ncol = p) # nvar x p

  if (method == "responses") {
    ms_perm <- matrix(nrow = nvar + 1L, ncol = p)
    for (var in seq_len(nvar + 1L)) {
      ms_perm[var, ] <- colSums(rbind(
        regr_perm$effects[index_vars[var, 1]:index_vars[var, 2], ]^2
      )) /
        df_vars[var]
    }
    t_part_row <- ms_perm[seq_len(nvar), , drop = FALSE] /
      matrix(ms_perm[nvar + 1L, ], nrow = nvar, ncol = p, byrow = TRUE)
  } else {
    residui_perm <- residui[, permutazioni, ]
    for (ii in seq_len(nvar)) {
      coeff_perm_ii <- fitted_part[ii, , ] + residui_perm[ii, , ]
      regr_perm_ii <- stats::lm.fit(design_matrix, coeff_perm_ii)
      ms_perm <- matrix(nrow = nvar + 1L, ncol = p)
      for (var in seq_len(nvar + 1L)) {
        ms_perm[var, ] <- colSums(rbind(
          regr_perm_ii$effects[index_vars[var, 1]:index_vars[var, 2], ]^2
        )) /
          df_vars[var]
      }
      ratio_mat <- ms_perm[seq_len(nvar), , drop = FALSE] /
        matrix(ms_perm[nvar + 1L, ], nrow = nvar, ncol = p, byrow = TRUE)
      t_part_row[ii, ] <- ratio_mat[ii, ]
    }
  }

  list(t_glob_row = t_glob_row, t_part_row = t_part_row)
}

# Internal helper: shared pointwise permutation test for AOV functions
# (global_aov, iwt_aov, twt_aov).
# Returns a list with all computed quantities needed by the combination step.
aov_permtest <- function(formula, dx, n_perm, method) {
  coeff <- formula2coeff(formula, dx = dx)

  formula_env <- new.env(parent = environment(formula))
  formula_env$coeff <- coeff

  dummynames_all <- colnames(attr(stats::terms(formula), "factors"))
  formula_const <- deparse(formula[[3]], width.cutoff = 500L)

  formula_discrete <- stats::as.formula(
    paste("coeff ~", formula_const),
    env = formula_env
  )
  design_matrix <- stats::model.matrix(formula_discrete)
  mf <- stats::model.frame(formula_discrete)

  n <- dim(coeff)[1]
  p <- dim(coeff)[2]

  coeffnames <- paste0("coeff[,", as.character(seq_len(p)), "]")
  formula_coeff <- paste(coeffnames, "~", formula_const)
  formula_coeff <- sapply(formula_coeff, stats::as.formula, env = formula_env)

  aovcoeff1 <- stats::aov(formula_coeff[[1]], data = mf)
  var_names <- rownames(summary(aovcoeff1)[[1]])
  df_vars <- summary(aovcoeff1)[[1]][, 1]
  var_names <- var_names[-length(var_names)]
  nvar <- length(var_names)
  for (ii in seq_len(nvar)) {
    var_names[ii] <- gsub(" ", "", var_names[ii])
  }

  index_vars <- cbind(
    c(2, (cumsum(df_vars) + 2)[-length(df_vars)]),
    cumsum(df_vars) + 1
  )
  regr0 <- stats::lm.fit(design_matrix, coeff)

  ms0 <- matrix(nrow = nvar + 1, ncol = p)
  for (var in seq_len(nvar + 1)) {
    ms0[var, ] <- colSums(rbind(
      regr0$effects[index_vars[var, 1]:index_vars[var, 2], ]^2
    )) /
      df_vars[var]
  }
  t0_part <- ms0[seq_len(nvar), , drop = FALSE] /
    matrix(ms0[nvar + 1, ], nrow = nvar, ncol = p, byrow = TRUE)

  resvar <- colSums(regr0$residuals^2) / regr0$df.residual

  if (nvar > 1) {
    mu_fit <- matrix(colMeans(regr0$fitted), nrow = n, ncol = p, byrow = TRUE)
    t0_glob <- colSums((regr0$fitted - mu_fit)^2) / (nvar * resvar)
  } else if (nvar == 1) {
    # only one factor: permuting residuals is equivalent to permuting responses
    method <- "responses"
    mu_fit <- matrix(colMeans(regr0$fitted), nrow = n, ncol = p, byrow = TRUE)
    t0_glob <- colSums((regr0$fitted - mu_fit)^2) / (nvar * resvar)
  } else {
    # intercept only: permute signs
    method <- "responses"
    t0_glob <- numeric(p)
  }

  residui <- fitted_part <- NULL
  if (method == "residuals") {
    var_names2 <- var_names
    if (length(grep("factor", formula_const)) > 0) {
      index_factor <- grep("factor", var_names)
      replace_names <- paste0("group", seq_along(index_factor))
      var_names2[index_factor] <- replace_names
    }

    residui <- array(dim = c(nvar, n, p))
    fitted_part <- array(dim = c(nvar, n, p))
    formula_coeff_part <- vector("list", nvar)
    regr0_part <- vector("list", nvar)
    dummy_interaz <- grep(":", dummynames_all)

    for (ii in seq_len(nvar)) {
      var_ii <- var_names2[ii]
      if (length(grep(":", var_ii)) > 0) {
        var12 <- strsplit(var_ii, ":")
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

      formula_temp <- if (length(dummynames_reduced) > 0) {
        paste(dummynames_reduced, collapse = " + ")
      } else {
        "1"
      }

      formula_coeff_temp <- paste(coeffnames, "~", formula_temp)
      formula_coeff_part[[ii]] <- sapply(
        formula_coeff_temp,
        FUN = stats::as.formula,
        env = formula_env
      )
      regr0_part[[ii]] <- lapply(formula_coeff_part[[ii]], stats::lm)
      residui[ii, , ] <- simplify2array(
        lapply(regr0_part[[ii]], extract_residuals)
      )
      fitted_part[ii, , ] <- simplify2array(
        lapply(regr0_part[[ii]], extract_fitted)
      )
    }
  }

  perm_args <- list(
    coeff = coeff,
    n = n,
    p = p,
    nvar = nvar,
    design_matrix = design_matrix,
    index_vars = index_vars,
    df_vars = df_vars,
    regr0_df_residual = regr0$df.residual,
    method = method,
    residui = residui,
    fitted_part = fitted_part
  )

  # Run permutations in parallel via mirai_map().
  if (mirai::daemons_set()) {
    perm_tasks <- mirai::mirai_map(
      seq_len(n_perm),
      function(.x) {
        rlang::inject(.aov_one_perm(!!!perm_args))
      }
    )
    perm_results <- perm_tasks[.progress]
  } else {
    mirai::daemons(1L, seed = 42)
    seeds <- lapply(seq_len(n_perm), \(.x) mirai::nextstream())
    mirai::daemons(0L)
    perm_results <- lapply(
      seq_len(n_perm),
      function(.x) {
        .Random.seed <<- seeds[[.x]]
        rlang::inject(.aov_one_perm(!!!perm_args))
      }
    )
  }

  t_glob <- do.call(
    rbind,
    lapply(perm_results, function(.x) .x$t_glob_row)
  )
  # t_glob <- do.call(rbind, lapply(perm_results, `[[`, "t_glob_row"))
  # t_part: array dim c(n_perm, nvar, p)
  t_part <- array(dim = c(n_perm, nvar, p), data = NA_real_)
  for (i in seq_len(n_perm)) {
    t_part[i, , ] <- perm_results[[i]]$t_part_row
  }

  pval_glob <- colSums(
    t_glob >= matrix(t0_glob, nrow = n_perm, ncol = p, byrow = TRUE)
  ) /
    n_perm
  pval_part <- matrix(nrow = nvar, ncol = p)
  for (i in seq_len(p)) {
    pval_part[, i] <- colSums(
      t_part[,, i] >=
        matrix(t0_part[, i], nrow = n_perm, ncol = nvar, byrow = TRUE)
    ) /
      n_perm
  }

  list(
    coeff = coeff,
    n = n,
    p = p,
    nvar = nvar,
    var_names = var_names,
    coeffnames = coeffnames,
    design_matrix = design_matrix,
    regr0 = regr0,
    t0_part = t0_part,
    t0_glob = t0_glob,
    t_glob = t_glob,
    t_part = t_part,
    pval_glob = pval_glob,
    pval_part = pval_part,
    method = method
  )
}

# Internal helper: one permutation for lm_permtest.
.lm_one_perm <- function(
  coeff,
  n,
  p,
  nvar,
  design_matrix,
  regr0_df_residual,
  method,
  residui,
  fitted_part
) {
  permutazioni <- sample(n)
  coeff_perm <- coeff[permutazioni, ]

  regr_perm <- stats::lm.fit(design_matrix, coeff_perm)
  sigma_perm <- chol2inv(regr_perm$qr$qr)
  resvar <- colSums(regr_perm$residuals^2) / regr_perm$df.residual

  t_glob_row <- if (nvar > 0) {
    colSums(
      (regr_perm$fitted -
        matrix(colMeans(regr_perm$fitted), nrow = n, ncol = p, byrow = TRUE))^2
    ) /
      (nvar * resvar)
  } else {
    numeric(p)
  }

  # (nvar+1) x p matrix
  t_part_row <- matrix(nrow = nvar + 1L, ncol = p)

  if (method == "responses") {
    se <- sqrt(
      matrix(diag(sigma_perm), nrow = nvar + 1L, ncol = p, byrow = FALSE) *
        matrix(resvar, nrow = nvar + 1L, ncol = p, byrow = TRUE)
    )
    t_part_row <- abs(regr_perm$coeff / se)^2
  } else {
    residui_perm <- residui[, permutazioni, ]
    for (ii in seq_len(nvar + 1L)) {
      coeff_perm_ii <- fitted_part[ii, , ] + residui_perm[ii, , ]
      regr_perm_ii <- stats::lm.fit(design_matrix, coeff_perm_ii)
      sigma_perm_ii <- chol2inv(regr_perm_ii$qr$qr)
      resvar_ii <- colSums(regr_perm_ii$residuals^2) / regr_perm_ii$df.residual
      se_ii <- sqrt(
        matrix(
          diag(sigma_perm_ii),
          nrow = nvar + 1L,
          ncol = p,
          byrow = FALSE
        ) *
          matrix(resvar_ii, nrow = nvar + 1L, ncol = p, byrow = TRUE)
      )
      t_part_row[ii, ] <- abs(regr_perm_ii$coeff / se_ii)[ii, ]^2
    }
  }

  list(t_glob_row = t_glob_row, t_part_row = t_part_row)
}

# Internal helper: shared pointwise permutation test for LM functions
# (global_lm, iwt_lm, twt_lm).
# Returns a list with all computed quantities needed by the combination step.
lm_permtest <- function(formula, dx, n_perm, method) {
  coeff <- formula2coeff(formula, dx = dx)
  design_matrix <- formula2design_matrix(formula, coeff)

  nvar <- dim(design_matrix)[2] - 1
  var_names <- colnames(design_matrix)
  p <- dim(coeff)[2]
  n <- dim(coeff)[1]

  regr0 <- stats::lm.fit(design_matrix, coeff)
  sigma <- chol2inv(regr0$qr$qr)
  resvar <- colSums(regr0$residuals^2) / regr0$df.residual
  se <- sqrt(
    matrix(diag(sigma), nrow = nvar + 1, ncol = p, byrow = FALSE) *
      matrix(resvar, nrow = nvar + 1, ncol = p, byrow = TRUE)
  )
  t0_part <- abs(regr0$coeff / se)^2

  if (nvar > 0) {
    mu_fit <- matrix(colMeans(regr0$fitted), nrow = n, ncol = p, byrow = TRUE)
    t0_glob <- colSums((regr0$fitted - mu_fit)^2) / (nvar * resvar)
  } else {
    method <- "responses"
    t0_glob <- numeric(p)
    t0_part <- matrix(t0_part, nrow = 1L, ncol = p)
  }

  residui <- fitted_part <- NULL
  if (method == "residuals") {
    formula_const <- deparse(formula[[3]], width.cutoff = 500L)
    var_names2 <- var_names
    coeffnames <- paste0("coeff[,", as.character(seq_len(p)), "]")
    formula_temp_dm <- coeff ~ design_matrix
    mf_temp_raw <- stats::model.frame(formula_temp_dm)[
      -((p + 1):(p + nvar + 1))
    ]
    mf_temp_cov <- as.data.frame(design_matrix[, -1, drop = FALSE])
    colnames(mf_temp_cov) <- var_names[-1]
    mf_temp <- cbind(mf_temp_raw, mf_temp_cov)

    design_matrix_names2 <- design_matrix
    if (length(grep("factor", formula_const)) > 0) {
      index_factor <- grep("factor", var_names)
      replace_names <- paste0("group", seq_along(index_factor))
      var_names2[index_factor] <- replace_names
      colnames(design_matrix_names2) <- var_names2
    }

    residui <- array(dim = c(nvar + 1, n, p))
    fitted_part <- array(dim = c(nvar + 1, n, p))
    formula_coeff_part <- vector("list", nvar + 1)
    regr0_part <- vector("list", nvar + 1)

    # Build mf_temp2 once: response columns + renamed predictor columns.
    # Using renamed column names (var_names2) avoids re-evaluating factor()
    # expressions against a data frame that lacks the original variable.
    formula_temp2 <- coeff ~ design_matrix_names2
    mf_temp2_raw <- stats::model.frame(formula_temp2)[
      -((p + 1):(p + nvar + 1))
    ]
    mf_temp2_cov <- as.data.frame(design_matrix_names2[, -1, drop = FALSE])
    colnames(mf_temp2_cov) <- var_names2[-1]
    mf_temp2 <- cbind(mf_temp2_raw, mf_temp2_cov)

    for (ii in seq(2L, nvar + 1L)) {
      var_ii <- var_names2[ii]
      variables_reduced <- var_names2[-c(1L, which(var_names2 == var_ii))]
      formula_temp_str <- if (nvar > 1) {
        paste(variables_reduced, collapse = " + ")
      } else {
        "1"
      }

      formula_coeff_temp <- paste(coeffnames, "~", formula_temp_str)
      formula_coeff_part[[ii]] <- sapply(
        formula_coeff_temp,
        stats::as.formula
      )
      regr0_part[[ii]] <- lapply(
        formula_coeff_part[[ii]],
        stats::lm,
        data = mf_temp2
      )
      residui[ii, , ] <- simplify2array(
        lapply(regr0_part[[ii]], extract_residuals)
      )
      fitted_part[ii, , ] <- simplify2array(
        lapply(regr0_part[[ii]], extract_fitted)
      )
    }

    ii <- 1L # intercept
    # Use renamed predictor column names so that factor() expressions in the
    # original formula are not re-evaluated against design-matrix columns.
    formula_temp_str <- paste0(paste(var_names2[-1], collapse = " + "), " - 1")
    formula_coeff_temp <- paste(coeffnames, "~", formula_temp_str)
    formula_coeff_part[[ii]] <- sapply(formula_coeff_temp, stats::as.formula)
    regr0_part[[ii]] <- lapply(
      formula_coeff_part[[ii]],
      stats::lm,
      data = mf_temp2
    )
    residui[ii, , ] <- simplify2array(
      lapply(regr0_part[[ii]], extract_residuals)
    )
    fitted_part[ii, , ] <- simplify2array(
      lapply(regr0_part[[ii]], extract_fitted)
    )
  }

  perm_args <- list(
    coeff = coeff,
    n = n,
    p = p,
    nvar = nvar,
    design_matrix = design_matrix,
    regr0_df_residual = regr0$df.residual,
    method = method,
    residui = residui,
    fitted_part = fitted_part
  )

  # Run permutations in parallel via mirai_map().
  if (mirai::daemons_set()) {
    perm_tasks <- mirai::mirai_map(
      seq_len(n_perm),
      function(.x) {
        rlang::inject(.lm_one_perm(!!!perm_args))
      }
    )
    perm_results <- perm_tasks[]
  } else {
    perm_results <- lapply(
      seq_len(n_perm),
      function(i) {
        rlang::inject(.lm_one_perm(!!!perm_args))
      }
    )
  }

  t_glob <- do.call(rbind, lapply(perm_results, `[[`, "t_glob_row"))
  # t_part: array dim c(n_perm, nvar+1, p)
  t_part <- array(dim = c(n_perm, nvar + 1L, p), data = NA_real_)
  for (i in seq_len(n_perm)) {
    t_part[i, , ] <- perm_results[[i]]$t_part_row
  }

  pval_glob <- colSums(
    t_glob >= matrix(t0_glob, nrow = n_perm, ncol = p, byrow = TRUE)
  ) /
    n_perm
  pval_part <- matrix(nrow = nvar + 1L, ncol = p)
  for (i in seq_len(p)) {
    pval_part[, i] <- colSums(
      t_part[,, i] >=
        matrix(t0_part[, i], nrow = n_perm, ncol = nvar + 1L, byrow = TRUE)
    ) /
      n_perm
  }

  list(
    coeff = coeff,
    n = n,
    p = p,
    nvar = nvar,
    var_names = var_names,
    design_matrix = design_matrix,
    regr0 = regr0,
    t0_part = t0_part,
    t0_glob = t0_glob,
    t_glob = t_glob,
    t_part = t_part,
    pval_glob = pval_glob,
    pval_part = pval_part,
    method = method
  )
}

optimize_order <- function(x) {
  n <- length(x)
  half <- ceiling(n / 2)
  x_half_1 <- x[1:half]
  x_half_2 <- rev(x[(half + 1):n])
  x_reordered <- numeric(n)
  x_reordered[seq(1, n, by = 2)] <- x_half_1
  x_reordered[seq(2, n, by = 2)] <- x_half_2
  x_reordered
}
