# Internal helper: compute the pointwise test statistics for AOV.
# Returns list(t_glob, t_part) where t_glob is a numeric vector of length p
# (global F-statistic) and t_part is a (nvar x p) matrix (per-factor
# F-statistics).
aov_pointwise_stat <- function(
  coeff_fit,
  n,
  p,
  nvar,
  index_vars,
  df_vars,
  resvar,
  effects,
  method_is_responses
) {
  t_glob <- if (nvar > 0) {
    colSums(
      (coeff_fit -
        matrix(colMeans(coeff_fit), nrow = n, ncol = p, byrow = TRUE))^2
    ) /
      (nvar * resvar)
  } else {
    numeric(p)
  }

  if (method_is_responses) {
    ms <- matrix(nrow = nvar + 1L, ncol = p)
    for (var in seq_len(nvar + 1L)) {
      ms[var, ] <- colSums(rbind(
        effects[index_vars[var, 1]:index_vars[var, 2], ]^2
      )) /
        df_vars[var]
    }
    t_part <- ms[seq_len(nvar), , drop = FALSE] /
      matrix(ms[nvar + 1L, ], nrow = nvar, ncol = p, byrow = TRUE)
  } else {
    t_part <- matrix(nrow = nvar, ncol = p)
  }

  list(t_glob = t_glob, t_part = t_part)
}

# Internal helper: one permutation iteration for aov_permtest.
# Returns list(t_glob, t_part) — same shape as aov_pointwise_stat().
aov_single_perm <- function(
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

  stat <- aov_pointwise_stat(
    coeff_fit = regr_perm$fitted,
    n = n,
    p = p,
    nvar = nvar,
    index_vars = index_vars,
    df_vars = df_vars,
    resvar = resvar,
    effects = regr_perm$effects,
    method_is_responses = (method == "responses")
  )

  if (method != "responses") {
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
      stat$t_part[ii, ] <- ratio_mat[ii, ]
    }
  }

  list(t_glob_row = stat$t_glob, t_part_row = stat$t_part)
}

# Internal helper: shared pointwise permutation test for AOV functions
# (global_aov, iwt_aov, twt_aov).
# Returns list(t0_glob, t0_part, t_glob, t_part, pval_glob, pval_part).
aov_permtest <- function(
  coeff,
  n,
  p,
  nvar,
  design_matrix,
  index_vars,
  df_vars,
  regr0,
  method,
  residui,
  fitted_part,
  n_perm
) {
  resvar <- colSums(regr0$residuals^2) / regr0$df.residual

  stat0 <- aov_pointwise_stat(
    coeff_fit = regr0$fitted,
    n = n,
    p = p,
    nvar = nvar,
    index_vars = index_vars,
    df_vars = df_vars,
    resvar = resvar,
    effects = regr0$effects,
    method_is_responses = (method == "responses")
  )
  t0_glob <- stat0$t_glob

  ms0 <- matrix(nrow = nvar + 1L, ncol = p)
  for (var in seq_len(nvar + 1L)) {
    ms0[var, ] <- colSums(rbind(
      regr0$effects[index_vars[var, 1]:index_vars[var, 2], ]^2
    )) /
      df_vars[var]
  }
  t0_part <- ms0[seq_len(nvar), , drop = FALSE] /
    matrix(ms0[nvar + 1L, ], nrow = nvar, ncol = p, byrow = TRUE)

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
  # Each task returns list(t_glob_row, t_part_row).
  if (mirai::daemons_set()) {
    perm_tasks <- mirai::mirai_map(seq_len(n_perm), function(.x) {
      rlang::inject(aov_single_perm(!!!perm_args))
    })
    perm_results <- perm_tasks[.progress]
  } else {
    perm_results <- lapply(seq_len(n_perm), function(.x) {
      rlang::inject(aov_single_perm(!!!perm_args))
    })
  }

  t_glob <- do.call(rbind, lapply(perm_results, `[[`, "t_glob_row"))
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
    t0_glob = t0_glob,
    t0_part = t0_part,
    t_glob = t_glob,
    t_part = t_part,
    pval_glob = pval_glob,
    pval_part = pval_part
  )
}

# Internal helper: data preparation + pointwise permutation test for
# AOV functions (global_aov, iwt_aov, twt_aov). Parses formula, builds
# design matrix and model quantities, then calls aov_permtest(). Returns
# a list with all computed quantities needed by the p-value adjustment step.
aov_prepare_data <- function(formula, dx, n_perm, method) {
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
    var_names[ii] <- gsub(" ", "", var_names[ii], fixed = TRUE)
  }

  index_vars <- cbind(
    c(2, (cumsum(df_vars) + 2)[-length(df_vars)]),
    cumsum(df_vars) + 1
  )
  regr0 <- stats::lm.fit(design_matrix, coeff)

  if (nvar > 1) {
    # keep requested method
  } else if (nvar == 1) {
    # only one factor: permuting residuals is equivalent to permuting responses
    method <- "responses"
  } else {
    # intercept only
    method <- "responses"
  }

  residui <- fitted_part <- NULL
  if (method == "residuals") {
    var_names2 <- var_names
    if (length(grep("factor", formula_const, fixed = TRUE)) > 0) {
      index_factor <- grep("factor", var_names, fixed = TRUE)
      replace_names <- paste0("group", seq_along(index_factor))
      var_names2[index_factor] <- replace_names
    }

    residui <- array(dim = c(nvar, n, p))
    fitted_part <- array(dim = c(nvar, n, p))
    formula_coeff_part <- vector("list", nvar)
    regr0_part <- vector("list", nvar)
    dummy_interaz <- grep(":", dummynames_all, fixed = TRUE)

    for (ii in seq_len(nvar)) {
      var_ii <- var_names2[ii]
      if (length(grep(":", var_ii, fixed = TRUE)) > 0) {
        var12 <- strsplit(var_ii, ":")
        var1 <- var12[[1]][1]
        var2 <- var12[[1]][2]
        dummy_test1 <- grep(var1, dummynames_all, fixed = TRUE)
        dummy_test2 <- grep(var2, dummynames_all, fixed = TRUE)
        dummy_test <- intersect(dummy_test1, dummy_test2)
        dummynames_reduced <- dummynames_all[-dummy_test]
      } else {
        dummy_test <- grep(var_ii, dummynames_all, fixed = TRUE)
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
        lapply(regr0_part[[ii]], function(.x) .x$residuals)
      )
      fitted_part[ii, , ] <- simplify2array(
        lapply(regr0_part[[ii]], function(.x) .x$fitted.values)
      )
    }
  }

  perm_out <- aov_permtest(
    coeff = coeff,
    n = n,
    p = p,
    nvar = nvar,
    design_matrix = design_matrix,
    index_vars = index_vars,
    df_vars = df_vars,
    regr0 = regr0,
    method = method,
    residui = residui,
    fitted_part = fitted_part,
    n_perm = n_perm
  )

  list(
    coeff = coeff,
    n = n,
    p = p,
    nvar = nvar,
    var_names = var_names,
    coeffnames = coeffnames,
    design_matrix = design_matrix,
    regr0 = regr0,
    t0_part = perm_out$t0_part,
    t0_glob = perm_out$t0_glob,
    t_glob = perm_out$t_glob,
    t_part = perm_out$t_part,
    pval_glob = perm_out$pval_glob,
    pval_part = perm_out$pval_part,
    method = method
  )
}
