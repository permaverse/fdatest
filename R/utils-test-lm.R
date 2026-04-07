# Internal helper: compute the pointwise test statistics for LM.
# Returns list(t_glob, t_part) where t_glob is a numeric vector of length p
# (global F-statistic) and t_part is a ((nvar+1) x p) matrix (per-coefficient
# t^2-statistics).
lm_pointwise_stat <- function(
  coeff_fit,
  coeff_regr,
  n,
  p,
  nvar,
  sigma,
  resvar
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

  se <- sqrt(
    matrix(diag(sigma), nrow = nvar + 1L, ncol = p, byrow = FALSE) *
      matrix(resvar, nrow = nvar + 1L, ncol = p, byrow = TRUE)
  )
  t_part <- abs(coeff_regr / se)^2

  list(t_glob = t_glob, t_part = t_part)
}

# Internal helper: one permutation iteration for lm_permtest.
# Returns list(t_glob_row, t_part_row) — same shape as lm_pointwise_stat().
lm_single_perm <- function(
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

  stat <- lm_pointwise_stat(
    coeff_fit = regr_perm$fitted,
    coeff_regr = regr_perm$coeff,
    n = n,
    p = p,
    nvar = nvar,
    sigma = sigma_perm,
    resvar = resvar
  )

  if (method != "responses") {
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
      stat$t_part[ii, ] <- abs(regr_perm_ii$coeff / se_ii)[ii, ]^2
    }
  }

  list(t_glob_row = stat$t_glob, t_part_row = stat$t_part)
}

# Internal helper: shared pointwise permutation test for LM functions
# (global_lm, iwt_lm, twt_lm).
# Returns list(t0_glob, t0_part, t_glob, t_part, pval_glob, pval_part).
lm_permtest <- function(
  coeff,
  n,
  p,
  nvar,
  design_matrix,
  regr0,
  method,
  residui,
  fitted_part,
  n_perm
) {
  sigma <- chol2inv(regr0$qr$qr)
  resvar <- colSums(regr0$residuals^2) / regr0$df.residual

  stat0 <- lm_pointwise_stat(
    coeff_fit = regr0$fitted,
    coeff_regr = regr0$coeff,
    n = n,
    p = p,
    nvar = nvar,
    sigma = sigma,
    resvar = resvar
  )
  t0_glob <- stat0$t_glob
  t0_part <- stat0$t_part

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
  # Each task returns list(t_glob_row, t_part_row).
  if (mirai::daemons_set()) {
    perm_tasks <- mirai::mirai_map(seq_len(n_perm), function(.x) {
      rlang::inject(lm_single_perm(!!!perm_args))
    })
    perm_results <- perm_tasks[.progress]
  } else {
    perm_results <- lapply(seq_len(n_perm), function(.x) {
      rlang::inject(lm_single_perm(!!!perm_args))
    })
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
    t0_glob = t0_glob,
    t0_part = t0_part,
    t_glob = t_glob,
    t_part = t_part,
    pval_glob = pval_glob,
    pval_part = pval_part
  )
}

# Internal helper: data preparation + pointwise permutation test for
# LM functions (global_lm, iwt_lm, twt_lm). Parses formula, builds
# design matrix and model quantities, then calls lm_permtest(). Returns
# a list with all computed quantities needed by the p-value adjustment step.
lm_prepare_data <- function(formula, dx, n_perm, method) {
  coeff <- formula2coeff(formula, dx = dx)
  design_matrix <- formula2design_matrix(formula, coeff)

  nvar <- dim(design_matrix)[2] - 1
  var_names <- colnames(design_matrix)
  p <- dim(coeff)[2]
  n <- dim(coeff)[1]

  regr0 <- stats::lm.fit(design_matrix, coeff)

  if (nvar == 0) {
    method <- "responses"
  }

  residui <- fitted_part <- NULL
  if (method == "residuals") {
    formula_const <- deparse(formula[[3]], width.cutoff = 500L)
    var_names2 <- var_names
    coeffnames <- paste0("coeff[,", as.character(seq_len(p)), "]")

    design_matrix_names2 <- design_matrix
    if (length(grep("factor", formula_const, fixed = TRUE)) > 0) {
      index_factor <- grep("factor", var_names, fixed = TRUE)
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
        lapply(regr0_part[[ii]], function(.x) .x$residuals)
      )
      fitted_part[ii, , ] <- simplify2array(
        lapply(regr0_part[[ii]], function(.x) .x$fitted.values)
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
      lapply(regr0_part[[ii]], function(.x) .x$residuals)
    )
    fitted_part[ii, , ] <- simplify2array(
      lapply(regr0_part[[ii]], function(.x) .x$fitted.values)
    )
  }

  perm_out <- lm_permtest(
    coeff = coeff,
    n = n,
    p = p,
    nvar = nvar,
    design_matrix = design_matrix,
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
