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
    perm_tasks <- mirai::mirai_map(seq_len(n_perm), function(i) {
      rlang::inject(.lm_one_perm(!!!perm_args))
    })
    perm_results <- perm_tasks[.progress]
  } else {
    perm_results <- lapply(seq_len(n_perm), function(i) {
      rlang::inject(.lm_one_perm(!!!perm_args))
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
