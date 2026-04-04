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

# Use index mapping (modulo) instead of building a 2x duplicated matrix.
pval_correct <- function(pval_matrix) {
  p <- ncol(pval_matrix)
  corrected <- matrix(NA_real_, nrow = p, ncol = p)

  corrected[p, ] <- pval_matrix[p, p:1]

  for (var in seq_len(p)) {
    k0 <- var
    orig0 <- p - ((k0 - 1L) %% p)
    pval_var <- pval_matrix[p, orig0]
    fine <- var

    for (riga in seq.int(p - 1L, 1L)) {
      fine <- fine + 1L
      k_seq <- seq.int(var, fine)
      orig_idx <- p - ((k_seq - 1L) %% p)
      seg_max <- max(pval_matrix[riga, orig_idx], na.rm = TRUE)
      pval_var <- max(pval_var, seg_max, na.rm = TRUE)
      corrected[riga, var] <- pval_var
    }
  }

  corrected[, p:1]
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
