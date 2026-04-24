#' Heatmap plot of the Interval Wise Testing Procedure results
#'
#' Plotting function creating a ggplot2 graphical output of the IWT: the
#' p-value heat-map, the adjusted p-value function, and the functional data,
#' assembled via \pkg{patchwork}.
#'
#' @param IWT_result Results of the IWT, as created by
#'   [`functional_one_sample_test()`], [`iwt1()`], [`functional_two_sample_test()`],
#'   [`iwt2()`], or the legacy functions [`IWT1()`] and [`IWTaov()`]. When using
#'   [`functional_two_sample_test()`] or [`iwt2()`], `correction` must be `"IWT"`.
#' @param alpha Threshold for the interval-wise error rate used for the
#'   hypothesis test. Regions where the adjusted p-value is below `alpha` are
#'   highlighted. The default is `alpha = 0.05`.
#' @param abscissa_range Range of the plot abscissa. The default is `c(0, 1)`.
#' @param nlevel Number of desired color levels for the p-value heatmap. The
#'   default is `nlevel = 20`.
#' @param plot_unadjusted Flag indicating if the unadjusted p-value function
#'   has to be overlaid (dashed line) on the adjusted p-value panel. The
#'   default is `FALSE`.
#'
#' @return An object of class `patchwork` containing the assembled ggplot2
#'   panels, returned invisibly. The plot is also printed as a side effect.
#'
#' @seealso See [`plot.fos()`], [`plot.fts()`], [`plot.flm()`] and
#'   [`plot.faov()`] for the plot method applied to the IWT results of one-
#'   and two-population tests, linear models, and ANOVA, respectively.
#'
#' @references
#' Pini, A., & Vantini, S. (2018). Interval-wise testing for functional data.
#' \emph{Journal of Nonparametric Statistics}, 29(2), 407-424.
#'
#' Pini, A., Vantini, S., Colosimo, B. M., & Grasso, M. (2018). Domain‐selective
#' functional analysis of variance for supervised statistical profile monitoring
#' of signal data. \emph{Journal of the Royal Statistical Society: Series C
#' (Applied Statistics)} 67(1), 55-81.
#'
#' Abramowicz, K., Hager, C. K., Pini, A., Schelin, L., Sjostedt de Luna, S., &
#' Vantini, S. (2018). Nonparametric inference for functional‐on‐scalar linear
#' models applied to knee kinematic hop data after injury of the anterior
#' cruciate ligament. \emph{Scandinavian Journal of Statistics} 45(4),
#' 1036-1061.
#'
#' @export
#' @importFrom rlang .data
#' @examples
#' # Performing the IWT for one population
#' IWT_result <- functional_one_sample_test(
#'   NASAtemp$paris, mu = 4, n_perm = 10L
#' )
#'
#' # Plotting the results of the IWT
#' IWTimage(IWT_result, abscissa_range = c(0, 12))
#'
#' # Selecting the significant components at 5% level
#' which(IWT_result$adjusted_pvalues < 0.05)
IWTimage <- function(
  IWT_result,
  alpha = 0.05,
  abscissa_range = c(0, 1),
  nlevel = 20L,
  plot_unadjusted = FALSE
) {
  out <- if (inherits(IWT_result, "fos")) {
    if (is.null(IWT_result$pvalue_matrix)) {
      cli::cli_abort(
        "IWTimage() requires an IWT result with a p-value matrix.
        Use {.arg correction = \"IWT\"} when calling {.fn functional_one_sample_test}."
      )
    }
    .iwt_image_one(
      pval_matrix = IWT_result$pvalue_matrix,
      adjusted_pval = IWT_result$adjusted_pvalues,
      unadjusted_pval = IWT_result$unadjusted_pvalues,
      data_eval = IWT_result$data,
      group_colors = rep(1L, nrow(IWT_result$data)),
      mu = IWT_result$mu,
      alpha = alpha,
      abscissa_range = abscissa_range,
      nlevel = nlevel,
      plot_unadjusted = plot_unadjusted
    )
  } else if (inherits(IWT_result, "fts")) {
    if (is.null(IWT_result$pvalue_matrix)) {
      cli::cli_abort(
        "IWTimage() requires an IWT result with a p-value matrix.
        Use {.arg correction = \"IWT\"} when calling {.fn functional_two_sample_test}."
      )
    }
    .iwt_image_one(
      pval_matrix = IWT_result$pvalue_matrix,
      adjusted_pval = IWT_result$adjusted_pvalues,
      unadjusted_pval = IWT_result$unadjusted_pvalues,
      data_eval = IWT_result$data,
      group_colors = IWT_result$group_labels,
      mu = NULL,
      alpha = alpha,
      abscissa_range = abscissa_range,
      nlevel = nlevel,
      plot_unadjusted = plot_unadjusted
    )
  } else if (inherits(IWT_result, "IWT1")) {
    .iwt_image_one(
      pval_matrix = IWT_result$pval_matrix,
      adjusted_pval = IWT_result$adjusted_pval,
      unadjusted_pval = IWT_result$unadjusted_pval,
      data_eval = IWT_result$data_eval,
      group_colors = rep(1L, nrow(IWT_result$data_eval)),
      mu = IWT_result$mu,
      alpha = alpha,
      abscissa_range = abscissa_range,
      nlevel = nlevel,
      plot_unadjusted = plot_unadjusted
    )
  } else if (inherits(IWT_result, "IWTaov")) {
    .iwt_image_aov(
      IWT_result = IWT_result,
      alpha = alpha,
      abscissa_range = abscissa_range,
      nlevel = nlevel,
      plot_unadjusted = plot_unadjusted
    )
  } else {
    cli::cli_abort(
      "Unsupported class {.cls {class(IWT_result)}}. \\
       {.fn IWTimage} accepts objects of class {.cls fos}, {.cls fts}, \\
       {.cls IWT1}, or {.cls IWTaov}."
    )
  }

  print(out)
  invisible(out)
}

# --------------------------------------------------------------------------
# Internal helpers
# --------------------------------------------------------------------------

# Convert a p×p pval_matrix into a long-format data frame for geom_raster(),
# reproducing the same triangular layout as the base-R image() call.
.iwt_heatmap_data <- function(pval_matrix, p, abscissa_range) {
  hm <- matrix(NA_real_, nrow = p, ncol = 4L * p)
  for (i in 0L:(p - 1L)) {
    for (j in seq_len(2L * p)) {
      hm[p - i, j + i + p] <- pval_matrix[p - i, (j + 1L) %/% 2L]
      if (j + i > 2L * p - i) {
        hm[p - i, j + i - p] <- pval_matrix[p - i, (j + 1L) %/% 2L]
      }
    }
  }
  matrice_quad <- hm[, (p + 1L):(3L * p), drop = FALSE]

  x_grid <- seq(abscissa_range[1], abscissa_range[2], length.out = 2L * p)
  y_grid <- seq(abscissa_range[1], abscissa_range[2], length.out = p) -
    abscissa_range[1]

  # image(x, y, t(matrice_quad[p:1,])) maps:
  #   z[xi, yi] = matrice_quad[p+1-yi, xi]  at (x_grid[xi], y_grid[yi])
  grid <- expand.grid(xi = seq_len(2L * p), yj = seq_len(p))
  grid$x <- x_grid[grid$xi]
  grid$y <- y_grid[grid$yj]
  grid$pval <- matrice_quad[cbind(p + 1L - grid$yj, grid$xi)]
  grid[!is.na(grid$pval), c("x", "y", "pval")]
}

# Build three ggplot2 panels (heatmap / adjusted p-values / functional data)
# for one test dimension (e.g., the global F-test or a single factor).
.iwt_image_panels <- function(
  pval_matrix,
  adjusted_pval,
  unadjusted_pval,
  data_eval,
  group_colors,
  abscissa_range,
  alpha,
  nlevel,
  plot_unadjusted,
  heat_title,
  mu = NULL
) {
  p <- length(adjusted_pval)
  abscissa_pts <- seq(abscissa_range[1], abscissa_range[2], length.out = p)
  step <- if (p > 1L) abscissa_pts[2L] - abscissa_pts[1L] else 0

  iwt_colors <- rev(grDevices::rainbow(nlevel, start = 0.15, end = 0.67))

  # --- Heatmap ---
  hm_data <- .iwt_heatmap_data(pval_matrix, p, abscissa_range)
  p_heat <- ggplot2::ggplot(
    hm_data,
    ggplot2::aes(x = .data$x, y = .data$y, fill = .data$pval)
  ) +
    ggplot2::geom_raster(interpolate = FALSE) +
    ggplot2::scale_fill_gradientn(
      colors = iwt_colors,
      limits = c(0, 1),
      name = "p-value",
      na.value = "white"
    ) +
    ggplot2::coord_equal(expand = FALSE) +
    ggplot2::labs(
      subtitle = heat_title,
      x = "Abscissa",
      y = "Interval length"
    ) +
    ggplot2::theme_minimal()

  # --- Shared significance band layer ---
  sig_idx <- which(adjusted_pval < alpha)
  sig_layer <- if (length(sig_idx) > 0L) {
    sig_df <- data.frame(
      xmin = abscissa_pts[sig_idx] - step / 2,
      xmax = abscissa_pts[sig_idx] + step / 2
    )
    ggplot2::geom_rect(
      data = sig_df,
      ggplot2::aes(
        xmin = .data$xmin,
        xmax = .data$xmax,
        ymin = -Inf,
        ymax = Inf
      ),
      inherit.aes = FALSE,
      fill = "gray90",
      alpha = 1
    )
  } else {
    NULL
  }

  # --- Adjusted p-value panel ---
  pval_df <- data.frame(
    x = abscissa_pts,
    adj = adjusted_pval,
    unadj = unadjusted_pval
  )
  p_pval <- ggplot2::ggplot(pval_df, ggplot2::aes(x = .data$x)) +
    sig_layer +
    ggplot2::geom_hline(
      yintercept = seq(0, 1, 0.1),
      color = "lightgray",
      linetype = "dotted"
    ) +
    ggplot2::geom_line(ggplot2::aes(y = .data$adj), linewidth = 1) +
    (if (plot_unadjusted) {
      ggplot2::geom_line(
        ggplot2::aes(y = .data$unadj),
        linewidth = 1,
        linetype = "dashed"
      )
    }) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(
      subtitle = "Adjusted p-value function",
      x = "Abscissa",
      y = "p-value"
    ) +
    ggplot2::theme_minimal()

  # --- Functional data panel ---
  data_long <- data.frame(
    x = rep(abscissa_pts, nrow(data_eval)),
    y = as.vector(t(data_eval)),
    id = as.factor(rep(seq_len(nrow(data_eval)), each = p)),
    grp = as.factor(rep(group_colors, each = p))
  )
  p_data <- ggplot2::ggplot(
    data_long,
    ggplot2::aes(
      x = .data$x,
      y = .data$y,
      group = .data$id,
      color = .data$grp
    )
  ) +
    sig_layer +
    ggplot2::geom_line(linewidth = 0.5) +
    ggplot2::scale_color_viridis_d(name = "Group") +
    ggplot2::labs(
      subtitle = "Functional data",
      x = "Abscissa",
      y = "Value"
    ) +
    ggplot2::theme_minimal()

  if (!is.null(mu)) {
    mu_vals <- if (length(mu) == 1L) rep(mu, p) else mu
    mu_df <- data.frame(x = abscissa_pts, y = mu_vals)
    p_data <- p_data +
      ggplot2::geom_line(
        data = mu_df,
        ggplot2::aes(x = .data$x, y = .data$y, group = NULL),
        color = "blue",
        linewidth = 1,
        linetype = "dashed",
        inherit.aes = FALSE
      )
  }

  list(heat = p_heat, pval = p_pval, data = p_data)
}

# Assemble a single-comparison IWT image (one-sample or two-sample).
.iwt_image_one <- function(
  pval_matrix,
  adjusted_pval,
  unadjusted_pval,
  data_eval,
  group_colors,
  mu,
  alpha,
  abscissa_range,
  nlevel,
  plot_unadjusted
) {
  panels <- .iwt_image_panels(
    pval_matrix = pval_matrix,
    adjusted_pval = adjusted_pval,
    unadjusted_pval = unadjusted_pval,
    data_eval = data_eval,
    group_colors = group_colors,
    abscissa_range = abscissa_range,
    alpha = alpha,
    nlevel = nlevel,
    plot_unadjusted = plot_unadjusted,
    heat_title = "p-value heatmap",
    mu = mu
  )
  patchwork::wrap_plots(panels$heat, panels$pval, panels$data, ncol = 1L)
}

# Assemble an IWTaov image: F-test column + one column per factor.
.iwt_image_aov <- function(
  IWT_result,
  alpha,
  abscissa_range,
  nlevel,
  plot_unadjusted
) {
  nvar <- dim(IWT_result$adjusted_pval_factor)[1L]
  factor_names <- rownames(IWT_result$unadjusted_pval_factor)
  all_names <- colnames(IWT_result$design_matrix)
  interaz <- grep(":", all_names, fixed = TRUE)

  # F-test column
  col_pop <- factor(apply(IWT_result$design_matrix, 1L, paste, collapse = ""))
  f_panels <- .iwt_image_panels(
    pval_matrix = IWT_result$pval_matrix_F,
    adjusted_pval = IWT_result$adjusted_pval_F,
    unadjusted_pval = IWT_result$unadjusted_pval_F,
    data_eval = IWT_result$data_eval,
    group_colors = col_pop,
    abscissa_range = abscissa_range,
    alpha = alpha,
    nlevel = nlevel,
    plot_unadjusted = plot_unadjusted,
    heat_title = "p-value heatmap, F test",
    mu = NULL
  )
  columns <- list(
    patchwork::wrap_plots(
      f_panels$heat,
      f_panels$pval,
      f_panels$data,
      ncol = 1L
    )
  )

  for (var in seq_len(nvar)) {
    var_name <- factor_names[var]
    if (length(grep(":", var_name, fixed = TRUE)) > 0L) {
      # Interaction term: colour by the intersection of the two factor columns
      var12 <- strsplit(var_name, ":", fixed = TRUE)[[1L]]
      dummy_test <- intersect(
        grep(var12[1L], all_names, fixed = TRUE),
        grep(var12[2L], all_names, fixed = TRUE)
      )
    } else {
      # Main effect: colour by matched design-matrix columns (excluding interactions)
      dummy_test <- setdiff(
        grep(var_name, all_names, fixed = TRUE),
        interaz
      )
    }
    colors <- IWT_result$design_matrix[, dummy_test]
    if (!is.null(dim(colors))) {
      colors <- apply(colors, 1L, paste, collapse = "")
    }
    colors <- as.factor(colors)

    var_panels <- .iwt_image_panels(
      pval_matrix = IWT_result$pval_matrix_factor[var, , ],
      adjusted_pval = IWT_result$adjusted_pval_factor[var, ],
      unadjusted_pval = IWT_result$unadjusted_pval_factor[var, ],
      data_eval = IWT_result$data_eval,
      group_colors = colors,
      abscissa_range = abscissa_range,
      alpha = alpha,
      nlevel = nlevel,
      plot_unadjusted = plot_unadjusted,
      heat_title = paste0("p-value heatmap, factor: ", var_name),
      mu = NULL
    )
    columns[[var + 1L]] <- patchwork::wrap_plots(
      var_panels$heat,
      var_panels$pval,
      var_panels$data,
      ncol = 1L
    )
  }

  patchwork::wrap_plots(columns, nrow = 1L)
}
