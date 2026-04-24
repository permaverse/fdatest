#' Plot for Functional One-Sample Test Results
#'
#' The `S3` methods `autoplot.fos()` and `plot.fos()` are methods for plotting
#' results of functional one-sample tests. They visualize the functional data
#' and the adjusted p-values obtained from the testing procedures for testing
#' the center of symmetry of a functional population. The plots highlight
#' significant effects at two levels of significance, `alpha1` and `alpha2`,
#' using shaded areas.
#'
#' @param object,x An object of class `fos`, usually a result of a call to
#'   [`functional_one_sample_test()`] or [`iwt1()`].
#' @param xrange A length-2 numeric vector specifying the range of the x-axis
#'   for the plots. Defaults to `c(0, 1)`. This should match the domain of the
#'   functional data.
#' @param alpha1 A numeric value specifying the first level of significance used
#'   to select and display significant effects. Defaults to `alpha1 = 0.05`.
#' @param alpha2 A numeric value specifying the second level of significance
#'   used to select and display significant effects. Defaults to `alpha2 =
#'   0.01`.
#' @param ylabel A string specifying the label of the y-axis of the functional
#'   data plot. Defaults to `"Functional Data"`.
#' @param title A string specifying the title of the plots. Defaults to `NULL`
#'   in which case no title is displayed.
#' @param linewidth A numeric value specifying the width of the line for the
#'   functional data plot. Note that the line width for the adjusted p-value
#'   plot will be twice this value. Defaults to `linewidth = 0.5`.
#' @param ... Other arguments passed to specific methods. Not used in this
#'   function.
#'
#' @returns The `autoplot.fos()` function creates a ggplot object that displays
#'   the functional data (with the null mean function `mu` overlaid as a dashed
#'   reference line) and the adjusted p-values. The significant intervals at
#'   levels `alpha1` and `alpha2` are highlighted in both panels. The
#'   `plot.fos()` function is a wrapper around `autoplot.fos()` that prints the
#'   plot directly.
#'
#' @seealso [`IWTimage()`] for the plot of p-value heatmaps (for IWT).
#'
#' @references
#' Abramowicz, K., Pini, A., Schelin, L., Stamm, A., & Vantini, S. (2022).
#' “Domain selection and familywise error rate for functional data: A unified
#' framework. *Biometrics* 79(2), 1119-1132.
#'
#'
#' Pini, A., Vantini, S., Colosimo, B. M., & Grasso, M. (2018). Domain‐selective
#' functional analysis of variance for supervised statistical profile monitoring
#' of signal data. *Journal of the Royal Statistical Society: Series C
#' (Applied Statistics)* 67(1), 55-81.
#'
#' Abramowicz, K., Hager, C. K., Pini, A., Schelin, L., Sjostedt de Luna, S., &
#' Vantini, S. (2018). Nonparametric inference for functional‐on‐scalar linear
#' models applied to knee kinematic hop data after injury of the anterior
#' cruciate ligament. *Scandinavian Journal of Statistics* 45(4),
#' 1036-1061.
#' 
#' Pini, A., & Vantini, S. (2017). Interval-wise testing for functional data.
#' *Journal of Nonparametric Statistics*, 29(2), 407-424.
#'
#' @name plot.fos
#'
#' @examples
#' # Performing the IWT for one population
#' IWT_result <- functional_one_sample_test(
#'   NASAtemp$paris, mu = 4, n_perm = 10L
#' )
#'
#' # Plotting the results
#' plot(IWT_result, xrange = c(0, 12), title = "Paris temperatures")
#'
#' # Selecting the significant components at 5% level
#' which(IWT_result$adjusted_pvalues < 0.05)
NULL

#' @rdname plot.fos
#' @importFrom ggplot2 autoplot
#' @importFrom rlang .data
#' @export
autoplot.fos <- function(
  object,
  xrange = c(0, 1),
  alpha1 = 0.05,
  alpha2 = 0.01,
  ylabel = "Functional Data",
  title = NULL,
  linewidth = 0.5,
  ...
) {
  if (length(alpha1) != 1L || length(alpha2) != 1L) {
    cli::cli_abort(
      "{.arg alpha1} and {.arg alpha2} must each be a single numeric value."
    )
  }
  if (alpha1 < alpha2) {
    temp <- alpha1
    alpha1 <- alpha2
    alpha2 <- temp
  }

  p <- length(object$adjusted_pvalues)
  abscissa_pval <- seq(xrange[1], xrange[2], length.out = p)
  step <- if (p > 1L) abscissa_pval[2L] - abscissa_pval[1L] else 0

  # Significance regions
  sig_idx1 <- which(object$adjusted_pvalues < alpha1)
  sig_regions <- data.frame(
    xmin = if (length(sig_idx1) > 0) {
      abscissa_pval[sig_idx1] - step / 2
    } else {
      numeric(0)
    },
    xmax = if (length(sig_idx1) > 0) {
      abscissa_pval[sig_idx1] + step / 2
    } else {
      numeric(0)
    },
    alpha_level = if (length(sig_idx1) > 0) "alpha1" else character(0)
  )

  sig_idx2 <- which(object$adjusted_pvalues < alpha2)
  if (length(sig_idx2) > 0) {
    sig_regions <- rbind(
      sig_regions,
      data.frame(
        xmin = abscissa_pval[sig_idx2] - step / 2,
        xmax = abscissa_pval[sig_idx2] + step / 2,
        alpha_level = "alpha2"
      )
    )
  }

  sig_layer <- if (nrow(sig_regions) > 0) {
    ggplot2::geom_rect(
      data = sig_regions,
      ggplot2::aes(
        xmin = .data$xmin,
        xmax = .data$xmax,
        ymin = -Inf,
        ymax = Inf,
        fill = .data$alpha_level
      ),
      inherit.aes = FALSE,
      alpha = 0.3
    )
  } else {
    NULL
  }

  sig_layer_scale <- if (nrow(sig_regions) > 0) {
    ggplot2::scale_fill_manual(
      values = c("alpha1" = "gray70", "alpha2" = "gray50"),
      labels = c(
        "alpha1" = paste0("p < ", alpha1),
        "alpha2" = paste0("p < ", alpha2)
      ),
      name = "Significance"
    )
  } else {
    NULL
  }

  # Functional data in long format
  data_long <- data.frame(
    x = rep(abscissa_pval, nrow(object$data)),
    y = as.vector(t(object$data)),
    id = as.factor(rep(seq_len(nrow(object$data)), each = p))
  )

  # mu reference line
  mu_eval <- if (length(object$mu) == 1L) rep(object$mu, p) else object$mu
  mu_data <- data.frame(x = abscissa_pval, y = mu_eval)

  # Functional data panel
  p1 <- ggplot2::ggplot(
    data_long,
    ggplot2::aes(x = .data$x, y = .data$y, group = .data$id)
  ) +
    sig_layer +
    sig_layer_scale +
    ggplot2::geom_line(linewidth = linewidth, color = "gray30", alpha = 0.6) +
    ggplot2::geom_line(
      data = mu_data,
      ggplot2::aes(x = .data$x, y = .data$y, group = NULL),
      linewidth = linewidth * 2,
      color = "steelblue",
      linetype = "dashed"
    ) +
    ggplot2::labs(subtitle = "Functional Data", x = "Domain", y = ylabel) +
    ggplot2::theme_minimal()

  # Adjusted p-values panel
  pval_data <- data.frame(
    x = abscissa_pval,
    pval = object$adjusted_pvalues
  )

  p2 <- ggplot2::ggplot(
    pval_data,
    ggplot2::aes(x = .data$x, y = .data$pval)
  ) +
    sig_layer +
    sig_layer_scale +
    ggplot2::geom_hline(
      yintercept = seq(0, 1, 0.1),
      color = "lightgray",
      linetype = "dotted"
    ) +
    ggplot2::geom_line(linewidth = 2 * linewidth) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(subtitle = "Adjusted p-values", x = "Domain", y = "p-value") +
    ggplot2::theme_minimal()

  patchwork::wrap_plots(
    p1,
    p2,
    ncol = 1,
    guides = "collect",
    axis_titles = "collect"
  ) +
    patchwork::plot_annotation(
      title = title,
      theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    )
}

#' @rdname plot.fos
#' @importFrom graphics plot
#' @export
plot.fos <- function(
  x,
  xrange = c(0, 1),
  alpha1 = 0.05,
  alpha2 = 0.01,
  ylabel = "Functional Data",
  title = NULL,
  linewidth = 0.5,
  ...
) {
  print(autoplot(
    x,
    xrange = xrange,
    alpha1 = alpha1,
    alpha2 = alpha2,
    ylabel = ylabel,
    title = title,
    linewidth = linewidth,
    ...
  ))
}

# --------------------------------------------------------------------------
# Backward-compatible plot.IWT1 — delegates to autoplot.fos().
# The legacy parameter names (ylab, main, lwd) are preserved so that
# existing user code continues to work unchanged.
# --------------------------------------------------------------------------
#' @rdname plot.fos
#' @param ylab Label of the y-axis (legacy alias for `ylabel`). Defaults to
#'   `"Functional Data"`.
#' @param main Plot title (legacy alias for `title`). Defaults to `NULL`.
#' @param lwd Line width (legacy alias for `linewidth`; divided by 2 for
#'   ggplot2 scaling). Defaults to `1`.
#' @param col,ylim,type Ignored; retained for backward compatibility only.
#' @export
plot.IWT1 <- function(
  x,
  xrange = c(0, 1),
  alpha1 = 0.05,
  alpha2 = 0.01,
  ylab = "Functional Data",
  main = NULL,
  lwd = 1,
  col = 1,
  ylim = NULL,
  type = "l",
  ...
) {
  fos_obj <- structure(
    list(
      data = x$data_eval,
      mu = x$mu,
      unadjusted_pvalues = x$unadjusted_pval,
      adjusted_pvalues = x$adjusted_pval,
      correction_method = "IWT",
      pvalue_matrix = x$pval_matrix
    ),
    class = "fos"
  )
  print(autoplot.fos(
    fos_obj,
    xrange = xrange,
    alpha1 = alpha1,
    alpha2 = alpha2,
    ylabel = ylab,
    title = main,
    linewidth = lwd / 2
  ))
}
