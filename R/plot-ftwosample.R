#' Plot for Functional Two-Sample Test Results
#'
#' The `S3` methods `autoplot.ftwosample()` and `plot.ftwosample()` are methods
#' for plotting results of functional two-sample tests. They visualize the
#' functional data and the adjusted p-values obtained from the testing
#' procedures for mean comparison of two groups. The plots highlight significant
#' effects at two levels of significance, `alpha1` and `alpha2`, using shaded
#' areas.
#'
#' @param object,x An object of class `ftwosample`, usually a result of a call
#'   to [`functional_two_sample_test()`], [`IWT2()`], [`TWT2()`], [`FDR2()`],
#'   [`PCT2()`] or [`Global2()`].
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
#' @param title A string specifying the title of the functional data plot.
#'   Defaults to `NULL` in which case no title is displayed.
#' @param linewidth A numeric value specifying the width of the line for the
#'   functional data plot. Note that the line width for the adjusted p-value
#'   plot will be twice this value. Defaults to `linewidth = 0.5`.
#' @param ... Other arguments passed to specific methods. Not used in this
#'   function.
#'
#' @returns The `autoplot.ftwosample()` function creates a ggplot object that
#'   displays the functional data and the adjusted p-values. The significant
#'   intervals at levels `alpha1` and `alpha2` are highlighted in the plots. The
#'   `plot.ftwosample()` function is a wrapper around `autoplot.ftwosample()`
#'   that prints the plot directly.
#'
#' @seealso [`IWTimage()`] for the plot of p-values heatmaps (for IWT).
#'
#' @references
#' Pini, A., & Vantini, S. (2017). Interval-wise testing for functional data.
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
#' @name plot.ftwosample
#'
#' @examples
#' # Performing the TWT for two populations
#' TWT_result <- functional_two_sample_test(
#'   NASAtemp$paris, NASAtemp$milan,
#'   correction = "TWT", B = 10L
#' )
#'
#' # Plotting the results of the TWT
#' plot(
#'   TWT_result,
#'   xrange = c(0, 12),
#'   title = 'TWT results for testing mean differences'
#' )
#'
#' # Selecting the significant components at 5% level
#' which(TWT_result$adjusted_pval < 0.05)
#'
#' # Performing the IWT for two populations
#' IWT_result <- functional_two_sample_test(
#'   NASAtemp$paris, NASAtemp$milan,
#'   correction = "IWT", B = 10L
#' )
#'
#' # Plotting the results of the IWT
#' plot(
#'   IWT_result,
#'   xrange = c(0, 12),
#'   title = 'IWT results for testing mean differences'
#' )
#'
#' # Selecting the significant components at 5% level
#' which(IWT_result$adjusted_pval < 0.05)
NULL

#' @rdname plot.ftwosample
#' @importFrom ggplot2 autoplot
#' @importFrom rlang .data
#' @export
autoplot.ftwosample <- function(
  object,
  xrange = c(0, 1),
  alpha1 = 0.05,
  alpha2 = 0.01,
  ylabel = "Functional Data",
  title = NULL,
  linewidth = 0.5,
  ...
) {
  if (alpha1 < alpha2) {
    temp <- alpha1
    alpha1 <- alpha2
    alpha2 <- temp
  }
  abscissa_pval <- seq(
    xrange[1],
    xrange[2],
    length.out = length(object$adjusted_pvalues)
  )

  # Create data frame for functional data plot
  data_long <- data.frame(
    x = rep(
      seq(xrange[1], xrange[2], length.out = ncol(object$data)),
      nrow(object$data)
    ),
    y = as.vector(t(object$data)),
    group = as.factor(rep(object$group_labels, each = ncol(object$data))),
    id = as.factor(rep(seq_len(nrow(object$data)), each = ncol(object$data)))
  )

  # Add significance regions
  sig_idx1 <- which(object$adjusted_pvalues < alpha1)
  sig_regions <- data.frame(
    xmin = if (length(sig_idx1) > 0) {
      abscissa_pval[sig_idx1] - (abscissa_pval[2] - abscissa_pval[1]) / 2
    } else {
      numeric(0)
    },
    xmax = if (length(sig_idx1) > 0) {
      abscissa_pval[sig_idx1] + (abscissa_pval[2] - abscissa_pval[1]) / 2
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
        xmin = abscissa_pval[sig_idx2] -
          (abscissa_pval[2] - abscissa_pval[1]) / 2,
        xmax = abscissa_pval[sig_idx2] +
          (abscissa_pval[2] - abscissa_pval[1]) / 2,
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

  # Functional data plot
  p1 <- ggplot2::ggplot(
    data_long,
    ggplot2::aes(
      x = .data$x,
      y = .data$y,
      group = .data$id,
      color = .data$group
    )
  ) +
    sig_layer +
    ggplot2::geom_line(linewidth = linewidth) +
    ggplot2::scale_color_viridis_d(name = "Group") +
    ggplot2::scale_fill_manual(
      values = c("alpha1" = "gray70", "alpha2" = "gray50"),
      labels = c(
        "alpha1" = paste0("p < ", alpha1),
        "alpha2" = paste0("p < ", alpha2)
      ),
      name = "Significance"
    ) +
    ggplot2::labs(subtitle = "Functional Data", x = "Domain", y = ylabel) +
    ggplot2::theme_minimal()

  # P-values plot
  pval_data <- data.frame(
    x = seq(xrange[1], xrange[2], length.out = length(object$adjusted_pvalues)),
    pval = object$adjusted_pvalues
  )

  p2 <- ggplot2::ggplot(pval_data, ggplot2::aes(x = .data$x, y = .data$pval)) +
    sig_layer +
    ggplot2::geom_hline(
      yintercept = seq(0, 1, 0.1),
      color = "lightgray",
      linetype = "dotted"
    ) +
    ggplot2::geom_line(linewidth = 2 * linewidth) +
    ggplot2::scale_fill_manual(
      values = c("alpha1" = "gray70", "alpha2" = "gray50"),
      labels = c(
        "alpha1" = paste0("p < ", alpha1),
        "alpha2" = paste0("p < ", alpha2)
      ),
      name = "Significance"
    ) +
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

#' @rdname plot.ftwosample
#' @importFrom graphics plot
#' @export
plot.ftwosample <- function(
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
