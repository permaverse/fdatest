#' Plot for Functional ANOVA Test Results
#'
#' The `S3` methods `autoplot.fanova()` and `plot.fanova()` are methods
#' for plotting results of functional analysis of variance tests. They visualize the
#' functional data and the adjusted p-values obtained from the testing
#' procedures for mean comparison of multiple groups. The plots highlight significant
#' effects at two levels of significance, `alpha1` and `alpha2`, using shaded
#' areas.
#'
#' @param object,x An object of class `fanova`, usually a result of a call
#'   to [`functional_anova_test()`], [`iwt_aov()`], [`twt_aov()`] or [`global_aov()`].
#' @param xrange A length-2 numeric vector specifying the range of the x-axis
#'   for the plots. Defaults to `c(0, 1)`. This should match the domain of the
#'   functional data.
#' @param alpha1 A numeric value specifying the first level of significance used
#'   to select and display significant effects. Defaults to `alpha1 = 0.05`.
#' @param alpha2 A numeric value specifying the second level of significance
#'   used to select and display significant effects. Defaults to `alpha2 =
#'   0.01`.
#' @param plot_adjpval A boolean value specifying whether the plots of adjusted
#'   p-values should be displayed. Defaults to `FALSE`.
#' @param ylim A 2-length numeric vector specifying the range of the y-axis.
#'   Defaults to `NULL`, which determines automatically the range from functional data.
#' @param col An integer specifying the color for the plot of functional data. Defaults
#'   to `1L`.
#' @param ylabel,ylab A string specifying the label of the y-axis of the functional
#'   data plot. Defaults to `"Functional Data"`.
#' @param title,main A string specifying the title of the functional data plot.
#'   Defaults to `NULL` in which case no title is displayed.
#' @param linewidth,lwd A numeric value specifying the width of the line for the
#'   functional data plot. Note that the line width for the adjusted p-value
#'   plot will be twice this value. Defaults to `0.5`.
#' @param type A string specifying the type of plot for the functional data. Defaults
#'   to `"l"` for lines.
#' @param ... Other arguments passed to specific methods. Not used in this
#'   function.
#'
#' @returns The `autoplot.fanova()` function creates a ggplot object that
#'   displays the functional data and the adjusted p-values. The significant
#'   intervals at levels `alpha1` and `alpha2` are highlighted in the plots.
#'   The `plot.fanova()` function is a wrapper around `autoplot.fanova()`
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
#' @name plot.fanova
#'
#' @examples
#' temperature <- rbind(NASAtemp$milan, NASAtemp$paris)
#' groups <- c(rep(0, 22), rep(1, 22))
#'
#' # Performing the TWT
#' TWT_result <- functional_anova_test(
#'   temperature ~ groups,
#'   correction = "TWT",
#'   B = 5L
#' )
#'
#' # Plotting the results of the TWT
#' plot(
#'   TWT_result,
#'   xrange = c(0, 12),
#'   main = 'TWT results for testing mean differences'
#' )
NULL

#' @rdname plot.fanova
#' @importFrom ggplot2 autoplot
#' @importFrom rlang .data
#' @export
autoplot.fanova <- function(
  object,
  xrange = c(0, 1),
  alpha1 = 0.05,
  alpha2 = 0.01,
  plot_adjpval = FALSE,
  ylim = NULL,
  col = 1,
  ylabel = "Functional Data",
  title = NULL,
  linewidth = 0.5,
  type = "l",
  ...
) {
  if (alpha1 < alpha2) {
    temp <- alpha1
    alpha1 <- alpha2
    alpha2 <- temp
  }

  nvar <- dim(object$adjusted_pval_factors)[1]
  p <- length(object$unadjusted_pval_F)
  xmin <- xrange[1]
  xmax <- xrange[2]
  abscissa_pval <- seq(xmin, xmax, len = p)
  dx <- abscissa_pval[2] - abscissa_pval[1]

  if (is.null(ylim)) {
    ylim <- range(object$data_eval)
  }

  # Helper: build a significance ribbon data frame
  make_sig_df <- function(pvals, abscissa, alpha1, alpha2) {
    sig1 <- which(pvals < alpha1)
    sig2 <- which(pvals < alpha2)
    dx <- abscissa[2] - abscissa[1]
    rows1 <- if (length(sig1) > 0) {
      data.frame(
        xmin = abscissa[sig1] - dx / 2,
        xmax = abscissa[sig1] + dx / 2,
        level = "alpha1"
      )
    } else {
      data.frame(xmin = numeric(0), xmax = numeric(0), level = character(0))
    }
    rows2 <- if (length(sig2) > 0) {
      data.frame(
        xmin = abscissa[sig2] - dx / 2,
        xmax = abscissa[sig2] + dx / 2,
        level = "alpha2"
      )
    } else {
      data.frame(xmin = numeric(0), xmax = numeric(0), level = character(0))
    }
    rbind(rows1, rows2)
  }

  # Helper: add significance ribbons to a ggplot
  add_sig_ribbons <- function(p, sig_df, ylim) {
    if (nrow(sig_df) == 0) {
      return(p)
    }
    p +
      ggplot2::geom_rect(
        data = sig_df,
        ggplot2::aes(
          xmin = .data$xmin,
          xmax = .data$xmax,
          ymin = ylim[1],
          ymax = ylim[2],
          fill = .data$level
        ),
        inherit.aes = FALSE,
        alpha = 1
      ) +
      ggplot2::scale_fill_manual(
        name = "Significance",
        values = c("alpha1" = "gray90", "alpha2" = "gray80"),
        labels = c(
          "alpha1" = paste0("p < ", alpha1),
          "alpha2" = paste0("p < ", alpha2)
        ),
        drop = FALSE
      )
  }

  # Build long-format data for functional data lines
  data_long <- as.data.frame(t(object$data_eval))
  data_long$x <- abscissa_pval
  data_long <- stats::reshape(
    data_long,
    varying = setdiff(names(data_long), "x"),
    v.names = "y",
    timevar = "curve",
    times = setdiff(names(data_long), "x"),
    direction = "long"
  )[, c("x", "curve", "y")]

  plots <- list()

  # --- F-test functional data plot ---
  if (nvar > 1) {
    sig_df_f <- make_sig_df(
      object$adjusted_pval_F,
      abscissa_pval,
      alpha1,
      alpha2
    )
    main_f <- if (!is.null(title)) {
      paste0(title, ": Functional Data and F-test")
    } else {
      "Functional Data and F-test"
    }

    p_f <- ggplot2::ggplot(
      data_long,
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$curve)
    )
    p_f <- add_sig_ribbons(p_f, sig_df_f, ylim)
    p_f <- p_f +
      ggplot2::geom_line(
        color = if (is.numeric(col)) "black" else col,
        linewidth = linewidth
      ) +
      ggplot2::coord_cartesian(ylim = ylim) +
      ggplot2::labs(title = main_f, x = NULL, y = ylabel) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "bottom")

    plots <- c(plots, list(p_f))
  }

  # --- Per-factor functional data plots ---
  names_all <- colnames(object$design_matrix)
  interaz <- grep(":", names_all)

  for (var in seq_len(nvar)) {
    var_name <- rownames(object$adjusted_pval_factors)[var]
    main_t <- if (!is.null(title)) {
      paste(title, ": factor", var_name)
    } else {
      paste("Factor", var_name)
    }

    if (length(grep(":", var_name)) > 0) {
      var12 <- strsplit(var_name, ":")
      var1 <- var12[[1]][1]
      var2 <- var12[[1]][2]
      dummy_test1 <- grep(var1, names_all)
      dummy_test2 <- grep(var2, names_all)
      dummy_test <- intersect(dummy_test1, dummy_test2)
      colors <- object$design_matrix[, dummy_test]
      if (length(dim(colors)) > 1) {
        colors <- apply(colors, 1, paste, collapse = "")
      }
      colors <- as.factor(colors)
    } else {
      dummy_test <- grep(var_name, names_all)
      dummy_test <- setdiff(dummy_test, interaz)
      colors <- object$design_matrix[, dummy_test]
      if (length(dim(colors)) > 1) {
        colors <- apply(colors, 1, paste, collapse = "")
      }
      colors <- as.factor(colors)
    }

    # Map colors to curves (rows of data_eval = observations)
    n_obs <- ncol(object$data_eval)
    curve_names <- colnames(as.data.frame(t(object$data_eval)))
    if (is.null(curve_names)) {
      curve_names <- paste0("V", seq_len(n_obs))
    }
    color_map <- data.frame(curve = curve_names, group = colors)
    data_long_var <- merge(data_long, color_map, by = "curve", sort = FALSE)

    sig_df_var <- make_sig_df(
      object$adjusted_pval_factors[var, ],
      abscissa_pval,
      alpha1,
      alpha2
    )

    p_var <- ggplot2::ggplot(
      data_long_var,
      ggplot2::aes(
        x = .data$x,
        y = .data$y,
        group = .data$curve,
        color = .data$group
      )
    )
    p_var <- add_sig_ribbons(p_var, sig_df_var, ylim)
    p_var <- p_var +
      ggplot2::geom_line(linewidth = linewidth) +
      ggplot2::geom_hline(
        yintercept = 0,
        linetype = "dashed",
        color = "black"
      ) +
      ggplot2::coord_cartesian(ylim = ylim) +
      ggplot2::labs(title = main_t, x = NULL, y = ylabel, color = var_name) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "bottom")

    plots <- c(plots, list(p_var))
  }

  # --- Adjusted p-value plots ---
  if (plot_adjpval) {
    # F-test p-value plot
    pval_df_f <- data.frame(x = abscissa_pval, pval = object$adjusted_pval_F)
    main_p <- if (!is.null(title)) {
      paste0(title, ": Adjusted p-values - F-test")
    } else {
      "Adjusted p-values - F-test"
    }
    sig_df_f_p <- make_sig_df(
      object$adjusted_pval_F,
      abscissa_pval,
      alpha1,
      alpha2
    )

    p_pval_f <- ggplot2::ggplot(
      pval_df_f,
      ggplot2::aes(x = .data$x, y = .data$pval)
    )
    p_pval_f <- add_sig_ribbons(p_pval_f, sig_df_f_p, c(0, 1))
    p_pval_f <- p_pval_f +
      ggplot2::geom_hline(
        yintercept = seq(0, 1, by = 0.1),
        color = "lightgray",
        linetype = "dotted"
      ) +
      ggplot2::geom_line(linewidth = linewidth * 2) +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::labs(title = main_p, x = NULL, y = "p-value") +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "bottom")

    plots <- c(plots, list(p_pval_f))

    # Per-factor p-value plots
    for (var in seq_len(nvar)) {
      var_name <- rownames(object$adjusted_pval_factors)[var]
      main_p <- if (!is.null(title)) {
        paste(title, ": Adjusted p-values - factor", var_name)
      } else {
        paste("Adjusted p-values - factor", var_name)
      }

      pval_df_var <- data.frame(
        x = abscissa_pval,
        pval = object$adjusted_pval_factors[var, ]
      )
      sig_df_var_p <- make_sig_df(
        object$adjusted_pval_factors[var, ],
        abscissa_pval,
        alpha1,
        alpha2
      )

      p_pval_var <- ggplot2::ggplot(
        pval_df_var,
        ggplot2::aes(x = .data$x, y = .data$pval)
      )
      p_pval_var <- add_sig_ribbons(p_pval_var, sig_df_var_p, c(0, 1))
      p_pval_var <- p_pval_var +
        ggplot2::geom_hline(
          yintercept = seq(0, 1, by = 0.1),
          color = "lightgray",
          linetype = "dotted"
        ) +
        ggplot2::geom_line(linewidth = linewidth * 2) +
        ggplot2::coord_cartesian(ylim = c(0, 1)) +
        ggplot2::labs(title = main_p, x = NULL, y = "p-value") +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "bottom")

      plots <- c(plots, list(p_pval_var))
    }
  }

  # Combine all plots with patchwork
  if (length(plots) == 1) {
    plots[[1]]
  } else {
    patchwork::wrap_plots(plots, ncol = 1)
  }
}

#' @rdname plot.fanova
#' @importFrom graphics plot
#' @export
plot.fanova <- function(
  x,
  xrange = c(0, 1),
  alpha1 = 0.05,
  alpha2 = 0.01,
  plot_adjpval = FALSE,
  ylim = NULL,
  col = 1,
  ylab = "Functional Data",
  main = NULL,
  lwd = 0.5,
  type = "l",
  ...
) {
  print(autoplot(
    x,
    xrange = xrange,
    alpha1 = alpha1,
    alpha2 = alpha2,
    plot_adjpval = plot_adjpval,
    ylim = ylim,
    col = col,
    ylabel = ylab,
    title = main,
    linewidth = lwd,
    type = type,
    ...
  ))
}
