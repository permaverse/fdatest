#' Plot for Functional Regression Test Results
#'
#' The `S3` methods `autoplot.fanova()` and `plot.fanova()` are methods
#' for plotting results of functional analysis of variance tests. They visualize the
#' functional data and the adjusted p-values obtained from the testing
#' procedures for mean comparison of multiple groups. The plots highlight significant
#' effects at two levels of significance, `alpha1` and `alpha2`, using shaded
#' areas.
#'
#' @param object,x  The object to be plotted. An object of class "\code{IWTlm}",
#'   usually, a result of a call to \code{\link{IWTlm}}.
#' @inherit plot.fanova params seealso
#'
#' @returns The `autoplot.flm()` function creates a ggplot object that
#'   displays the functional data and the adjusted p-values. The significant
#'   intervals at levels `alpha1` and `alpha2` are highlighted in the plots.
#'   The `plot.flm()` function is a wrapper around `autoplot.flm()`
#'   that prints the plot directly.
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
#' @name plot.flm
#'
#' @examples
#' temperature <- rbind(NASAtemp$milan, NASAtemp$paris)
#' groups <- c(rep(0, 22), rep(1, 22))
#'
#' # Performing the IWT
#' IWT.result <- IWTlm(temperature ~ groups, B = 2L)
#'
#' # Summary of the IWT results
#' summary(IWT.result)
#'
#' # Plot of the IWT results
#' plot(IWT.result)
#'
#' plot(
#'   IWT.result,
#'   main = 'NASA data',
#'   plot_adjpval = TRUE,
#'   xlab = 'Day',
#'   xrange = c(1, 365)
#' )
NULL

#' @rdname plot.flm
#' @importFrom ggplot2 autoplot
#' @importFrom rlang .data
#' @export
autoplot.flm <- function(
  object,
  xrange = c(0, 1),
  alpha1 = 0.05,
  alpha2 = 0.01,
  plot_adjpval = FALSE,
  col = c(1, grDevices::rainbow(dim(object$adjusted_pval_part)[1])),
  ylim = NULL,
  ylabel = "Functional Data",
  title = NULL,
  linewidth = 1,
  type = "l",
  ...
) {
  if (alpha1 < alpha2) {
    temp <- alpha1
    alpha1 <- alpha2
    alpha2 <- temp
  }

  p <- length(object$unadjusted_pval_F)
  J <- p
  n <- dim(object$data.eval)[1]
  xmin <- xrange[1]
  xmax <- xrange[2]
  abscissa_pval <- seq(xmin, xmax, len = p)
  abscissa_smooth <- seq(xmin, xmax, len = J)
  half_step <- (abscissa_pval[2] - abscissa_pval[1]) / 2

  # Helper: build significance rectangles data frame
  make_sig_df <- function(pvals, alpha1, alpha2, abscissa) {
    hs <- (abscissa[2] - abscissa[1]) / 2
    rows <- list()
    idx1 <- which(pvals < alpha1 & pvals >= alpha2)
    if (length(idx1) > 0) {
      rows[[length(rows) + 1]] <- data.frame(
        xmin = abscissa[idx1] - hs,
        xmax = abscissa[idx1] + hs,
        level = paste0("< ", alpha1)
      )
    }
    idx2 <- which(pvals < alpha2)
    if (length(idx2) > 0) {
      rows[[length(rows) + 1]] <- data.frame(
        xmin = abscissa[idx2] - hs,
        xmax = abscissa[idx2] + hs,
        level = paste0("< ", alpha2)
      )
    }
    if (length(rows) > 0) {
      df <- do.call(rbind, rows)
      df$level <- factor(
        df$level,
        levels = c(paste0("< ", alpha1), paste0("< ", alpha2))
      )
      df
    } else {
      NULL
    }
  }

  plots <- list()

  # --- Panel 1: Functional Data and F-test ---
  main_F <- if (!is.null(title)) {
    paste(title, ": Functional Data and F-test")
  } else {
    "Functional Data and F-test"
  }

  data_long <- data.frame(
    x = rep(abscissa_smooth, each = n),
    y = as.vector(t(object$data.eval)),
    id = rep(seq_len(n), times = J)
  )

  p_F <- ggplot2::ggplot()

  sig_df_F <- make_sig_df(object$adjusted_pval_F, alpha1, alpha2, abscissa_pval)
  if (!is.null(sig_df_F)) {
    p_F <- p_F +
      ggplot2::geom_rect(
        data = sig_df_F,
        ggplot2::aes(
          xmin = .data$xmin,
          xmax = .data$xmax,
          ymin = -Inf,
          ymax = Inf,
          fill = .data$level
        ),
        alpha = 1
      ) +
      ggplot2::scale_fill_manual(
        name = "Significance",
        values = stats::setNames(
          c("gray90", "gray80"),
          c(paste0("< ", alpha1), paste0("< ", alpha2))
        ),
        drop = FALSE
      )
  }

  p_F <- p_F +
    ggplot2::geom_line(
      data = data_long,
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$id),
      color = col[1],
      linewidth = linewidth
    ) +
    ggplot2::labs(title = main_F, x = NULL, y = ylabel) +
    ggplot2::theme_bw()

  if (!is.null(ylim)) {
    p_F <- p_F + ggplot2::coord_cartesian(ylim = ylim)
  }

  plots[[length(plots) + 1]] <- p_F

  # --- Panels for each variable: regression coefficient + t-test ---
  n_vars <- dim(object$adjusted_pval_part)[1]
  for (var in seq_len(n_vars)) {
    var_name <- rownames(object$adjusted_pval_part)[var]
    main_t <- if (!is.null(title)) {
      paste(title, ": t-test -", var_name)
    } else {
      paste("t-test -", var_name)
    }

    coeff_df <- data.frame(
      x = abscissa_smooth,
      y = object$coeff.regr.eval[var, ]
    )

    p_t <- ggplot2::ggplot()

    sig_df_t <- make_sig_df(
      object$adjusted_pval_part[var, ],
      alpha1,
      alpha2,
      abscissa_pval
    )
    if (!is.null(sig_df_t)) {
      p_t <- p_t +
        ggplot2::geom_rect(
          data = sig_df_t,
          ggplot2::aes(
            xmin = .data$xmin,
            xmax = .data$xmax,
            ymin = -Inf,
            ymax = Inf,
            fill = .data$level
          ),
          alpha = 1
        ) +
        ggplot2::scale_fill_manual(
          name = "Significance",
          values = stats::setNames(
            c("gray90", "gray80"),
            c(paste0("< ", alpha1), paste0("< ", alpha2))
          ),
          drop = FALSE
        )
    }

    p_t <- p_t +
      ggplot2::geom_hline(
        yintercept = 0,
        linetype = "dashed",
        color = "black"
      ) +
      ggplot2::geom_line(
        data = coeff_df,
        ggplot2::aes(x = .data$x, y = .data$y),
        color = col[var + 1],
        linewidth = linewidth
      ) +
      ggplot2::labs(title = main_t, x = NULL, y = "Regression Coefficient") +
      ggplot2::theme_bw()

    plots[[length(plots) + 1]] <- p_t
  }

  # --- Adjusted p-value plots ---
  if (plot_adjpval) {
    # F-test p-values
    main_pF <- if (!is.null(title)) {
      paste(title, ": Adjusted p-values - F-test")
    } else {
      "Adjusted p-values - F-test"
    }

    pval_F_df <- data.frame(x = abscissa_pval, y = object$adjusted_pval_F)

    p_pF <- ggplot2::ggplot()

    sig_df_pF <- make_sig_df(
      object$adjusted_pval_F,
      alpha1,
      alpha2,
      abscissa_pval
    )
    if (!is.null(sig_df_pF)) {
      p_pF <- p_pF +
        ggplot2::geom_rect(
          data = sig_df_pF,
          ggplot2::aes(
            xmin = .data$xmin,
            xmax = .data$xmax,
            ymin = -Inf,
            ymax = Inf,
            fill = .data$level
          ),
          alpha = 1
        ) +
        ggplot2::scale_fill_manual(
          name = "Significance",
          values = stats::setNames(
            c("gray90", "gray80"),
            c(paste0("< ", alpha1), paste0("< ", alpha2))
          ),
          drop = FALSE
        )
    }

    p_pF <- p_pF +
      ggplot2::geom_hline(
        yintercept = seq(0, 1, by = 0.1),
        color = "lightgray",
        linetype = "dotted"
      ) +
      ggplot2::geom_line(
        data = pval_F_df,
        ggplot2::aes(x = .data$x, y = .data$y),
        linewidth = linewidth
      ) +
      ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::labs(title = main_pF, x = NULL, y = "p-value") +
      ggplot2::theme_bw()

    plots[[length(plots) + 1]] <- p_pF

    # t-test p-values per variable
    for (var in seq_len(n_vars)) {
      var_name <- rownames(object$adjusted_pval_part)[var]
      main_pt <- if (!is.null(title)) {
        paste(title, ": Adjusted p-values - t-test -", var_name)
      } else {
        paste("Adjusted p-values - t-test -", var_name)
      }

      pval_t_df <- data.frame(
        x = abscissa_pval,
        y = object$adjusted_pval_part[var, ]
      )

      p_pt <- ggplot2::ggplot()

      sig_df_pt <- make_sig_df(
        object$adjusted_pval_part[var, ],
        alpha1,
        alpha2,
        abscissa_pval
      )
      if (!is.null(sig_df_pt)) {
        p_pt <- p_pt +
          ggplot2::geom_rect(
            data = sig_df_pt,
            ggplot2::aes(
              xmin = .data$xmin,
              xmax = .data$xmax,
              ymin = -Inf,
              ymax = Inf,
              fill = .data$level
            ),
            alpha = 1
          ) +
          ggplot2::scale_fill_manual(
            name = "Significance",
            values = stats::setNames(
              c("gray90", "gray80"),
              c(paste0("< ", alpha1), paste0("< ", alpha2))
            ),
            drop = FALSE
          )
      }

      p_pt <- p_pt +
        ggplot2::geom_hline(
          yintercept = seq(0, 1, by = 0.1),
          color = "lightgray",
          linetype = "dotted"
        ) +
        ggplot2::geom_line(
          data = pval_t_df,
          ggplot2::aes(x = .data$x, y = .data$y),
          linewidth = linewidth
        ) +
        ggplot2::scale_y_continuous(limits = c(0, 1)) +
        ggplot2::labs(title = main_pt, x = NULL, y = "p-value") +
        ggplot2::theme_bw()

      plots[[length(plots) + 1]] <- p_pt
    }
  }

  patchwork::wrap_plots(plots, ncol = 1)
}

#' @rdname plot.flm
#' @importFrom graphics plot
#' @export
plot.flm <- function(
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
