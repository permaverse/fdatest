# Tests for plot and autoplot methods:
# plot-ftwosample.R: autoplot.ftwosample, plot.ftwosample
# plot-fanova.R:     autoplot.fanova, plot.fanova
# plot-flm.R:        autoplot.flm, plot.flm
# plot.IWT1.R:       plot.IWT1
library(fdatest)
library(ggplot2)
library(tinysnapshot)
using(tinysnapshot)

# Use PNG device (no rsvg dependency) for snapshot comparisons
options(tinysnapshot_device = "png")
# Snapshots were produced on macOS; skip comparisons on other OSes
options(tinysnapshot_os = "Darwin")

data("NASAtemp", package = "fdatest")
d1 <- NASAtemp$milan[1:4, 1:8]
d2 <- NASAtemp$paris[1:4, 1:8]
temperature <- rbind(d1, d2)
groups <- c(rep(0L, 4L), rep(1L, 4L))

# Compute once to reuse across multiple plot tests
set.seed(42)
res_iwt2 <- IWT2(d1, d2, B = 5L, verbose = FALSE)
set.seed(42)
res_iwt1 <- IWT1(d1, B = 5L)
set.seed(42)
res_fanova <- IWTaov(temperature ~ groups, B = 5L)
set.seed(42)
res_flm <- IWTlm(temperature ~ groups, B = 5L)

# Also need TWT results (no pval_matrix) to exercise alternate branches
set.seed(42)
res_twt2 <- TWT2(d1, d2, B = 5L, verbose = FALSE)
set.seed(42)
res_twtaov <- TWTaov(temperature ~ groups, B = 5L)
set.seed(42)
res_twtlm <- TWTlm(temperature ~ groups, B = 5L)

# Two-factor design: exercises if(nvar > 1) F-test panel and per-factor loop
grpA <- c(0L, 0L, 1L, 1L, 0L, 0L, 1L, 1L)
grpB <- c(0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L)
set.seed(42)
res_fanova2 <- IWTaov(temperature ~ grpA + grpB, B = 5L)
set.seed(42)
res_twtaov2 <- TWTaov(temperature ~ grpA + grpB, B = 5L)
set.seed(42)
res_flm2 <- IWTlm(temperature ~ grpA + grpB, B = 5L)
set.seed(42)
res_twtlm2 <- TWTlm(temperature ~ grpA + grpB, B = 5L)

# ===========================================================================
# autoplot.ftwosample, plot.ftwosample
# ===========================================================================
# autoplot returns a ggplot / patchwork object
p_fts <- autoplot(res_iwt2)
expect_true(inherits(p_fts, "gg") || inherits(p_fts, "patchwork"))

# Snapshot the default plot (IWT — has pvalue_matrix so heatmap is shown)
expect_snapshot_plot(p_fts, "autoplot_ftwosample_iwt")

# plot.ftwosample is an alias for autoplot
p_fts2 <- plot(res_iwt2)
expect_true(inherits(p_fts2, "gg") || inherits(p_fts2, "patchwork"))

# TWT result (no pvalue_matrix) — different branch inside autoplot.ftwosample
p_fts_twt <- autoplot(res_twt2)
expect_true(inherits(p_fts_twt, "gg") || inherits(p_fts_twt, "patchwork"))
expect_snapshot_plot(p_fts_twt, "autoplot_ftwosample_twt")

# Custom arguments
p_fts_custom <- autoplot(
  res_iwt2,
  alpha1 = 0.1,
  alpha2 = 0.05,
  xrange = c(0, 1),
  ylabel = "Temperature",
  title = "IWT two-sample",
  linewidth = 0.5
)
expect_true(inherits(p_fts_custom, "gg") || inherits(p_fts_custom, "patchwork"))

# Swapped alpha (alpha1 < alpha2 triggers internal swap)
p_fts_swap <- autoplot(res_iwt2, alpha1 = 0.01, alpha2 = 0.05)
expect_true(inherits(p_fts_swap, "gg") || inherits(p_fts_swap, "patchwork"))

# alpha1 > all p-values triggers significance ribbon drawing (geom_rect path)
p_fts_sig <- autoplot(res_iwt2, alpha1 = 1.1, alpha2 = 0.9)
expect_true(inherits(p_fts_sig, "gg") || inherits(p_fts_sig, "patchwork"))

# alpha1 / alpha2 of length > 1 must error
expect_error(
  autoplot(res_iwt2, alpha1 = c(0.05, 0.1), alpha2 = 0.01),
  "single numeric value"
)
expect_error(
  autoplot(res_iwt2, alpha1 = 0.05, alpha2 = c(0.01, 0.02)),
  "single numeric value"
)

# ===========================================================================
# autoplot.fanova, plot.fanova
# ===========================================================================
p_fa <- autoplot(res_fanova)
expect_true(inherits(p_fa, "gg") || inherits(p_fa, "patchwork"))
expect_snapshot_plot(p_fa, "autoplot_fanova_default")

# plot_adjpval = TRUE (additional panels)
p_fa_adj <- autoplot(res_fanova, plot_adjpval = TRUE)
expect_true(inherits(p_fa_adj, "gg") || inherits(p_fa_adj, "patchwork"))
expect_snapshot_plot(p_fa_adj, "autoplot_fanova_adjpval")

# plot.fanova alias
p_fa2 <- plot(res_fanova)
expect_true(inherits(p_fa2, "gg") || inherits(p_fa2, "patchwork"))

# TWTaov result (no pval_matrix_F) — exercises null-matrix branch
p_fa_twt <- autoplot(res_twtaov)
expect_true(inherits(p_fa_twt, "gg") || inherits(p_fa_twt, "patchwork"))
expect_snapshot_plot(p_fa_twt, "autoplot_fanova_twt")

# plot_adjpval = TRUE with TWT (no heatmap available)
p_fa_twt_adj <- autoplot(res_twtaov, plot_adjpval = TRUE)
expect_true(inherits(p_fa_twt_adj, "gg") || inherits(p_fa_twt_adj, "patchwork"))

# Custom visual options
p_fa_custom <- autoplot(
  res_fanova,
  xrange = c(0, 1),
  alpha1 = 0.1,
  alpha2 = 0.05,
  ylabel = "Temp",
  title = "ANOVA",
  linewidth = 0.8,
  col = c("red", "blue")
)
expect_true(inherits(p_fa_custom, "gg") || inherits(p_fa_custom, "patchwork"))

# Two-factor ANOVA — exercises if(nvar > 1) F-test panel + per-factor loop
p_fa_2f <- autoplot(res_fanova2)
expect_true(inherits(p_fa_2f, "gg") || inherits(p_fa_2f, "patchwork"))

p_fa_2f_adj <- autoplot(res_fanova2, plot_adjpval = TRUE)
expect_true(inherits(p_fa_2f_adj, "gg") || inherits(p_fa_2f_adj, "patchwork"))

# 3-level factor x binary interaction ANOVA — exercises:
#   * nvar > 1 + !is.null(title): paste0(title, ": Functional Data...") (line 198)
#   * interaction branch of per-factor loop (lines 234–244)
#   * multi-column intersect in interaction: apply(colors,...) (line 242)
#   * multi-column non-interaction factor: apply(colors,...) (line 250)
groups3 <- factor(c(1L, 2L, 3L, 1L, 2L, 3L, 1L, 2L))
set.seed(42)
res_fanova_3l <- IWTaov(
  temperature ~ groups3 * grpB,
  B = 5L,
  method = "responses"
)
# title branch (line 198) + interaction/multi-column branches (234–244, 250)
p_fa_3l <- autoplot(res_fanova_3l, title = "3L")
expect_true(inherits(p_fa_3l, "gg") || inherits(p_fa_3l, "patchwork"))

# plot_adjpval = TRUE + title: covers paste0(title, ": Adjusted p-values - F-test") (line 301)
# and paste(title, ": Adjusted p-values - factor", ...) (line 335)
p_fa_3l_adj <- autoplot(res_fanova_3l, plot_adjpval = TRUE, title = "3L")
expect_true(inherits(p_fa_3l_adj, "gg") || inherits(p_fa_3l_adj, "patchwork"))

p_fa_twt_2f <- autoplot(res_twtaov2)
expect_true(inherits(p_fa_twt_2f, "gg") || inherits(p_fa_twt_2f, "patchwork"))

p_fa_twt_2f_adj <- autoplot(res_twtaov2, plot_adjpval = TRUE)
expect_true(inherits(p_fa_twt_2f_adj, "gg") || inherits(p_fa_twt_2f_adj, "patchwork"))

# Swapped alpha (alpha1 < alpha2 triggers internal swap)
p_fa_swap <- autoplot(res_fanova, alpha1 = 0.01, alpha2 = 0.05)
expect_true(inherits(p_fa_swap, "gg") || inherits(p_fa_swap, "patchwork"))

# alpha1 > all p-values triggers significance ribbon drawing (geom_rect path)
p_fa_sig <- autoplot(res_fanova, alpha1 = 1.1, alpha2 = 0.9)
expect_true(inherits(p_fa_sig, "gg") || inherits(p_fa_sig, "patchwork"))

# alpha1 / alpha2 of length > 1 must error
expect_error(
  autoplot(res_fanova, alpha1 = c(0.05, 0.1), alpha2 = 0.01),
  "single numeric value"
)
expect_error(
  autoplot(res_fanova, alpha1 = 0.05, alpha2 = c(0.01, 0.02)),
  "single numeric value"
)

# ===========================================================================
# autoplot.flm, plot.flm
# ===========================================================================
p_flm <- autoplot(res_flm)
expect_true(inherits(p_flm, "gg") || inherits(p_flm, "patchwork"))
expect_snapshot_plot(p_flm, "autoplot_flm_default")

# plot_adjpval TRUE
p_flm_adj <- autoplot(res_flm, plot_adjpval = TRUE)
expect_true(inherits(p_flm_adj, "gg") || inherits(p_flm_adj, "patchwork"))
expect_snapshot_plot(p_flm_adj, "autoplot_flm_adjpval")

# plot_adjpval = TRUE + title + ylim: covers title branches inside adjpval
# panels (lines 248, 307) and the coord_cartesian(ylim) branch (line 175)
p_flm_adj_t <- autoplot(
  res_flm,
  plot_adjpval = TRUE,
  title = "LM",
  ylim = c(-5, 5)
)
expect_true(inherits(p_flm_adj_t, "gg") || inherits(p_flm_adj_t, "patchwork"))

# plot.flm alias
p_flm2 <- plot(res_flm)
expect_true(inherits(p_flm2, "gg") || inherits(p_flm2, "patchwork"))

# TWT result (no pval_matrix_F)
p_flm_twt <- autoplot(res_twtlm)
expect_true(inherits(p_flm_twt, "gg") || inherits(p_flm_twt, "patchwork"))
expect_snapshot_plot(p_flm_twt, "autoplot_flm_twt")

# plot_adjpval = TRUE with TWT
p_flm_twt_adj <- autoplot(res_twtlm, plot_adjpval = TRUE)
expect_true(
  inherits(p_flm_twt_adj, "gg") || inherits(p_flm_twt_adj, "patchwork")
)

# Custom visual options
p_flm_custom <- autoplot(
  res_flm,
  xrange = c(0, 1),
  alpha1 = 0.1,
  alpha2 = 0.05,
  ylabel = "Temp",
  title = "LM",
  linewidth = 0.8
)
expect_true(inherits(p_flm_custom, "gg") || inherits(p_flm_custom, "patchwork"))

# Two-predictor LM — exercises longer per-variable loop (intercept + grpA + grpB)
p_flm_2v <- autoplot(res_flm2)
expect_true(inherits(p_flm_2v, "gg") || inherits(p_flm_2v, "patchwork"))

p_flm_2v_adj <- autoplot(res_flm2, plot_adjpval = TRUE)
expect_true(inherits(p_flm_2v_adj, "gg") || inherits(p_flm_2v_adj, "patchwork"))

p_flm_twt_2v <- autoplot(res_twtlm2)
expect_true(inherits(p_flm_twt_2v, "gg") || inherits(p_flm_twt_2v, "patchwork"))

p_flm_twt_2v_adj <- autoplot(res_twtlm2, plot_adjpval = TRUE)
expect_true(inherits(p_flm_twt_2v_adj, "gg") || inherits(p_flm_twt_2v_adj, "patchwork"))

# Swapped alpha (alpha1 < alpha2 triggers internal swap)
p_flm_swap <- autoplot(res_flm, alpha1 = 0.01, alpha2 = 0.05)
expect_true(inherits(p_flm_swap, "gg") || inherits(p_flm_swap, "patchwork"))

# alpha1 > all p-values triggers significance ribbon drawing (geom_rect path)
p_flm_sig <- autoplot(res_flm, alpha1 = 1.1, alpha2 = 0.9)
expect_true(inherits(p_flm_sig, "gg") || inherits(p_flm_sig, "patchwork"))

# alpha1 / alpha2 of length > 1 must error
expect_error(
  autoplot(res_flm, alpha1 = c(0.05, 0.1), alpha2 = 0.01),
  "single numeric value"
)
expect_error(
  autoplot(res_flm, alpha1 = 0.05, alpha2 = c(0.01, 0.02)),
  "single numeric value"
)

# ===========================================================================
# plot.IWT1 — base graphics; redirect to null device to avoid interactive prompt
# ===========================================================================
tmp_pdf <- tempfile(fileext = ".pdf")
grDevices::pdf(tmp_pdf)
plot(res_iwt1, xrange = c(0, 1)) # default: mu is zero
grDevices::dev.off()
expect_true(file.exists(tmp_pdf))

# vector mu — exercises the else branch (lines 184-185) of the mu length check
set.seed(42)
res_iwt1_vmu <- IWT1(d1, mu = colMeans(d1), B = 5L)
grDevices::pdf(tempfile(fileext = ".pdf"))
plot(res_iwt1_vmu, xrange = c(0, 1))
grDevices::dev.off()

# Swapped alpha (alpha1 < alpha2 triggers internal swap)
grDevices::pdf(tempfile(fileext = ".pdf"))
plot(res_iwt1, xrange = c(0, 1), alpha1 = 0.01, alpha2 = 0.05)
grDevices::dev.off()

# alpha1 > all p-values triggers difference1/difference2 rectangle drawing
grDevices::pdf(tempfile(fileext = ".pdf"))
plot(res_iwt1, xrange = c(0, 1), alpha1 = 1.1, alpha2 = 0.9)
grDevices::dev.off()

# alpha1 / alpha2 of length > 1 must error
expect_error(
  plot(res_iwt1, xrange = c(0, 1), alpha1 = c(0.05, 0.1), alpha2 = 0.01),
  "single numeric value"
)
expect_error(
  plot(res_iwt1, xrange = c(0, 1), alpha1 = 0.05, alpha2 = c(0.01, 0.02)),
  "single numeric value"
)

# Clean up
unlink(tmp_pdf)
