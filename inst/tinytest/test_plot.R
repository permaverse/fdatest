# Tests for plot and autoplot methods:
# R/plot-ftwosample.R  (autoplot.ftwosample / plot.ftwosample)
# R/plot-fanova.R      (autoplot.fanova    / plot.fanova)
# R/plot-flm.R         (autoplot.flm       / plot.flm)
# R/plot.IWT1.R        (plot.IWT1)
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

# ===========================================================================
# autoplot.ftwosample + plot.ftwosample
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

# ===========================================================================
# autoplot.fanova + plot.fanova
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

# ===========================================================================
# autoplot.flm + plot.flm
# ===========================================================================
p_flm <- autoplot(res_flm)
expect_true(inherits(p_flm, "gg") || inherits(p_flm, "patchwork"))
expect_snapshot_plot(p_flm, "autoplot_flm_default")

# plot_adjpval = TRUE
p_flm_adj <- autoplot(res_flm, plot_adjpval = TRUE)
expect_true(inherits(p_flm_adj, "gg") || inherits(p_flm_adj, "patchwork"))
expect_snapshot_plot(p_flm_adj, "autoplot_flm_adjpval")

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

# ===========================================================================
# plot.IWT1 — base graphics; redirect to null device to avoid interactive prompt
# ===========================================================================
tmp_pdf <- tempfile(fileext = ".pdf")
grDevices::pdf(tmp_pdf)
plot(res_iwt1, xrange = c(0, 1)) # default: mu = 0 (scalar)
grDevices::dev.off()
expect_true(file.exists(tmp_pdf))

# Clean up
unlink(tmp_pdf)
