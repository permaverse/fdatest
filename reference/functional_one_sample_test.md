# Local testing procedures for the functional one-sample test

The function implements local testing procedures for testing the center
of symmetry of a functional population. Functional data are tested
locally and unadjusted and adjusted p-value functions are provided. The
unadjusted p-value function controls the point-wise error rate. The
adjusted p-value function can be computed according to the following
methods:

- interval-wise testing (controlling the interval-wise error rate)

## Usage

``` r
functional_one_sample_test(
  data,
  correction = c("IWT"),
  mu = 0,
  dx = NULL,
  n_perm = 1000L,
  verbose = FALSE,
  aggregation_strategy = c("integral", "max"),
  recycle = TRUE
)
```

## Arguments

- data:

  Either a numeric matrix or an object of class
  [`fda::fd`](https://rdrr.io/pkg/fda/man/fd.html) specifying the sample
  data. If the data is provided within a matrix, it should be of shape
  \\n \times J\\ and it should contain in each row one of the \\n\\
  functions in the sample and in columns the evaluation of each function
  on a **same** uniform grid of size \\J\\.

- correction:

  A string specifying the correction method to perform the local
  functional testing procedure and adjust the p-value function.
  Currently only `"IWT"` is available.

- mu:

  Either a numeric value or a numeric vector or an object of class
  [`fda::fd`](https://rdrr.io/pkg/fda/man/fd.html) specifying the
  functional center of symmetry under the null hypothesis. If `mu` is a
  constant, then a constant function is used. If `mu` is a numeric
  vector, it must correspond to evaluations of the center function on
  the **same** grid that has been used to evaluate the data sample.
  Defaults to `0`.

- dx:

  A numeric value specifying the step of the uniform grid on which the
  data are evaluated. If `NULL`, the step is automatically inferred from
  the data. Defaults to `NULL`.

- n_perm:

  An integer value specifying the number of permutations to use for the
  local testing procedure. Defaults to `1000L`.

- verbose:

  A boolean value specifying whether to print the progress of the
  computation. Defaults to `FALSE`.

- aggregation_strategy:

  A string specifying the strategy to aggregate the point-wise test
  statistics for the correction procedure. Possible values are
  `"integral"` and `"max"`. Defaults to `"integral"`.

- recycle:

  A boolean value specifying whether to recycle the test statistic
  values across permutations for the IWT procedure. Defaults to `TRUE`.

## Value

An object of class `fos` containing the following components:

- `data`: A numeric matrix of shape \\n \times J\\ containing the
  evaluation of the \\n\\ functions on a uniform grid of size \\J\\.

- `mu`: A numeric vector of shape \\J\\ containing the evaluation of the
  functional center of symmetry under the null hypothesis on the same
  uniform grid used to evaluate the functional sample.

- `unadjusted_pvalues`: A numeric vector of size \\J\\ containing the
  evaluation of the unadjusted p-value function on the same uniform grid
  used to evaluate the functional sample.

- `adjusted_pvalues`: A numeric vector of size \\J\\ containing the
  evaluation of the adjusted p-value function on the same uniform grid
  used to evaluate the functional sample.

- `correction_method`: A string containing the correction method used to
  compute the adjusted p-value function.

Optionally, the list may contain the following component:

- `pvalue_matrix`: A numeric matrix of shape \\p \times p\\ containing
  the p-values of the interval-wise tests. Element \\i, j\\ contains the
  p-value of the test performed on the interval indexed by \\j, j+1,
  \dots, j+(p-i)\\. Only present if the `correction` argument is set to
  `"IWT"`.

## References

For the interval-wise testing procedure:

- Pini, Alessia, and Simone Vantini. 2016. "The interval testing
  procedure: a general framework for inference in functional data
  analysis." Biometrics 72 (3): 835–845.

- Pini, Alessia, and Simone Vantini. 2017. "Interval-Wise Testing for
  Functional Data." Journal of Nonparametric Statistics 29 (2): 407–24.

## See also

[`iwt1()`](https://permaverse.github.io/fdatest/reference/IWT1.md) for
calling directly the IWT test and
[`plot.fos()`](https://permaverse.github.io/fdatest/reference/plot.fos.md)
and
[`autoplot.fos()`](https://permaverse.github.io/fdatest/reference/plot.fos.md)
for plotting the results.

## Examples

``` r
# Performing the IWT for one population
IWT_result <- functional_one_sample_test(
  NASAtemp$paris, mu = 4, n_perm = 10L
)

# Plotting the results of the IWT
IWTimage(IWT_result, abscissa_range = c(0, 12))


# Selecting the significant components at 5% level
which(IWT_result$adjusted_pvalues < 0.05)
#>   [1]   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18
#>  [19]  19  20  21  22  23  24  25  26  27  28  29  30  31  32  33  34  35  36
#>  [37]  37  38  39  40  41  42  43  44  45  46  47  48  49  50  51  52  53  54
#>  [55]  55  56  57  58  59  60  61  62  63  64  65  66  67  68  69  70  71  72
#>  [73]  73  74  75  76  77  78  79  80  81  82  83  84  85  86  87  88  89  90
#>  [91]  91  92  93  94  95  96  97  98  99 100 101 102 103 104 105 106 107 108
#> [109] 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126
#> [127] 127 128 129 130 131 132 133 262 263 264 265 266 267 268 269 270 271 272
#> [145] 273 274 275 276 277 278 279 280 281 282 283 284 285 286 287 288 289 290
#> [163] 291 292 293 294 295 296 297 298 299 300 301 302 303 304 305 306 307 308
#> [181] 309 310 311 312 313 314 315 316 317 318 319 320 321 322 323 324 325 326
#> [199] 327 328 329 330 331 332 333 334 335 336 337 338 339 340 341 342 343 344
#> [217] 345 346 347 348 349 350 351 352 353 354 355 356 357 358 359 360 361 362
#> [235] 363 364 365
```
