# Two population Interval Wise Testing procedure

The function implements the Interval Wise Testing procedure for testing
mean differences between two functional populations. Functional data are
tested locally and unadjusted and adjusted p-value functions are
provided. The unadjusted p-value function controls the point-wise error
rate. The adjusted p-value function controls the interval-wise error
rate.

## Usage

``` r
ITP2bspline(
  data1,
  data2,
  mu = 0,
  B = 1000,
  paired = FALSE,
  order = 2,
  nknots = dim(data1)[2]
)

ITP2fourier(
  data1,
  data2,
  mu = 0,
  B = 1000,
  paired = FALSE,
  maxfrequency = floor(dim(data1)[2]/2)
)

ITP2pafourier(
  data1,
  data2,
  mu = 0,
  B = 1000,
  paired = FALSE,
  maxfrequency = floor(dim(data1)[2]/2)
)

IWT2(
  data1,
  data2,
  mu = 0,
  dx = NULL,
  B = 1000L,
  paired = FALSE,
  alternative = c("two.sided", "less", "greater"),
  verbose = FALSE,
  recycle = TRUE
)
```

## Arguments

- data1:

  Either a numeric matrix or an object of class
  [`fda::fd`](https://rdrr.io/pkg/fda/man/fd.html) specifying the data
  in the first sample. If the data is provided within a matrix, it
  should be of shape \\n_1 \times J\\ and it should contain in each row
  one of the \\n_1\\ functions in the sample and in columns the
  evaluation of each function on a **same** uniform grid of size \\J\\.

- data2:

  Either a numeric matrix or an object of class
  [`fda::fd`](https://rdrr.io/pkg/fda/man/fd.html) specifying the data
  in the second sample. If the data is provided within a matrix, it
  should be of shape \\n_2 \times J\\ and it should contain in each row
  one of the \\n_2\\ functions in the sample and in columns the
  evaluation of each function on a **same** uniform grid of size \\J\\.

- mu:

  Either a numeric value or a numeric vector or an object of class
  [`fda::fd`](https://rdrr.io/pkg/fda/man/fd.html) specifying the
  functional mean difference under the null hypothesis. If `mu` is a
  constant, then a constant function is used. If `mu` is a numeric
  vector, it must correspond to evaluation of the mean difference
  function on the **same** grid that has been used to evaluate the data
  samples. Defaults to `0`.

- B:

  An integer value specifying the number of iterations of the MC
  algorithm to evaluate the p-value of the permutation tests. Defaults
  to `1000L`.

- paired:

  A boolean value specifying whether a paired test should be performed.
  Defaults to `FALSE`.

- order:

  Order of the B-spline basis expansion. Defaults to `2L`.

- nknots:

  Number of knots of the B-spline basis expansion. Defaults to
  `dim(data1)[2]`.

- maxfrequency:

  The maximum frequency to be used in the Fourier basis expansion of
  data. Defaults to `floor(dim(data1)[2] / 2)`, leading to an
  interpolating expansion.

- dx:

  A numeric value specifying the discretization step of the grid used to
  evaluate functional data when it is provided as objects of class
  [`fda::fd`](https://rdrr.io/pkg/fda/man/fd.html). Defaults to `NULL`,
  in which case a default value of `0.01` is used which corresponds to a
  grid of size `100L`. Unused if functional data is provided in the form
  of matrices.

- alternative:

  A string specifying the type of alternative hypothesis. Choices are
  `"two.sided"`, `"less"` or `"greater"`. Defaults to `"two.sided"`.

- verbose:

  A boolean value specifying whether to print the progress of the
  computation. Defaults to `FALSE`.

- recycle:

  A boolean value specifying whether the recycled version of the
  interval-wise testing procedure should be used. See Pini and
  Vantini (2017) for details. Defaults to `TRUE`.

## Value

An object of class `ftwosample` containing the following components:

- `data`: A numeric matrix of shape \\n \times J\\ containing the
  evaluation of the \\n = n_1 + n_2\\ functions on a **common** uniform
  grid of size \\p\\.

- `group_labels`: An integer vector of size \\n = n_1 + n_2\\ containing
  the group membership of each function.

- `mu`: A numeric vector of shape \\J\\ containing the evaluation of the
  functional mean difference under the null hypothesis on the same
  uniform grid used to evaluate the functional samples.

- `unadjusted_pvalues`: A numeric vector of size \\J\\ containing the
  evaluation of the unadjusted p-value function on the **same** uniform
  grid used to evaluate the functional samples.

- `adjusted_pvalues`: A numeric vector of size \\J\\ containing the
  evaluation of the adjusted p-value functione on the **same** uniform
  grid used to evaluate the functional samples.

Optionally, the list may contain the following components:

- `global_pvalue`: A numeric value containing the global p-value. Only
  present if the `correction` argument is set to `"Global"`.

- `pvalue_matrix`: A numeric matrix of shape \\p \times p\\ containing
  the p-values of the interval-wise tests. Element \\i, j\\ contains the
  p-value of the test performed on the interval indexed by \\j, j+1 ,
  \dots, j+(p-i)\\. Only present if the `correction` argument is set to
  `"IWT"`.

## References

A. Pini and S. Vantini (2017). The Interval Testing Procedure: Inference
for Functional Data Controlling the Family Wise Error Rate on Intervals.
*Biometrics*, 73(3): 835–845.

A. Pini and S. Vantini (2017). Interval-wise testing for functional
data. *Journal of Nonparametric Statistics*, 29(2), 407-424.

## See also

See also
[`plot.ftwosample()`](https://permaverse.github.io/fdatest/reference/plot.ftwosample.md)
for plotting the results.

## Examples

``` r
# Performing the IWT for two populations
IWT_result <- IWT2(NASAtemp$paris, NASAtemp$milan, B = 10L)

# Plotting the results of the IWT
plot(
  IWT_result,
  xrange = c(0, 12),
  title = 'IWT results for testing mean differences'
)


# Plotting the p-value heatmap
IWTimage(IWT_result, abscissa_range = c(0, 12))

# Selecting the significant components at 5% level
which(IWT_result$adjusted_pvalues < 0.05)
#>   [1]  81  87  88  89  90  91  92  93  94  95  96  97  98 100 101 102 103 104
#>  [19] 105 106 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126
#>  [37] 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144
#>  [55] 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161 162
#>  [73] 163 164 165 166 167 168 169 170 171 172 173 174 175 176 177 178 179 180
#>  [91] 181 182 183 184 185 186 187 188 189 190 191 192 193 194 195 196 197 198
#> [109] 199 200 201 202 203 204 205 206 207 208 209 210 211 212 213 214 215 216
#> [127] 217 218 219 220 221 222 223 224 225 226 227 228 229 230 231 232 233 234
#> [145] 235 236 237 238 240 241 242 243 253 254 255 257 258 259 260 261 262 263
#> [163] 264 265 266 267 268 269 270 272 273 274 275 276 277
```
