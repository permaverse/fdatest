# Two population Threshold Wise Testing procedure

The function implements the Threshold Wise Testing procedure for testing
mean differences between two functional populations. Functional data are
tested locally and unadjusted and adjusted p-value functions are
provided. The unadjusted p-value function controls the point-wise error
rate. The adjusted p-value function controls the family-wise error rate
asymptotically.

## Usage

``` r
TWT2(
  data1,
  data2,
  mu = 0,
  dx = NULL,
  B = 1000L,
  paired = FALSE,
  alternative = c("two.sided", "less", "greater"),
  verbose = FALSE
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

- dx:

  A numeric value specifying the discretization step of the grid used to
  evaluate functional data when it is provided as objects of class
  [`fda::fd`](https://rdrr.io/pkg/fda/man/fd.html). Defaults to `NULL`,
  in which case a default value of `0.01` is used which corresponds to a
  grid of size `100L`. Unused if functional data is provided in the form
  of matrices.

- B:

  An integer value specifying the number of iterations of the MC
  algorithm to evaluate the p-value of the permutation tests. Defaults
  to `1000L`.

- paired:

  A boolean value specifying whether a paired test should be performed.
  Defaults to `FALSE`.

- alternative:

  A string specifying the type of alternative hypothesis. Choices are
  `"two.sided"`, `"less"` or `"greater"`. Defaults to `"two.sided"`.

- verbose:

  A boolean value specifying whether to print the progress of the
  computation. Defaults to `FALSE`.

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

Abramowicz, K., Pini, A., Schelin, L., Stamm, A., & Vantini, S. (2022).
“Domain selection and familywise error rate for functional data: A
unified framework. *Biometrics* 79(2), 1119-1132.

Pini, A., & Vantini, S. (2017). Interval-wise testing for functional
data. *Journal of Nonparametric Statistics*, 29(2), 407-424.

## See also

See also
[`plot.ftwosample()`](https://permaverse.github.io/fdatest/reference/plot.ftwosample.md)
for plotting the results.

## Examples

``` r
# Performing the TWT for two populations
TWT_result <- TWT2(NASAtemp$paris, NASAtemp$milan)

# Plotting the results of the TWT
plot(
  TWT_result,
  xrange = c(0, 12),
  title = 'TWT results for testing mean differences'
)


# Selecting the significant components at 5% level
which(TWT_result$adjusted_pvalues < 0.05)
#>   [1]  49  50  61  64  65  69  70  71  72  73  74  88  89  90  91  92  93  94
#>  [19]  95  96  97 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115
#>  [37] 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133
#>  [55] 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 151
#>  [73] 152 153 154 155 156 157 158 159 160 161 162 163 164 165 166 167 168 169
#>  [91] 170 171 172 173 174 175 176 177 178 179 180 181 182 183 184 185 186 187
#> [109] 188 189 190 191 192 193 194 195 196 197 198 199 200 201 202 203 204 205
#> [127] 206 207 208 209 210 211 212 213 214 215 216 217 218 219 220 221 222 223
#> [145] 224 225 226 227 228 229 230 231 232 233 234 235 236 237 238 239 240 241
#> [163] 242 243 244 245 246 247 248 249 250 251 252 253 254 255 256 257 258 259
#> [181] 260 261 262 264 265 266 267 270 271 272 273 274 275 276 281 286 287 288
#> [199] 289 290 291 299
```
