# Joint Test via Contrast Statement

Joint Test via Contrast Statement

## Usage

``` r
joint_test(model = NULL, coef_df_input = NULL, cov_mtx_input = NULL, L, C)
```

## Arguments

- model:

  Model object that contrast statement applies to

- coef_df_input:

  Coefficients for contrast

- cov_mtx_input:

  Covariance matrix for contrast

- L:

  Vector or matrix corresponding to contrast coefficients

- C:

  Vector of contrasts

## Value

Tibble of contrast statement results containing the point estimate, test
statistic, and p-value
