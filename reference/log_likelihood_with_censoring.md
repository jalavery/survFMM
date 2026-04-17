# Log-likelihood computation for use within the survFMM function

Log-likelihood computation for use within the survFMM function

## Usage

``` r
log_likelihood_with_censoring(input_df, distribution)
```

## Arguments

- input_df:

  Input data frame containing each observation's survival status
  (0=censored, 1=event) and time

- distribution:

  Acceptable values are currently Weibull or log-normal (not case
  sensitive)

## Value

Log-likelihood value
