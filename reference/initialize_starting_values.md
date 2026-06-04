# Initialize starting values for each initial partition of the EM algorithm call

Initialize starting values for each initial partition of the EM
algorithm call

## Usage

``` r
initialize_starting_values(
  n_inits,
  k,
  starting_values_type,
  starting_values_df = NULL,
  starting_values_window,
  input_df,
  outc_model_formula,
  weights_input,
  outc_distribution,
  starting_scale_logn = "exp"
)
```

## Arguments

- n_inits:

  Number of initial partitions

- k:

  Number of subgroups

- starting_values_type:

  One of "single_survreg", "uniform_pct", or "non_random_start".
  "single_survreg" fits a single AFT model to all of the data and then
  generates random starting values based on the
  \`starting_values_window\` parameter for each initial partition. If
  not supplying starting values, they are randomly generated for each
  initial partition. Be sure to set a seed at the top of your script to
  ensure reproducibility.

- starting_values_df:

  Optional input dataset with starting values for algorithm

- starting_values_window:

  The percent margin around the starting values. For example,
  starting_values_type = 'single_survreg' and starting_values_window =
  0.5 means that starting values are randomly generated uniformly +/- 50
  fit. For starting_values_type = 'uniform' and starting_values_df is
  supplied, the starting values are generated uniformly

- input_df:

  Input data frame containing 1 row/observation, along with each
  observation's event status (0=censored, 1=event) and event time
  variables

- outc_model_formula:

  Formula object for subgroup-specific outcome models

- weights_input:

  Variable for inverse probability of treatment weights, if applicable.
  For IPCW-FMM, supply the IPCW. For IPCW-FMM that also uses IPTW, the
  product of the ITPW and IPCW may be supplied.

- outc_distribution:

  Outcome distribution for subgroup-specific outcome models. Currently
  allowed values are "Weibull" and "Log-Normal" (not case-sensitive)

## Value

List of starting values, with 1 list element per set of starting values
