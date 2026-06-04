# EM Algorithm for Finite Mixture Models with Survival Endpoints

EM Algorithm for Finite Mixture Models with Survival Endpoints

## Usage

``` r
fmm_em_algorithm(
  input_df,
  weights_input,
  starting_values_input_df,
  k,
  outc_model_formula,
  outc_model_time,
  outc_model_status,
  outc_model_covars,
  outc_distribution,
  covariates_subgroup_model,
  n_inits = 5,
  tolerance = 0.001,
  conv_pct_criteria = -1,
  max_iter = 200
)
```

## Arguments

- input_df:

  Input data frame containing 1 row/observation, along with each
  observation's event status (0=censored, 1=event) and event time
  variables

- weights_input:

  Variable for inverse probability of treatment weights, if applicable.
  For IPCW-FMM, supply the IPCW. For IPCW-FMM that also uses IPTW, the
  product of the ITPW and IPCW may be supplied.

- starting_values_input_df:

  Input dataset with starting values for algorithm

- k:

  Number of subgroups

- outc_model_formula:

  Formula object for subgroup-specific outcome models

- outc_model_time:

  Variable indicating the time to event or censoring for each
  observation.

- outc_model_status:

  Variable indicating the event status for the time-to-event outcome. It
  is assumed that 0 = censored, 1 = event.

- outc_model_covars:

  Names of covariates to include in the outcome models for each subgroup

- outc_distribution:

  Outcome distribution for subgroup-specific outcome models. Currently
  allowed values are "Weibull" and "Log-Normal" (not case-sensitive)

- covariates_subgroup_model:

  Vector of covariates to include in subgroup membership model

- n_inits:

  Number of initial partitions for the EM algorithm. Default is 5. A
  higher number of initial partitions may result in greater stability of
  estimates.

- tolerance:

  Convergence criteria for the change in the log-likelihood for the EM
  algorithm.

- conv_pct_criteria:

  Convergence criteria for the percentage of observations changing
  subgroup. Specify -1 to only use the log-likelihood as the convergence
  criteria.

- max_iter:

  Maximum number of iterations. Default is 200.

## Value

List
