# Simulated Dataset

A dataset containing outcome data for 1500 simulated observations
belonging to three underlying latent subgroups

## Usage

``` r
sim_data
```

## Format

A data table of simulated data with:

- record_id:

  Record ID

- covariate_sim_normal:

  Simulated covariate following a N(0, 1) distribution

- covariate_sim_binary1:

  Simulated binary covariate with p=0.5

- covariate_sim_binary2:

  Simulated binary covariate with p=0.4

- tx:

  Binary treatment indicator

- latent_subgroup:

  Simulated latent subgroup

- event_status:

  Event indicator, 0 = censored, 1 = event

- time_to_event_days:

  Simulated time to event or censoring

- iptw_trim97:

  Inverse probability of treatment weights (IPTW), trimmed at the 97th
  percentile

- ipcw_trim97:

  Inverse probability of censoring weights (IPCW), trimmed at the 97th
  percentile; only populated for observations that are not censored

- iptw_ipcw_trim97:

  Product of IPTW and inverse probability of censoring weights (IPCW),
  trimmed at the 97th percentile; only populated for observations that
  are not censored
