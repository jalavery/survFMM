#' Simulated Dataset
#'
#' A dataset containing outcome data for 1500 simulated observations
#' belonging to three underlying latent subgroups
#'
#' @format A data table of simulated data with:
#' \describe{
#'   \item{record_id}{Record ID}
#'   \item{covariate_sim_normal}{Simulated covariate following a N(0, 1)
#'   distribution}
#'   \item{covariate_sim_binary1}{Simulated binary covariate with p=0.5}
#'   \item{covariate_sim_binary2}{Simulated binary covariate with p=0.4}
#'   \item{tx}{Binary treatment indicator}
#'   \item{latent_subgroup}{Simulated latent subgroup}
#'   \item{event_status}{Event indicator, 0 = censored, 1 = event}
#'   \item{time_to_event_days}{Simulated time to event or censoring}
#'   \item{iptw_trim97}{Inverse probability of treatment weights (IPTW), trimmed
#'   at the 97th percentile}
#'   \item{ipcw_trim97}{Inverse probability of censoring weights (IPCW), trimmed
#'   at the 97th percentile; only populated for observations that are not censored}
#'   \item{iptw_ipcw_trim97}{Product of IPTW and inverse probability of censoring weights (IPCW), trimmed
#'   at the 97th percentile; only populated for observations that are not censored}

#'   ...
#' }
"sim_data"
