#' Log-likelihood computation for use within the survFMM function
#'
#' @param input_df Input data frame containing each observation's survival status (0=censored, 1=event) and time
#' @param distribution Acceptable values are currently Weibull or log-normal (not case sensitive)
#'
#' @returns Log-likelihood value
#'
log_likelihood_with_censoring <- function(input_df,
                                          distribution) {
  if (str_to_lower(distribution) == "weibull") {
    input_df %>%
      mutate(
        survival_distribution = status * dweibull(time,
          shape = shape,
          scale = scale * exp(beta_variable)
        ),
        censoring_distribution = (1 - status) * (1 - pweibull(time,
          shape = shape,
          scale = scale * exp(beta_variable)
        )),
        prior_prob_survival_cens = prior_probability * (survival_distribution + censoring_distribution)
      ) %>%
      # sum across components w/in records
      group_by(record_id) %>%
      summarize(
        sum_prior_dists = sum(prior_prob_survival_cens),
        .groups = "drop"
      ) %>%
      mutate(
        log_sum_prior_dists = log(sum_prior_dists)
        # weights_log_sum_prior_dists = weights_input_loglik*log_sum_prior_dists
      ) %>%
      # sum across records
      summarize(
        log_likelihood_values = sum(log_sum_prior_dists)
        # log_likelihood_values_iptw = sum(weights_log_sum_prior_dists)
      )
  } else if (str_to_lower(distribution) == "lognormal") {
    input_df %>%
      mutate(
        survival_distribution = status * dlnorm(time,
          sdlog = shape,
          meanlog = scale + beta_variable
        ),
        censoring_distribution = (1 - status) * (1 - plnorm(time,
          sdlog = shape,
          meanlog = scale + beta_variable
        )),
        prior_prob_survival_cens = prior_probability * (survival_distribution + censoring_distribution)
      ) %>%
      group_by(record_id) %>%
      summarize(
        sum_prior_dists = sum(prior_prob_survival_cens),
        .groups = "drop"
      ) %>%
      mutate(
        log_sum_prior_dists = log(sum_prior_dists)
        # weights_log_sum_prior_dists = weights_input_loglik*log_sum_prior_dists
      ) %>%
      summarize(
        log_likelihood_values = sum(log_sum_prior_dists)
        # log_likelihood_values_iptw = sum(weights_log_sum_prior_dists)
      )
  } # end of else if for lognormal
} # end of function
