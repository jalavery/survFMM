#' Log-likelihood computation for use within the survFMM function
#'
#' @param input_df Input data frame containing each observation's survival status (0=censored, 1=event) and time
#' @param distribution Acceptable values are currently Weibull or log-normal (not case sensitive)
#'
#' @returns Log-likelihood value
#'
#' @keywords internal
log_likelihood_with_censoring <- function(input_df,
                                          distribution) {
  if (stringr::str_to_lower(distribution) == "weibull") {
    input_df %>%
      dplyr::mutate(
        survival_distribution = .data$status * stats::dweibull(.data$time,
          shape = .data$shape,
          scale = .data$scale * exp(.data$beta_variable)
        ),
        censoring_distribution = (1 - .data$status) * (1 - stats::pweibull(.data$time,
          shape = .data$shape,
          scale = .data$scale * exp(.data$beta_variable)
        )),
        prior_prob_survival_cens = .data$prior_probability * (.data$survival_distribution + .data$censoring_distribution)
      ) %>%
      # sum across components w/in records
      dplyr::group_by(.data$record_id) %>%
      dplyr::summarize(
        sum_prior_dists = sum(.data$prior_prob_survival_cens),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        log_sum_prior_dists = log(.data$sum_prior_dists)
        # weights_log_sum_prior_dists = weights_input_loglik*log_sum_prior_dists
      ) %>%
      # sum across records
      dplyr::summarize(
        log_likelihood_values = sum(.data$log_sum_prior_dists)
        # log_likelihood_values_iptw = sum(weights_log_sum_prior_dists)
      )
  } else if (stringr::str_to_lower(distribution) == "lognormal") {
    input_df %>%
      dplyr::mutate(
        survival_distribution = .data$status * stats::dlnorm(.data$time,
          sdlog = .data$shape,
          meanlog = .data$scale + .data$beta_variable
        ),
        censoring_distribution = (1 - .data$status) * (1 - stats::plnorm(.data$time,
          sdlog = .data$shape,
          meanlog = .data$scale + .data$beta_variable
        )),
        prior_prob_survival_cens = .data$prior_probability * (.data$survival_distribution + .data$censoring_distribution)
      ) %>%
      dplyr::group_by(.data$record_id) %>%
      dplyr::summarize(
        sum_prior_dists = sum(.data$prior_prob_survival_cens),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        log_sum_prior_dists = log(.data$sum_prior_dists)
        # weights_log_sum_prior_dists = weights_input_loglik*log_sum_prior_dists
      ) %>%
      dplyr::summarize(
        log_likelihood_values = sum(.data$log_sum_prior_dists)
        # log_likelihood_values_iptw = sum(weights_log_sum_prior_dists)
      )
  } # end of else if for lognormal
} # end of function
