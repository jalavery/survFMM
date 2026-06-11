#' Initialize starting values for each initial partition of the EM algorithm
#' call
#'
#' @param n_inits Number of initial partitions
#' @param k Number of subgroups
#' @param starting_values_type One of "single_survreg", "uniform_pct", or
#'   "non_random_start". "single_survreg" fits a single AFT model to all of the
#'   data and then generates random starting values based on the
#'   `starting_values_window` parameter for each initial partition. If not
#'   supplying starting values, they are randomly generated for each initial
#'   partition. Be sure to set a seed at the top of your script to ensure
#'   reproducibility.
#' @param starting_values_df Optional input dataset with starting values for algorithm
#' @param starting_values_window The percent margin around the starting values.
#'   For example, starting_values_type = 'single_survreg' and
#'   starting_values_window = 0.5 means that starting values are randomly
#'   generated uniformly +/- 50% of the estimates from a single survreg model
#'   fit. For starting_values_type = 'uniform' and starting_values_df is supplied, the starting values are generated uniformly
#' @param input_df Input data frame containing 1 row/observation, along with
#'   each observation's event status (0=censored, 1=event) and event time
#'   variables
#' @param outc_model_formula Formula object for subgroup-specific outcome models
#' @param weights_input Variable for inverse probability of treatment weights,
#'   if applicable. For IPCW-FMM, supply the IPCW. For IPCW-FMM that also uses
#'   IPTW, the product of the ITPW and IPCW may be supplied.
#' @param outc_distribution Outcome distribution for subgroup-specific outcome
#'   models. Currently allowed values are "Weibull" and "Log-Normal" (not
#'   case-sensitive)
#'
#' @returns List of starting values, with 1 list element per set of starting values
#' @export
initialize_starting_values <- function(n_inits,
                                       k,
                                       starting_values_type,
                                       starting_values_df = NULL,
                                       starting_values_window,
                                       input_df,
                                       outc_model_formula,
                                       weights_input,
                                       outc_distribution,
                                       starting_scale_logn = "exp") {
  # check inputs ---------------------------------------------------------
  if (starting_values_window == 0 & n_inits > 1){
    message("Note: Since no noise is added to the starting values, the only variation across initial partition is random subgroup initialization.")
  }
  if (!is.null(starting_values_df) & starting_values_type == "single_survreg"){
    message("Note: `starting_values_df` is ignored when `starting_values_type` = 'single_survreg'. To use the supplied starting values, change `starting_values_type` to `uniform_pct`.")
  }
  if (is.null(starting_values_df) & starting_values_type == "uniform_pct"){
    stop("`starting_values_df` is required when `starting_values_type` = 'uniform_pct'.")
  }

  # starting values ---------------------------------------------------------
  # generate n_inits different starting values
  # set up random starting values based on the true data
  if (starting_values_type == "single_survreg") {
    survreg_for_start <- survival::survreg(outc_model_formula,
      data = input_df,
      weights = get(weights_input),
      robust = TRUE,
      dist = outc_distribution
    )

    # browser()

    survreg_for_start_tidy <- survreg_for_start %>%
      broom::tidy() %>%
      # with flexsurvreg, shape and scale parameters are returned in a parameterization consistent with rweibull
      # with survreg, the log(scale) term returned is 1/rweibull shape
      # with survreg, the exp(survreg intercept) = rweib scale
      # this is documented in the survreg example
      #   survreg's scale  =    1/(rweibull shape)
      #   survreg's intercept = log(rweibull scale)
      dplyr::rename(
        estimate_original = .data$estimate,
        term_original = .data$term
      ) %>%
      dplyr::mutate(
        estimate = dplyr::case_when(
          outc_distribution == "weibull" & term_original == "Log(scale)" ~ 1 / exp(estimate_original),
          outc_distribution == "weibull" & term_original == "(Intercept)" ~ exp(estimate_original),
          outc_distribution == "lognormal" & term_original == "Log(scale)" ~ exp(estimate_original),
          outc_distribution == "lognormal" & term_original == "(Intercept)" ~ estimate_original,
          TRUE ~ estimate_original
        ),
        # term is correct regardless of weibull or lognormal distribution
        term = dplyr::case_when(
          term_original == "Log(scale)" ~ "shape",
          term_original == "(Intercept)" ~ "scale",
          TRUE ~ paste0("beta", stringr::str_remove_all(string = term_original,
                                               pattern = "_|:"))
        )
      ) %>%
      dplyr::select(term, estimate) %>%
      replicate(
        n = k, .,
        simplify = FALSE
      ) %>%
      dplyr::bind_rows(.id = "subgroup") %>%
      dplyr::mutate(name = paste0(.data$term, "_", .data$subgroup)) %>%
      dplyr::select(name, estimate)

    # replicate k times for the 3 latent subgroups (below starting_values_df_list will replicate for the number of initial partitions)
    starting_values_df <- tidyr::pivot_wider(survreg_for_start_tidy,
        names_from = .data$name,
        values_from = .data$estimate
      )
  }

  if (starting_values_type %in% c("single_survreg", "uniform_pct")) {
    # browser()

    # don't want seed here or the starting values will be the same for each repetition!
    starting_values_df_list <- replicate(n = n_inits, expr = starting_values_df, simplify = FALSE) %>%
      purrr::map(., ~
        tidyr::pivot_longer(starting_values_df,
          cols = c(
            dplyr::starts_with("shape"),
            dplyr::starts_with("scale"),
            dplyr::starts_with("beta")
          )
        ) %>%
          # so that runif isn't just 0 to 0
          dplyr::mutate(value = dplyr::case_when(
            value == 0 ~ 0.0001,
            TRUE ~ value
          )) %>%
          dplyr::rowwise() %>%
          dplyr::mutate(
            min_for_runif = dplyr::case_when(
              # if weibull, shape and scale can't be negative
              outc_distribution == "weibull" & !grepl("shape|scale", name, ignore.case = TRUE) ~ min(
                value * (1 - starting_values_window),
                value * (1 + starting_values_window)
              ),
              outc_distribution == "weibull" & grepl("shape|scale", name, ignore.case = TRUE) ~ max(value * (1 - starting_values_window), 0),
              # if lognormal, only shape can't be negative
              outc_distribution == "lognormal" & grepl("scale", name, ignore.case = TRUE) & starting_scale_logn == "exp" ~ min(
                exp(value) * (1 - starting_values_window),
                exp(value) * (1 + starting_values_window)
              ),
              outc_distribution == "lognormal" & grepl("scale", name, ignore.case = TRUE) & starting_scale_logn == "log" ~ min(
                value * (1 - starting_values_window),
                value * (1 + starting_values_window)
              ),
              outc_distribution == "lognormal" & !grepl("shape|scale", name, ignore.case = TRUE) ~ min(
                value * (1 - starting_values_window),
                value * (1 + starting_values_window)
              ),
              outc_distribution == "lognormal" & grepl("shape", name, ignore.case = TRUE) ~ max(value * (1 - starting_values_window), 0)
            ),
            max_for_runif = dplyr::case_when(
              # if weibull, shape and scale can't be negative
              outc_distribution == "weibull" & !grepl("shape|scale", name, ignore.case = TRUE) ~ max(
                value * (1 + starting_values_window),
                value * (1 - starting_values_window)
              ),
              outc_distribution == "weibull" & grepl("shape|scale", name, ignore.case = TRUE) ~ value * (1 + starting_values_window),
              # if lognormal, only shape can't be negative)
              outc_distribution == "lognormal" & grepl("scale", name, ignore.case = TRUE) & starting_scale_logn == "exp" ~ max(
                exp(value) * (1 + starting_values_window),
                exp(value) * (1 - starting_values_window)
              ),
              outc_distribution == "lognormal" & grepl("scale", name, ignore.case = TRUE) & starting_scale_logn == "log" ~ max(
                value * (1 + starting_values_window),
                value * (1 - starting_values_window)
              ),
              outc_distribution == "lognormal" & !grepl("shape", name, ignore.case = TRUE) ~ max(
                value * (1 + starting_values_window),
                value * (1 - starting_values_window)
              ),
              outc_distribution == "lognormal" & grepl("shape", name, ignore.case = TRUE) ~ value * (1 + starting_values_window)
            ),
            hat1 = runif(n = 1, min_for_runif, max_for_runif),
            hat = ifelse(outc_distribution == "lognormal" & grepl("scale", name, ignore.case = TRUE) & starting_scale_logn == "exp",
                         log(hat1),
                         hat1),
            # name2a = str_replace_all(
            #   string = name,
            #   pattern = "beta_tx",
            #   replacement = "beta"
            # ),
            name2 = paste0(stringr::str_remove_all(
              string = name,
              pattern = "_"
            ), "_hat"),
          ) %>%
          dplyr::select(-value, -min_for_runif, -max_for_runif, -name, -hat1) %>%
          tidyr::pivot_wider(
            names_from = name2,
            values_from = hat
          ))
  }

  return(starting_values_df_list)
}
