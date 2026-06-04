#' EM Algorithm for Finite Mixture Models with Survival Endpoints
#'
#' @param input_df Input data frame containing 1 row/observation, along with
#'   each observation's event status (0=censored, 1=event) and event time
#'   variables
#' @param weights_input Variable for inverse probability of treatment weights,
#'   if applicable. For IPCW-FMM, supply the IPCW. For IPCW-FMM that also uses
#'   IPTW, the product of the ITPW and IPCW may be supplied.
#' @param starting_values_input_df Input dataset with starting values for algorithm
#' @param k Number of subgroups
#' @param outc_model_time Variable indicating the time to event or censoring for
#'   each observation.
#' @param outc_model_status Variable indicating the event status for the
#'   time-to-event outcome. It is assumed that 0 = censored, 1 = event.
#' @param outc_model_covars Names of covariates to include in the outcome models
#'   for each subgroup
#' @param outc_model_formula Formula object for subgroup-specific outcome models
#' @param outc_distribution Outcome distribution for subgroup-specific outcome
#'   models. Currently allowed values are "Weibull" and "Log-Normal" (not
#'   case-sensitive)
#' @param covariates_subgroup_model Vector of covariates to include in subgroup membership model
#' @param n_inits Number of initial partitions for the EM algorithm. Default is
#'   5. A higher number of initial partitions may result in greater stability of
#'   estimates.
#' @param tolerance Convergence criteria for the change in the log-likelihood
#'   for the EM algorithm.
#' @param conv_pct_criteria Convergence criteria for the percentage of
#'   observations changing subgroup. Specify -1 to only use the log-likelihood
#'   as the convergence criteria.
#' @param max_iter Maximum number of iterations. Default is 200.
#'
#' @returns List
fmm_em_algorithm <- function(input_df,
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
                             tolerance = 1e-3,
                             conv_pct_criteria = -1,
                             max_iter = 200) {

  # initialize empty vectors to store results
  log_likelihood_values <- c()
  convergence_iter <- 1
  convergence_status <- 0
  convergence_message <- NULL

  # null dataframes
  dfs_to_create <- c(
    "ests",
    "pi_hat_ests",
    paste0("final_outcome_model_tidy_", 1:k),
    paste0("final_outcome_model_cov_mtx", 1:k),
    "final_subgroup_model_tidy",
    "final_subgroup_model_cov_mtx",
    "subgroup_assn",
    "subgroup_assn_all",
    "n_subgroups_assigned"
  )

  empty_dfs <- purrr::map(dfs_to_create, ~ assign(.x, data.frame()))

  names(empty_dfs) <- dfs_to_create

  list2env(empty_dfs,
           envir = environment()
  )

  # browser()

  # for each iteration of algorithm
  for (i in 1:max_iter) {
    # print(paste0("iteration ", i))

    # initial partition -------------------------------------------------------
    # if first iteration, estimate pi-hat based on distribution in initial partition
    if (i == 1) {
      # browser()
      # define starting values for parameters (only at first iteration)
      # create objects based on [parameter]_hat variables in starting_values_df_list
      # starting values
      ests_long <- starting_values_input_df %>%
        dplyr::select(tidyselect::ends_with("_hat")) %>%
        tidyr::pivot_longer(
          cols = tidyselect::everything(),
          names_to = "term",
          values_to = "estimate"
        ) %>%
        dplyr::mutate(
          iter = 0,
          # k = stringr::str_extract(pattern = "[0-9]+", string = .data$term),
          # term = stringr::str_remove_all(pattern = "[0-9]+|", string = .data$term)
          k = stringr::str_remove_all(string = .data$term,
                                           pattern = paste0(paste0("beta",
                                                            str_remove_all(
                                                              string = outc_model_covars[order(nchar(outc_model_covars), decreasing = TRUE)],
                                                              pattern = "_|:|\\*"),
                                                            collapse = "|"),
                                           "|shape|scale|_hat")),
          term = paste0(stringr::str_extract(string = term,
                                       pattern = paste0("shape|scale|",
                                              paste0("beta",
                                                     str_remove_all(
                                                       string = outc_model_covars[order(nchar(outc_model_covars), decreasing = TRUE)],
                                                       pattern = "_|:|\\*"),
                                                     collapse = "|"))),
                                       "_hat")
        )

      # 1 row/k, 1 col/term
      ests_wide <- ests_long %>%
        tidyr::pivot_wider(
          id_cols = c(.data$iter, .data$k),
          names_from = .data$term,
          values_from = .data$estimate
        )

      # repetition 250 of S1257_2257_3257_TxIPTW2Null_LNul_E30_N500 has initial single survreg that says '2 not defined because of singularities'
      if (nrow(ests_wide %>%
               dplyr::filter(.data$iter == 0) %>%
               tidyr::drop_na()) == 0) {
        convergence_status <- 0
        convergence_iter <- i
        convergence_message <- paste0("Algorithm did not converge after ", i, " iterations (survreg error - singularities).")
        # message(convergence_message)
        break
      }

      # add variable for assn_subgroup, randomly assigned
      x <- input_df %>%
        dplyr::mutate(assn_subgroup = sample(seq(1, k),
                                      size = nrow(.),
                                      replace = TRUE
        ))

      # also replicate k times for processing through algorithm
      x_k <- replicate(x, n = k, simplify = FALSE) %>%
        dplyr::bind_rows(.id = "k") %>%
        dplyr::arrange(.data$record_id, .data$k)

      # browser()
      # now model P(assn subgroup)
      pi_hat_logistic <- (nnet::multinom(assn_subgroup ~ .,
                                         data = x %>%
                                           dplyr::select(
                                             assn_subgroup,
                                             tidyr::all_of(covariates_subgroup_model)
                                           ), trace = FALSE
      ))

      # 1 row/record_id
      prior_probability_wide <- stats::predict(pi_hat_logistic,
                                        newdata = x,
                                        type = "probs"
      ) %>%
        # add name_repair since name is different w/ <3 vs >=3 subgroups
        tibble::as_tibble(.name_repair = custom_name_repair) %>%
        dplyr::rename_with(~ paste0("priorprob_subgroup", .x))

      if (k == 2) {
        prior_probability_wide <- prior_probability_wide %>%
          dplyr::mutate(priorprob_subgroup1 = 1 - priorprob_subgroup2) %>%
          # order in df matters for tidyr::pivot_wider below
          dplyr::select(priorprob_subgroup1, priorprob_subgroup2)
      }

      # 1 row/record_id/subgroup
      prior_probability <- cbind(
        x_k %>%
          dplyr::select(record_id),
        prior_probability_wide %>%
          # get 1 row/subgroup
          tidyr::pivot_longer(
            cols = tidyr::starts_with("priorprob"),
            names_to = "k",
            names_prefix = "priorprob_subgroup",
            values_to = "prior_probability"
          )
      )

      # need to get current parameter estimates * covariate values for each
      # patient to supply in E-step calculation of posterior probability
      # need to remove intercept from formula
      # don't need to do this at each iteration
      mm <- stats::model.matrix(stats::update(outc_model_formula, . ~ . - 1),
                                data = input_df)

      # for first iteration use ests_wide
      # for subsequent use aft_output_tidy
      ests_coefs <- t(starting_values_input_df %>% dplyr::select(tidyr::starts_with("beta")))
      # ests_coefs <- as.matrix(ests_wide %>% dplyr::select(tidyr::starts_with("beta"))) # doesn't work w >1 covar

      # split coefficients by subgroup
      ests_coefs_by_k <- lapply(split(ests_coefs, rep(c(1:k), each=length(outc_model_covars))),
                                matrix,
                                nrow = length(outc_model_covars))

      # for each subgroup, get product of model matrix and current estimates
      # stacks each model matrix * estimates for each subgroup on top of each other
      mm_ests <- purrr::map(ests_coefs_by_k, ~mm %*% .x) %>%
        purrr::map_df(., as.data.frame, row.names = FALSE, .id = "k") %>%
        dplyr::rename(beta_variable = V1) %>%
        dplyr::group_by(.data$k) %>%
        # use input_df$record_id instead of recreating here with 1:n()
        dplyr::mutate(record_id = dplyr::pull(input_df, record_id)) %>%
        dplyr::ungroup()
    } # end iter 1

    # browser()

    # E-Step ------------------------------------------------------------------
    # E-Step: Compute the posterior probabilities
    # to incorporate covariates, use logistic regression to estimate posterior probability
    # get predicted probability of observed subgroup
    x1 <- dplyr::full_join(
      x_k %>%
        # remove prior probabilities from the previous iteration
        dplyr::select(-tidyr::starts_with("priorprob")),
      prior_probability,
      by = c("record_id", "k")
    ) %>%
      # merge on shape/scale/tx parameters
      dplyr::left_join(ests_wide %>%
                  dplyr::select(k, tidyr::contains("shape"), tidyr::contains("scale")),
                by = "k") %>%
      dplyr::left_join(mm_ests,
                by = c("k", "record_id")) %>%
      {
        if (outc_distribution == "weibull") {
          dplyr::mutate(.,
                 tau = dplyr::case_when(
                   # weibull distribution
                   outc_distribution == "weibull" & get(outc_model_status) %in% c(1) ~ prior_probability * pmax(
                     0.00000000000000000001,
                     stats::dweibull(get(outc_model_time),
                              shape = shape_hat,
                              # scale = scale_hat * exp(betatx_hat*tx)
                              scale = scale_hat * exp(beta_variable)
                     )
                   ),
                   outc_distribution == "weibull" & get(outc_model_status) == 0 ~ prior_probability * pmax(
                     0.00000000000000000001,
                     (1 - stats::pweibull(get(outc_model_time),
                                   shape = shape_hat,
                                   # scale = scale_hat * exp(betatx_hat*tx)
                                   scale = scale_hat * exp(beta_variable)
                     ))
                   )
                 )
          )
        } else if (outc_distribution == "lognormal") {
          dplyr::mutate(.,
                 tau = dplyr::case_when(
                   # lognormal distribution
                   outc_distribution == "lognormal" & get(outc_model_status) %in% c(1) ~ prior_probability * pmax(
                     0.00000000000000000001,
                     stats::dlnorm(get(outc_model_time),
                            sdlog = shape_hat,
                            meanlog = scale_hat + beta_variable
                            # meanlog = scale_hat + betatx_hat * tx +
                            #   betacovariatesimtmbzscore_hat * covariate_sim_tmb_zscore +
                            #   betatxcovariatesimtmbzscore_hat * tx * covariate_sim_tmb_zscore
                     )
                   ),
                   outc_distribution == "lognormal" & get(outc_model_status) == 0 ~ prior_probability * pmax(
                     0.00000000000000000001,
                     (1 - stats::plnorm(get(outc_model_time),
                                 sdlog = shape_hat,
                                 meanlog = scale_hat + beta_variable
                                 # meanlog = scale_hat + betatx_hat * tx +
                                 #   betacovariatesimtmbzscore_hat * covariate_sim_tmb_zscore +
                                 #   betatxcovariatesimtmbzscore_hat * tx * covariate_sim_tmb_zscore
                     ))
                   )
                 )
          )
        }
      } %>%
      dplyr::group_by(.data$record_id) %>%
      dplyr::mutate(
        sum_tau = sum(.data$tau),
        # posterior probabilities
        posterior_prob = .data$tau / .data$sum_tau,
        # assign subgroup based on highest posterior probability
        assn_subgroup = which.max(.data$posterior_prob),
        # multiply weights by IPTW for weighting in AFT model
        posterior_prob_weights_input = get(weights_input) * .data$posterior_prob,
      ) %>%
      dplyr::ungroup()

    # check convergence ------------------------------------------------------------------
    # final models ------------------------------------------------------------
    # if the algorithm converged, pull the final model info
    if (convergence_status == 1) {
      # browser()
      # final outcome models
      # tidy final outcome models
      final_outcome_model_tidy <- purrr::map(est_survreg, broom::tidy, conf.int = TRUE)
      names(final_outcome_model_tidy) <- paste0("final_outcome_model_tidy_", names(final_outcome_model_tidy))

      # need matrices for joint test
      final_outcome_model_coefs <- purrr::map(est_survreg, broom::tidy, conf.int = TRUE) %>%
        purrr::map(., purrr::pluck, "estimate") %>%
        purrr::map(., ~ matrix(.x))

      names(final_outcome_model_coefs) <- paste0("final_outcome_model_coefs_", names(final_outcome_model_coefs))

      # also need covariance matrices
      final_outcome_model_cov_mtx <- purrr::map(est_survreg, stats::vcov) %>%
        purrr::map(., tibble::as_tibble) %>%
        purrr::map(., janitor::clean_names)

      names(final_outcome_model_cov_mtx) <- paste0("final_outcome_model_cov_mtx_", names(final_outcome_model_cov_mtx))

      # outcome model objects (has to go below other objects)
      names(est_survreg) <- paste0("final_outcome_model_", names(est_survreg))

      # export as individual objects to environment
      list2env(c(est_survreg, final_outcome_model_tidy, final_outcome_model_coefs,
                 final_outcome_model_cov_mtx),
               envir = environment()
      )

      # final subgroup model
      final_subgroup_model_tidy <- broom::tidy(pi_hat_logistic2, conf.int = TRUE)

      # need matrices for joint test
      # DO NOT USE THIS ROW, re-orders as int1, int2, tmb1, tmb2, instead of int1, tmb1, etc.
      # final_subgroup_model_coefs <- matrix(coef(pi_hat_logistic2))
      final_subgroup_model_coefs <- matrix(broom::tidy(pi_hat_logistic2)$estimate)

      # also need covariance matrix
      final_subgroup_model_cov_mtx <- unname(stats::vcov(pi_hat_logistic2))

      final_subgroup_model_cov_mtx_tidy <- tibble::as_tibble(stats::vcov(pi_hat_logistic2)) %>%
        janitor::clean_names()

      # final subgroup assignment from posterior probability
      subgroup_assn <- x1 %>%
        tidyr::pivot_wider(
          id_cols = c(.data$record_id, tidyr::contains("latent_subgroup"), .data$assn_subgroup),
          names_from = k,
          names_prefix = "posterior_prob",
          values_from = .data$posterior_prob
        )

      break
    }

    # M-Step ------------------------------------------------------------------
    # browser()

    ## re-estimate prior probabilities -----------------------------------------
    # number of subgroups assigned
    n_subgroups_assigned_i <- x1 %>%
      dplyr::mutate(iter = i) %>%
      dplyr::distinct(.data$iter, .data$assn_subgroup) %>%
      dplyr::count(.data$iter, name = "n_subgroups")

    n_subgroups_assigned <- dplyr::bind_rows(
      n_subgroups_assigned,
      n_subgroups_assigned_i
    )

    # check number of subgroups
    if (n_subgroups_assigned_i$n_subgroups == 1) {
      convergence_status <- 0
      convergence_message <- paste0("Algorithm did not converge (patients assigned to only 1 subgroup).")
      convergence_iter <- i
      break
    }

    # M-Step: Update the parameter estimates
    # logistic regression with outcome = cluster assignment, covariate
    # this will be individual specific, depending on their covariates
    # need to set up data in particular structure
    # browser()

    pi_hat_logistic2 <- (nnet::multinom(k ~ . - posterior_prob,
                                        data = x1 %>%
                                          dplyr::select(
                                            k,
                                            tidyr::all_of(covariates_subgroup_model),
                                            posterior_prob
                                          ),
                                        weights = posterior_prob,
                                        trace = FALSE
    ))

    # 1 rec/record ID and subgroup w/ a var for K
    # also have to add var for k in prior_probability from starting values
    # 1 row/record_id
    # when using syntax of . - posterior_prob above,
    # for some reason predict requires the posterior prob
    # variable on the newdata object
    # set to 1
    # (checked by comparing predict using old syntax where covariates are explicitly listed and using predict w/o posterior_prob = 1)
    prior_probability_wide <- stats::predict(pi_hat_logistic2,
                                      newdata = x %>%
                                        dplyr::mutate(posterior_prob = 1),
                                      type = "probs"
    ) %>%
      tibble::as_tibble(.name_repair = custom_name_repair) %>%
      dplyr::rename_with(~ paste0("priorprob_subgroup", .x))

    # browser()

    if (k == 2) {
      prior_probability_wide <- prior_probability_wide %>%
        dplyr::mutate(priorprob_subgroup1 = 1 - .data$priorprob_subgroup2) %>%
        # order in df matters for tidyr::pivot_wider below
        dplyr::select(priorprob_subgroup1, priorprob_subgroup2)
    }

    # 1 row/record_id/subgroup
    prior_probability <- cbind(
      x1 %>%
        dplyr::select(record_id),
      prior_probability_wide %>%
        # get 1 row/subgroup
        tidyr::pivot_longer(
          cols = tidyr::starts_with("priorprob"),
          names_to = "k",
          names_prefix = "priorprob_subgroup",
          values_to = "prior_probability"
        )
    )

    ## estimate outcome model -----------------------------------------
    # to add left truncation
    # to incorporate covariates
    # reset at each repetition to pick up when it errors out due to 0 weights
    est_survreg <- NULL

    # split the data by k subgroups
    # this is keeping 1 full copy of the dataset with corresponding posterior probability values
    est_survreg <- x1 %>%
      split(.$k) %>%
      purrr::map(., ~ tryCatch(
        {
          survival::survreg(outc_model_formula,
                            data = .x,
                            weights = posterior_prob_weights_input,
                            robust = TRUE,
                            dist = outc_distribution
          )
        }, # end tryCatch expression
        # error = function(e) NA,
        # warning = function(w) NA
        error = function(e) {
          message("Caught an error!")
          print(e)
        },
        warning = function(w) {
          message(paste0("Caught a warning in iter ", i))
          print(w)
        }
      ))

    # browser()

    # if weights become 0, survreg doesn't run and the number of rows in the model becomes 0
    # want to exit that repetition
    # error message is: Error in survreg.fit(X, Y, weights, offset, init =
    # init, controlvals = control,  : Invalid weights, must be >0
    if (nrow(purrr::map_chr(est_survreg,
      ~ class(.x)[1]
    ) %>%
    tibble::as_tibble() %>%
    dplyr::filter(.data$value != "survreg")) > 0) {
      convergence_status <- 0
      convergence_iter <- i
      convergence_message <- paste0("Algorithm did not converge after ", i, " iterations (survreg error).")
      # message(convergence_message)
      break
    }

    # tidy AFT model output
    aft_output_tidy <- purrr::map_df(est_survreg,
                              ~ broom::tidy(.x, conf.int = TRUE) %>%
                                # with flexsurvreg, shape and scale parameters are returned in a parameterization consistent with rweibull
                                # with survreg, the log(scale) term returned is 1/rweibull shape
                                # with survreg, the exp(survreg intercept) = rweib scale
                                # this is documented in the survreg example
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
                                  term = dplyr::case_when(
                                    term_original == "Log(scale)" ~ "shape",
                                    term_original == "(Intercept)" ~ "scale",
                                    TRUE ~ term_original
                                  ),
                                  term_for_loglik = dplyr::case_when(
                                    # term == "tx" ~ "beta",
                                    term == "tx:covariate_sim_tmb_zscore" ~ "betatxcovariatesimtmbzscore",
                                    grepl("tx|tmb", term) ~ paste0("beta", stringr::str_remove_all(
                                      string = term,
                                      pattern = "_"
                                    )),
                                    TRUE ~ term
                                  )
                                ),
                              .id = "k"
    )

    ests_coefs <- aft_output_tidy %>%
        dplyr::filter(!grepl("shape|scale", .data$term, ignore.case = TRUE)) %>%
        dplyr::select(estimate) %>%
        as.matrix()

    # split coefficients by subgroup
    ests_coefs_by_k <- lapply(split(ests_coefs, rep(c(1:k), each=length(outc_model_covars))),
                              matrix,
                              nrow = length(outc_model_covars))

    # for each subgroup, get product of model matrix and current estimates
    # stacks each model matrix * estimates for each subgroup on top of each other
    mm_ests <- purrr::map(ests_coefs_by_k, ~mm %*% .x) %>%
      purrr::map_df(., as.data.frame, row.names = FALSE, .id = "k") %>%
      dplyr::rename(beta_variable = V1) %>%
      dplyr::group_by(.data$k) %>%
      # use input_df$record_id instead of recreating here with 1:n()
      dplyr::mutate(record_id = dplyr::pull(input_df, record_id)) %>%
      dplyr::ungroup()

    # extract log-likelihood
    aft_loglik <- purrr::map(est_survreg, "loglik") %>%
      purrr::map(., purrr::pluck, 2) %>%
      purrr::map_df(., tibble::as_tibble, .id = "k") %>%
      dplyr::rename(estimate = .data$value)

    # combine estimates
    ests <- dplyr::bind_rows(
      ests,
      aft_output_tidy %>% dplyr::mutate(
        model = "Outcome Model",
        assn_subgroup = as.numeric(.data$k),
        iter = i
      ),
      # for log likelihood, want model$loglik[2],
      # which corresponds to the loglik(model) object
      aft_loglik %>%
        dplyr::mutate(
          term = "loglik",
          model = "Outcome Model",
          assn_subgroup = as.numeric(.data$k),
          iter = i
        )
    )

    # browser()

    # set up with 1 col/estimate and 1 row/subgroup to merge onto dataset
    # at next iteration of E-step
    ests_wide <- aft_output_tidy %>%
        dplyr::mutate(iter = i) %>%
        tidyr::pivot_wider(
          id_cols = c(.data$k, .data$iter),
          names_from = term_for_loglik,
          names_glue = "{term_for_loglik}_hat",
          values_from = .data$estimate
        )

    ests <- dplyr::bind_rows(
      ests,
      # coefficient for the covariate in the model of P(mixture 1)
      broom::tidy(pi_hat_logistic2, conf.int = TRUE) %>%
        dplyr::mutate(
          model = "Latent Subgroup Model",
          iter = i,
          assn_subgroup = as.numeric(.data$y.level)
        ) %>%
        dplyr::select(-y.level)
    )

    # pi is now a vector, need to pull for each subject
    # only save for last iteration
    pi_hat_ests <- dplyr::bind_rows(
      pi_hat_ests, # uncomment if we want across iterations
      cbind(
        x %>%
          dplyr::select(-tidyr::starts_with("priorprob_")),
        prior_probability_wide
      ) %>%
        dplyr::mutate( # model = "Latent Subgroup Model",
          term = "pi_hat",
          iter = i
        ) %>%
        dplyr::select(
          term, iter,
          # will only have latent_subgroup variable if running on the simulated data
          record_id, tidyr::any_of("latent_subgroup"), assn_subgroup,
          tidyr::starts_with("latent_vprob"),
          tidyr::starts_with("priorprob_"),
          tidyr::starts_with("posterior")
        )
    )

    # repeat but w/o reference to latent_subgroup
    subgroup_assn_all <- dplyr::bind_rows(
      subgroup_assn_all,
      x1 %>%
        dplyr::distinct(.data$record_id, .data$assn_subgroup) %>%
        dplyr::mutate(iter = i) %>%
        dplyr::select(iter, record_id, assn_subgroup)
    )

    # log likelihood ----------------------------------------------------------
    # Compute log-likelihood
    # browser()
    # df with prior probabilities from E-step (1 row/record ID/k)
    # merge on original df (1 row/record ID)
    # then merge on AFT outcome model estimates (1 row/term/k)
    for_loglik <- dplyr::left_join(
      dplyr::left_join(x1 %>%
                         dplyr::select(record_id, k, prior_probability),
                       mm_ests,
                       by = c("k", "record_id")),
      x %>%
        dplyr::select(record_id,
               time = outc_model_time, status = outc_model_status,
               # tx, covariate_sim_tmb_zscore,
               tidyr::all_of(weights_input)
        ),
      by = "record_id"
    ) %>%
      dplyr::left_join(
        aft_output_tidy %>%
          tidyr::pivot_wider(
            id_cols = k,
            names_from = term_for_loglik,
            values_from = .data$estimate
          ) %>%
          dplyr::select(k, shape, scale),
        by = "k"
      )

    log_likelihood_df <- log_likelihood_with_censoring(
      input_df = for_loglik,
      distribution = outc_distribution
    )

    log_likelihood_values <- c(
      log_likelihood_values,
      log_likelihood_df$log_likelihood_values
    )

    # browser()
    # check convergence ----------------------------------------------------------
    # if beyond the 1st iteration, check for convergence
    # if (i > 1) {
    # browser()
    ## computations relevant to checking convergence
    # want these outside of the check convergence function so they update at each iteration
    # compute % change subgroup assignment
    subgroup_assn_all_chg <- subgroup_assn_all %>%
      dplyr::filter(.data$iter %in% c(i, i - 1)) %>%
      dplyr::group_by(.data$record_id) %>%
      dplyr::arrange(.data$record_id, .data$iter) %>%
      dplyr::mutate(chg_assn_subgroup = .data$assn_subgroup != dplyr::lag(.data$assn_subgroup)) %>%
      tidyr::drop_na() %>%
      dplyr::ungroup() %>%
      dplyr::summarize(pct_chg_assn_subgroup = mean(.data$chg_assn_subgroup)) %>%
      dplyr::pull(.data$pct_chg_assn_subgroup)


    ## non-convergence
    # if log-likelihood = NA ()
    if (is.nan(log_likelihood_values[i]) | is.na(log_likelihood_values[i])) {
      convergence_status <- 0
      convergence_message <- paste0("Algorithm did not converge (log likelihood NA)")
      convergence_iter <- i
      # message(convergence_message)
      # break
    }

    if (n_subgroups_assigned_i$n_subgroups == 1) {
      convergence_status <- 0
      convergence_message <- paste0("Algorithm did not converge (patients assigned to only 1 subgroup).")
      convergence_iter <- i
      # break
    }

    # browser()
    ## convergence
    # % records changing assigned subgroup
    # browser()

    if (i > 1 && subgroup_assn_all_chg <= conv_pct_criteria) {
      convergence_status <- 1
      convergence_message <- paste0("Algorithm converged (<=", 100 * conv_pct_criteria, "% records changed assigned subgroup).")
      convergence_iter <- i
      # break
    }

    # algorithm converged based on change in log-likelhiood
    if (i > 1 && (100 * abs((log_likelihood_values[i] - log_likelihood_values[i - 1]) / log_likelihood_values[i - 1]) < tolerance)) {
      # if (i > 1 && abs(log_likelihood_values[i] - log_likelihood_values[i - 1]) < tolerance) {
      convergence_status <- 1
      convergence_message <- paste0("Algorithm converged (percent change in log likelihood).")
      convergence_iter <- i
      # message(convergence_message)
      # break
    }

    # algorithm went to maximum number of specified iterations without converging
    if (i > 1 && i == max_iter && ((100 * abs((log_likelihood_values[i] - log_likelihood_values[i - 1]) / log_likelihood_values[i - 1]) >= tolerance))) {
      # if (i > 1 && i == max_iter && abs(log_likelihood_values[i] - log_likelihood_values[i - 1]) >= tolerance){
      convergence_status <- 0
      convergence_message <- paste0("Algorithm did not converge (maxiter reached)")
      convergence_iter <- i
      # message(convergence_message)
      # break
    }
    # browser()
  } # end loop over max_iter

  # return objects ----------------------------------------------------------
  return(c(
    purrr::list_flatten(list(
      "starting_values" = dplyr::bind_rows(starting_values_input_df) %>%
        dplyr::select(tidyr::contains("hat")),
      # outcome models
      mget(ls(pattern = "^final_outcome_model_\\d$")),
      mget(ls(pattern = "final_outcome_model_tidy_")),
      mget(ls(pattern = "final_outcome_model_cov_mtx_")),
      # subgroup membership model
      "final_subgroup_model" = pi_hat_logistic2,
      "final_subgroup_model_tidy" = final_subgroup_model_tidy,
      # "final_subgroup_model_coefs" = final_subgroup_model_coefs,
      "final_subgroup_model_cov_mtx" = final_subgroup_model_cov_mtx,
      # "final_subgroup_model_cov_mtx_tidy" = final_subgroup_model_cov_mtx_tidy,
      "subgroup_assn" = subgroup_assn,
      "pi_hat_ests" = pi_hat_ests,
      # return final dataset
      "final_df" = x1,
      # log likelihood values and convergence
      "log_likelihood_values" = log_likelihood_values,
      "convergence_status" = convergence_status,
      "convergence_message" = convergence_message,
      "convergence_iter" = convergence_iter,
      "ests" = ests # estimates across iterations of algorithm
    )
  )))
} # end em_function
