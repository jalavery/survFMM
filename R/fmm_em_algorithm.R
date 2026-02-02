#' EM Algorithm for Finite Mixture Models with Survival Endpoints
#'
#' @param starting_values_input_df Input dataset with starting values for algorithm
#' @param input_df Input data frame containing 1 row/observation, along with
#'   each observation's event status (0=censored, 1=event) and event time
#'   variables
#' @param weights_input Variable for inverse probability of treatment weights,
#'   if applicable. For IPCW-FMM, supply the IPCW. For IPCW-FMM that also uses
#'   IPTW, the product of the ITPW and IPCW may be supplied.
#' @param k Number of subgroups
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
#'
#' @examples
fmm_em_algorithm <- function(input_df,
                             starting_values_input_df,
                             k,
                             outc_model_formula,
                             weights_input,
                             outc_distribution,
                             outc_model_covars,
                             covariates_subgroup_model,
                             n_inits = 5,
                             tolerance = 1e-3,
                             conv_pct_criteria = -1,
                             max_iter = 200) {

  # initialize empty vectors to store results
  log_likelihood_values <- c()
  convergence_message <- NULL

  # null dataframes
  dfs_to_create <- c(
    "ests",
    "pi_hat_ests",
    "final_outcome_model_tidy_1",
    "final_outcome_model_tidy_2",
    "final_outcome_model_tidy_3",
    "final_outcome_model_cov_mtx_1",
    "final_outcome_model_cov_mtx_2",
    "final_outcome_model_cov_mtx_3",
    "final_subgroup_model_tidy",
    "final_subgroup_model_cov_mtx",
    "subgroup_assn",
    "subgroup_assn_all",
    "n_subgroups_assigned"
  )

  empty_dfs <- map(dfs_to_create, ~ assign(.x, data.table::data.table()))

  names(empty_dfs) <- dfs_to_create

  list2env(empty_dfs,
           envir = environment()
  )

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
        select(ends_with("_hat")) %>%
        pivot_longer(
          cols = everything(),
          names_to = "term",
          values_to = "estimate"
        ) %>%
        mutate(
          iter = 0,
          k = str_extract(pattern = "[0-9]+", string = term),
          term = str_remove_all(pattern = "[0-9]+|", string = term)
        )

      # 1 row/k, 1 col/term
      ests_wide <- ests_long %>%
        pivot_wider(
          id_cols = c(iter, k),
          names_from = term,
          values_from = estimate
        )

      # repetition 250 of S1257_2257_3257_TxIPTW2Null_LNul_E30_N500 has initial single survreg that says '2 not defined because of singularities'
      if (nrow(ests_wide %>%
               filter(iter == 0) %>%
               drop_na()) == 0) {
        convergence_status <- 0
        convergence_iter <- i
        convergence_message <- paste0("Algorithm did not converge after ", i, " iterations (survreg error - singularities).")
        # message(convergence_message)
        break
      }

      # add variable for assn_subgroup, randomly assigned
      x <- input_df %>%
        mutate(assn_subgroup = sample(seq(1, k),
                                      size = nrow(.),
                                      replace = TRUE
        ))

      # also replicate k times for processing through algorithm
      x_k <- replicate(x, n = k, simplify = FALSE) %>%
        bind_rows(.id = "k") %>%
        arrange(record_id, k)


      # browser()
      # now model P(assn subgroup)
      pi_hat_logistic <- (nnet::multinom(assn_subgroup ~ .,
                                         data = x %>%
                                           select(
                                             assn_subgroup,
                                             all_of(covariates_subgroup_model)
                                           ), trace = FALSE
      ))

      # 1 row/record_id
      prior_probability_wide <- predict(pi_hat_logistic,
                                        newdata = x,
                                        type = "probs"
      ) %>%
        # add name_repair since name is different w/ <3 vs >=3 subgroups
        as_tibble(.name_repair = survFMM:::custom_name_repair) %>%
        rename_with(~ paste0("priorprob_subgroup", .x))

      if (k == 2) {
        prior_probability_wide <- prior_probability_wide %>%
          mutate(priorprob_subgroup1 = 1 - priorprob_subgroup2) %>%
          # order in df matters for pivot_wider below
          select(priorprob_subgroup1, priorprob_subgroup2)
      }

      # 1 row/record_id/subgroup
      prior_probability <- cbind(
        x_k %>%
          select(record_id),
        prior_probability_wide %>%
          # get 1 row/subgroup
          pivot_longer(
            cols = starts_with("priorprob"),
            names_to = "k",
            names_prefix = "priorprob_subgroup",
            values_to = "prior_probability"
          )
      )

      # need to get current parameter estimates * covariate values for each
      # patient to supply in E-step calculation of posterior probability
      # need to remove intercept from formula
      # don't need to do this at each iteration
      mm <- model.matrix(update(outc_model_formula, . ~ . - 1), data = input_df)

      # for first iteration use ests_wide
      # for subsequent use aft_output_tidy
      ests_coefs <- t(starting_values_input_df %>% select(starts_with("beta")))
      # ests_coefs <- as.matrix(ests_wide %>% select(starts_with("beta"))) # doesn't work w >1 covar
    } # end iter 1

    # E-Step ------------------------------------------------------------------
    # browser()


    if (i>1){
      ests_coefs <- aft_output_tidy %>%
      filter(!grepl("shape|scale", term, ignore.case = TRUE)) %>%
      select(estimate) %>%
      as.matrix()
    }

    # split coefficients by subgroup
    ests_coefs_by_k <- lapply(split(ests_coefs, rep(c(1:k), each=length(outc_model_covars))),
                              matrix,
                              nrow = length(outc_model_covars))

    # for each subgroup, get product of model matrix and current estimates
    # stacks each model matrix * estimates for each subgroup on top of each other
    mm_ests <- map(ests_coefs_by_k, ~mm %*% .x) %>%
      map_df(., as.data.frame, row.names = FALSE, .id = "k") %>%
      rename(beta_variable = V1) %>%
      group_by(k) %>%
      mutate(record_id = 1:n()) %>%
      ungroup()

    # browser()

    # E-Step: Compute the posterior probabilities
    # to incorporate covariates, use logistic regression to estimate posterior probability
    # get predicted probability of observed subgroup
    x1 <- full_join(
      x_k %>%
        # remove prior probabilities from the previous iteration
        select(-starts_with("priorprob")),
      prior_probability,
      by = c("record_id", "k")
    ) %>%
      # merge on shape/scale/tx parameters
      left_join(ests_wide %>%
                  select(k, contains("shape"), contains("scale")),
                by = "k") %>%
      left_join(mm_ests,
                by = c("k", "record_id")) %>%
      {
        if (outc_distribution == "weibull") {
          mutate(.,
                 tau = case_when(
                   # weibull distribution
                   outc_distribution == "weibull" & pfs_status %in% c(1) ~ prior_probability * pmax(
                     0.00001,
                     dweibull(tt_pfs_days,
                              shape = shape_hat,
                              # scale = scale_hat * exp(betatx_hat*tx)
                              scale = scale_hat * exp(beta_variable)
                     )
                   ),
                   outc_distribution == "weibull" & pfs_status == 0 ~ prior_probability * pmax(
                     0.00001,
                     (1 - pweibull(tt_pfs_days,
                                   shape = shape_hat,
                                   # scale = scale_hat * exp(betatx_hat*tx)
                                   scale = scale_hat * exp(beta_variable)
                     ))
                   )
                 )
          )
        } else if (outc_distribution == "lognormal") {
          mutate(.,
                 tau = case_when(
                   # lognormal distribution
                   outc_distribution == "lognormal" & pfs_status %in% c(1) ~ prior_probability * pmax(
                     0.00001,
                     dlnorm(tt_pfs_days,
                            sdlog = shape_hat,
                            meanlog = scale_hat + beta_variable
                            # meanlog = scale_hat + betatx_hat * tx +
                            #   betacovariatesimtmbzscore_hat * covariate_sim_tmb_zscore +
                            #   betatxcovariatesimtmbzscore_hat * tx * covariate_sim_tmb_zscore
                     )
                   ),
                   outc_distribution == "lognormal" & pfs_status == 0 ~ prior_probability * pmax(
                     0.00001,
                     (1 - plnorm(tt_pfs_days,
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
      group_by(record_id) %>%
      mutate(
        sum_tau = sum(tau),
        # posterior probabilities
        posterior_prob = tau / sum_tau,
        # assign subgroup based on highest posterior probability
        assn_subgroup = which.max(posterior_prob),
        # multiply weights by IPTW for weighting in AFT model
        posterior_prob_weights_input = get(weights_input) * posterior_prob,
      ) %>%
      ungroup()

    # M-Step ------------------------------------------------------------------
    # browser()

    ## re-estimate prior probabilities -----------------------------------------
    # number of subgroups assigned
    n_subgroups_assigned_i <- x1 %>%
      mutate(iter = i) %>%
      distinct(iter, assn_subgroup) %>%
      count(iter, name = "n_subgroups")

    n_subgroups_assigned <- bind_rows(
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
                                          select(
                                            k,
                                            all_of(covariates_subgroup_model),
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
    prior_probability_wide <- predict(pi_hat_logistic2,
                                      newdata = x %>%
                                        mutate(posterior_prob = 1),
                                      type = "probs"
    ) %>%
      as_tibble(.name_repair = survFMM:::custom_name_repair) %>%
      rename_with(~ paste0("priorprob_subgroup", .x))

    if (k == 2) {
      prior_probability_wide <- prior_probability_wide %>%
        mutate(priorprob_subgroup1 = 1 - priorprob_subgroup2) %>%
        # order in df matters for pivot_wider below
        select(priorprob_subgroup1, priorprob_subgroup2)
    }

    # 1 row/record_id/subgroup
    prior_probability <- cbind(
      x1 %>%
        select(record_id),
      prior_probability_wide %>%
        # get 1 row/subgroup
        pivot_longer(
          cols = starts_with("priorprob"),
          names_to = "k",
          names_prefix = "priorprob_subgroup",
          values_to = "prior_probability"
        )
    )

    # browser()

    ## estimate outcome model -----------------------------------------
    # to add left truncation
    # to incorporate covariates
    # reset at each repetition to pick up when it errors out due to 0 weights
    est_survreg <- NULL

    # split the data by k subgroups
    # this is keeping 1 full copy of the dataset with corresponding posterior probability values
    est_survreg <- x1 %>%
      split(.$k) %>%
      map(., ~ tryCatch(
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
    if (nrow(map_chr(est_survreg,
      ~ class(.x)[1]
    ) %>%
    as_tibble() %>%
    filter(value != "survreg")) > 0) {
      convergence_status <- 0
      convergence_iter <- i
      convergence_message <- paste0("Algorithm did not converge after ", i, " iterations (survreg error).")
      # message(convergence_message)
      break
    }

    # tidy AFT model output
    aft_output_tidy <- map_df(est_survreg,
                              ~ broom::tidy(.x, conf.int = TRUE) %>%
                                # with flexsurvreg, shape and scale parameters are returned in a parameterization consistent with rweibull
                                # with survreg, the log(scale) term returned is 1/rweibull shape
                                # with survreg, the exp(survreg intercept) = rweib scale
                                # this is documented in the survreg example
                                rename(
                                  estimate_original = estimate,
                                  term_original = term
                                ) %>%
                                mutate(
                                  estimate = case_when(
                                    outc_distribution == "weibull" & term_original == "Log(scale)" ~ 1 / exp(estimate_original),
                                    outc_distribution == "weibull" & term_original == "(Intercept)" ~ exp(estimate_original),
                                    outc_distribution == "lognormal" & term_original == "Log(scale)" ~ exp(estimate_original),
                                    outc_distribution == "lognormal" & term_original == "(Intercept)" ~ estimate_original,
                                    TRUE ~ estimate_original
                                  ),
                                  term = case_when(
                                    term_original == "Log(scale)" ~ "shape",
                                    term_original == "(Intercept)" ~ "scale",
                                    TRUE ~ term_original
                                  ),
                                  term_for_loglik = case_when(
                                    # term == "tx" ~ "beta",
                                    term == "tx:covariate_sim_tmb_zscore" ~ "betatxcovariatesimtmbzscore",
                                    grepl("tx|tmb", term) ~ paste0("beta", str_remove_all(
                                      string = term,
                                      pattern = "_"
                                    )),
                                    TRUE ~ term
                                  )
                                ),
                              .id = "k"
    )

    # extract log-likelihood
    aft_loglik <- map(est_survreg, "loglik") %>%
      map(., pluck, 2) %>%
      map_df(., as_tibble, .id = "k") %>%
      rename(estimate = value)

    # combine estimates
    ests <- bind_rows(
      # ests,
      aft_output_tidy %>% mutate(
        model = "Outcome Model",
        assn_subgroup = as.numeric(k),
        iter = i
      ),
      # for log likelihood, want model$loglik[2],
      # which corresponds to the loglik(model) object
      aft_loglik %>%
        mutate(
          term = "loglik",
          model = "Outcome Model",
          assn_subgroup = as.numeric(k),
          iter = i
        )
    )

    # browser()

    # set up with 1 col/estimate and 1 row/subgroup to merge onto dataset
    # at next iteration of E-step
    ests_wide <- aft_output_tidy %>%
        mutate(iter = i) %>%
        pivot_wider(
          id_cols = c(k, iter),
          names_from = term_for_loglik,
          names_glue = "{term_for_loglik}_hat",
          values_from = estimate
        )

    ests <- bind_rows(
      ests,
      # coefficient for the covariate in the model of P(mixture 1)
      broom::tidy(pi_hat_logistic2, conf.int = TRUE) %>%
        mutate(
          model = "Latent Subgroup Model",
          iter = i,
          assn_subgroup = as.numeric(y.level)
        ) %>%
        select(-y.level)
    )

    # pi is now a vector, need to pull for each subject
    # only save for last iteration
    pi_hat_ests <- bind_rows(
      # pi_hat_ests, # uncomment if we want across iterations
      cbind(
        x %>%
          select(-starts_with("priorprob_")),
        prior_probability_wide
      ) %>%
        mutate( # model = "Latent Subgroup Model",
          term = "pi_hat",
          iter = i
        ) %>%
        select(
          term, iter,
          # will only have latent_subgroup variable if running on the simulated data
          record_id, any_of("latent_subgroup"), assn_subgroup,
          starts_with("latent_vprob"),
          starts_with("priorprob_"),
          starts_with("posterior")
        )
    )

    # repeat but w/o reference to latent_subgroup
    subgroup_assn_all <- bind_rows(
      subgroup_assn_all,
      x1 %>%
        distinct(record_id, assn_subgroup) %>%
        mutate(iter = i) %>%
        select(iter, record_id, assn_subgroup)
    )

    # log likelihood ----------------------------------------------------------
    # Compute log-likelihood
    # browser()
    # df with prior probabilities from E-step (1 row/record ID/k)
    # merge on original df (1 row/record ID)
    # then merge on AFT outcome model estimates (1 row/term/k)
    for_loglik <- left_join(
      x1 %>%
        select(record_id, k, prior_probability, beta_variable),
      x %>%
        select(record_id,
               time = tt_pfs_days, status = pfs_status,
               # tx, covariate_sim_tmb_zscore,
               all_of(weights_input)
        ),
      by = "record_id"
    ) %>%
      left_join(
        aft_output_tidy %>%
          pivot_wider(
            id_cols = k,
            names_from = term_for_loglik,
            values_from = estimate
          ),
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
      filter(iter %in% c(i, i - 1)) %>%
      group_by(record_id) %>%
      arrange(record_id, iter) %>%
      mutate(chg_assn_subgroup = assn_subgroup != lag(assn_subgroup)) %>%
      drop_na() %>%
      ungroup() %>%
      summarize(pct_chg_assn_subgroup = mean(chg_assn_subgroup)) %>%
      pull(pct_chg_assn_subgroup)


    ## non-convergence
    # if log-likelihood = NA ()
    if (is.nan(log_likelihood_values[i]) | is.na(log_likelihood_values[i])) {
      convergence_status <- 0
      convergence_message <- paste0("Algorithm did not converge (log likelihood NA)")
      convergence_iter <- i
      # message(convergence_message)
      break
    }

    if (n_subgroups_assigned_i$n_subgroups == 1) {
      convergence_status <- 0
      convergence_message <- paste0("Algorithm did not converge (patients assigned to only 1 subgroup).")
      convergence_iter <- i
      break
    }

    # browser()
    ## convergence
    # % records changing assigned subgroup
    # browser()

    if (i > 1 && subgroup_assn_all_chg <= conv_pct_criteria) {
      convergence_status <- 1
      convergence_message <- paste0("Algorithm converged (<=", 100 * conv_pct_criteria, "% records changed assigned subgroup).")
      convergence_iter <- i
      break
    }

    # algorithm converged based on change in log-likelhiood
    if (i > 1 && (100 * abs((log_likelihood_values[i] - log_likelihood_values[i - 1]) / log_likelihood_values[i - 1]) < tolerance)) {
      # if (i > 1 && abs(log_likelihood_values[i] - log_likelihood_values[i - 1]) < tolerance) {
      convergence_status <- 1
      convergence_message <- paste0("Algorithm converged (percent change in log likelihood).")
      convergence_iter <- i
      # message(convergence_message)
      break
    }

    # algorithm went to maximum number of specified iterations without converging
    if (i > 1 && i == max_iter && ((100 * abs((log_likelihood_values[i] - log_likelihood_values[i - 1]) / log_likelihood_values[i - 1]) >= tolerance))) {
      # if (i > 1 && i == max_iter && abs(log_likelihood_values[i] - log_likelihood_values[i - 1]) >= tolerance){
      convergence_status <- 0
      convergence_message <- paste0("Algorithm did not converge (maxiter reached)")
      convergence_iter <- i
      # message(convergence_message)
      break
    }
    # browser()
  } # end loop over max_iter
  #
  # browser()

  # final models ------------------------------------------------------------
  # if the algorithm converged, pull the final model info
  if (convergence_status == 1) {
    # final outcome models
    final_outcome_model_tidy <- map(est_survreg, broom::tidy, conf.int = TRUE)
    names(final_outcome_model_tidy) <- paste0("final_outcome_model_tidy_", names(final_outcome_model_tidy))

    # need matrices for joint test
    final_outcome_model_coefs <- map(est_survreg, broom::tidy, conf.int = TRUE) %>%
      map(., pluck, "estimate") %>%
      map(., ~ matrix(.x))

    names(final_outcome_model_coefs) <- paste0("final_outcome_model_coefs_", names(final_outcome_model_coefs))

    # also need covariance matrices
    final_outcome_model_cov_mtx <- map(est_survreg, vcov) %>%
      map(., as_tibble) %>%
      map(., janitor::clean_names)

    names(final_outcome_model_cov_mtx) <- paste0("final_outcome_model_cov_mtx_", names(final_outcome_model_cov_mtx))

    # export as individual objects to environment
    list2env(c(final_outcome_model_tidy, final_outcome_model_coefs, final_outcome_model_cov_mtx),
             envir = environment()
    )

    # for final subgroup model: need dfs instead of model object in order to unnest the em_results variable
    final_subgroup_model_tidy <- broom::tidy(pi_hat_logistic2, conf.int = TRUE)

    # need matrices for joint test
    # DO NOT USE THIS ROW, re-orders as int1, int2, tmb1, tmb2, instead of int1, tmb1, etc.
    # final_subgroup_model_coefs <- matrix(coef(pi_hat_logistic2))
    final_subgroup_model_coefs <- matrix(broom::tidy(pi_hat_logistic2)$estimate)

    # also need covariance matrix
    final_subgroup_model_cov_mtx <- unname(vcov(pi_hat_logistic2))

    final_subgroup_model_cov_mtx_tidy <- as_tibble(vcov(pi_hat_logistic2)) %>%
      janitor::clean_names()

    # final subgroup assignment from posterior probability
    subgroup_assn <- x1 %>%
      pivot_wider(
        id_cols = c(record_id, contains("latent_subgroup"), assn_subgroup),
        names_from = k,
        names_prefix = "posterior_prob",
        values_from = posterior_prob
      )
  }

  # return objects ----------------------------------------------------------
  return(c(
    list_flatten(list(
      "starting_values" = bind_rows(starting_values_input_df) %>%
        select(contains("hat")),
      # outcome models
      mget(ls(pattern = "final_outcome_model_tidy_")),
      mget(ls(pattern = "final_outcome_model_cov_mtx_")),
      # subgroup membership model
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
      "convergence_iter" = convergence_iter
      # "ests" = ests, # estimates across iterations of algorithm
    )
  )))
} # end em_function
