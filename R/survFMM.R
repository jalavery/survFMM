#' Fit a finite mixture model for survival endpoints
#'
#' This function fits adaptations of the finite mixture model for time-to-event
#' endpoints. One adaptation is a finite mixture of continuous outcomes, fit
#' among individuals with an event and correcting for censoring bias via
#' supplied inverse probability of censoring (IPCW) weights. This approach and
#' its properties are described in
#' \href{https://pubmed.ncbi.nlm.nih.gov/40637678/}{Unveiling non-small cell
#' lung cancer treatment effect heterogeneity: a comparative analysis of
#' statistical methods (2025)}.
#' The second adaptation is a mixture of accelerated failure time models
#' (AFT-FMM), in which censored observations are directly incorporated into the
#' subgroup-specific outcome models.
#'
#' Due to the random initialization of starting values and initial
#' partitions of patients into latent subgroups as part of the EM algorithm,
#' setting a seed prior to calling survFMM is strongly recommended to ensure
#' full reproducibility.
#'
#' @param input_df Input data frame containing 1 row/observation, along with
#'   each observation's event status (0=censored, 1=event) and event time
#'   variables
#' @param weights_input Variable for inverse probability of treatment weights,
#'   if applicable. If nothing is supplied for an AFT-FMM, IPTW will not be
#'   computed. For IPCW-FMM, supply the IPCW. For IPCW-FMM that also uses IPTW,
#'   the product of the ITPW and IPCW may be supplied.
#' @param outc_model_time Variable indicating the time to event or censoring for
#'   each observation.
#' @param outc_model_status Variable indicating the event status for the
#'   time-to-event outcome. It is assumed that 0 = censored, 1 = event.
#' @param outc_model_covars Names of covariates to include in the outcome models
#'   for each subgroup. Note that covariates should be continuous or in the form
#'   of numeric dummy variables.
#' @param outc_distribution Outcome distribution for subgroup-specific outcome
#'   models. Currently allowed values are "Weibull" and "Log-Normal" (not
#'   case-sensitive)
#' @param covariates_subgroup_model Names of covariates to include in subgroup
#'   membership model.
#' @param model "AFT-FMM" for a finite mixture of accelerated failure time (AFT)
#'   models, "IPCW-FMM" for a finite mixture of a continuous distribution
#'   weighted by inverse probability of censoring weights (IPCW). Input is not
#'   case sensitive. Default is AFT-FMM.
#' @param k Number of subgroups. Default is 2.
#' @param starting_values_window The percent margin around the starting values.
#'   For example, starting_values_type = 'single_survreg' and
#'   starting_values_window = 0.5 means that starting values are randomly
#'   generated uniformly around +\/- 50\% of the estimates from a single survreg
#'   model fit. Default is 1 (i.e., starting values +\/- 100\%).
#' @param starting_values_type One of "single_survreg", "uniform_pct", or
#'   "non_random_start". "single_survreg" fits a single AFT model to all of the
#'   data and then generates random starting values based on the
#'   `starting_values_window` parameter for each initial partition. If not
#'   supplying starting values, they are randomly generated for each initial
#'   partition. Be sure to set a seed at the top of your script to ensure
#'   reproducibility. Default is "single_survreg."
#' @param starting_values_df Input dataset with starting values for algorithm
#' @param n_inits Number of initial partitions for the EM algorithm. Default is
#'   5. A higher number of initial partitions may result in greater stability of
#'   estimates.
#' @param tolerance Convergence criteria for the change in the log-likelihood
#'   for the EM algorithm. Default is 0.001.
#' @param conv_pct_criteria Convergence criteria for the percentage of
#'   observations changing subgroup. Specify -1 to only use the log-likelihood
#'   as the convergence criteria.
#' @param max_iter Maximum number of iterations. Default is 200.
#' @param save_all_init Whether results for all initial partitions are saved.
#'   Default is FALSE.
#' @return List of results with the following components:
#' \itemize{
#'    \item starting_values: Dataframe with starting values used for algorithm initialization
#'    \item final_outcome_model_1-final_outcome_model_k: Outcome model objects
#'    for each latent subgroup, 1 to k
#'    \item final_outcome_model_tidy_1-final_outcome_model_tidy_k: Tidy dataframe
#'   corresponding to the outcome model for each latent subgroup, 1 to k
#'    \item final_outcome_model_cov_mtx_1-final_outcome_model_cov_mtx_k: Tibble of
#'   the covariance matrix corresponding to the outcome model for each latent
#'   subgroup, 1 to k
#'    \item  subgroup_assn: Dataframe containing the posterior probability of
#'   subgroup membership and corresponding assigned subgroup (based on maximum
#'   posterior probability) for each observation
#'    \item final_subgroup_model: Final model for subgroup membership
#'    \item final_subgroup_model_tidy: Tidy dataframe corresponding to the final
#'   model for subgroup membership
#'    \item  final_df: One observation per record, per subgroup containing the
#'   input dataset and corresponding prior and posterior probabilities.
#'    \item  log_likelihood_values: Vector of log-likelihood values across
#'   algorithm iterations
#'    \item  convergence_status: Numeric convergence status. 0 = did not
#'   converge, 1 = algorithm converged
#'    \item  convergence_message: Message indicating whether algorithm converged
#'    \item  convergence_iter: Final iteration of the algorithm. Either the
#'   iteration that the algorithm achieved convergence, the final iteration
#'   following an error, or the maximum number of iterations (non-convergence)
#'}
#' @export
#'
#' @examples
#' # Examples using package test data
#' # Example 1 ----------------------------------
#' # Fit a mixture of accelerated failure time models
#' ex_aft_fmm <- survFMM(
#'                  model = "AFT-FMM",
#'                  input_df = sim_data,
#'                  weights = "iptw_trim97",
#'                  outc_distribution = "weibull",
#'                  outc_model_time = "time_to_event_days",
#'                  outc_model_status = "event_status",
#'                  outc_model_covars = "tx",
#'                  covariates_subgroup_model = "covariate_sim_normal",
#'                  n_inits = 5)
#' #' # Example 2 ----------------------------------
#' # Fit a mixture of Weibull models, weighted by the inverse probability of
#' # censoring
#' ex_ipcw_fmm <- survFMM(
#'                  model = "IPCW-FMM",
#'                  input_df = sim_data,
#'                  weights = "iptw_ipcw_trim97",
#'                  outc_distribution = "weibull",
#'                  outc_model_time = "time_to_event_days",
#'                  outc_model_status = "event_status",
#'                  outc_model_covars = "tx",
#'                  covariates_subgroup_model = "covariate_sim_normal",
#'                  n_inits = 5)
survFMM <- function(input_df,
                    weights_input = NULL,
                    outc_model_time = NULL,
                    outc_model_status = NULL,
                    outc_model_covars = NULL,
                    outc_distribution = "weibull",
                    covariates_subgroup_model = NULL,
                    model = "AFT-FMM",
                    k = 2,
                    starting_values_type = "single_survreg",
                    starting_values_window = 1,
                    starting_values_df = NULL,
                    n_inits = 5,
                    tolerance = 1e-3,
                    conv_pct_criteria = -1,
                    max_iter = 200,
                    save_all_init = FALSE) {
  # initialization ----------------------------------------------------------
  # take any casing of the model input and distribution parameters
  model_input <- stringr::str_to_lower(model)
  outc_distribution <- stringr::str_to_lower(outc_distribution)

  # error checking
  if (!(model_input %in% c("aft-fmm", "ipcw-fmm"))){
    stop("`model_input` must be one of 'AFT-FMM' or 'IPCW-FMM'")
  }

  if (!(outc_distribution %in% c("weibull", "lognormal"))){
    stop("Currently only the Weibull and lognormal distributions are supported.")
  }

  if (!(starting_values_type %in% c("single_survreg", "uniform_pct"))){
    stop("Input parameter `starting_values_type` must be one of: single_survreg or uniform_pct.")
  }

  if (is.null(input_df)){
    stop("Input dataframe is required. Please specify the `input_df` parameter.")
  }

  if (model_input == "aft-fmm" & is.null(weights_input)){
    message("Note: Inverse probability of treatment weights were not supplied and will not be computed.")

    # assign weights to 1
    input_df <- input_df %>%
      dplyr::mutate(weights = 1)

    weights_input <- "weights"
  }

  if (model_input == "ipcw-fmm" & is.null(weights_input)){
    stop("Inverse probability of censoring weights are required to be supplied when `method_input` = 'IPCW-FMM'.")
  }

  if (is.null(covariates_subgroup_model)){
    stop("Covariates for the latent subgroup membership model must be specified in the `covariates_subgroup_model` parameter.")
  }

  if (is.null(outc_model_time) | is.null(outc_model_status) | is.null(outc_model_covars)){
    stop("All outcome model terms (`outc_model_time`, `outc_model_status`, `outc_model_covars`) must be specified")
  }

  if (!any(grepl("record_id", names(input_df)))){
    # assign a record_id variable if not present in the input dataset
    input_df <- input_df %>%
      dplyr::mutate(record_id = 1:dplyr::n())
  } else {
    input_df <- input_df %>%
      dplyr::arrange(record_id)
  }

  # define outcome model formula for use in survreg when calling
  # starting_values_df_list (for single_survreg start) and when calling
  # fmm_em_algorithm for fitting outcome models in E-step
  outc_model_formula <- stats::as.formula(
    paste0(
      "survival::Surv(", outc_model_time, ", ", outc_model_status, ") ~ ",
      paste(outc_model_covars, collapse = " + ")
    )
  )

  # initialize starting values ----------------------------------------------------
  starting_values_df_list <- initialize_starting_values(n_inits = n_inits,
                             k = k,
                             starting_values_df = starting_values_df,
                             starting_values_type = starting_values_type,
                             starting_values_window = starting_values_window,
                             input_df = input_df,
                             outc_model_formula = outc_model_formula,
                             weights_input = weights_input,
                             outc_distribution = outc_distribution)


  # subset to events if ipcw-fmm
  # do this after setting up starting values
  if (model_input %in% c("ipcw-fmm", "ipcw_fmm")) {
    input_df <- input_df %>%
      dplyr::filter(.data[[outc_model_status]] == 1)
  }

  # input_df_list <- replicate(n_inits, input_df, simplify = FALSE)

  # loop over iterations of EM algorithm ----------------------------------------------------
  # browser()
  # call em_function
  em_diff_start <- purrr::map(starting_values_df_list,
    ~ try(fmm_em_algorithm(input_df = input_df,
                           starting_values_input_df = .x,
                           k = k,
                           outc_model_formula = outc_model_formula,
                           outc_model_time = outc_model_time,
                           outc_model_status = outc_model_status,
                           outc_model_covars = outc_model_covars,
                           weights_input = weights_input,
                           outc_distribution = outc_distribution,
                           covariates_subgroup_model = covariates_subgroup_model,
                           n_inits = n_inits,
                           tolerance = tolerance,
                           conv_pct_criteria = conv_pct_criteria,
                           max_iter = max_iter),
      silent = FALSE
    ),
    .progress = TRUE
  ) # end purrr::map over each initial partition

  # browser()

  # list of initial partitions that converged
  init_partit_converged_index <- purrr::map(em_diff_start, purrr::pluck, "convergence_status") %>%
    tibble::enframe(
      name = "initial_partition",
      value = "convergence_status"
    ) %>%
    tidyr::unnest(cols = convergence_status) %>%
    dplyr::filter(.data$convergence_status == 1)

  # only keep initial partitions that converged
  # (keep all that converged for now; will select best based on log-likelihood in next step)
  em_diff_start_converged <- em_diff_start[init_partit_converged_index$initial_partition]

  # get initial start with highest log-likelihood
  max_log_l <- purrr::map(em_diff_start_converged, purrr::pluck, "log_likelihood_values") %>%
    # get final log-likelihood from each initial split (might not be the max if
    # log-lik decreases across iterations)
    purrr::map(., utils::tail, 1) %>%
    # get max final log-lik from each initial split
    # purrr::map(., max) %>%
    # select the initial split with the largest log-lik
    which.max()

  # extract initial partition resulting in highest log-likelihood for each repetition
  if (length(max_log_l)) {
    best_initial_partition <- purrr::pluck(em_diff_start_converged, max_log_l)
  } else {
    # if no initial splits converged
    best_initial_partition <- purrr::pluck(em_diff_start, 1) # [c("convergence_status",
    # "convergence_message")]
  }

  # name initializations for em_diff_start
  names(em_diff_start) <- paste0("ninit_", 1:n_inits)

  # browser()

  # return all initial partitions AND the best initial partition (for now)
  # save pred_time_to_event as its own column, don't save it for all iterations
  # or as part of em_results
  # https://stackoverflow.com/questions/48467884/remove-an-element-of-a-list-by-name
  if (save_all_init == TRUE) {
    return(tibble::lst(
      "all_init_partitions" = em_diff_start,
      best_initial_partition
    ))
  } else {
    # return(best_initial_partition)
    return(purrr::discard_at(best_initial_partition, at = c("pi_hat_ests")))
  }
}
