#' Clean subgroup labels for survFMM objects
#'
#' A common problem in finite mixture modeling is that subgroup labels are
#' superficial, e.g., the group labeled 'subgroup 1' does not always correspond
#' to the subgroup with the smallest treatment effect. The purpose of this
#' function is to update subgroup labels in order of increasing magnitude of the
#' treatment term such that subgroup 1 always has the smallest treatment term
#' and subgroup 2 always has the largest treatment estimate. Note: this function
#' is currently only compatible when `save_all_init` = FALSE on the `survFMM()`
#' call.
#'
#' @param survFMM_object Object returned from the survFMM function (may be
#'   either AFT-FMM or IPCW-FMM model)
#' @param tx_term The name of the variable corresponding to the treatment
#'
#' @returns aft_fmm_clean A survFMM object with the same list elements as the
#'   input object, but with subgroup labels ordered by increasing magnitude of
#'   the treatment term
#'
#' @export
clean_subgroup_labels <- function(survFMM_object,
                                  tx_term){
  # get list of tidied outcome models
  tidy_outcome_models <- names(survFMM_object)[grepl("final_outcome_model_tidy_", names(survFMM_object))]

  # if tx term doesn't exist
  if (!any(grepl(tx_term, purrr::pluck(survFMM_object, tidy_outcome_models[1])$term))){
    stop("Check that the `tx_term` input exists in the data and matches the tx term supplied for the survFMM call.")
  }

  # browser()

  # determine which subgroup had the higher tx term
  subgroup_order <- survFMM_object %>%
    # keep results objects corresponding to the outcome models
    purrr::keep_at(tidy_outcome_models) %>%
    # identify which subgroup
    purrr::imap(., ~dplyr::mutate(.x, k = as.numeric(stringr::str_remove(pattern = "final_outcome_model_tidy_",
                                                  string = .y)))) %>%
    dplyr::bind_rows() %>%
    dplyr::filter(.data$term == tx_term) %>%
    dplyr::arrange(.data$estimate) %>%
    dplyr::mutate(k_clean = 1:dplyr::n())

  # determine if re-ordering is required ---------------------------------
  reorder_flag <- subgroup_order %>%
    dplyr::filter(.data$k != .data$k_clean) %>%
    nrow() > 0

  # correct label switching ---------------------------------
  if (reorder_flag == TRUE){

    # flip probabilities in final_df
    final_df_list <- dplyr::bind_rows(
      purrr::pluck(survFMM_object, "final_df") %>% dplyr::mutate(iter = "final"),
      purrr::pluck(survFMM_object, "x1_prev_iter") %>% dplyr::mutate(iter = "final-1")) %>%
      dplyr::rename(assn_subgroup_orig = assn_subgroup,
                    k_orig = k) %>%
      # merge on updated subgroup labels for k
      dplyr::left_join(.,
                       subgroup_order %>%
                         dplyr::select(k, k_clean) %>%
                         dplyr::mutate(k = as.character(k)),
                       by = c("k_orig" = "k")) %>%
      dplyr::mutate(k = as.character(k_clean)) %>%
      # now merge updated subgroup labels on for assn_subgroup
      dplyr::left_join(.,
                       subgroup_order %>%
                         dplyr::select(k, k_clean),
                       by = c("assn_subgroup_orig" = "k")) %>%
      dplyr::mutate(assn_subgroup = k_clean.y) %>%
      dplyr::select(-contains("k_clean"), -k_orig, -assn_subgroup_orig) %>%
      split(.$iter)

    final_df <- final_df_list$`final` %>% dplyr::select(-iter)

    # get list of old vs new subgroup names based on k and k_clean
    rename_list <- tidyr::tibble(
      original = paste0("posterior_prob", subgroup_order$k),
      new = paste0("posterior_prob", subgroup_order$k_clean))

    # subgroup assignment
    subgroup_assn <- purrr::pluck(survFMM_object, "subgroup_assn") %>%
      dplyr::rename_with(.fn = ~rename_list$new, .cols = rename_list$original) %>%
      dplyr::rename(assn_subgroup_original = assn_subgroup) %>%
      # merge on clean subgroup assignment
      dplyr::left_join(subgroup_order %>% dplyr::select(k, k_clean),
                by = c("assn_subgroup_original" = "k")) %>%
      dplyr::rename(assn_subgroup = k_clean) %>%
      # latent_subgroup only known in simulations, but want to keep it if it's available
      dplyr::select(record_id, dplyr::starts_with("latent_subgroup"),
                    assn_subgroup, paste0("posterior_prob", 1:max(subgroup_order$k)))

    # subgroup model: use final_df created above with corrected levels
    final_subgroup_model <- (nnet::multinom(k ~ . - posterior_prob,
                                            data = final_df_list$`final-1` %>%
                                              dplyr::select(
                                                k,
                                                tidyr::all_of(
                                                  labels(terms(purrr::pluck(survFMM_object, "final_subgroup_model")))
                                                ),
                                                posterior_prob
                                              ),
                                            weights = posterior_prob,
                                            trace = FALSE
    ))

    # tidy updated subgroup model
    final_subgroup_model_tidy <- broom::tidy(final_subgroup_model, conf.int = TRUE)

    # subgroup model covariance matrix
    final_subgroup_model_cov_mtx <- unname(stats::vcov(final_subgroup_model))

    # outcome model objects
    final_outcome_model_rename_list <- tidyr::tibble(
      original = paste0("final_outcome_model_tidy_", subgroup_order$k),
      new = paste0("final_outcome_model_tidy_", subgroup_order$k_clean)) %>%
      dplyr::arrange(original)

    final_outcome_model_list_clean <- survFMM_object %>%
      purrr::keep_at(names(survFMM_object)[grepl("final_outcome_model_tidy", names(survFMM_object))]) %>%
      setNames(., final_outcome_model_rename_list$new) %>%
      .[order(names(.))]

    # covariance matrix (of outcome models) objects
    final_outcome_model_covar_mtx_rename_list <- tidyr::tibble(
      original = paste0("final_outcome_model_cov_mtx_", subgroup_order$k),
      new = paste0("final_outcome_model_cov_mtx_", subgroup_order$k_clean)) %>%
      dplyr::arrange(original)

    final_outcome_model_cov_mtx_clean <- survFMM_object %>%
      purrr::keep_at(names(survFMM_object)[grepl("final_outcome_model_cov_mtx", names(survFMM_object))]) %>%
      setNames(., final_outcome_model_covar_mtx_rename_list$new) %>%
      .[order(names(.))]

    # also clean ests df
    # can only return final iter b/c ow would have to update all latent subgroup model ests
    ests <- bind_rows(
      purrr::pluck(survFMM_object, "ests") %>%
        dplyr::filter(model == "Outcome Model") %>%
        dplyr::slice_max(iter) %>%
        dplyr::rename(assn_subgroup_orig = assn_subgroup,
                      k_orig = k) %>%
        # merge on updated subgroup labels for k
        dplyr::left_join(.,
                         subgroup_order %>%
                           dplyr::select(k, k_clean) %>%
                           dplyr::mutate(k = as.character(k)),
                         by = c("k_orig" = "k")) %>%
        dplyr::mutate(k = as.character(k_clean)) %>%
        # now merge updated subgroup labels on for assn_subgroup
        dplyr::left_join(.,
                         subgroup_order %>%
                           dplyr::select(k, k_clean),
                         by = c("assn_subgroup_orig" = "k")) %>%
        dplyr::mutate(assn_subgroup = k_clean.y) %>%
        dplyr::select(-contains("k_clean"), -k_orig, -assn_subgroup_orig),
      # subgroup labels already corrected above
      final_subgroup_model_tidy %>%
        dplyr::mutate(
          model = "Latent Subgroup Model",
          assn_subgroup = as.numeric(.data$y.level)
        ) %>%
        dplyr::select(-y.level)
    ) %>%
      tidyr::fill(iter, .direction = "down") %>%
      dplyr::arrange(desc(model), k)

    # new results object
    # combine original list + corrected objects
    aft_fmm_clean <-
      # original list
      c(survFMM_object %>%
      # drop list elements that required correction
      purrr::discard_at(names(survFMM_object)[grepl("final_outcome_model_tidy_|final_outcome_model_cov_mtx_|subgroup_assn|final_subgroup_model|final_subgroup_model_tidy|final_subgroup_model_cov_mtx|final_df|ests", names(survFMM_object))]),
      # corrected objects
      tibble::lst(subgroup_assn, final_subgroup_model, final_subgroup_model_tidy, final_subgroup_model_cov_mtx,
                  final_df, ests),
      # these objects are fully flipped, no manipulation required, only re-naming, done above
      final_outcome_model_list_clean,
      final_outcome_model_cov_mtx_clean)

    # order list elements in original order
    aft_fmm_clean <- aft_fmm_clean[names(survFMM_object)]
    # end if for correcting label switching
  } else {
    aft_fmm_clean <- c(survFMM_object %>%
                         purrr::discard_at(~str_detect(.x, "ests")),
                       list("ests" = purrr::pluck(survFMM_object, "ests") %>%
                         dplyr::slice_max(iter)))

    # order list elements in original order
    aft_fmm_clean <- aft_fmm_clean[names(survFMM_object)]
  }

  return(aft_fmm_clean)
}
