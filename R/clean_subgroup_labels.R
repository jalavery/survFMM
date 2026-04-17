#' Clean subgroup labels for survFMM objects
#' A common problem in finite mixture modeling is that subgroup labels are
#' superficial, e.g., the group labeled 'subgroup 1' does not always correspond
#' to the subgroup with the smallest treatment effect. The purpose of this
#' function is to update subgroup labels in order of increasing magnitude of the
#' treatment term such that subgroup 1 always has the smallest treatment term
#' and subgroup 2 always has the largest treatment estimate.
#'
#' @param survFMM_object Object returned from the survFMM function (may be
#'   either AFT-FMM or IPCW-FMM model)
#' @param tx_term The name of the variable corresponding to the treatment
#'
#' @returns aft_fmm_clean A survFMM object with the same list elements as the
#'   input object, but with subgroup labels ordered by increasing magnitude of
#'   the treatment term
#' @export
#'
clean_subgroup_labels <- function(survFMM_object,
                                  tx_term){
  # get list of tidied outcome models
  tidy_outcome_models <- names(survFMM_object)[grepl("final_outcome_model_tidy_", names(survFMM_object))]

  # currently only available for k = 2
  if (max(purrr::pluck(survFMM_object, "final_df")$k)>2){
    stop("Subgroup label alignment is currently only available for 2 subgroups.")
  }

  # if tx term doesn't exist
  if (!any(grepl(tx_term, purrr::pluck(survFMM_object, tidy_outcome_models[1])$term))){
    stop("Check that the `tx_term` input exists in the data and matches the tx term supplied for the survFMM call.")
  }

  # determine which subgroup had the higher tx term
  subgroup_order <- survFMM_object %>%
    # keep results objects corresponding to the outcome models
    purrr::keep_at(tidy_outcome_models) %>%
    # identify which subgroup
    purrr::imap(., ~dplyr::mutate(.x, k = as.numeric(str_remove(pattern = "final_outcome_model_tidy_",
                                                  string = .y)))) %>%
    bind_rows() %>%
    dplyr::filter(term == tx_term) %>%
    dplyr::arrange(estimate) %>%
    dplyr::mutate(k_clean = 1:n())

  # determine if re-ordering is required ---------------------------------
  reorder_flag <- subgroup_order %>%
    dplyr::filter(k != k_clean) %>%
    nrow() > 0

  # correct label switching ---------------------------------
  if (reorder_flag == TRUE){

    # get list of old vs new subgroup names based on k and k_clean
    rename_list <- tidyr::tibble(original = paste0("posterior_prob", subgroup_order$k),
                          new = paste0("posterior_prob", subgroup_order$k_clean))

    # subgroup assignment
    subgroup_assn <- purrr:::pluck(survFMM_object, "subgroup_assn") %>%
      dplyr::rename_with(.fn = ~rename_list$new, .cols = rename_list$original) %>%
      dplyr::rename(assn_subgroup_original = assn_subgroup) %>%
      # merge on clean subgroup assignment
      dplyr::left_join(subgroup_order %>% dplyr::select(k, k_clean),
                by = c("assn_subgroup_original" = "k")) %>%
      dplyr::rename(assn_subgroup = k_clean) %>%
      # latent_sugroup only known in simulations, but want to keep it if it's available
      dplyr::select(record_id, starts_with("latent_subgroup"), assn_subgroup, paste0("posterior_prob", 1:max(subgroup_order$k)))

    # subgroup model
    final_subgroup_model_tidy <- purrr::pluck(survFMM_object, "final_subgroup_model_tidy") %>%
      dplyr::rename(conf.low_orig = conf.low,
             conf.high_orig = conf.high) %>%
      dplyr::mutate(estimate = -estimate,
             conf.low = -conf.high_orig,
             conf.high = -conf.low_orig,
             y.level = 2)

    # flip probabilities in final_df
    final_df <- purrr::pluck(survFMM_object, "final_df") %>%
      dplyr::rename(assn_subgroup_orig = assn_subgroup,
             k_orig = k) %>%
      dplyr::mutate(assn_subgroup = case_when(
        assn_subgroup_orig == 1 ~ 2,
        assn_subgroup_orig == 2 ~ 1
      ),
      k = case_when(
        k_orig == 1 ~ as.character(2),
        k_orig == 2 ~ as.character(1)
      ))

    # new results object
    # combine original list + corrected objects
    aft_fmm_clean <-
      # original list
      c(survFMM_object %>%
      # drop list elements that required correction
      purrr::discard_at(names(survFMM_object)[grepl("final_outcome_model_tidy_|final_outcome_model_cov_mtx_|subgroup_assn|final_subgroup_model_tidy|final_df", names(survFMM_object))]),
      # corrected objects
      tibble::lst(subgroup_assn, final_subgroup_model_tidy, final_df),
      list("final_outcome_model_tidy_1" = pluck(survFMM_object, "final_outcome_model_tidy_2"),
           "final_outcome_model_tidy_2" = pluck(survFMM_object, "final_outcome_model_tidy_1"),
           "final_outcome_model_cov_mtx_1" = pluck(survFMM_object, "final_outcome_model_cov_mtx_2"),
           "final_outcome_model_cov_mtx_2" = pluck(survFMM_object, "final_outcome_model_cov_mtx_1")))

    # order list elements in original order
    aft_fmm_clean <- aft_fmm_clean[names(survFMM_object)]
    # end if for correcting label switching
  } else {
    aft_fmm_clean <- survFMM_object
  }

  return(aft_fmm_clean)
}
