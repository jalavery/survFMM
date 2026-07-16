test_that("clean_subgroup_labels reorders subgroups correctly", {
  cleaned_aft_fmm_obj <- clean_subgroup_labels(aft_fmm_obj_k2_req_label_switch, "tx")

  # check that subgroup 2 flipped to subgroup 1 in the outcome models
  expect_equal(cleaned_aft_fmm_obj$final_outcome_model_tidy_1 %>%
                 dplyr::filter(term == "tx") %>%
                 dplyr::pull(estimate),
               aft_fmm_obj_k2_req_label_switch$final_outcome_model_tidy_2 %>%
                 dplyr::filter(term == "tx") %>%
                 dplyr::pull(estimate))

  # check that subgroup 2 flipped to subgroup 1 in the covariance matrices
  expect_equal(cleaned_aft_fmm_obj$final_outcome_model_cov_mtx_1,
               aft_fmm_obj_k2_req_label_switch$final_outcome_model_cov_mtx_2)

  # check that subgroup 1 flipped to subgroup 2 in the outcome models
  expect_equal(cleaned_aft_fmm_obj$final_outcome_model_tidy_2 %>%
                 dplyr::filter(term == "tx") %>%
                 dplyr::pull(estimate),
               aft_fmm_obj_k2_req_label_switch$final_outcome_model_tidy_1 %>%
                 dplyr::filter(term == "tx") %>%
                 dplyr::pull(estimate))

  # check that subgroup 1 flipped to subgroup 2 in the covariance matrices
  expect_equal(cleaned_aft_fmm_obj$final_outcome_model_cov_mtx_2,
               aft_fmm_obj_k2_req_label_switch$final_outcome_model_cov_mtx_1)

  # check that assn_subgroup is updated in final_df
  expect_equal(dplyr::count(cleaned_aft_fmm_obj$final_df, assn_subgroup) %>%
                 dplyr::arrange(desc(assn_subgroup)) %>%
                 dplyr::pull(n),
               dplyr::count(aft_fmm_obj_k2_req_label_switch$final_df, assn_subgroup) %>%
                 dplyr::pull(n))

  # check subgroup_assn
  expect_equal(cleaned_aft_fmm_obj$subgroup_assn$posterior_prob2,
               aft_fmm_obj_k2_req_label_switch$subgroup_assn$posterior_prob1)

  expect_equal(cleaned_aft_fmm_obj$subgroup_assn$posterior_prob1,
               aft_fmm_obj_k2_req_label_switch$subgroup_assn$posterior_prob2)

  expect_equal(dplyr::count(cleaned_aft_fmm_obj$subgroup_assn, assn_subgroup) %>%
                 dplyr::arrange(desc(assn_subgroup)) %>%
                 dplyr::pull(n),
               dplyr::count(aft_fmm_obj_k2_req_label_switch$subgroup_assn, assn_subgroup) %>%
                 dplyr::pull(n))

  # check final_subgroup_model_tidy
  expect_equal(cleaned_aft_fmm_obj$final_subgroup_model_tidy$estimate,
               -aft_fmm_obj_k2_req_label_switch$final_subgroup_model_tidy$estimate)

  expect_equal(cleaned_aft_fmm_obj$final_subgroup_model_tidy$conf.low,
               -aft_fmm_obj_k2_req_label_switch$final_subgroup_model_tidy$conf.high)

  expect_equal(cleaned_aft_fmm_obj$final_subgroup_model_tidy$conf.high,
               -aft_fmm_obj_k2_req_label_switch$final_subgroup_model_tidy$conf.low)
})

test_that("clean_subgroup_labels does not reorder if already correct", {
  cleaned_aft_fmm_obj <- clean_subgroup_labels(aft_fmm_obj_k2_labels_no_switch, "tx")

  # when passing object to cleaning fcn, even if cleaning labels isn't required,
  # ests df gets subset on last iter since if cleaning is required only the last
  # iter is cleaned in ests
  aft_fmm_obj_k2_labels_no_switch$ests <- aft_fmm_obj_k2_labels_no_switch$ests %>%
    dplyr::slice_max(iter)

  expect_equal(aft_fmm_obj_k2_labels_no_switch, cleaned_aft_fmm_obj)
})

test_that("clean_subgroup_labels handles missing tx_term", {
  set.seed(1213)

  expect_error(clean_subgroup_labels(aft_fmm_obj_k2, "abc"), "^Check that the `tx_term` input")
})

