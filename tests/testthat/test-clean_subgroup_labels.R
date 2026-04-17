test_that("clean_subgroup_labels reorders subgroups correctly", {
  set.seed(1213) # with this seed, subgroup labels require switching
  aft_fmm_obj <- survFMM(
    model = "aft-fmm",
    input_df = sim_data,
    weights = "iptw_trim97",
    outc_model_time = "time_to_event_days",
    outc_model_status = "event_status",
    outc_model_covars = "tx",
    covariates_subgroup_model = "covariate_sim_normal",
    n_inits = 1
  )
  cleaned_aft_fmm_obj <- clean_subgroup_labels(aft_fmm_obj, "tx")

  # check that subgroup 2 flipped to subgroup 1 in the outcome models
  expect_equal(cleaned_aft_fmm_obj$final_outcome_model_tidy_1 %>%
                 dplyr::filter(term == "tx") %>%
                 pull(estimate),
               aft_fmm_obj$final_outcome_model_tidy_2 %>%
                 dplyr::filter(term == "tx") %>%
                 pull(estimate))

  # check that subgroup 2 flipped to subgroup 1 in the covariance matrices
  expect_equal(cleaned_aft_fmm_obj$final_outcome_model_cov_mtx_1,
               aft_fmm_obj$final_outcome_model_cov_mtx_2)

  # check that subgroup 1 flipped to subgroup 2 in the outcome models
  expect_equal(cleaned_aft_fmm_obj$final_outcome_model_tidy_2 %>%
                 dplyr::filter(term == "tx") %>%
                 pull(estimate),
               aft_fmm_obj$final_outcome_model_tidy_1 %>%
                 dplyr::filter(term == "tx") %>%
                 pull(estimate))

  # check that subgroup 1 flipped to subgroup 2 in the covariance matrices
  expect_equal(cleaned_aft_fmm_obj$final_outcome_model_cov_mtx_2,
               aft_fmm_obj$final_outcome_model_cov_mtx_1)

  # check that assn_subgroup is updated in final_df
  expect_equal(count(cleaned_aft_fmm_obj$final_df, assn_subgroup) %>%
                 arrange(desc(assn_subgroup)) %>%
                 pull(n),
               count(aft_fmm_obj$final_df, assn_subgroup) %>%
                 pull(n))

  # check subgroup_assn
  expect_equal(cleaned_aft_fmm_obj$subgroup_assn$posterior_prob2,
               aft_fmm_obj$subgroup_assn$posterior_prob1)

  expect_equal(cleaned_aft_fmm_obj$subgroup_assn$posterior_prob1,
               aft_fmm_obj$subgroup_assn$posterior_prob2)

  expect_equal(count(cleaned_aft_fmm_obj$subgroup_assn, assn_subgroup) %>%
                 arrange(desc(assn_subgroup)) %>%
                 pull(n),
               count(aft_fmm_obj$subgroup_assn, assn_subgroup) %>%
                 pull(n))

  # check final_subgroup_model_tidy
  expect_equal(cleaned_aft_fmm_obj$final_subgroup_model_tidy$estimate,
               -aft_fmm_obj$final_subgroup_model_tidy$estimate)

  expect_equal(cleaned_aft_fmm_obj$final_subgroup_model_tidy$conf.low,
               -aft_fmm_obj$final_subgroup_model_tidy$conf.high)

  expect_equal(cleaned_aft_fmm_obj$final_subgroup_model_tidy$conf.high,
               -aft_fmm_obj$final_subgroup_model_tidy$conf.low)
})

test_that("clean_subgroup_labels does not reorder if already correct", {
  set.seed(007) # with this seed, subgroup labels DO NOT require switching
  aft_fmm_obj <- survFMM(
    model = "aft-fmm",
    input_df = sim_data,
    weights = "iptw_trim97",
    outc_model_time = "time_to_event_days",
    outc_model_status = "event_status",
    outc_model_covars = "tx",
    covariates_subgroup_model = "covariate_sim_normal",
    n_inits = 1
  )
  cleaned_aft_fmm_obj <- clean_subgroup_labels(aft_fmm_obj, "tx")

  expect_equal(aft_fmm_obj, cleaned_aft_fmm_obj)
})

test_that("clean_subgroup_labels throws error for more than 2 subgroups", {
  set.seed(1213)
  aft_fmm_obj_k3 <- survFMM(
    model = "aft-fmm",
    input_df = sim_data,
    weights = "iptw_trim97",
    outc_model_time = "time_to_event_days",
    outc_model_status = "event_status",
    outc_model_covars = "tx",
    k = 3,
    covariates_subgroup_model = "covariate_sim_normal",
    n_inits = 1
  )

  expect_error(clean_subgroup_labels(aft_fmm_obj_k3, "tx"),
               "Subgroup label alignment is currently only available for 2 subgroups.")
})

test_that("clean_subgroup_labels handles missing tx_term", {
  set.seed(1213)
  aft_fmm_obj <- survFMM(
    model = "aft-fmm",
    input_df = sim_data,
    weights = "iptw_trim97",
    outc_model_time = "time_to_event_days",
    outc_model_status = "event_status",
    outc_model_covars = "tx",
    covariates_subgroup_model = "covariate_sim_normal",
    n_inits = 1
  )

  expect_error(clean_subgroup_labels(aft_fmm_obj, "abc"), "^Check that the `tx_term` input")
})
