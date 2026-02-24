## check input parameters ---------------------
test_that("input_model", {
  expect_error(survFMM(model = "apples"),
               "`model_input` must be one of 'AFT-FMM' or 'IPCW-FMM'")

  # check that casing doesn't matter
  set.seed(1213)
  aft_fmm_lower <- survFMM(model = "aft-fmm",
                           input_df = sim_data,
                           weights = "iptw_trim97",
                           outc_model_time = "tt_pfs_days",
                           outc_model_status = "pfs_status",
                           outc_model_covars = "tx",
                           covariates_subgroup_model = "covariate_sim_tmb_zscore",
                           n_inits = 1)

  set.seed(1213)
  aft_fmm_mixed <- survFMM(model = "afT-fmM",
                           input_df = sim_data,
                           weights = "iptw_trim97",
                           outc_model_time = "tt_pfs_days",
                           outc_model_status = "pfs_status",
                           outc_model_covars = "tx",
                           covariates_subgroup_model = "covariate_sim_tmb_zscore",
                           n_inits = 1)

  set.seed(1213)
  aft_fmm_upper <- survFMM(model = "AFT-FMM",
                           input_df = sim_data,
                           weights = "iptw_trim97",
                           outc_model_time = "tt_pfs_days",
                           outc_model_status = "pfs_status",
                           outc_model_covars = "tx",
                           covariates_subgroup_model = "covariate_sim_tmb_zscore",
                           n_inits = 1)

  expect_equal(aft_fmm_lower,
               aft_fmm_mixed)

  expect_equal(aft_fmm_lower,
               aft_fmm_upper)
})

test_that("weights", {
  # error for IPCW-FMM w/o weights
  expect_error(survFMM(model = "IPCW-FMM",
                       input_df = sim_data,
                       # weights = "iptw_trim97",
                       outc_model_time = "tt_pfs_days",
                       outc_model_status = "pfs_status",
                       outc_model_covars = "tx",
                       covariates_subgroup_model = "covariate_sim_tmb_zscore"),
               "^Inverse probability of censoring weights are required")


  # message for AFT-FMM w/o weights
  expect_message(survFMM(model = "AFT-FMM",
                         input_df = sim_data,
                         # weights = "iptw_trim97",
                         outc_model_time = "tt_pfs_days",
                         outc_model_status = "pfs_status",
                         outc_model_covars = "tx",
                         covariates_subgroup_model = "covariate_sim_tmb_zscore",
                         n_inits = 1),
               "^Note: Inverse probability of treatment weights were not supplied")

  # equal if weights of 1 are supplied vs not supplied
  set.seed(0804)
  aft_fmm_no_wts <- survFMM(model = "AFT-FMM",
                            input_df = sim_data,
                            # weights = "iptw_trim97",
                            outc_model_time = "tt_pfs_days",
                            outc_model_status = "pfs_status",
                            outc_model_covars = "tx",
                            covariates_subgroup_model = "covariate_sim_tmb_zscore",
                            n_inits = 1)

  set.seed(0804)
  expect_equal(aft_fmm_no_wts, survFMM(model = "AFT-FMM",
                         input_df = sim_data %>% mutate(weights = 1),
                         weights = "weights",
                         outc_model_time = "tt_pfs_days",
                         outc_model_status = "pfs_status",
                         outc_model_covars = "tx",
                         covariates_subgroup_model = "covariate_sim_tmb_zscore",
                         n_inits = 1))
})

test_that("model terms specified", {
  # missing outcome time
  expect_error(survFMM(model = "AFT-FMM",
                       input_df = sim_data,
                       weights = "iptw_trim97",
                       # outc_model_time = "tt_pfs_days",
                       outc_model_status = "pfs_status",
                       outc_model_covars = "tx",
                       covariates_subgroup_model = "covariate_sim_tmb_zscore",
                       n_inits = 1),
               "^All outcome model terms")

  # missing outcome status
  expect_error(survFMM(model = "AFT-FMM",
                       input_df = sim_data,
                       weights = "iptw_trim97",
                       outc_model_time = "tt_pfs_days",
                       # outc_model_status = "pfs_status",
                       outc_model_covars = "tx",
                       covariates_subgroup_model = "covariate_sim_tmb_zscore",
                       n_inits = 1),
               "^All outcome model terms")

  # missing outcome covariates
  expect_error(survFMM(model = "AFT-FMM",
                       input_df = sim_data,
                       weights = "iptw_trim97",
                       outc_model_time = "tt_pfs_days",
                       outc_model_status = "pfs_status",
                       # outc_model_covars = "tx",
                       covariates_subgroup_model = "covariate_sim_tmb_zscore",
                       n_inits = 1),
               "^All outcome model terms")

  # missing covariates subgroup model
  expect_error(survFMM(model = "AFT-FMM",
                       input_df = sim_data,
                       weights = "iptw_trim97",
                       outc_model_time = "tt_pfs_days",
                       outc_model_status = "pfs_status",
                       outc_model_covars = "tx",
                       # covariates_subgroup_model = "covariate_sim_tmb_zscore",
                       n_inits = 1),
               "^Covariates for the latent subgroup membership model")
})
