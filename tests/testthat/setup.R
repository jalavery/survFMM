set.seed(1213) # with this seed, subgroup labels require switching
aft_fmm_obj_k2 <- survFMM(
  model = "aft-fmm",
  input_df = sim_data,
  weights = "iptw_trim97",
  outc_model_time = "time_to_event_days",
  outc_model_status = "event_status",
  outc_model_covars = "tx",
  covariates_subgroup_model = "covariate_sim_normal",
  n_inits = 1
)

aft_fmm_obj_k2_req_label_switch <- aft_fmm_obj_k2


set.seed(007) # with this seed, subgroup labels DO NOT require switching
aft_fmm_obj_k2_labels_no_switch <- survFMM(
  model = "aft-fmm",
  input_df = sim_data,
  weights = "iptw_trim97",
  outc_model_time = "time_to_event_days",
  outc_model_status = "event_status",
  outc_model_covars = "tx",
  covariates_subgroup_model = "covariate_sim_normal",
  n_inits = 1
)

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
