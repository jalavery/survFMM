## code to prepare `sim_data` dataset goes here
library(tidyverse)
load(here::here("data-raw/S1172_2156_32112_TxIPTW2Opp4_LStr_E30_N1500.RData"))

sim_data <- S1172_2156_32112_TxIPTW2Opp4_LStr_E30_N1500 %>%
  dplyr::filter(repetition == 1) %>%
  tidyr::unnest(data) %>%
  dplyr::rename(covariate_sim_normal = covariate_sim_tmb_zscore,
         covariate_sim_binary1 = covariate_sim_sex_binary,
         covariate_sim_binary2 = covariate_sim_histology_binary,
         event_status = pfs_status,
         time_to_event_days = tt_pfs_days) %>%
  dplyr::select(-repetition, -setting_name) %>%
  dplyr::select(record_id, latent_subgroup,
         covariate_sim_normal, covariate_sim_binary1, covariate_sim_binary2,
         tx,
         time_to_event_days, event_status,
         iptw_trim97, ipcw_trim97,
         iptw_ipcw_trim97)

usethis::use_data(sim_data, overwrite = TRUE)
