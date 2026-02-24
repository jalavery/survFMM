## code to prepare `sim_data` dataset goes here

load("H:/Biostatistics/Jessica L/Lavery/dissertation/data/Project 2/Simulation/simulated_data/S1172_2156_32112_TxIPTW2Opp4_LStr_E30_N1500.RData")

sim_data <- S1172_2156_32112_TxIPTW2Opp4_LStr_E30_N1500 %>%
  filter(repetition ==1) %>%
  unnest(data)

usethis::use_data(sim_data, overwrite = TRUE)
