# IPM -- Starkey Elk
# Kenneth Loonam
# April 2021

#Variables======================================================================

# Specify the model
model_file <- "models//ipm//ipm_nimble.R"
save_file <- "results//ipm_result_09may2024.Rdata"

# Loop dimension parameters
n_year <- 34

# JAGS control parameters
n_i <- 500
n_a <- 50
n_b <- 100
n_c <- 3
n_t <- 1

#Environment====================================================================

require(tidyverse); require(nimble); require(mcmcplots)
load("data//elk_ipm_data_21apr2023.Rdata")

#Data_prep======================================================================

nimble_data <- list(
  s_cjs      = ipm_data$s_cjs,
  r_ratio    = ipm_data$r_ratio,
  n_sight_ca = ipm_data$n_sight_ca,
  n_sight_am = ipm_data$n_sight_am,
  n_sight_af = ipm_data$n_sight_af,
  n_a_mov    = ipm_data$n_ad_add - abs(ipm_data$n_ad_rem),
  n_c_mov    = ipm_data$n_ca_add - abs(ipm_data$n_ca_rem),
  n_har      = ipm_data$n_hnt,
  min_ad     = ipm_data$min_ad,
  min_ca     = ipm_data$min_ca,
  prcp       = ipm_data$summer_precip,
  temp       = ipm_data$summer_temp,
  af_count   = ipm_data$n_f_p_count,
  am_count   = ipm_data$n_m_p_count,
  cdens      = ipm_data$cougar_density,
  n_adj      = ipm_data$af_density,
  min_n1     = ipm_data$min_n1,
  est_n1     = ipm_data$est_n1,
  clim       = ipm_data$palmer_index
)

nimble_constants <- list(
  n_year     = n_year,
  nn_ca      = nrow(ipm_data$n_sight_ca),
  nn_af      = nrow(ipm_data$n_sight_af),
  nn_am      = nrow(ipm_data$n_sight_am),
  ns         = nrow(ipm_data$s_cjs),
  nr         = nrow(ipm_data$r_ratio),
  nn_fc      = nrow(ipm_data$n_f_p_count),
  nn_mc      = nrow(ipm_data$n_m_p_count)
)


nimble_inits <- function(){
  # Starting abundance
  N <- array(data = NA, dim = c(4,2))
  N[,] <- 500
  
  # Recruitment
  R <- 0.9
  
  # Survival
  S_C_B0_ps <- 0.99
  S_Y_F_B0_ps <- 0.99
  S_Y_M_B0_ps <- 0.99
  S_A_F_B0_ps <- 0.99
  S_A_M_B0_ps <- 0.99
  
  out <- list(
    init_N = N,
    R_B0_ps = R,
    S_C_B0_ps = S_C_B0_ps,
    S_Y_F_B0_ps = S_Y_F_B0_ps,
    S_Y_M_B0_ps = S_Y_M_B0_ps,
    S_A_F_B0_ps = S_A_F_B0_ps,
    S_A_M_B0_ps = S_A_M_B0_ps
  )
  return(out)
}

initial_values <- inits()

params = c(
  "N_tot",
  "survival_ca",
  "survival_af",
  "survival_am",
  "R",
  "N_f",
  "N_c",
  "N_m",
  "N_lam",
  "R_B0",
  "R_wt",
  "R_wm",
  "R_dd",
  "R_cg",
  "S_C_0",
  "S_wt",
  "S_wm",
  "S_dd",
  "S_cg",
  "lambda",
  "LAMBDA_mean",
  "pop_r",
  "mean_pop_r",
  "R_mean",
  "S_F_mean",
  "S_M_mean",
  "S_C_mean"
)

source(model_file)

#Model==========================================================================

rslt <- nimbleMCMC(
  code = nimble_code,
  constants = nimble_constants,
  data = nimble_data,
  inits = nimble_inits,
  monitors = params,
  thin = n_t,
  nburnin = n_b,
  niter = n_i,
  nchains = n_c
)

mcmcplots::mcmcplot(rslt)

# data <- summary(rslt)
# 
# q <- data$quantiles
# 
# n_tot <- q[grep("N_tot", dimnames(q)[[1]]), c(1,3,5)] %>%
#   as_tibble() %>%
#   mutate(year = 1988:2021)


save(rslt, file = save_file)

