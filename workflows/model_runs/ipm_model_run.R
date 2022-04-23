# IPM -- Starkey Elk
# Kenneth Loonam
# April 2021

#Variables======================================================================

# Specify the model
model_file <- "models//ipm//ipm_elk_null_in_progress.txt"
save_file <- "results//ipm_result_22apr2022.Rdata"

# Loop dimension parameters
n_year <- 34

# JAGS control parameters
n_i <- 500000
n_a <- 10000
n_b <- 500000
n_c <- 3
n_t <- 100

#Environment====================================================================

require(tidyverse); require(rjags); require(mcmcplots)
# load("data//elk_ipm_data.Rdata")
load("data//elk_ipm_data_12apr2022.Rdata")

#Data_prep======================================================================

jags_data <- list(
  s_cjs = ipm_data$s_cjs,
  r_ratio = ipm_data$r_ratio,
  n_sight_ca = ipm_data$n_sight_ca,
  n_sight_am = ipm_data$n_sight_am,
  n_sight_af = ipm_data$n_sight_af,
  n_a_mov = ipm_data$n_ad_add - abs(ipm_data$n_ad_rem),
  n_c_mov = ipm_data$n_ca_add - abs(ipm_data$n_ca_rem),
  n_year = n_year,
  nn_ca = nrow(ipm_data$n_sight_ca),
  nn_af = nrow(ipm_data$n_sight_af),
  nn_am = nrow(ipm_data$n_sight_am),
  ns = nrow(ipm_data$s_cjs),
  nr = nrow(ipm_data$r_ratio),
  n_har = ipm_data$n_hnt,
  min_ad = ipm_data$min_ad,
  min_ca = ipm_data$min_ca,
  x1 = ipm_data$annual_temp
)

inits <- function(){
  # Starting abundance
  N <- array(data = NA, dim = c(4,2))
  N[,] <- 500
  
  # Recruitment
  R <- 0.9
  
  # # Harvest
  # p_har <- array(data = 0, dim = c(3,2,33))
  # p_har[,,1] <- NA
  
  # # Removals
  # p_rem <- array(data = 0.01, dim = c(2,2,33))
  # p_rem[,,1] <- NA
  
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
    # p_har = p_har,
    # p_rem = p_rem
  )
  return(out)
}

initial_values <- inits()

params = c(
  "N_tot",
  "survival_ca",
  "survival_af",
  "survival_yf",
  "survival_am",
  "survival_ym",
  "R",
  "N_f",
  "N_c",
  "N_yf",
  "R_B1"
)

#Model==========================================================================

jgs_mdl <- jags.model(
  file = model_file,
  data = jags_data,
  inits = inits,
  n.chains = n_c,
  n.adapt = n_a
)

update(jgs_mdl, n.iter = n_b)

rslt <- coda.samples(
  jgs_mdl,
  variable.names = params,
  n.iter = n_i,
  thin = n_t
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
