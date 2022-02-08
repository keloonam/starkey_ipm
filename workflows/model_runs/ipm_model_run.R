# IPM -- Starkey Elk
# Kenneth Loonam
# April 2021

#Variables======================================================================

# Specify the model
model_file <- "models//ipm//ipm_elk_null_in_progress.txt"
save_file <- paste0("results//ipm_result_", date(), ".Rdata")

# Loop dimension parameters
n_year <- 33

# JAGS control parameters
n_i <- 5000
n_a <- 1000
n_b <- 10000
n_c <- 3
n_t <- 10

#Environment====================================================================

require(tidyverse); require(rjags); require(mcmcplots)
# load("data//elk_ipm_data.Rdata")
load("data//elk_ipm_data_07jan2022.Rdata")

#Data_prep======================================================================

# Build data for jags

jags_data <- list(
  s_cjs = ipm_data$s_cjs,
  r_ratio = ipm_data$r_ratio,
  n_sight_ca = ipm_data$n_sight_ca,
  n_sight_am = ipm_data$n_sight_am,
  n_sight_af = ipm_data$n_sight_af,
  n_fg_c = ipm_data$n_fg_c,
  nn_fg_ca = nrow(ipm_data$n_fg_c),
  n_fg_f = ipm_data$n_fg_f,
  nn_fg_f = nrow(ipm_data$n_fg_f),
  n_fg_m = ipm_data$n_fg_m,
  nn_fg_m = nrow(ipm_data$n_fg_m),
  n_ad_add = ipm_data$n_ad_add,
  n_ca_add = ipm_data$n_ca_add,
  n_ad_rem = abs(ipm_data$n_ad_rem),
  n_ca_rem = abs(ipm_data$n_ca_rem),
  n_year = n_year,
  nn_ca = nrow(ipm_data$n_sight_ca),
  nn_af = nrow(ipm_data$n_sight_af),
  nn_am = nrow(ipm_data$n_sight_am),
  ns = nrow(ipm_data$s_cjs),
  nr = nrow(ipm_data$r_ratio),
  n_har = ipm_data$n_hnt[2:4,,]
)

inits <- function(){
  # Starting abundance
  N <- array(data = NA, dim = c(4,2))
  N[,] <- 500
  
  # Recruitment
  R <- rep(0.99, 33)
  
  # Harvest
  p_har <- array(data = 0.01, dim = c(3,2,33))
  p_har[,,1] <- NA
  
  # Removals
  p_rem <- array(data = 0.01, dim = c(2,2,33))
  p_rem[,,1] <- NA
  
  # Survival
  S_C_B0_ps <- 0.99
  S_Y_F_B0_ps <- 0.99
  S_Y_M_B0_ps <- 0.99
  S_A_F_B0_ps <- 0.99
  S_A_M_B0_ps <- 0.99
  
  out <- list(
    init_N = N,
    R = R,
    S_C_B0_ps = S_C_B0_ps,
    S_Y_F_B0_ps = S_Y_F_B0_ps,
    S_Y_M_B0_ps = S_Y_M_B0_ps,
    S_A_F_B0_ps = S_A_F_B0_ps,
    S_A_M_B0_ps = S_A_M_B0_ps,
    p_har = p_har,
    p_rem = p_rem
  )
  return(out)
}

initial_values <- inits()

params = c(
  "N_tot",
  "survival_ca",
  "survival_af",
  "R"
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


save(rslt, file = save_file)

# plot(rslt$BUGSoutput$summary[grep("N_tot", dimnames(rslt$BUGSoutput$summary)[[1]]),1], type = "l")
# plot(rslt$BUGSoutput$summary[grep("R", dimnames(rslt$BUGSoutput$summary)[[1]]),1], type = "l")
# plot(rslt$BUGSoutput$summary[grep("survival_af", dimnames(rslt$BUGSoutput$summary)[[1]]),1], type = "l")
# plot(rslt$BUGSoutput$summary[grep("survival_ca", dimnames(rslt$BUGSoutput$summary)[[1]]),1], type = "l")
# 
# 
# jags_data$n_fg_c[,4] / jags_data$n_fg_f[,4]
# jags_data$r_ratio[,4]
# jags_data$s_cjs

