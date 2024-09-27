# IPM -- Starkey Elk
# Kenneth Loonam
# April 2021

#Variables======================================================================

# Specify the model
model_file <- "s2//ipm.txt"

climate_variable <- "spei12"
 # pdsi, precip, temp, ndvi, spei3, spei6, spei12
save_file <- paste0("s2//ipmrs_25sep2025_", climate_variable, ".rds")

# Loop dimension parameters
n_year <- 36

# JAGS control parameters
n_i <- 1000000
n_a <- 100000
n_b <- 500000
n_c <- 3
n_t <- 100

#Environment====================================================================

require(tidyverse); require(rjags); require(mcmcplots)
ipm_data <- readRDS("s2//ipm_data_25sep2024.rds")

#Data_prep======================================================================
jags_data <- list(
  s_cjs      = ipm_data$s_cjs,
  r_ratio    = ipm_data$r_ratio,
  n_sight_ca = ipm_data$n_sight_ca,
  n_sight_am = ipm_data$n_sight_am,
  n_sight_af = ipm_data$n_sight_af,
  n_a_mov    = ipm_data$n_a_mov,
  n_c_mov    = ipm_data$n_c_mov,
  n_year     = ipm_data$n_year,
  nn_ca      = ipm_data$nn_ca,
  nn_af      = ipm_data$nn_af,
  nn_am      = ipm_data$nn_am,
  ns         = nrow(ipm_data$s_cjs),
  nr         = nrow(ipm_data$r_ratio),
  n_har      = ipm_data$n_har,
  min_ad     = ipm_data$min_ad,
  min_ca     = ipm_data$min_ca,
  af_count   = ipm_data$af_count,
  nn_fc      = ipm_data$nn_fc,
  am_count   = ipm_data$am_count,
  nn_mc      = ipm_data$nn_mc,
  cdens      = ipm_data$cdens,
  n_adj      = ipm_data$nelk,
  min_n1     = ipm_data$min_n1,
  est_n1     = ipm_data$est_n1
)


inits <- function(){
  # Starting abundance
  N <- array(data = NA, dim = c(4,2))
  N[,] <- 500
  
  # Recruitment
  R <- 2
  
  # Survival
  S_C_B0_ps <- 0.99
  S_F_B0_ps <- 0.99
  S_M_B0_ps <- 0.99
  
  out <- list(
    init_N = N,
    R_B0 = R,
    S_C_B0_ps = S_C_B0_ps,
    S_F_B0_ps = S_F_B0_ps,
    S_M_B0_ps = S_M_B0_ps
  )
  return(out)
}

initial_values <- inits()

params = c(
  "N_tot",
  "N_f",
  "N_m",
  "N_c",
  "R",
  "lambda",
  "sd_afcount",
  "sd_amcount",
  "sd_R",
  "R_B0",
  "R_wt",
  "R_wm",
  "R_cg",
  "R_dd",
  "R_yr",
  "S_C_B0_ps",
  "S_Y_F_B0_ps",
  "S_Y_M_B0_ps",
  "test_ps",
  "S_wt",
  "S_wm",
  "S_cg",
  "S_dd",
  "sd_S_C",
  "sd_S_Y_F",
  "sd_S_Y_M",
  "S_C_yr",
  "S_Y_F_yr",
  "S_Y_M_yr",
  "survival_af",
  "survival_am",
  "survival_ca"
)

if(climate_variable != "null"){
  jags_data$clim <- ipm_data[[which(names(ipm_data) == climate_variable)]]
  model_file <- "s2//ipm.txt"
}


#Model==========================================================================
cat("This run was started at", as.character(Sys.time()), "\n")
tictoc::tic()

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

cat("From start to end of this model run,")
tictoc::toc()

mcmcplots::mcmcplot(rslt)


saveRDS(rslt, file = save_file)

