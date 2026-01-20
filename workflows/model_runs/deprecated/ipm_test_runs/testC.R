# IPM -- Starkey Elk
# Kenneth Loonam
# April 2021

#Variables======================================================================

# Specify the model
model_file <- "models//ipm//tests//21june2024_testC.txt"
save_file <- "results//ipm//tests//21june2024_rsltC.rds"

# Loop dimension parameters
n_year <- 34

# JAGS control parameters
n_i <- 25000
n_a <- 5000
n_b <- 10000
n_c <- 3
n_t <- 5


#Environment====================================================================

require(tidyverse); require(rjags); require(mcmcplots)
# load("data//elk_ipm_data.Rdata")
load("data//elk_ipm_data_14jun2024.Rdata")

#Data_prep======================================================================

jags_data <- list(
  # s_cjs      = ipm_data$s_cjs,
  # p_cjs      = p_cjs[-c(95:96),], # drop the rows with counts of zero
  # r_ratio    = ipm_data$r_ratio,
  n_sight_ca = ipm_data$n_sight_ca,
  n_sight_am = ipm_data$n_sight_am,
  n_sight_af = ipm_data$n_sight_af,
  n_a_mov    = ipm_data$n_ad_add - abs(ipm_data$n_ad_rem),
  n_c_mov    = ipm_data$n_ca_add - abs(ipm_data$n_ca_rem),
  n_year     = n_year,
  nn_ca      = nrow(ipm_data$n_sight_ca),
  nn_af      = nrow(ipm_data$n_sight_af),
  nn_am      = nrow(ipm_data$n_sight_am),
  nn_count   = nrow(ipm_data$p_cjs),
  # ns         = nrow(ipm_data$s_cjs),
  # nr         = nrow(ipm_data$r_ratio),
  n_har      = ipm_data$n_hnt,
  min_ad     = ipm_data$min_ad,
  min_ca     = ipm_data$min_ca,
  af_count   = ipm_data$n_f_p_count,
  nn_fc      = nrow(ipm_data$n_f_p_count),
  am_count   = ipm_data$n_m_p_count,
  nn_mc      = nrow(ipm_data$n_m_p_count),
  cdens      = ipm_data$puma_derived,
  n_adj      = ipm_data$af_density,
  min_n1     = ipm_data$min_n1,
  est_n1     = ipm_data$est_n1,
  clim       = ipm_data$palmer_index,
  y          = ipm_data$y,
  z          = ipm_data$z,
  f          = ipm_data$f,
  l          = ipm_data$l,
  calf       = ipm_data$calf,
  male       = ipm_data$male,
  female     = ipm_data$female,
  herd       = ipm_data$herd,
  n_ind      = ipm_data$n_ind,
  n_R_data   = nrow(ipm_data$r_data),
  r_data     = ipm_data$r_data
  # NF_ct      = ipm_data$cjs_n_female,
  # NM_ct      = ipm_data$cjs_n_male
)



inits <- function(){
  # Starting abundance
  N <- array(data = NA, dim = c(4,2))
  N[,] <- 500
  
  # Recruitment
  R <- 2
  
  # Survival
  S_C_B0 <- 3
  S_F_B0 <- 3
  S_M_B0 <- 3
  
  out <- list(
    init_N = N,
    R_B0 = R,
    S_C_B0 = S_C_B0,
    S_F_B0 = S_F_B0,
    S_M_B0 = S_M_B0
  )
  return(out)
}

initial_values <- inits()

params = c(
  "N_tot",
  "R",
  "S_f",
  "S_m",
  "S_c",
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
  "S_C_mean",
  # "sd_afcount",
  # "sd_amcount",
  "sd_R",
  "S_C_B0",
  "S_F_B0",
  "S_M_B0",
  "sd_S_C",
  "sd_S_M",
  "sd_S_F",
  "PF",
  "PM",
  "P",
  "exp_ct",
  "tau_ct",
  "Ne",
  "Te",
  "Naug"
)

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

# update(jgs_mdl, n.iter = n_b)

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

