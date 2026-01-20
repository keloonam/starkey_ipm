# IPM -- Starkey Elk
# Kenneth Loonam
# April 2021

#Variables======================================================================

# Specify the model
model_file <- "models//ipm//ipm_first_submission_version.txt"
save_file <- "results//amiinsane.Rdata"

# Loop dimension parameters
n_year <- 34

# JAGS control parameters
n_i <- 500000
n_a <- 1000
n_b <- 100000
n_c <- 2
n_t <- 5

#Environment====================================================================

require(tidyverse); require(rjags); require(mcmcplots)
# load("data//elk_ipm_data.Rdata")
load("data//elk_ipm_data_14jun2024.Rdata")

#Data_prep======================================================================

jags_data <- list(
  s_cjs      = ipm_data$s_cjs,
  p_cjs      = ipm_data$p_cjs,
  r_ratio    = ipm_data$r_ratio,
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
  ns         = nrow(ipm_data$s_cjs),
  nr         = nrow(ipm_data$r_ratio),
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
  clim       = ipm_data$palmer_index
  # y          = ipm_data$y,
  # z          = ipm_data$z,
  # f          = ipm_data$f,
  # l          = ipm_data$l,
  # calf       = ipm_data$calf,
  # male       = ipm_data$male,
  # female     = ipm_data$female,
  # herd       = ipm_data$herd,
  # n_ind      = ipm_data$n_ind,
  # n_R_data   = nrow(ipm_data$r_data),
  # r_data     = ipm_data$r_data
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
    S_C_B0_ps = 0.9,
    S_F_B0_ps = 0.9,
    S_M_B0_ps = 0.9
  )
  return(out)
}

initial_values <- inits()

params = c(
  "N_tot",
  "R",
  "survival_af",
  "aurvival_am",
  "survival_ca",
  "N_f",
  "N_c",
  "N_m",
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
  # "sd_afcount",
  # "sd_amcount",
  "sd_R"
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

