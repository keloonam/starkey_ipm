# IPM -- Starkey Elk
# Kenneth Loonam
# December 2025

#Variables======================================================================
# JAGS control parameters
n_i <- 100000
n_a <- 100000
n_b <- 100000
n_c <- 3
n_t <- 50
run_date <- "05jan2025_"
### Set cougar covariate ### (un-comment one line)
# cgvar <- "cgnorm"; cgname <- "norm"
# cgvar <- "cgunif"; cgname <- "unif"
# cgvar <- "cg_logistic_growth"; cgname <- "lgst"
# cgvar <- "cg_odfw_estimate"; cgname <- "odfw"
# cgvar <- "cg_reconstruction"; cgname <- "rcns"
# cgvar <- "cg_mean_estimate"; cgname <- "mean"

### Specify which demographic rate covariate set to use ###
cov_set <- "all"; cvname <- "all"
# cov_set <- "original_set"; cvname <- "org"
# cov_set <- "bst"; cvname <- "bst"
# cov_set <- "test"; cvname <- "test"

save_file <- paste0("results//ipm_rslt_", run_date, cvname, "_", cgname, ".rds")
# Specify the model
if(cgvar %in% c("cgnorm", "cgunif")){
  model_file <- paste0("models//ipm_", cov_set, "_", cgvar, ".txt")
}else{
  model_file <- paste0("models//ipm_", cov_set, ".txt")
}

# Loop dimension parameters
n_year <- 36
params = c(
  "N_f", "N_m", "N_c", "Ntot", "N_fy", "N_fp", "N_fo", "R", "Pbar",
  "sd_afcount", "sd_amcount", "cdrng",
  "P",  "P_B0",  "P_sd",            "PB_wm",  "PB_cg",  "PB_dd",  "P_Byr",
  "SN", "SN_B0", "SN_sd", "SNB_wt",           "SNB_cg", "SNB_dd", "SN_Byr",
  "SC", "SC_B0", "SC_sd", "SCB_wt", "SCB_wm", "SCB_cg", "SCB_dd", "SC_Byr",
  "SF", "SF_B0", "SF_sd", "SFB_wt", "SFB_wm", "SFB_cg", "SFB_dd", "SF_Byr",
  "SM", "SM_B0", "SM_sd", "SMB_wt", "SMB_wm", "SMB_cg", "SMB_dd", "SM_Byr",
  "H", "H_B0", "H_sd"
)

#Environment====================================================================

require(tidyverse); require(rjags); require(mcmcplots); require(furrr)
ipm_data <- readRDS("data//ipm_data_30dec2025.rds")

#Data_prep======================================================================
jags_data <- list(
  s_cjs      = ipm_data$s_cjs,
  h_cjs      = ipm_data$h_cjs,
  r_ratio    = ipm_data$r_ratio,
  p_ratio    = ipm_data$p_ratio,
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
  nh         = ipm_data$nh,
  np         = ipm_data$np,
  na         = ipm_data$na,
  n_har      = ipm_data$n_har,
  n_har_hack = ipm_data$n_har,
  min_ad     = ipm_data$min_ad,
  min_ca     = ipm_data$min_ca,
  af_count   = ipm_data$af_count,
  nn_fc      = ipm_data$nn_fc,
  am_count   = ipm_data$am_count,
  nn_mc      = ipm_data$nn_mc,
  # cdens      = ipm_data$cg_logistic_growth,
  cd_mean    = ipm_data$cg_mean_estimate,
  cd_tau     = 1/ipm_data$cg_sd_estiamte^2,
  cd_min     = ipm_data$cg_min_estimate,
  cd_max     = ipm_data$cg_max_estimate,
  n_adj      = ipm_data$nelk,
  har_i      = ipm_data$harvest_indicator,
  clim       = ipm_data$spei12,
  min_n1     = ipm_data$min_n1,
  est_n1     = ipm_data$est_n1
)

if(cgname %in% c("lgst", "odfw", "rcns", "mean")){
  jags_data$cdens <- ipm_data[[cgvar]]
}

#Initial Values=================================================================
inits <- function(chain){
  # Variable Assignment
  nyr = 36
  gsd = 0.5
  ppb0 = 1
  snb0 = 1
  scb0 = 1
  sfb0 = 1
  smb0 = 1
  
  # Starting abundance
  N <- array(data = NA, dim = c(ipm_data$na,2,ipm_data$n_year))
  N[,,1] <- 1000
  H_Byr <- array(-1, dim = c(2,2,ipm_data$n_year))
  H_Byr[,,1] <- NA
  
  P_B0 <- ppb0
  P_sd <- gsd
  P_Byr <- rep(ppb0, nyr)
  P_Byr[1] <- NA
  SN_B0 <- snb0
  SN_sd <- gsd
  SN_Byr <- rep(snb0, nyr)
  SN_Byr[1] <- NA
  SC_B0 <- scb0
  SC_sd <- gsd
  SC_Byr <- rep(scb0, nyr)
  SC_Byr[1] <- NA
  SF_B0 <- sfb0
  SF_sd <- gsd
  SF_Byr <- rep(sfb0, nyr)
  SF_Byr[1] <- NA
  SM_B0 <- smb0
  SM_sd <- gsd
  SM_Byr <- rep(smb0, nyr)
  SM_Byr[1] <- NA
  
  name_options <- c(
    "base::Wichmann-Hill",
    "base::Marsaglia-Multicarry",
    "base::Super-Duper",
    "base::Mersenne-Twister"
  )
  
  out <- list(
    N = N,
    H_Byr = H_Byr,
    P_B0 = rep(P_B0, 3),
    P_sd = P_sd,
    P_Byr = P_Byr,
    SN_B0 = SN_B0,
    SN_sd = SN_sd,
    SN_Byr = SN_Byr,
    SC_B0 = SC_B0,
    SC_sd = SC_sd,
    SC_Byr = SC_Byr,
    SF_B0 = SF_B0,
    SF_sd = SF_sd,
    SF_Byr = SF_Byr,
    SM_B0 = SM_B0,
    SM_sd = SM_sd,
    SM_Byr = SM_Byr,
    .RNG.name = name_options[chain%%length(name_options)+1],
    .RNG.seed = chain
  )
  return(out)
}
#Model Run Function=============================================================
run_jags_model <- function(
    ipm_inits,
    mf, 
    ipm_dt, 
    params, 
    nc = 1, 
    na = 1000, 
    nb = 1000, 
    ni = 1000, 
    nt = 1){
  
  jgs_mdl <- jags.model(
    file = mf,
    data = ipm_dt,
    inits = ipm_inits,
    n.chains = nc,
    n.adapt = na
  )
  
  update(jgs_mdl, n.iter = nb)
  
  rslt <- coda.samples(
    jgs_mdl,
    variable.names = params,
    n.iter = ni,
    thin = nt
  )
  return(rslt)
}
# Run Model=====================================================================

jgs_ivs <- map(1:n_c, inits)

plan(multisession, workers = n_c)
rslt <- future_map(
  .x = jgs_ivs, 
  .f = run_jags_model, 
  .options = furrr_options(seed = FALSE),
  mf = model_file,
  ipm_dt = jags_data,
  params = params,
  na = n_a,
  nb = n_b,
  ni = n_i,
  nt = n_t
)

# mcmcplots::mcmcplot(rslt)


saveRDS(rslt, file = save_file)


