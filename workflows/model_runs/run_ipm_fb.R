# IPM -- Starkey Elk
# Kenneth Loonam
# December 2025

#Variables======================================================================
# JAGS control parameters
n_i <- 100000
n_a <- 10000
n_b <- 50000
n_c <- 3
n_t <- 50
run_date <- "06feb2026_"

### Set cougar covariate ### (un-comment one line)
# cgvar <- "cg_logistic_growth"; cgname <- "lgst"
# cgvar <- "cg_mean_estimate"; cgname <- "mean"
cgvar <- "cgnorm"; cgname <- "norm"
# cgvar <- "cg_odfw_estimate"; cgname <- "odfw"
# cgvar <- "cg_reconstruction"; cgname <- "rcns"

### Specify which demographic rate covariate set to use ###
# cov_set <- "full"; cvname <- "full"
cov_set <- "best"; cvname <- "best"
# cov_set <- "null"; cvname <- "null"

save_file <- paste0("results//fbipm_rslt_", run_date, cvname, "_", cgname, ".rds")
if(cvname == "null"){
  save_file <- paste0("results//fbipm_rslt_", run_date, cvname, ".rds")
}
# Specify the model
if(cgvar == "cgnorm"){
  model_file <- paste0("models//fbipm_", cov_set, "_", cgvar, ".txt")
}else{
  model_file <- paste0("models//fbipm_", cov_set, ".txt")
}
if(cvname == "null"){
  model_file <- "models//fbipm_null.txt"
}

# Loop dimension parameters
n_year <- 36
params = c(
  "N_f", "N_m", "N_c", "Ntot", "N_fy", "N_fp", "N_fo", "R", "Pbar",
  "sd_afcount", "sd_amcount", "cdrng",
  "P",  "P_B0",  "P_sd",            "PB_wm",  "PB_cg",  "PB_dd",  "P_Byr",
  "SN", "SN_B0", "SN_sd", "SNB_wt", "SNB_wm", "SNB_cg", "SNB_dd", "SN_Byr",
  "SC", "SC_B0", "SC_sd", "SCB_wt", "SCB_wm", "SCB_cg", "SCB_dd", "SC_Byr",
  "SF", "SF_B0", "SF_sd", "SFB_wt", "SFB_wm", "SFB_cg", "SFB_dd", "SF_Byr",
  "SM", "SM_B0", "SM_sd",                                         "SM_Byr",
  "H", "H_B0", "H_sd", "H_Byr", "Pf", "Pm"
)
#Environment====================================================================

require(tidyverse); require(rjags); require(mcmcplots); require(furrr)
ipm_data <- readRDS("data//ipm_data_03feb2026.rds")
prg_dt <- read.csv("data//pregnancy_rates.csv") %>%
  as_tibble() %>%
  mutate(preg = round(pregnancy_rate * n_observations)) %>%
  group_by(yr, age_class) %>%
  summarise(no = sum(n_observations), np = sum(preg)) %>%
  mutate(yr = yr - 1986) %>%
  mutate(age = case_when(
    age_class == "young" ~ 1,
    age_class == "prime" ~ 2,
    age_class == "old" ~ 3
  )) %>%
  select(yr, age, np, no) %>%
  as.matrix()
rec_dt_raw <- readRDS("data//recruitment_data.rds")
rec_dt <- tibble(
  yr = rec_dt_raw$constants$years - 1987,
  nca = rec_dt_raw$data$n_calf,
  naf = rec_dt_raw$data$n_cow
) %>% as.matrix()
cjs_dt_raw <- readRDS("data//cjs_data.rds")
sfyrdt <- readRDS("data//sfyr_informative_prior_information.rds")
scyrdt <- readRDS("data//scyr_informative_prior_information.rds")
#Data_prep======================================================================
keep_sf_yr <- ipm_data$s_cjs %>%
  as_tibble() %>% filter(age == 3 & sex == 1) %>% pull(year)
keep_sc_yr <- ipm_data$s_cjs %>%
  as_tibble() %>% filter(age == 2 & sex == 1) %>% pull(year)
keep_sm_yr <- ipm_data$s_cjs %>%
  as_tibble() %>% filter(age == 3 & sex == 2) %>% pull(year)

jags_data <- list(
  r.af       = cjs_dt_raw$data$r.af,
  r.am       = cjs_dt_raw$data$r.am,
  r.jf       = cjs_dt_raw$data$r.jf,
  r.jm       = cjs_dt_raw$data$r.jm,
  ma.af      = cjs_dt_raw$data$ma.af,
  ma.am      = cjs_dt_raw$data$ma.am,
  ma.jf      = cjs_dt_raw$data$ma.jf,
  ma.jm      = cjs_dt_raw$data$ma.jm,
  sf_yri     = as.numeric(1:ipm_data$n_year %in% keep_sf_yr),
  sm_yri     = as.numeric(1:ipm_data$n_year %in% keep_sm_yr),
  sc_yri     = as.numeric(1:ipm_data$n_year %in% keep_sc_yr),
  nt         = length(cjs_dt_raw$data$r.af) + 1,
  rec_dt     = rec_dt,
  nr         = nrow(rec_dt),
  prg_dt     = prg_dt,
  np         = nrow(prg_dt),
  n_a_mov    = ipm_data$n_a_mov,
  n_c_mov    = ipm_data$n_c_mov,
  n_year     = ipm_data$n_year,
  na         = ipm_data$na,
  n_har      = ipm_data$n_har,
  n_har_hack = ipm_data$n_har,
  min_ad     = ipm_data$min_ad,
  min_ca     = ipm_data$min_ca,
  af_count   = ipm_data$af_count,
  nn_fc      = ipm_data$nn_fc,
  am_count   = ipm_data$am_count,
  nn_mc      = ipm_data$nn_mc,
  har_i      = ipm_data$harvest_indicator,
  min_nf1    = ipm_data$min_nf1,
  est_nf1    = ipm_data$est_nf1,
  min_nm1    = ipm_data$min_nm1,
  est_nm1    = ipm_data$est_nm1,
  cd_mean    = ipm_data$cg_mean_estimate,
  cd_tau     = 1/ipm_data$cg_sd_estiamte^2,
  clim       = ipm_data$spei12,
  n_adj      = ipm_data$nelk
)
if(cgname %in% c("lgst", "odfw", "rcns", "mean")){
  jags_data$cdens <- ipm_data[[cgvar]]
}

#Initial Values=================================================================
inits <- function(chain){
  # Variable Assignment
  nyr = 36
  gsd = 0.5
  ppb0 = 3
  snb0 = 3
  scb0 = 3
  sfb0 = 3
  smb0 = 3
  
  # Starting abundance
  N <- array(data = NA, dim = c(ipm_data$na,2,ipm_data$n_year))
  N[,,1] <- 250
  H_Byr <- array(-2, dim = c(2,2,ipm_data$n_year))
  H_Byr[,,1] <- NA
  
  P_Byr <- rep(0, nyr)
  P_Byr[1] <- NA
  SN_Byr <- rep(0, nyr)
  SN_Byr[1] <- NA
  SC_Byr <- rep(0, nyr)
  SC_Byr[1] <- NA
  SF_Byr <- rep(0, nyr)
  SF_Byr[1] <- NA
  SM_Byr <- rep(0, nyr)
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
    P_B0 = c(1.26, 1.79, 0.55),
    P_sd = 0.54,
    P_Byr = P_Byr,
    SN_B0 = 0.2462668,
    SN_sd = 0.7573277,
    SN_Byr = SN_Byr,
    SC_B0 = 0.8352736,
    SC_sd = 0.6027709,
    SC_Byr = SC_Byr,
    SF_B0 = 2.588588,
    SF_sd = 0.637042,
    SF_Byr = SF_Byr,
    SM_B0 = 0.7176843,
    SM_sd = 0.2750366,
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

# tst_rslt <- run_jags_model(
#   ipm_inits = jgs_ivs[[1]],
#   mf = model_file,
#   ipm_dt = jags_data,
#   params = params
# )

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

mcmcplots::mcmcplot(rslt)


saveRDS(rslt, file = save_file)


