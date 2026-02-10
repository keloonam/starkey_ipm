# IPM -- Starkey Elk
# Kenneth Loonam
# December 2025

#Variables======================================================================
# JAGS control parameters
n_i <- 25000
n_a <- 1000
n_b <- 10000
n_c <- 7
n_t <- 50
run_date <- "05feb2026_"
model_file <- "models//mscr_fb_test.txt"
save_file <- paste0("results//mscr_fb_test_rslt_", run_date, "null_fb", ".rds")

# Loop dimension parameters
n_year <- 36
params = c(
  "SC", "SC_B0", "SC_sd", "SC_Byr", "SC_adj", "Snca", "Shjf", "Shjm",
  "SF", "SF_B0", "SF_sd", "SF_Byr", "SF_adj", "Snaf", "Shaf",
  "SM", "SM_B0", "SM_sd", "SM_Byr", "SM_adj", "Snam", "Sham",
  "H", "H_B0", "H_sd", "Pf", "Pm"
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
  n_year     = ipm_data$n_year,
  har_i      = ipm_data$harvest_indicator
)

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


