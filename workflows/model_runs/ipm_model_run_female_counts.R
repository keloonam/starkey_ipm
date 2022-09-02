# IPM -- Starkey Elk
# Kenneth Loonam
# April 2021

#Variables======================================================================

# Specify the model
model_file <- "models//ipm//ipm.txt"
save_file <- "results//ipm_result_01sep2022.Rdata"

# Loop dimension parameters
n_year <- 34

# JAGS control parameters
n_i <- 5000
n_a <- 1000
n_b <- 5000
n_c <- 3
n_t <- 100

#Environment====================================================================

require(tidyverse); require(rjags); require(mcmcplots)
# load("data//elk_ipm_data.Rdata")
load("data//elk_ipm_data_02sep2022.Rdata")

#Data_prep======================================================================

jags_data <- list(
  s_cjs = ipm_data$s_cjs,
  # r_ratio = ipm_data$r_ratio,
  n_sight_ca = ipm_data$n_sight_ca,
  n_sight_am = ipm_data$n_sight_am,
  n_sight_af = ipm_data$n_sight_af,
  Na_mov = ipm_data$n_ad_add - abs(ipm_data$n_ad_rem),
  Nc_mov = ipm_data$n_ca_add - abs(ipm_data$n_ca_rem),
  nyr = n_year,
  nn_ca = nrow(ipm_data$n_sight_ca),
  nn_af = nrow(ipm_data$n_sight_af),
  nn_am = nrow(ipm_data$n_sight_am),
  ns = nrow(ipm_data$s_cjs),
  # nr = nrow(ipm_data$r_ratio),
  Nhar = ipm_data$n_hnt,
  Na_obs = ipm_data$min_ad,
  Nc_obs = ipm_data$min_ca,
  # st = ipm_data$summer_temp,
  # wt = ipm_data$winter_temp,
  # sp = ipm_data$summer_precip,
  # wp = ipm_data$winter_precip,
  # est_mean_n = 450,
  # est_sd_n = 126,
  af_count = ipm_data$n_f_p_count,
  nn_fc = nrow(ipm_data$n_f_p_count),
  r_NfNc = ipm_data$ratio_counts,
  n_rc = nrow(ipm_data$ratio_counts)
)

inits <- function(){
  # Starting abundance
  N <- array(data = NA, dim = c(4,2))
  N[,] <- 500
  
  # Recruitment
  Rb0 <- 10
  
  # Survival
  Sc0 <- 4
  Sf0 <- 4
  Sm0 <- 4
  
  out <- list(
    init_N = N,
    Rb0    = Rb0,
    Sc0    = Sc0,
    Sf0    = Sf0,
    Sm0    = Sm0
  )
  return(out)
}

initial_values <- inits()

params = c(
  "Nt",
  "Nf",
  "Nm",
  "Nc",
  "Sf",
  "Sc",
  "Sm",
  "R",
  "Rdd"
)

#Model==========================================================================

jgs_mdl <- jags.model(
  file     = model_file,
  data     = jags_data,
  inits    = inits,
  n.chains = n_c,
  n.adapt  = n_a
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
