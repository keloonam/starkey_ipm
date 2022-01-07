# IPM -- Starkey Elk
# Kenneth Loonam
# April 2021

#Variables======================================================================

# Specify the model
model_file <- "models//ipm//ipm_elk_null_07jan2022.txt"

# Loop dimension parameters
n_year <- 33

# Priors for starting abundance
n1_calf <- 200
n1_f <- 400
n1_m <- 100
n1_sd <- 500

# JAGS control parameters
n_i <- 60000
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
  # n_ad_rem = abs(ipm_data$n_ad_rem),
  # n_ca_rem = abs(ipm_data$n_ca_rem),
  n_year = n_year,
  nn_ca = nrow(ipm_data$n_sight_ca),
  nn_af = nrow(ipm_data$n_sight_af),
  nn_am = nrow(ipm_data$n_sight_am),
  ns = nrow(ipm_data$s_cjs),
  nr = nrow(ipm_data$r_ratio)
  # n_hnt = ipm_data$n_hnt
)

inits <- function(){
  N <- array(data = NA, dim = c(4,2))
  N[,] <- 5000
  out <- list(init_N = N)
  return(out)
}

params = c(
  "N_tot",
  "survival_ca",
  "survival_af",
  "R"
)

#Model==========================================================================

# jgs_mdl <- jags.model(
#   file = model_file,
#   data = jags_data,
#   inits = inits,
#   n.chains = n_c
# )
# 
# update(jgs_mdl, n.iter = n_b)
# 
# rslt <- coda.samples(
#   jgs_mdl,
#   variable.names = params,
#   n.iter = n_i,
#   thin = n_t
# )

rslt <- R2jags::jags(
  data = jags_data,
  parameters.to.save = params,
  model.file = model_file,
  n.chains = n_c,
  n.iter = n_i,
  n.burnin = n_b,
  n.thin = n_t,
  inits = inits
)

rslt

mcmcplots::mcmcplot(rslt)


save(rslt, file = "results//ipm_result.Rdata")

plot(rslt$BUGSoutput$summary[grep("N_tot", dimnames(rslt$BUGSoutput$summary)[[1]]),1], type = "l")
plot(rslt$BUGSoutput$summary[grep("R", dimnames(rslt$BUGSoutput$summary)[[1]]),1], type = "l")
plot(rslt$BUGSoutput$summary[grep("survival_af", dimnames(rslt$BUGSoutput$summary)[[1]]),1], type = "l")
plot(rslt$BUGSoutput$summary[grep("survival_ca", dimnames(rslt$BUGSoutput$summary)[[1]]),1], type = "l")


jags_data$n_fg_c[,4] / jags_data$n_fg_f[,4]
jags_data$r_ratio[,4]
jags_data$s_cjs

