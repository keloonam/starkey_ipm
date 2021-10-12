# IPM -- Starkey Elk
# Kenneth Loonam
# April 2021

#Variables======================================================================

# Specify the model
model_file <- "models//ipm//stky_elk_ipm_norm_approx.txt"

# Loop dimension parameters
n_year <- 34

# Priors for logistically transformed parameters
logit_tau_prior <- 0.5
logit_mean_prior <- 0

# Priors for starting abundance
n1_calf <- 100
n1_f <- 500
n1_m <- 500
n1_sa_mf <- 500
n1_sd <- 500

# JAGS control parameters
n_i <- 2500
n_b <- 1000
n_c <- 3
n_t <- 10

#Environment====================================================================

require(tidyverse); require(rjags); require(mcmcplots)
load("data//elk_ipm_data.Rdata")

#Data_prep======================================================================

# Priors on recruitment -- only 2.5+ females breed
r_prior <- array(
  data = NA,
  dim = c(3, 2, 2)
)
r_prior[3,1,1] <- logit_mean_prior
r_prior[3,1,2] <- logit_tau_prior

# Priors on survival -- same prior for sex and age classes 
s_prior <- array(
  data = NA,
  dim = c(3,2,2)
)
s_prior[2:3,,1] <- logit_mean_prior
s_prior[2:3,,2] <- logit_tau_prior

# Priors on starting abundance -- same variance on all, initial am and af pops
n1_prior <- array(
  data = NA,
  dim = c(3,2,2)
)
n1_prior[1,1:2,1] <- n1_calf
n1_prior[2,1:2,1] <- n1_sa_mf
n1_prior[3,1,1] <- n1_f
n1_prior[3,2,1] <- n1_m
n1_prior[,,2] <- n1_sd

# Build data for jags

jags_data <- list(
  s_cjs = ipm_data$s_cjs,
  r_ratio = ipm_data$r_ratio,
  n_sight_ca = ipm_data$n_sight_ca,
  n_sight_am = ipm_data$n_sight_am,
  n_sight_af = ipm_data$n_sight_af,
  n_year = n_year,
  r_prior = r_prior,
  s_prior = s_prior,
  n1_pr = n1_prior,
  nn_ca = nrow(ipm_data$n_sight_ca),
  nn_af = nrow(ipm_data$n_sight_af),
  nn_am = nrow(ipm_data$n_sight_am),
  ns = nrow(ipm_data$s_cjs),
  nr = nrow(ipm_data$r_ratio),
  mindat = ipm_data$min_n_live,
  n_rem = ipm_data$net_remove,
  n_mindat = nrow(ipm_data$min_n_live)
)

inits <- function(){
  list(
    N = array(data = rpois(3*2*n_year, 500), dim = c(3,2,n_year))
  )
}

params = c(
  # "lambda",
  "N_tot",
  "N_ca",
  "N_af",
  "N_am",
  "survival_af",
  "survival_am",
  "survival_ca",
  "recruitment"
)

#Model==========================================================================

jgs_mdl <- jags.model(
  file = model_file,
  data = jags_data,
  # inits = inits,
  n.chains = n_c
)

update(jgs_mdl, n.iter = n_b)

rslt <- coda.samples(
  jgs_mdl,
  variable.names = params,
  n.iter = n_i,
  thin = n_t
)

summary(rslt)


save(rslt, file = "results//ipm_result.Rdata")

plot(rslt$BUGSoutput$summary[grep("N_af", dimnames(rslt$BUGSoutput$summary)[[1]]),1], type = "l")
plot(rslt$BUGSoutput$summary[grep("recruitment", dimnames(rslt$BUGSoutput$summary)[[1]]),1], type = "l")
plot(rslt$BUGSoutput$summary[grep("survival_ca", dimnames(rslt$BUGSoutput$summary)[[1]]),1], type = "l")
