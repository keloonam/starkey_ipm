# IPM -- Starkey Elk
# Kenneth Loonam
# April 2021

#Variables======================================================================

# Specify the model
model_file <- "models//ipm//ipm_elk_null_dec2021.txt"

# Loop dimension parameters
n_year <- 34

# Priors for logistically transformed parameters
logit_tau_prior <- 0.5
logit_mean_prior <- 0

# Priors for starting abundance
n1_calf <- 200
n1_f <- 400
n1_m <- 100
n1_sd <- 500

# JAGS control parameters
n_i <- 6000
n_b <- 1000
n_c <- 3
n_t <- 10

#Environment====================================================================

require(tidyverse); require(rjags); require(mcmcplots)
# load("data//elk_ipm_data.Rdata")
load("data//elk_ipm_data.Rdata")

#Temporary_section_bypassing_ipm_data_prep_workflow=============================

# Fix columns for new indexing and remove "bad" rows in survival data
censor_rows <- c(9,11,12,13,27,29,31,42,62,66,70,71,73,74,75,76,89,91,93)
ipm_data$s_cjs[,2] <- ipm_data$s_cjs[,2] - 1
ipm_data$s_cjs[,1] <- ipm_data$s_cjs[,1] + 1
ipm_data$s_cjs <- ipm_data$s_cjs[-censor_rows,]

# Build 3 dim n_hnt array [age,sex,year]
load("data//elk_harvest_data.Rdata")
hntdat$age[hntdat$age == 3] <- 2
n_hnt <- array(data = 0, dim = c(2,2,n_year))
hntdat <- as.matrix(hntdat)
for(i in 1:nrow(hntdat)){
  n_hnt[hntdat[i,2], hntdat[i,3], hntdat[i,1]] <- n_hnt[hntdat[i,2], hntdat[i,3], hntdat[i,1]] + hntdat[i,4]
}

# Build 3 dim n_mov array [age,sex,year]
mov_dat <- read_csv("data//mov_data_handle_summaries.csv") %>%
  as.matrix()
mov_dat[,1] <- mov_dat[,1] - 1987
n_mov <- array(0, dim = c(2,2,n_year))
for(i in 1:nrow(mov_dat)){
  n_mov[1,1,mov_dat[i,1]] <- mov_dat[i,5]
  n_mov[1,2,mov_dat[i,1]] <- mov_dat[i,4]
  n_mov[2,1,mov_dat[i,1]] <- mov_dat[i,2]
  n_mov[2,2,mov_dat[i,1]] <- mov_dat[i,3]
}

# Build 3 dim min_n array [age,sex,year]
min_dat_han <- read_csv("data//min_n_handle_summaries.csv") %>%
  as.matrix()
min_dat_han[,1] <- min_dat_han[,1] - 1987
n_min <- array(0, dim = c(2,2,n_year))
for(i in 1:nrow(min_dat_han)){
  n_min[1,1,min_dat_han[i,1]] <- floor(min_dat_han[i,4] / 2)
  n_min[1,2,min_dat_han[i,1]] <- floor(min_dat_han[i,4] / 2)
  n_min[2,1,min_dat_han[i,1]] <- min_dat_han[i,2]
  n_min[2,2,min_dat_han[i,1]] <- min_dat_han[i,3]
}


#Data_prep======================================================================

# Priors on recruitment -- only 2.5+ females breed
r_prior <- array(
  data = NA,
  dim = c(2, 2, 2)
)
r_prior[2,1,1] <- logit_mean_prior
r_prior[2,1,2] <- logit_tau_prior

# Priors on survival -- same prior for sex and age classes 
s_prior <- array(
  data = NA,
  dim = c(2,2,2)
)
s_prior[1:2,,1] <- logit_mean_prior
s_prior[1:2,,2] <- logit_tau_prior

# Priors on starting abundance -- same variance on all, initial am and af pops
n1_prior <- array(
  data = NA,
  dim = c(2,2,2)
)
n1_prior[1,1:2,1] <- n1_calf
n1_prior[2,1,1] <- n1_f
n1_prior[2,2,1] <- n1_m
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
  # mindat = ipm_data$min_n_live,
  # n_rem = ipm_data$net_remove,
  # n_mindat = nrow(ipm_data$min_n_live)
  n_min = n_min,
  n_mov = n_mov,
  n_hnt = n_hnt
)

# inits <- function(){
#   list(
#     N = array(data = rpois(3*2*n_year, 500), dim = c(3,2,n_year))
#   )
# }

params = c(
  # "lambda",
  "S_C_B0",
  "sd_S_C",
  "S_A_F_B0",
  "sd_S_AF",
  "S_A_M_B0",
  "sd_S_AM",
  "R_B0",
  "sd_R",
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
  n.thin = n_t
)

rslt

mcmcplots::mcmcplot(rslt)


save(rslt, file = "results//ipm_result.Rdata")

plot(rslt$BUGSoutput$summary[grep("N_tot", dimnames(rslt$BUGSoutput$summary)[[1]]),1], type = "l")
plot(rslt$BUGSoutput$summary[grep("recruitment", dimnames(rslt$BUGSoutput$summary)[[1]]),1], type = "l")
plot(rslt$BUGSoutput$summary[grep("survival_af", dimnames(rslt$BUGSoutput$summary)[[1]]),1], type = "l")


dim(jags_data$s_cjs)
dim(jags_data$r_ratio)
jags_data$s_cjs

