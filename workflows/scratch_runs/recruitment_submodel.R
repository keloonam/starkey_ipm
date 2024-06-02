require(nimble); require(mcmcplots); require(dplyr)

full_data <- readRDS("data//the_ipm_data.rds")

load("data//elk_ipm_data_21apr2023.Rdata")

nimble_code <- nimbleCode({
  #Priors and GLMs==============================================================
  ##### Recruitment #####
  # Random intercept setup
  R_B0_mean ~ T(dlogis(0, 1), -5, 5)
  R_B0_sd ~ dunif(0, 5)
  
  # Covariates on recruitment
  R_Bvy ~ T(dlogis(0, 1), -5, 5) # veg proxy in year t
  R_Bvm ~ T(dlogis(0, 1), -5, 5) # veg proxy in year t-1
  R_Bpu ~ T(dlogis(0, 1), -5, 5) # puma index
  R_Bdd ~ T(dlogis(0, 1), -5, 5) # elk density
  
  # Annual recruitment assignment
  for(t in 2:n_year){
    R_B0[t] ~ T(dnorm(R_B0_mean, sd = R_B0_sd), -5, 5)
    logit(R[t]) <- R_B0[t] + 
      R_Bvy * veg[t] +
      R_Bvm * veg[t-1] +
      R_Bpu * pum[t] +
      R_Bdd * elk[t]
  }
  ##### Recruitment #####
  for(i in 1:nR){ # loop over number of observations of recruitment
    r_dt[i,1] ~ dbinom(R[r_dt_t[i]], r_dt[i,2])
  }
})
full_data <- readRDS("data//the_ipm_data.rds")
load("data//elk_ipm_data_21apr2023.Rdata")

dtf <- list(
  r_dt = full_data$r_dt[,2:3]
)

#Constants======================================================================
cnst <- list(
  n_year = full_data$n_year,
  nR = full_data$nR,
  veg = full_data$sep_pdi,
  pum = full_data$puma_composit,
  elk = full_data$elk_density,
  r_dt_t = full_data$r_dt[,1]
)

#Initial values=================================================================

inits <- readRDS("data//the_ipm_inits.rds")
inits <- inits[c(1,8,9,10,11,12,13,14)]
#Model setup====================================================================
mons <- c(
  "R", "R_Bvy", "R_Bvm", "R_B0_mean", 
  "R_Bpu", "R_Bdd", "R_B0_sd"
)
# "PF", "PM",

# nimbleMCMC(code = model_code, data = dtf, constants = cnst, inits = inits,
#            monitors = mons, niter = 1000, nburnin = 100, nchains = 1)

ipm <- nimbleModel(
  code = nimble_code,
  constants = cnst,
  data = dtf,
  inits = inits
)

ipm_conf <- configureMCMC(
  model = ipm, 
  monitors = mons,
  onlySlice = T
)

ipm_mcmc <- buildMCMC(ipm_conf)
c_ipm_mcmc <- compileNimble(ipm_mcmc, ipm)
# c_ipm_mcmc$ipm_mcmc$run(10000)
# x <- c_ipm_mcmc$ipm_mcmc$mvSamples %>% as.matrix()

rslt <- runMCMC(
  c_ipm_mcmc$ipm_mcmc, 
  niter = 25000, 
  nburnin = 0, 
  nchains = 3, 
  inits = inits,
  thin = 10
)
mcmcplot(rslt)
