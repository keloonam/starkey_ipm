nimble_code <- nimbleCode({
  ##### Survival #####
  # Random intercept setup -- Survival
  SC_B0_mean ~ T(dlogis(0, 1), -5, 5)
  # SF_mean ~ dbeta(3, 1.2)
  # SM_mean ~ dbeta(3, 1.2)
  SF_B0_mean ~ T(dlogis(1.3, .8), -5, 5)
  SM_B0_mean ~ T(dlogis(1.3, .8), -5, 5)
  SC_B0_sd ~ dunif(0, 5)
  SF_B0_sd ~ dunif(0, 5)
  SM_B0_sd ~ dunif(0, 5)
  
  # Random intercept setup -- Detection probability
  # PF_B0_mean ~ dlogis(0, 1)
  # PM_B0_mean ~ dlogis(0, 1)
  # PF_B0_sd ~ dunif(0, 5)
  # PM_B0_sd ~ dunif(0, 5)
  
  # Covariates on calf survival
  SC_Bvy ~ T(dlogis(0, 1), -5, 5) # veg proxy in year t
  SC_Bvm ~ T(dlogis(0, 1), -5, 5) # veg proxy in year t-1
  SC_Bpu ~ T(dlogis(0, 1), -5, 5) # puma index
  SC_Bdd ~ T(dlogis(0, 1), -5, 5) # elk density
  
  # Fixed effects for non-Main Study elk (some switch back and forth)
  # S__Bhe ~ dlogis(0, 1)
  # P__Bhe ~ dlogis(0, 1)
  
  # Annual survival assignment
  for(t in 2:n_year){
    SF_B0[t] ~ T(dnorm(SF_B0_mean, sd = SF_B0_sd), -5, 5)
    SM_B0[t] ~ T(dnorm(SM_B0_mean, sd = SM_B0_sd), -5, 5)
    SC_B0[t] ~ T(dnorm(SC_B0_mean, sd = SC_B0_sd), -5, 5)
    
    # PM_B0[t] ~ dnorm(PM_B0_mean, sd = PM_B0_sd)
    # PF_B0[t] ~ dnorm(PF_B0_mean, sd = PF_B0_sd)
    
    logit(SF[t]) <- SF_B0[t]
    logit(SM[t]) <- SM_B0[t]
    # logit(PF[t]) <- PF_B0[t]
    # logit(PM[t]) <- PM_B0[t]
    SC_Byr[t]    <- SC_B0[t] +
      SC_Bvy * veg[t] +
      SC_Bvm * veg[t-1] +
      SC_Bpu * pum[t] +
      SC_Bdd * elk[t]
    logit(SC[t]) <- SC_Byr[t]
  }
  
  ##### Survival #####
  # CJS - Feedgrounds
  for(i in 1:nsc){
    sc_cjs[i,1] ~ T(dnorm(SC[sc_cjs_t[i]], sc_cjs[i,2]), 0, 1)
  }
  for(i in 1:nsf){
    sf_cjs[i,1] ~ T(dnorm(SF[sf_cjs_t[i]], sf_cjs[i,2]), 0, 1)
  }
  for(i in 1:nsm){
    sm_cjs[i,1] ~ T(dnorm(SM[sm_cjs_t[i]], sm_cjs[i,2]), 0, 1)
  }
})

require(nimble); require(mcmcplots); require(dplyr)

full_data <- readRDS("data//the_ipm_data.rds")

load("data//elk_ipm_data_21apr2023.Rdata")

#Data===========================================================================

dtf <- list(
  sc_cjs = ipm_data$s_cjs[62:92,4:5],
  sf_cjs = ipm_data$s_cjs[1:30, 4:5],
  sm_cjs = ipm_data$s_cjs[31:61,4:5]
)

#Constants======================================================================
cnst <- list(
  n_year = full_data$n_year,
  veg = full_data$sep_pdi,
  pum = full_data$puma_composit,
  nsc = nrow(dtf$sc_cjs),
  nsf = nrow(dtf$sf_cjs),
  nsm = nrow(dtf$sm_cjs),
  elk = full_data$elk_density,
  sc_cjs_t = ipm_data$s_cjs[62:92,1],
  sf_cjs_t = ipm_data$s_cjs[1:30, 1],
  sm_cjs_t = ipm_data$s_cjs[31:61,1]
)

#Initial values=================================================================

inits <- readRDS("data//the_ipm_inits.rds")
inits <- inits[c(15:31)]

#Model setup====================================================================
mons <- c(
  "SC", "SF", "SM", "SC_B0_mean", "SC_B0_sd", "SC_Bvy", "SC_Bvm", 
  "SC_Bpu", "SC_Bdd", "SM_B0_mean", "SF_B0_mean", "SM_B0_sd",
  "SF_B0_sd"
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
