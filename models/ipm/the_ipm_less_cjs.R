nimble_code <- nimbleCode({
  # IPM -- Built for the Starkey elk herd in Main Study
  # Kenneth Loonam
  # Parameters to estimate:
  # R: recruitment
  # S: survival (C, F, and M stand for calf, female, and male)
  # N: abundance (C, F, and M stand for calf, female, and male)
  
  # Estimated parameters are CAPITALIZED
  # Data parameters are lowercase
  # Index information is the exception, e.g. N_e
  
  # N_e = expected number
  # N_man = net change from management
  # N_S = N ____ surviving
  # N__est = estimates of N_ passed as data
    # indexed by observation (row), and columns = {(1)year (2)mean (3)sd}
  # r_dt = recruitment data
    # indexed by observation (row), and columns = {(1)year (2)calves (3)cows}
  
  #Process Model================================================================
  for(t in 2:n_year){
    ##### Recruitment #####
    # Expected recruitment for normal approximation
    NCe_mean[t] <- R[t] * NF[t]
    NCe_sd[t]   <- sqrt(R[t] * NF[t] * (1 - R[t]))
    
    # Recruitment process
    NC[t] ~ T(dnorm(NCe_mean[t], sd = NCe_sd[t]), NCmin[t], )
    
    ##### Survival #####
    # Expected survival for normal approximation
    ## Calf survival - expected
    NCSe_mean[t] <- (NC[t-1] + NCman[t]) * SC[t]
    NCSe_sd[t] <- sqrt(NCSe_mean[t] * (1 - SC[t]))
    ## Female survival - expected
    NFSe_mean[t] <- (NF[t-1] + NFman[t]) * SF[t]
    NFSe_sd[t] <- sqrt(NFSe_mean[t] * (1 - SF[t]))
    ## Male survival - expected
    NMSe_mean[t] <- (NM[t-1] + NMman[t]) * SM[t]
    NMSe_sd[t] <- sqrt(NMSe_mean[t] * (1 - SM[t]))
    
    # Survival process - make it to August 1 (start of hunting season)
    NCaug[t] ~ T(dnorm(NCSe_mean[t], sd = NCSe_sd[t]), NCmin[t], )
    NFaug[t] ~ T(dnorm(NFSe_mean[t], sd = NFSe_sd[t]), NFmin[t], )
    NMaug[t] ~ T(dnorm(NMSe_mean[t], sd = NMSe_sd[t]), NMmin[t], )
    
    # Harvest process - make it to November 1 (biological year start)
    NF[t] <- NFaug[t] + NCaug[t]/2 - NFhar[t]
    NM[t] <- NMaug[t] + NCaug[t]/2 - NMhar[t]
    
    ##### Abundance #####
    # Expected abundance ignoring management & harvest - for calculating lambda
    NCSe_lambda[t] <- NC[t-1] * SC[t]
    NFe_lambda[t] <- NF[t-1] * SF[t] + NCSe_lambda[t] / 2
    NMe_lambda[t] <- NM[t-1] * SM[t] + NCSe_lambda[t] / 2
    NCe_lambda[t] <- (NFe_lambda[t] + NCSe_lambda[t]) * R[t]
    N_lambda[t] <- NCe_lambda[t] + NMe_lambda[t] + NFe_lambda[t]
    LAMBDA[t] <- N_lambda[t] / Ntot[t-1]
    
    # Total abundance
    Ntot[t] <- NM[t] + NF[t] + NC[t]
  }
  #Priors and GLMs==============================================================
  ##### Recruitment #####
  # Random intercept setup
  R_B0_mean ~ dlogis(0, 1)
  R_B0_sd ~ dunif(0, 5)
  
  # Covariates on recruitment
  R_Bvy ~ dlogis(0, 1) # veg proxy in year t
  R_Bvm ~ dlogis(0, 1) # veg proxy in year t-1
  R_Bpu ~ dlogis(0, 1) # puma index
  R_Bdd ~ dlogis(0, 1) # elk density
  
  # Annual recruitment assignment
  for(t in 2:n_year){
    R_B0[t] ~ dnorm(R_B0_mean, sd = R_B0_sd)
    logit(R[t]) <- R_B0[t] + 
      R_Bvy * veg[t] +
      R_Bvm * veg[t-1] +
      R_Bpu * pum[t] +
      R_Bdd * elk[t]
  }
  
  ##### Survival #####
  # Random intercept setup -- Survival
  SC_B0_mean ~ dlogis(0, 1)
  SF_mean ~ dbeta(3, 1.2)
  SM_mean ~ dbeta(3, 1.2)
  SF_B0_mean <- logit(SF_mean)
  SM_B0_mean <- logit(SM_mean)
  SC_B0_sd ~ dunif(0, 5)
  SF_B0_sd ~ dunif(0, 5)
  SM_B0_sd ~ dunif(0, 5)
  
  # Random intercept setup -- Detection probability
  PF_B0_mean ~ dlogis(0, 1)
  PM_B0_mean ~ dlogis(0, 1)
  PF_B0_sd ~ dunif(0, 5)
  PM_B0_sd ~ dunif(0, 5)
  
  # Covariates on calf survival
  SC_Bvy ~ dlogis(0, 1) # veg proxy in year t
  SC_Bvm ~ dlogis(0, 1) # veg proxy in year t-1
  SC_Bpu ~ dlogis(0, 1) # puma index
  SC_Bdd ~ dlogis(0, 1) # elk density
  
  # Fixed effects for non-Main Study elk (some switch back and forth)
  S__Bhe ~ dlogis(0, 1)
  P__Bhe ~ dlogis(0, 1)
  
  # Annual survival assignment
  for(t in 2:n_year){
    SF_B0[t] ~ dnorm(SF_B0_mean, sd = SF_B0_sd)
    SM_B0[t] ~ dnorm(SM_B0_mean, sd = SM_B0_sd)
    SC_B0[t] ~ dnorm(SC_B0_mean, sd = SC_B0_sd)
    
    PM_B0[t] ~ dnorm(PM_B0_mean, sd = PM_B0_sd)
    PF_B0[t] ~ dnorm(PF_B0_mean, sd = PF_B0_sd)
    
    logit(SF[t]) <- SF_B0[t]
    logit(SM[t]) <- SM_B0[t]
    logit(PF[t]) <- PF_B0[t]
    logit(PM[t]) <- PM_B0[t]
    SC_Byr[t]    <- SC_B0[t] +
      SC_Bvy * veg[t] +
      SC_Bvm * veg[t-1] +
      SC_Bpu * pum[t] +
      SC_Bdd * elk[t]
    logit(SC[t]) <- SC_Byr[t]
  }
  
  ##### Starting Abundance #####
  NC[1] ~ T(dnorm(nc1e, sd = 100), nc1e_min, )
  NF[1] ~ T(dnorm(nf1e, sd = 100), nf1e_min, )
  NM[1] ~ T(dnorm(nm1e, sd = 100), nm1e_min, )
  
  #Observation Models===========================================================
  ##### Abundance #####
  for(i in 1:nNC){ 
    NC_est[i,2] ~ T(dnorm(NC[NC_est[i,1]], sd = NC_est[i,3]), 0, )
  }
  for(i in 1:nNF){
    NF_est[i,2] ~ T(dnorm(NF[NF_est[i,1]], sd = NF_est[i,3]), 0, )
  }
  for(i in 1:nNM){
    NM_est[i,2] ~ T(dnorm(NM[NM_est[i,1]], sd = NM_est[i,3]), 0, )
  }
  
  ##### Counts #####
  sd_afcount ~ dunif(0,50)
  tau_afcount <- 1/sd_afcount^2
  for(i in 1:nn_fc){
    af_count[i,4] ~ T(dnorm(NFaug[af_count[i,1]], tau_afcount), 0, )
  }
  
  sd_amcount ~ dunif(0,50)
  tau_amcount <- 1/sd_amcount^2
  for(i in 1:nn_mc){
    am_count[i,4] ~ T(dnorm(NMaug[am_count[i,1]], tau_amcount), 0, )
  }
  # for(t in 2:n_year){ # loop over number of observations - no P in year 1
  #   # prep normal approximation of binomial based on cjs p and NF
  #   NFct_sd[t] <- sqrt(PF[t] * NF[t] * (1 - PF[t]))
  #   NFct_ex[t] <- PF[t] * NF[t]
  #   NF_ct[t] ~ dnorm(NFct_ex[t], sd = NFct_sd[t])
  # }
  # for(t in 2:n_year){ # loop over number of observations - no P in year 1
  #   # prep normal approximation of binomial based on cjs p and NF
  #   NMct_sd[t] <- sqrt(PM[t] * NM[t] * (1 - PM[t]))
  #   NMct_ex[t] <- PM[t] * NM[t]
  #   NM_ct[t] ~ dnorm(NMct_ex[t], sd = NMct_sd[t])
  # }
  
  ##### Recruitment #####
  for(i in 1:nR){ # loop over number of observations of recruitment
    r_dt[i,2] ~ dbinom(R[r_dt[i,1]], r_dt[i,3])
  }
  
  ##### Survival #####
  # CJS - Feedgrounds
  # S[1,1:3,1:2] <- 0
  # S[2:n_year,1,1:2] <- 0
  for(t in 2:n_year){
    S[t,2,1:2] <- SC[t]
    S[t,3,1] <- SF[t]
    S[t,3,2] <- SM[t]
  }
  for(i in 1:ns){
    s_cjs[i,4] ~ T(dnorm(S[s_cjs[i,1], s_cjs[i,2], s_cjs[i,3]], s_cjs[i,5]),0,1)
  }
  
  
  # for(i in 1:nind){
  #   for(t in (f[i]+1):l[i]){
  #     # Setting Probabilities
  #     logit(p[i,t]) <- PF_B0[t] * female[i,t-1] + 
  #       PM_B0[t] * male[i,t-1] + 
  #       P__Bhe * herd[i,t-1]
  #     logit(s[i,t]) <- SF_B0[t] * female[i,t-1] +
  #       SM_B0[t] * male[i,t-1] +
  #       SC_Byr[t] * calf[i,t-1] +
  #       S__Bhe * herd[i,t-1]
  #     
  #     # State Process ------- alive >>> z = 1
  #     z[i,t] ~ dbern(s[i,t] * z[i,t-1])
  #     # Observation Process ------- seen >>> y = 1
  #     y[i,t] ~ dbern(p[i,t] * z[i,t])
  #   }
  # }
})






