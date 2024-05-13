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
  
  
  #Process Model================================================================
  for(t in 2:n_year){
    ##### Recruitment #####
    # Expected recruitment for normal approximation
    NCe_mean[t] <- R[t] * NF[t]
    NCe_sd[t]   <- sqrt(R[t] * NF[t] * (1 - R[t]))
    
    # Recruitment process
    NC[t] ~ T(dnorm(NCe_mean[t], sd = NCe_sd[t]), NC_min[t], )
    
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
    NF[t] <- NF_aug[t] + NC_aug[t]/2 - NFhar[t]
    NM[t] <- NM_aug[t] + NC_aug[t]/2 - NMhar[t]
    
    ##### Abundance #####
    # Expected abundance ignoring management & harvest - for calculating lambda
    NCSe_lambda[t] <- NC[t-1] * SC[t]
    NFe_lambda[t] <- NF[t-1] * SF[t] + NCSe_lambda[t] / 2
    NMe_lambda[t] <- NM[t-1] * SM[t] + NCSe+lambda[t] / 2
    NCe_lambda[t] <- (NFe_lambda[t] + NCSe_lambda[t]) * R[t]
    N_lambda[t] <- NCe_lambda[t] + NMe_lambda[t] + NFe_lambda[t]
    LAMBDA <- N_lambda[t] / N_tot[t-1]
    
    # Total abundance
    N_tot[t] <- NM[t] + NF[t] + NC[t]
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
    R_B0[t] ~ dnorm(R_B0_mean, sd = R_B0_shape)
    logit(R[t]) <- R_B0 + 
      R_Bvy * veg[t] +
      R_Bvm * veg[t-1] +
      R_Bpu * pum[t] +
      R_Bdd * elk[t]
  }
  
  ##### Survival #####
  # Random intercept setup
  SC_B0_mean ~ dlogis(0, 1)
  SF_mean ~ dbeta(3, 1.2)
  SM_mean ~ dbeta(3, 1.2)
  SF_B0_mean <- logit(SF_mean)
  SM_B0_mean <- logit(SM_mean)
  SC_B0_sd ~ dunif(0, 5)
  SF_B0_sd ~ dunif(0, 5)
  SM_B0_sd ~ dunif(0, 5)
  
  # Covariates on calf survival
  SC_Bvy ~ dlogis(0, 1) # veg proxy in year t
  SC_Bvm ~ dlogis(0, 1) # veg proxy in year t-1
  SC_Bpu ~ dlogis(0, 1) # puma index
  SC_Bdd ~ dlogis(0, 1) # elk density
  
  # Annual survival assignment
  for(t in 2:n_year){
    SF_B0[t] ~ dnorm(SF_B0_mean, SF_B0_sd)
    SM_B0[t] ~ dnorm(SM_B0_mean, SM_B0_sd)
    SC_B0[t] ~ dnorm(SC_B0_mean, SC_B0_sd)
    logit(SF[t]) <- SF_B0[t]
    logit(SM[t]) <- SM_B0[t]
    logit(SC[t]) <- SC_B0[t] +
      SC_Bvy * veg[t] +
      SC_Bvm * veg[t-1] +
      SC_Bpu * pum[t] +
      SC_Bdd * elk[t]
  }
  
  ##### Starting Abundance #####
  NC[1] T(dnorm(nc1e, sd = 100), nc1e_min, )
  NF[1] T(dnorm(nf1e, sd = 100), nf1e_min, )
  NM[1] T(dnorm(nm1e, sd = 100), nm1e_min, )
  
  #Observation Models===========================================================
  ##### Abundance #####
  
})






