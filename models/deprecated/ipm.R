code <- nimbleCode({
  # IPM -- Built for the Starkey elk population
  # Kenneth Loonam
  # Clean re-write as of December 2025
  
  # Demographic parameters:
    # R: recruitment (PP*SC)
    # PP: pregnancy probability
    # SR: 6-month neonate survival (survive to recruitment)
    # SN: natural survival
    # SH: harvest survival
    # NS: abundance surviving the year (August)
    # NR: abundance recruited (6-month-olds from each age class)
    # NH: number harvested (known exactly)
    # NM: number moved during management activities (winter)
    # N: abundance (November)
  
  #####################################
  ########### Process Model ###########
  #####################################
  #======= Recruitment Calculations =======#
  for(t in 2:nt){ # t = time (years)
    for(a in 1:na){ # a = age (years - 0.5)
      R[a,t] <- PP[a,t-1] * SR[t]
      NRm[a,t] <- R[a,t] * N[a,1,t]
      NRv[a,t] <- NRm[a,t] * (1 - R[a,t])
    } # a
    for(x in 1:2){ # x = sex (1=female; 2=male)
      Nm[1,x,t] <- sum(NRm[1:na,t]) / 2
      Nv[1,x,t] <- sum(NRv[1:na,t]) / 2
      N[1,x,t] ~ dnorm(Nm[1,x,t], 1/Nv[1,x,t])
    }
  } # t
  #======= Survival Calculations =======#
  for(t in 2:nt){ 
    for(x in 1:2){
      #----- Natural Survival -----#
      # Survive from November of t-1 to August of t
      for(a in 2:(na-1)){
        NSm[a,x,t] <- (N[a-1,x,t-1] + NM[a-1,x,t-1]) * SN[a-1,x,t-1]
        NSv[a,x,t] <- NSm[a,x,t] * (1 - SN[a-1,x,t-1])
      } # a
      NSm[na,x,t] <- (N[na-1,x,t-1] + NM[na-1,x,t-1]) * SN[na-1,x,t-1] + 
        (N[na,x,t-1] + NM[na,x,t-1]) * SN[na,x,t-1]
      NSv[na,x,t] <- (N[na-1,x,t-1] + NM[na-1,x,t-1]) * SN[na-1,x,t-1] * 
        (1 - SN[na-1,x,t-1]) + 
        (N[na,x,t-1] + NM[na,x,t-1]) * SN[na,x,t-1] * (1 - SN[na,x,t-1])
      for(a in 2:na){
        NS[a,x,t] ~ dnorm(NSm[a,x,t], 1/NSv[a,x,t])
        
        #----- Harvest Survival -----#
        # Survive harvest pressure from August of t to November of t
        Nm[a,x,t] <- NS[a,x,t] * SH[a,x,t]
        Nv[a,x,t] <- Nm[a,x,t] * (1 - SH[a,x,t])
        N[a,x,t] ~ dnorm(Nm[a,x,t], 1/Nv[a,x,t])
        
        NH[a,x,t] <- NS[a,x,t] - N[a,x,t]
      } # a
    } # x
  } # t
  #======= Initial Abundance =======#
  for(a in 1:na){
    for(x in 1:2){
      N[a,x,1] ~ dnorm(N1m[a,x], 0.0001)
    } # x
  } # a
  #======= Parameter Tracking/Organization =======#
  for(t in 1:nt){
    Naf[t] <- sum(N[2:na,1,t])
    Nam[t] <- sum(N[2:na,2,t])
    Nca[t] <- sum(N[1,1:2,t])
  }
  
  #############################################
  ########### Pregnancy Module (PP) ###########
  #############################################
  #======= Process Model Tie-in =======#
  # Must account for variable lactation rates in each age class
  for(t in 1:(nt-1)){
    PP[1,t] <- 0
    PP[2:3,t] <- Py[t] # lactation not observed in young 
    PP[4:(na-1),t] <- (Ppl[t] * Npl[t] + Ppd[t] * Npd[t]) / (Npl[t] + Npd[t])
    PP[na,t] <- (Pol[t] * Nol[t] + Pod[t] * Nod[t]) / (Nol[t] + Nod[t])
  }
  #======= Priors and GLMM =======#
  BPC ~ dlogis(0, 1)
  for(c in 1:5){
    BP0m[c] ~ dlogis(0, 1)
    BP0v[c] ~ dunif(0, 5)
    BPW[c] ~ dlogis(0, 1)
    BPE[c] ~ dlogis(0, 1)
    for(t in 1:(nt-1)){
      BP0[c,t] ~ dnorm(BP0m[c], 1/BP0v[c])
      BP[c,t] <- BP0[c,t] + 
        BPW[c] * clim[t] +
        BPE[c] * nelk[t] +
        BPC * puma[t]
    }
  }
  #======= Observation Model =======#
  for(i in 1:nPo){
    # p_dt: pregnancy data -- each row is an observation with columns:
      # 1 - class (age/lactation); 1-5
      # 2 - year (integer in 1:nt)
      # 3 - number pregnant
      # 4 - number observed
    logit(preg_p[i]) <- BP[p_dt[i,1],p_dt[i,2]]
    p_dt[i,3] ~ dbinom(preg_p[i], p_dt[i,4])
  }
  
  #======= Parameter Tracking/Organization =======#
  for(t in 1:(nt-1)){
    logit(Py[t])  <- BP[1,t]
    logit(Ppl[t]) <- BP[2,t]
    logit(Ppd[t]) <- BP[3,t]
    logit(Pol[t]) <- BP[4,t]
    logit(Pod[t]) <- BP[5,t]
  }
  beta_preg_ynga_clim <- BPW[1]
  beta_preg_prml_clim <- BPW[2]
  beta_preg_prmd_clim <- BPW[3]
  beta_preg_oldl_clim <- BPW[4]
  beta_preg_oldd_clim <- BPW[5]
  beta_preg_ynga_nelk <- BPE[1]
  beta_preg_prml_nelk <- BPE[2]
  beta_preg_prmd_nelk <- BPE[3]
  beta_preg_oldl_nelk <- BPE[4]
  beta_preg_oldd_nelk <- BPE[5]
  beta_preg_puma <- BPC
  
  ####################################################
  ########### Neonate Survival Module (SR) ###########
  ####################################################
  #======= Priors and GLMM =======#
  BSR0m ~ dlogis(0, 1)
  BSR0v ~ dunif(0, 5)
  BSRW ~ dlogis(0, 1)
  BSRE ~ dlogis(0, 1)
  BSRC ~ dlogis(0, 1)
  for(t in 1:(nt-1)){
    BSR0[t] ~ dnorm(BSR0m, 1/BSR0v)
    logit(SR[t]) <- BSR0[t] + 
      BSRW * clim[t] +
      BSRE * nelk[t] +
      BSRC * puma[t]
  }
  #======= Parameter Tracking/Organization =======#
  beta_neos_clim <- BSRW
  beta_neos_puma <- BSRC
  beta_neos_nelk <- BSRE
  
  #########################################################
  ########### Capture-Recapture Module (SH, SN) ###########
  #########################################################
  #======= Process Model Connections =======#
  for(t in 1:(nt-1)){
    SH[1,1,t] <- Shjf[t]
    SH[2:na,1,t] <- Shaf[t]
    SH[1,2,t] <- Shjm[t]
    SH[2:na,2,t] <- Sham[t]
    
    logit(SN[2:na,1,t]) <- Shaf[t] <- logit(Snaf[t]) - 
      error_correction[1,t] * needs_correction[t]
    logit(SN[2:na,2,t]) <- logit(Snam[t]) - 
      error_correction[2,t] * needs_correction[t]
    logit(SN[1,1:2,t]) <- logit(Snca[t]) - 
      error_correction[3,t] * needs_correction[t]
  }
  #======= Observation Model (m-arrays) =======#
  for(t in 1:(nt-1)){
    ma.af[t,(2*t-1):(2*nt-1)] ~ dmulti(p.af[t,(2*t-1):(2*nt-1)], r.af[t])
    ma.am[t,(2*t-1):(2*nt-1)] ~ dmulti(p.am[t,(2*t-1):(2*nt-1)], r.am[t])
    ma.jf[t,(2*t-1):(2*nt-1)] ~ dmulti(p.jf[t,(2*t-1):(2*nt-1)], r.jf[t])
    ma.jm[t,(2*t-1):(2*nt-1)] ~ dmulti(p.jm[t,(2*t-1):(2*nt-1)], r.jm[t])
  } # t
  
  #======= Priors and GLMM =======#
  #----- Hyperpriors -----#
  BSF0m ~ dlogis(0, 1)
  BSM0m ~ dlogis(0, 1)
  BSC0m ~ dlogis(0, 1)
  BSF0v ~ dunif(0, 5)
  BSM0v ~ dunif(0, 5)
  BSC0v ~ dunif(0, 5)
  
  BSFWm ~ dlogis(0, 1)
  BSFWp ~ dlogis(0, 1)
  BSFE  ~ dlogis(0, 1)
  BSFC  ~ dlogis(0, 1)
  BSCWm ~ dlogis(0, 1)
  BSCWp ~ dlogis(0, 1)
  BSCE  ~ dlogis(0, 1)
  BSCC  ~ dlogis(0, 1)
  
  for(t in 1:(nt-1)){
    #----- Detection Probability -----#
    Pf[t] ~ dunif(0, 1)
    Pm[t] ~ dunif(0, 1)
    Qf[t] <- 1 - Pf[t]
    Qm[t] <- 1 - Pm[t]
    
    #----- Natural survival -----#
    for(i in 1:3){
      error_correction[i,t] ~ dnorm(0, 0.001)
    }
    B0Snaf[t] ~ dnorm(BSF0m, 1/BSF0v)
    B0Snam[t] ~ dnorm(BSF0m, 1/BSF0v)
    B0Snca[t] ~ dnorm(BSF0m, 1/BSF0v)
    logit(Snaf[t]) <- B0Snaf[t] +
      BSFWm * clim[t] +
      BSFWp * clim[t+1] +
      BSFE * nelk[t] +
      BSFC * puma[t] +
      error_correction[1,t] * needs_correction[t]
    logit(Snam[t]) <- B0Snam[t] + error_correction[2,t] * needs_correction[t]
    logit(Snca[t]) <- B0Snca[t] +
      BSCWm * clim[t] +
      BSCWp * clim[t+1] +
      BSCE * nelk[t] +
      BSCC * puma[t] +
      error_correction[3,t] * needs_correction[t]
    
    #----- Harvest Probability -----#
    Haf[t] ~ dunif(0, 1)
    Ham[t] ~ dunif(0, 1)
    Hjf[t] ~ dunif(0, 1)
    Hjm[t] ~ dunif(0, 1)
    Shaf[t] <- (1 - Haf[t] * afh[t])
    Sham[t] <- (1 - Ham[t] * amh[t])
    Shjf[t] <- (1 - Hjf[t] * jfh[t])
    Shjm[t] <- (1 - Hjm[t] * jmh[t])
    
    #----- Combined Survival -----#
    Saf[t] <- Snaf[t] * Shaf[t]
    Sam[t] <- Snam[t] * Sham[t]
    Sjf[t] <- Snca[t] * Shjf[t]
    Sjm[t] <- Snca[t] * Shjm[t]
  } # t
  
  #======= Cell Observation Probabilities =======#
  for(t in 1:(nt-1)){
    #----- Seen Alive Next Session -----#
    # Adult
    p.af[t,t*2-1] <- Saf[t] * Pf[t]
    p.am[t,t*2-1] <- Sam[t] * Pm[t]
    # Calf 
    p.jf[t,t*2-1] <- Sjf[t] * Pf[t]
    p.jm[t,t*2-1] <- Sjm[t] * Pm[t]
    
    #----- Harvested Next Session -----#
    # Adult
    p.af[t,t*2] <- Snaf[t] * (1 - Shaf[t])
    p.am[t,t*2] <- Snam[t] * (1 - Sham[t])
    # Calf 
    p.jf[t,t*2] <- Snca[t] * (1 - Shjf[t])
    p.jm[t,t*2] <- Snca[t] * (1 - Shjm[t])
  } # t
  
  for(t in 1:(nt-2)){
    for(u in (t+1):(nt-1)){
      #----- Seen Alive In `x` Sessions -----#
      # Adult
      p.af[t,u*2-1] <- prod(Saf[t:u]) * prod(Qf[t:(u-1)]) * Pf[u]
      p.am[t,u*2-1] <- prod(Sam[t:u]) * prod(Qm[t:(u-1)]) * Pm[u]
      # Calf 
      p.jf[t,u*2-1] <- Sjf[t] * prod(Saf[(t+1):u]) * prod(Qf[t:(u-1)]) * Pf[u]
      p.jm[t,u*2-1] <- Sjm[t] * prod(Sam[(t+1):u]) * prod(Qm[t:(u-1)]) * Pm[u]
      
      #----- Harvested In `x` Sessions -----#
      # Adult 
      p.af[t,u*2] <- prod(Snaf[t:u]) * 
        prod(Qf[t:(u-1)]) * 
        prod(Shaf[t:(u-1)]) * 
        (1 - Shaf[u])
      p.am[t,u*2] <- prod(Snam[t:u]) * 
        prod(Qm[t:(u-1)]) * 
        prod(Sham[t:(u-1)]) * 
        (1 - Sham[u])
      # Calf 
      p.jf[t,u*2] <- Sjf[t] * 
        prod(Snaf[(t+1):u]) * 
        prod(Qf[t:(u-1)]) * 
        prod(Shaf[(t+1):(u-1)]) * 
        (1 - Shaf[u])
      p.jm[t,u*2] <- Sjm[t] * 
        prod(Snam[(t+1):u]) * 
        prod(Qm[t:(u-1)]) * 
        prod(Sham[(t+1):(u-1)]) * 
        (1 - Sham[u])
    } # u
  } # t
  
  for(t in 2:(nt-1)){
    #----- Structural Zeroes -----#
    # Adult 
    p.af[t,1:((t-1)*2)] <- 0
    p.am[t,1:((t-1)*2)] <- 0
    # Calf 
    p.jf[t,1:((t-1)*2)] <- 0
    p.jm[t,1:((t-1)*2)] <- 0
  } # t
  
  for(t in 1:(nt-1)){
    #----- Not Observed Again -----#
    # Adult 
    p.af[t,nt*2-1] <- 1 - sum(p.af[t,1:(2*nt-2)])
    p.am[t,nt*2-1] <- 1 - sum(p.am[t,1:(2*nt-2)])
    # Calf 
    p.jf[t,nt*2-1] <- 1 - sum(p.jf[t,1:(2*nt-2)])
    p.jm[t,nt*2-1] <- 1 - sum(p.jm[t,1:(2*nt-2)])
  } # t
})





