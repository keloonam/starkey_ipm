code <- nimbleCode({
  # Code for annual, multi-state capture recapture model of Starkey elk
  # Kenneth Loonam
  # November 2025
  
  ######## Likelihood ########
  for(t in 1:(nt-1)){
    ma.af[t,(2*t-1):(2*nt-1)] ~ dmulti(p.af[t,(2*t-1):(2*nt-1)], r.af[t])
    ma.am[t,(2*t-1):(2*nt-1)] ~ dmulti(p.am[t,(2*t-1):(2*nt-1)], r.am[t])
    ma.jf[t,(2*t-1):(2*nt-1)] ~ dmulti(p.jf[t,(2*t-1):(2*nt-1)], r.jf[t])
    ma.jm[t,(2*t-1):(2*nt-1)] ~ dmulti(p.jm[t,(2*t-1):(2*nt-1)], r.jm[t])
  } # t
  
  ######## Priors ########
  for(t in 1:(nt-1)){
    #===== Detection Probability =====#
    Pf[t] ~ dunif(0, 1)
    Pm[t] ~ dunif(0, 1)
    Qf[t] <- 1 - Pf[t]
    Qm[t] <- 1 - Pm[t]
    
    #===== Natural survival =====#
    Snaf[t] ~ dunif(0, 1)
    Snam[t] ~ dunif(0, 1)
    Snca[t] ~ dunif(0, 1)
    
    #===== Harvest Survival =====#
    Shaf[t] ~ dunif(0, 1)
    Sham[t] ~ dunif(0, 1)
    Shjf[t] ~ dunif(0, 1)
    Shjm[t] ~ dunif(0, 1)
    
    #===== Combined Survival =====#
    Saf[t] <- Snaf[t] * Shaf[t]
    Sam[t] <- Snam[t] * Sham[t]
    Sjf[t] <- Snca[t] * Shjf[t]
    Sjm[t] <- Snca[t] * Shjm[t]
  } # t
  
  ######## Cell Observation Probabilities ########
  for(t in 1:(nt-1)){
    #===== Seen Alive Next Session =====#
    #--- Adult ---#
    p.af[t,t*2-1] <- Saf[t] * Pf[t]
    p.am[t,t*2-1] <- Sam[t] * Pm[t]
    #--- Calf ---#
    p.jf[t,t*2-1] <- Sjf[t] * Pf[t]
    p.jm[t,t*2-1] <- Sjm[t] * Pm[t]
    
    #===== Harvested Next Session =====#
    #--- Adult ---#
    p.af[t,t*2] <- Snaf[t] * (1 - Shaf[t])
    p.am[t,t*2] <- Snam[t] * (1 - Sham[t])
    #--- Calf ---#
    p.jf[t,t*2] <- Snca[t] * (1 - Shjf[t])
    p.jm[t,t*2] <- Snca[t] * (1 - Shjm[t])
  } # t
  
  for(t in 1:(nt-2)){
    for(u in (t+1):(nt-1)){
      #===== Seen Alive In `x` Sessions =====#
      #--- Adult ---#
      p.af[t,u*2-1] <- prod(Saf[t:u]) * prod(Qf[t:(u-1)]) * Pf[u]
      p.am[t,u*2-1] <- prod(Sam[t:u]) * prod(Qm[t:(u-1)]) * Pm[u]
      #--- Calf ---#
      p.jf[t,u*2-1] <- Sjf[t] * prod(Saf[(t+1):u]) * prod(Qf[t:(u-1)]) * Pf[u]
      p.jm[t,u*2-1] <- Sjm[t] * prod(Sam[(t+1):u]) * prod(Qm[t:(u-1)]) * Pm[u]
      
      #===== Harvested In `x` Sessions =====#
      #--- Adult ---#
      p.af[t,u*2] <- prod(Snaf[t:u]) * 
        prod(Qf[t:(u-1)]) * 
        prod(Shaf[t:(u-1)]) * 
        (1 - Shaf[u])
      p.am[t,u*2] <- prod(Snam[t:u]) * 
        prod(Qm[t:(u-1)]) * 
        prod(Sham[t:(u-1)]) * 
        (1 - Sham[u])
      #--- Calf ---#
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
    #===== Structural Zeroes =====#
    #--- Adult ---#
    p.af[t,1:((t-1)*2)] <- 0
    p.am[t,1:((t-1)*2)] <- 0
    #--- Calf ---#
    p.jf[t,1:((t-1)*2)] <- 0
    p.jm[t,1:((t-1)*2)] <- 0
  } # t
  
  for(t in 1:(nt-1)){
    #===== Not Observed Again =====#
    #--- Adult ---#
    p.af[t,nt*2-1] <- 1 - sum(p.af[t,1:(2*nt-2)])
    p.am[t,nt*2-1] <- 1 - sum(p.am[t,1:(2*nt-2)])
    #--- Calf ---#
    p.jf[t,nt*2-1] <- 1 - sum(p.jf[t,1:(2*nt-2)])
    p.jm[t,nt*2-1] <- 1 - sum(p.jm[t,1:(2*nt-2)])
  } # t
  
  ######## Goodness-of-Fit ########
  for(t in 1:(nt-1)){
    #===== Replicate Data =====#
    mn.af[t,1:(2*nt-1)] ~ dmulti(p.af[t,1:(2*nt-1)], r.af[t])
    mn.am[t,1:(2*nt-1)] ~ dmulti(p.am[t,1:(2*nt-1)], r.am[t])
    mn.jf[t,1:(2*nt-1)] ~ dmulti(p.jf[t,1:(2*nt-1)], r.jf[t])
    mn.jm[t,1:(2*nt-1)] ~ dmulti(p.jm[t,1:(2*nt-1)], r.jm[t])
    
    for(u in 1:(2*nt-1)){
      #===== Calculate Expected Values =====#
      ex.af[t,u] <- p.af[t,u] * r.af[t]
      ex.am[t,u] <- p.am[t,u] * r.am[t]
      ex.jf[t,u] <- p.jf[t,u] * r.jf[t]
      ex.jm[t,u] <- p.jm[t,u] * r.jm[t]
      
      #===== Freeman-Tukey Statistic =====#
      #--- Original Data --- #
      E.af.org[t,u] <- (sqrt(ma.af[t,u]) - sqrt(ex.af[t,u]))^2
      E.am.org[t,u] <- (sqrt(ma.am[t,u]) - sqrt(ex.am[t,u]))^2
      E.jf.org[t,u] <- (sqrt(ma.jf[t,u]) - sqrt(ex.jf[t,u]))^2
      E.jm.org[t,u] <- (sqrt(ma.jm[t,u]) - sqrt(ex.jm[t,u]))^2
      #--- New Data --- #
      E.af.new[t,u] <- (sqrt(mn.af[t,u]) - sqrt(ex.af[t,u]))^2
      E.am.new[t,u] <- (sqrt(mn.am[t,u]) - sqrt(ex.am[t,u]))^2
      E.jf.new[t,u] <- (sqrt(mn.jf[t,u]) - sqrt(ex.jf[t,u]))^2
      E.jm.new[t,u] <- (sqrt(mn.jm[t,u]) - sqrt(ex.jm[t,u]))^2
    }
  }
  #===== Calculate Final Statistic =====#
  #--- New Data --- #
  fit.af.new <- sum(E.af.new[1:(nt-1), 1:(2*nt-1)])
  fit.am.new <- sum(E.am.new[1:(nt-1), 1:(2*nt-1)])
  fit.jf.new <- sum(E.jf.new[1:(nt-1), 1:(2*nt-1)])
  fit.jm.new <- sum(E.jm.new[1:(nt-1), 1:(2*nt-1)])
  fit.new <- fit.af.new + fit.am.new + fit.jf.new + fit.jm.new
  #--- Original Data --- #
  fit.af.org <- sum(E.af.org[1:(nt-1), 1:(2*nt-1)])
  fit.am.org <- sum(E.am.org[1:(nt-1), 1:(2*nt-1)])
  fit.jf.org <- sum(E.jf.org[1:(nt-1), 1:(2*nt-1)])
  fit.jm.org <- sum(E.jm.org[1:(nt-1), 1:(2*nt-1)])
  fit.org <- fit.af.org + fit.am.org + fit.jf.org + fit.jm.org
})
