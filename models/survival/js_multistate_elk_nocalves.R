code <- nimbleCode({
  # 24-August-2022
  
  #Priors=======================================================================
  for(t in 1:(nocc-1)){
    # time varying effects on
    # survival (s), entry (e), and detection (p) probabilities
    # including females (f), males (m), and calves (c)
    sf[t] ~ dunif(0, 1)
    sm[t] ~ dunif(0, 1)
    sc[t] ~ dunif(0, 1)
    pf[t] ~ dunif(0, 1)
    pm[t] ~ dunif(0, 1)
    pc[t] <- 1 # should be no calves during the initial all 0 capture session
    ep[t] ~ dunif(0, 1)
  }

  # Fixed effects of non-main study herd (h) on p and s
  sh ~ dunif(0, 1)
  ph ~ dunif(0, 1)

  #Probabilities================================================================
  for(t in 1:(nocc-1)){
    for(i in 1:nind){
      s[i,t] <- (1 - h[i,t]) * sm[t] *      m[i]  * (1 - c[i,t]) + 
                (1 - h[i,t]) * sf[t] * (1 - m[i]) * (1 - c[i,t]) + 
                (1 - h[i,t]) * sc[t] * c[i,t] +
                     h[i,t]  * sh
      
      p[i,t] <- (1 - h[i,t]) * pm[t] *      m[i]  * (1 - c[i,t]) + 
                (1 - h[i,t]) * pf[t] * (1 - m[i]) * (1 - c[i,t]) + 
                (1 - h[i,t]) * pc[t] * c[i,t] +
                     h[i,t]  * ph
      
      # State transition probabilities. States are:
      # 1 -- does not exist yet
      # 2 -- alive
      # 3 -- dead
      ps[1,i,t,1] <- 1 - ep[t]
      ps[1,i,t,2] <- ep[t]
      ps[1,i,t,3] <- 0
      ps[2,i,t,1] <- 0
      ps[2,i,t,2] <- s[i,t]
      ps[2,i,t,3] <- 1 - s[i,t]
      ps[3,i,t,1] <- 0
      ps[3,i,t,2] <- 0
      ps[3,i,t,3] <- 1
      
      # State detection probabilities.
      po[1,i,t,1] <- 0
      po[1,i,t,2] <- 1
      po[2,i,t,1] <- p[i,t]
      po[2,i,t,2] <- 1 - p[i,t]
      po[3,i,t,1] <- 0
      po[3,i,t,2] <- 1
    } #i
  } #t
  
  #Model========================================================================
  for(i in 1:nind){
    m[i] ~ dbern(0.5)
    z[i,1] <- 1
    for(t in 2:nocc){
      # State Process -------- z ~ c(1, 2, 3)
      z[i,t] ~ dcat(ps[ z[i,t-1], i, t-1, 1:3])
      
      # Observation Process -- seen >>> y = 2
      y[i,t] ~ dcat(po[ z[i,t],   i, t-1, 1:2])
    } #t
  } #i
  
  #Derived======================================================================
  # Abundance
  for(t in 1:nocc){
    for(i in 1:nind){
      a.al[i,t]  <- equals(z[i,t], 2)
      m.al[i,t]  <- ad.al[i,t] *      m[i]  * (1 - c[i,t])
      f.al[i,t]  <- ad.al[i,t] * (1 - m[i]) * (1 - c[i,t])
    } #i
    Nm[t] <- sum(m.al[1:nind,t])
    Nf[t] <- sum(f.al[1:nind,t])
    Nt[t] <- sum(a.al[1:nind,t])
  } #t
  
  # Superpopulation
  for(i in 1:nind){
    nt.al[i] <- sum(a.al[i,1:nocc])
    w[i]     <- 1 - equals(nt.al[i], 0)
  }
  Ns <- sum(w[1:nind])
})