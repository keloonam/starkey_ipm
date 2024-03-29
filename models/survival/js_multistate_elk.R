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
    pc[t] <- pf[t]
    ec[t] ~ dunif(0, 1)
  }
  ea[1] ~ dunif(0, 1 - ec[1]) # can only enter as an adult from t = 1
  for(t in 2:(nocc - 1)){
    ea[t] <- 0
  }
  # Fixed effects of non-main study herd (h) on p and s
  sh ~ dunif(0, 1)
  ph ~ dunif(0, 1)
  for(i in 1:nind){
    pi[i]
  }
  
  #Probabilities================================================================
  for(t in 1:(nocc-1)){
    for(i in 1:nind){
      s[i,t] <- (1 - h[i,t]) * sm[t] *      m[i]   + 
        (1 - h[i,t]) * sf[t] * (1 - m[i])  + 
        sh * h[i,t]
      p[i,t] <- (1 - h[i,t]) * pm[t] *      m[i]   + 
        (1 - h[i,t]) * pf[t] * (1 - m[i])  + 
        ph * h[i,t]
      # State transition probabilities. States are:
      # 1 -- does not exist yet
      # 2 -- live calf
      # 3 -- live adult
      # 4 -- dead
      ps[1,i,t,1] <- 1 - ec[t] - ea[t]
      ps[1,i,t,2] <- ec[t]
      ps[1,i,t,3] <- ea[t]
      ps[1,i,t,4] <- 0
      ps[2,i,t,1] <- 0
      ps[2,i,t,2] <- 0
      ps[2,i,t,3] <- sc[t]
      ps[2,i,t,4] <- 1 - sc[t]
      ps[3,i,t,1] <- 0
      ps[3,i,t,2] <- 0
      ps[3,i,t,3] <- s[i,t]
      ps[3,i,t,4] <- 1 - s[i,t]
      ps[4,i,t,1] <- 0
      ps[4,i,t,2] <- 0
      ps[4,i,t,3] <- 0
      ps[4,i,t,4] <- 1
      
      # State detection probabilities.
      po[1,i,t,1] <- 0
      po[1,i,t,2] <- 1
      po[2,i,t,1] <- p[i,t]
      po[2,i,t,2] <- 1 - p[i,t]
      po[3,i,t,1] <- p[i,t]
      po[3,i,t,2] <- 1 - p[i,t]
      po[4,i,t,1] <- 0
      po[4,i,t,2] <- 1
    } #i
  } #t
  
  #Model========================================================================
  for(i in 1:nind){
    m[i] ~ dbern(0.5)
    z[i,1] <- 1
    for(t in 2:nocc){
      # State Process -------- z ~ c(1, 2, 3, 4)
      z[i,t] ~ dcat(ps[ z[i,t-1], i, t-1, 1:4])
      
      # Observation Process -- seen >>> y = 2
      y[i,t] ~ dcat(po[ z[i,t],   i, t-1, 1:2])
    } #t
  } #i
  
  #Derived======================================================================
  # Abundance
  for(t in 1:nocc){
    for(i in 1:nind){
      ad.al[i,t] <- equals(z[i,t], 3)
      ca.al[i,t] <- equals(z[i,t], 2)
      m.al[i,t]  <- ad.al[i,t] *      m[i]
      f.al[i,t]  <- ad.al[i,t] * (1 - m[i])
      an.al[i,t] <- ad.al[i,t] + ca.al[i,t]
    } #i
    Nm[t] <- sum(m.al[1:nind,t])
    Nf[t] <- sum(f.al[1:nind,t])
    Nc[t] <- sum(ca.al[1:nind,t])
    Nt[t] <- sum(an.al[1:nind,t])
  } #t
  
  # Superpopulation
  for(i in 1:nind){
    nt.al[i] <- sum(an.al[i,1:nocc])
    w[i]     <- 1 - equals(nt.al[i], 0)
  }
  Ns <- sum(w[1:nind])
})