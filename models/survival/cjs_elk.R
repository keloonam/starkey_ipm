nimble_code <- nimbleCode({
  # cjs code for Starkey elk with js abundance estimate added
  # 18-August-2022
  
  #Priors=======================================================================
  for(t in 1:(nocc-1)){
    # time varying effects on
    # detection probability (p), and survival probability (s)
    # including intercepts (0), males (m), and calves (c)
    p0[t] ~ dlogis(0, 1)
    pm[t] ~ dlogis(0, 1)
    s0[t] ~ dlogis(0, 1)
    sm[t] ~ dlogis(0, 1)
    sc[t] ~ dlogis(0, 1)
  }
  
  # Fixed effects of non-main study herd (h) on p and s
  sh ~ dlogis(0, 1)
  ph ~ dlogis(0, 1)
  
  for(t in f[i]:(nocc - 1)){
    for(i in 1:nind){
      # probabilities to actually use
      logit(s[i,t]) <- 
        s0[t] + 
        sc[t] * c[i,t] + 
        sm[t] * m[i]   * (1 - c[i,t]) + 
        sh    * h[i,t]
      
      logit(p[i,t]) <- 
        p0[t] + 
        pm[t] * m[i]   * (1 - c[i,t]) + 
        ph    * h[i,t]
    } #i
    
    # Tracking values
    logit(s.m[t]) <- s0[t] + sm[t]
    logit(p.m[t]) <- p0[t] + pm[t]
    logit(s.f[t]) <- s0[t]
    logit(p.f[t]) <- p0[t]
    logit(s.c[t]) <- s0[t] + sc[t]
  } #t
  
  for(i in 1:nind){
    for(t in 1:nocc){
      # Constrain z to be a number, not NA
      z_con[i,t] ~ dconstraint(z[i,t] >= 0)
    }
  }
  
  #Survival=====================================================================
  for(i in 1:nind){
    z[i,f[i]] <- 1
    for(t in (f[i]+1):l[i]){
      # State Process -------- alive >>> z = 1
      z[i,t] ~ dbern(mu1[i,t])
      mu1[i,t] <- s[i,t-1] * z[i,t-1]
      
      # Observation Process -- seen >>> y = 1
      y[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- p[i,t-1] * z[i,t]
    } #t
  } #i
})