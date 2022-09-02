code <- nimbleCode({
  # cjs code for Starkey elk with js abundance estimate added
  # 18-August-2022
  
  #Priors=======================================================================
  for(t in 1:(nocc)){
    # time varying effects on
    # detection probability (p), and survival probability (s)
    # including intercepts (0), males (m), and calves (c)
    pf[t] ~ dunif(0, 1)
    pm[t] ~ dunif(0, 1)
    sf[t] ~ dunif(0, 1)
    sm[t] ~ dunif(0, 1)
    sc[t] ~ dunif(0, 1)
  }
  
  # Fixed effects of non-main study herd (h) on p and s
  sh ~ dunif(0, 1)
  ph ~ dunif(0, 1)
  
  for(i in 1:nind){
    for(t in f[i]:l[i]){
      # probabilities to actually use
      s[i,t] <- 0 +
        sf[t] * (1 - m[i]) * (1 - c[i,t]) * (1 - h[i,t]) + 
        sc[t] *                   c[i,t]  * (1 - h[i,t]) + 
        sm[t] *      m[i]  * (1 - c[i,t]) * (1 - h[i,t]) + 
        sh    *                                  h[i,t]
      
      p[i,t] <- 0 +
        pf[t] * (1 - m[i]) * (1 - c[i,t]) * (1 - h[i,t]) + 
        pm[t] *      m[i]  * (1 - c[i,t]) * (1 - h[i,t]) + 
        ph    *                                  h[i,t]
    } #i
  } #t
  
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