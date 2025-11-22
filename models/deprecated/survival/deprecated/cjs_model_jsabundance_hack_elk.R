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
  
  #Abundance====================================================================
  for(i in 1:nind){
    for(t in f[i]:l[i]){
      # Record identities: adult male = m.a; adult female = f.a
      m.a[i,t] <- z[i,t] *      m[i]  * (1 - c[i,t])
      f.a[i,t] <- z[i,t] * (1 - m[i]) * (1 - c[i,t])
    } #t
  } #i
  
  for(t in 2:nocc){
    # Number of adult males/females observed
    Nm.a[t] <- sum(m.a[1:nind,t])
    Nf.a[t] <- sum(f.a[1:nind,t])
    
    # Expected value of N previously unobserved males/females
    expNm.un[t] <- newNm[t] / p.m[t]
    expNf.un[t] <- newNf[t] / p.f[t]
    
    # SD of expected value (normal approx. of binomial)
    sdNm.un[t] <- sqrt(expNm.un[t] * p.m[t] * (1 - p.m[t]))
    sdNf.un[t] <- sqrt(expNf.un[t] * p.f[t] * (1 - p.f[t]))
    
    # Estimate of previously unobserved males/females
    Nm.un[t] ~ T(dnorm(expNm.un[t], sd = sdNm.un[t]), 0, 1)
    Nf.un[t] ~ T(dnorm(expNf.un[t], sd = sdNf.un[t]), 0, 1)
    
    # N = previously unobserved + know alive - 2*newly observed
    # Newly observed included in previously unobserved and n.alive
    Nm[t] <- Nm.un[t] + Nm.a[t] - 2*newNm[t]
    Nf[t] <- Nf.un[t] + Nm.a[t] - 2*newNm[t]
  } #t
})