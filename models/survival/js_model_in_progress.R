js_model_code <- nimbleCode({
  # Nimble version of jolly-seber super population formulation with age data
  #Priors=======================================================================
  e    ~ dbeta(1, 1)
  p.b0 ~ dlogis(0, 1)
  s.b0 ~ dlogis(0, 1)
  r_n[1:nocc] ~ ddirch(rep(1, nocc))
  r[1] <- r_n[1]
  
  for(t in 1:nocc){
    logit(p0[t]) <- p.b0
    logit(s0[t]) <- s.b0
  }
  for(t in 2:nocc){
    r[t] <- r_n[t] / (1 - sum(r_n[1:(t-1)]))
  }
  
  #Likelihood===================================================================
  for (i in 1:naug){
    ##### First Occasion
    ### Probabilities
    logit(p[i,1]) <- p.b0
    
    ### State Process
    w[i] ~ dbern(e)
    z[i,1] <- dbern(r[1] * w[i])
    
    ### Observation Process
    y[i,1] ~ dbern(z[i,1] * p[i,1])
    
    for (t in 2:(L[i])){
      ##### Subsequent Occasions
      ### Probabilities
      logit(s[i,t]) <- s.b0
      logit(p[i,t]) <- p.b0
      
      ### State Process
      q[i,t-1] <- 1 - z[i,1:t-1]
      z[i,t] ~ dbern(z[i,t-1] * s[i,t] + prod(q[i,1:(t-1)] * r[t] * w[i])
                     
      ### Observation Process
      y[i,t] ~ dbern(z[i,t] * p[i,t])
    } #t
  } #i
  
  #Derived_Values===============================================================
  for(t in 1:nocc){
    N[t] <- sum(z[1:naug,t])
  }
  logit(survival_af) <- s.b0
  logit(detection_f) <- p.b0
  N_super <- sum(w[1:naug])
  
})# end model