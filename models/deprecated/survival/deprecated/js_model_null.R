js_nimble_code <- nimbleCode({
  # Nimble version of jolly-seber super population formulation with age data
  #Priors=======================================================================
  e    ~ dbeta(1, 1)
  p.b0 ~ dlogis(0, 1)
  s.b0 ~ dlogis(0, 1)
  r_n[1:nocc] ~ ddirch(R[1:nocc])
  r[1] <- r_n[1]
  
  for(t in 2:nocc){
    r[t] <- r_n[t] / (1 - sum(r_n[1:(t-1)]))
  }
  for(i in 1:naug){
    for(t in 1:(nocc-1)){
      logit(s[i,t]) <- s.b0
    }
    for(t in 1:nocc){
      logit(p[i,t]) <- p.b0
    }
  }
  #Likelihood===================================================================
  for(i in 1:naug){
    ##### First Occasion
    ### State Process
    w[i] ~ dbern(e)
    z[i,1] ~ dbern(r[1])

    ### Observation Process
    y[i,1] ~ dbern(z[i,1] * p[i,1] * w[i])

    for(t in 2:nocc){
      ##### Subsequent Occasions
      ### State Process
      q[i,(t-1)] <- 1 - z[i,t-1]
      z[i,t] ~ dbern(z[i,t-1] * s[i,t-1] + prod(q[i,1:(t-1)]) * r[t])

      ### Observation Process
      y[i,t] ~ dbern(z[i,t] * p[i,t] * w[i])
    }
  }
  #Derived_Values===============================================================
  for(i in 1:naug){
    for(t in 1:nocc){
      u[i,t] <- z[i,t] * w[i]
    }
  }
  
  for(t in 1:nocc){
    N[t] <- sum(u[1:naug,t])
  }
  logit(survival_af) <- s.b0
  logit(detection_f) <- p.b0
  N_super <- sum(w[1:naug])
})