js_nimble_code <- nimbleCode({
  # cjs code for Starkey elk with js abundance estimate added
  
  #Priors=======================================================================
  s.b0 ~ dlogis(0, 1)
  p.b0 ~ dlogis(0, 1)
  
  for(i in 1:nind){
    for(t in f[i]:(nocc - 1)){
      logit(s[i,t]) <- s.b0
      logit(p[i,t]) <- p.b0
    } #t
    for(t in 1:nocc){
      z.constrain[i,t] ~ dconstraint(z[i,t] >= 0)
    }
  } #i
  for(t in 1:nocc){
    N.new.bin[t] ~ T(dnorm(0, sd = 1000), 0, 1)
  }
  
  #Likelihood===================================================================
  for(i in 1:nind){
    z[i,f[i]] <- 1
    for(t in (f[i]+1):nocc){
      # State Process
      z[i,t] ~ dbern(mu1[i,t])
      mu1[i,t] <- s[i,t-1] * z[i,t-1]
      # new.ind[i,t] <- z[i,t] * (1 - z[i,t-1])
      # old.ind[i,t] <- z[i,t] * z[i,t-1]
      
      # Observation Process
      y[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- p[i,t-1] * z[i,t]
    } #t
  } #i
  
  # Abundance Process
  for(t in 2:nocc){
    # N.seen[t] <- sum(y[1:nind,t])
    # N[t] <- N.seen[t] / mean.p
    # N.old[t] <- sum(old.ind[1:nind,t])
    N.alive[t]  <- sum(z[1:nind, t])
    
    expN.new[t] <- newN[t] / mean.p
    sdN.new[t]  <- sqrt(expN.new[t] * mean.p * (1 - mean.p))
    N.new[t]    ~ dnorm(expN.new[t], sd = sdN.new[t])
    
    newN[t] ~ dbinom(N.new.bin[t], mean.p)
    
    N.bin[t]  <- N.new.bin[t] + N.alive[t] - newN[t]
    N[t]        <- N.new[t] + N.alive[t] - newN[t]
  }
  
  #Values=======================================================================
  logit(mean.s) <- s.b0
  logit(mean.p) <- p.b0
  
})