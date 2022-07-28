js_superpop_model <- nimbleCode({
  # Nimble version of jolly-seber super population formulation with age data
  # Constants and data are CAPITALIZED
  # Estimated parameters are lowercase
  # Basically, if it needs to be passed to the model it is in all caps
  # The exceptions to this are y and z. Just because it would bother me. :D
  
  ##### Priors #####
  e0 ~ dbeta(1, 1) # probability animal is in super-population (i.e. exists)
  
  for(t in 1:K){
    p0[t]   ~ dbeta(1, 1) # capture probability
    p.b0[t] <- logit(p0[t])
    p.ca[t] ~ dbeta(1, 1) # calf adjustment
    p.bc[t] <- logit(p.ca[t])
    p.ma[t] ~ dbeta(1, 1) # male adjustment
    p.bm[t] <- logit(p.ma[t])
    p.he[t] ~ dbeta(1, 1) # herd adjustment
    p.bh[t] <- logit(p.he[t])
  }
  
  for(t in 1:(K)){
    s0[t]   ~ dbeta(1, 1) # survival probability
    s.b0[t] <- logit(s0[t])
    s.ca[t] ~ dbeta(1, 1) # calf adjustment
    s.bc[t] <- logit(s.ca[t])
    s.ma[t] ~ dbeta(1, 1) # male adjustment
    s.bm[t] <- logit(s.ma[t])
    s.he[t] ~ dbeta(1, 1) # herd adjustment
    s.bh[t] <- logit(s.he[t])
  }
  
  # age distribution - given alive, what is the probability of being age 1:max?
  a0[1:MAX_A] ~ ddirch(A[1:MAX_A])
  
  # age dist including unentered individuals, c(p(!entered), p(age(1:max))) == 1
  # aka, age probabilities at start inclusive (i) of not yet entered age (0)
  a0_i[1:(MAX_A + 1)] <- c((1 - eta[1]), eta[1] * a0[1:MAX_A])
  
  # recruitment rate from 'not entered population' at t 
  # r_n is a vector of length n_occasions that sums to 1
    # it represents the naive probability of individuals entering. It is not 
    # adjusted by time-step/n ind remaining
  # r is the probability of the remaining true individuals to enter the pop at
    # each time step. The final value of r is 1.
  r_n[1:K] ~ ddirch(B[1:K])
  r[1] <- r_n[1]
  for(k in 2:K){
    r[k] <- r_n[k]/(1 - sum(r_n[1:(k-1)]))
  }
  
  ##### Model #####
  ### Session 1 ###
  for (i in 1:M){
    # is individual i real?
    w[i] ~ dbern(e0)
    
    # initial ages
    A_P1[i] ~ dcat(a0_i[1:(MAX_A + 1)]) # A_P1 = age plust one
    AGE[i,1] <- (A_P1[i] - 1) 
    c[i,1] <- AGE[i,1] == 1 # is i a calf?
    # dcat gives integers from 1 to the length of the vector passed
    # age 0 can't be given, so add one to all initial ages.
    # AGE is also data, but has NAs for augmented ind and ind w/o age
    
    # state process
    u[i,1] <- step(AGE[i,1]-.1) # sets u to zero if age is zero, 1 otherwise
    z[i,1] <- u[i,1] * w[i] # z is the "real" state
    
    # Observation probability
    logit(p[i,1]) <- p.b0[1] + p.bc[1]*c[i,1] + p.bm[1]*M[i] + p.bh[1]*H[i,t]
    
    # Observation process
    y[i,1] ~ dbern(z[i,1] * p[i,1])
    
    # derived stuff
    avail[i,1] <- 1 - u[i,1] # still available -- i.e. not yet recruited
    
    ### Sessions 2:K ###     
    for (t in 2:K[L[i]]){ # 2 to last session i was available
      # Survival probabilities
      logit(s[i,t]) <- s.b0[t] + s.bc[t]*c[i,t] + s.bm[t]*M[i] + s.bh[t]*H[i,t]
      
      # Detection probabilities
      logit(p[i,t]) <- p.b0[t] + p.bc[t]*c[i,t] + p.bm[t]*M[i] + p.bh[t]*H[i,t]
      
      # State process
      u[i,t] ~ dbern(u[i,t-1] * s[i,t] + avail[i,t-1] * r[t])   
      z[i,t] <- u[i,t] * w[i]
      
      # Age process
      AGE[i,t] <- AGE[i,t-1] + max(u[i,1:t]) 
      # ages by one year after recruitment (NIMBLE allows this syntax?)
      c[i,t] <- AGE[i,t] == 1 # is i a calf?
      
      # Observation process
      y[i,t] ~ dbern(z[i,t] * p[i,t])
      
      # derived stuff
      avail[i,t] <- 1- max(u[i,1:t]) # still available -- i.e. not yet recruited
    } #t
  } #i
  
  ##### Derived/Tracked Values #####
  for(t in 1:K){
    logit(survival_af[t]) <- s.b0[t]
    logit(survival_am[t]) <- s.b0[t] + s.bm[t]
    logit(survival_ca[t]) <- s.b0[t] + s.bc[t]
    logit(detection_f[t]) <- p.b0[t]
    logit(detection_m[t]) <- b.b0[t] + p.bm[t]
    logit(detection_c[t]) <- b.b0[t] + p.bc[t]
  }
  # Annual abundance
  for(t in 1:K){
    N[t] <- sum(z[1:M,t])               
  } #t
  
  # Annual growth rate
  for(t in 1:(K-1)){
    lambda[t] <- N[t+1]/N[t]               
  } #t
  
  # Super-population size
  N_super <- sum(w[1:M])       
  
})# end model