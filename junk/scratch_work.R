js_superpop_model <- nimbleCode({
  # Nimble version of jolly-seber super population formulation with age data
  # Probabilities are capitalized
  # Data,  and constants are lowercase
  
  # Priors
  E0 ~ dbeta(1, 1)    # probability animal is in superpopulation (i.e. exists)
  P0   ~ dbeta(1, 1)    # mean capture probability
  S0   ~ dbeta(1, 1)    # survival
  s.b0 <- logit(S0)
  p.b0 <- logit(P0)
  
  # age distribution - given alive, what is the probability of being age 1:max?
  piAGE[1:max.age] ~ ddirch(a[1:max.age])
  
  # age dist including unentered individuals, c(p(!entered), p(age(1:max))) == 1
  piAGEuncond[1:(max.age + 1)] <- c((1 - eta[1]), eta[1] * piAGE[1:max.age])
  
  # recruitment rate from 'not entered population' at t 
  beta[1:K] ~ ddirch(b[1:K])
  eta[1] <- beta[1]
  for(k in 2:K){
    eta[k] <- beta[k]/(1-sum(beta[1:(k-1)]))
  }
  
  # Likelihoods 
  for (i in 1:M){
    # is individual i real?
    w[i] ~ dbern(psi)
    
    # initial ages
    agePlusOne[i] ~ dcat(piAGEuncond[1:(max.age+1)]) # where agePlusOne are data
    age[i,1] <- (agePlusOne[i]-1) 
    # I think age+1 is solving the age = zero indexing issue??
    
    # state process
    u[i,1] <- step(age[i,1]-.1) # sets u to zero if age is zero, 1 otherwise
    z[i,1] <- u[i,1]*w[i] # z is the "real" state
    
    # Observation process
    y[i,1] ~ dbern(z[i,1]*p)
    
    # derived stuff
    avail[i,1] <- 1 - u[i,1] # still available -- i.e. not yet recruited
    
    # for occasions > 1     
    for (t in 2:K){
      
      # State process
      u[i,t] ~ dbern(u[i,t-1]*phi[i,t] + avail[i,t-1]*eta[t])   
      logit(phi[i,t]) <- alpha0 
      z[i,t] <- u[i,t]*w[i]
      
      # Age process
      age[i,t] <- age[i,t-1] + max(u[i,1:t]) 
      # ages by one year after recruitment (NIMBLE allows this syntax)
      
      # Observation process
      y[i,t] ~ dbern(z[i,t]*p)
      
      # derived stuff
      avail[i,t] <- 1- max(u[i,1:t]) # still available -- i.e. not yet recruited
    } #t
  } #i
  
  ## Derived population level stuff
  # Annual abundance
  for (t in 1:K){
    N[t] <- sum(z[1:M,t])               
  } #t
  
  # Annual growth rate
  for (t in 1:(K-1)){
    lambda[t] <- N[t+1]/N[t]               
  } #t
  
  # Superpopulation size
  Nsuper <- sum(w[1:M])       
  
})# end model