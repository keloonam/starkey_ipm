# PARALLEL NIMBLE ON CGRB
# SGE_Batch -c "R CMD BATCH js_cqls_speed_test.R js_speed_out" -r js_speed_test -q otter -P 3 -M kenneth.loonam@oregonstate.edu


rm(list=ls())

n_cores <- 3

run_js_par <- function(seed){
  
  nb <- 1000
  ni <- 10000
  nc <- 1
  
  load("js_nimble_data_cqls.Rdata")
  
  data      <- jsnm$data
  constants <- jsnm$constants
  params    <- c(
    "survival_af", 
    "survival_am", 
    "survival_ca",
    "detection_f",
    "detection_m",
    "N_super",
    "N",
    "N_females",
    "N_males",
    "N_calves",
    "lambda"
  )
  
  inits <- list(
    s0   = runif(ncol(data$y), 0.5, 0.9),
    p0   = runif(ncol(data$y), 0, 1),
    e0   = runif(1, 0, 1),
    a0   = rep(1/constants$MAX_A, constants$MAX_A),
    s.ma = rnorm(ncol(data$y)),
    s.ca = rnorm(ncol(data$y)),
    s.he = rnorm(ncol(data$y)),
    p.ma = rnorm(ncol(data$y)),
    p.he = rnorm(ncol(data$y))
  )
  
  library(nimble)
  
  js_model_code <- nimbleCode({
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
    a0_i[1:(MAX_A + 1)] <- c((1 - r[1]), r[1] * a0[1:MAX_A])
    
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
    for (i in 1:NAUG){
      # is individual i real?
      w[i] ~ dbern(e0)
      
      # initial ages
      A_P1[i] ~ dcat(a0_i[1:(MAX_A + 1)]) # A_P1 = age plus one
      AGE[i,1] <- (A_P1[i] - 1) 
      # c[i,1] <- c_1[i] # is i a calf?
      # dcat gives integers from 1 to the length of the vector passed
      # age 0 can't be given, so add one to all initial ages.
      # AGE is also data, but has NAs for augmented ind and ind w/o age
      
      # state process
      u[i,1] <- AGE[i,1] > 0 # sets u to zero if age is zero, 1 otherwise
      z[i,1] <- u[i,1] * w[i] # z is the "real" state
      
      # Observation probability
      logit(p[i,1]) <- p.b0[1] + p.bm[1]*M[i] + p.bh[1]*H[i,1]
      
      # Observation process
      y[i,1] ~ dbern(z[i,1] * p[i,1])
      
      # derived stuff
      avail[i,1] <- 1 - u[i,1] # still available -- i.e. not yet recruited
      
      ### Sessions 2:K ###     
      for (t in 2:L[i]){ # 2 to last session i was available
        # Survival probabilities
        logit(s[i,t]) <- s.b0[t] + s.bm[t]*M[i] + s.bh[t]*H[i,t] + s.bc[t]*c[i,t]*avail[i,t-1]
        
        # Detection probabilities
        logit(p[i,t]) <- p.b0[t] + p.bm[t]*M[i] + p.bh[t]*H[i,t]
        
        # State process
        u[i,t] ~ dbern(u[i,t-1] * s[i,t] + avail[i,t-1] * r[t])   
        z[i,t] <- u[i,t] * w[i]
        
        # Age process
        AGE[i,t] <- AGE[i,t-1] + max(u[i,1:t]) 
        # ages by one year after recruitment (NIMBLE allows this syntax?)
        # c[i,t] <- AGE[i,t] * avail[i,t-1] # is i a calf?
        
        # Observation process
        y[i,t] ~ dbern(z[i,t] * p[i,t])
        
        # derived stuff
        avail[i,t] <- 1 - max(u[i,1:t]) # still available -- i.e. not recruited
      } #t
    } #i
    
    ##### Derived/Tracked Values #####
    for(t in 1:K){
      logit(survival_af[t]) <- s.b0[t]
      logit(survival_am[t]) <- s.b0[t] + s.bm[t]
      logit(survival_ca[t]) <- s.b0[t] + s.bc[t]
      logit(detection_f[t]) <- p.b0[t] # also calf detection
      logit(detection_m[t]) <- b.b0[t] + p.bm[t]
    }
    # Annual abundance
    for(t in 1:K){
      N[t]         <- sum(z[1:NAUG,t])
      N_males[t]   <- sum(AM[1:NAUG,t])
      N_females[t] <- sum(AF[1:NAUG,t])
      N_calves[t]  <- sum(CA[1:NAUG,t])
    } #t
    
    # Identities for abundance
    for(i in 1:NAUG){
      for(t in 1:K){
        AM[i,t] <- z[i,t] * M[i] * (AGE[i,t] > 1)
        AF[i,t] <- z[i,t] * (M[i] == 0) * (AGE[i,t] > 1)
        CA[i,t] <- z[i,t] * (AGE[i,t] == 1)
      }
    }
    
    # Annual growth rate
    for(t in 1:(K-1)){
      lambda[t] <- N[t+1]/N[t]
    } #t
    
    # Super-population size
    N_super <- sum(w[1:NAUG])
    
  })# end model
  
  out <- nimbleMCMC(
    code              = js_model_code,
    constants         = constants,
    data              = data,
    monitors          = params,
    inits             = inits,
    niter             = ni,
    nburnin           = nb,
    nchains           = nc,
    progressBar       = T,
    samplesAsCodaMCMC = T,
    summary           = F,
    setSeed           = seed,
    check             = F
  )  
  
  return(out)
}

kl_cluster <- parallel::makeCluster(n_cores, type = "SOCK")

rslt <- parallel::parLapply(
  cl  = kl_cluster, 
  X   = 1:n_cores, 
  fun = run_js_par
)

parallel::stopCluster(kl_cluster)

library(coda)
results <- as.mcmc.list(rslt)
summary(results)
gelman.diag(bear_nimble_results, multivariate = F)

save(results, file = "js_results_2aug2022_speed_test.RData")