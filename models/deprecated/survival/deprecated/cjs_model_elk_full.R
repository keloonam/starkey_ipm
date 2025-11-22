cjs_code <- nimbleCode({
  # cjs
  # phi -- year, sex-year (adults), age-year
  # p   -- year, sex-year, random individual effect
  # harvest is a terminal detection; the individual was alive to be harvested
    # (e.g. survived the natural mortality period), but unavailable after
  
  ##### Priors
  # Indexing follows [age, sex, time]
  # [a,x,h,t]
  # age  (a): 1 = calf;   2 = adult
  # sex  (x): 1 = female; 2 = male
  # herd (h): 0 = main;   1 = not in main
  
  # Time effects
  for(t in 1:n_occ){
    s0_ps[t] ~ dunif(0,1)
    s0[t] <- logit(s0_ps[t])
    sm[t] ~ dnorm(0, sd = 50)
    sc[t] ~ dnorm(0, sd = 50)
    p0_ps[t] ~ dunif(0,1)
    p0[t] <- logit(p0_ps[t])
  }
  
  # Constant effects (herds that aren't main)
  sh ~ T(dnorm(0, sd = 50), -10, 10)
  ph ~ T(dnorm(0, sd = 50), -10, 10)
  
  # Individual effects -- mean = 0 because there shouldn't be a direction bias
  sd_ind ~ dunif(0,50)
  for(i in 1:n_ind){
    b_ind[i] ~ dnorm(0, sd = sd_ind)
  } #i
  
  # Probabilities
  for(i in 1:n_ind){
    for(t in f[i]:(n_occ-1)){
      logit(s[i,t])   <- s0[t] + sc[t]*c[i,t] + sm[t]*m[i]*abs((c[i,t]-1))   + sh*h[i,t]
      logit(p[i,t])   <- p0[t] + pm[t]*m[i]   + b_ind[i]     + ph*h[i,t]
    } #t
  } #i

  ##### Likelihood
  # s and p don't have values for occ = 1, p/s-1 shifts the vector
  for(i in 1:n_ind){
    for(t in (f[i]+1):(l[i])){
      
      # state -- did ind i survive?
      z[i,t] ~ dbern(mu1[i,t])
      mu1[i,t] <- s[i,t-1] * z[i,t-1]
      
      # observation -- did we see ind i?
      y[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- p[i,t-1] * z[i,t]
    } #t
  } #i
  
  ##### Derived quantities
  for(t in 1:n_occ){
    logit(survival_af[t]) <- s0[t]
    logit(survival_am[t]) <- s0[t] + sm[t]
    logit(survival_ca[t]) <- s0[t] + sc[t]
  }
}
)